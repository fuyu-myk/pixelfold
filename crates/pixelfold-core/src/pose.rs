//! Comparing poses of a ligand by the contacts it makes and the clashes it
//! suffers.
//!
//! Co-folding and docking return many poses of one ligand in one pocket, and
//! what tells a plausible pose from an implausible one is whether it makes the
//! contacts a trusted reference makes. A pose reduces to its fingerprint: the
//! residues the ligand touches, each with the kind of interaction, so three
//! hydrogen bonds to one aspartate count once. Two fingerprints compare by
//! Tanimoto similarity, which charges a pose for the contacts it adds as well as
//! the ones it loses, and by recovery, the share of the reference's contacts it
//! keeps.
//!
//! Residues match by chain, number and insertion code, so a pose must number
//! them as the reference does; [`first_renumbered`] finds where one does not.

use std::collections::{BTreeSet, HashMap, HashSet};

use glam::Vec3;

use crate::fixed_str::FixedStr;
use crate::interactions::{Interaction, InteractionKind, classify, connectivity, topology};
use crate::sasa::get_vdw_radius;
use crate::select::AtomSet;
use crate::spatial::Grid;
use crate::structure::{Atom, ChainId, Protein};

/// The overlap that makes a clash: a ligand atom and a protein atom closer than
/// this fraction of the sum of their van der Waals radii. PoseBusters' default
/// (Buttenschoen, Morris and Deane, Chem. Sci. 15:3130, 2024).
pub const CLASH_RATIO: f32 = 0.75;

/// The farthest a clash can reach: the ratio times twice the largest radius
/// the van der Waals table holds (iodine, 1.98 A).
const CLASH_REACH: f32 = 3.0;

/// A residue's identity: chain, number and insertion code.
pub type ResidueKey = (ChainId, i32, Option<char>);

/// One element of a fingerprint: a residue the ligand touches, and how.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Contact {
    pub residue: ResidueKey,
    pub residue_name: FixedStr<6>,
    pub kind: InteractionKind,
}

/// The contacts a ligand makes, as a set.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Fingerprint {
    contacts: BTreeSet<Contact>,
}

impl FromIterator<Contact> for Fingerprint {
    fn from_iter<I: IntoIterator<Item = Contact>>(contacts: I) -> Self {
        Self {
            contacts: contacts.into_iter().collect(),
        }
    }
}

impl Fingerprint {
    /// How many contacts there are.
    pub fn len(&self) -> usize {
        self.contacts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.contacts.is_empty()
    }

    /// The contacts, ordered by residue and then kind.
    pub fn contacts(&self) -> impl Iterator<Item = &Contact> {
        self.contacts.iter()
    }

    /// How many contacts this fingerprint and `other` both hold.
    pub fn shared(&self, other: &Fingerprint) -> usize {
        self.contacts.intersection(&other.contacts).count()
    }

    /// The Tanimoto similarity to `other`: shared contacts over all contacts
    /// either holds. Two empty fingerprints are identical.
    pub fn similarity(&self, other: &Fingerprint) -> f32 {
        let either = self.contacts.union(&other.contacts).count();
        if either == 0 {
            return 1.0;
        }

        self.shared(other) as f32 / either as f32
    }

    /// The share of `reference`'s contacts this fingerprint keeps, 1 when the
    /// reference has none.
    pub fn recovery(&self, reference: &Fingerprint) -> f32 {
        if reference.is_empty() {
            return 1.0;
        }

        self.shared(reference) as f32 / reference.len() as f32
    }

    /// The contacts of `reference` this fingerprint lacks.
    pub fn lost<'a>(&'a self, reference: &'a Fingerprint) -> impl Iterator<Item = &'a Contact> {
        reference.contacts.difference(&self.contacts)
    }
}

/// The fingerprint of the `ligand` atoms in `interactions`.
///
/// An interaction is a contact when exactly one side is ligand; the other side
/// names the residue. Water is left out: a water bridge, or a hydrogen bond to a
/// water, is not a contact, since a predicted pose carries no water and would
/// lose every one of them.
pub fn fingerprint(
    protein: &Protein,
    interactions: &[Interaction],
    ligand: &AtomSet,
) -> Fingerprint {
    let is_ligand = |atoms: &[usize]| atoms.iter().any(|&i| ligand.contains(i));

    interactions
        .iter()
        .filter(|found| found.kind != InteractionKind::WaterBridge)
        .filter_map(|found| {
            let partner = match (is_ligand(&found.atoms_a), is_ligand(&found.atoms_b)) {
                (true, false) => &found.atoms_b,
                (false, true) => &found.atoms_a,
                _ => return None,
            };
            let atom = &protein.atoms[partner[0]];
            if classify::is_water(atom.residue_name.as_str()) {
                return None;
            }

            Some(Contact {
                residue: residue_key(atom),
                residue_name: atom.residue_name,
                kind: found.kind,
            })
        })
        .collect()
}

/// A protein atom's identity across structures: residue and atom name.
type AtomIdentity = (ChainId, i32, Option<char>, FixedStr<4>);

/// A ligand atom's identity across structures: component and atom name. A pose
/// may put its ligand in another chain, under another number, than the
/// reference does, as co-folding output usually does.
type LigandAtom = (FixedStr<6>, FixedStr<4>);

/// The ligand-protein atom pairs a reference bonds covalently, and the pairs one
/// bond from each, by identity.
///
/// A covalent ligand sits at bonding distance from its anchor, well inside the
/// clash threshold, as do the atoms one bond further. The attachment is read off
/// the reference rather than each pose because in a pose a severe clash and a
/// bond are the same distance, and so is the bond graph that finds the pairs one
/// bond further: a pose may not define its ligand at all.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Attachments {
    pairs: Vec<(LigandAtom, AtomIdentity)>,
}

/// The covalent attachments of the `ligand` atoms in `reference`.
pub fn attachments(reference: &Protein, ligand: &AtomSet) -> Attachments {
    let bonds = connectivity::bonds(reference);
    let atoms = &reference.atoms;
    let inside = |atom: usize| {
        bonds
            .of(atom)
            .iter()
            .copied()
            .filter(|&other| ligand.contains(other))
    };
    let outside = |atom: usize| {
        bonds
            .of(atom)
            .iter()
            .copied()
            .filter(|&other| !ligand.contains(other))
    };

    let mut pairs = BTreeSet::new();
    for atom in ligand.iter() {
        for anchor in outside(atom) {
            for near in std::iter::once(atom).chain(inside(atom)) {
                pairs.insert((ligand_atom(&atoms[near]), identity(&atoms[anchor])));
            }
            for near in outside(anchor) {
                pairs.insert((ligand_atom(&atoms[atom]), identity(&atoms[near])));
            }
        }
    }

    Attachments {
        pairs: pairs.into_iter().collect(),
    }
}

/// How many ligand-protein heavy-atom pairs clash.
///
/// Hydrogens, water and metal ions are left out, a metal because coordination
/// sits inside the threshold, and so are the pairs `attached` names.
pub fn clashes(protein: &Protein, ligand: &AtomSet, attached: &Attachments) -> usize {
    let counts = |atom: &Atom| {
        !is_hydrogen(atom)
            && !classify::is_water(atom.residue_name.as_str())
            && !classify::is_metal(atom)
    };
    let excluded = attached_pairs(protein, ligand, attached);
    let positions: Vec<Vec3> = protein.atoms.iter().map(|atom| atom.position).collect();
    let grid = Grid::build(&positions, CLASH_REACH);

    let mut found = 0;
    for atom in ligand.iter().filter(|&atom| counts(&protein.atoms[atom])) {
        let radius = get_vdw_radius(protein.atoms[atom].element.as_str());
        grid.for_each_within(positions[atom], CLASH_REACH, |other| {
            let partner = &protein.atoms[other];
            if ligand.contains(other) || !counts(partner) || excluded.contains(&(atom, other)) {
                return;
            }

            let limit = CLASH_RATIO * (radius + get_vdw_radius(partner.element.as_str()));
            if (positions[atom] - positions[other]).length() < limit {
                found += 1;
            }
        });
    }

    found
}

/// The ligand-protein index pairs `attached` names in `protein`, every
/// conformer of an atom included.
fn attached_pairs(
    protein: &Protein,
    ligand: &AtomSet,
    attached: &Attachments,
) -> HashSet<(usize, usize)> {
    if attached.pairs.is_empty() {
        return HashSet::new();
    }

    let mut ligand_atoms: HashMap<LigandAtom, Vec<usize>> = HashMap::new();
    let mut others: HashMap<AtomIdentity, Vec<usize>> = HashMap::new();
    for (index, atom) in protein.atoms.iter().enumerate() {
        if ligand.contains(index) {
            ligand_atoms
                .entry(ligand_atom(atom))
                .or_default()
                .push(index);
        } else {
            others.entry(identity(atom)).or_default().push(index);
        }
    }

    attached
        .pairs
        .iter()
        .filter_map(|(near, anchor)| Some((ligand_atoms.get(near)?, others.get(anchor)?)))
        .flat_map(|(nears, anchors)| {
            nears
                .iter()
                .flat_map(move |&a| anchors.iter().map(move |&b| (a, b)))
        })
        .collect()
}

/// The ligand residues whose chemistry the structure does not define for every
/// heavy atom, each once with its name. Their contacts are incomplete, since an
/// atom without a definition has no bonds to classify it by.
pub fn undefined_ligand_residues(
    protein: &Protein,
    ligand: &AtomSet,
) -> Vec<(ResidueKey, FixedStr<6>)> {
    let undefined: BTreeSet<(ResidueKey, FixedStr<6>)> = ligand
        .iter()
        .map(|atom| &protein.atoms[atom])
        .filter(|atom| {
            let residue = atom.residue_name.as_str();
            let defined = topology::is_standard_residue(residue)
                || protein
                    .components
                    .get(residue)
                    .is_some_and(|component| component.element_of(atom.name.as_str()).is_some());

            !defined
                && !is_hydrogen(atom)
                && !classify::is_water(residue)
                && !classify::is_metal(atom)
        })
        .map(|atom| (residue_key(atom), atom.residue_name))
        .collect();

    undefined.into_iter().collect()
}

/// A residue the reference makes a contact with that a pose numbers
/// differently: another residue under the same number, or none.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Renumbered {
    pub residue: ResidueKey,
    pub expected: FixedStr<6>,
    pub found: Option<FixedStr<6>>,
}

/// The first residue of `reference` that `pose` does not hold under the same
/// number and name.
pub fn first_renumbered(reference: &Fingerprint, pose: &Protein) -> Option<Renumbered> {
    let names: HashMap<ResidueKey, FixedStr<6>> = pose
        .atoms
        .iter()
        .map(|atom| (residue_key(atom), atom.residue_name))
        .collect();

    reference.contacts().find_map(|contact| {
        let found = names.get(&contact.residue).copied();
        (found != Some(contact.residue_name)).then_some(Renumbered {
            residue: contact.residue,
            expected: contact.residue_name,
            found,
        })
    })
}

fn residue_key(atom: &Atom) -> ResidueKey {
    (atom.chain_id, atom.residue_seq, atom.insertion_code)
}

fn identity(atom: &Atom) -> AtomIdentity {
    (
        atom.chain_id,
        atom.residue_seq,
        atom.insertion_code,
        atom.name,
    )
}

fn ligand_atom(atom: &Atom) -> LigandAtom {
    (atom.residue_name, atom.name)
}

fn is_hydrogen(atom: &Atom) -> bool {
    matches!(
        FixedStr::<2>::uppercased(atom.element.as_str()).as_str(),
        "H" | "D"
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::components::Dictionary;
    use crate::select;
    use crate::structure::{Molecule, SecondaryStructure};
    use InteractionKind::{HydrogenBond, Hydrophobic, SaltBridge, WaterBridge};
    use glam::Vec3;

    fn atom(seq: i32, residue: &str, name: &str, element: &str, position: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(element),
            residue_name: FixedStr::new(residue),
            residue_seq: seq,
            chain_id: FixedStr::new("A"),
            is_hetatm: matches!(residue, "LIG" | "HOH" | "ZN"),
            molecule: Molecule::Unknown,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    fn protein(atoms: Vec<Atom>) -> Protein {
        Protein {
            atoms,
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        }
    }

    fn ligand_of(structure: &Protein) -> AtomSet {
        select("resn LIG", structure).expect("a valid query")
    }

    fn between(kind: InteractionKind, a: usize, b: usize) -> Interaction {
        Interaction {
            kind,
            atoms_a: vec![a],
            atoms_b: vec![b],
            distance: 3.0,
            angle: None,
            bridge: None,
        }
    }

    fn contact(seq: i32, residue: &str, kind: InteractionKind) -> Contact {
        Contact {
            residue: (FixedStr::new("A"), seq, None),
            residue_name: FixedStr::new(residue),
            kind,
        }
    }

    /// A ligand (atoms 0 and 1), a lysine (2), an aspartate (3), a second
    /// ligand residue (4) and a water (5). A fingerprint reads no positions.
    fn site() -> Protein {
        protein(vec![
            atom(900, "LIG", "O1", "O", Vec3::ZERO),
            atom(900, "LIG", "C2", "C", Vec3::ZERO),
            atom(10, "LYS", "NZ", "N", Vec3::ZERO),
            atom(20, "ASP", "OD1", "O", Vec3::ZERO),
            atom(901, "LIG", "C1", "C", Vec3::ZERO),
            atom(950, "HOH", "O", "O", Vec3::ZERO),
        ])
    }

    #[test]
    fn a_fingerprint_holds_each_residue_and_kind_the_ligand_makes_once() {
        let structure = site();
        let found = [
            between(HydrogenBond, 2, 0),
            between(HydrogenBond, 2, 1), // a second bond to the same lysine
            between(Hydrophobic, 1, 3),
            between(SaltBridge, 2, 3), // touches no ligand atom
        ];

        let expected: Fingerprint = [
            contact(10, "LYS", HydrogenBond),
            contact(20, "ASP", Hydrophobic),
        ]
        .into_iter()
        .collect();
        assert_eq!(
            fingerprint(&structure, &found, &ligand_of(&structure)),
            expected
        );
    }

    #[test]
    fn water_and_the_ligand_itself_are_not_contacts() {
        let structure = site();
        let found = [
            between(WaterBridge, 2, 0),
            between(HydrogenBond, 5, 0), // to a crystal water
            between(Hydrophobic, 1, 4),  // between the two ligand residues
        ];

        assert!(fingerprint(&structure, &found, &ligand_of(&structure)).is_empty());
    }

    #[test]
    fn similarity_charges_added_and_lost_contacts_and_recovery_only_lost_ones() {
        let [a, b, c, d, e] = [10, 20, 30, 40, 50].map(|seq| contact(seq, "ALA", Hydrophobic));
        let reference: Fingerprint = [a, b, c, d].into_iter().collect();
        let pose: Fingerprint = [a, b, e].into_iter().collect();

        assert_eq!(pose.shared(&reference), 2);
        assert!((pose.similarity(&reference) - 2.0 / 5.0).abs() < 1e-6);
        assert!((pose.recovery(&reference) - 2.0 / 4.0).abs() < 1e-6);
        assert_eq!(
            pose.lost(&reference).copied().collect::<Vec<_>>(),
            vec![c, d]
        );
        assert_eq!(reference.similarity(&reference), 1.0);
        assert_eq!(
            Fingerprint::default().similarity(&Fingerprint::default()),
            1.0
        );
    }

    /// A ligand carbon at the origin and one other atom `distance` along x.
    fn pair(residue: &str, name: &str, element: &str, distance: f32) -> Protein {
        protein(vec![
            atom(900, "LIG", "C1", "C", Vec3::ZERO),
            atom(10, residue, name, element, Vec3::X * distance),
        ])
    }

    #[test]
    fn a_clash_is_a_pair_closer_than_three_quarters_of_its_radii() {
        // Two carbons: 0.75 x (1.70 + 1.70) = 2.55 A. At 1.4 A the pair is at
        // bonding distance, which in a pose is a clash, not a bond.
        for (distance, expected) in [(2.5, 1), (2.6, 0), (1.4, 1)] {
            let structure = pair("LEU", "CD1", "C", distance);
            let found = clashes(&structure, &ligand_of(&structure), &Attachments::default());

            assert_eq!(found, expected, "at {distance} A");
        }
    }

    #[test]
    fn hydrogens_water_and_metals_never_clash() {
        for (residue, name, element) in
            [("LEU", "HD11", "H"), ("HOH", "O", "O"), ("ZN", "ZN", "ZN")]
        {
            let structure = pair(residue, name, element, 1.5);
            let found = clashes(&structure, &ligand_of(&structure), &Attachments::default());

            assert_eq!(found, 0, "{residue} {name}");
        }
    }

    const LIG: &str = "\
data_LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_ordinal
LIG C1 C N 1
LIG C2 C N 2
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_ordinal
LIG C1 C2 sing N 1
#
";

    /// A two-carbon ligand bonded through C1 to a serine OG, as a covalent
    /// inhibitor is. C1-CB and C2-OG, one bond further, sit inside the clash
    /// threshold too; C2-CB, two bonds further, does not.
    fn attached() -> Protein {
        let mut structure = protein(vec![
            atom(900, "LIG", "C1", "C", Vec3::ZERO),
            atom(900, "LIG", "C2", "C", Vec3::new(-0.495, 1.437, 0.0)),
            atom(195, "SER", "OG", "O", Vec3::new(1.42, 0.0, 0.0)),
            atom(195, "SER", "CB", "C", Vec3::new(2.024, -1.296, 0.0)),
        ]);
        structure.components = Dictionary::parse(LIG);
        structure
    }

    #[test]
    fn a_covalent_attachment_the_reference_makes_is_not_a_clash() {
        let structure = attached();
        let ligand = ligand_of(&structure);

        assert_eq!(clashes(&structure, &ligand, &Attachments::default()), 3);
        assert_eq!(
            clashes(&structure, &ligand, &attachments(&structure, &ligand)),
            0
        );
    }

    #[test]
    fn an_attachment_carries_to_a_pose_by_atom_names() {
        let reference = attached();
        let attached_in_reference = attachments(&reference, &ligand_of(&reference));

        // The same pose with its atoms in another order and its ligand in a chain
        // of its own, numbered from 1, as co-folding output writes it.
        let mut pose = attached();
        pose.atoms.rotate_left(2);
        for atom in pose
            .atoms
            .iter_mut()
            .filter(|atom| atom.residue_name == "LIG")
        {
            atom.chain_id = FixedStr::new("B");
            atom.residue_seq = 1;
        }

        assert_eq!(clashes(&pose, &ligand_of(&pose), &attached_in_reference), 0);
    }

    #[test]
    fn an_attachment_holds_in_a_pose_without_ligand_chemistry() {
        let reference = attached();
        let attached_in_reference = attachments(&reference, &ligand_of(&reference));

        let mut pose = attached();
        pose.components = Dictionary::default();

        assert_eq!(clashes(&pose, &ligand_of(&pose), &attached_in_reference), 0);
    }

    #[test]
    fn a_close_pair_beside_an_attachment_still_clashes() {
        let reference = attached();
        let attached_in_reference = attachments(&reference, &ligand_of(&reference));

        // A leucine carbon jammed 1.43 A from C2, which the reference bonds to
        // nothing outside the ligand.
        let mut pose = attached();
        pose.atoms
            .push(atom(50, "LEU", "CD1", "C", Vec3::new(-1.495, 2.337, 0.5)));

        assert_eq!(clashes(&pose, &ligand_of(&pose), &attached_in_reference), 1);
    }

    #[test]
    fn a_ligand_atom_its_definition_does_not_name_leaves_the_residue_undefined() {
        let mut structure = attached();
        assert!(undefined_ligand_residues(&structure, &ligand_of(&structure)).is_empty());

        structure
            .atoms
            .push(atom(900, "LIG", "N3", "N", Vec3::new(0.0, 0.0, 4.0)));
        let undefined = vec![((FixedStr::new("A"), 900, None), FixedStr::new("LIG"))];
        assert_eq!(
            undefined_ligand_residues(&structure, &ligand_of(&structure)),
            undefined
        );

        structure.components = Dictionary::default();
        assert_eq!(
            undefined_ligand_residues(&structure, &ligand_of(&structure)),
            undefined
        );
    }

    #[test]
    fn a_pose_that_numbers_a_contact_residue_differently_is_caught() {
        let structure = site();
        let fingerprint_of = |seq, residue| -> Fingerprint {
            [contact(seq, residue, HydrogenBond)].into_iter().collect()
        };

        assert_eq!(
            first_renumbered(&fingerprint_of(10, "LYS"), &structure),
            None
        );
        assert_eq!(
            first_renumbered(&fingerprint_of(20, "LYS"), &structure),
            Some(Renumbered {
                residue: (FixedStr::new("A"), 20, None),
                expected: FixedStr::new("LYS"),
                found: Some(FixedStr::new("ASP")),
            })
        );
        assert_eq!(
            first_renumbered(&fingerprint_of(30, "GLY"), &structure).map(|missing| missing.found),
            Some(None)
        );
    }
}
