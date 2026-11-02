//! Checks of a ligand pose that need no reference structure, as PoseBusters
//! makes them of a docked pose (Buttenschoen, Morris and Deane, Chem. Sci.
//! 15:3130, 2024): its geometry against the bounds RDKit's distance geometry
//! gives its bonds, and its stereochemistry against the ideal conformer of its
//! chemical component.
//!
//! The geometry needs the ligand's bond table: a structure file's own
//! component definitions, a dictionary entry borrowed for them, or an SD
//! record's bonds. The stereochemistry needs more, since a pose's own
//! coordinates say only which way its centres turn and not which way they
//! should: the Chemical Component Dictionary's ideal model of the component,
//! whose R and S centres and E and Z double bonds it labels. A centre or a
//! double bond the dictionary leaves unlabelled is unspecified and goes
//! unchecked.

use std::collections::{BTreeMap, HashMap};

use crate::components::{BondOrder, Component, Dictionary};
use crate::fixed_str::FixedStr;
use crate::ideal;
use crate::interactions::classify::is_metal;
use crate::interactions::connectivity;
use crate::molecule::{Molecule, inferred_charge};
use crate::numbering::{ResidueKey, residue_key};
use crate::pose::is_hydrogen;
use crate::rmsd::{AtomId, RmsdReference};
use crate::select::AtomSet;
use crate::structure::{Atom, Protein};

pub use crate::geometry::Geometry;

/// The share of a component's bonds between a residue's atoms that must lie
/// at bonding distance for the residue to be taken for that component: enough
/// that atoms named alike by chance are not, few enough that a badly
/// distorted pose of the component is, and fails its geometry instead.
const BONDED_SHARE: f32 = 0.5;

/// What of a pose's ligand turns the other way from its chemical component's
/// ideal conformer.
#[derive(Clone, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct Stereo {
    /// The labelled stereocentres the pose inverts.
    pub inverted: Vec<AtomId>,
    /// The labelled E or Z double bonds the pose turns the other way, each by
    /// its two atoms.
    pub flipped: Vec<[AtomId; 2]>,
    /// The ligand's components with no ideal conformer its atoms match, whose
    /// stereochemistry goes unchecked.
    pub unchecked: Vec<FixedStr<6>>,
    /// How many of the ligand's residues were checked.
    pub checked: usize,
}

/// The `ligand` atoms of `pose` checked as PoseBusters checks a docked pose's
/// own geometry: its bond lengths, its bond angles, its atoms clashing with
/// one another, and the flatness of its aromatic rings and its double bonds
/// between carbons, against the bounds RDKit gives the bonds `pose`'s
/// definitions say it has. Residues of the ligand bond one another where
/// their atoms lie at bonding distance. `Err` names a residue no definition
/// covers, whose chemistry is unknown.
pub fn geometry(pose: &Protein, ligand: &AtomSet) -> Result<Geometry, FixedStr<6>> {
    let heavy: Vec<usize> = ligand
        .iter()
        .filter(|&atom| !is_hydrogen(&pose.atoms[atom]))
        .collect();
    let molecule = molecule(pose, &heavy)?;
    let positions: Vec<glam::Vec3> = heavy
        .iter()
        .map(|&atom| pose.atoms[atom].position)
        .collect();

    Ok(ideal::template(&molecule, &positions).check(&positions, true))
}

/// The `heavy` atoms of `pose` as a molecule, each residue's bonds, orders,
/// hydrogens and charges from its definition, a bond between residues
/// wherever two atoms lie at bonding distance, each taking one of its atom's
/// hydrogens. A definition without charges, as a structure file's own are,
/// leaves them to the valences.
fn molecule(pose: &Protein, heavy: &[usize]) -> Result<Molecule, FixedStr<6>> {
    let atoms: Vec<&Atom> = heavy.iter().map(|&atom| &pose.atoms[atom]).collect();
    let definitions: Vec<&Component> = atoms
        .iter()
        .map(|atom| {
            pose.components
                .get(atom.residue_name.as_str())
                .filter(|component| component.element_of(atom.name.as_str()).is_some())
                .ok_or(atom.residue_name)
        })
        .collect::<Result<_, _>>()?;
    let mut index: HashMap<(ResidueKey, FixedStr<4>), usize> = HashMap::new();
    for (at, atom) in atoms.iter().enumerate() {
        index.entry((residue_key(atom), atom.name)).or_insert(at);
    }

    let mut bonds: Vec<(usize, usize, BondOrder)> = Vec::new();
    for (at, atom) in atoms.iter().enumerate() {
        for (other, order) in definitions[at].bonded(atom.name.as_str()) {
            if let Some(&partner) = index.get(&(residue_key(atom), FixedStr::new(other)))
                && partner > at
            {
                bonds.push((at, partner, order));
            }
        }
    }
    let mut links = vec![0usize; atoms.len()];
    for (a, first) in atoms.iter().enumerate() {
        for (b, second) in atoms.iter().enumerate().skip(a + 1) {
            if residue_key(first) != residue_key(second)
                && connectivity::at_bonding_distance(first, second)
            {
                bonds.push((a, b, BondOrder::Single));
                links[a] += 1;
                links[b] += 1;
            }
        }
    }

    let hydrogens: Vec<u8> = atoms
        .iter()
        .enumerate()
        .map(|(at, atom)| {
            let defined = definitions[at].hydrogen_count(atom.name.as_str());
            defined.saturating_sub(links[at]).min(usize::from(u8::MAX)) as u8
        })
        .collect();
    let mut valences: Vec<u8> = hydrogens.clone();
    for &(a, b, order) in &bonds {
        let count = match order {
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            _ => 1,
        };
        valences[a] = valences[a].saturating_add(count);
        valences[b] = valences[b].saturating_add(count);
    }
    let charges: Vec<i8> = atoms
        .iter()
        .enumerate()
        .map(|(at, atom)| {
            definitions[at]
                .charge(atom.name.as_str())
                .unwrap_or_else(|| inferred_charge(atom.element.as_str(), valences[at]))
        })
        .collect();

    Ok(Molecule::new(
        atoms
            .iter()
            .map(|atom| FixedStr::uppercased(atom.element.as_str()))
            .collect(),
        &bonds,
        hydrogens,
        charges,
    ))
}

/// The `ligand` residues of `pose`, each with its atoms, in order.
fn residues(pose: &Protein, ligand: &AtomSet) -> BTreeMap<ResidueKey, (FixedStr<6>, Vec<usize>)> {
    let mut residues: BTreeMap<ResidueKey, (FixedStr<6>, Vec<usize>)> = BTreeMap::new();
    for atom in ligand.iter() {
        let found = &pose.atoms[atom];
        residues
            .entry(residue_key(found))
            .or_insert_with(|| (found.residue_name, Vec::new()))
            .1
            .push(atom);
    }

    residues
}

/// Whether the `atoms` of one residue of `pose` are `component`'s: each heavy
/// atom named as one of its heavy atoms, of the same element, and at least
/// half of its bonds between them at bonding distance, bonds to a metal left
/// out since no covalent radius bonds one. A component's code says little
/// alone: co-folding output names a ligand given as a SMILES string `LIG`,
/// which the dictionary defines as another molecule.
fn matches(pose: &Protein, atoms: &[usize], component: &Component) -> bool {
    let mut by_name: HashMap<&str, &Atom> = HashMap::new();
    for atom in atoms
        .iter()
        .map(|&atom| &pose.atoms[atom])
        .filter(|atom| !is_hydrogen(atom))
    {
        let same = component
            .element_of(atom.name.as_str())
            .is_some_and(|element| element.eq_ignore_ascii_case(atom.element.as_str()));
        if !same || by_name.insert(atom.name.as_str(), atom).is_some() {
            return false;
        }
    }
    let (mut bonds, mut bonded) = (0usize, 0usize);
    for (&name, &atom) in &by_name {
        for other in component.neighbours(name) {
            if let Some(&partner) = by_name.get(other)
                && name < other
                && !is_metal(atom)
                && !is_metal(partner)
            {
                bonds += 1;
                bonded += usize::from(connectivity::at_bonding_distance(atom, partner));
            }
        }
    }

    !by_name.is_empty() && bonded as f32 >= BONDED_SHARE * bonds as f32
}

/// The entries of `dictionary` that every residue of their component among
/// the `ligand` atoms of `pose` matches, as [`stereo`] would check them, and
/// the codes of the components it holds whose residues do not.
pub fn matching(
    pose: &Protein,
    ligand: &AtomSet,
    dictionary: &Dictionary,
) -> (Dictionary, Vec<FixedStr<6>>) {
    let mut fits: BTreeMap<FixedStr<6>, bool> = BTreeMap::new();
    for (code, atoms) in residues(pose, ligand).into_values() {
        if let Some(component) = dictionary.get(code.as_str()) {
            let fit = fits.entry(code).or_insert(true);
            *fit = *fit && matches(pose, &atoms, component);
        }
    }

    let mut matched = Dictionary::default();
    let mut refused = Vec::new();
    for (code, fit) in fits {
        match dictionary.get(code.as_str()) {
            Some(component) if fit => matched = matched.with(code.as_str(), component.clone()),
            _ => refused.push(code),
        }
    }

    (matched, refused)
}

/// The labelled stereocentres and E or Z double bonds that the `ligand` atoms
/// of `pose` turn the other way from the ideal conformers `ideals` holds,
/// residue by residue, each residue's atoms paired with its component's by
/// name. A residue whose component has no conformer there, or whose atoms it
/// lacks, goes unchecked.
pub fn stereo(pose: &Protein, ligand: &AtomSet, ideals: &Dictionary) -> Stereo {
    let mut stereo = Stereo::default();
    for (code, atoms) in residues(pose, ligand).into_values() {
        let checked = ideals
            .get(code.as_str())
            .and_then(|component| residue_stereo(pose, &atoms, code, component));
        match checked {
            Some((inverted, flipped)) => {
                stereo.inverted.extend(inverted);
                stereo.flipped.extend(flipped);
                stereo.checked += 1;
            }
            None if !stereo.unchecked.contains(&code) => stereo.unchecked.push(code),
            None => {}
        }
    }

    stereo
}

/// Component `code`'s ideal conformer as a reference ligand, which an SD
/// record's atoms can be named after ([`crate::sdf::Record::defined`]):
/// `None` when the dictionary gives some heavy atom no ideal position.
pub fn ideal_reference(code: &str, component: &Component) -> Option<RmsdReference> {
    let atoms: Vec<Atom> = component
        .atoms_in_order()
        .enumerate()
        .map(|(at, name)| {
            Some(Atom {
                serial: at as u32 + 1,
                name: FixedStr::new(name),
                element: FixedStr::new(component.element_of(name)?),
                residue_name: FixedStr::new(code),
                residue_seq: 1,
                chain_id: crate::structure::ChainId::new("L"),
                is_hetatm: true,
                molecule: crate::structure::Molecule::NonPolymer,
                altloc: None,
                occupancy: 1.0,
                insertion_code: None,
                model_number: 1,
                position: component.ideal(name)?,
                b_factor: 0.0,
                secondary_structure: crate::structure::SecondaryStructure::Coil,
            })
        })
        .collect::<Option<_>>()?;
    let mut all = AtomSet::empty(atoms.len());
    for atom in 0..atoms.len() {
        all.insert(atom);
    }
    let conformer = Protein {
        atoms,
        title: String::new(),
        hbonds: Vec::new(),
        assembly: None,
        components: Dictionary::default().with(code, component.clone()),
    };

    Some(RmsdReference::new(&conformer, &all))
}

/// One residue's stereochemistry against `component`'s ideal conformer, the
/// conformer made of the atoms the residue has (a residue linked to others
/// lacks the leaving atoms the conformer holds, and its centre at the link is
/// not checked).
fn residue_stereo(
    pose: &Protein,
    atoms: &[usize],
    code: FixedStr<6>,
    component: &Component,
) -> Option<(Vec<AtomId>, Vec<[AtomId; 2]>)> {
    if !matches(pose, atoms, component) {
        return None;
    }
    let ideal_atoms: Vec<Atom> = atoms
        .iter()
        .map(|&atom| &pose.atoms[atom])
        .filter(|atom| !is_hydrogen(atom))
        .map(|atom| {
            Some(Atom {
                position: component.ideal(atom.name.as_str())?,
                ..atom.clone()
            })
        })
        .collect::<Option<_>>()?;
    let mut all = AtomSet::empty(ideal_atoms.len());
    for atom in 0..ideal_atoms.len() {
        all.insert(atom);
    }
    let conformer = Protein {
        atoms: ideal_atoms,
        title: String::new(),
        hbonds: Vec::new(),
        assembly: None,
        components: Dictionary::default().with(code.as_str(), component.clone()),
    };
    let reference = RmsdReference::new(&conformer, &all).labelled(
        |name| component.is_stereocentre(name),
        |first, second| component.has_stereo_bond(first, second),
    );
    let mut residue = AtomSet::empty(pose.atoms.len());
    for &atom in atoms {
        residue.insert(atom);
    }

    reference.stereo(pose, &residue)
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use glam::Vec3;

    use super::*;
    use crate::structure::{ChainId, Molecule as Kind, SecondaryStructure};

    const NAG: &str = include_str!("../tests/data/NAG.cif");
    const GPP: &str = include_str!("../tests/data/GPP.cif");

    /// The dictionary entry `code` of `text` as a pose: its heavy atoms at
    /// their ideal coordinates, one residue, defined by the entry.
    fn ideal_pose(text: &str, code: &str) -> Protein {
        let dictionary = Dictionary::parse(text);
        let component = dictionary.get(code).unwrap();
        let atoms = component
            .atoms_in_order()
            .enumerate()
            .map(|(at, name)| Atom {
                serial: at as u32 + 1,
                name: FixedStr::new(name),
                element: FixedStr::new(component.element_of(name).unwrap()),
                residue_name: FixedStr::new(code),
                residue_seq: 1,
                chain_id: ChainId::new("A"),
                is_hetatm: true,
                molecule: Kind::NonPolymer,
                altloc: None,
                occupancy: 1.0,
                insertion_code: None,
                model_number: 1,
                position: component.ideal(name).unwrap(),
                b_factor: 0.0,
                secondary_structure: SecondaryStructure::Coil,
            })
            .collect();

        Protein {
            atoms,
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: dictionary,
        }
    }

    fn everything(pose: &Protein) -> AtomSet {
        let mut all = AtomSet::empty(pose.atoms.len());
        for atom in 0..pose.atoms.len() {
            all.insert(atom);
        }
        all
    }

    fn named(pose: &Protein, name: &str) -> usize {
        pose.atoms
            .iter()
            .position(|atom| atom.name.as_str() == name)
            .unwrap()
    }

    fn names(ids: &[AtomId]) -> Vec<&str> {
        let mut names: Vec<&str> = ids.iter().map(|(_, _, name)| name.as_str()).collect();
        names.sort_unstable();
        names
    }

    /// The pose with every atom on `start`'s side of its bond to `across`
    /// turned half a turn about that bond.
    fn turned(pose: &Protein, start: &str, across: &str) -> Protein {
        let component = pose
            .components
            .get(pose.atoms[0].residue_name.as_str())
            .unwrap();
        let mut side: HashSet<&str> = HashSet::from([start]);
        let mut queue = vec![start];
        while let Some(name) = queue.pop() {
            for other in component.neighbours(name) {
                if other != across && side.insert(other) {
                    queue.push(other);
                }
            }
        }
        let (origin, end) = (
            pose.atoms[named(pose, start)].position,
            pose.atoms[named(pose, across)].position,
        );
        let turn = glam::Quat::from_axis_angle((end - origin).normalize(), std::f32::consts::PI);
        let mut turned = pose.clone();
        for atom in &mut turned.atoms {
            if atom.name.as_str() != start && side.contains(atom.name.as_str()) {
                atom.position = origin + turn * (atom.position - origin);
            }
        }
        turned
    }

    #[test]
    fn a_dictionary_ligand_at_its_ideal_coordinates_passes_every_check() {
        for (text, code) in [(NAG, "NAG"), (GPP, "GPP")] {
            let pose = ideal_pose(text, code);
            let ligand = everything(&pose);

            let geometry = geometry(&pose, &ligand).unwrap();
            assert_eq!(
                (geometry.bonds, geometry.angles, geometry.clashes),
                (0, 0, 0),
                "{code}"
            );
            assert_eq!(geometry.unflat, Some(0), "{code}");

            let stereo = stereo(&pose, &ligand, &pose.components);
            assert!(
                stereo.inverted.is_empty() && stereo.flipped.is_empty(),
                "{code}: {stereo:?}"
            );
            assert!(stereo.unchecked.is_empty(), "{code}");
        }
    }

    #[test]
    fn the_mirror_image_of_n_acetylglucosamine_inverts_its_five_labelled_centres() {
        let pose = ideal_pose(NAG, "NAG");
        let mut mirrored = pose.clone();
        for atom in &mut mirrored.atoms {
            atom.position.x = -atom.position.x;
        }
        let ligand = everything(&mirrored);

        let stereo = stereo(&mirrored, &ligand, &pose.components);
        assert_eq!(names(&stereo.inverted), ["C1", "C2", "C3", "C4", "C5"]);
        // A mirror image keeps its bonds and angles.
        let geometry = geometry(&mirrored, &ligand).unwrap();
        assert_eq!((geometry.bonds, geometry.angles), (0, 0));
    }

    #[test]
    fn geranyl_diphosphate_turned_about_its_e_bond_is_flipped() {
        let pose = ideal_pose(GPP, "GPP");
        let flipped = turned(&pose, "C2", "C3");
        let ligand = everything(&flipped);

        let stereo = stereo(&flipped, &ligand, &pose.components);
        let bonds: Vec<[String; 2]> = stereo
            .flipped
            .iter()
            .map(|ends| ends.map(|(_, _, name)| name.to_string()))
            .collect();
        assert_eq!(bonds, [["C2", "C3"]]);
        assert!(stereo.inverted.is_empty());
    }

    #[test]
    fn a_stretched_bond_fails_and_an_undefined_residue_is_named() {
        let mut pose = ideal_pose(NAG, "NAG");
        let oxygen = named(&pose, "O7");
        pose.atoms[oxygen].position += Vec3::new(0.0, 0.0, 0.8);
        let ligand = everything(&pose);
        assert!(geometry(&pose, &ligand).unwrap().bonds >= 1);

        let undefined = Protein {
            components: Dictionary::default(),
            ..pose
        };
        assert_eq!(geometry(&undefined, &ligand), Err(FixedStr::new("NAG")));
    }

    #[test]
    fn a_badly_distorted_residue_still_matches_its_entry_and_fails_its_bonds() {
        // Every third heavy atom moved 0.9 A, breaking more than a tenth of
        // the bonds: the pose is NAG's, and what it fails is its geometry.
        let mut pose = ideal_pose(NAG, "NAG");
        for atom in pose.atoms.iter_mut().step_by(3) {
            atom.position.z -= 0.9;
        }
        let ligand = everything(&pose);
        let file = std::mem::take(&mut pose.components);

        let (matched, refused) = matching(&pose, &ligand, &file);
        assert!(refused.is_empty());
        let pose = Protein {
            components: matched,
            ..pose
        };
        assert!(geometry(&pose, &ligand).unwrap().bonds > 0);
    }

    #[test]
    fn a_cluster_bonded_only_to_metals_matches_by_its_atoms() {
        // An iron-sulfur pair: no covalent radius bonds iron, so its bonds
        // cannot be read from distances and are not asked for.
        let entry = Dictionary::parse(
            "data_FES\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
             _chem_comp_atom.type_symbol\nFES FE1 FE\nFES FE2 FE\nFES S1 S\nFES S2 S\n\
             loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n\
             _chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n\
             FES FE1 S1 SING\nFES FE1 S2 SING\nFES FE2 S1 SING\nFES FE2 S2 SING\n",
        );
        let mut pose = ideal_pose(GPP, "GPP");
        pose.atoms.truncate(4);
        for (atom, (name, element, at)) in pose.atoms.iter_mut().zip([
            ("FE1", "FE", Vec3::new(1.35, 0.0, 0.0)),
            ("FE2", "FE", Vec3::new(-1.35, 0.0, 0.0)),
            ("S1", "S", Vec3::new(0.0, 1.8, 0.0)),
            ("S2", "S", Vec3::new(0.0, -1.8, 0.0)),
        ]) {
            atom.name = FixedStr::new(name);
            atom.element = FixedStr::new(element);
            atom.residue_name = FixedStr::new("FES");
            atom.position = at;
        }
        let ligand = everything(&pose);

        let (matched, refused) = matching(&pose, &ligand, &entry);
        assert!(refused.is_empty(), "{refused:?}");
        assert!(matched.get("FES").is_some());
    }

    #[test]
    fn a_residue_named_as_a_component_but_built_otherwise_is_not_matched() {
        // GPP's atoms renamed as NAG's would be: no NAG entry matches.
        let mut pose = ideal_pose(GPP, "GPP");
        for atom in &mut pose.atoms {
            atom.residue_name = FixedStr::new("NAG");
        }
        let ligand = everything(&pose);
        let nag = Dictionary::parse(NAG);

        let (matched, refused) = matching(&pose, &ligand, &nag);
        assert!(matched.is_empty());
        assert_eq!(refused, [FixedStr::new("NAG")]);
        assert_eq!(
            stereo(&pose, &ligand, &nag).unchecked,
            [FixedStr::new("NAG")]
        );
    }
}
