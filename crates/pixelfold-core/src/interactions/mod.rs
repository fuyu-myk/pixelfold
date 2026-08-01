//! Non-covalent interaction detection.
//!
//! One engine over the uniform spatial grid, pointed at different atom subsets:
//! protein-protein, protein-ligand, and docking all fall out of the same
//! detectors rather than being separate features.
//!
//! Geometric thresholds live in [`params`] and are taken from PLIP's source;
//! atom roles live in [`classify`] and follow explicit residue tables.

pub mod aromatic;
pub mod classify;
pub mod connectivity;
mod detectors;
mod energy;
pub mod hydrogens;
pub mod params;
pub mod topology;

use crate::structure::Protein;

pub use classify::Charge;

/// The kind of a detected interaction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum InteractionKind {
    HydrogenBond,
    SaltBridge,
    Hydrophobic,
    PiStacking,
    PiCation,
    HalogenBond,
    WaterBridge,
    MetalCoordination,
    Disulfide,
}

impl InteractionKind {
    /// A short, stable label for reports and exports.
    pub fn label(&self) -> &'static str {
        match self {
            InteractionKind::HydrogenBond => "hydrogen-bond",
            InteractionKind::SaltBridge => "salt-bridge",
            InteractionKind::Hydrophobic => "hydrophobic",
            InteractionKind::PiStacking => "pi-stacking",
            InteractionKind::PiCation => "pi-cation",
            InteractionKind::HalogenBond => "halogen-bond",
            InteractionKind::WaterBridge => "water-bridge",
            InteractionKind::MetalCoordination => "metal-coordination",
            InteractionKind::Disulfide => "disulfide",
        }
    }
}

/// One detected interaction between two atom groups.
///
/// `atoms_a` and `atoms_b` hold the participating atom indices on each side: a
/// single atom for pairwise types (disulfide, metal coordination), or the whole
/// group for centroid types (a salt bridge's charged group).
#[derive(Clone, Debug)]
pub struct Interaction {
    pub kind: InteractionKind,
    pub atoms_a: Vec<usize>,
    pub atoms_b: Vec<usize>,
    /// The value the type's cutoff was applied to, e.g. centroid separation.
    pub distance: f32,
    pub angle: Option<f32>,
    /// The atom mediating an indirect interaction, like the bridging water
    /// oxygen of a water bridge. `None` where the two sides meet directly.
    pub bridge: Option<usize>,
}

/// Detect every supported interaction in a structure.
pub fn detect(protein: &Protein) -> Vec<Interaction> {
    let bonds = connectivity::bonds(protein);

    let mut interactions = Vec::new();
    interactions.extend(detectors::hydrogen_bonds(protein, &bonds));
    interactions.extend(detectors::disulfides(protein));
    interactions.extend(detectors::metal_coordination(protein));
    interactions.extend(detectors::salt_bridges(protein, &bonds));
    interactions.extend(detectors::pi_stacking(protein));
    interactions.extend(detectors::pi_cation(protein, &bonds));
    interactions.extend(detectors::hydrophobic_contacts(protein, &bonds));
    interactions.extend(detectors::water_bridges(protein, &bonds));
    interactions.extend(detectors::halogen_bonds(protein, &bonds));

    dedup_conformers(protein, interactions)
}

/// An atom's residue-and-name identity, ignoring its alternate-location label.
type AtomIdentity = (
    crate::fixed_str::FixedStr<4>,
    u32,
    Option<char>,
    crate::fixed_str::FixedStr<4>,
);

/// Collapse interactions that are identical once alternate-location labels are
/// ignored.
///
/// Under `--altloc all` every conformer of an atom is a separate index, so one
/// hydrogen bond from a split hydroxyl is detected once per conformer; keyed on
/// residue-and-name identity (with each side's role preserved, so a reciprocal
/// bond is not folded into its partner) those duplicates fold to a single edge,
/// keeping the first modelled. With one conformer kept per residue (the default)
/// every identity is already unique, so this changes nothing.
fn dedup_conformers(protein: &Protein, interactions: Vec<Interaction>) -> Vec<Interaction> {
    let identity = |idx: usize| -> AtomIdentity {
        let atom = &protein.atoms[idx];
        (
            atom.chain_id,
            atom.residue_seq,
            atom.insertion_code,
            atom.name,
        )
    };
    let group = |atoms: &[usize]| {
        let mut ids: Vec<AtomIdentity> = atoms.iter().map(|&i| identity(i)).collect();
        ids.sort_unstable();
        ids
    };

    let mut seen = std::collections::HashSet::new();
    interactions
        .into_iter()
        .filter(|it| {
            seen.insert((
                it.kind,
                group(&it.atoms_a),
                group(&it.atoms_b),
                it.bridge.map(identity),
            ))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixed_str::FixedStr;
    use crate::structure::{Atom, SecondaryStructure};
    use glam::Vec3;

    fn atom(seq: u32, name: &str, altloc: Option<char>) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(&name[..1]),
            residue_name: FixedStr::new("SER"),
            residue_seq: seq,
            chain_id: FixedStr::new("A"),
            is_hetatm: false,
            altloc,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: Vec3::ZERO,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    fn protein(atoms: Vec<Atom>) -> Protein {
        Protein {
            atoms,
            title: String::new(),
            surface_points: Vec::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        }
    }

    fn hbond(atoms_a: Vec<usize>, atoms_b: Vec<usize>) -> Interaction {
        Interaction {
            kind: InteractionKind::HydrogenBond,
            atoms_a,
            atoms_b,
            distance: 3.0,
            angle: None,
            bridge: None,
        }
    }

    #[test]
    fn conformer_duplicate_edges_collapse_but_reciprocals_and_distinct_survive() {
        // A split hydroxyl: OG conformer A (0) and B (1) donate to the same
        // acceptor O (2). The two are one interaction reported per conformer.
        let structure = protein(vec![
            atom(5, "OG", Some('A')),
            atom(5, "OG", Some('B')),
            atom(10, "O", None),
        ]);
        let deduped = dedup_conformers(
            &structure,
            vec![
                hbond(vec![0], vec![2]), // OG(A) -> O
                hbond(vec![1], vec![2]), // OG(B) -> O: conformer duplicate, dropped
                hbond(vec![2], vec![0]), // reciprocal (roles swapped): distinct, kept
            ],
        );

        assert_eq!(deduped.len(), 2);
        assert_eq!(
            deduped[0].atoms_a,
            vec![0],
            "the first-modelled conformer is kept"
        );
        assert_eq!(
            deduped[1].atoms_a,
            vec![2],
            "the reciprocal edge is not folded away"
        );
    }

    #[test]
    fn single_conformer_edges_are_left_untouched() {
        // Distinct atom identities (the default, one conformer per residue): the
        // dedup must remove nothing.
        let structure = protein(vec![
            atom(5, "OG", None),
            atom(10, "O", None),
            atom(11, "N", None),
        ]);
        let deduped = dedup_conformers(
            &structure,
            vec![hbond(vec![0], vec![1]), hbond(vec![0], vec![2])],
        );
        assert_eq!(deduped.len(), 2);
    }
}
