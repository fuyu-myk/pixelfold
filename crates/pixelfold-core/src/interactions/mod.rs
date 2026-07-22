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
    interactions.extend(detectors::salt_bridges(protein));
    interactions.extend(detectors::pi_stacking(protein));
    interactions.extend(detectors::pi_cation(protein));
    interactions.extend(detectors::hydrophobic_contacts(protein, &bonds));
    interactions.extend(detectors::water_bridges(protein, &bonds));
    interactions.extend(detectors::halogen_bonds(protein, &bonds));

    interactions
}
