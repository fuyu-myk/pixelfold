//! Pixelfold core: structure model, spatial index, selection, secondary structure,
//! SASA, interactions, and residue interaction network.

mod assembly;
pub mod components;
mod dssp;
pub mod fixed_str;
pub mod interactions;
pub mod parser;
pub mod rin;
pub mod sasa;
pub mod select;
mod spatial;
pub mod structure;

pub use components::{Component, Dictionary};
pub use dssp::HBond;
pub use fixed_str::FixedStr;
pub use interactions::{Interaction, InteractionKind};
pub use rin::Network;
pub use select::{AtomSet, SelectError, select};
pub use structure::{
    AltlocPolicy, Atom, BiologicalAssembly, ChainId, Molecule, Protein, SecondaryStructure,
    calculate_bfactor_range, get_calpha_connections, get_calpha_indices,
};
