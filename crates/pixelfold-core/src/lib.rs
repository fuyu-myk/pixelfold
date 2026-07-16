//! Pixelfold core: structure model, spatial index, selection, secondary structure,
//! SASA, interactions, and residue interaction network.

pub mod assembly;
pub mod dssp;
pub mod fixed_str;
pub mod interactions;
pub mod parser;
pub mod rin;
pub mod sasa;
pub mod spatial;
pub mod structure;

pub use dssp::HBond;
pub use fixed_str::FixedStr;
pub use interactions::{Interaction, InteractionKind};
pub use sasa::SurfacePoint;
pub use structure::{
    AltlocPolicy, Atom, BiologicalAssembly, DisplayMode, Protein, SecondaryStructure,
    calculate_bfactor_range, get_calpha_connections, get_calpha_indices,
};
