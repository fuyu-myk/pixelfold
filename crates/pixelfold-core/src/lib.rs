//! Pixelfold core: structure model, spatial index, selection, secondary structure,
//! SASA, interactions, and residue interaction network.

pub mod parser;
pub mod rin;
pub mod sasa;
pub mod structure;

pub use sasa::SurfacePoint;
pub use structure::{
    AltlocPolicy, Atom, DisplayMode, Protein, SecondaryStructure, calculate_bfactor_range,
    get_calpha_connections, get_calpha_indices,
};
