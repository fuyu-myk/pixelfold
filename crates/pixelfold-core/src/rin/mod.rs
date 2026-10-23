//! Residue interaction networks: one node per residue, one edge per typed
//! non-covalent interaction, built over the whole [`crate::interactions`]
//! engine.

mod analysis;
mod layout;
mod network;

pub use analysis::{RinAnalysis, analyze, shortest_path};
pub use layout::{Layout, force_directed};
pub use network::{Edge, Network, Node, build};
