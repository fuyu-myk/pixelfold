//! Residue interaction networks.
//!
//! Two graphs live here. [`hbond`] is the backbone hydrogen-bond graph the DSSP
//! pass feeds, used by the viewer's network overlay. [`network`] is the residue
//! interaction network proper: one node per residue, one edge per typed
//! non-covalent interaction, built over the whole [`crate::interactions`] engine.

mod analysis;
mod hbond;
mod layout;
mod network;

pub use analysis::{RinAnalysis, analyze};
pub use hbond::{BondType, HBondEdge, HBondGraph, HBondMotifs, NetworkAnalysis, ResidueNode};
pub use layout::{Layout, force_directed};
pub use network::{Edge, Network, Node, build};
