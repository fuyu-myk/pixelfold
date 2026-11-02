//! The linked residue-interaction-network pane.
//!
//! A force-directed 2D view of the residue interaction network, drawn on a
//! braille canvas beside the 3D structure. Clicking a node selects that residue,
//! which the structure pane highlights; picking a residue in the structure
//! highlights its node here. One `selected_atom_idx` drives both panes, so the
//! two views stay in sync.
//!
//! For a multi-model file the pane can show the network of the whole ensemble
//! instead: every edge any model has, brightened by the fraction of models it
//! is present in, laid over the atoms of whichever model the 3D view shows.

mod build;
mod draw;
mod palette;
#[cfg(test)]
mod tests;

use std::collections::HashMap;

use glam::Vec3;

use pixelfold_core::interactions::InteractionKind;
use pixelfold_core::rin::Network;

/// How the network's nodes are placed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LayoutMode {
    /// Force-directed: positions from the graph's topology, independent of the
    /// structure's orientation.
    #[default]
    Force,
    /// Spatial: each residue at its 3D position projected through the camera, so
    /// the graph rotates with the structure.
    Spatial,
}

/// How the network's nodes are coloured.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NodeColorMode {
    /// By DSSP secondary structure, the same palette the 3D view uses.
    #[default]
    Structure,
    /// By burial, or relative solvent accessibility.
    Burial,
}

/// The fractions of an ensemble's models an edge must be present in to draw,
/// stepped through in turn: every edge, then those in half the models, then
/// those in nearly all.
const FREQUENCY_FLOORS: [f32; 3] = [0.0, 0.5, 0.9];

/// What the pane shows: the network of the model the 3D view shows, or of
/// every model folded into one with the edges in fewer than a floor of them
/// hidden. Kept by the app rather than the view, so a key pressed while a view
/// is still being built never loses it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NetworkScope {
    #[default]
    Model,
    Ensemble {
        /// The index of the floor in [`FREQUENCY_FLOORS`].
        floor: usize,
    },
}

impl NetworkScope {
    /// The next scope: the ensemble after the model, each floor in turn, and
    /// the model again after the last.
    pub fn next(self) -> Self {
        match self {
            NetworkScope::Model => NetworkScope::Ensemble { floor: 0 },
            NetworkScope::Ensemble { floor } if floor + 1 < FREQUENCY_FLOORS.len() => {
                NetworkScope::Ensemble { floor: floor + 1 }
            }
            NetworkScope::Ensemble { .. } => NetworkScope::Model,
        }
    }

    pub fn is_ensemble(self) -> bool {
        matches!(self, NetworkScope::Ensemble { .. })
    }

    /// The least fraction of the models an edge drawn is present in.
    pub fn min_frequency(self) -> f32 {
        match self {
            NetworkScope::Model => 0.0,
            NetworkScope::Ensemble { floor } => FREQUENCY_FLOORS[floor],
        }
    }
}

/// One edge as the pane draws it: its two node indices, kind (which sets the
/// hue), summed energy (which sets the brightness of one model's network),
/// whether it crosses chains (for the chain filter), and the fraction of the
/// models it is present in (the brightness of an ensemble's network).
struct EdgeVis {
    u: usize,
    v: usize,
    kind: InteractionKind,
    energy: f32,
    interchain: bool,
    frequency: f32,
}

/// Which edges the chain filter keeps.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum ChainFilter {
    /// Every edge.
    #[default]
    All,
    /// Only edges within a single chain.
    IntraOnly,
    /// Only edges between different chains, the inter-chain interface.
    InterOnly,
}

pub struct NetworkView {
    network: Network,
    mode: LayoutMode,
    color_mode: NodeColorMode,
    /// Force-directed layout position of each node in the unit square, y up.
    positions: Vec<(f32, f32)>,
    /// The 3D position of each node's representative atom, projected through the
    /// camera in spatial mode.
    node_pos3: Vec<Vec3>,
    /// Each edge for drawing: node indices, kind (hue), and energy (brightness).
    edges: Vec<EdgeVis>,
    /// The interaction kinds actually present, in display order, for cycling the
    /// kind filter over only the kinds this network has.
    present_kinds: Vec<InteractionKind>,
    /// Show only this interaction kind's edges, or every kind when `None`.
    kind_filter: Option<InteractionKind>,
    /// Which edges the chain filter keeps.
    chain_filter: ChainFilter,
    /// The least fraction of an ensemble's models an edge drawn is present in.
    min_frequency: f32,
    /// How many of the models each pair of nodes (the lower index first)
    /// interacts in, by any kind.
    touches: HashMap<(usize, usize), usize>,
    /// Normalised betweenness centrality per node (`sqrt` scaled to `0..1`), for
    /// hub sizing.
    hub: Vec<f32>,
    /// Relative solvent accessibility per node, `None` for residues with no
    /// reference maximum (ligands, non-standard). Drives the burial colour mode.
    rsa: Vec<Option<f32>>,
    /// A representative atom for each node, for selecting into the 3D view;
    /// `None` for a residue of an ensemble's network the shown model lacks.
    repr_atom: Vec<Option<usize>>,
    /// The node each interacting atom belongs to, for reading the selection back.
    node_of_atom: HashMap<usize, usize>,
    /// Where each node last drew, in terminal cells, for click picking.
    node_cell: Vec<(u16, u16)>,
}

impl NetworkView {
    pub fn network(&self) -> &Network {
        &self.network
    }

    /// Whether the network is an ensemble's, drawn by how often each edge is
    /// present rather than by its energy.
    fn is_ensemble(&self) -> bool {
        self.network.models > 1
    }

    /// Hide the edges in fewer than `fraction` of the models.
    pub fn set_min_frequency(&mut self, fraction: f32) {
        self.min_frequency = fraction;
    }

    #[cfg(test)]
    pub fn min_frequency(&self) -> f32 {
        self.min_frequency
    }

    /// The residues `node` interacts with and how many of the models it does
    /// in, by any kind of interaction, the most often first.
    pub fn partners(&self, node: usize) -> Vec<(&str, usize)> {
        let mut partners: Vec<(usize, usize)> = self
            .touches
            .iter()
            .filter_map(|(&(u, v), &models)| match (u == node, v == node) {
                (true, _) => Some((v, models)),
                (_, true) => Some((u, models)),
                _ => None,
            })
            .collect();
        partners.sort_by_key(|&(partner, models)| (std::cmp::Reverse(models), partner));

        partners
            .into_iter()
            .map(|(partner, models)| (self.network.nodes[partner].id.as_str(), models))
            .collect()
    }

    /// Switch between the force-directed and spatial layouts.
    pub fn toggle_mode(&mut self) {
        self.mode = match self.mode {
            LayoutMode::Force => LayoutMode::Spatial,
            LayoutMode::Spatial => LayoutMode::Force,
        };
    }

    /// Switch node colouring between secondary structure and burial.
    pub fn toggle_color(&mut self) {
        self.color_mode = match self.color_mode {
            NodeColorMode::Structure => NodeColorMode::Burial,
            NodeColorMode::Burial => NodeColorMode::Structure,
        };
    }

    /// Step the kind filter: every kind, then each present kind in turn, and
    /// back to every kind.
    pub fn cycle_kind_filter(&mut self) {
        self.kind_filter = match self.kind_filter {
            None => self.present_kinds.first().copied(),
            Some(current) => match self.present_kinds.iter().position(|&k| k == current) {
                Some(i) => self.present_kinds.get(i + 1).copied(),
                None => None,
            },
        };
    }

    /// Step the chain filter: all edges, then intra-chain only, then the
    /// inter-chain interface only, and back.
    pub fn cycle_chain_filter(&mut self) {
        self.chain_filter = match self.chain_filter {
            ChainFilter::All => ChainFilter::IntraOnly,
            ChainFilter::IntraOnly => ChainFilter::InterOnly,
            ChainFilter::InterOnly => ChainFilter::All,
        };
    }

    /// A representative atom of a node, to select into the 3D view.
    pub fn atom_of_node(&self, node: usize) -> Option<usize> {
        self.repr_atom[node]
    }

    /// The node an atom belongs to, if the atom is in the network.
    pub fn node_of_atom(&self, atom: usize) -> Option<usize> {
        self.node_of_atom.get(&atom).copied()
    }
}

/// Whether an edge passes the kind, chain and frequency filters and should be
/// drawn.
fn edge_passes(
    edge: &EdgeVis,
    kind_filter: Option<InteractionKind>,
    chain_filter: ChainFilter,
    min_frequency: f32,
) -> bool {
    let kind_ok = kind_filter.is_none_or(|k| k == edge.kind);
    let chain_ok = match chain_filter {
        ChainFilter::All => true,
        ChainFilter::IntraOnly => !edge.interchain,
        ChainFilter::InterOnly => edge.interchain,
    };

    kind_ok && chain_ok && edge.frequency >= min_frequency
}
