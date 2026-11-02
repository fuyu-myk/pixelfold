//! Building the pane's network: one model's, or an ensemble's, laid out once
//! and its residues mapped to the atoms of the model the 3D view shows.

use std::collections::HashMap;

use anyhow::Result;
use glam::Vec3;
use pixelfold_core::interactions::InteractionKind;
use pixelfold_core::parser::{Load, ParsedStructure};
use pixelfold_core::rin::{self, Ensemble, Network};
use pixelfold_core::sasa::{SurfaceCalculator, relative_sasa};
use pixelfold_core::{Protein, interactions};

use super::palette::KIND_ORDER;
use super::{ChainFilter, EdgeVis, FREQUENCY_FLOORS, LayoutMode, NetworkView, NodeColorMode};

/// Fruchterman-Reingold steps run once when the pane is opened.
const LAYOUT_ITERATIONS: usize = 300;

impl NetworkView {
    /// Detect the interactions, build and lay out the network, and map residues
    /// to atoms both ways.
    pub fn build(protein: &Protein) -> Self {
        let detected = interactions::detect(protein);
        let network = rin::build(protein, &detected);
        Self::from_network(protein, network)
    }

    /// Detect the interactions of every one of `models` of `source`, built as
    /// `load` says, and lay the network of them all over `protein`, the model
    /// shown.
    pub fn build_ensemble(
        protein: &Protein,
        source: &ParsedStructure,
        load: &Load,
        models: &[usize],
    ) -> Result<Self> {
        let mut ensemble = Ensemble::default();
        for &model in models {
            let built = source.protein(&Load {
                model: Some(model),
                ..load.clone()
            })?;
            ensemble.add(&built, &interactions::detect(&built));
        }

        Ok(Self::from_ensemble(protein, &ensemble))
    }

    /// Lay `ensemble`'s network over `protein`, knowing how many models each
    /// residue pair interacts in by any kind, which its edges, one per kind,
    /// cannot say.
    pub(super) fn from_ensemble(protein: &Protein, ensemble: &Ensemble) -> Self {
        let mut view = Self::from_network(protein, ensemble.network(0.0));
        let index: HashMap<&str, usize> = view
            .network
            .nodes
            .iter()
            .enumerate()
            .map(|(i, node)| (node.id.as_str(), i))
            .collect();
        view.touches = ensemble
            .contacts()
            .map(|(source, target, models)| ((index[source], index[target]), models))
            .collect();

        view
    }

    /// Lay out an already-built network and map its residues to `protein`'s
    /// atoms both ways. Split from [`build`] so the mapping can be tested without
    /// running the interaction engine.
    pub(super) fn from_network(protein: &Protein, network: Network) -> Self {
        let positions: Vec<(f32, f32)> = rin::force_directed(&network, LAYOUT_ITERATIONS)
            .positions
            .iter()
            .map(|p| (p.x, p.y))
            .collect();

        let edges: Vec<EdgeVis> = {
            let index: HashMap<&str, usize> = network
                .nodes
                .iter()
                .enumerate()
                .map(|(i, node)| (node.id.as_str(), i))
                .collect();
            network
                .edges
                .iter()
                .filter_map(|edge| {
                    Some(EdgeVis {
                        u: *index.get(edge.source.as_str())?,
                        v: *index.get(edge.target.as_str())?,
                        kind: edge.kind,
                        energy: edge.energy,
                        interchain: edge.interchain,
                        frequency: network.frequency(edge),
                    })
                })
                .collect()
        };

        // One model's pairs interact in it or not; an ensemble's are counted by
        // `from_ensemble`.
        let touches = edges.iter().map(|e| ((e.u, e.v), 1)).collect();

        // The kinds actually present, in the legend's fixed order.
        let present: std::collections::BTreeSet<InteractionKind> =
            edges.iter().map(|e| e.kind).collect();
        let present_kinds: Vec<InteractionKind> = KIND_ORDER
            .into_iter()
            .filter(|k| present.contains(k))
            .collect();

        // Size hubs by betweenness centrality; the square root lifts the many
        // mid-ranked residues out of the shadow of the few very high ones.
        let betweenness = rin::analyze(&network).betweenness;
        let max = betweenness.iter().copied().fold(0.0f64, f64::max);
        let hub: Vec<f32> = betweenness
            .iter()
            .map(|&b| {
                if max > 0.0 {
                    (b / max).sqrt() as f32
                } else {
                    0.0
                }
            })
            .collect();

        let atoms = AtomMap::of(protein, &network);

        Self {
            network,
            mode: LayoutMode::default(),
            color_mode: NodeColorMode::default(),
            positions,
            edges,
            present_kinds,
            kind_filter: None,
            chain_filter: ChainFilter::default(),
            min_frequency: FREQUENCY_FLOORS[0],
            touches,
            hub,
            rsa: atoms.rsa,
            repr_atom: atoms.repr_atom,
            node_of_atom: atoms.node_of_atom,
            node_pos3: atoms.node_pos3,
            node_cell: Vec::new(),
        }
    }

    /// This view laid over `protein`, another model of the same file: the
    /// network, its layout and every setting kept, the residues mapped to the
    /// new model's atoms and placed and buried as it has them.
    pub fn remapped(self, protein: &Protein) -> Self {
        let atoms = AtomMap::of(protein, &self.network);

        Self {
            rsa: atoms.rsa,
            repr_atom: atoms.repr_atom,
            node_of_atom: atoms.node_of_atom,
            node_pos3: atoms.node_pos3,
            ..self
        }
    }
}

/// A network's residues mapped to the atoms of one model, and what that model
/// says of each: where it sits and how buried it is.
struct AtomMap {
    /// A representative atom for each node, `None` for a residue of an
    /// ensemble's network the model lacks.
    repr_atom: Vec<Option<usize>>,
    /// The node each atom of a residue in the network belongs to.
    node_of_atom: HashMap<usize, usize>,
    /// Each node's representative atom's position; a residue the model lacks
    /// sits at the middle of those it has.
    node_pos3: Vec<Vec3>,
    /// Each node's relative solvent accessibility, `None` for residues with no
    /// reference maximum (ligands, non-standard) or none in the model.
    rsa: Vec<Option<f32>>,
}

impl AtomMap {
    fn of(protein: &Protein, network: &Network) -> Self {
        let key_index: HashMap<(String, i32, Option<char>), usize> = network
            .nodes
            .iter()
            .enumerate()
            .map(|(i, node)| ((node.chain.clone(), node.resi, node.icode), i))
            .collect();

        let mut repr_atom = vec![None; network.nodes.len()];
        let mut node_of_atom = HashMap::new();
        for (atom_idx, atom) in protein.atoms.iter().enumerate() {
            let key = (
                atom.chain_id.as_str().to_string(),
                atom.residue_seq,
                atom.insertion_code,
            );
            if let Some(&node) = key_index.get(&key) {
                node_of_atom.insert(atom_idx, node);
                repr_atom[node].get_or_insert(atom_idx);
            }
        }

        let present: Vec<Vec3> = repr_atom
            .iter()
            .flatten()
            .map(|&atom| protein.atoms[atom].position)
            .collect();
        let middle = present.iter().sum::<Vec3>() / present.len().max(1) as f32;
        let node_pos3 = repr_atom
            .iter()
            .map(|atom| atom.map_or(middle, |atom| protein.atoms[atom].position))
            .collect();

        let atom_sasa = SurfaceCalculator::default().calculate_atom_sasa(&protein.atoms);
        let mut node_sasa = vec![0.0f32; network.nodes.len()];
        for (&atom_idx, &node) in &node_of_atom {
            node_sasa[node] += atom_sasa[atom_idx];
        }
        let rsa: Vec<Option<f32>> = network
            .nodes
            .iter()
            .enumerate()
            .map(|(i, node)| repr_atom[i].and_then(|_| relative_sasa(node_sasa[i], &node.resn)))
            .collect();

        Self {
            repr_atom,
            node_of_atom,
            node_pos3,
            rsa,
        }
    }
}
