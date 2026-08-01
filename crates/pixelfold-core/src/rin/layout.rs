//! A force-directed layout of the residue interaction network.
//!
//! Fruchterman-Reingold: every node repels every other, edges pull their two
//! endpoints together, and the whole thing cools over a fixed number of steps.
//! Communities that interact densely settle close, so the layout reads as the
//! structure's contact map rather than its sequence. It is deterministic (a
//! fixed circular start, no randomness), so the same network always lays out the
//! same way and it can be tested.

use std::collections::{BTreeSet, HashMap};

use glam::Vec2;

use crate::rin::network::Network;

/// Node positions in the unit square `[0, 1] x [0, 1]`, aligned with
/// `network.nodes` by index. The consumer maps them to its own viewport.
#[derive(Clone, Debug, PartialEq)]
pub struct Layout {
    pub positions: Vec<Vec2>,
}

/// The distance below which two nodes are treated as coincident, to keep the
/// repulsion force finite.
const MIN_DISTANCE: f32 = 1e-4;

/// Lay `network` out over `iterations` steps of Fruchterman-Reingold.
pub fn force_directed(network: &Network, iterations: usize) -> Layout {
    let n = network.nodes.len();
    if n == 0 {
        return Layout {
            positions: Vec::new(),
        };
    }
    if n == 1 {
        return Layout {
            positions: vec![Vec2::splat(0.5)],
        };
    }

    let edges = edge_indices(network);
    let mut pos = circle_start(n);

    // Ideal edge length for a unit-area frame, and a cooling cap on how far a
    // node may move in one step.
    let k = (1.0 / n as f32).sqrt();
    let mut temperature = 0.1;

    for _ in 0..iterations {
        let mut disp = vec![Vec2::ZERO; n];

        for i in 0..n {
            for j in (i + 1)..n {
                let delta = pos[i] - pos[j];
                let dist = delta.length().max(MIN_DISTANCE);
                let repulsion = (k * k / dist) * (delta / dist);
                disp[i] += repulsion;
                disp[j] -= repulsion;
            }
        }

        for &(u, v) in &edges {
            let delta = pos[u] - pos[v];
            let dist = delta.length().max(MIN_DISTANCE);
            let attraction = (dist * dist / k) * (delta / dist);
            disp[u] -= attraction;
            disp[v] += attraction;
        }

        for (p, d) in pos.iter_mut().zip(&disp) {
            let dist = d.length().max(MIN_DISTANCE);
            *p += (*d / dist) * dist.min(temperature);
        }

        temperature *= 0.95;
    }

    normalize(&mut pos);
    Layout { positions: pos }
}

/// Distinct undirected edges as node-index pairs. Several interaction kinds
/// between one residue pair collapse to one spring.
fn edge_indices(network: &Network) -> Vec<(usize, usize)> {
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
            let u = *index.get(edge.source.as_str())?;
            let v = *index.get(edge.target.as_str())?;
            (u != v).then_some((u.min(v), u.max(v)))
        })
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

/// A deterministic starting layout: the nodes evenly spaced on a circle. Any
/// non-degenerate start breaks the symmetry the forces would otherwise preserve.
fn circle_start(n: usize) -> Vec<Vec2> {
    (0..n)
        .map(|i| {
            let angle = (i as f32) / (n as f32) * std::f32::consts::TAU;
            Vec2::new(0.5 + 0.4 * angle.cos(), 0.5 + 0.4 * angle.sin())
        })
        .collect()
}

/// Rescale the positions to fill the unit square, preserving aspect so the
/// layout is not stretched.
fn normalize(pos: &mut [Vec2]) {
    let mut min = Vec2::splat(f32::INFINITY);
    let mut max = Vec2::splat(f32::NEG_INFINITY);
    for &p in pos.iter() {
        min = min.min(p);
        max = max.max(p);
    }

    let extent = (max - min).max_element().max(MIN_DISTANCE);
    let center = (min + max) * 0.5;
    for p in pos.iter_mut() {
        // Center on 0.5 and scale the larger axis to span the frame.
        *p = Vec2::splat(0.5) + (*p - center) / extent;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixed_str::FixedStr;
    use crate::interactions::{Interaction, InteractionKind};
    use crate::rin::network::build;
    use crate::structure::{Atom, Protein, SecondaryStructure};
    use glam::Vec3;

    fn atom(chain: &str, seq: u32, name: &str) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(&name[..1]),
            residue_name: FixedStr::new("ALA"),
            residue_seq: seq,
            chain_id: FixedStr::new(chain),
            is_hetatm: false,
            altloc: None,
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
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        }
    }

    /// A chain of residues each hydrogen-bonded to the next.
    fn chain_network(len: u32) -> Network {
        let atoms: Vec<Atom> = (0..len).map(|i| atom("A", i + 1, "N")).collect();
        let interactions: Vec<Interaction> = (0..len - 1)
            .map(|i| Interaction {
                kind: InteractionKind::HydrogenBond,
                atoms_a: vec![i as usize],
                atoms_b: vec![(i + 1) as usize],
                distance: 3.0,
                angle: None,
                bridge: None,
            })
            .collect();
        build(&protein(atoms), &interactions)
    }

    #[test]
    fn an_empty_network_has_no_positions() {
        let layout = force_directed(&Network::default(), 50);
        assert!(layout.positions.is_empty());
    }

    #[test]
    fn positions_are_finite_normalised_and_one_per_node() {
        let network = chain_network(6);
        let layout = force_directed(&network, 200);

        assert_eq!(layout.positions.len(), network.nodes.len());
        for p in &layout.positions {
            assert!(p.x.is_finite() && p.y.is_finite(), "non-finite position");
            assert!(
                (-0.01..=1.01).contains(&p.x) && (-0.01..=1.01).contains(&p.y),
                "position {p:?} outside the unit square"
            );
        }
    }

    #[test]
    fn the_layout_is_deterministic() {
        let network = chain_network(8);
        assert_eq!(force_directed(&network, 150), force_directed(&network, 150));
    }

    #[test]
    fn bonded_neighbours_settle_closer_than_the_chain_ends() {
        // Down a hydrogen-bonded chain, adjacent residues should end up nearer
        // each other than the two ends, which the springs pull apart the most.
        let network = chain_network(6);
        let p = force_directed(&network, 400).positions;

        let id_pos = |id: &str| {
            let i = network.nodes.iter().position(|n| n.id == id).unwrap();
            p[i]
        };
        let neighbour = id_pos("A/1").distance(id_pos("A/2"));
        let ends = id_pos("A/1").distance(id_pos("A/6"));
        assert!(
            neighbour < ends,
            "adjacent {neighbour} should be closer than the ends {ends}"
        );
    }
}
