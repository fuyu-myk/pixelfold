//! The residue interaction network.
//!
//! Nodes are residues; edges are the non-covalent interactions between them,
//! one per residue pair per interaction kind. The atom-level interactions the
//! engine reports are aggregated to the residue level, which is the scale a
//! network is read at: a lysine and a glutamate that make three hydrogen bonds
//! and a salt bridge are two nodes joined by two edges, not eight.
//!
//! The graph is undirected. A hydrogen bond has a donor and an acceptor, but a
//! residue pair joined by several interaction types has no single direction.
//! The atom-level direction is kept by the interaction list itself.

use std::collections::BTreeMap;

use crate::fixed_str::FixedStr;
use crate::interactions::{Interaction, InteractionKind};
use crate::structure::{Atom, Protein};

type ResidueKey = (FixedStr<4>, u32, Option<char>);

fn residue_key(atom: &Atom) -> ResidueKey {
    (atom.chain_id, atom.residue_seq, atom.insertion_code)
}

/// A residue, identified as `chain/seq` with the insertion code appended when it
/// has one (`A/60B`).
fn node_id(atom: &Atom) -> String {
    match atom.insertion_code {
        Some(code) => format!("{}/{}{}", atom.chain_id.as_str(), atom.residue_seq, code),
        None => format!("{}/{}", atom.chain_id.as_str(), atom.residue_seq),
    }
}

/// A residue in the network.
#[derive(Clone, Debug, PartialEq)]
pub struct Node {
    pub id: String,
    pub chain: String,
    pub resi: u32,
    pub icode: Option<char>,
    pub resn: String,
    /// The DSSP-style secondary-structure letter.
    pub ss: char,
    /// How many other residues this one interacts with.
    pub degree: usize,
}

/// One residue pair joined by interactions of a single kind.
#[derive(Clone, Debug, PartialEq)]
pub struct Edge {
    pub source: String,
    pub target: String,
    pub kind: InteractionKind,
    /// How many atom-level interactions of this kind join the pair.
    pub count: usize,
    /// The shortest of those interactions' gate distances.
    pub min_distance: f32,
    /// Whether the two residues are on different chains.
    pub interchain: bool,
    /// The summed interaction energy of this edge (kJ/mol): the per-contact weight
    /// of `kind` times `count`. See [`InteractionKind::contact_energy_kj`].
    pub energy: f32,
}

/// A residue interaction network: the residues that interact and the typed
/// edges between them.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Network {
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
}

/// What an edge accumulates before it is finished: the count and the running
/// minimum distance of the interactions folded into it.
struct EdgeStats {
    count: usize,
    min_distance: f32,
    interchain: bool,
}

/// Build the network of `interactions` over `protein`.
///
/// Every interaction names a residue on each side through its first atom. A
/// self-interaction is dropped. The residue pair is ordered so the two
/// directions of one interaction collapse to one edge.
pub fn build(protein: &Protein, interactions: &[Interaction]) -> Network {
    // Keyed by (lower residue, higher residue, kind).
    let mut edges: BTreeMap<(ResidueKey, ResidueKey, InteractionKind), EdgeStats> = BTreeMap::new();
    // The atom that stands for each residue, for its attributes.
    let mut representative: BTreeMap<ResidueKey, usize> = BTreeMap::new();

    for interaction in interactions {
        let (Some(&a), Some(&b)) = (interaction.atoms_a.first(), interaction.atoms_b.first())
        else {
            continue;
        };
        let (key_a, key_b) = (
            residue_key(&protein.atoms[a]),
            residue_key(&protein.atoms[b]),
        );
        if key_a == key_b {
            continue;
        }

        representative.entry(key_a).or_insert(a);
        representative.entry(key_b).or_insert(b);

        let interchain = protein.atoms[a].chain_id != protein.atoms[b].chain_id;
        let (low, high) = if key_a <= key_b {
            (key_a, key_b)
        } else {
            (key_b, key_a)
        };

        edges
            .entry((low, high, interaction.kind))
            .and_modify(|stats| {
                stats.count += 1;
                stats.min_distance = stats.min_distance.min(interaction.distance);
            })
            .or_insert(EdgeStats {
                count: 1,
                min_distance: interaction.distance,
                interchain,
            });
    }

    let id_of = |key: &ResidueKey| node_id(&protein.atoms[representative[key]]);
    let edges: Vec<Edge> = edges
        .into_iter()
        .map(|((low, high, kind), stats)| Edge {
            source: id_of(&low),
            target: id_of(&high),
            kind,
            count: stats.count,
            min_distance: stats.min_distance,
            interchain: stats.interchain,
            energy: kind.contact_energy_kj() * stats.count as f32,
        })
        .collect();

    Network {
        nodes: nodes(protein, &representative, &edges),
        edges,
    }
}

/// One node per residue that interacts, with its degree, the number of distinct
/// residues it touches.
fn nodes(
    protein: &Protein,
    representative: &BTreeMap<ResidueKey, usize>,
    edges: &[Edge],
) -> Vec<Node> {
    let mut partners: BTreeMap<&str, std::collections::BTreeSet<&str>> = BTreeMap::new();
    for edge in edges {
        partners
            .entry(&edge.source)
            .or_default()
            .insert(&edge.target);
        partners
            .entry(&edge.target)
            .or_default()
            .insert(&edge.source);
    }

    representative
        .values()
        .map(|&index| {
            let atom = &protein.atoms[index];
            let id = node_id(atom);
            let degree = partners.get(id.as_str()).map_or(0, |set| set.len());

            Node {
                chain: atom.chain_id.as_str().to_string(),
                resi: atom.residue_seq,
                icode: atom.insertion_code,
                resn: atom.residue_name.as_str().to_string(),
                ss: atom.secondary_structure.code(),
                degree,
                id,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::SecondaryStructure;
    use glam::Vec3;

    fn atom(chain: &str, seq: u32, name: &str, residue: &str) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(&name[..1]),
            residue_name: FixedStr::new(residue),
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

    fn interaction(kind: InteractionKind, a: usize, b: usize, distance: f32) -> Interaction {
        Interaction {
            kind,
            atoms_a: vec![a],
            atoms_b: vec![b],
            distance,
            angle: None,
            bridge: None,
        }
    }

    /// Three residues: A/1 and A/2 joined by two hydrogen bonds and a salt
    /// bridge, A/1 and B/5 by one hydrophobic contact.
    fn sample() -> (Protein, Vec<Interaction>) {
        let structure = protein(vec![
            atom("A", 1, "N", "LYS"),   // 0
            atom("A", 1, "NZ", "LYS"),  // 1
            atom("A", 2, "O", "GLU"),   // 2
            atom("A", 2, "OE1", "GLU"), // 3
            atom("B", 5, "CB", "LEU"),  // 4
        ]);
        let interactions = vec![
            interaction(InteractionKind::HydrogenBond, 0, 2, 3.0),
            interaction(InteractionKind::HydrogenBond, 1, 2, 2.8),
            interaction(InteractionKind::SaltBridge, 1, 3, 3.5),
            interaction(InteractionKind::Hydrophobic, 0, 4, 3.9),
        ];

        (structure, interactions)
    }

    #[test]
    fn interactions_of_one_kind_between_a_pair_collapse_to_one_edge() {
        let (structure, interactions) = sample();
        let network = build(&structure, &interactions);

        let hbond = network
            .edges
            .iter()
            .find(|e| e.kind == InteractionKind::HydrogenBond && !e.interchain)
            .expect("the intra-chain hydrogen bond edge");
        assert_eq!(hbond.count, 2, "two hydrogen bonds fold into one edge");
        assert!(
            (hbond.min_distance - 2.8).abs() < 1e-6,
            "the shorter is kept"
        );
    }

    #[test]
    fn a_pair_with_two_kinds_keeps_an_edge_for_each() {
        let (structure, interactions) = sample();
        let network = build(&structure, &interactions);

        let between: Vec<InteractionKind> = network
            .edges
            .iter()
            .filter(|e| {
                (e.source == "A/1" && e.target == "A/2") || (e.source == "A/2" && e.target == "A/1")
            })
            .map(|e| e.kind)
            .collect();

        assert!(between.contains(&InteractionKind::HydrogenBond));
        assert!(between.contains(&InteractionKind::SaltBridge));
        assert_eq!(between.len(), 2);
    }

    #[test]
    fn edge_energy_is_the_per_contact_weight_summed_over_the_count() {
        let (structure, interactions) = sample();
        let network = build(&structure, &interactions);

        // Two hydrogen bonds fold into one edge, so its energy is twice one bond.
        let hbond = network
            .edges
            .iter()
            .find(|e| e.kind == InteractionKind::HydrogenBond && !e.interchain)
            .expect("the intra-chain hydrogen bond edge");
        assert_eq!(hbond.count, 2);
        assert!(
            (hbond.energy - 2.0 * InteractionKind::HydrogenBond.contact_energy_kj()).abs() < 1e-4
        );

        // The single hydrophobic contact carries one per-contact weight.
        let hydrophobic = network
            .edges
            .iter()
            .find(|e| e.kind == InteractionKind::Hydrophobic)
            .expect("the hydrophobic edge");
        assert_eq!(hydrophobic.count, 1);
        assert!(
            (hydrophobic.energy - InteractionKind::Hydrophobic.contact_energy_kj()).abs() < 1e-4
        );
    }

    #[test]
    fn degree_counts_distinct_partner_residues() {
        let (structure, interactions) = sample();
        let network = build(&structure, &interactions);

        // A/1 touches A/2 (hydrogen bond and salt bridge) and B/5 (hydrophobic):
        // two partners, not three edges.
        let a1 = network.nodes.iter().find(|n| n.id == "A/1").expect("A/1");
        assert_eq!(a1.degree, 2);

        let b5 = network.nodes.iter().find(|n| n.id == "B/5").expect("B/5");
        assert_eq!(b5.degree, 1);
    }

    #[test]
    fn only_interacting_residues_are_nodes() {
        let mut atoms = sample().0.atoms;
        atoms.push(atom("A", 3, "CA", "GLY")); // interacts with nothing
        let structure = protein(atoms);

        let network = build(&structure, &sample().1);
        assert!(
            network.nodes.iter().all(|n| n.id != "A/3"),
            "an isolated residue is not in the network"
        );
        assert_eq!(network.nodes.len(), 3);
    }

    #[test]
    fn an_edge_across_chains_is_marked_interchain() {
        let (structure, interactions) = sample();
        let network = build(&structure, &interactions);

        let hydrophobic = network
            .edges
            .iter()
            .find(|e| e.kind == InteractionKind::Hydrophobic)
            .expect("the hydrophobic edge");
        assert!(hydrophobic.interchain, "A and B are different chains");
    }

    #[test]
    fn the_network_is_deterministic() {
        let (structure, interactions) = sample();

        assert_eq!(
            build(&structure, &interactions),
            build(&structure, &interactions)
        );
    }
}
