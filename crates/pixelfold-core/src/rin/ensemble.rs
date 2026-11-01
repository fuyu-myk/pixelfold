//! The residue interaction network of an ensemble: an NMR bundle, frames of a
//! trajectory, samples of a co-folding model.
//!
//! The models' networks fold into one, each edge weighted by the fraction of
//! models it is present in, the probabilistic network RING 3.0 builds
//! (Clementel et al., Nucleic Acids Res. 50:W651, 2022). Residues are matched
//! across models by chain, number and insertion code, so the models must
//! number them alike, as the models of one file do.

use std::collections::{BTreeMap, HashMap, HashSet};

use crate::fixed_str::FixedStr;
use crate::interactions::{Interaction, InteractionKind};
use crate::rin::network::{Edge, Network, Node, ResidueKey, build, degrees, residue_key};
use crate::structure::Protein;

/// The networks of an ensemble's models, folded together one model at a time
/// so the ensemble is never held whole.
#[derive(Clone, Debug, Default)]
pub struct Ensemble {
    models: usize,
    /// Each residue that interacts in any model, as the first such model has
    /// it.
    residues: BTreeMap<ResidueKey, Node>,
    /// How many models give each residue each secondary-structure letter, in
    /// the order they first do. Every model counts, not only those the residue
    /// interacts in, so the letter does not hang on the interactions asked for.
    ss: HashMap<ResidueKey, Vec<(char, usize)>>,
    edges: BTreeMap<(ResidueKey, ResidueKey, InteractionKind), Tally>,
    /// How many models each residue pair interacts in, by any kind.
    pairs: BTreeMap<(ResidueKey, ResidueKey), usize>,
}

/// What an edge gathers across models.
#[derive(Clone, Debug)]
struct Tally {
    models: usize,
    count: usize,
    min_distance: f32,
    interchain: bool,
}

impl Ensemble {
    /// Fold in one more model: its structure and the interactions found in it.
    pub fn add(&mut self, protein: &Protein, interactions: &[Interaction]) {
        self.models += 1;

        let mut seen = HashSet::new();
        for atom in &protein.atoms {
            let key = residue_key(atom);
            if seen.insert(key) {
                let letters = self.ss.entry(key).or_default();
                let ss = atom.secondary_structure.code();
                match letters.iter_mut().find(|(letter, _)| *letter == ss) {
                    Some((_, models)) => *models += 1,
                    None => letters.push((ss, 1)),
                }
            }
        }

        let network = build(protein, interactions);
        let key: HashMap<&str, ResidueKey> = network
            .nodes
            .iter()
            .map(|node| (node.id.as_str(), key_of(node)))
            .collect();
        for node in &network.nodes {
            self.residues
                .entry(key_of(node))
                .or_insert_with(|| node.clone());
        }
        let mut paired = HashSet::new();
        for edge in &network.edges {
            let pair = (key[edge.source.as_str()], key[edge.target.as_str()]);
            if paired.insert(pair) {
                *self.pairs.entry(pair).or_default() += 1;
            }
            self.edges
                .entry((pair.0, pair.1, edge.kind))
                .and_modify(|tally| {
                    tally.models += 1;
                    tally.count += edge.count;
                    tally.min_distance = tally.min_distance.min(edge.min_distance);
                })
                .or_insert(Tally {
                    models: 1,
                    count: edge.count,
                    min_distance: edge.min_distance,
                    interchain: edge.interchain,
                });
        }
    }

    /// The ensemble's network: every edge present in at least `at_least` of
    /// the models (the exact fraction, from 0, which keeps every edge, to 1),
    /// and the residues those edges join. An edge's count and energy are
    /// summed over the models, its distance the shortest in any; a residue's
    /// secondary structure is the one most models give it.
    pub fn network(&self, at_least: f32) -> Network {
        let edges: Vec<Edge> = self
            .edges
            .iter()
            .filter(|(_, tally)| tally.models as f32 / self.models as f32 >= at_least)
            .map(|((source, target, kind), tally)| Edge {
                source: self.residues[source].id.clone(),
                target: self.residues[target].id.clone(),
                kind: *kind,
                count: tally.count,
                min_distance: tally.min_distance,
                interchain: tally.interchain,
                energy: kind.contact_energy_kj() * tally.count as f32,
                models: tally.models,
            })
            .collect();

        let degrees = degrees(&edges);
        let nodes = self
            .residues
            .iter()
            .filter_map(|(key, node)| {
                let degree = *degrees.get(node.id.as_str())?;
                // The first of the most common letters: `max_by_key` keeps
                // the last of equals, so the letters are walked backwards.
                let ss = self.ss[key]
                    .iter()
                    .rev()
                    .max_by_key(|(_, models)| *models)
                    .map_or(node.ss, |&(ss, _)| ss);

                Some(Node {
                    ss,
                    degree,
                    ..node.clone()
                })
            })
            .collect();

        Network {
            nodes,
            edges,
            models: self.models,
        }
    }
}

impl Ensemble {
    /// Each residue pair that interacts in any model, as the ids of its two
    /// residues, and how many models it interacts in by any kind: a pair
    /// hydrogen-bonded in some models and in hydrophobic contact in the rest
    /// interacts in every one, which no single edge of the network says.
    pub fn contacts(&self) -> impl Iterator<Item = (&str, &str, usize)> {
        self.pairs.iter().map(|((source, target), &models)| {
            (
                self.residues[source].id.as_str(),
                self.residues[target].id.as_str(),
                models,
            )
        })
    }
}

/// The key [`crate::rin::build`] gives the residue of `node`; its chain is
/// already one an atom holds, so nothing is truncated.
fn key_of(node: &Node) -> ResidueKey {
    (FixedStr::new(&node.chain), node.resi, node.icode)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Atom, SecondaryStructure};
    use glam::Vec3;

    use InteractionKind::{HydrogenBond, Hydrophobic};
    use SecondaryStructure::{Coil, Helix, Sheet};

    fn atom(chain: &str, seq: i32, ss: SecondaryStructure) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new("CA"),
            element: FixedStr::new("C"),
            residue_name: FixedStr::new("ALA"),
            residue_seq: seq,
            chain_id: FixedStr::new(chain),
            is_hetatm: false,
            molecule: crate::structure::Molecule::Unknown,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: Vec3::ZERO,
            b_factor: 0.0,
            secondary_structure: ss,
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

    /// A/1, A/2 and B/5, one atom each, with A/1's secondary structure `ss`.
    fn residues(ss: SecondaryStructure) -> Protein {
        protein(vec![
            atom("A", 1, ss),
            atom("A", 2, Helix),
            atom("B", 5, Sheet),
        ])
    }

    /// Two models: A/1 and A/2 hydrogen-bonded in both, twice in the first;
    /// A/1 and B/5 in hydrophobic contact in the first alone.
    fn two_models() -> Ensemble {
        let mut ensemble = Ensemble::default();
        ensemble.add(
            &residues(Helix),
            &[
                interaction(HydrogenBond, 0, 1, 2.8),
                interaction(HydrogenBond, 0, 1, 3.1),
                interaction(Hydrophobic, 0, 2, 3.9),
            ],
        );
        ensemble.add(&residues(Helix), &[interaction(HydrogenBond, 0, 1, 3.0)]);

        ensemble
    }

    fn find(network: &Network, kind: InteractionKind) -> &Edge {
        network
            .edges
            .iter()
            .find(|edge| edge.kind == kind)
            .expect("an edge of the kind")
    }

    #[test]
    fn an_edge_is_weighted_by_the_fraction_of_models_it_is_present_in() {
        let network = two_models().network(0.0);

        assert_eq!(network.models, 2);
        let hbond = find(&network, HydrogenBond);
        assert_eq!(hbond.models, 2);
        assert_eq!(network.frequency(hbond), 1.0);
        assert_eq!(hbond.count, 3, "summed over the models");
        assert_eq!(hbond.min_distance, 2.8, "the shortest in any model");
        assert!((hbond.energy - 3.0 * HydrogenBond.contact_energy_kj()).abs() < 1e-4);
        let contact = find(&network, Hydrophobic);
        assert_eq!(network.frequency(contact), 0.5);
        assert!(contact.interchain);
        assert_eq!(network.nodes.len(), 3);
    }

    #[test]
    fn a_model_without_the_edge_counts_against_it() {
        let mut ensemble = two_models();
        ensemble.add(&residues(Helix), &[]);

        let network = ensemble.network(0.0);

        assert_eq!(network.models, 3);
        assert_eq!(find(&network, HydrogenBond).models, 2);
        assert!((network.frequency(find(&network, Hydrophobic)) - 1.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn a_threshold_keeps_the_frequent_edges_and_the_residues_they_join() {
        let network = two_models().network(1.0);

        assert_eq!(network.edges.len(), 1);
        assert_eq!(network.edges[0].kind, HydrogenBond);
        let ids: Vec<&str> = network.nodes.iter().map(|node| node.id.as_str()).collect();
        assert_eq!(ids, ["A/1", "A/2"]);
        assert!(network.nodes.iter().all(|node| node.degree == 1));
    }

    #[test]
    fn residues_sort_by_number() {
        let mut ensemble = Ensemble::default();
        ensemble.add(
            &protein(vec![atom("A", 9, Helix), atom("A", 10, Helix)]),
            &[interaction(HydrogenBond, 0, 1, 2.9)],
        );

        let network = ensemble.network(0.0);

        let ids: Vec<&str> = network.nodes.iter().map(|node| node.id.as_str()).collect();
        assert_eq!(ids, ["A/9", "A/10"]);
    }

    // A/1 interacts only in the first model, as coil; the two models after
    // give it a helix without its interactions, and the helix still wins.
    #[test]
    fn a_residue_takes_the_secondary_structure_most_models_give_it() {
        let mut ensemble = Ensemble::default();
        ensemble.add(&residues(Coil), &[interaction(Hydrophobic, 0, 2, 3.9)]);
        ensemble.add(&residues(Helix), &[]);
        ensemble.add(&residues(Helix), &[]);

        let network = ensemble.network(0.0);

        assert_eq!(network.nodes[0].id, "A/1");
        assert_eq!(network.nodes[0].ss, 'H');
    }

    #[test]
    fn of_secondary_structures_as_common_the_first_seen_wins() {
        let mut ensemble = Ensemble::default();
        ensemble.add(&residues(Coil), &[interaction(Hydrophobic, 0, 2, 3.9)]);
        ensemble.add(&residues(Helix), &[interaction(Hydrophobic, 0, 2, 3.9)]);

        assert_eq!(ensemble.network(0.0).nodes[0].ss, 'C');
    }

    #[test]
    fn a_pair_touching_by_one_kind_or_another_touches_in_every_model() {
        let mut ensemble = Ensemble::default();
        ensemble.add(&residues(Helix), &[interaction(HydrogenBond, 0, 1, 2.9)]);
        ensemble.add(
            &residues(Helix),
            &[
                interaction(HydrogenBond, 0, 1, 2.9),
                interaction(Hydrophobic, 0, 1, 3.9),
            ],
        );
        ensemble.add(&residues(Helix), &[interaction(Hydrophobic, 0, 1, 3.9)]);
        ensemble.add(&residues(Helix), &[]);

        let contacts: Vec<_> = ensemble.contacts().collect();

        assert_eq!(contacts, [("A/1", "A/2", 3)]);
        let network = ensemble.network(0.0);
        assert!(network.edges.iter().all(|edge| edge.models == 2));
    }

    #[test]
    fn an_ensemble_of_one_model_is_that_model_s_network() {
        let structure = residues(Helix);
        let interactions = [
            interaction(HydrogenBond, 0, 1, 2.8),
            interaction(Hydrophobic, 0, 2, 3.9),
        ];
        let mut ensemble = Ensemble::default();
        ensemble.add(&structure, &interactions);

        assert_eq!(ensemble.network(0.0), build(&structure, &interactions));
    }
}
