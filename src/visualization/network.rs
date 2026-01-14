use crate::{Protein, SecondaryStructure};
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;
use std::collections::HashMap;


/// Residue node in the H-bond network graph
#[derive(Clone, Debug)]
pub struct ResidueNode {
    pub residue_seq: u32,
    pub chain_id: String,
    pub residue_name: String,
    pub secondary_structure: SecondaryStructure,
    pub atom_indices: Vec<usize>,
}

/// H-bond edge in the network graph
#[derive(Clone, Debug)]
pub struct HBondEdge {
    pub energy: f32,  // kcal/mol
    pub bond_type: BondType,
    pub donor_atom_idx: usize,
    pub acceptor_atom_idx: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondType {
    IntraChain,
    InterChain,
}

/// H-bond network graph with analysis capabilities
pub struct HBondGraph {
    pub graph: DiGraph<ResidueNode, HBondEdge>,
    residue_to_node: HashMap<(String, u32), NodeIndex>,  // (chain_id, residue_seq) -> NodeIndex
}

impl HBondGraph {
    /// Builds a H-bond network graph from protein structure
    pub fn build(protein: &Protein) -> Self {
        let mut graph = DiGraph::new();
        let mut residue_to_node = HashMap::new();
        
        // Group atoms by residue
        let mut residues: HashMap<(String, u32), Vec<usize>> = HashMap::new();
        for (idx, atom) in protein.atoms.iter().enumerate() {
            residues
                .entry((atom.chain_id.clone(), atom.residue_seq))
                .or_default()
                .push(idx);
        }
        
        // Create nodes for each residue
        for ((chain_id, residue_seq), atom_indices) in residues.iter() {
            let first_atom = &protein.atoms[atom_indices[0]];
            let node = ResidueNode {
                residue_seq: *residue_seq,
                chain_id: chain_id.clone(),
                residue_name: first_atom.residue_name.clone(),
                secondary_structure: first_atom.secondary_structure,
                atom_indices: atom_indices.clone(),
            };
            
            let node_idx = graph.add_node(node);
            residue_to_node.insert((chain_id.clone(), *residue_seq), node_idx);
        }
        
        // Add edges for H-bonds
        for hbond in &protein.hbonds {
            let donor_atom = &protein.atoms[hbond.donor_atom_idx];
            let acceptor_atom = &protein.atoms[hbond.acceptor_atom_idx];
            
            let donor_key = (donor_atom.chain_id.clone(), donor_atom.residue_seq);
            let acceptor_key = (acceptor_atom.chain_id.clone(), acceptor_atom.residue_seq);
            
            if let (Some(&donor_node), Some(&acceptor_node)) = 
                (residue_to_node.get(&donor_key), residue_to_node.get(&acceptor_key)) {
                
                let bond_type = if donor_atom.chain_id == acceptor_atom.chain_id {
                    BondType::IntraChain
                } else {
                    BondType::InterChain
                };
                
                let edge = HBondEdge {
                    energy: hbond.energy,
                    bond_type,
                    donor_atom_idx: hbond.donor_atom_idx,
                    acceptor_atom_idx: hbond.acceptor_atom_idx,
                };
                
                graph.add_edge(donor_node, acceptor_node, edge);
            }
        }
        
        HBondGraph {
            graph,
            residue_to_node,
        }
    }
    
    /// Filter graph by energy threshold (returns new filtered graph)
    pub fn filter_by_energy(&self, threshold: f32) -> Self {
        let mut new_graph = DiGraph::new();
        let mut old_to_new: HashMap<NodeIndex, NodeIndex> = HashMap::new();
        let mut new_residue_to_node = HashMap::new();
        
        // Copy all nodes
        for node_idx in self.graph.node_indices() {
            let node = &self.graph[node_idx];
            let new_idx = new_graph.add_node(node.clone());
            old_to_new.insert(node_idx, new_idx);
            new_residue_to_node.insert((node.chain_id.clone(), node.residue_seq), new_idx);
        }
        
        // Copy edges that pass the threshold
        for edge in self.graph.edge_references() {
            if edge.weight().energy < threshold {
                let source = old_to_new[&edge.source()];
                let target = old_to_new[&edge.target()];
                new_graph.add_edge(source, target, edge.weight().clone());
            }
        }
        
        HBondGraph {
            graph: new_graph,
            residue_to_node: new_residue_to_node,
        }
    }
    
    /// Get node index for a residue
    pub fn get_node(&self, chain_id: &str, residue_seq: u32) -> Option<NodeIndex> {
        self.residue_to_node.get(&(chain_id.to_string(), residue_seq)).copied()
    }
}

/// Network analysis results
#[derive(Clone, Debug)]
pub struct NetworkAnalysis {
    pub degree_centrality: HashMap<NodeIndex, usize>,
    pub connected_components: Vec<Vec<NodeIndex>>,
    pub motifs: HBondMotifs,
}

/// Detected H-bond motifs (secondary structure patterns)
#[derive(Clone, Debug, Default)]
pub struct HBondMotifs {
    pub helix_ladders: Vec<Vec<NodeIndex>>,      // i -> i+4 patterns
    pub sheet_ladders: Vec<Vec<NodeIndex>>,      // Non-local bidirectional patterns
    pub turns: Vec<NodeIndex>,                   // Single residue turns
}

impl HBondGraph {
    /// Compute degree centrality for all nodes
    pub fn compute_degree_centrality(&self) -> HashMap<NodeIndex, usize> {
        let mut centrality = HashMap::new();
        
        for node_idx in self.graph.node_indices() {
            let in_degree = self.graph.neighbors_directed(node_idx, petgraph::Direction::Incoming).count();
            let out_degree = self.graph.neighbors_directed(node_idx, petgraph::Direction::Outgoing).count();
            centrality.insert(node_idx, in_degree + out_degree);
        }
        
        centrality
    }
    
    /// Find connected components in the H-bond network
    /// 
    /// Uses depth-first search to find weakly connected components
    pub fn find_connected_components(&self) -> Vec<Vec<NodeIndex>> {
        use petgraph::visit::Visitable;
        use petgraph::visit::depth_first_search;
        use petgraph::visit::DfsEvent;
        
        let mut visited = self.graph.visit_map();
        let mut components = Vec::new();
        
        for start in self.graph.node_indices() {
            let start_idx = start.index();
            if visited.contains(start_idx) {
                continue;
            }
            
            let mut component = Vec::new();
            depth_first_search(&self.graph, Some(start), |event| {
                if let DfsEvent::Discover(n, _) = event {
                    let n_idx = n.index();
                    visited.insert(n_idx);
                    component.push(n);
                }
            });
            
            if !component.is_empty() {
                components.push(component);
            }
        }
        
        components
    }
    
    /// Detect H-bond motifs (helix ladders, sheet patterns, turns)
    pub fn detect_hbond_motifs(&self) -> HBondMotifs {
        let mut motifs = HBondMotifs::default();
        
        // Detect α-helix ladders (i → i+4 patterns)
        for edge in self.graph.edge_references() {
            let source = &self.graph[edge.source()];
            let target = &self.graph[edge.target()];
            
            if source.chain_id == target.chain_id {
                let offset = (target.residue_seq as i32 - source.residue_seq as i32).abs();
                
                if offset == 4 {
                    let mut ladder = vec![edge.source(), edge.target()];
                    
                    // Try to extend forward
                    let mut current = edge.target();
                    while let Some(next_edge) = self.graph.edges(current).find(|e| {
                        let next_target = &self.graph[e.target()];
                        next_target.chain_id == source.chain_id &&
                        (next_target.residue_seq as i32 - self.graph[current].residue_seq as i32).abs() == 4
                    }) {
                        ladder.push(next_edge.target());
                        current = next_edge.target();
                    }
                    
                    if ladder.len() >= 3 {
                        motifs.helix_ladders.push(ladder);
                    }
                }
                
                // Sheet patterns (non-local, offset > 4)
                if offset > 4 {
                    // Check for bidirectional H-bonds (antiparallel sheet)
                    let has_reverse = self.graph.edges(edge.target()).any(|e| e.target() == edge.source());
                    
                    if has_reverse {
                        let ladder = vec![edge.source(), edge.target()];
                        motifs.sheet_ladders.push(ladder);
                    }
                }
            }
        }
        
        // Detect turns (high connectivity in small sequence regions)
        let centrality = self.compute_degree_centrality();
        for (node_idx, &degree) in &centrality {
            if degree >= 3 {
                let node = &self.graph[*node_idx];
                if node.secondary_structure == SecondaryStructure::Turn {
                    motifs.turns.push(*node_idx);
                }
            }
        }
        
        motifs
    }
    
    /// Perform full network analysis
    pub fn analyze(&self) -> NetworkAnalysis {
        NetworkAnalysis {
            degree_centrality: self.compute_degree_centrality(),
            connected_components: self.find_connected_components(),
            motifs: self.detect_hbond_motifs(),
        }
    }
}
