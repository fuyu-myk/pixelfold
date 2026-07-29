//! Structural analysis of a residue interaction network.
//!
//! Connected components, betweenness centrality, and articulation points read
//! the network as an undirected simple graph with the interaction kinds joining a
//! residue pair collapsing to one connection.

use std::collections::{BTreeSet, HashMap, VecDeque};

use crate::rin::network::Network;

/// The structural analysis of a [`Network`]. Every per-node vector is aligned
/// with `network.nodes` by index.
#[derive(Clone, Debug, PartialEq)]
pub struct RinAnalysis {
    /// Connected components as lists of node indices, largest first.
    pub components: Vec<Vec<usize>>,
    /// Betweenness centrality per node: the number of shortest paths between
    /// other residues that pass through this one, summed over all pairs
    /// (undirected, so each pair is counted once).
    pub betweenness: Vec<f64>,
    /// Node indices whose removal splits their connected component (cut
    /// residues).
    pub articulation_points: Vec<usize>,
}

/// Analyse the structure of `network`.
pub fn analyze(network: &Network) -> RinAnalysis {
    let adj = adjacency(network);
    let mut components = connected_components(&adj);
    components.sort_by_key(|component| std::cmp::Reverse(component.len()));

    RinAnalysis {
        betweenness: betweenness(&adj),
        articulation_points: articulation_points(&adj),
        components,
    }
}

/// The undirected simple-graph adjacency of the network, indexed by node
/// position. Parallel edges (a residue pair joined by several interaction kinds)
/// collapse to one neighbour, and a self-edge is ignored.
fn adjacency(network: &Network) -> Vec<Vec<usize>> {
    let index: HashMap<&str, usize> = network
        .nodes
        .iter()
        .enumerate()
        .map(|(i, node)| (node.id.as_str(), i))
        .collect();

    let mut sets = vec![BTreeSet::new(); network.nodes.len()];
    for edge in &network.edges {
        if let (Some(&u), Some(&v)) = (
            index.get(edge.source.as_str()),
            index.get(edge.target.as_str()),
        ) && u != v
        {
            sets[u].insert(v);
            sets[v].insert(u);
        }
    }

    sets.into_iter()
        .map(|set| set.into_iter().collect())
        .collect()
}

/// Label each node with its connected component through breadth-first search.
fn connected_components(adj: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let n = adj.len();
    let mut seen = vec![false; n];
    let mut components = Vec::new();

    for start in 0..n {
        if seen[start] {
            continue;
        }

        let mut component = Vec::new();
        let mut queue = VecDeque::from([start]);
        seen[start] = true;
        while let Some(node) = queue.pop_front() {
            component.push(node);
            for &next in &adj[node] {
                if !seen[next] {
                    seen[next] = true;
                    queue.push_back(next);
                }
            }
        }

        components.push(component);
    }

    components
}

/// Betweenness centrality by Brandes' algorithm on the unweighted graph. The
/// raw accumulation counts each shortest path in both directions, so the result
/// is halved for the undirected network.
fn betweenness(adj: &[Vec<usize>]) -> Vec<f64> {
    let n = adj.len();
    let mut centrality = vec![0.0f64; n];

    for source in 0..n {
        let mut stack = Vec::new();
        let mut predecessors: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut sigma = vec![0.0f64; n];
        let mut distance = vec![-1i64; n];
        sigma[source] = 1.0;
        distance[source] = 0;

        let mut queue = VecDeque::from([source]);
        while let Some(v) = queue.pop_front() {
            stack.push(v);
            for &w in &adj[v] {
                if distance[w] < 0 {
                    distance[w] = distance[v] + 1;
                    queue.push_back(w);
                }
                if distance[w] == distance[v] + 1 {
                    sigma[w] += sigma[v];
                    predecessors[w].push(v);
                }
            }
        }

        let mut delta = vec![0.0f64; n];
        while let Some(w) = stack.pop() {
            for &v in &predecessors[w] {
                delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
            }
            if w != source {
                centrality[w] += delta[w];
            }
        }
    }

    for value in &mut centrality {
        *value /= 2.0;
    }

    centrality
}

/// The cut vertices of the graph, by the depth-first low-link method (Tarjan).
/// Recursion depth is bounded by the size of a connected component.
fn articulation_points(adj: &[Vec<usize>]) -> Vec<usize> {
    let n = adj.len();
    let mut visited = vec![false; n];
    let mut discovery = vec![0usize; n];
    let mut low = vec![0usize; n];
    let mut is_cut = vec![false; n];
    let mut timer = 0usize;

    for start in 0..n {
        if !visited[start] {
            dfs_cut(
                start,
                None,
                adj,
                &mut visited,
                &mut discovery,
                &mut low,
                &mut is_cut,
                &mut timer,
            );
        }
    }

    (0..n).filter(|&i| is_cut[i]).collect()
}

#[allow(clippy::too_many_arguments)]
fn dfs_cut(
    u: usize,
    parent: Option<usize>,
    adj: &[Vec<usize>],
    visited: &mut [bool],
    discovery: &mut [usize],
    low: &mut [usize],
    is_cut: &mut [bool],
    timer: &mut usize,
) {
    visited[u] = true;
    *timer += 1;
    discovery[u] = *timer;
    low[u] = *timer;
    let mut children = 0;

    for &v in &adj[u] {
        if !visited[v] {
            children += 1;
            dfs_cut(v, Some(u), adj, visited, discovery, low, is_cut, timer);
            low[u] = low[u].min(low[v]);
            // A non-root u is a cut vertex if a child cannot reach an ancestor
            // of u without going through u.
            if parent.is_some() && low[v] >= discovery[u] {
                is_cut[u] = true;
            }
        } else if Some(v) != parent {
            low[u] = low[u].min(discovery[v]);
        }
    }

    // The root is a cut vertex only if it roots more than one subtree.
    if parent.is_none() && children > 1 {
        is_cut[u] = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Undirected adjacency from an edge list.
    fn graph(n: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
        let mut sets = vec![BTreeSet::new(); n];
        for &(u, v) in edges {
            sets[u].insert(v);
            sets[v].insert(u);
        }
        sets.into_iter().map(|s| s.into_iter().collect()).collect()
    }

    #[test]
    fn betweenness_of_a_path_peaks_in_the_middle() {
        // 0 - 1 - 2 - 3: the two interior nodes each lie on the shortest path of
        // two residue pairs, the endpoints on none.
        let g = graph(4, &[(0, 1), (1, 2), (2, 3)]);
        let bc = betweenness(&g);
        assert_eq!(bc, vec![0.0, 2.0, 2.0, 0.0]);
    }

    #[test]
    fn betweenness_of_a_star_centre_is_every_leaf_pair() {
        // Centre 0, leaves 1..3: the centre lies on all three leaf-to-leaf paths.
        let g = graph(4, &[(0, 1), (0, 2), (0, 3)]);
        let bc = betweenness(&g);
        assert_eq!(bc, vec![3.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn components_separate_disconnected_subgraphs() {
        // Two triangles with no edge between them.
        let g = graph(6, &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]);
        let mut components = connected_components(&g);
        components.iter_mut().for_each(|c| c.sort_unstable());
        components.sort_by_key(|c| c[0]);
        assert_eq!(components, vec![vec![0, 1, 2], vec![3, 4, 5]]);
    }

    #[test]
    fn articulation_point_is_the_bridging_residue() {
        // A triangle 0-1-2 with a tail 2-3: removing 2 strands node 3.
        let g = graph(4, &[(0, 1), (1, 2), (2, 0), (2, 3)]);
        assert_eq!(articulation_points(&g), vec![2]);
    }

    #[test]
    fn a_path_cuts_at_every_interior_node() {
        let g = graph(4, &[(0, 1), (1, 2), (2, 3)]);
        assert_eq!(articulation_points(&g), vec![1, 2]);
    }

    #[test]
    fn a_cycle_has_no_cut_vertices() {
        let g = graph(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]);
        assert!(articulation_points(&g).is_empty());
    }

    #[test]
    fn an_empty_graph_analyses_to_nothing() {
        let g: Vec<Vec<usize>> = Vec::new();
        assert!(connected_components(&g).is_empty());
        assert!(betweenness(&g).is_empty());
        assert!(articulation_points(&g).is_empty());
    }
}
