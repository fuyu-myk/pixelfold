//! The rings of a molecule's bond graph, as RDKit perceives them: the
//! smallest set of smallest rings, made symmetric, so a cage none of whose
//! faces is smaller than another keeps every face (cubane's six, where a
//! smallest set holds five).
//!
//! A ring belongs when no set of strictly smaller rings sums to it, bond by
//! bond (the relevant cycles). Every such ring is a cycle of Horton's set: from
//! each atom, the shortest paths to the two ends of a bond, closed by the bond.

use std::collections::{HashMap, HashSet, VecDeque};

/// The rings of a bond graph and which hold each atom and bond.
#[derive(Clone, Debug, Default)]
pub(crate) struct Rings {
    /// Each ring's atoms in order around it, smallest rings first.
    rings: Vec<Vec<usize>>,
    /// The rings holding each atom, by index into `rings`.
    by_atom: Vec<Vec<usize>>,
}

impl Rings {
    /// The rings of the graph whose atoms are bonded as `neighbours` says.
    pub(crate) fn of(neighbours: &[Vec<usize>]) -> Rings {
        let edges = edges(neighbours);
        let mut candidates: Vec<Vec<usize>> = cyclic_pieces(neighbours)
            .iter()
            .flat_map(|piece| {
                let (local, atoms) = induced(neighbours, piece);
                horton(&local, &self::edges(&local))
                    .into_iter()
                    .map(move |ring| ring.into_iter().map(|at| atoms[at]).collect::<Vec<usize>>())
            })
            .collect();
        candidates.sort_by_key(Vec::len);

        let mut basis: Vec<Vec<u64>> = Vec::new();
        let mut rings = Vec::new();
        let mut start = 0;
        while start < candidates.len() {
            let size = candidates[start].len();
            let end = candidates[start..]
                .iter()
                .position(|ring| ring.len() != size)
                .map_or(candidates.len(), |offset| start + offset);
            let shorter = basis.clone();
            for ring in &candidates[start..end] {
                let bits = bond_bits(ring, &edges);
                if independent(&shorter, &bits) {
                    rings.push(ring.clone());
                }
                if independent(&basis, &bits) {
                    insert(&mut basis, bits);
                }
            }
            start = end;
        }

        let mut by_atom = vec![Vec::new(); neighbours.len()];
        for (index, ring) in rings.iter().enumerate() {
            for &atom in ring {
                by_atom[atom].push(index);
            }
        }

        Rings { rings, by_atom }
    }

    /// Each ring's atoms in order around it, smallest rings first.
    pub(crate) fn all(&self) -> &[Vec<usize>] {
        &self.rings
    }

    /// The rings holding `atom`, by index into [`Rings::all`].
    pub(crate) fn of_atom(&self, atom: usize) -> &[usize] {
        &self.by_atom[atom]
    }

    /// The size of the smallest ring holding `atom`, if any does.
    pub(crate) fn smallest_of_atom(&self, atom: usize) -> Option<usize> {
        self.by_atom[atom]
            .iter()
            .map(|&ring| self.rings[ring].len())
            .min()
    }

    /// The rings holding the bond between `a` and `b`.
    pub(crate) fn of_bond(&self, a: usize, b: usize) -> impl Iterator<Item = usize> + '_ {
        self.by_atom[a]
            .iter()
            .copied()
            .filter(move |&ring| holds_bond(&self.rings[ring], a, b))
    }

    /// Whether some ring of `size` atoms holds `atom`.
    pub(crate) fn atom_in_ring_of(&self, atom: usize, size: usize) -> bool {
        self.by_atom[atom]
            .iter()
            .any(|&ring| self.rings[ring].len() == size)
    }

    /// Whether some ring of `size` atoms holds the bond between `a` and `b`.
    pub(crate) fn bond_in_ring_of(&self, a: usize, b: usize, size: usize) -> bool {
        self.of_bond(a, b)
            .any(|ring| self.rings[ring].len() == size)
    }
}

/// Whether `ring`, its atoms in order around it, holds the bond `a`-`b`.
pub(crate) fn holds_bond(ring: &[usize], a: usize, b: usize) -> bool {
    (0..ring.len()).any(|at| {
        let next = ring[(at + 1) % ring.len()];
        (ring[at] == a && next == b) || (ring[at] == b && next == a)
    })
}

/// The pieces a ring can lie in: the graph without its bridges (a bond no
/// ring holds, whose removal splits its piece) falls into pieces, and each of
/// three atoms or more holds rings. A protein's chain so leaves only its side
/// chains' rings to search.
fn cyclic_pieces(neighbours: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let bridges = bridges(neighbours);
    let mut piece = vec![usize::MAX; neighbours.len()];
    let mut pieces = Vec::new();
    for start in 0..neighbours.len() {
        if piece[start] != usize::MAX {
            continue;
        }
        let id = pieces.len();
        piece[start] = id;
        let mut members = vec![start];
        let mut at = 0;
        while at < members.len() {
            let atom = members[at];
            for &next in &neighbours[atom] {
                if piece[next] == usize::MAX && !bridges.contains(&(atom.min(next), atom.max(next)))
                {
                    piece[next] = id;
                    members.push(next);
                }
            }
            at += 1;
        }
        pieces.push(members);
    }

    pieces
        .into_iter()
        .filter(|members| members.len() >= 3)
        .collect()
}

/// The graph's bridges, lesser atom first, by Tarjan's low links, walked
/// without recursion so a long chain cannot overflow the stack.
fn bridges(neighbours: &[Vec<usize>]) -> HashSet<(usize, usize)> {
    let mut order = vec![usize::MAX; neighbours.len()];
    let mut low = vec![0; neighbours.len()];
    let mut found = HashSet::new();
    let mut next_order = 0;
    for root in 0..neighbours.len() {
        if order[root] != usize::MAX {
            continue;
        }
        order[root] = next_order;
        low[root] = next_order;
        next_order += 1;
        // Each frame: an atom, the atom it was reached from, and how many of
        // its neighbours are done.
        let mut stack = vec![(root, usize::MAX, 0)];
        while let Some(&mut (atom, parent, ref mut done)) = stack.last_mut() {
            if let Some(&next) = neighbours[atom].get(*done) {
                *done += 1;
                if next == parent {
                    continue;
                }
                if order[next] == usize::MAX {
                    order[next] = next_order;
                    low[next] = next_order;
                    next_order += 1;
                    stack.push((next, atom, 0));
                } else {
                    low[atom] = low[atom].min(order[next]);
                }
            } else {
                stack.pop();
                if parent != usize::MAX {
                    low[parent] = low[parent].min(low[atom]);
                    if low[atom] > order[parent] {
                        found.insert((atom.min(parent), atom.max(parent)));
                    }
                }
            }
        }
    }

    found
}

/// The graph induced on `atoms`, renumbered in their order, with the atoms.
fn induced(neighbours: &[Vec<usize>], atoms: &[usize]) -> (Vec<Vec<usize>>, Vec<usize>) {
    let index: HashMap<usize, usize> = atoms
        .iter()
        .enumerate()
        .map(|(at, &atom)| (atom, at))
        .collect();
    let local = atoms
        .iter()
        .map(|atom| {
            neighbours[*atom]
                .iter()
                .filter_map(|other| index.get(other).copied())
                .collect()
        })
        .collect();

    (local, atoms.to_vec())
}

/// Every bond once, lesser atom first, in order.
fn edges(neighbours: &[Vec<usize>]) -> Vec<(usize, usize)> {
    let mut edges: Vec<(usize, usize)> = neighbours
        .iter()
        .enumerate()
        .flat_map(|(a, bonded)| bonded.iter().filter(move |&&b| b > a).map(move |&b| (a, b)))
        .collect();
    edges.sort_unstable();
    edges.dedup();
    edges
}

/// Horton's candidate rings, each once: for each root and each bond, the
/// shortest paths from the root to the bond's ends, closed by the bond, where
/// the two paths meet only at the root.
fn horton(neighbours: &[Vec<usize>], edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut seen: HashSet<Vec<usize>> = HashSet::new();
    let mut found = Vec::new();
    for root in 0..neighbours.len() {
        let parents = shortest_paths(root, neighbours);
        for &(x, y) in edges {
            let (Some(to_x), Some(to_y)) = (path(&parents, root, x), path(&parents, root, y))
            else {
                continue;
            };
            if to_x.contains(&y) || to_y.contains(&x) {
                continue;
            }
            let on_x: HashSet<usize> = to_x.iter().copied().collect();
            if to_y.iter().filter(|atom| on_x.contains(atom)).count() != 1 {
                continue;
            }
            let ring: Vec<usize> = to_x
                .iter()
                .copied()
                .chain(to_y.iter().rev().copied())
                .collect();
            if ring.len() < 3 {
                continue;
            }
            let ring = ring[..ring.len() - 1].to_vec();
            let mut key = ring.clone();
            key.sort_unstable();
            if seen.insert(key) {
                found.push(ring);
            }
        }
    }

    found
}

/// The breadth-first parent of each atom on a shortest path from `root`.
fn shortest_paths(root: usize, neighbours: &[Vec<usize>]) -> Vec<Option<usize>> {
    let mut parents = vec![None; neighbours.len()];
    let mut visited = vec![false; neighbours.len()];
    visited[root] = true;
    let mut queue = VecDeque::from([root]);
    while let Some(atom) = queue.pop_front() {
        for &next in &neighbours[atom] {
            if !visited[next] {
                visited[next] = true;
                parents[next] = Some(atom);
                queue.push_back(next);
            }
        }
    }

    parents
}

/// The path from `root` to `atom` through `parents`, root first, or `None`
/// when `atom` is in another piece.
fn path(parents: &[Option<usize>], root: usize, atom: usize) -> Option<Vec<usize>> {
    let mut path = vec![atom];
    let mut at = atom;
    while at != root {
        at = parents[at]?;
        path.push(at);
    }
    path.reverse();

    Some(path)
}

/// The bonds of `ring` as bits over `edges`.
fn bond_bits(ring: &[usize], edges: &[(usize, usize)]) -> Vec<u64> {
    let mut bits = vec![0u64; edges.len().div_ceil(64)];
    for at in 0..ring.len() {
        let (a, b) = (ring[at], ring[(at + 1) % ring.len()]);
        let edge = (a.min(b), a.max(b));
        let index = edges
            .binary_search(&edge)
            .expect("a ring's bonds are bonds of the graph");
        bits[index / 64] |= 1 << (index % 64);
    }

    bits
}

/// The lowest bit set, if any.
fn leading(bits: &[u64]) -> Option<usize> {
    bits.iter()
        .enumerate()
        .find(|(_, word)| **word != 0)
        .map(|(at, word)| at * 64 + word.trailing_zeros() as usize)
}

/// Whether `bits` is no sum of the `basis` rows, each kept reduced so that its
/// leading bit is set in no other row.
fn independent(basis: &[Vec<u64>], bits: &[u64]) -> bool {
    leading(&reduced(basis, bits)).is_some()
}

fn reduced(basis: &[Vec<u64>], bits: &[u64]) -> Vec<u64> {
    let mut rest = bits.to_vec();
    for row in basis {
        let lead = leading(row).expect("a basis row is never empty");
        if rest[lead / 64] & (1 << (lead % 64)) != 0 {
            for (word, other) in rest.iter_mut().zip(row) {
                *word ^= other;
            }
        }
    }

    rest
}

/// Add `bits`, independent of `basis`, keeping every row's leading bit its
/// own.
fn insert(basis: &mut Vec<Vec<u64>>, bits: Vec<u64>) {
    let row = reduced(basis, &bits);
    let lead = leading(&row).expect("an independent row reduces to something");
    for other in basis.iter_mut() {
        if other[lead / 64] & (1 << (lead % 64)) != 0 {
            for (word, bit) in other.iter_mut().zip(&row) {
                *word ^= bit;
            }
        }
    }
    basis.push(row);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A bond graph from its bonds.
    fn graph(atoms: usize, bonds: &[(usize, usize)]) -> Vec<Vec<usize>> {
        let mut neighbours = vec![Vec::new(); atoms];
        for &(a, b) in bonds {
            neighbours[a].push(b);
            neighbours[b].push(a);
        }
        neighbours
    }

    fn sizes(rings: &Rings) -> Vec<usize> {
        rings.all().iter().map(Vec::len).collect()
    }

    #[test]
    fn a_chain_has_no_rings() {
        let rings = Rings::of(&graph(4, &[(0, 1), (1, 2), (2, 3)]));

        assert!(rings.all().is_empty());
        assert_eq!(rings.smallest_of_atom(1), None);
    }

    #[test]
    fn benzene_is_one_ring_in_order() {
        let rings = Rings::of(&graph(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]));

        assert_eq!(sizes(&rings), [6]);
        let ring = &rings.all()[0];
        for at in 0..6 {
            let (a, b) = (ring[at], ring[(at + 1) % 6]);
            assert!(
                a.abs_diff(b) == 1 || a.abs_diff(b) == 5,
                "{ring:?} is out of order"
            );
        }
    }

    #[test]
    fn naphthalene_is_two_rings_and_not_its_rim() {
        // Fused at atoms 0 and 5.
        let bonds = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 0),
        ];
        let rings = Rings::of(&graph(10, &bonds));

        assert_eq!(sizes(&rings), [6, 6]);
        assert_eq!(rings.of_atom(0).len(), 2);
        assert_eq!(rings.of_bond(0, 5).count(), 2);
        assert_eq!(rings.of_bond(0, 1).count(), 1);
    }

    #[test]
    fn indole_s_fused_atoms_lie_in_a_five_membered_ring_first() {
        let bonds = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 0),
        ];
        let rings = Rings::of(&graph(9, &bonds));

        assert_eq!(sizes(&rings), [5, 6]);
        assert_eq!(rings.smallest_of_atom(0), Some(5));
        assert_eq!(rings.smallest_of_atom(2), Some(6));
        assert!(rings.bond_in_ring_of(0, 5, 5));
        assert!(rings.atom_in_ring_of(0, 6));
    }

    #[test]
    fn cubane_keeps_all_six_faces() {
        let bonds = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        let rings = Rings::of(&graph(8, &bonds));

        assert_eq!(sizes(&rings), [4; 6]);
    }

    #[test]
    fn norbornane_is_two_five_membered_rings() {
        // Bridgeheads 0 and 3; bridges 1-2, 4-5 and 6.
        let bonds = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (0, 6),
            (6, 3),
        ];
        let rings = Rings::of(&graph(7, &bonds));

        assert_eq!(sizes(&rings), [5, 5]);
    }

    #[test]
    fn bicyclooctane_keeps_its_three_equal_rings() {
        // Bridgeheads 0 and 3, three bridges of two.
        let bonds = [
            (0, 1),
            (1, 2),
            (2, 3),
            (0, 4),
            (4, 5),
            (5, 3),
            (0, 6),
            (6, 7),
            (7, 3),
        ];
        let rings = Rings::of(&graph(8, &bonds));

        assert_eq!(sizes(&rings), [6, 6, 6]);
    }

    #[test]
    fn rings_hanging_off_a_long_chain_are_each_found_once() {
        // A chain of 3000 atoms with a benzene ring on every tenth: only the
        // rings' atoms can lie in a ring, so only they are searched.
        let mut bonds: Vec<(usize, usize)> = (0..2999).map(|at| (at, at + 1)).collect();
        let mut atoms = 3000;
        for anchor in (0..3000).step_by(10) {
            let ring: Vec<usize> = (atoms..atoms + 6).collect();
            atoms += 6;
            bonds.push((anchor, ring[0]));
            for at in 0..6 {
                bonds.push((ring[at], ring[(at + 1) % 6]));
            }
        }
        let rings = Rings::of(&graph(atoms, &bonds));

        assert_eq!(sizes(&rings), [6; 300]);
        assert_eq!(rings.smallest_of_atom(0), None);
    }

    #[test]
    fn a_spiro_atom_lies_in_both_rings() {
        let bonds = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 5), (5, 0)];
        let rings = Rings::of(&graph(6, &bonds));

        assert_eq!(sizes(&rings), [3, 4]);
        assert_eq!(rings.of_atom(0).len(), 2);
    }
}
