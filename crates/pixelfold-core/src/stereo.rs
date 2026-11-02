//! The tetrahedral stereocentres and E/Z double bonds of a ligand, and
//! whether a pose keeps their handedness.
//!
//! A co-folding model can return a ligand with a stereocentre inverted, which
//! is another molecule however well it sits; PoseBusters counts this among its
//! checks (Buttenschoen, Morris and Deane, Chem. Sci. 15:3130, 2024). Here the
//! reference ligand's own geometry says which atoms are stereocentres and which
//! way each turns, and a pose whose atoms are paired with the reference's
//! inverts a centre when its substituents clearly turn the other way. A centre
//! the pose leaves flat is not counted: a covalent warhead is tetrahedral in the
//! reference only through its bond to the protein, and a pose without that bond
//! holds it planar. A double bond turns from E to Z, or back, when the pose's
//! substituents on its ends clearly sit on the other sides of it, which
//! PoseBusters also checks.

use std::collections::{BTreeSet, HashSet, VecDeque};

use glam::Vec3;

use crate::fixed_str::FixedStr;

/// The least volume the unit bonds from a centre to three of its neighbours
/// may span for the centre to be tetrahedral, about 0.19 A out of their plane
/// for bonds of 1.5 A: a regular tetrahedron spans 0.77, a planar centre none.
const TETRAHEDRAL: f32 = 0.35;

/// The elements that hold four substituents tetrahedrally. An oxygen or a
/// halogen with three neighbours is a broken bond graph, such as conformers
/// modelled as overlapping residues.
const CENTRES: &[&str] = &["C", "N", "P", "S", "B", "SI", "SE", "AS"];

/// The least cosine of the dihedral angle across a double bond for its
/// substituents to be cis (or, negated, trans), 60 degrees from either.
const PLANAR: f32 = 0.5;

/// The elements a double bond with E or Z joins.
const ENDS: &[&str] = &["C", "N"];

/// The fewest atoms a ring may hold for a double bond in it to be E or Z: in a
/// smaller ring a double bond can only be cis.
const RING_WITH_TRANS: usize = 8;

/// A stereocentre of a ligand: the atom, three of its neighbours, the sign of
/// the volume they span from it, and the atom's class under colour refinement.
/// For a double bond, one end, then a substituent on it, the other end and a
/// substituent on that, with the sign of the cosine of their dihedral angle,
/// and the lesser class of its ends.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Stereocentre {
    pub(crate) centre: usize,
    around: [usize; 3],
    turn: f32,
    pub(crate) class: usize,
    kind: Kind,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Tetrahedral,
    DoubleBond,
}

impl Stereocentre {
    /// The centre and the three neighbours its turn is read from; for a double
    /// bond, the four atoms of its dihedral, the first end first.
    pub(crate) fn atoms(&self) -> [usize; 4] {
        let [a, b, c] = self.around;
        [self.centre, a, b, c]
    }

    /// The ends of the double bond this is, or `None` for a tetrahedral centre.
    pub(crate) fn bond(&self) -> Option<[usize; 2]> {
        (self.kind == Kind::DoubleBond).then_some([self.centre, self.around[1]])
    }

    /// Whether atoms at `positions`, those of [`Stereocentre::atoms`] in
    /// order, clearly turn the other way; a centre they leave flat, or a
    /// double bond they leave twisted across, is not.
    pub(crate) fn inverted_at(&self, positions: [Vec3; 4]) -> bool {
        match self.kind {
            Kind::Tetrahedral => volume(positions) * self.turn <= -TETRAHEDRAL,
            Kind::DoubleBond => dihedral_cosine(positions) * self.turn <= -PLANAR,
        }
    }
}

/// The stereocentres of a ligand of heavy atoms with these `elements`,
/// `positions` and bonded `neighbours`: each tetrahedral atom with three or
/// four neighbours no symmetry of the bond graph can exchange. A nitrogen with
/// three neighbours is left out, since it inverts freely.
pub(crate) fn stereocentres(
    elements: &[FixedStr<2>],
    positions: &[Vec3],
    neighbours: &[Vec<usize>],
) -> Vec<Stereocentre> {
    let classes = classes(elements, neighbours);

    (0..elements.len())
        .filter_map(|centre| {
            let bonded = &neighbours[centre];
            let element = elements[centre].as_str();
            let inverts = bonded.len() == 3 && element == "N";
            if !(3..=4).contains(&bonded.len()) || inverts || !CENTRES.contains(&element) {
                return None;
            }
            let distinct: HashSet<usize> = bonded.iter().map(|&atom| classes[atom]).collect();
            if distinct.len() < bonded.len() {
                return None;
            }

            let around = [bonded[0], bonded[1], bonded[2]];
            let spanned =
                volume([centre, around[0], around[1], around[2]].map(|atom| positions[atom]));
            (spanned.abs() >= TETRAHEDRAL).then_some(Stereocentre {
                centre,
                around,
                turn: spanned.signum(),
                class: classes[centre],
                kind: Kind::Tetrahedral,
            })
        })
        .collect()
}

/// The `doubles` bonds of a ligand that are E or Z, its heavy atoms with these
/// `elements`, `positions` and bonded `neighbours`: each joining carbons or
/// nitrogens outside any ring of fewer than eight atoms, cis or trans in the
/// reference, with a substituent on either end that no symmetry of the bond
/// graph can exchange with another there.
pub(crate) fn double_bonds(
    elements: &[FixedStr<2>],
    positions: &[Vec3],
    neighbours: &[Vec<usize>],
    doubles: &[(usize, usize)],
) -> Vec<Stereocentre> {
    let classes = classes(elements, neighbours);
    let substituent = |end: usize, other: usize| {
        let around: Vec<usize> = neighbours[end]
            .iter()
            .copied()
            .filter(|&atom| atom != other)
            .collect();
        match around[..] {
            [only] => Some(only),
            [first, second] if classes[first] != classes[second] => Some(first),
            _ => None,
        }
    };

    doubles
        .iter()
        .filter_map(|&(first, second)| {
            if !ENDS.contains(&elements[first].as_str())
                || !ENDS.contains(&elements[second].as_str())
                || in_small_ring(neighbours, first, second)
            {
                return None;
            }
            let around = [
                substituent(first, second)?,
                second,
                substituent(second, first)?,
            ];
            let cosine = dihedral_cosine(
                [first, around[0], around[1], around[2]].map(|atom| positions[atom]),
            );

            (cosine.abs() >= PLANAR).then_some(Stereocentre {
                centre: first,
                around,
                turn: cosine.signum(),
                class: classes[first].min(classes[second]),
                kind: Kind::DoubleBond,
            })
        })
        .collect()
}

/// Whether the bond from `first` to `second` closes a ring of fewer than
/// [`RING_WITH_TRANS`] atoms.
fn in_small_ring(neighbours: &[Vec<usize>], first: usize, second: usize) -> bool {
    let mut hops = vec![None; neighbours.len()];
    hops[first] = Some(0);
    let mut queue = VecDeque::from([first]);
    while let Some(atom) = queue.pop_front() {
        let next = hops[atom].map_or(0, |count: usize| count + 1);
        if next >= RING_WITH_TRANS - 1 {
            continue;
        }
        for &bonded in &neighbours[atom] {
            if atom == first && bonded == second {
                continue;
            }
            if bonded == second {
                return true;
            }
            if hops[bonded].is_none() {
                hops[bonded] = Some(next);
                queue.push_back(bonded);
            }
        }
    }

    false
}

/// The cosine of the dihedral angle of `a` and `d` about the bond from
/// `first` to `second`, at positions `[first, a, second, d]`.
fn dihedral_cosine([first, a, second, d]: [Vec3; 4]) -> f32 {
    let axis = (second - first).normalize_or_zero();
    let across = |from: Vec3, to: Vec3| {
        let bond = to - from;
        (bond - axis * bond.dot(axis)).normalize_or_zero()
    };

    across(first, a).dot(across(second, d))
}

/// The volume the unit bonds from a centre to three neighbours span, signed by
/// the way they turn.
fn volume([centre, a, b, c]: [Vec3; 4]) -> f32 {
    let bond = |atom: Vec3| (atom - centre).normalize_or_zero();
    bond(a).dot(bond(b).cross(bond(c)))
}

/// Each atom's class under colour refinement of the bond graph, from element
/// and number of bonds. Atoms of different classes are never exchanged by a
/// symmetry of the graph, so a centre whose neighbours all differ in class is
/// a stereocentre; atoms of one class nearly always are exchanged, and a
/// centre the refinement cannot tell apart is passed over rather than
/// misjudged.
pub(crate) fn classes(elements: &[FixedStr<2>], neighbours: &[Vec<usize>]) -> Vec<usize> {
    let mut classes = relabelled(
        &elements
            .iter()
            .zip(neighbours)
            .map(|(&element, bonded)| (element, bonded.len()))
            .collect::<Vec<_>>(),
    );
    loop {
        let signatures: Vec<(usize, Vec<usize>)> = neighbours
            .iter()
            .zip(&classes)
            .map(|(bonded, &class)| {
                let mut around: Vec<usize> = bonded.iter().map(|&atom| classes[atom]).collect();
                around.sort_unstable();
                (class, around)
            })
            .collect();
        let refined = relabelled(&signatures);
        if count(&refined) == count(&classes) {
            return refined;
        }
        classes = refined;
    }
}

/// Each key's rank among the distinct keys.
fn relabelled<T: Ord>(keys: &[T]) -> Vec<usize> {
    let distinct: Vec<&T> = keys.iter().collect::<BTreeSet<_>>().into_iter().collect();
    keys.iter()
        .map(|key| {
            distinct
                .binary_search(&key)
                .expect("every key is among them")
        })
        .collect()
}

fn count(classes: &[usize]) -> usize {
    classes.iter().max().map_or(0, |&most| most + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The four corners of a regular tetrahedron around the origin, 1.5 A out.
    fn corners() -> [Vec3; 4] {
        [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
            Vec3::new(-1.0, -1.0, 1.0),
        ]
        .map(|corner| corner.normalize() * 1.5)
    }

    /// A centre (atom 0) bonded to each of `substituents`, which take the
    /// corners in turn; a substituent named twice is bonded to atom 0 in each
    /// place, and `extra` atoms hang off the substituent at the given index.
    fn star(centre: &str, substituents: &[&str], extra: &[(usize, &str)]) -> Graph {
        let mut elements = vec![FixedStr::new(centre)];
        let mut positions = vec![Vec3::ZERO];
        let mut neighbours = vec![Vec::new()];
        for (&element, corner) in substituents.iter().zip(corners()) {
            let atom = elements.len();
            elements.push(FixedStr::new(element));
            positions.push(corner);
            neighbours.push(vec![0]);
            neighbours[0].push(atom);
        }
        for &(on, element) in extra {
            let atom = elements.len();
            elements.push(FixedStr::new(element));
            positions.push(positions[on] * 2.0);
            neighbours.push(vec![on]);
            neighbours[on].push(atom);
        }

        Graph {
            elements,
            positions,
            neighbours,
        }
    }

    struct Graph {
        elements: Vec<FixedStr<2>>,
        positions: Vec<Vec3>,
        neighbours: Vec<Vec<usize>>,
    }

    impl Graph {
        fn stereocentres(&self) -> Vec<Stereocentre> {
            stereocentres(&self.elements, &self.positions, &self.neighbours)
        }
    }

    #[test]
    fn a_carbon_with_four_different_neighbours_is_a_stereocentre() {
        let graph = star("C", &["C", "N", "O", "F"], &[]);

        let found = graph.stereocentres();

        assert_eq!(found.len(), 1);
        assert_eq!(found[0].centre, 0);
        let at = |positions: &[Vec3]| found[0].atoms().map(|atom| positions[atom]);
        assert!(!found[0].inverted_at(at(&graph.positions)));
        let mirrored: Vec<Vec3> = graph
            .positions
            .iter()
            .map(|&p| Vec3::new(-p.x, p.y, p.z))
            .collect();
        assert!(found[0].inverted_at(at(&mirrored)));
    }

    #[test]
    fn a_centre_a_pose_leaves_flat_is_not_inverted() {
        let graph = star("C", &["C", "N", "O"], &[]);
        let found = graph.stereocentres();
        // The centre pushed through the plane of its neighbours by a little,
        // as a warhead carbon without its bond to the protein sits.
        let mut flat = graph.positions.clone();
        let plane = (flat[1] + flat[2] + flat[3]) / 3.0;
        flat[0] = plane + (plane - graph.positions[0]).normalize() * 0.05;

        assert!(!found[0].inverted_at(found[0].atoms().map(|atom| flat[atom])));
    }

    #[test]
    fn a_carbon_with_one_hydrogen_left_implicit_is_a_stereocentre() {
        let graph = star("C", &["C", "N", "O"], &[]);
        assert_eq!(graph.stereocentres().len(), 1);
    }

    #[test]
    fn neighbours_the_graph_can_exchange_make_no_stereocentre() {
        // Two methyls on one carbon, and two ethyls alike all the way out.
        let methyls = star("C", &["C", "C", "N", "O"], &[]);
        let ethyls = star("C", &["C", "C", "N", "O"], &[(1, "C"), (2, "C")]);
        // An ethyl beside a methyl differs one bond out.
        let ethyl_and_methyl = star("C", &["C", "C", "N", "O"], &[(1, "C")]);

        assert!(methyls.stereocentres().is_empty());
        assert!(ethyls.stereocentres().is_empty());
        assert_eq!(ethyl_and_methyl.stereocentres().len(), 1);
    }

    #[test]
    fn a_planar_centre_a_three_bonded_nitrogen_and_an_oxygen_are_not_stereocentres() {
        let mut planar = star("C", &["C", "N", "O"], &[]);
        for position in &mut planar.positions {
            position.z = 0.0;
        }
        let amine = star("N", &["C", "S", "O"], &[]);
        let broken = star("O", &["C", "S", "N"], &[]);

        assert!(planar.stereocentres().is_empty());
        assert!(amine.stereocentres().is_empty());
        assert!(broken.stereocentres().is_empty());
    }

    /// But-2-ene laid flat: atom 0 on C1, C1=C2, atom 3 on C2, trans as
    /// written; `d` places atom 3.
    fn butene(d: Vec3) -> Graph {
        Graph {
            elements: ["C", "C", "C", "C"].map(FixedStr::new).to_vec(),
            positions: vec![
                Vec3::new(-0.75, 1.3, 0.0),
                Vec3::ZERO,
                Vec3::new(1.34, 0.0, 0.0),
                d,
            ],
            neighbours: vec![vec![1], vec![0, 2], vec![1, 3], vec![2]],
        }
    }

    const TRANS: Vec3 = Vec3::new(2.09, -1.3, 0.0);

    #[test]
    fn a_double_bond_turned_from_e_to_z_is_inverted() {
        let graph = butene(TRANS);
        let found = double_bonds(
            &graph.elements,
            &graph.positions,
            &graph.neighbours,
            &[(1, 2)],
        );
        let at = |positions: &[Vec3]| found[0].atoms().map(|atom| positions[atom]);
        let mut cis = graph.positions.clone();
        cis[3] = Vec3::new(2.09, 1.3, 0.0);
        let mut twisted = graph.positions.clone();
        twisted[3] = Vec3::new(2.09, 0.0, 1.3);

        assert_eq!(found.len(), 1);
        assert_eq!(found[0].bond(), Some([1, 2]));
        assert!(!found[0].inverted_at(at(&graph.positions)));
        assert!(found[0].inverted_at(at(&cis)));
        assert!(!found[0].inverted_at(at(&twisted)), "neither E nor Z");
    }

    #[test]
    fn a_double_bond_end_holding_two_alike_substituents_or_none_has_no_e_or_z() {
        let mut isobutene = butene(TRANS);
        isobutene.elements.push(FixedStr::new("C"));
        isobutene.positions.push(Vec3::new(-0.75, -1.3, 0.0));
        isobutene.neighbours[1].push(4);
        isobutene.neighbours.push(vec![1]);
        let mut propene = butene(TRANS);
        propene.neighbours[2].retain(|&atom| atom != 3);
        propene.neighbours[3].clear();

        for graph in [isobutene, propene] {
            let found = double_bonds(
                &graph.elements,
                &graph.positions,
                &graph.neighbours,
                &[(1, 2)],
            );
            assert!(found.is_empty());
        }
    }

    #[test]
    fn a_double_bond_in_a_small_ring_or_twisted_in_the_reference_has_no_e_or_z() {
        // Cyclohexene: the ring holds C1=C2 cis.
        let ring: Vec<Vec3> = (0..6)
            .map(|at| {
                let angle = at as f32 * std::f32::consts::TAU / 6.0;
                Vec3::new(angle.cos(), angle.sin(), 0.0) * 1.4
            })
            .collect();
        let cyclohexene = Graph {
            elements: vec![FixedStr::new("C"); 6],
            positions: ring,
            neighbours: (0..6).map(|at| vec![(at + 5) % 6, (at + 1) % 6]).collect(),
        };
        let twisted = butene(Vec3::new(2.09, 0.0, 1.3));

        for (graph, double) in [(cyclohexene, (0, 1)), (twisted, (1, 2))] {
            let found = double_bonds(
                &graph.elements,
                &graph.positions,
                &graph.neighbours,
                &[double],
            );
            assert!(found.is_empty());
        }
    }
}
