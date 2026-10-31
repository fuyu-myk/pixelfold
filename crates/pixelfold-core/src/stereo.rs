//! The tetrahedral stereocentres of a ligand, and whether a pose keeps their
//! handedness.
//!
//! A co-folding model can return a ligand with a stereocentre inverted, which
//! is another molecule however well it sits; PoseBusters counts this among its
//! checks (Buttenschoen, Morris and Deane, Chem. Sci. 15:3130, 2024). Here the
//! reference ligand's own geometry says which atoms are stereocentres and which
//! way each turns, and a pose whose atoms are paired with the reference's
//! inverts a centre when its substituents clearly turn the other way. A centre
//! the pose leaves flat is not counted: a covalent warhead is tetrahedral in the
//! reference only through its bond to the protein, and a pose without that bond
//! holds it planar.

use std::collections::{BTreeSet, HashSet};

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

/// A stereocentre of a ligand: the atom, three of its neighbours, the sign of
/// the volume they span from it, and the atom's class under colour refinement.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Stereocentre {
    pub(crate) centre: usize,
    around: [usize; 3],
    turn: f32,
    pub(crate) class: usize,
}

impl Stereocentre {
    /// The centre and the three neighbours its turn is read from.
    pub(crate) fn atoms(&self) -> [usize; 4] {
        let [a, b, c] = self.around;
        [self.centre, a, b, c]
    }

    /// Whether atoms at `positions`, those of [`Stereocentre::atoms`] in
    /// order, clearly turn the other way; a centre they leave flat is not.
    pub(crate) fn inverted_at(&self, positions: [Vec3; 4]) -> bool {
        volume(positions) * self.turn <= -TETRAHEDRAL
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
            })
        })
        .collect()
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
fn classes(elements: &[FixedStr<2>], neighbours: &[Vec<usize>]) -> Vec<usize> {
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
}
