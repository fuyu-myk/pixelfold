//! The checks PoseBusters makes of a ligand pose's own geometry
//! (Buttenschoen, Morris and Deane, Chem. Sci. 15:3130, 2024): its bond
//! lengths, its bond angles, its atoms clashing with one another, and the
//! flatness of its aromatic rings and double bonds.
//!
//! PoseBusters bounds each distance from RDKit's distance geometry, which
//! starts from a force field's ideal lengths and angles. Here the reference
//! ligand's own geometry stands in for the ideal one, under RDKit's margins
//! and PoseBusters' tolerances (its `redock` settings), so a pose is held to
//! the reference ligand's shape. Atoms three or more bonds apart are bounded as
//! RDKit bounds them before triangle smoothing, but never further apart than
//! the reference holds them, and a group must lie flat only where the
//! reference's does: a graph read from distances bonds no metal, and a
//! definition can call a bond double that the structure makes single.
//! Hydrogens are left out, as PoseBusters leaves them.

use std::collections::{HashMap, HashSet, VecDeque};

use glam::Vec3;
use nalgebra::{Matrix3, SymmetricEigen};

use crate::components::{BondOrder, Component};
use crate::fixed_str::FixedStr;
use crate::numbering::residue_key;
use crate::structure::Protein;

/// RDKit's margin around a bond length, in angstroms (`DIST12_DELTA`).
const BOND_MARGIN: f32 = 0.01;

/// RDKit's margin around a bond angle, in radians (`ANGLE_DELTA`).
const ANGLE_MARGIN: f32 = 0.035;

/// How far outside its bounds a bond length or the distance across a bond
/// angle may stray, as a fraction of the bound crossed.
const TOLERANCE: f32 = 0.25;

/// How far inside their lower bound two atoms may come, as a fraction of it.
const CLASH_TOLERANCE: f32 = 0.3;

/// How far from their shared plane the atoms of an aromatic ring or a double
/// bond may stray, in angstroms.
const FLATNESS: f32 = 0.25;

/// The share of their van der Waals radii summed that RDKit lets atoms four
/// bonds apart come within (`VDW_SCALE_15`); five bonds apart, halfway from it
/// to the whole sum.
const VDW_SCALE_15: f32 = 0.7;

/// What a pose's ligand fails of the checks of its own geometry.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct Geometry {
    /// Bonds more than a quarter longer or shorter than the reference's.
    pub bonds: usize,
    /// Bond angles stretching or closing the distance across them by more
    /// than a quarter.
    pub angles: usize,
    /// Pairs of atoms neither bonded nor sharing a neighbour that come within
    /// seven tenths of their lower bound.
    pub clashes: usize,
    /// Aromatic rings and double bonds whose atoms stray more than 0.25 A from
    /// a shared plane, or `None` when the reference ligand's definitions do not
    /// say which atoms are aromatic and which bonds double.
    pub unflat: Option<usize>,
}

/// What a ligand's definitions say of its heavy atoms, in order: which are
/// aromatic, and which of the bonds between them are double outside an
/// aromatic ring.
#[derive(Clone, Debug)]
pub(crate) struct Chemistry {
    pub(crate) aromatic: Vec<bool>,
    pub(crate) doubles: Vec<(usize, usize)>,
}

impl Chemistry {
    /// The chemistry of the `heavy` atoms of `protein`, bonded as `neighbours`
    /// says, or `None` when a definition is missing or lacks one of them.
    pub(crate) fn of(
        protein: &Protein,
        heavy: &[usize],
        neighbours: &[Vec<usize>],
    ) -> Option<Chemistry> {
        let atom = |at: usize| &protein.atoms[heavy[at]];
        let definitions: Vec<&Component> = (0..heavy.len())
            .map(|at| {
                let component = protein.components.get(atom(at).residue_name.as_str())?;
                component.element_of(atom(at).name.as_str())?;
                Some(component)
            })
            .collect::<Option<_>>()?;
        let aromatic: Vec<bool> = (0..heavy.len())
            .map(|at| definitions[at].is_aromatic(atom(at).name.as_str()))
            .collect();
        let mut doubles = Vec::new();
        for (a, bonded) in neighbours.iter().enumerate() {
            for &b in bonded.iter().filter(|&&b| b > a) {
                let double = residue_key(atom(a)) == residue_key(atom(b))
                    && definitions[a]
                        .bonded(atom(a).name.as_str())
                        .any(|(name, order)| {
                            name == atom(b).name.as_str() && matches!(order, BondOrder::Double)
                        });
                if double && !(aromatic[a] && aromatic[b]) {
                    doubles.push((a, b));
                }
            }
        }

        Some(Chemistry { aromatic, doubles })
    }

    /// Whether each atom shares its class under colour refinement, `classes`,
    /// with an atom of other chemistry (aromatic or not, or with another
    /// number of double bonds): a symmetry of the bond graph alone may then
    /// exchange the two, as one turns pent-2-ene-1,5-diol end for end, so a
    /// pairing by bond graph cannot say which is which.
    pub(crate) fn exchangeable(&self, classes: &[usize]) -> Vec<bool> {
        let kind = |atom: usize| {
            let doubled = self
                .doubles
                .iter()
                .filter(|&&(a, b)| a == atom || b == atom)
                .count();
            (self.aromatic[atom], doubled)
        };
        let mut kinds: HashMap<usize, HashSet<(bool, usize)>> = HashMap::new();
        for (atom, &class) in classes.iter().enumerate() {
            kinds.entry(class).or_default().insert(kind(atom));
        }

        classes.iter().map(|class| kinds[class].len() > 1).collect()
    }

    /// The groups of atoms that must lie flat: each aromatic ring, and each
    /// double bond joining carbons, neither of them a cumulene's, with the
    /// atoms bonded to it, the bonds PoseBusters' pattern matches; each
    /// `exchangeable` when one of its atoms is.
    pub(crate) fn flat(
        &self,
        elements: &[FixedStr<2>],
        neighbours: &[Vec<usize>],
        exchangeable: &[bool],
    ) -> Vec<Flat> {
        let doubled = |atom: usize| {
            self.doubles
                .iter()
                .filter(|&&(a, b)| a == atom || b == atom)
                .count()
        };
        let carbons: Vec<(usize, usize)> = self
            .doubles
            .iter()
            .copied()
            .filter(|&(a, b)| {
                [a, b]
                    .iter()
                    .all(|&atom| elements[atom].as_str() == "C" && doubled(atom) == 1)
            })
            .collect();

        flat_groups(neighbours, &self.aromatic, &carbons)
            .into_iter()
            .map(|atoms| Flat {
                exchangeable: atoms.iter().any(|&atom| exchangeable[atom]),
                atoms,
            })
            .collect()
    }
}

/// A group of atoms that must lie flat, and whether a symmetry of the bond
/// graph may carry it onto atoms of other chemistry.
#[derive(Clone, Debug)]
pub(crate) struct Flat {
    atoms: Vec<usize>,
    exchangeable: bool,
}

/// The bounds a pose's ligand is held to, taken from the reference ligand.
#[derive(Clone, Debug)]
pub(crate) struct Template {
    bonds: Vec<Bound>,
    angles: Vec<Bound>,
    /// Pairs three or more bonds apart, or in separate pieces, bounded below
    /// and never above their distance in the reference.
    apart: Vec<Lower>,
    flat: Option<Vec<Flat>>,
}

/// The least and most two atoms may be apart.
#[derive(Clone, Debug)]
struct Bound {
    atoms: (usize, usize),
    lower: f32,
    upper: f32,
}

/// The least two atoms may be apart, kept small since a ligand selection of
/// a few thousand atoms has millions of pairs.
#[derive(Clone, Debug)]
struct Lower {
    atoms: [u32; 2],
    lower: f32,
}

impl Template {
    /// The bounds of a ligand of heavy atoms with these `elements`,
    /// `positions` and bonded `neighbours`, and of the groups of atoms its
    /// definitions say are `flat`, those it holds flat.
    pub(crate) fn new(
        elements: &[FixedStr<2>],
        positions: &[Vec3],
        neighbours: &[Vec<usize>],
        flat: Option<Vec<Flat>>,
    ) -> Template {
        let length = |a: usize, b: usize| positions[a].distance(positions[b]);
        let angle = |a: usize, at: usize, b: usize| {
            (positions[a] - positions[at]).angle_between(positions[b] - positions[at])
        };
        let mut template = Template {
            bonds: Vec::new(),
            angles: Vec::new(),
            apart: Vec::new(),
            flat: flat.map(|groups| {
                groups
                    .into_iter()
                    .filter(|group| off_plane(&group.atoms, positions) <= FLATNESS)
                    .collect()
            }),
        };

        for a in 0..elements.len() {
            let hops = hops_from(a, neighbours);
            for b in a + 1..elements.len() {
                let lower = match hops[b] {
                    Some(1) => {
                        template.bonds.push(Bound {
                            atoms: (a, b),
                            lower: length(a, b) - BOND_MARGIN,
                            upper: length(a, b) + BOND_MARGIN,
                        });
                        continue;
                    }
                    Some(2) => {
                        let at = *neighbours[a]
                            .iter()
                            .find(|&&at| neighbours[at].contains(&b))
                            .expect("a pair two bonds apart shares a neighbour");
                        let (narrow, wide) = angle_range(angle(a, at, b));
                        let (first, second) = (length(a, at), length(at, b));
                        template.angles.push(Bound {
                            atoms: (a, b),
                            lower: across(first - BOND_MARGIN, second - BOND_MARGIN, narrow),
                            upper: across(first + BOND_MARGIN, second + BOND_MARGIN, wide),
                        });
                        continue;
                    }
                    Some(3) => {
                        let (next, last) = neighbours[a]
                            .iter()
                            .find_map(|&next| {
                                neighbours[next]
                                    .iter()
                                    .find(|&&last| {
                                        hops[last] == Some(2) && neighbours[last].contains(&b)
                                    })
                                    .map(|&last| (next, last))
                            })
                            .expect("a pair three bonds apart has a path");
                        let shortened = |x: usize, y: usize| length(x, y) - BOND_MARGIN;
                        cis(
                            [
                                shortened(a, next),
                                shortened(next, last),
                                shortened(last, b),
                            ],
                            [
                                angle_range(angle(a, next, last)).0,
                                angle_range(angle(next, last, b)).0,
                            ],
                        )
                    }
                    apart => {
                        let radii = vdw_radius(elements[a]) + vdw_radius(elements[b]);
                        let scale = match apart {
                            Some(4) => VDW_SCALE_15,
                            Some(5) => VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15),
                            _ => 1.0,
                        };
                        scale * radii
                    }
                };
                template.apart.push(Lower {
                    atoms: [a, b].map(|atom| atom as u32),
                    lower: lower.min(length(a, b)),
                });
            }
        }

        template
    }

    /// What a pose's ligand fails, its atoms at `positions` in the
    /// reference's order. Unless they were put there `by_name`, a flat group
    /// a symmetry of the bond graph could exchange for atoms of other
    /// chemistry is passed over, since the pairing may have done so.
    pub(crate) fn check(&self, positions: &[Vec3], by_name: bool) -> Geometry {
        let strays = |bounds: &[Bound]| {
            bounds
                .iter()
                .filter(|bound| bound.error(positions).abs() > TOLERANCE)
                .count()
        };

        Geometry {
            bonds: strays(&self.bonds),
            angles: strays(&self.angles),
            clashes: self
                .apart
                .iter()
                .filter(|bound| bound.error(positions) < -CLASH_TOLERANCE)
                .count(),
            unflat: self.flat.as_ref().map(|groups| {
                groups
                    .iter()
                    .filter(|group| by_name || !group.exchangeable)
                    .filter(|group| off_plane(&group.atoms, positions) > FLATNESS)
                    .count()
            }),
        }
    }
}

impl Lower {
    /// How far inside its bound the pair is at `positions`, as a fraction of
    /// the bound: negative inside, zero beyond.
    fn error(&self, positions: &[Vec3]) -> f32 {
        let [a, b] = self.atoms.map(|atom| atom as usize);
        let distance = positions[a].distance(positions[b]);

        ((distance - self.lower) / self.lower).min(0.0)
    }
}

impl Bound {
    /// How far outside its bounds the pair is at `positions`, as a fraction of
    /// the bound crossed: negative below, positive above, zero within.
    fn error(&self, positions: &[Vec3]) -> f32 {
        let (a, b) = self.atoms;
        let distance = positions[a].distance(positions[b]);
        if distance < self.lower {
            (distance - self.lower) / self.lower
        } else if distance > self.upper {
            (distance - self.upper) / self.upper
        } else {
            0.0
        }
    }
}

/// The groups of atoms that must lie flat in a ligand with these bonded
/// `neighbours`: each five- or six-membered ring of `aromatic` atoms, and each
/// of the `doubles` bonds with the atoms bonded to its ends.
fn flat_groups(
    neighbours: &[Vec<usize>],
    aromatic: &[bool],
    doubles: &[(usize, usize)],
) -> Vec<Vec<usize>> {
    let mut groups = rings(neighbours, aromatic);
    for &(a, b) in doubles {
        let mut group: Vec<usize> = [a, b]
            .into_iter()
            .chain(neighbours[a].iter().copied())
            .chain(neighbours[b].iter().copied())
            .collect();
        group.sort_unstable();
        group.dedup();
        if group.len() >= 4 {
            groups.push(group);
        }
    }

    groups
}

/// Each five- or six-membered cycle of `aromatic` atoms, once, from its
/// lowest atom.
fn rings(neighbours: &[Vec<usize>], aromatic: &[bool]) -> Vec<Vec<usize>> {
    fn walk(
        path: &mut Vec<usize>,
        neighbours: &[Vec<usize>],
        aromatic: &[bool],
        found: &mut Vec<Vec<usize>>,
    ) {
        let (start, last) = (path[0], path[path.len() - 1]);
        for &next in &neighbours[last] {
            // Each ring is walked both ways round from its lowest atom; the
            // way whose second atom is lower is kept.
            if next == start && path.len() >= 5 && path[1] < last {
                found.push(path.clone());
            } else if next > start && aromatic[next] && !path.contains(&next) && path.len() < 6 {
                path.push(next);
                walk(path, neighbours, aromatic, found);
                path.pop();
            }
        }
    }

    let mut found = Vec::new();
    for start in (0..neighbours.len()).filter(|&atom| aromatic[atom]) {
        walk(&mut vec![start], neighbours, aromatic, &mut found);
    }

    found
}

/// How many bonds each atom is from `start`, `None` for another piece.
fn hops_from(start: usize, neighbours: &[Vec<usize>]) -> Vec<Option<usize>> {
    let mut hops = vec![None; neighbours.len()];
    hops[start] = Some(0);
    let mut queue = VecDeque::from([start]);
    while let Some(atom) = queue.pop_front() {
        let next = hops[atom].map(|count| count + 1);
        for &bonded in &neighbours[atom] {
            if hops[bonded].is_none() {
                hops[bonded] = next;
                queue.push_back(bonded);
            }
        }
    }

    hops
}

/// The narrowest and widest a bond angle may be, as RDKit's `getAngleRange`
/// gives them.
fn angle_range(angle: f32) -> (f32, f32) {
    use std::f32::consts::PI;
    if angle + ANGLE_MARGIN >= PI {
        (PI - ANGLE_MARGIN, PI)
    } else if angle - ANGLE_MARGIN <= 0.0 {
        (0.0, ANGLE_MARGIN)
    } else {
        (angle - ANGLE_MARGIN, angle + ANGLE_MARGIN)
    }
}

/// The distance across a bond angle from its two bond lengths.
fn across(first: f32, second: f32, angle: f32) -> f32 {
    (first * first + second * second - 2.0 * first * second * angle.cos()).sqrt()
}

/// The distance between the ends of three bonds held cis, the least a torsion
/// lets them come, as RDKit's `compute14DistCis` gives it.
fn cis([first, second, third]: [f32; 3], [at_second, at_third]: [f32; 2]) -> f32 {
    let along = second - third * at_third.cos() - first * at_second.cos();
    let aside = third * at_third.sin() - first * at_second.sin();

    (along * along + aside * aside).sqrt()
}

/// The farthest any of the `group` atoms at `positions` lies from the plane
/// that fits them best. The plane is found in double precision: in single,
/// nalgebra can return eigenvectors that do not match their eigenvalues for a
/// scatter close to diagonal.
fn off_plane(group: &[usize], positions: &[Vec3]) -> f32 {
    let centre = group.iter().map(|&atom| positions[atom]).sum::<Vec3>() / group.len() as f32;
    let mut scatter = Matrix3::<f64>::zeros();
    for &atom in group {
        let offset = (positions[atom] - centre).as_dvec3();
        let column = nalgebra::Vector3::new(offset.x, offset.y, offset.z);
        scatter += column * column.transpose();
    }
    let eigen = SymmetricEigen::new(scatter);
    let least = eigen.eigenvalues.imin();
    let normal = eigen.eigenvectors.column(least);
    let normal = Vec3::new(normal[0] as f32, normal[1] as f32, normal[2] as f32);

    group
        .iter()
        .map(|&atom| (positions[atom] - centre).dot(normal).abs())
        .fold(0.0, f32::max)
}

/// The van der Waals radius RDKit gives an element (from the Blue Obelisk
/// data repository), in angstroms.
fn vdw_radius(element: FixedStr<2>) -> f32 {
    match element.as_str() {
        "C" => 1.7,
        "N" => 1.6,
        "O" => 1.55,
        "F" => 1.5,
        "P" => 1.95,
        "S" => 1.8,
        "CL" => 1.8,
        "BR" => 1.9,
        "I" => 2.1,
        "B" => 1.8,
        "SI" => 2.1,
        "SE" => 1.9,
        "AS" => 2.05,
        "NA" => 2.4,
        "MG" => 2.2,
        "K" => 2.8,
        "CA" => 2.4,
        "MN" | "FE" => 2.05,
        "CO" | "NI" | "CU" => 2.0,
        "ZN" => 2.1,
        _ => 2.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    impl Template {
        fn check_by_name(&self, positions: &[Vec3]) -> Geometry {
            self.check(positions, true)
        }
    }

    /// Groups no symmetry exchanges.
    fn fixed(groups: Vec<Vec<usize>>) -> Vec<Flat> {
        groups
            .into_iter()
            .map(|atoms| Flat {
                atoms,
                exchangeable: false,
            })
            .collect()
    }

    /// A chain of `atoms` carbons in a flat zigzag: bonds of 1.53 A, angles of
    /// 111.7 degrees, every torsion trans.
    fn zigzag(atoms: usize) -> (Vec<FixedStr<2>>, Vec<Vec3>, Vec<Vec<usize>>) {
        let positions = (0..atoms)
            .map(|at| Vec3::new(at as f32 * 1.27, (at % 2) as f32 * 0.86, 0.0))
            .collect();
        let neighbours = (0..atoms)
            .map(|at| {
                [at.checked_sub(1), Some(at + 1).filter(|&next| next < atoms)]
                    .into_iter()
                    .flatten()
                    .collect()
            })
            .collect();

        (vec![FixedStr::new("C"); atoms], positions, neighbours)
    }

    /// Benzene, or naphthalene when `fused`, flat in z = 0, every atom
    /// aromatic.
    fn rings(fused: bool) -> (Vec<FixedStr<2>>, Vec<Vec3>, Vec<Vec<usize>>) {
        let hexagon = |centre: Vec3| -> Vec<Vec3> {
            (0..6)
                .map(|at| {
                    let angle = (at as f32 + 0.5) * std::f32::consts::TAU / 6.0;
                    centre + Vec3::new(angle.cos(), angle.sin(), 0.0) * 1.39
                })
                .collect()
        };
        let mut positions = hexagon(Vec3::ZERO);
        let mut neighbours: Vec<Vec<usize>> =
            (0..6).map(|at| vec![(at + 5) % 6, (at + 1) % 6]).collect();
        if fused {
            // The second ring shares atoms 0 and 5, the hexagon's right edge.
            let second = hexagon(Vec3::new(1.39 * 3f32.sqrt(), 0.0, 0.0));
            positions.extend([second[1], second[2], second[3], second[4]]);
            neighbours[0].push(6);
            neighbours[5].push(9);
            neighbours.extend([vec![0, 7], vec![6, 8], vec![7, 9], vec![8, 5]]);
        }

        (
            vec![FixedStr::new("C"); positions.len()],
            positions,
            neighbours,
        )
    }

    #[test]
    fn a_pose_shaped_like_the_reference_passes() {
        let (elements, positions, neighbours) = zigzag(8);
        let template = Template::new(&elements, &positions, &neighbours, None);

        assert_eq!(
            template.check_by_name(&positions),
            Geometry {
                bonds: 0,
                angles: 0,
                clashes: 0,
                unflat: None,
            }
        );
    }

    #[test]
    fn a_bond_a_quarter_longer_than_the_reference_s_fails() {
        let (elements, positions, neighbours) = zigzag(4);
        let template = Template::new(&elements, &positions, &neighbours, None);
        let stretched = |by: f32| {
            let mut moved = positions.clone();
            moved[3] = positions[2] + (positions[3] - positions[2]) * by;
            template.check_by_name(&moved).bonds
        };

        assert_eq!(stretched(1.2), 0, "within the tolerance");
        assert_eq!(stretched(1.35), 1);
        assert_eq!(stretched(0.7), 1);
    }

    #[test]
    fn an_angle_bent_far_from_the_reference_s_fails() {
        let (elements, positions, neighbours) = zigzag(4);
        let template = Template::new(&elements, &positions, &neighbours, None);
        // Atom 3 swung about atom 2 until C1-C2-C3 closes to 70 degrees.
        let mut bent = positions.clone();
        let back = (positions[1] - positions[2]).normalize();
        let normal = Vec3::Z;
        let turned = glam::Quat::from_axis_angle(normal, 70f32.to_radians()) * back;
        bent[3] = positions[2] + turned * 1.534;

        let geometry = template.check_by_name(&bent);

        assert_eq!(geometry.bonds, 0);
        assert_eq!(geometry.angles, 1);
    }

    #[test]
    fn a_pose_folding_onto_itself_clashes() {
        let (elements, positions, neighbours) = zigzag(8);
        let template = Template::new(&elements, &positions, &neighbours, None);
        let mut folded = positions.clone();
        folded[7] = positions[0] + Vec3::new(0.0, 0.0, 1.2);

        assert!(template.check_by_name(&folded).clashes >= 1);
    }

    #[test]
    fn atoms_further_apart_in_the_graph_must_stay_further_apart() {
        let (elements, positions, neighbours) = zigzag(8);
        let template = Template::new(&elements, &positions, &neighbours, None);
        let lower = |a: usize, b: usize| {
            template
                .apart
                .iter()
                .find(|bound| bound.atoms == [a, b].map(|atom| atom as u32))
                .map(|bound| bound.lower)
        };

        assert_eq!(lower(0, 1), None, "bonded");
        assert_eq!(lower(0, 2), None, "an angle");
        let cis = lower(0, 3).expect("a 1-4 pair");
        assert!((2.4..2.6).contains(&cis), "the cis distance, not {cis}");
        let close = |found: Option<f32>, expected: f32| {
            found.is_some_and(|found| (found - expected).abs() < 1e-5)
        };
        assert!(close(lower(0, 4), 0.7 * 3.4));
        assert!(close(lower(0, 5), 0.85 * 3.4));
        assert!(close(lower(0, 6), 3.4));
    }

    #[test]
    fn aromatic_rings_are_found_one_ring_at_a_time() {
        let (_, _, benzene) = rings(false);
        let (_, _, naphthalene) = rings(true);

        let one = flat_groups(&benzene, &[true; 6], &[]);
        let two = flat_groups(&naphthalene, &[true; 10], &[]);

        assert_eq!(one.len(), 1);
        assert_eq!(two.len(), 2, "not the ten-membered rim");
        assert!(two.iter().all(|ring| ring.len() == 6));
        assert!(flat_groups(&benzene, &[false; 6], &[]).is_empty());
    }

    #[test]
    fn a_ring_lifted_out_of_its_plane_is_not_flat() {
        let (elements, positions, neighbours) = rings(true);
        let flat = flat_groups(&neighbours, &[true; 10], &[]);
        let template = Template::new(&elements, &positions, &neighbours, Some(fixed(flat)));
        let mut lifted = positions.clone();
        lifted[8].z += 0.6;

        assert_eq!(template.check_by_name(&positions).unflat, Some(0));
        assert_eq!(template.check_by_name(&lifted).unflat, Some(1));
    }

    #[test]
    fn a_double_bond_and_its_neighbours_lie_flat() {
        // Atoms 1 and 2 of a zigzag joined by a double bond, with an atom out
        // of plane on either end.
        let (elements, positions, neighbours) = zigzag(4);
        let groups = flat_groups(&neighbours, &[false; 4], &[(1, 2)]);
        let template = Template::new(
            &elements,
            &positions,
            &neighbours,
            Some(fixed(groups.clone())),
        );
        let mut twisted = positions.clone();
        twisted[3].z += 1.0;

        assert_eq!(groups, [vec![0, 1, 2, 3]]);
        assert_eq!(template.check_by_name(&positions).unflat, Some(0));
        assert_eq!(template.check_by_name(&twisted).unflat, Some(1));
    }

    // An iron-sulfur cluster or a heme's iron is bonded to nothing in a graph
    // read from distances, and a pose is not faulted for what the reference
    // does.
    #[test]
    fn atoms_as_close_as_in_the_reference_do_not_clash() {
        let (elements, mut positions, neighbours) = zigzag(8);
        positions[7] = positions[0] + Vec3::new(0.0, 0.0, 2.0);
        let template = Template::new(&elements, &positions, &neighbours, None);
        let mut closer = positions.clone();
        closer[7] = positions[0] + Vec3::new(0.0, 0.0, 1.3);

        assert_eq!(template.check_by_name(&positions).clashes, 0);
        assert_eq!(template.check_by_name(&closer).clashes, 1);
    }

    // HEC defines its vinyl carbons doubly bonded, though in cytochrome c each
    // is bonded to a cysteine and tetrahedral.
    #[test]
    fn a_group_the_reference_holds_out_of_plane_need_not_be_flat() {
        let (elements, mut positions, neighbours) = zigzag(4);
        positions[3].z += 1.0;
        let groups = flat_groups(&neighbours, &[false; 4], &[(1, 2)]);
        let template = Template::new(&elements, &positions, &neighbours, Some(fixed(groups)));

        assert_eq!(template.check_by_name(&positions).unflat, Some(0));
    }

    #[test]
    fn a_flat_group_lying_along_the_axes_is_flat() {
        let groups = [
            [
                (-1.42, -1.299, 0.0),
                (-0.67, 0.0, 0.0),
                (0.67, 0.0, 0.0),
                (1.42, -1.299, 0.0),
            ],
            [
                (-1.42, 1.299, 0.0),
                (-0.67, 0.0, 0.0),
                (0.67, 0.0, 0.0),
                (1.42, -1.299, 0.0),
            ],
            [
                (0.0, -1.299, -1.42),
                (0.0, 0.0, -0.67),
                (0.0, 0.0, 0.67),
                (0.0, -1.299, 1.42),
            ],
        ];
        for group in groups {
            let positions = group.map(|(x, y, z)| Vec3::new(x, y, z));

            let off = off_plane(&[0, 1, 2, 3], &positions);

            assert!(off < 0.01, "{off} A off: {positions:?}");
        }
    }
}
