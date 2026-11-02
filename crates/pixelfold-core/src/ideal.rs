//! The distance bounds RDKit's distance geometry gives a molecule from its
//! bonds alone, which PoseBusters holds a pose to when nothing else says what
//! the molecule should look like (Buttenschoen, Morris and Deane, Chem. Sci.
//! 15:3130, 2024).
//!
//! This follows RDKit 2024.09's `BoundsMatrixBuilder.cpp` under PoseBusters'
//! settings (bounds five bonds out, scaled van der Waals floors, triangle
//! smoothing): a bond's length from UFF's radii and electronegativities, a
//! bond angle from the hybridisation of its middle atom and the rings it lies
//! in, a torsion's two ends cis or trans where the rings, double bonds and
//! amides about it fix them, five atoms in a row from the two torsions, and
//! everything further apart kept out by van der Waals radii; then the
//! triangle inequality tightens every pair through every third atom. A double
//! bond's E or Z is read from the pose itself, as RDKit reads a 3D molecule's.

use std::collections::{HashMap, HashSet, VecDeque};
use std::f64::consts::PI;

use glam::Vec3;

use crate::components::BondOrder;
use crate::geometry::{Template, puckered_rings, vdw_radius};
use crate::molecule::{Hybrid, Molecule, key};
use crate::stereo;

/// RDKit's margin around a bond length, in angstroms.
const DIST12_DELTA: f64 = 0.01;

/// RDKit's margin around the distance across a bond angle.
const DIST13_TOL: f64 = 0.04;

/// RDKit's margin around a distance a torsion fixes.
const GEN_DIST_TOL: f64 = 0.06;

/// RDKit's margin around a distance two fixed torsions fix.
const DIST15_TOL: f64 = 0.08;

/// The share of their van der Waals radii summed that atoms four bonds apart
/// may come within; five apart, halfway from it to the whole sum.
const VDW_SCALE_15: f64 = 0.7;

/// RDKit's upper bound for a pair nothing bounds.
const MAX_UPPER: f64 = 1000.0;

/// UFF's scaling of the Pauling bond-order correction.
const LAMBDA: f64 = 0.1332;

/// The most atoms whose bounds are triangle smoothed: the pass is cubic, and a
/// selection of whole chains as a ligand would take minutes.
const MAX_SMOOTHED: usize = 512;

/// How a torsion's two ends sit, as RDKit records each path of four atoms for
/// the paths of five.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Path {
    Cis,
    Trans,
    Other,
}

/// A molecule's heavy-atom pair bounds, lower below the diagonal and upper
/// above it, as RDKit keeps them.
struct Matrix {
    size: usize,
    values: Vec<f64>,
}

impl Matrix {
    fn new(size: usize) -> Matrix {
        let mut values = vec![0.0; size * size];
        for i in 0..size {
            for j in i + 1..size {
                values[i * size + j] = MAX_UPPER;
            }
        }
        Matrix { size, values }
    }

    fn lower(&self, a: usize, b: usize) -> f64 {
        self.values[a.max(b) * self.size + a.min(b)]
    }

    fn upper(&self, a: usize, b: usize) -> f64 {
        self.values[a.min(b) * self.size + a.max(b)]
    }

    fn set_lower(&mut self, a: usize, b: usize, value: f64) {
        self.values[a.max(b) * self.size + a.min(b)] = value;
    }

    fn set_upper(&mut self, a: usize, b: usize, value: f64) {
        self.values[a.min(b) * self.size + a.max(b)] = value;
    }

    /// RDKit's `_checkAndSetBounds`: an unset bound takes the new one, a set
    /// one only loosens.
    fn widen(&mut self, a: usize, b: usize, lower: f64, upper: f64) {
        let (current_lower, current_upper) = (self.lower(a, b), self.upper(a, b));
        if current_lower <= DIST12_DELTA || (lower < current_lower && lower > DIST12_DELTA) {
            self.set_lower(a, b, lower);
        }
        if current_upper >= MAX_UPPER || (upper > current_upper && upper < MAX_UPPER) {
            self.set_upper(a, b, upper);
        }
    }

    /// RDKit's `triangleSmoothBounds` with no tolerance: every pair's bounds
    /// tightened through every third atom, stopping where they cross.
    fn smooth(&mut self) {
        let n = self.size;
        let v = &mut self.values;
        for k in 0..n {
            for i in 0..n.saturating_sub(1) {
                if i == k {
                    continue;
                }
                let (ii, ik) = (i.min(k), i.max(k));
                let upper_ik = v[ii * n + ik];
                let lower_ik = v[ik * n + ii];
                for j in i + 1..n {
                    if j == k {
                        continue;
                    }
                    let (jj, jk) = (j.min(k), j.max(k));
                    let upper_kj = v[jj * n + jk];
                    let sum = upper_ik + upper_kj;
                    if v[i * n + j] > sum {
                        v[i * n + j] = sum;
                    }
                    let through_i = lower_ik - upper_kj;
                    let through_j = v[jk * n + jj] - upper_ik;
                    if v[j * n + i] < through_i {
                        v[j * n + i] = through_i;
                    } else if v[j * n + i] < through_j {
                        v[j * n + i] = through_j;
                    }
                    if v[j * n + i] - v[i * n + j] > 0.0 {
                        return;
                    }
                }
            }
        }
    }
}

/// UFF's bond radius and electronegativity of an atom type, by RDKit's label.
fn uff(label: &str) -> Option<(f64, f64)> {
    Some(match label {
        "H_" => (0.354, 4.528),
        "B_3" => (0.838, 5.11),
        "B_2" => (0.828, 5.11),
        "C_3" => (0.757, 5.343),
        "C_R" => (0.729, 5.343),
        "C_2" => (0.732, 5.343),
        "C_1" => (0.706, 5.343),
        "N_3" => (0.7, 6.899),
        "N_R" => (0.699, 6.899),
        "N_2" => (0.685, 6.899),
        "N_1" => (0.656, 6.899),
        "O_3" => (0.658, 8.741),
        "O_R" => (0.68, 8.741),
        "O_2" => (0.634, 8.741),
        "O_1" => (0.639, 8.741),
        "F_" => (0.668, 10.874),
        "Si3" => (1.117, 4.168),
        "P_3+3" => (1.101, 5.463),
        "P_3+5" => (1.056, 5.463),
        "S_3+2" => (1.064, 6.928),
        "S_3+4" => (1.049, 6.928),
        "S_3+6" => (1.027, 6.928),
        "S_R" => (1.077, 6.928),
        "S_2" => (0.854, 6.928),
        "Cl" => (1.044, 8.564),
        "As3+3" => (1.211, 5.188),
        "Se3+2" => (1.19, 6.428),
        "Br" => (1.192, 7.79),
        "I_" => (1.382, 6.822),
        _ => return None,
    })
}

/// RDKit's UFF atom label: the symbol, its hybridisation (`R` for an
/// aromatic or conjugated sp2 C, N, O or S), and for P and S the valence.
fn uff_label(molecule: &Molecule, atom: usize) -> String {
    let element = molecule.element(atom);
    let symbol = molecule.elements[atom].as_str();
    let mut label = match symbol.len() {
        1 => format!("{symbol}_"),
        _ => {
            let mut chars = symbol.chars();
            let first = chars.next().unwrap_or(' ');
            format!("{first}{}", chars.as_str().to_ascii_lowercase())
        }
    };
    let hybrid = molecule.hybrids[atom];
    if element.default_valence().is_none() || (element.outer != 1 && element.outer != 7) {
        match element.number {
            12..=15 => label.push('3'),
            _ => match hybrid {
                Hybrid::Sp => label.push('1'),
                Hybrid::Sp2 => {
                    let resonant = molecule.aromatic[atom] || molecule.has_conjugated_bond(atom);
                    let organic = matches!(element.number, 6 | 7 | 8 | 16);
                    label.push(if resonant && organic { 'R' } else { '2' });
                }
                Hybrid::Sp3 => label.push('3'),
                Hybrid::Sp3d => label.push('5'),
                Hybrid::Sp3d2 => label.push('6'),
                Hybrid::S | Hybrid::Unspecified => {}
            },
        }
    }
    let valence = molecule.valence(atom);
    match (element.number, valence) {
        (15, 3) => label.push_str("+3"),
        (15, 5) => label.push_str("+5"),
        (16, 2 | 4 | 6) if hybrid != Hybrid::Sp2 => label.push_str(&format!("+{valence}")),
        (33, 3) => label.push_str("+3"),
        (34, 2) => label.push_str("+2"),
        (30, 2) => label.push_str("+2"),
        _ => {}
    }

    label
}

/// UFF's rest length of a bond of `order` between atoms of these radii and
/// electronegativities.
fn rest_length(order: f64, (ri, xi): (f64, f64), (rj, xj): (f64, f64)) -> f64 {
    let pauling = -LAMBDA * (ri + rj) * order.ln();
    let electronegativity = ri * rj * (xi.sqrt() - xj.sqrt()).powi(2) / (xi * ri + xj * rj);

    ri + rj + pauling - electronegativity
}

fn across(first: f64, second: f64, angle: f64) -> f64 {
    (first * first + second * second - 2.0 * first * second * angle.cos()).sqrt()
}

fn cis(first: f64, second: f64, third: f64, at_second: f64, at_third: f64) -> f64 {
    let dx = second - third * at_third.cos() - first * at_second.cos();
    let dy = third * at_third.sin() - first * at_second.sin();
    (dx * dx + dy * dy).sqrt()
}

fn trans(first: f64, second: f64, third: f64, at_second: f64, at_third: f64) -> f64 {
    let dx = second - third * at_third.cos() - first * at_second.cos();
    let dy = third * at_third.sin() + first * at_second.sin();
    (dx * dx + dy * dy).sqrt()
}

/// The distance between the ends of three bonds held at `torsion`.
fn at_torsion(
    first: f64,
    second: f64,
    third: f64,
    at_second: f64,
    at_third: f64,
    torsion: f64,
) -> f64 {
    let p1 = (first * at_second.cos(), first * at_second.sin(), 0.0);
    let (x, y) = (second - third * at_third.cos(), third * at_third.sin());
    let p4 = (x, y * torsion.cos(), y * torsion.sin());
    ((p4.0 - p1.0).powi(2) + (p4.1 - p1.1).powi(2) + (p4.2 - p1.2).powi(2)).sqrt()
}

/// RDKit's distances between the ends of four bonds, the first three held
/// `first` and the last three `second` (cis or trans), from the bond lengths
/// `d` and angles `angles`.
fn five(d: [f64; 4], angles: [f64; 3], first: Path, second: Path) -> f64 {
    let [d1, d2, d3, d4] = d;
    let [a12, a23, a34] = angles;
    let (dy14, cval) = match first {
        Path::Cis => (
            d3 * a23.sin() - d1 * a12.sin(),
            d3 - d2 * a23.cos() + d1 * (a12 + a23).cos(),
        ),
        _ => (
            d3 * a23.sin() + d1 * a12.sin(),
            d3 - d2 * a23.cos() + d1 * (a12 - a23).cos(),
        ),
    };
    let dx14 = d2 - d3 * a23.cos() - d1 * a12.cos();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let a143 = (cval / d14).clamp(-1.0, 1.0).acos();
    let a145 = match (first, second) {
        (Path::Cis, Path::Cis) | (Path::Trans, Path::Cis) => a34 - a143,
        _ => a34 + a143,
    };

    across(d14, d4, a145)
}

/// The geometry checks a pose of `molecule` at `positions` is held to, as
/// PoseBusters makes them with no reference.
pub(crate) fn template(molecule: &Molecule, positions: &[Vec3]) -> Template {
    let mut builder = Builder::new(molecule, positions);
    let matrix = builder.bounds();
    builder.template(&matrix)
}

/// What the bounds are built from, and what they record along the way.
struct Builder<'a> {
    molecule: &'a Molecule,
    positions: &'a [Vec3],
    bonds: Vec<(usize, usize)>,
    index: HashMap<(usize, usize), usize>,
    /// Topological distances.
    hops: Vec<Vec<usize>>,
    lengths: Vec<f64>,
    angles: HashMap<(usize, usize), f64>,
}

impl<'a> Builder<'a> {
    fn new(molecule: &'a Molecule, positions: &'a [Vec3]) -> Builder<'a> {
        let mut bonds: Vec<(usize, usize)> = molecule
            .neighbours
            .iter()
            .enumerate()
            .flat_map(|(a, bonded)| bonded.iter().filter(move |&&b| b > a).map(move |&b| (a, b)))
            .collect();
        bonds.sort_unstable();
        let index = bonds
            .iter()
            .enumerate()
            .map(|(at, &bond)| (bond, at))
            .collect();
        let hops = (0..molecule.len())
            .map(|atom| hops_from(atom, &molecule.neighbours))
            .collect();

        Builder {
            molecule,
            positions,
            bonds,
            index,
            hops,
            lengths: Vec::new(),
            angles: HashMap::new(),
        }
    }

    fn bond(&self, a: usize, b: usize) -> usize {
        self.index[&key(a, b)]
    }

    /// The atom two bonds share.
    fn shared(&self, first: usize, second: usize) -> usize {
        let (a, b) = self.bonds[first];
        let (c, d) = self.bonds[second];
        if a == c || a == d { a } else { b }
    }

    fn other(&self, bond: usize, atom: usize) -> usize {
        let (a, b) = self.bonds[bond];
        if a == atom { b } else { a }
    }

    fn incident(&self, atom: usize) -> impl Iterator<Item = usize> + '_ {
        self.molecule.neighbours[atom]
            .iter()
            .map(move |&other| self.bond(atom, other))
    }

    fn angle(&self, first: usize, second: usize) -> f64 {
        self.angles.get(&(first, second)).copied().unwrap_or(0.0)
    }

    fn bounds(&mut self) -> Matrix {
        let mut matrix = Matrix::new(self.molecule.len());
        self.set12(&mut matrix);
        self.set13(&mut matrix);
        let paths = self.set14(&mut matrix);
        self.set15(&mut matrix, &paths);
        self.set_vdw(&mut matrix);
        if matrix.size <= MAX_SMOOTHED {
            matrix.smooth();
        }

        matrix
    }

    fn set12(&mut self, matrix: &mut Matrix) {
        let molecule = self.molecule;
        let params: Vec<Option<(f64, f64)>> = (0..molecule.len())
            .map(|atom| uff(&uff_label(molecule, atom)))
            .collect();
        let mut squished = vec![false; molecule.len()];
        for &(a, b) in &self.bonds {
            let heavy = molecule.element(a).number > 10 || molecule.element(b).number > 10;
            if molecule.is_conjugated(a, b) && heavy && molecule.rings.bond_in_ring_of(a, b, 5) {
                squished[a] = true;
                squished[b] = true;
            }
        }

        self.lengths = self
            .bonds
            .iter()
            .map(|&(a, b)| {
                let order = match molecule.order(a, b) {
                    Some(BondOrder::Double) => 2.0,
                    Some(BondOrder::Triple) => 3.0,
                    Some(BondOrder::Aromatic) => 1.5,
                    _ => 1.0,
                };
                if let (Some(first), Some(second)) = (params[a], params[b]) {
                    let length = rest_length(order, first, second);
                    let squish = if squished[a] || squished[b] { 0.2 } else { 0.0 };
                    matrix.set_upper(a, b, length + squish + DIST12_DELTA);
                    matrix.set_lower(a, b, length - squish - DIST12_DELTA);
                    length
                } else {
                    let length = f64::from(
                        vdw_radius(molecule.elements[a]) + vdw_radius(molecule.elements[b]),
                    ) / 2.0;
                    matrix.set_upper(a, b, 1.5 * length);
                    matrix.set_lower(a, b, 0.5 * length);
                    length
                }
            })
            .collect();
    }

    fn set13_pair(&self, matrix: &mut Matrix, first: usize, second: usize, angle: f64) {
        let molecule = self.molecule;
        let centre = self.shared(first, second);
        let (a, b) = (self.other(first, centre), self.other(second, centre));
        let distance = across(self.lengths[first], self.lengths[second], angle);
        let larger_sp2 = |atom: usize| {
            molecule.element(atom).number > 13
                && molecule.hybrids[atom] == Hybrid::Sp2
                && !molecule.rings.of_atom(atom).is_empty()
        };
        let tolerance = [a, centre, b]
            .into_iter()
            .filter(|&atom| larger_sp2(atom))
            .fold(DIST13_TOL, |tolerance, _| tolerance * 2.0);
        matrix.widen(a, b, distance - tolerance, distance + tolerance);
    }

    fn set13(&mut self, matrix: &mut Matrix) {
        let molecule = self.molecule;
        let n = molecule.len();
        let mut visited = vec![0usize; n];
        let mut taken = vec![0.0f64; n];
        let mut done: HashSet<(usize, usize)> = HashSet::new();

        for ring in molecule.rings.all() {
            let size = ring.len();
            let mut previous = ring[size - 1];
            for at in 0..size {
                let (centre, next) = (ring[at], ring[(at + 1) % size]);
                let first = self.bond(previous, centre);
                let second = self.bond(centre, next);
                if !done.contains(&(first, second)) && !done.contains(&(second, first)) {
                    let angle = ring_angle(molecule.hybrids[centre], size);
                    self.set13_pair(matrix, first, second, angle);
                    self.angles.insert((first, second), angle);
                    self.angles.insert((second, first), angle);
                    visited[centre] += 1;
                    taken[centre] += angle;
                    done.insert((first, second));
                    done.insert((second, first));
                }
                previous = centre;
            }
        }

        for centre in 0..n {
            let degree = molecule.total_degree(centre);
            let heavy = molecule.neighbours[centre].len();
            let pairs = heavy * heavy.saturating_sub(1) / 2;
            if pairs == visited[centre] {
                continue;
            }
            let hybrid = molecule.hybrids[centre];
            let incident: Vec<usize> = self.incident(centre).collect();
            let in_ring = visited[centre] >= 1;
            for (at, &first) in incident.iter().enumerate() {
                for &second in &incident[..at] {
                    if self.angles.contains_key(&(first, second)) {
                        continue;
                    }
                    let angle = if in_ring {
                        match hybrid {
                            Hybrid::Sp2 => {
                                (2.0 * PI - taken[centre]) / (pairs - visited[centre]) as f64
                            }
                            Hybrid::Sp3 if molecule.rings.atom_in_ring_of(centre, 3) => {
                                116.0_f64.to_radians()
                            }
                            Hybrid::Sp3 if molecule.rings.atom_in_ring_of(centre, 4) => {
                                112.0_f64.to_radians()
                            }
                            Hybrid::Sp3 => 109.5_f64.to_radians(),
                            _ if degree == 5 => 105.0_f64.to_radians(),
                            _ if degree == 6 => 135.0_f64.to_radians(),
                            _ => 120.0_f64.to_radians(),
                        }
                    } else {
                        match hybrid {
                            Hybrid::Sp => PI,
                            Hybrid::Sp2 => 2.0 * PI / 3.0,
                            Hybrid::Sp3 => 109.5_f64.to_radians(),
                            Hybrid::Sp3d => 105.0_f64.to_radians(),
                            Hybrid::Sp3d2 => 135.0_f64.to_radians(),
                            _ => 120.0_f64.to_radians(),
                        }
                    };
                    if in_ring || molecule.total_degree(centre) <= 4 {
                        self.set13_pair(matrix, first, second, angle);
                    } else {
                        let (a, b) = (self.other(first, centre), self.other(second, centre));
                        let most = self.lengths[first] + self.lengths[second];
                        matrix.widen(a, b, 1.0, most * 1.2);
                    }
                    self.angles.insert((first, second), angle);
                    self.angles.insert((second, first), angle);
                    taken[centre] += angle;
                    visited[centre] += 1;
                }
            }
        }
    }

    /// Whether the pose holds the ends of the path `a`-`b`=`c`-`d` cis, as
    /// RDKit reads a stereogenic double bond's E or Z from 3D, or `None` when
    /// the bond is not one.
    fn stereo(&self, a: usize, b: usize, c: usize, d: usize, classes: &[usize]) -> Option<bool> {
        let molecule = self.molecule;
        if !matches!(molecule.order(b, c), Some(BondOrder::Double))
            || molecule
                .rings
                .of_bond(b, c)
                .any(|ring| molecule.rings.all()[ring].len() < 8)
        {
            return None;
        }
        let stereogenic = |end: usize, partner: usize| {
            let others: Vec<usize> = molecule.neighbours[end]
                .iter()
                .copied()
                .filter(|&other| other != partner)
                .collect();
            match (others.len(), molecule.hydrogens(end)) {
                (1, _) => true,
                (2, _) => classes[others[0]] != classes[others[1]],
                _ => false,
            }
        };
        if !stereogenic(b, c) || !stereogenic(c, b) {
            return None;
        }
        let [pa, pb, pc, pd] = [a, b, c, d].map(|atom| self.positions[atom]);
        let (first, axis, last) = (pa - pb, pc - pb, pd - pc);
        let normal_a = axis.cross(first);
        let normal_d = axis.cross(last);
        if normal_a.length_squared() < 1e-6 || normal_d.length_squared() < 1e-6 {
            return None;
        }

        Some(normal_a.dot(normal_d) > 0.0)
    }

    /// RDKit's `set14Bounds`, returning every path of four atoms it records
    /// with how its ends sit.
    fn set14(&self, matrix: &mut Matrix) -> Vec<([usize; 3], Path)> {
        let molecule = self.molecule;
        let classes = stereo::classes(&molecule.elements, &molecule.neighbours);
        let mut paths: Vec<([usize; 3], Path)> = Vec::new();
        let mut ring_pairs: HashSet<(usize, usize)> = HashSet::new();
        let mut done: HashSet<[usize; 3]> = HashSet::new();

        for ring in molecule.rings.all() {
            let size = ring.len();
            let bonds: Vec<usize> = (0..size)
                .map(|at| self.bond(ring[at], ring[(at + 1) % size]))
                .collect();
            let mut first = bonds[size - 1];
            for at in 0..size {
                let (second, third) = (bonds[at], bonds[(at + 1) % size]);
                ring_pairs.insert((first, second));
                ring_pairs.insert((second, first));
                done.insert([first, second, third]);
                done.insert([third, second, first]);
                if size > 5 {
                    if let Some(path) =
                        self.in_ring14(matrix, [first, second, third], size, &classes)
                    {
                        paths.push(path);
                    }
                } else {
                    let (b, c) = (self.shared(first, second), self.shared(second, third));
                    let both_sp2 =
                        molecule.hybrids[b] == Hybrid::Sp2 && molecule.hybrids[c] == Hybrid::Sp2;
                    paths.push((
                        [first, second, third],
                        if both_sp2 { Path::Cis } else { Path::Other },
                    ));
                }
                first = second;
            }
        }

        for (second, &(b, c)) in self.bonds.iter().enumerate() {
            for first in self.incident(b).filter(|&bond| bond != second) {
                for third in self.incident(c).filter(|&bond| bond != second) {
                    if done.contains(&[first, second, third])
                        || done.contains(&[third, second, first])
                    {
                        continue;
                    }
                    let triple = [first, second, third];
                    let in_ring = |bond: usize| {
                        let (x, y) = self.bonds[bond];
                        molecule.rings.of_bond(x, y).next().is_some()
                    };
                    let path = if ring_pairs.contains(&(first, second))
                        || ring_pairs.contains(&(second, third))
                    {
                        self.two_in_same_ring14(matrix, triple)
                    } else if (in_ring(first) && in_ring(second))
                        || (in_ring(second) && in_ring(third))
                        || in_ring(second)
                    {
                        self.in_ring14(matrix, triple, 0, &classes)
                    } else {
                        self.chain14(matrix, triple, &classes)
                    };
                    paths.extend(path);
                }
            }
        }

        paths
    }

    /// The path's atoms, bond lengths and angles, or `None` when its ends are
    /// fewer than three bonds apart.
    fn path14(
        &self,
        [first, second, third]: [usize; 3],
    ) -> Option<([usize; 4], [f64; 3], [f64; 2])> {
        let b = self.shared(first, second);
        let c = self.shared(second, third);
        let a = self.other(first, b);
        let d = self.other(third, c);
        if a == d || self.hops[a][d] < 3 {
            return None;
        }

        Some((
            [a, b, c, d],
            [
                self.lengths[first],
                self.lengths[second],
                self.lengths[third],
            ],
            [self.angle(first, second), self.angle(second, third)],
        ))
    }

    fn set_path(
        &self,
        matrix: &mut Matrix,
        a: usize,
        d: usize,
        lengths: [f64; 3],
        angles: [f64; 2],
        path: Path,
    ) -> Path {
        let [l1, l2, l3] = lengths;
        let [a12, a23] = angles;
        let (lower, upper) = match path {
            Path::Cis => {
                let lower = cis(l1, l2, l3, a12, a23) - GEN_DIST_TOL;
                (lower, lower + 2.0 * GEN_DIST_TOL)
            }
            Path::Trans => {
                let lower = trans(l1, l2, l3, a12, a23) - GEN_DIST_TOL;
                (lower, lower + 2.0 * GEN_DIST_TOL)
            }
            Path::Other => {
                let (mut lower, mut upper) =
                    (cis(l1, l2, l3, a12, a23), trans(l1, l2, l3, a12, a23));
                if upper < lower {
                    std::mem::swap(&mut lower, &mut upper);
                }
                if (upper - lower).abs() < DIST12_DELTA {
                    (lower - GEN_DIST_TOL, upper + GEN_DIST_TOL)
                } else {
                    (lower, upper)
                }
            }
        };
        matrix.widen(a, d, lower, upper);

        path
    }

    /// RDKit's `_setInRing14Bounds`; `size` is the ring's, or zero when the
    /// bonds lie in no one ring.
    fn in_ring14(
        &self,
        matrix: &mut Matrix,
        triple: [usize; 3],
        size: usize,
        classes: &[usize],
    ) -> Option<([usize; 3], Path)> {
        let molecule = self.molecule;
        let ([a, b, c, d], lengths, angles) = self.path14(triple)?;
        let stereo = self.stereo(a, b, c, d, classes);
        let both_sp2 = molecule.hybrids[b] == Hybrid::Sp2 && molecule.hybrids[c] == Hybrid::Sp2;
        let [first, _, third] = triple;
        let rings_of = |bond: usize| {
            let (x, y) = self.bonds[bond];
            molecule.rings.of_bond(x, y).collect::<Vec<usize>>()
        };
        let path = if size <= 8 && both_sp2 && stereo != Some(false) {
            if molecule.rings.of_bond(b, c).count() > 1 {
                let (around_first, around_third) = (rings_of(first), rings_of(third));
                if around_first.len() == 1
                    && around_third.len() == 1
                    && around_first[0] == around_third[0]
                {
                    Path::Cis
                } else {
                    Path::Other
                }
            } else {
                Path::Cis
            }
        } else {
            match stereo {
                Some(true) => Path::Cis,
                Some(false) => Path::Trans,
                None => Path::Other,
            }
        };

        Some((triple, self.set_path(matrix, a, d, lengths, angles, path)))
    }

    /// RDKit's `_setTwoInSameRing14Bounds`.
    fn two_in_same_ring14(
        &self,
        matrix: &mut Matrix,
        triple: [usize; 3],
    ) -> Option<([usize; 3], Path)> {
        let molecule = self.molecule;
        let ([a, b, c, d], lengths, angles) = self.path14(triple)?;
        if molecule.order(a, c).is_some() || molecule.order(d, b).is_some() {
            return None;
        }
        let both_sp2 = molecule.hybrids[b] == Hybrid::Sp2 && molecule.hybrids[c] == Hybrid::Sp2;
        let path = if both_sp2 { Path::Trans } else { Path::Other };

        Some((triple, self.set_path(matrix, a, d, lengths, angles, path)))
    }

    /// RDKit's `_setChain14Bounds`, amides and esters held trans.
    fn chain14(
        &self,
        matrix: &mut Matrix,
        triple: [usize; 3],
        classes: &[usize],
    ) -> Option<([usize; 3], Path)> {
        let molecule = self.molecule;
        let ([a, b, c, d], lengths, angles) = self.path14(triple)?;
        let order = |x: usize, y: usize| molecule.order(x, y).unwrap_or(BondOrder::Single);
        let path = match order(b, c) {
            BondOrder::Double => {
                if order(a, b) == BondOrder::Double || order(c, d) == BondOrder::Double {
                    Path::Cis
                } else {
                    match self.stereo(a, b, c, d, classes) {
                        Some(true) => Path::Cis,
                        Some(false) => Path::Trans,
                        None => Path::Other,
                    }
                }
            }
            BondOrder::Single => {
                let number = |atom: usize| molecule.element(atom).number;
                if number(b) == 16 && number(c) == 16 {
                    let [l1, l2, l3] = lengths;
                    let lower =
                        at_torsion(l1, l2, l3, angles[0], angles[1], PI / 2.0) - GEN_DIST_TOL;
                    matrix.widen(a, d, lower, lower + 2.0 * GEN_DIST_TOL);
                    return Some((triple, Path::Other));
                } else if self.amide_ester14(a, b, c, d) || self.amide_ester14(d, c, b, a) {
                    Path::Cis
                } else if self.amide_ester15(a, b, c) || self.amide_ester15(d, c, b) {
                    Path::Trans
                } else {
                    Path::Other
                }
            }
            _ => Path::Other,
        };

        Some((triple, self.set_path(matrix, a, d, lengths, angles, path)))
    }

    /// RDKit's `_checkAmideEster14`: `a`-`b` single, `c`=`d` with `c` a
    /// carbon and `d` an oxygen or nitrogen, `b` an oxygen or an NH.
    fn amide_ester14(&self, a: usize, b: usize, c: usize, d: usize) -> bool {
        let molecule = self.molecule;
        let number = |atom: usize| molecule.element(atom).number;
        number(c) == 6
            && molecule.order(c, d) == Some(BondOrder::Double)
            && matches!(number(d), 7 | 8)
            && molecule.order(a, b) == Some(BondOrder::Single)
            && (number(b) == 8 || (number(b) == 7 && molecule.hydrogens(b) == 1))
    }

    /// RDKit's `_checkAmideEster15`: `a`-`b` single with `b` an oxygen or an
    /// NH, `b`-`c` single with `c` a carbonyl carbon.
    fn amide_ester15(&self, a: usize, b: usize, c: usize) -> bool {
        let molecule = self.molecule;
        let number = |atom: usize| molecule.element(atom).number;
        let carbonyl = number(c) == 6
            && molecule.neighbours[c].len() > 2
            && molecule.neighbours[c].iter().any(|&other| {
                matches!(number(other), 7 | 8)
                    && molecule.order(c, other) == Some(BondOrder::Double)
            });
        (number(b) == 8 || (number(b) == 7 && molecule.hydrogens(b) == 1))
            && molecule.order(a, b) == Some(BondOrder::Single)
            && molecule.order(b, c) == Some(BondOrder::Single)
            && carbonyl
    }

    /// RDKit's `set15Bounds`.
    fn set15(&self, matrix: &mut Matrix, paths: &[([usize; 3], Path)]) {
        let molecule = self.molecule;
        let n = molecule.len();
        let mut kinds: HashMap<[usize; 3], Path> = HashMap::new();
        for &([first, second, third], path) in paths {
            if path != Path::Other {
                kinds.insert([first, second, third], path);
                kinds.insert([third, second, first], path);
            }
        }
        let mut set = vec![false; n * n];
        for &([first, second, third], path) in paths {
            for (first, third) in [(first, third), (third, first)] {
                let b = self.shared(first, second);
                let a = self.other(first, b);
                let c = self.shared(second, third);
                let d = self.other(third, c);
                for fourth in self.incident(d).filter(|&bond| bond != third) {
                    let e = self.other(fourth, d);
                    if a == e || self.hops[a][e] < 4 {
                        continue;
                    }
                    if !(matrix.lower(a, e) < DIST12_DELTA || set[a * n + e]) {
                        continue;
                    }
                    let lengths = [
                        self.lengths[first],
                        self.lengths[second],
                        self.lengths[third],
                        self.lengths[fourth],
                    ];
                    let angles = [
                        self.angle(first, second),
                        self.angle(second, third),
                        self.angle(third, fourth),
                    ];
                    let next = kinds
                        .get(&[second, third, fourth])
                        .copied()
                        .unwrap_or(Path::Other);
                    let (lower, upper) = match (path, next) {
                        (Path::Other, Path::Other) => {
                            let radii = f64::from(
                                vdw_radius(molecule.elements[a]) + vdw_radius(molecule.elements[e]),
                            );
                            (VDW_SCALE_15 * radii, MAX_UPPER)
                        }
                        (Path::Other, fixed) => {
                            let [d1, d2, d3, d4] = lengths;
                            let [a12, a23, a34] = angles;
                            let (lengths, angles) = ([d4, d3, d2, d1], [a34, a23, a12]);
                            (
                                five(lengths, angles, fixed, Path::Cis) - DIST15_TOL,
                                five(lengths, angles, fixed, Path::Trans) + DIST15_TOL,
                            )
                        }
                        (_, Path::Other) => (
                            five(lengths, angles, path, Path::Cis) - DIST15_TOL,
                            five(lengths, angles, path, Path::Trans) + DIST15_TOL,
                        ),
                        (_, _) => {
                            let distance = five(lengths, angles, path, next);
                            (distance - DIST15_TOL, distance + DIST15_TOL)
                        }
                    };
                    matrix.widen(a, e, lower, upper);
                    set[a * n + e] = true;
                    set[e * n + a] = true;
                }
            }
        }
    }

    /// RDKit's `setLowerBoundVDW`: every pair nothing else bounds below kept
    /// apart by its van der Waals radii, scaled four and five bonds apart.
    fn set_vdw(&self, matrix: &mut Matrix) {
        let elements = &self.molecule.elements;
        for i in 1..elements.len() {
            for j in 0..i {
                if matrix.lower(i, j) < DIST12_DELTA {
                    let radii = f64::from(vdw_radius(elements[i]) + vdw_radius(elements[j]));
                    let scale = match self.hops[i][j] {
                        4 => VDW_SCALE_15,
                        5 => VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15),
                        _ => 1.0,
                    };
                    matrix.set_lower(i, j, scale * radii);
                }
            }
        }
    }

    /// The bounds as PoseBusters reads them: bonds, the pairs sharing a
    /// neighbour, and every other pair's lower bound, with the bonds of
    /// conjugated terminal groups (a carboxylate's two oxygens) bounded alike.
    fn template(&self, matrix: &Matrix) -> Template {
        let molecule = self.molecule;
        let terminal = self.terminal_groups();
        let mut pooled: HashMap<(String, String), (f64, f64)> = HashMap::new();
        for &(a, b) in &terminal {
            let entry = pooled
                .entry(element_pair(molecule, a, b))
                .or_insert((f64::MAX, f64::MIN));
            entry.0 = entry.0.min(matrix.lower(a, b));
            entry.1 = entry.1.max(matrix.upper(a, b));
        }
        let bonds = self
            .bonds
            .iter()
            .map(|&(a, b)| {
                let (lower, upper) = if terminal.contains(&(a, b)) {
                    pooled[&element_pair(molecule, a, b)]
                } else {
                    (matrix.lower(a, b), matrix.upper(a, b))
                };
                ((a, b), lower as f32, upper as f32)
            })
            .collect();
        let mut angles = Vec::new();
        let mut apart = Vec::new();
        for a in 0..molecule.len() {
            for b in a + 1..molecule.len() {
                match self.hops[a][b] {
                    1 => {}
                    2 => {
                        angles.push(((a, b), matrix.lower(a, b) as f32, matrix.upper(a, b) as f32))
                    }
                    _ => apart.push(((a, b), matrix.lower(a, b) as f32)),
                }
            }
        }

        let multiple =
            |a: usize, b: usize| !matches!(molecule.order(a, b), Some(BondOrder::Single) | None);
        let puckered = puckered_rings(
            &molecule.elements,
            &molecule.neighbours,
            &molecule.aromatic,
            multiple,
        );

        Template::of_bounds(bonds, angles, apart, self.flat_groups(), puckered)
    }

    /// PoseBusters' conjugated terminal bonds,
    /// `[O,N;D1;$([O,N;D1]-[*]=[O,N;D1]),$([O,N;D1]=[*]-[O,N;D1])]~[*]`, each
    /// lesser atom first.
    fn terminal_groups(&self) -> HashSet<(usize, usize)> {
        let molecule = self.molecule;
        let lone = |atom: usize| {
            matches!(molecule.element(atom).number, 7 | 8)
                && molecule.neighbours[atom].len() == 1
                && molecule.hydrogens(atom) == 0
        };
        let mut found = HashSet::new();
        for end in (0..molecule.len()).filter(|&atom| lone(atom)) {
            let centre = molecule.neighbours[end][0];
            let own = molecule.order(end, centre);
            let partner = molecule.neighbours[centre].iter().any(|&other| {
                other != end
                    && lone(other)
                    && match own {
                        Some(BondOrder::Single) => {
                            molecule.order(centre, other) == Some(BondOrder::Double)
                        }
                        Some(BondOrder::Double) => {
                            molecule.order(centre, other) == Some(BondOrder::Single)
                        }
                        _ => false,
                    }
            });
            if partner {
                found.insert(key(end, centre));
            }
        }

        found
    }

    /// PoseBusters' flat groups: each five- or six-membered ring of aromatic
    /// sp2 atoms no smaller ring holds, and each double bond between sp2
    /// carbons of three substituents with the heavy atoms bonded to it.
    fn flat_groups(&self) -> Vec<Vec<usize>> {
        let molecule = self.molecule;
        let rings = &molecule.rings;
        let mut groups: Vec<Vec<usize>> = rings
            .all()
            .iter()
            .filter(|ring| {
                let size = ring.len();
                (5..=6).contains(&size)
                    && ring.iter().all(|&atom| {
                        molecule.aromatic[atom]
                            && molecule.hybrids[atom] == Hybrid::Sp2
                            && rings.smallest_of_atom(atom) == Some(size)
                    })
            })
            .cloned()
            .collect();
        for &(a, b) in &self.bonds {
            let trigonal = |atom: usize| {
                molecule.element(atom).number == 6
                    && molecule.total_degree(atom) == 3
                    && molecule.hybrids[atom] == Hybrid::Sp2
            };
            if molecule.order(a, b) != Some(BondOrder::Double) || !trigonal(a) || !trigonal(b) {
                continue;
            }
            let mut group: Vec<usize> = [a, b]
                .into_iter()
                .chain(molecule.neighbours[a].iter().copied())
                .chain(molecule.neighbours[b].iter().copied())
                .collect();
            group.sort_unstable();
            group.dedup();
            if group.len() >= 4 {
                groups.push(group);
            }
        }

        groups
    }
}

fn element_pair(molecule: &Molecule, a: usize, b: usize) -> (String, String) {
    let (x, y) = (
        molecule.elements[a].to_string(),
        molecule.elements[b].to_string(),
    );
    if x <= y { (x, y) } else { (y, x) }
}

/// RDKit's `_setRingAngle`: a ring's angle at an atom of this hybridisation.
fn ring_angle(hybrid: Hybrid, size: usize) -> f64 {
    if (hybrid == Hybrid::Sp2 && size <= 8) || size == 3 || size == 4 {
        PI * (1.0 - 2.0 / size as f64)
    } else {
        match hybrid {
            Hybrid::Sp3 if size == 5 => 104.0_f64.to_radians(),
            Hybrid::Sp3 => 109.5_f64.to_radians(),
            Hybrid::Sp3d => 105.0_f64.to_radians(),
            Hybrid::Sp3d2 => 90.0_f64.to_radians(),
            _ => 120.0_f64.to_radians(),
        }
    }
}

/// How many bonds each atom is from `start`, `usize::MAX` for another piece.
fn hops_from(start: usize, neighbours: &[Vec<usize>]) -> Vec<usize> {
    let mut hops = vec![usize::MAX; neighbours.len()];
    hops[start] = 0;
    let mut queue = VecDeque::from([start]);
    while let Some(atom) = queue.pop_front() {
        for &next in &neighbours[atom] {
            if hops[next] == usize::MAX {
                hops[next] = hops[atom] + 1;
                queue.push_back(next);
            }
        }
    }

    hops
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::components::{Component, Dictionary};
    use crate::fixed_str::FixedStr;

    /// A dictionary component as a molecule, its atoms in the file's order,
    /// at its ideal coordinates.
    fn molecule_of(component: &Component) -> (Molecule, Vec<Vec3>) {
        let names: Vec<&str> = component.atoms_in_order().collect();
        let index: HashMap<&str, usize> = names
            .iter()
            .enumerate()
            .map(|(at, &name)| (name, at))
            .collect();
        let bonds: Vec<(usize, usize, BondOrder)> = names
            .iter()
            .enumerate()
            .flat_map(|(at, name)| {
                let index = &index;
                component
                    .bonded(name)
                    .map(move |(other, order)| (at, index[other], order))
                    .filter(|&(a, b, _)| b > a)
            })
            .collect();
        let molecule = Molecule::new(
            names
                .iter()
                .map(|name| FixedStr::uppercased(component.element_of(name).unwrap()))
                .collect(),
            &bonds,
            names
                .iter()
                .map(|name| component.hydrogen_count(name) as u8)
                .collect(),
            names
                .iter()
                .map(|name| component.charge(name).unwrap_or(0))
                .collect(),
        );
        let positions = names
            .iter()
            .map(|name| component.ideal(name).unwrap())
            .collect();

        (molecule, positions)
    }

    #[test]
    fn bounds_match_rdkit_s_for_dictionary_ligands() {
        let fixtures = [
            ("BEN", include_str!("../tests/data/BEN.cif")),
            ("NAG", include_str!("../tests/data/NAG.cif")),
            ("GPP", include_str!("../tests/data/GPP.cif")),
        ];
        let expected = include_str!("../tests/data/rdkit-bounds.tsv");
        for (code, text) in fixtures {
            let dictionary = Dictionary::parse(text);
            let (molecule, positions) = molecule_of(dictionary.get(code).unwrap());
            let matrix = Builder::new(&molecule, &positions).bounds();
            let rows = expected
                .lines()
                .filter(|line| line.starts_with(code))
                .map(|line| line.split('\t').collect::<Vec<&str>>());
            let mut checked = 0;
            for row in rows {
                let [_, a, b, apart, lower, upper] = row[..] else {
                    panic!("a row of six fields");
                };
                let (a, b): (usize, usize) = (a.parse().unwrap(), b.parse().unwrap());
                let lower: f64 = lower.parse().unwrap();
                assert!(
                    (matrix.lower(a, b) - lower).abs() < 0.005,
                    "{code} {a}-{b} lower {} against RDKit's {lower}",
                    matrix.lower(a, b)
                );
                if apart.parse::<usize>().unwrap() <= 2 {
                    let upper: f64 = upper.parse().unwrap();
                    assert!(
                        (matrix.upper(a, b) - upper).abs() < 0.005,
                        "{code} {a}-{b} upper {} against RDKit's {upper}",
                        matrix.upper(a, b)
                    );
                }
                checked += 1;
            }
            assert_eq!(checked, molecule.len() * (molecule.len() - 1) / 2, "{code}");
        }
    }

    #[test]
    fn an_aromatic_bond_takes_uff_s_resonant_length() {
        let length = rest_length(1.5, uff("C_R").unwrap(), uff("C_R").unwrap());

        // 2 x 0.729 less the Pauling correction for order 1.5.
        assert!((length - 1.3793).abs() < 1e-4, "{length}");
    }

    /// A six-membered ring of `elements` as a chair, its atoms 0.25 A above
    /// and below their mean plane in turn, or flat with `pucker` 0, with
    /// `extra` atoms bonded on as given.
    fn six_ring(
        elements: &[&str],
        pucker: f32,
        bonds: &[(usize, usize, BondOrder)],
        hydrogens: Vec<u8>,
    ) -> (Molecule, Vec<Vec3>) {
        let mut positions: Vec<Vec3> = (0..6)
            .map(|at| {
                let turn = at as f32 * std::f32::consts::PI / 3.0;
                let side = if at % 2 == 0 { pucker } else { -pucker };
                Vec3::new(1.45 * turn.cos(), 1.45 * turn.sin(), side)
            })
            .collect();
        for at in 6..elements.len() {
            positions.push(positions[at - 6] * 1.9);
        }
        let mut all: Vec<(usize, usize, BondOrder)> = (0..6)
            .map(|at| (at, (at + 1) % 6, BondOrder::Single))
            .collect();
        for &(a, b, order) in bonds {
            match all.iter_mut().find(|bond| (bond.0, bond.1) == (a, b)) {
                Some(bond) => bond.2 = order,
                None => all.push((a, b, order)),
            }
        }
        let molecule = Molecule::new(
            elements
                .iter()
                .map(|&element| FixedStr::new(element))
                .collect(),
            &all,
            hydrogens,
            vec![0; elements.len()],
        );

        (molecule, positions)
    }

    #[test]
    fn a_saturated_six_ring_held_flat_fails_and_a_chair_passes() {
        for (elements, hydrogens) in [
            (vec!["C"; 6], vec![2; 6]),
            (vec!["N", "C", "C", "O", "C", "C"], vec![1, 2, 2, 0, 2, 2]),
        ] {
            let (chair, positions) = six_ring(&elements, 0.25, &[], hydrogens.clone());
            let checked = template(&chair, &positions);
            assert_eq!(checked.check(&positions, true).flat_rings, Some(0));

            let (_, flat) = six_ring(&elements, 0.0, &[], hydrogens);
            assert_eq!(
                checked.check(&flat, true).flat_rings,
                Some(1),
                "{elements:?}"
            );
        }
    }

    #[test]
    fn a_ring_with_one_double_bond_must_pucker_but_an_aromatic_one_need_not() {
        let (cyclohexene, puckered) = six_ring(
            &["C"; 6],
            0.25,
            &[(0, 1, BondOrder::Double)],
            vec![1, 1, 2, 2, 2, 2],
        );
        let (_, flat) = six_ring(&["C"; 6], 0.0, &[], vec![2; 6]);
        assert_eq!(
            template(&cyclohexene, &puckered)
                .check(&flat, true)
                .flat_rings,
            Some(1)
        );

        let kekule = [
            (0, 1, BondOrder::Double),
            (2, 3, BondOrder::Double),
            (4, 5, BondOrder::Double),
        ];
        let (benzene, flat) = six_ring(&["C"; 6], 0.0, &kekule, vec![1; 6]);
        assert_eq!(
            template(&benzene, &flat).check(&flat, true).flat_rings,
            Some(0)
        );
    }

    #[test]
    fn a_ring_fused_to_another_is_not_held_to_pucker() {
        // Decalin flattened: its two fused atoms each lie in both rings.
        let (decalin, flat) = six_ring(
            &["C"; 10],
            0.0,
            &[
                (0, 6, BondOrder::Single),
                (6, 7, BondOrder::Single),
                (7, 8, BondOrder::Single),
                (8, 9, BondOrder::Single),
                (9, 1, BondOrder::Single),
            ],
            vec![1, 1, 2, 2, 2, 2, 2, 2, 2, 2],
        );

        assert_eq!(
            template(&decalin, &flat).check(&flat, true).flat_rings,
            Some(0)
        );
    }

    #[test]
    fn a_stretched_bond_and_a_folded_chain_fail_their_checks() {
        // Hexane in a flat zigzag: bonds of 1.53 A and angles of 111.7 degrees.
        let zigzag: Vec<Vec3> = (0..6)
            .map(|at| Vec3::new(at as f32 * 1.27, (at % 2) as f32 * 0.86, 0.0))
            .collect();
        let hexane = Molecule::new(
            vec![FixedStr::new("C"); 6],
            &(0..5)
                .map(|at| (at, at + 1, BondOrder::Single))
                .collect::<Vec<_>>(),
            vec![3, 2, 2, 2, 2, 3],
            vec![0; 6],
        );
        let ideal = template(&hexane, &zigzag);
        let clean = crate::geometry::Geometry {
            unflat: Some(0),
            flat_rings: Some(0),
            ..Default::default()
        };
        assert_eq!(ideal.check(&zigzag, true), clean);

        let mut stretched = zigzag.clone();
        for position in &mut stretched[3..] {
            *position += Vec3::X * 0.8;
        }
        let failed = ideal.check(&stretched, true);
        assert_eq!(failed.bonds, 1);

        // The last atom folded back onto the first.
        let mut folded = zigzag.clone();
        folded[5] = zigzag[0] + Vec3::new(0.3, 0.4, 0.9);
        assert!(ideal.check(&folded, true).clashes > 0);
    }
}
