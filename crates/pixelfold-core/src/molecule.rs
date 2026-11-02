//! A ligand's chemistry as RDKit models it once it has sanitized a molecule:
//! which rings are aromatic, which bonds conjugated, and how each atom is
//! hybridised. PoseBusters bounds a pose's geometry from these, so the rules
//! follow RDKit 2024.09 (`Aromaticity.cpp`, `ConjugHybrid.cpp`) rather than a
//! textbook: an amide nitrogen or an ester oxygen next to a double bond is sp2,
//! a sulfone sulfur sp3, and 2-pyridone aromatic.
//!
//! Hydrogens are counts on the heavy atoms, not atoms. Bond orders are those
//! of a Kekule structure, as the Chemical Component Dictionary and most SD
//! files write them; a bond already written aromatic stays so.

use std::collections::{HashMap, HashSet};

use crate::components::BondOrder;
use crate::fixed_str::FixedStr;
use crate::rings::{Rings, holds_bond};

/// What RDKit's periodic table says of an element: its atomic number, the
/// electrons in its outer shell, and its valences, the first the default.
/// A metal has no valences: any will do.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Element {
    pub(crate) number: u8,
    pub(crate) outer: u8,
    valences: &'static [u8],
}

impl Element {
    /// The element of `symbol`, read without regard to case.
    pub(crate) fn of(symbol: &str) -> Element {
        let (number, outer, valences): (u8, u8, &'static [u8]) =
            match symbol.to_ascii_uppercase().as_str() {
                "H" | "D" | "T" => (1, 1, &[1]),
                "B" => (5, 3, &[3]),
                "C" => (6, 4, &[4]),
                "N" => (7, 5, &[3]),
                "O" => (8, 6, &[2]),
                "F" => (9, 7, &[1]),
                "NA" => (11, 1, &[1]),
                "MG" => (12, 2, &[2]),
                "AL" => (13, 3, &[3]),
                "SI" => (14, 4, &[4]),
                "P" => (15, 5, &[3, 5, 7]),
                "S" => (16, 6, &[2, 4, 6]),
                "CL" => (17, 7, &[1]),
                "K" => (19, 1, &[1]),
                "CA" => (20, 2, &[2]),
                "AS" => (33, 5, &[3, 5, 7]),
                "SE" => (34, 6, &[2, 4, 6]),
                "BR" => (35, 7, &[1]),
                "I" => (53, 7, &[1, 3, 5]),
                "MN" => (25, 7, &[]),
                "FE" => (26, 8, &[]),
                "CO" => (27, 9, &[]),
                "NI" => (28, 10, &[]),
                "CU" => (29, 11, &[]),
                "ZN" => (30, 2, &[]),
                _ => (0, 0, &[]),
            };
        Element {
            number,
            outer,
            valences,
        }
    }

    /// The default valence, or `None` for an element any valence suits.
    pub(crate) fn default_valence(&self) -> Option<u8> {
        self.valences.first().copied()
    }

    /// Whether this element pulls electrons harder than `other`, as RDKit
    /// orders them: by outer electrons, then the lighter.
    fn more_electronegative(&self, other: &Element) -> bool {
        self.outer > other.outer || (self.outer == other.outer && self.number < other.number)
    }
}

/// How RDKit hybridises an atom.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Hybrid {
    S,
    Sp,
    Sp2,
    Sp3,
    Sp3d,
    Sp3d2,
    Unspecified,
}

/// A molecule's heavy atoms and bonds with what RDKit perceives of them.
#[derive(Clone, Debug)]
pub(crate) struct Molecule {
    pub(crate) elements: Vec<FixedStr<2>>,
    pub(crate) neighbours: Vec<Vec<usize>>,
    /// Each bond's order once aromaticity is perceived, keyed lesser atom
    /// first.
    orders: HashMap<(usize, usize), BondOrder>,
    hydrogens: Vec<u8>,
    charges: Vec<i8>,
    /// Each atom's valence in the Kekule structure, hydrogens included.
    valences: Vec<u8>,
    pub(crate) rings: Rings,
    pub(crate) aromatic: Vec<bool>,
    conjugated: HashSet<(usize, usize)>,
    pub(crate) hybrids: Vec<Hybrid>,
}

/// A bond keyed lesser atom first.
pub(crate) fn key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

impl Molecule {
    /// The molecule of heavy atoms with these `elements`, `bonds` between
    /// them, `hydrogens` on each and formal `charges`, perceived as RDKit
    /// sanitizes it.
    pub(crate) fn new(
        elements: Vec<FixedStr<2>>,
        bonds: &[(usize, usize, BondOrder)],
        hydrogens: Vec<u8>,
        charges: Vec<i8>,
    ) -> Molecule {
        let mut neighbours = vec![Vec::new(); elements.len()];
        let mut orders = HashMap::new();
        for &(a, b, order) in bonds {
            if a != b && orders.insert(key(a, b), order).is_none() {
                neighbours[a].push(b);
                neighbours[b].push(a);
            }
        }
        let rings = Rings::of(&neighbours);
        let mut molecule = Molecule {
            valences: Vec::new(),
            aromatic: vec![false; elements.len()],
            conjugated: HashSet::new(),
            hybrids: Vec::new(),
            elements,
            neighbours,
            orders,
            hydrogens,
            charges,
            rings,
        };
        molecule.valences = kekule_valences(&molecule.elements, bonds, &molecule.hydrogens);
        molecule.perceive_aromaticity();
        molecule.conjugated = molecule.conjugation();
        molecule.hybrids = (0..molecule.len())
            .map(|atom| molecule.hybridisation(atom))
            .collect();

        molecule
    }

    pub(crate) fn len(&self) -> usize {
        self.elements.len()
    }

    pub(crate) fn element(&self, atom: usize) -> Element {
        Element::of(self.elements[atom].as_str())
    }

    pub(crate) fn order(&self, a: usize, b: usize) -> Option<BondOrder> {
        self.orders.get(&key(a, b)).copied()
    }

    pub(crate) fn hydrogens(&self, atom: usize) -> u8 {
        self.hydrogens[atom]
    }

    pub(crate) fn valence(&self, atom: usize) -> u8 {
        self.valences[atom]
    }

    pub(crate) fn is_conjugated(&self, a: usize, b: usize) -> bool {
        self.conjugated.contains(&key(a, b))
    }

    pub(crate) fn has_conjugated_bond(&self, atom: usize) -> bool {
        self.neighbours[atom]
            .iter()
            .any(|&other| self.is_conjugated(atom, other))
    }

    /// Heavy neighbours and hydrogens.
    pub(crate) fn total_degree(&self, atom: usize) -> usize {
        self.neighbours[atom].len() + self.hydrogens[atom] as usize
    }

    /// The bond's share of an atom's valence: half bonds summed are three for
    /// an aromatic bond.
    fn contribution(&self, a: usize, b: usize) -> f32 {
        match self.order(a, b) {
            Some(BondOrder::Single) | None => 1.0,
            Some(BondOrder::Double) => 2.0,
            Some(BondOrder::Triple) => 3.0,
            Some(BondOrder::Aromatic) => 1.5,
        }
    }

    /// The bond orders summed, less the bonds: an atom's unsaturations.
    fn unsaturations(&self, atom: usize) -> u8 {
        self.valences[atom]
            .saturating_sub(self.hydrogens[atom])
            .saturating_sub(self.neighbours[atom].len() as u8)
    }

    /// RDKit's `countAtomElec`: the electrons an atom could give a pi system,
    /// or `None` when it can be neither aromatic nor conjugated.
    fn pi_electrons(&self, atom: usize) -> Option<i32> {
        let element = self.element(atom);
        let default = i32::from(element.default_valence()?);
        if default <= 1 {
            return None;
        }
        let degree = self.total_degree(atom) as i32;
        if degree > 3 {
            return None;
        }
        let lone = (i32::from(element.outer) - default - i32::from(self.charges[atom])).max(0);
        let electrons = default - degree + lone;
        if electrons > 1 && self.unsaturations(atom) > 1 {
            return Some(1);
        }

        Some(electrons)
    }

    /// RDKit's aromaticity model: each ring, or set of up to six fused rings,
    /// of candidate atoms whose pi electrons, counted at their perimeter, can
    /// make 4n + 2.
    fn perceive_aromaticity(&mut self) {
        let donors: Vec<Option<Donor>> = (0..self.len()).map(|atom| self.donor(atom)).collect();
        let candidates: Vec<Vec<usize>> = self
            .rings
            .all()
            .iter()
            .filter(|ring| {
                ring.iter()
                    .all(|&atom| donors[atom].is_some() && self.aromatic_candidate(atom))
            })
            .cloned()
            .collect();
        let mut aromatic_bonds: HashSet<(usize, usize)> = self
            .orders
            .iter()
            .filter(|(_, order)| matches!(order, BondOrder::Aromatic))
            .map(|(&bond, _)| bond)
            .collect();
        for (a, b) in aromatic_bonds.clone() {
            self.aromatic[a] = true;
            self.aromatic[b] = true;
        }

        for system in fused_systems(&candidates) {
            let mut done: HashSet<(usize, usize)> = HashSet::new();
            let bonds_in_system: HashSet<(usize, usize)> = system
                .iter()
                .flat_map(|&ring| ring_bonds(&candidates[ring]))
                .collect();
            'sizes: for size in 1..=system.len().min(6) {
                for combination in combinations(system.len(), size) {
                    if done.len() >= bonds_in_system.len() {
                        break 'sizes;
                    }
                    let chosen: Vec<&Vec<usize>> = combination
                        .iter()
                        .map(|&at| &candidates[system[at]])
                        .collect();
                    if !connected(&chosen) {
                        continue;
                    }
                    let mut counts: HashMap<usize, usize> = HashMap::new();
                    for ring in &chosen {
                        for &atom in ring.iter() {
                            *counts.entry(atom).or_default() += 1;
                        }
                    }
                    let perimeter: Vec<usize> = counts
                        .iter()
                        .filter(|&(_, &count)| count <= 2)
                        .map(|(&atom, _)| atom)
                        .collect();
                    if !huckel(
                        perimeter
                            .iter()
                            .map(|&atom| donors[atom].expect("a candidate")),
                    ) {
                        continue;
                    }
                    let mut bond_counts: HashMap<(usize, usize), usize> = HashMap::new();
                    for ring in &chosen {
                        for &atom in ring.iter() {
                            self.aromatic[atom] = true;
                        }
                        for bond in ring_bonds(ring) {
                            *bond_counts.entry(bond).or_default() += 1;
                        }
                    }
                    for (bond, count) in bond_counts {
                        if count == 1 {
                            aromatic_bonds.insert(bond);
                            done.insert(bond);
                        }
                    }
                }
            }
        }

        for bond in aromatic_bonds {
            if let Some(order) = self.orders.get_mut(&bond)
                && matches!(order, BondOrder::Single | BondOrder::Double)
            {
                *order = BondOrder::Aromatic;
            }
        }
    }

    /// RDKit's `isAtomCandForArom`.
    fn aromatic_candidate(&self, atom: usize) -> bool {
        let element = self.element(atom);
        if element.number > 18 && element.number != 34 && element.number != 52 {
            return false;
        }
        let charge = i32::from(self.charges[atom]);
        let shifted = Element::of(symbol_of(i32::from(element.number) - charge));
        if element.default_valence().is_some()
            && shifted
                .default_valence()
                .is_some_and(|default| self.valences[atom] > default)
        {
            return false;
        }
        if self.unsaturations(atom) > 1 {
            let multiple = self.neighbours[atom]
                .iter()
                .filter(|&&other| {
                    matches!(
                        self.order(atom, other),
                        Some(BondOrder::Double | BondOrder::Triple)
                    )
                })
                .count();
            if multiple > 1 {
                return false;
            }
        }

        true
    }

    /// RDKit's `getAtomDonorTypeArom`, or `None` for an atom that donates
    /// nothing to a ring and so cannot be aromatic.
    fn donor(&self, atom: usize) -> Option<Donor> {
        let electrons = self.pi_electrons(atom)?;
        let exocyclic = self.neighbours[atom].iter().copied().find(|&other| {
            self.rings.of_bond(atom, other).next().is_none()
                && matches!(
                    self.order(atom, other),
                    Some(BondOrder::Double | BondOrder::Triple)
                )
        });
        let cyclic_multiple = self.neighbours[atom].iter().any(|&other| {
            self.rings.of_bond(atom, other).next().is_some()
                && matches!(
                    self.order(atom, other),
                    Some(BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic)
                )
        });
        let steals = |other: usize| {
            self.element(other)
                .more_electronegative(&self.element(atom))
        };

        match electrons {
            0 => {
                if exocyclic.is_some() {
                    Some(Donor::Vacant)
                } else if cyclic_multiple {
                    Some(Donor::One)
                } else {
                    None
                }
            }
            1 => {
                if let Some(other) = exocyclic {
                    Some(if steals(other) {
                        Donor::Vacant
                    } else {
                        Donor::One
                    })
                } else if self.unsaturations(atom) > 0 || cyclic_multiple {
                    Some(Donor::One)
                } else if self.charges[atom] == 1 {
                    Some(Donor::Vacant)
                } else {
                    None
                }
            }
            _ => {
                let electrons = match exocyclic {
                    Some(other) if steals(other) => electrons - 1,
                    _ => electrons,
                };
                Some(if electrons % 2 == 1 {
                    Donor::One
                } else {
                    Donor::Two
                })
            }
        }
    }

    /// RDKit's `isAtomConjugCand`.
    fn conjugation_candidate(&self, atom: usize) -> bool {
        let element = self.element(atom);
        if self.charges[atom] == 0
            && element
                .default_valence()
                .is_some_and(|default| self.valences[atom] > default)
        {
            return false;
        }
        let outer = element.outer;
        (element.number <= 10
            || (outer != 5 && outer != 6)
            || (outer == 6 && self.total_degree(atom) < 2))
            && self
                .pi_electrons(atom)
                .is_some_and(|electrons| electrons > 0)
    }

    /// RDKit's `setConjugation`: aromatic bonds, and around each atom with two
    /// or three substituents and a multiple bond to a candidate, that bond and
    /// every other bond to a candidate with at most three substituents.
    fn conjugation(&self) -> HashSet<(usize, usize)> {
        let mut conjugated: HashSet<(usize, usize)> = self
            .orders
            .iter()
            .filter(|(_, order)| matches!(order, BondOrder::Aromatic))
            .map(|(&bond, _)| bond)
            .collect();
        for atom in 0..self.len() {
            if !self.conjugation_candidate(atom) || !(2..=3).contains(&self.total_degree(atom)) {
                continue;
            }
            for &first in &self.neighbours[atom] {
                if self.contribution(atom, first) < 1.5 || !self.conjugation_candidate(first) {
                    continue;
                }
                for &second in &self.neighbours[atom] {
                    if second != first
                        && self.total_degree(second) <= 3
                        && self.conjugation_candidate(second)
                    {
                        conjugated.insert(key(atom, first));
                        conjugated.insert(key(atom, second));
                    }
                }
            }
        }

        conjugated
    }

    /// RDKit's `setHybridization`, from the bonds and lone pairs about an
    /// atom; four of them make sp2 rather than sp3 where the atom has a
    /// conjugated bond and at most three substituents (the second oxygen of
    /// O=CO, an amide nitrogen).
    fn hybridisation(&self, atom: usize) -> Hybrid {
        let element = self.element(atom);
        if element.number == 0 {
            return Hybrid::Unspecified;
        }
        let degree = self.total_degree(atom) as i32;
        let orbitals = if element.number <= 1 {
            degree
        } else {
            let outer = i32::from(element.outer);
            let valence = i32::from(self.valences[atom]);
            let charge = i32::from(self.charges[atom]);
            degree + (outer - (valence + charge)) / 2
        };

        match orbitals {
            0 | 1 => Hybrid::S,
            2 => Hybrid::Sp,
            3 => Hybrid::Sp2,
            4 if degree <= 3 && self.has_conjugated_bond(atom) => Hybrid::Sp2,
            4 => Hybrid::Sp3,
            5 => Hybrid::Sp3d,
            6 => Hybrid::Sp3d2,
            _ => Hybrid::Unspecified,
        }
    }
}

/// How many electrons an atom gives an aromatic ring.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Donor {
    Vacant,
    One,
    Two,
}

/// Whether atoms giving these electrons make an aromatic ring: some count
/// they allow is 4n + 2, or, with fewer than six, exactly two.
fn huckel(donors: impl Iterator<Item = Donor>) -> bool {
    let electrons: i32 = donors
        .map(|donor| match donor {
            Donor::Vacant => 0,
            Donor::One => 1,
            Donor::Two => 2,
        })
        .sum();

    if electrons >= 6 {
        (electrons - 2) % 4 == 0
    } else {
        electrons == 2
    }
}

/// A ring's bonds, lesser atom first.
fn ring_bonds(ring: &[usize]) -> Vec<(usize, usize)> {
    (0..ring.len())
        .map(|at| key(ring[at], ring[(at + 1) % ring.len()]))
        .collect()
}

/// The candidate rings grouped into systems of rings that share a bond, each
/// system as indices into `rings`.
fn fused_systems(rings: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut system_of = vec![usize::MAX; rings.len()];
    let mut systems: Vec<Vec<usize>> = Vec::new();
    for start in 0..rings.len() {
        if system_of[start] != usize::MAX {
            continue;
        }
        let id = systems.len();
        let mut members = vec![start];
        system_of[start] = id;
        let mut at = 0;
        while at < members.len() {
            let ring = members[at];
            for other in 0..rings.len() {
                if system_of[other] == usize::MAX && share_bond(&rings[ring], &rings[other]) {
                    system_of[other] = id;
                    members.push(other);
                }
            }
            at += 1;
        }
        members.sort_unstable();
        systems.push(members);
    }

    systems
}

fn share_bond(first: &[usize], second: &[usize]) -> bool {
    ring_bonds(first)
        .into_iter()
        .any(|(a, b)| holds_bond(second, a, b))
}

/// Whether the `rings` form one fused system.
fn connected(rings: &[&Vec<usize>]) -> bool {
    let mut reached = vec![false; rings.len()];
    reached[0] = true;
    let mut queue = vec![0];
    while let Some(ring) = queue.pop() {
        for other in 0..rings.len() {
            if !reached[other] && share_bond(rings[ring], rings[other]) {
                reached[other] = true;
                queue.push(other);
            }
        }
    }

    reached.into_iter().all(|reached| reached)
}

/// Every choice of `size` of `count` items, each in increasing order, in
/// lexicographic order.
fn combinations(count: usize, size: usize) -> Vec<Vec<usize>> {
    fn extend(
        start: usize,
        count: usize,
        left: usize,
        current: &mut Vec<usize>,
        out: &mut Vec<Vec<usize>>,
    ) {
        if left == 0 {
            out.push(current.clone());
            return;
        }
        for next in start..=count - left {
            current.push(next);
            extend(next + 1, count, left - 1, current, out);
            current.pop();
        }
    }

    let mut out = Vec::new();
    if size <= count {
        extend(0, count, size, &mut Vec::new(), &mut out);
    }
    out
}

/// Each atom's valence in the Kekule structure, its `hydrogens` included,
/// without perceiving anything else. A bond written aromatic counts one, and
/// its atom gains a double bond's worth more unless it gives the ring its
/// lone pair (a pyrrole nitrogen, a furan oxygen).
pub(crate) fn kekule_valences(
    elements: &[FixedStr<2>],
    bonds: &[(usize, usize, BondOrder)],
    hydrogens: &[u8],
) -> Vec<u8> {
    let mut whole = vec![0u8; elements.len()];
    let mut aromatic = vec![0u8; elements.len()];
    let mut degree = vec![0usize; elements.len()];
    let mut seen = HashSet::new();
    for &(a, b, order) in bonds {
        if a == b || !seen.insert(key(a, b)) {
            continue;
        }
        for atom in [a, b] {
            degree[atom] += 1;
            match order {
                BondOrder::Aromatic => aromatic[atom] += 1,
                BondOrder::Double => whole[atom] += 2,
                BondOrder::Triple => whole[atom] += 3,
                BondOrder::Single => whole[atom] += 1,
            }
        }
    }

    (0..elements.len())
        .map(|atom| {
            let gives_lone_pair = match Element::of(elements[atom].as_str()).number {
                7 | 15 => degree[atom] + usize::from(hydrogens[atom]) == 3,
                8 | 16 | 34 => true,
                _ => false,
            };
            let pi = u8::from(aromatic[atom] > 0 && !gives_lone_pair);
            whole[atom]
                .saturating_add(aromatic[atom])
                .saturating_add(pi)
                .saturating_add(hydrogens[atom])
        })
        .collect()
}

/// The hydrogens RDKit gives an atom of `symbol` and `charge` whose bonds,
/// its hydrogens already present included, sum to `valence`: up to the least
/// valence it allows at or above that, a charged atom allowing those of the
/// element it is isoelectronic with (an N+ is a carbon, an O- a fluorine).
pub(crate) fn implicit_hydrogens(symbol: &str, charge: i8, valence: u8) -> u8 {
    let element = Element::of(symbol);
    let allowed = match charge {
        0 => element.valences,
        _ => Element::of(symbol_of(i32::from(element.number) - i32::from(charge))).valences,
    };

    allowed
        .iter()
        .find(|&&allowed| allowed >= valence)
        .map_or(0, |&allowed| allowed - valence)
}

/// The formal charge of an atom of `symbol` whose bonds and hydrogens sum to
/// `valence`, for a definition that gives none: an ammonium nitrogen is +1,
/// an oxygen with one bond and no hydrogen -1.
pub(crate) fn inferred_charge(symbol: &str, valence: u8) -> i8 {
    match (symbol.to_ascii_uppercase().as_str(), valence) {
        ("N" | "P", 4) | ("O" | "S", 3) => 1,
        ("O" | "S", 1) | ("N", 2) => -1,
        _ => 0,
    }
}

/// The symbol of an element by atomic number, for the ones a charged atom of
/// a ligand can be isoelectronic with.
fn symbol_of(number: i32) -> &'static str {
    match number {
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        14 => "SI",
        15 => "P",
        16 => "S",
        17 => "CL",
        _ => "",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use BondOrder::{Double as D, Single as S};

    fn molecule(
        elements: &[&str],
        bonds: &[(usize, usize, BondOrder)],
        hydrogens: &[u8],
    ) -> Molecule {
        Molecule::new(
            elements
                .iter()
                .map(|&element| FixedStr::new(element))
                .collect(),
            bonds,
            hydrogens.to_vec(),
            vec![0; elements.len()],
        )
    }

    fn benzene() -> Molecule {
        molecule(
            &["C"; 6],
            &[
                (0, 1, D),
                (1, 2, S),
                (2, 3, D),
                (3, 4, S),
                (4, 5, D),
                (5, 0, S),
            ],
            &[1; 6],
        )
    }

    #[test]
    fn a_kekule_benzene_is_aromatic_and_sp2() {
        let benzene = benzene();

        assert!(benzene.aromatic.iter().all(|&aromatic| aromatic));
        assert_eq!(benzene.order(0, 1), Some(BondOrder::Aromatic));
        assert_eq!(benzene.order(1, 2), Some(BondOrder::Aromatic));
        assert!(benzene.hybrids.iter().all(|&hybrid| hybrid == Hybrid::Sp2));
        assert!(benzene.is_conjugated(0, 5));
    }

    #[test]
    fn pyrrole_furan_and_thiophene_are_aromatic_by_their_lone_pairs() {
        for (heteroatom, hydrogens) in [("N", 1), ("O", 0), ("S", 0)] {
            let ring = molecule(
                &[heteroatom, "C", "C", "C", "C"],
                &[(0, 1, S), (1, 2, D), (2, 3, S), (3, 4, D), (4, 0, S)],
                &[hydrogens, 1, 1, 1, 1],
            );
            assert!(
                ring.aromatic.iter().all(|&aromatic| aromatic),
                "{heteroatom}"
            );
            assert_eq!(ring.hybrids[0], Hybrid::Sp2, "{heteroatom}");
        }
    }

    #[test]
    fn cyclohexadiene_is_not_aromatic() {
        let ring = molecule(
            &["C"; 6],
            &[
                (0, 1, D),
                (1, 2, S),
                (2, 3, D),
                (3, 4, S),
                (4, 5, S),
                (5, 0, S),
            ],
            &[1, 1, 1, 1, 2, 2],
        );

        assert!(ring.aromatic.iter().all(|&aromatic| !aromatic));
        assert_eq!(ring.hybrids[4], Hybrid::Sp3);
    }

    #[test]
    fn two_pyridone_is_aromatic_with_its_carbonyl_outside() {
        // N1 C2(=O7) C3=C4 C5=C6
        let pyridone = molecule(
            &["N", "C", "C", "C", "C", "C", "O"],
            &[
                (0, 1, S),
                (1, 2, S),
                (2, 3, D),
                (3, 4, S),
                (4, 5, D),
                (5, 0, S),
                (1, 6, D),
            ],
            &[1, 0, 1, 1, 1, 1, 0],
        );

        assert!(pyridone.aromatic[..6].iter().all(|&aromatic| aromatic));
        assert!(!pyridone.aromatic[6]);
    }

    #[test]
    fn an_amide_nitrogen_and_an_ester_oxygen_are_sp2_and_a_sulfone_sulfur_sp3() {
        // Acetamide: C0 C1(=O2) N3H2
        let amide = molecule(
            &["C", "C", "O", "N"],
            &[(0, 1, S), (1, 2, D), (1, 3, S)],
            &[3, 0, 0, 2],
        );
        assert_eq!(amide.hybrids[3], Hybrid::Sp2);
        assert_eq!(amide.hybrids[0], Hybrid::Sp3);
        assert!(amide.is_conjugated(1, 3));

        // Methyl acetate: C0 C1(=O2) O3 C4
        let ester = molecule(
            &["C", "C", "O", "O", "C"],
            &[(0, 1, S), (1, 2, D), (1, 3, S), (3, 4, S)],
            &[3, 0, 0, 0, 3],
        );
        assert_eq!(ester.hybrids[3], Hybrid::Sp2);

        // Dimethyl sulfone: C0 S1(=O2)(=O3) C4
        let sulfone = molecule(
            &["C", "S", "O", "O", "C"],
            &[(0, 1, S), (1, 2, D), (1, 3, D), (1, 4, S)],
            &[3, 0, 0, 0, 3],
        );
        assert_eq!(sulfone.hybrids[1], Hybrid::Sp3);
        assert_eq!(sulfone.valence(1), 6);
    }

    #[test]
    fn a_ketone_flanked_by_methyls_is_not_conjugated() {
        let acetone = molecule(
            &["C", "C", "O", "C"],
            &[(0, 1, S), (1, 2, D), (1, 3, S)],
            &[3, 0, 0, 3],
        );

        assert!(!acetone.is_conjugated(1, 2));
        assert_eq!(acetone.hybrids[2], Hybrid::Sp2);
        assert_eq!(acetone.hybrids[1], Hybrid::Sp2);
    }

    #[test]
    fn indole_is_aromatic_in_both_rings() {
        // Benzene 0-5 fused at 0 and 5 with the pyrrole 5,6,7(N),8,0.
        let indole = molecule(
            &["C", "C", "C", "C", "C", "C", "C", "N", "C"],
            &[
                (0, 1, D),
                (1, 2, S),
                (2, 3, D),
                (3, 4, S),
                (4, 5, D),
                (5, 0, S),
                (5, 6, S),
                (6, 8, D),
                (8, 7, S),
                (7, 0, S),
            ],
            &[0, 1, 1, 1, 1, 0, 1, 1, 1],
        );

        assert!(indole.aromatic.iter().all(|&aromatic| aromatic));
    }

    #[test]
    fn a_ring_already_written_aromatic_stays_so() {
        let ring = molecule(
            &["C"; 6],
            &[
                (0, 1, BondOrder::Aromatic),
                (1, 2, BondOrder::Aromatic),
                (2, 3, BondOrder::Aromatic),
                (3, 4, BondOrder::Aromatic),
                (4, 5, BondOrder::Aromatic),
                (5, 0, BondOrder::Aromatic),
            ],
            &[1; 6],
        );

        assert!(ring.aromatic.iter().all(|&aromatic| aromatic));
        assert_eq!(ring.valence(0), 4);
        assert!(ring.hybrids.iter().all(|&hybrid| hybrid == Hybrid::Sp2));
    }
}
