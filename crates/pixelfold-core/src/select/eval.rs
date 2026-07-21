//! Evaluate a [`Selection`] against a structure to an [`AtomSet`].

use std::collections::HashSet;

use fixedbitset::FixedBitSet;
use glam::Vec3;

use super::Selection;
use crate::fixed_str::FixedStr;
use crate::interactions::classify::{is_metal, is_water};
use crate::interactions::topology::is_standard_residue;
use crate::spatial::Grid;
use crate::structure::{Atom, Protein};

/// A set of atom indices into a structure, held as a bitset so the boolean
/// operators are word-parallel over large structures.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomSet {
    bits: FixedBitSet,
}

impl AtomSet {
    /// Evaluate `selection` against `protein`.
    pub fn of(selection: &Selection, protein: &Protein) -> Self {
        evaluate(selection, protein)
    }

    /// Whether atom `index` is in the set.
    pub fn contains(&self, index: usize) -> bool {
        self.bits.contains(index)
    }

    /// How many atoms the set holds.
    pub fn count(&self) -> usize {
        self.bits.count_ones(..)
    }

    pub fn is_empty(&self) -> bool {
        self.count() == 0
    }

    /// The selected atom indices, ascending.
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.bits.ones()
    }

    fn empty(atoms: usize) -> Self {
        Self {
            bits: FixedBitSet::with_capacity(atoms),
        }
    }

    fn full(atoms: usize) -> Self {
        let mut bits = FixedBitSet::with_capacity(atoms);
        bits.insert_range(..);
        Self { bits }
    }

    fn insert(&mut self, index: usize) {
        self.bits.insert(index);
    }

    fn union(mut self, other: &Self) -> Self {
        self.bits.union_with(&other.bits);
        self
    }

    fn intersect(mut self, other: &Self) -> Self {
        self.bits.intersect_with(&other.bits);
        self
    }

    fn complement(mut self) -> Self {
        self.bits.toggle_range(..);
        self
    }
}

type ResidueKey = (FixedStr<4>, u32, Option<char>);

fn residue_key(atom: &Atom) -> ResidueKey {
    (atom.chain_id, atom.residue_seq, atom.insertion_code)
}

fn evaluate(selection: &Selection, protein: &Protein) -> AtomSet {
    let atoms = protein.atoms.len();
    match selection {
        Selection::All => AtomSet::full(atoms),
        Selection::None => AtomSet::empty(atoms),
        Selection::Not(inner) => evaluate(inner, protein).complement(),
        Selection::And(a, b) => evaluate(a, protein).intersect(&evaluate(b, protein)),
        Selection::Or(a, b) => evaluate(a, protein).union(&evaluate(b, protein)),
        Selection::ByRes(inner) => by_residue(&evaluate(inner, protein), protein),
        Selection::Within(distance, inner) => within(*distance, &evaluate(inner, protein), protein),
        Selection::Chain(ids) => select(protein, |atom| {
            ids.iter().any(|id| atom.chain_id == id.as_str())
        }),
        Selection::Resn(names) => select(protein, |atom| {
            names
                .iter()
                .any(|name| atom.residue_name.as_str().eq_ignore_ascii_case(name))
        }),
        Selection::Name(names) => select(protein, |atom| {
            names
                .iter()
                .any(|name| atom.name.as_str().eq_ignore_ascii_case(name))
        }),
        Selection::Element(symbols) => select(protein, |atom| {
            symbols
                .iter()
                .any(|symbol| atom.element.as_str().eq_ignore_ascii_case(symbol))
        }),
        Selection::Resi(ranges) => select(protein, |atom| {
            ranges.iter().any(|range| range.contains(atom.residue_seq))
        }),
        Selection::Model(number) => select(protein, |atom| atom.model_number == *number),
        Selection::Ss(ss) => select(protein, |atom| atom.secondary_structure == *ss),
        Selection::BFactor(cmp, bound) | Selection::Plddt(cmp, bound) => {
            select(protein, |atom| cmp.test(atom.b_factor, *bound))
        }
        Selection::Protein => select(protein, |atom| {
            is_standard_residue(atom.residue_name.as_str())
        }),
        Selection::Nucleic => select(protein, |atom| is_nucleic(atom.residue_name.as_str())),
        Selection::Water => select(protein, |atom| is_water(atom.residue_name.as_str())),
        Selection::Metal => select(protein, is_metal),
        Selection::Hetatm => select(protein, |atom| atom.is_hetatm),
        Selection::Ligand => select(protein, |atom| {
            atom.is_hetatm && !is_water(atom.residue_name.as_str()) && !is_metal(atom)
        }),
    }
}

/// The atoms satisfying a per-atom predicate.
fn select(protein: &Protein, keep: impl Fn(&Atom) -> bool) -> AtomSet {
    let mut set = AtomSet::empty(protein.atoms.len());
    for (index, atom) in protein.atoms.iter().enumerate() {
        if keep(atom) {
            set.insert(index);
        }
    }

    set
}

/// Widen a set to every atom of any residue it touches.
fn by_residue(seed: &AtomSet, protein: &Protein) -> AtomSet {
    let residues: HashSet<ResidueKey> = seed
        .iter()
        .map(|index| residue_key(&protein.atoms[index]))
        .collect();

    select(protein, |atom| residues.contains(&residue_key(atom)))
}

/// Every atom within `distance` of any atom in the seed, the seed included since
/// its own atoms are zero away. The grid makes this O(N) rather than seed times
/// atoms.
fn within(distance: f32, seed: &AtomSet, protein: &Protein) -> AtomSet {
    let mut out = AtomSet::empty(protein.atoms.len());
    if seed.is_empty() {
        return out;
    }

    let positions: Vec<Vec3> = protein.atoms.iter().map(|atom| atom.position).collect();
    let cell = if distance > 0.0 { distance } else { 1.0 };
    let grid = Grid::build(&positions, cell);

    for index in seed.iter() {
        grid.for_each_within(positions[index], distance, |neighbour| {
            out.insert(neighbour)
        });
    }

    out
}

/// The standard nucleotide residue names, DNA and RNA.
fn is_nucleic(residue_name: &str) -> bool {
    matches!(
        residue_name,
        "DA" | "DC" | "DG" | "DT" | "DU" | "DI" | "A" | "C" | "G" | "U" | "I"
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::select::select as run;
    use crate::structure::SecondaryStructure;

    fn atom(residue: &str, seq: u32, name: &str, element: &str, position: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(element),
            residue_name: FixedStr::new(residue),
            residue_seq: seq,
            chain_id: FixedStr::new("A"),
            is_hetatm: false,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    /// A small structure: an alanine, a water, and a zinc ion, all on a line.
    fn sample() -> Protein {
        let mut ala_n = atom("ALA", 1, "N", "N", Vec3::new(0.0, 0.0, 0.0));
        ala_n.b_factor = 20.0;
        let mut ala_ca = atom("ALA", 1, "CA", "C", Vec3::new(1.5, 0.0, 0.0));
        ala_ca.b_factor = 80.0;
        ala_ca.secondary_structure = SecondaryStructure::Helix;
        let mut water = atom("HOH", 2, "O", "O", Vec3::new(4.0, 0.0, 0.0));
        water.is_hetatm = true;
        let mut zinc = atom("ZN", 3, "ZN", "Zn", Vec3::new(10.0, 0.0, 0.0));
        zinc.is_hetatm = true;

        Protein {
            atoms: vec![ala_n, ala_ca, water, zinc],
            title: String::new(),
            surface_points: Vec::new(),
            hbonds: Vec::new(),
            assembly: None,
        }
    }

    fn selected(query: &str) -> Vec<usize> {
        run(query, &sample()).unwrap().iter().collect()
    }

    #[test]
    fn class_predicates_split_the_structure() {
        assert_eq!(selected("protein"), vec![0, 1]);
        assert_eq!(selected("water"), vec![2]);
        assert_eq!(selected("metal"), vec![3]);
        assert_eq!(selected("hetatm"), vec![2, 3]);
        assert_eq!(selected("ligand"), Vec::<usize>::new()); // water and metal are not ligand
        assert_eq!(selected("all"), vec![0, 1, 2, 3]);
        assert_eq!(selected("none"), Vec::<usize>::new());
    }

    #[test]
    fn element_matches_whole_symbols_not_prefixes() {
        // A carbon query must not catch the calcium-like zinc, and vice versa.
        assert_eq!(selected("element C"), vec![1]);
        assert_eq!(selected("element Zn"), vec![3]);
        assert_eq!(selected("element zn"), vec![3]); // case-insensitive
    }

    #[test]
    fn boolean_operators_combine_sets() {
        assert_eq!(selected("protein and name CA"), vec![1]);
        assert_eq!(selected("water or metal"), vec![2, 3]);
        assert_eq!(selected("not protein"), vec![2, 3]);
        assert_eq!(selected("protein and not name CA"), vec![0]);
    }

    #[test]
    fn numeric_predicates_read_the_right_fields() {
        assert_eq!(selected("bfactor > 50"), vec![1]);
        assert_eq!(selected("bfactor <= 50"), vec![0, 2, 3]);
        assert_eq!(selected("plddt > 50"), vec![1]); // same column as bfactor
        assert_eq!(selected("ss H"), vec![1]);
        assert_eq!(selected("resi 1"), vec![0, 1]);
        assert_eq!(selected("model 1"), vec![0, 1, 2, 3]);
    }

    #[test]
    fn byres_widens_to_whole_residues() {
        // The alanine CA alone widens to both alanine atoms.
        assert_eq!(selected("byres name CA"), vec![0, 1]);
    }

    #[test]
    fn within_gathers_neighbours_and_includes_the_seed() {
        // Within 2 A of the water (at x=4): nothing else is that close, so only
        // the water itself.
        assert_eq!(selected("within 2 of water"), vec![2]);
        // Within 3 A of the alanine CA (x=1.5): the N at x=0, the CA itself, and
        // the water at x=4 (2.5 A away); the zinc at x=10 is too far.
        assert_eq!(selected("within 3 of name CA"), vec![0, 1, 2]);
    }

    #[test]
    fn an_empty_seed_within_is_empty() {
        assert!(run("within 5 of resn XXX", &sample()).unwrap().is_empty());
    }
}
