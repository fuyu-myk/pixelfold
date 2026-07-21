//! The selection language.
//!
//! A query string names a set of atoms, which is what makes the viewer, the
//! network filter, the CLI, and every export composable: each takes the same
//! selection. A query parses to a [`Selection`] tree and evaluates against a
//! [`Protein`] to an [`AtomSet`], a bitset over atom indices.
//!
//! ```text
//! selection := or_expr
//! or_expr   := and_expr ("or" and_expr)*
//! and_expr  := not_expr ("and" not_expr)*
//! not_expr  := "not" not_expr | primary
//! primary   := "(" selection ")"
//!            | "byres" primary
//!            | "within" FLOAT "of" primary
//!            | predicate
//! predicate := "chain" NAME+   | "resi" RANGE+   | "resn" NAME+
//!            | "name" NAME+     | "element" NAME+ | "model" INT
//!            | "ss" (H|E|C|T)   | "bfactor" CMP FLOAT
//!            | "plddt" CMP FLOAT
//!            | "protein" | "nucleic" | "ligand" | "water" | "metal" | "hetatm"
//!            | "all" | "none"
//! ```
//!
//! `byres` and `within` bind to a single primary, so a compound operand is
//! parenthesised: `within 5 of (chain A or chain B)`. Operators bind not tighter
//! than and tighter than or, and `not` is prefix. The keywords above are
//! reserved, so a chain or residue literally named one of them cannot be
//! selected bare; no structure in practice carries such a name.

mod eval;
mod parser;

use crate::structure::{Protein, SecondaryStructure};

pub use eval::AtomSet;

/// A parsed selection, evaluated against a structure by [`AtomSet::of`].
#[derive(Clone, Debug, PartialEq)]
pub enum Selection {
    All,
    None,
    Chain(Vec<String>),
    Resi(Vec<ResiRange>),
    Resn(Vec<String>),
    Name(Vec<String>),
    Element(Vec<String>),
    Model(usize),
    Ss(SecondaryStructure),
    BFactor(Cmp, f32),
    /// The predicted local distance difference, read from the same column as the
    /// B-factor, which is where a predicted structure stores it.
    Plddt(Cmp, f32),
    Protein,
    Nucleic,
    Ligand,
    Water,
    Metal,
    Hetatm,
    /// Every atom within the distance of any atom in the operand.
    Within(f32, Box<Selection>),
    /// Every atom of any residue with an atom in the operand.
    ByRes(Box<Selection>),
    Not(Box<Selection>),
    And(Box<Selection>, Box<Selection>),
    Or(Box<Selection>, Box<Selection>),
}

/// A numeric comparison for the `bfactor` and `plddt` predicates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cmp {
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
}

impl Cmp {
    /// Whether `value` satisfies the comparison against `bound`.
    ///
    /// Equality is exact: a user wanting a band around a value uses `<` and `>`,
    /// since a stored B-factor is rarely the round number typed.
    #[allow(clippy::float_cmp)]
    pub fn test(self, value: f32, bound: f32) -> bool {
        match self {
            Cmp::Lt => value < bound,
            Cmp::Le => value <= bound,
            Cmp::Gt => value > bound,
            Cmp::Ge => value >= bound,
            Cmp::Eq => value == bound,
        }
    }
}

/// An inclusive residue-number range; a single number has equal ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResiRange {
    pub lo: u32,
    pub hi: u32,
}

impl ResiRange {
    fn contains(self, seq: u32) -> bool {
        self.lo <= seq && seq <= self.hi
    }
}

/// A selection that could not be parsed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SelectError {
    pub message: String,
}

impl std::fmt::Display for SelectError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for SelectError {}

/// Parse `query` and evaluate it against `protein`.
pub fn select(query: &str, protein: &Protein) -> Result<AtomSet, SelectError> {
    Ok(AtomSet::of(&parser::parse(query)?, protein))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_range_includes_its_ends() {
        let range = ResiRange { lo: 5, hi: 10 };
        assert!(range.contains(5));
        assert!(range.contains(10));
        assert!(!range.contains(4));
        assert!(!range.contains(11));
    }

    #[test]
    fn comparisons_read_as_written() {
        assert!(Cmp::Lt.test(1.0, 2.0));
        assert!(!Cmp::Lt.test(2.0, 2.0));
        assert!(Cmp::Le.test(2.0, 2.0));
        assert!(Cmp::Ge.test(2.0, 2.0));
        assert!(Cmp::Eq.test(2.0, 2.0));
        assert!(!Cmp::Eq.test(2.0, 2.5));
    }

    #[test]
    fn a_bad_query_reports_rather_than_panics() {
        let protein = Protein {
            atoms: Vec::new(),
            title: String::new(),
            surface_points: Vec::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        };
        assert!(select("chain", &protein).is_err());
        assert!(select("", &protein).is_err());
        assert!(select("resn ALA and", &protein).is_err());
    }
}
