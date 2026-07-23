//! Golden reference outputs from the external tools, keyed by residue.
//!
//! These are the JSON files under `benchmarks/golden/`, one per entry per tool,
//! produced by the dockerised reference tools. The harness deserialises them and
//! aligns each residue against pixelfold's own output by [`ResidueKey`].

use std::path::Path;

use anyhow::{Context, Result};
use serde::Deserialize;

/// A residue identity shared by predictions and golden references.
///
/// Ordered so an interaction's two residues can be put in a canonical order,
/// collapsing the two directions of one edge to a single unordered pair.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Deserialize)]
pub struct ResidueKey {
    pub chain: String,
    pub seq: i64,
    #[serde(default)]
    pub icode: Option<String>,
}

/// DSSP (`mkdssp`) reference: per-residue secondary structure and backbone
/// hydrogen bonds.
#[derive(Clone, Debug, Deserialize)]
pub struct DsspGolden {
    pub residues: Vec<DsspResidue>,
    /// `None` when the reference omits hydrogen bonds; `Some([])` means it
    /// genuinely found none. The distinction decides whether the F1 is scored.
    #[serde(default)]
    pub hbonds: Option<Vec<HBondEdge>>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct DsspResidue {
    #[serde(flatten)]
    pub key: ResidueKey,
    /// The DSSP single-letter code (`H G I E B T S`, blank for coil).
    pub ss: String,
}

/// A backbone hydrogen bond from the N-H of `donor` to the C=O of `acceptor`.
#[derive(Clone, Debug, Deserialize)]
pub struct HBondEdge {
    pub donor: ResidueKey,
    pub acceptor: ResidueKey,
}

/// PLIP or RING reference: residue-residue non-covalent interactions by type.
///
/// The generator normalises each tool's own type names and residue numbering to
/// pixelfold's, so the comparison here is a plain set intersection per type.
#[derive(Clone, Debug, Deserialize)]
pub struct InteractionGolden {
    /// The interaction kinds this reference reports, by pixelfold's label.
    pub kinds: Vec<String>,
    #[serde(default)]
    pub interactions: Vec<InteractionEdge>,
}

/// One residue-residue interaction of a given kind. The pair is unordered; the
/// loader puts it in a canonical order before comparison.
#[derive(Clone, Debug, Deserialize)]
pub struct InteractionEdge {
    pub kind: String,
    pub a: ResidueKey,
    pub b: ResidueKey,
}

/// FreeSASA reference: per-residue solvent-accessible surface area (A^2).
#[derive(Clone, Debug, Deserialize)]
pub struct SasaGolden {
    pub residues: Vec<SasaResidue>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct SasaResidue {
    #[serde(flatten)]
    pub key: ResidueKey,
    pub sasa: f64,
}

/// Load and parse a JSON golden file, returning `Ok(None)` if it does not exist
/// (the entry simply has no reference for that tool yet).
pub fn load<T: for<'de> Deserialize<'de>>(path: &Path) -> Result<Option<T>> {
    if !path.exists() {
        return Ok(None);
    }

    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading golden file {}", path.display()))?;
    let value = serde_json::from_str(&text)
        .with_context(|| format!("parsing golden file {}", path.display()))?;

    Ok(Some(value))
}
