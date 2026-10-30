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
/// pixelfold's, so the comparison here is a plain set intersection per type,
/// over the part of the structure the reference looked at.
#[derive(Clone, Debug, Deserialize)]
pub struct InteractionGolden {
    /// The interaction kinds this reference reports, by pixelfold's label.
    pub kinds: Vec<String>,
    /// Whether the reference only covers interactions that touch a ligand, as
    /// PLIP does by profiling binding sites.
    #[serde(default)]
    pub ligand_scoped: bool,
    /// The residues a ligand-scoped reference profiled as ligands. Only these
    /// count: the reference decides what a ligand is (PLIP skips modified
    /// residues and some crystallisation additives), not the HETATM record.
    #[serde(default)]
    pub ligands: Vec<ResidueKey>,
    /// Whether the reference left water out of the structure, as RING does by
    /// default, so it never reports an interaction with a water.
    #[serde(default)]
    pub water_excluded: bool,
    /// Standard amino acids of one chain fewer than this many positions apart
    /// are never paired by the reference (RING's sequence gap, 3 by default);
    /// other residues pair at any distance. Zero pairs any two.
    #[serde(default)]
    pub sequence_separation: u32,
    #[serde(default)]
    pub interactions: Vec<InteractionEdge>,
}

impl InteractionGolden {
    /// Errors when a ligand-scoped reference has an interaction touching none of
    /// its listed ligands: a stale or broken golden that would otherwise
    /// silently narrow what is scored.
    pub fn check_scope(&self) -> Result<()> {
        if !self.ligand_scoped {
            return Ok(());
        }

        let unlisted = self
            .interactions
            .iter()
            .find(|edge| !self.ligands.contains(&edge.a) && !self.ligands.contains(&edge.b));
        match unlisted {
            Some(edge) => anyhow::bail!(
                "{} interaction {}-{} touches no listed ligand; regenerate the golden file",
                edge.kind,
                edge.a,
                edge.b
            ),
            None => Ok(()),
        }
    }
}

impl std::fmt::Display for ResidueKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}/{}{}",
            self.chain,
            self.seq,
            self.icode.as_deref().unwrap_or("")
        )
    }
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

/// Every golden reference for one entry, each `None` when its tool has no file.
pub struct Goldens {
    pub dssp: Option<DsspGolden>,
    pub sasa: Option<SasaGolden>,
    pub plip: Option<InteractionGolden>,
    pub ring: Option<InteractionGolden>,
}

impl Goldens {
    pub fn load(golden_dir: &Path, id: &str) -> Result<Self> {
        let path = |tool: &str| golden_dir.join(format!("{id}.{tool}.json"));

        Ok(Self {
            dssp: load(&path("dssp"))?,
            sasa: load(&path("freesasa"))?,
            plip: load_interactions(&path("plip"))?,
            ring: load_interactions(&path("ring"))?,
        })
    }
}

/// [`load`] for an interaction golden, rejecting one whose scope is
/// inconsistent (see [`InteractionGolden::check_scope`]).
pub fn load_interactions(path: &Path) -> Result<Option<InteractionGolden>> {
    let golden: Option<InteractionGolden> = load(path)?;
    if let Some(golden) = &golden {
        golden
            .check_scope()
            .with_context(|| format!("checking golden file {}", path.display()))?;
    }

    Ok(golden)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn golden(json: &str) -> InteractionGolden {
        serde_json::from_str(json).expect("test golden is valid JSON")
    }

    const EDGE: &str = r#"{"kind": "hydrogen-bond", "a": {"chain": "A", "seq": 1}, "b": {"chain": "A", "seq": 9}}"#;

    // A ligand-scoped golden written before ligands were listed, or by a
    // generator that dropped one, would silently narrow what is scored.
    #[test]
    fn a_ligand_scoped_golden_must_list_the_ligand_of_every_interaction() {
        let listed = golden(&format!(
            r#"{{"kinds": ["hydrogen-bond"], "ligand_scoped": true,
                "ligands": [{{"chain": "A", "seq": 9}}], "interactions": [{EDGE}]}}"#
        ));
        assert!(listed.check_scope().is_ok());

        let unlisted = golden(&format!(
            r#"{{"kinds": ["hydrogen-bond"], "ligand_scoped": true, "interactions": [{EDGE}]}}"#
        ));
        let error = unlisted.check_scope().expect_err("ligand 9 is not listed");
        assert!(
            error.to_string().contains("A/1-A/9"),
            "the error names the edge: {error}"
        );
    }

    #[test]
    fn a_whole_structure_golden_lists_no_ligands() {
        let whole = golden(&format!(
            r#"{{"kinds": ["hydrogen-bond"], "ligand_scoped": false, "interactions": [{EDGE}]}}"#
        ));
        assert!(whole.check_scope().is_ok());
    }
}
