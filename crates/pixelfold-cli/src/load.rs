//! Getting from what the user typed to a loaded structure and a chosen subset.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use pixelfold_core::parser::load_protein_with_options;
use pixelfold_core::{AltlocPolicy, AtomSet, Protein, select};

use crate::resolve::{self, Resolution};

/// Resolve `input` to a file on disk, downloading it if it names a PDB entry
/// that is not cached yet.
pub fn structure_path(input: &str, cache_dir: &Path) -> Result<PathBuf> {
    match resolve::resolve(input, cache_dir)? {
        Resolution::Found(path) => Ok(path),
        Resolution::Fetch { id, dest } => {
            std::fs::create_dir_all(cache_dir)?;
            eprintln!("Fetching {id} from RCSB...");
            pixelfold_fetch::fetch_cif(&id, cache_dir)?;

            Ok(dest)
        }
    }
}

/// Load a structure for analysis, without the surface, i.e. nothing here renders.
pub fn analysis_structure(input: &str, cache_dir: &Path, altloc: AltlocPolicy) -> Result<Protein> {
    let path = structure_path(input, cache_dir)?;

    load_protein_with_options(&path, true, altloc)
        .with_context(|| format!("failed to load {}", path.display()))
}

/// The atoms a query chooses, or every atom when there is no query.
///
/// An empty result throws an error.
pub fn chosen(protein: &Protein, query: Option<&str>) -> Result<AtomSet> {
    let query = query.unwrap_or("all");

    let set = select(query, protein).map_err(|error| anyhow::anyhow!("{error}"))?;
    if set.is_empty() {
        anyhow::bail!("the selection `{query}` matched no atoms");
    }

    Ok(set)
}
