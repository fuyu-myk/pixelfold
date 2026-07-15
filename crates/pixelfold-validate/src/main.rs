//! Pixelfold validation harness.
//!
//! For each entry in the benchmark manifest, it loads the structure, runs
//! pixelfold's analyses, aligns them against the golden reference outputs under
//! `benchmarks/golden/`, and prints a markdown summary. Entries without golden
//! files are reported as pending.

mod analysis;
mod compare;
mod golden;
mod metrics;
mod report;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use serde::Deserialize;

use pixelfold_core::AltlocPolicy;
use pixelfold_core::parser::load_protein_with_options;

use crate::compare::{EntryMetrics, compare};
use crate::golden::{DsspGolden, SasaGolden};
use crate::report::render_markdown;

/// Validate pixelfold against reference implementations (DSSP, FreeSASA).
#[derive(Parser)]
#[command(name = "pixelfold-validate", version)]
struct Cli {
    /// Benchmark manifest (TOML) listing PDB ids per stratum.
    #[arg(long, default_value = "benchmarks/manifest.toml")]
    manifest: PathBuf,
    /// Directory of golden reference files (`<ID>.dssp.json`, `<ID>.freesasa.json`).
    #[arg(long, default_value = "benchmarks/golden")]
    golden_dir: PathBuf,
    /// Directory holding and caching downloaded structures.
    #[arg(long, default_value = "benchmarks/cache")]
    cache_dir: PathBuf,
    /// Only evaluate structures already cached; never reach the network.
    #[arg(long)]
    offline: bool,
}

#[derive(Deserialize)]
struct Manifest {
    strata: BTreeMap<String, Vec<String>>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    std::fs::create_dir_all(&cli.cache_dir)?;

    let manifest_text = std::fs::read_to_string(&cli.manifest)
        .with_context(|| format!("reading manifest {}", cli.manifest.display()))?;
    let manifest: Manifest = toml::from_str(&manifest_text).context("parsing manifest")?;

    let ids: Vec<String> = manifest.strata.values().flatten().cloned().collect();
    if ids.is_empty() {
        println!(
            "Manifest has no entries yet; add PDB ids to {}.",
            cli.manifest.display()
        );
        return Ok(());
    }

    let mut metrics = Vec::new();
    for id in &ids {
        match evaluate(id, &cli.cache_dir, &cli.golden_dir, cli.offline) {
            Ok(Some(entry)) => metrics.push(entry),
            Ok(None) => eprintln!("skip {id}: not cached (offline)"),
            Err(e) => eprintln!("skip {id}: {e:#}"),
        }
    }

    print!("{}", render_markdown(&metrics));
    Ok(())
}

/// Load one structure, run the analyses, and compare against its golden files.
/// Returns `None` when the structure is absent and offline mode forbids fetching.
fn evaluate(
    id: &str,
    cache_dir: &Path,
    golden_dir: &Path,
    offline: bool,
) -> Result<Option<EntryMetrics>> {
    let id = id.to_uppercase();
    let path = cache_dir.join(format!("{id}.cif"));
    if !path.exists() {
        if offline {
            return Ok(None);
        }

        pixelfold_fetch::fetch_cif(&id, cache_dir).with_context(|| format!("fetching {id}"))?;
    }

    let protein = load_protein_with_options(&path, true, AltlocPolicy::default())
        .with_context(|| format!("loading {}", path.display()))?;
    let prediction = analysis::predict(&protein);

    let dssp: Option<DsspGolden> = golden::load(&golden_dir.join(format!("{id}.dssp.json")))?;
    let sasa: Option<SasaGolden> = golden::load(&golden_dir.join(format!("{id}.freesasa.json")))?;

    Ok(Some(compare(
        &id,
        &prediction,
        dssp.as_ref(),
        sasa.as_ref(),
    )))
}
