//! Pixelfold validation harness.
//!
//! For each entry in the benchmark manifest, it loads the structure, runs
//! pixelfold's analyses, aligns them against the golden reference outputs under
//! `benchmarks/golden/`, and prints a markdown summary. Entries without golden
//! files are reported as pending.

mod analysis;
mod check;
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
use crate::golden::Goldens;
use crate::report::{render_markdown, summarize};

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
    /// Fail (exit 1) if an aggregate metric violates the thresholds file.
    #[arg(long)]
    check: bool,
    /// Regression thresholds for `--check`.
    #[arg(long, default_value = "benchmarks/thresholds.toml")]
    thresholds: PathBuf,
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
        let id = id.to_uppercase();
        // A golden file that fails to parse or contradicts itself is a defect in
        // the repository, not a missing input: skipping its entry would quietly
        // switch off the gates it feeds.
        let goldens = Goldens::load(&cli.golden_dir, &id)?;
        match evaluate(&id, &cli.cache_dir, cli.offline, &goldens) {
            Ok(Some(entry)) => metrics.push(entry),
            Ok(None) => eprintln!("skip {id}: not cached (offline)"),
            Err(e) => eprintln!("skip {id}: {e:#}"),
        }
    }

    print!("{}", render_markdown(&metrics));

    if cli.check {
        let thresholds = check::load_thresholds(&cli.thresholds)?;
        let failures = check::evaluate(&summarize(&metrics), &thresholds);
        if !failures.is_empty() {
            eprintln!("\nRegression gate failed:");
            for failure in &failures {
                eprintln!("  - {failure}");
            }

            std::process::exit(1);
        }
    }

    Ok(())
}

/// Load one structure, run the analyses, and compare against its golden files.
/// Returns `None` when the structure is absent and offline mode forbids fetching.
fn evaluate(
    id: &str,
    cache_dir: &Path,
    offline: bool,
    goldens: &Goldens,
) -> Result<Option<EntryMetrics>> {
    let path = cache_dir.join(format!("{id}.cif"));
    if !path.exists() {
        if offline {
            return Ok(None);
        }

        pixelfold_fetch::fetch_cif(id, cache_dir).with_context(|| format!("fetching {id}"))?;
    }

    let protein = load_protein_with_options(&path, AltlocPolicy::default())
        .with_context(|| format!("loading {}", path.display()))?;
    let prediction = analysis::predict(&protein);

    Ok(Some(compare(
        id,
        &prediction,
        goldens.dssp.as_ref(),
        goldens.sasa.as_ref(),
        goldens.plip.as_ref(),
        goldens.ring.as_ref(),
    )))
}
