//! The confidence a co-folding model writes beside each structure it predicts:
//! AlphaFold 3's summary and Boltz's, found by the names those tools give them.
//!
//! AlphaFold 3 writes `X_model.cif` with `X_summary_confidences.json` (each
//! sample's `model.cif` with `summary_confidences.json` in older releases);
//! AlphaFold Server `X_model_N.cif` with `X_summary_confidences_N.json`;
//! Boltz `X_model_N.cif` with `confidence_X_model_N.json`.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

/// A model's confidence as its co-folding tool gives it.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Confidence {
    /// The tool that wrote it: `alphafold3` or `boltz`.
    pub source: &'static str,
    /// The score the tool ranks its samples by: AlphaFold 3's
    /// `ranking_score`, Boltz's `confidence_score`.
    pub score: f32,
    /// The confidence in the ligand's interface: AlphaFold 3's `chain_iptm`
    /// of the ligand's chain, Boltz's `ligand_iptm`.
    pub ligand_iptm: Option<f32>,
}

#[derive(Deserialize)]
struct AlphaFold {
    ranking_score: f32,
    #[serde(default)]
    chain_ids: Vec<String>,
    #[serde(default)]
    chain_iptm: Vec<Option<f32>>,
}

#[derive(Deserialize)]
struct Boltz {
    confidence_score: f32,
    ligand_iptm: Option<f32>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Source {
    AlphaFold,
    Boltz,
}

/// The confidence files a model file `pose` could have beside it, each with
/// the tool that would write it, in the order they are tried.
fn candidates(pose: &Path) -> Vec<(PathBuf, Source)> {
    let Some(name) = pose.file_name().and_then(|name| name.to_str()) else {
        return Vec::new();
    };
    let name = name.strip_suffix(".gz").unwrap_or(name);
    let stem = [".cif", ".mmcif", ".bcif", ".pdb"]
        .iter()
        .find_map(|extension| name.strip_suffix(extension))
        .unwrap_or(name);
    let beside = |file: String| pose.with_file_name(file);

    let mut found = Vec::new();
    if let Some(prefix) = stem.strip_suffix("model") {
        found.push((
            beside(format!("{prefix}summary_confidences.json")),
            Source::AlphaFold,
        ));
    }
    if let Some((base, number)) = stem.rsplit_once("_model_")
        && !number.is_empty()
        && number.chars().all(|c| c.is_ascii_digit())
    {
        found.push((
            beside(format!("{base}_summary_confidences_{number}.json")),
            Source::AlphaFold,
        ));
        found.push((beside(format!("confidence_{stem}.json")), Source::Boltz));
    }

    found
}

/// The confidence written beside model file `pose`, the ligand's interface
/// read for chain `chain`, or `None` when no such file is there. AlphaFold 3
/// lists each chain's ipTM in the order the model holds its chains,
/// `chains`, unless the summary names them (`chain_ids`).
pub fn beside(pose: &Path, chain: &str, chains: &[String]) -> Result<Option<Confidence>> {
    let Some((path, source)) = candidates(pose)
        .into_iter()
        .find(|(path, _)| path.is_file())
    else {
        return Ok(None);
    };
    let text = std::fs::read_to_string(&path)
        .with_context(|| format!("failed to read {}", path.display()))?;
    let confidence = match source {
        Source::AlphaFold => {
            let summary: AlphaFold = serde_json::from_str(&text)
                .with_context(|| format!("{} is no AlphaFold 3 summary", path.display()))?;
            let order = if summary.chain_ids.is_empty() {
                chains
            } else {
                &summary.chain_ids
            };
            let ligand_iptm = order
                .iter()
                .position(|id| id == chain)
                .and_then(|at| summary.chain_iptm.get(at).copied().flatten());
            Confidence {
                source: "alphafold3",
                score: summary.ranking_score,
                ligand_iptm,
            }
        }
        Source::Boltz => {
            let summary: Boltz = serde_json::from_str(&text)
                .with_context(|| format!("{} is no Boltz confidence file", path.display()))?;
            Confidence {
                source: "boltz",
                score: summary.confidence_score,
                ligand_iptm: summary.ligand_iptm,
            }
        }
    };

    Ok(Some(confidence))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(pose: &str) -> Vec<String> {
        candidates(Path::new(pose))
            .into_iter()
            .map(|(path, _)| path.to_string_lossy().into_owned())
            .collect()
    }

    #[test]
    fn each_tool_s_confidence_file_is_named_after_its_model() {
        assert_eq!(
            names("out/job_model.cif"),
            ["out/job_summary_confidences.json"]
        );
        assert_eq!(
            names("seed-1_sample-0/model.cif"),
            ["seed-1_sample-0/summary_confidences.json"]
        );
        assert_eq!(
            names("fold_job_model_3.cif"),
            [
                "fold_job_summary_confidences_3.json",
                "confidence_fold_job_model_3.json"
            ]
        );
        assert_eq!(
            names("job_model_0.pdb"),
            [
                "job_summary_confidences_0.json",
                "confidence_job_model_0.json"
            ]
        );
        assert!(names("crystal.cif").is_empty());
        assert!(names("job_model_x.cif").is_empty());
    }

    #[test]
    fn an_alphafold_summary_lists_its_chains_in_the_model_s_order() {
        let dir = tempfile::tempdir().unwrap();
        let pose = dir.path().join("job_model.cif");
        std::fs::write(
            dir.path().join("job_summary_confidences.json"),
            r#"{"ranking_score": 0.7, "chain_iptm": [0.9, 0.8, 0.4]}"#,
        )
        .unwrap();
        let chains = ["A", "B", "C"].map(String::from);

        let confidence = beside(&pose, "C", &chains).unwrap().unwrap();
        assert_eq!(confidence.ligand_iptm, Some(0.4));
        assert_eq!(
            beside(&pose, "Z", &chains).unwrap().unwrap().ligand_iptm,
            None
        );
    }

    #[test]
    fn a_malformed_file_beside_a_model_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let pose = dir.path().join("job_model_0.cif");
        std::fs::write(
            dir.path().join("confidence_job_model_0.json"),
            "{\"ptm\": 0.8}",
        )
        .unwrap();

        assert!(beside(&pose, "A", &[]).is_err());
        assert_eq!(
            beside(&dir.path().join("crystal.cif"), "A", &[]).unwrap(),
            None
        );
    }
}
