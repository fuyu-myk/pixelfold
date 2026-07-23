//! Regression gate: fail when an aggregate metric falls below (or, for error,
//! above) a committed threshold. Used by CI via `--check`.

use std::path::Path;

use anyhow::{Context, Result};
use serde::Deserialize;

use crate::report::Summary;

/// Per-metric floors and ceilings. A metric is only enforced when both its
/// threshold is set here and its golden data exists.
#[derive(Deserialize, Default)]
pub struct Thresholds {
    pub min_q3: Option<f64>,
    pub min_hbond_f1: Option<f64>,
    pub max_sasa_mae: Option<f64>,
}

/// Load thresholds from TOML, or defaults if absent.
pub fn load_thresholds(path: &Path) -> Result<Thresholds> {
    if !path.exists() {
        return Ok(Thresholds::default());
    }

    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading thresholds {}", path.display()))?;
    toml::from_str(&text).with_context(|| format!("parsing thresholds {}", path.display()))
}

/// A message for each present metric that violates its threshold. Pending
/// metrics are not checked.
pub fn evaluate(summary: &Summary, thresholds: &Thresholds) -> Vec<String> {
    let mut failures = Vec::new();

    if let (Some(min), Some(value)) = (thresholds.min_q3, summary.mean_q3)
        && value < min
    {
        failures.push(format!("Q3 {value:.3} below minimum {min:.3}"));
    }
    if let (Some(min), Some(value)) = (thresholds.min_hbond_f1, summary.mean_hbond_f1)
        && value < min
    {
        failures.push(format!("H-bond F1 {value:.3} below minimum {min:.3}"));
    }
    if let (Some(max), Some(value)) = (thresholds.max_sasa_mae, summary.mean_sasa_mae)
        && value > max
    {
        failures.push(format!("SASA MAE {value:.3} above maximum {max:.3}"));
    }

    failures
}

#[cfg(test)]
mod tests {
    use super::*;

    fn summary(q3: Option<f64>) -> Summary {
        Summary {
            entries: 1,
            dssp_entries: 1,
            mean_q3: q3,
            hbond_entries: 0,
            mean_hbond_f1: None,
            sasa_entries: 0,
            mean_sasa_mae: None,
            mean_sasa_median: None,
            mean_sasa_pearson: None,
            interactions: Vec::new(),
        }
    }

    #[test]
    fn below_threshold_fails_and_above_passes() {
        let thresholds = Thresholds {
            min_q3: Some(0.8),
            ..Default::default()
        };
        assert_eq!(evaluate(&summary(Some(0.7)), &thresholds).len(), 1);
        assert!(evaluate(&summary(Some(0.9)), &thresholds).is_empty());
    }

    #[test]
    fn pending_metric_is_not_enforced() {
        let thresholds = Thresholds {
            min_q3: Some(0.8),
            ..Default::default()
        };
        assert!(evaluate(&summary(None), &thresholds).is_empty());
    }
}
