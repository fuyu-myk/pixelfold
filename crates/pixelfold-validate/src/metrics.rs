//! Comparison metrics against reference implementations. Pure functions over
//! already-aligned data; the harness handles the alignment.

use std::collections::HashSet;
use std::hash::Hash;

/// Three-state secondary structure, the target of Q3 agreement.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Ss3 {
    Helix,
    Strand,
    Coil,
}

/// Fraction of aligned positions where two secondary-structure labelings agree
/// (the Q3 score). `None` when the series are empty or differ in length.
pub fn q3_agreement(predicted: &[Ss3], reference: &[Ss3]) -> Option<f64> {
    if predicted.is_empty() || predicted.len() != reference.len() {
        return None;
    }

    let matches = predicted
        .iter()
        .zip(reference)
        .filter(|(a, b)| a == b)
        .count();

    Some(matches as f64 / predicted.len() as f64)
}

/// Precision, recall, and F1 for a predicted edge set against a reference set.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrF1 {
    pub precision: f64,
    pub recall: f64,
    pub f1: f64,
}

/// Precision / recall / F1 treating each set as a collection of true edges. An
/// empty predicted (or reference) set scores 1.0 on precision (or recall)
/// vacuously, so two empty sets agree perfectly.
pub fn set_prf1<T: Eq + Hash>(predicted: &HashSet<T>, reference: &HashSet<T>) -> PrF1 {
    let true_positives = predicted.intersection(reference).count() as f64;
    let precision = if predicted.is_empty() {
        1.0
    } else {
        true_positives / predicted.len() as f64
    };
    let recall = if reference.is_empty() {
        1.0
    } else {
        true_positives / reference.len() as f64
    };
    let f1 = if precision + recall == 0.0 {
        0.0
    } else {
        2.0 * precision * recall / (precision + recall)
    };

    PrF1 {
        precision,
        recall,
        f1,
    }
}

/// Error statistics between two aligned value series (for example per-residue
/// SASA in square Angstroms).
#[derive(Clone, Copy, Debug)]
pub struct ErrorStats {
    pub n: usize,
    pub mae: f64,
    pub median_abs: f64,
    pub pearson: f64,
}

/// Mean absolute error, median absolute error, and Pearson correlation between
/// `predicted` and `reference`. `None` when empty or mismatched in length.
pub fn error_stats(predicted: &[f64], reference: &[f64]) -> Option<ErrorStats> {
    if predicted.is_empty() || predicted.len() != reference.len() {
        return None;
    }

    let n = predicted.len();

    let mut abs_errors: Vec<f64> = predicted
        .iter()
        .zip(reference)
        .map(|(p, r)| (p - r).abs())
        .collect();
    let mae = abs_errors.iter().sum::<f64>() / n as f64;

    abs_errors.sort_by(f64::total_cmp);
    let median_abs = if n % 2 == 1 {
        abs_errors[n / 2]
    } else {
        (abs_errors[n / 2 - 1] + abs_errors[n / 2]) / 2.0
    };

    Some(ErrorStats {
        n,
        mae,
        median_abs,
        pearson: pearson_r(predicted, reference),
    })
}

/// Pearson correlation coefficient; 0.0 when either series has no variance.
fn pearson_r(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;

    let (mut cov, mut var_x, mut var_y) = (0.0, 0.0, 0.0);
    for (&xi, &yi) in x.iter().zip(y) {
        let (dx, dy) = (xi - mean_x, yi - mean_y);
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    let denom = (var_x * var_y).sqrt();
    if denom == 0.0 { 0.0 } else { cov / denom }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn q3_scores_matches_over_length() {
        use Ss3::*;
        assert_eq!(
            q3_agreement(&[Helix, Helix, Coil], &[Helix, Helix, Coil]),
            Some(1.0)
        );
        assert_eq!(q3_agreement(&[Helix, Strand], &[Helix, Coil]), Some(0.5));
        assert_eq!(q3_agreement(&[], &[]), None);
        assert_eq!(q3_agreement(&[Helix], &[Helix, Coil]), None);
    }

    #[test]
    fn prf1_on_partial_overlap() {
        let predicted: HashSet<i32> = [1, 2, 3].into_iter().collect();
        let reference: HashSet<i32> = [2, 3, 4].into_iter().collect();
        let r = set_prf1(&predicted, &reference);
        assert!((r.precision - 2.0 / 3.0).abs() < 1e-9);
        assert!((r.recall - 2.0 / 3.0).abs() < 1e-9);
        assert!((r.f1 - 2.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn prf1_identical_sets_are_perfect_and_empty_sets_agree() {
        let s: HashSet<i32> = [1, 2].into_iter().collect();
        let perfect = set_prf1(&s, &s);
        assert_eq!(
            (perfect.precision, perfect.recall, perfect.f1),
            (1.0, 1.0, 1.0)
        );

        let empty: HashSet<i32> = HashSet::new();
        assert_eq!(set_prf1(&empty, &empty).f1, 1.0);
    }

    #[test]
    fn error_stats_on_identical_and_offset_series() {
        let a = [1.0, 2.0, 3.0, 4.0];
        let identical = error_stats(&a, &a).unwrap();
        assert_eq!(identical.mae, 0.0);
        assert_eq!(identical.median_abs, 0.0);
        assert!((identical.pearson - 1.0).abs() < 1e-9);

        // A constant +1 offset: perfect correlation, unit error.
        let b = [2.0, 3.0, 4.0, 5.0];
        let offset = error_stats(&a, &b).unwrap();
        assert!((offset.mae - 1.0).abs() < 1e-9);
        assert!((offset.median_abs - 1.0).abs() < 1e-9);
        assert!((offset.pearson - 1.0).abs() < 1e-9);
    }

    #[test]
    fn error_stats_rejects_mismatched_input() {
        assert!(error_stats(&[1.0], &[1.0, 2.0]).is_none());
        assert!(error_stats(&[], &[]).is_none());
    }
}
