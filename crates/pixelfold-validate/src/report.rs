//! Aggregate per-entry metrics into a markdown summary.

use crate::compare::EntryMetrics;

/// Means across the entries that carry each kind of golden reference.
pub struct Summary {
    pub entries: usize,
    pub dssp_entries: usize,
    pub mean_q3: Option<f64>,
    pub hbond_entries: usize,
    pub mean_hbond_f1: Option<f64>,
    pub sasa_entries: usize,
    pub mean_sasa_mae: Option<f64>,
    pub mean_sasa_median: Option<f64>,
    pub mean_sasa_pearson: Option<f64>,
}

/// Average each metric over the entries that reported it.
pub fn summarize(entries: &[EntryMetrics]) -> Summary {
    let q3: Vec<f64> = entries.iter().filter_map(|e| e.q3).collect();
    let f1: Vec<f64> = entries
        .iter()
        .filter_map(|e| e.hbond_f1.map(|p| p.f1))
        .collect();
    let mae: Vec<f64> = entries
        .iter()
        .filter_map(|e| e.sasa.map(|s| s.mae))
        .collect();
    let median: Vec<f64> = entries
        .iter()
        .filter_map(|e| e.sasa.map(|s| s.median_abs))
        .collect();
    let pearson: Vec<f64> = entries
        .iter()
        .filter_map(|e| e.sasa.map(|s| s.pearson))
        .collect();

    Summary {
        entries: entries.len(),
        dssp_entries: q3.len(),
        mean_q3: mean(&q3),
        hbond_entries: f1.len(),
        mean_hbond_f1: mean(&f1),
        sasa_entries: mae.len(),
        mean_sasa_mae: mean(&mae),
        mean_sasa_median: mean(&median),
        mean_sasa_pearson: mean(&pearson),
    }
}

fn mean(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        None
    } else {
        Some(values.iter().sum::<f64>() / values.len() as f64)
    }
}

/// Render a per-entry breakdown (for spotting disagreements) followed by the
/// aggregate summary, as markdown for the terminal and the README.
pub fn render_markdown(entries: &[EntryMetrics]) -> String {
    let pct = |v: Option<f64>| v.map_or("pending".to_string(), |x| format!("{:.1}%", x * 100.0));
    let num =
        |v: Option<f64>, unit: &str| v.map_or("pending".to_string(), |x| format!("{x:.2}{unit}"));

    let mut out = String::new();

    let scored: Vec<&EntryMetrics> = entries
        .iter()
        .filter(|e| e.q3.is_some() || e.sasa.is_some())
        .collect();
    if !scored.is_empty() {
        out.push_str("| Entry | Residues | Q3 | H-bond F1 | SASA MAE | SASA n |\n");
        out.push_str("| --- | --- | --- | --- | --- | --- |\n");
        for e in scored {
            out.push_str(&format!(
                "| {} | {} | {} | {} | {} | {} |\n",
                e.id,
                e.aligned_residues,
                pct(e.q3),
                num(e.hbond_f1.map(|p| p.f1), ""),
                num(e.sasa.map(|s| s.mae), ""),
                e.sasa.map_or("-".to_string(), |s| s.n.to_string()),
            ));
        }

        out.push('\n');
    }

    let summary = summarize(entries);
    out.push_str(&format!(
        "Validated {} benchmark entries against DSSP 4 and FreeSASA.\n\n",
        summary.entries
    ));
    out.push_str("| Metric | Reference | Entries | Result |\n");
    out.push_str("| --- | --- | --- | --- |\n");
    out.push_str(&format!(
        "| Secondary structure (Q3) | DSSP 4 | {} | {} |\n",
        summary.dssp_entries,
        pct(summary.mean_q3),
    ));
    out.push_str(&format!(
        "| Backbone H-bond edges (F1) | DSSP 4 | {} | {} |\n",
        summary.hbond_entries,
        num(summary.mean_hbond_f1, ""),
    ));
    out.push_str(&format!(
        "| SASA mean abs error | FreeSASA | {} | {} |\n",
        summary.sasa_entries,
        num(summary.mean_sasa_mae, " A^2"),
    ));
    out.push_str(&format!(
        "| SASA median abs error | FreeSASA | {} | {} |\n",
        summary.sasa_entries,
        num(summary.mean_sasa_median, " A^2"),
    ));
    out.push_str(&format!(
        "| SASA correlation (Pearson r) | FreeSASA | {} | {} |\n",
        summary.sasa_entries,
        num(summary.mean_sasa_pearson, ""),
    ));

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::{ErrorStats, PrF1};

    fn entry(id: &str, q3: Option<f64>) -> EntryMetrics {
        EntryMetrics {
            id: id.to_string(),
            q3,
            hbond_f1: q3.map(|_| PrF1 {
                precision: 1.0,
                recall: 1.0,
                f1: 0.8,
            }),
            sasa: Some(ErrorStats {
                n: 10,
                mae: 2.0,
                median_abs: 1.0,
                pearson: 0.95,
            }),
            aligned_residues: 10,
        }
    }

    #[test]
    fn averages_only_over_entries_with_golden() {
        let entries = vec![
            entry("a", Some(0.8)),
            entry("b", Some(0.9)),
            entry("c", None),
        ];
        let s = summarize(&entries);
        assert_eq!(s.entries, 3);
        assert_eq!(s.dssp_entries, 2); // only a, b had Q3
        assert!((s.mean_q3.unwrap() - 0.85).abs() < 1e-9);
        assert_eq!(s.sasa_entries, 3);
        assert!((s.mean_sasa_mae.unwrap() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn missing_metrics_render_as_pending() {
        let md = render_markdown(&[]);
        assert!(md.contains("pending"), "{md}");
    }
}
