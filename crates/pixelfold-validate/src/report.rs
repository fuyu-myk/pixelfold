//! Aggregate per-entry metrics into a markdown summary.

use std::collections::BTreeMap;

use crate::compare::{EntryMetrics, Tool};

/// Mean interaction precision/recall/F1 for one tool and one type, over the
/// entries that carry that tool's golden and report the type.
pub struct InteractionSummary {
    pub tool: &'static str,
    pub kind: String,
    pub entries: usize,
    /// Total predicted and reference edges of this type summed over the entries,
    /// the scale the averaged precision and recall are measured on.
    pub pred_edges: usize,
    pub ref_edges: usize,
    pub precision: f64,
    pub recall: f64,
    pub f1: f64,
}

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
    /// Per-tool, per-type interaction agreement, in a stable order.
    pub interactions: Vec<InteractionSummary>,
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
        interactions: interaction_summary(entries),
    }
}

/// Running sums for one interaction type over the scored entries.
#[derive(Default)]
struct KindTotals {
    precision: f64,
    recall: f64,
    f1: f64,
    pred_edges: usize,
    ref_edges: usize,
    entries: usize,
}

/// Average each interaction type's precision/recall/F1 over the entries scored
/// against each tool. A type is macro-averaged: every scored entry counts once,
/// so a large complex does not dominate a small one.
fn interaction_summary(entries: &[EntryMetrics]) -> Vec<InteractionSummary> {
    let mut out = Vec::new();

    for tool in [Tool::Plip, Tool::Ring] {
        let per_entry = entries.iter().filter_map(|e| match tool {
            Tool::Plip => e.plip.as_ref(),
            Tool::Ring => e.ring.as_ref(),
        });

        // Kept in the type's first-seen order across the entries.
        let mut order: Vec<String> = Vec::new();
        let mut sums: BTreeMap<String, KindTotals> = BTreeMap::new();
        for metrics in per_entry {
            for metric in metrics {
                let slot = sums.entry(metric.kind.clone()).or_insert_with(|| {
                    order.push(metric.kind.clone());
                    KindTotals::default()
                });
                slot.precision += metric.prf1.precision;
                slot.recall += metric.prf1.recall;
                slot.f1 += metric.prf1.f1;
                slot.pred_edges += metric.predicted;
                slot.ref_edges += metric.reference;
                slot.entries += 1;
            }
        }

        for kind in order {
            let totals = &sums[&kind];
            let n = totals.entries as f64;
            out.push(InteractionSummary {
                tool: tool.label(),
                kind,
                entries: totals.entries,
                pred_edges: totals.pred_edges,
                ref_edges: totals.ref_edges,
                precision: totals.precision / n,
                recall: totals.recall / n,
                f1: totals.f1 / n,
            });
        }
    }

    out
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

    if !summary.interactions.is_empty() {
        out.push_str("\nPer-type interaction agreement (precision / recall / F1):\n\n");
        out.push_str(
            "| Interaction type | Reference | Entries | Pred edges | Ref edges | Precision | Recall | F1 |\n",
        );
        out.push_str("| --- | --- | --- | --- | --- | --- | --- | --- |\n");
        for row in &summary.interactions {
            out.push_str(&format!(
                "| {} | {} | {} | {} | {} | {:.2} | {:.2} | {:.2} |\n",
                row.kind,
                row.tool,
                row.entries,
                row.pred_edges,
                row.ref_edges,
                row.precision,
                row.recall,
                row.f1,
            ));
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare::KindMetric;
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
            plip: None,
            ring: None,
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

    fn kind_metric(kind: &str, f1: f64, reference: usize) -> KindMetric {
        KindMetric {
            kind: kind.to_string(),
            prf1: PrF1 {
                precision: f1,
                recall: f1,
                f1,
            },
            predicted: reference,
            reference,
        }
    }

    fn interaction_entry(id: &str, plip: Vec<KindMetric>) -> EntryMetrics {
        EntryMetrics {
            id: id.to_string(),
            q3: None,
            hbond_f1: None,
            sasa: None,
            aligned_residues: 0,
            plip: Some(plip),
            ring: None,
        }
    }

    /// Each type is macro-averaged over the entries that report it, and only the
    /// tool with a golden appears.
    #[test]
    fn interaction_types_macro_average_over_entries() {
        let entries = vec![
            interaction_entry("a", vec![kind_metric("hydrogen-bond", 0.8, 10)]),
            interaction_entry("b", vec![kind_metric("hydrogen-bond", 1.0, 30)]),
        ];
        let summary = summarize(&entries);

        assert_eq!(summary.interactions.len(), 1); // PLIP only, one type
        let hbond = &summary.interactions[0];
        assert_eq!(hbond.tool, "PLIP");
        assert_eq!(hbond.entries, 2);
        assert_eq!(hbond.ref_edges, 40);
        assert!((hbond.f1 - 0.9).abs() < 1e-9); // (0.8 + 1.0) / 2, not edge-weighted

        let md = render_markdown(&entries);
        assert!(md.contains("Per-type interaction agreement"), "{md}");
        assert!(md.contains("| hydrogen-bond | PLIP |"), "{md}");
    }
}
