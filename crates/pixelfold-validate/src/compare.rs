//! Compare a prediction against golden references, aligning by residue.

use std::collections::{HashMap, HashSet};

use crate::analysis::{Prediction, ordered_pair, reduce_dssp_letter};
use crate::golden::{DsspGolden, InteractionGolden, ResidueKey, SasaGolden};
use crate::metrics::{ErrorStats, PrF1, Ss3, error_stats, q3_agreement, set_prf1};

/// Which reference tool an interaction comparison is against.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tool {
    Plip,
    Ring,
}

impl Tool {
    pub fn label(self) -> &'static str {
        match self {
            Tool::Plip => "PLIP",
            Tool::Ring => "RING",
        }
    }
}

/// Precision, recall, and F1 for one interaction type against one reference,
/// with the edge counts the score was computed over.
pub struct KindMetric {
    pub kind: String,
    pub prf1: PrF1,
    pub predicted: usize,
    pub reference: usize,
}

/// Metrics for one benchmark entry; each field is `None` when its golden is
/// absent (or, for H-bonds, when the reference does not report them).
pub struct EntryMetrics {
    pub id: String,
    pub q3: Option<f64>,
    pub hbond_f1: Option<PrF1>,
    pub sasa: Option<ErrorStats>,
    pub aligned_residues: usize,
    /// Per-type interaction metrics against PLIP, absent when no PLIP golden.
    pub plip: Option<Vec<KindMetric>>,
    /// Per-type interaction metrics against RING, absent when no RING golden.
    pub ring: Option<Vec<KindMetric>>,
}

/// Compare a prediction against whichever golden references are available.
pub fn compare(
    id: &str,
    prediction: &Prediction,
    dssp: Option<&DsspGolden>,
    sasa: Option<&SasaGolden>,
    plip: Option<&InteractionGolden>,
    ring: Option<&InteractionGolden>,
) -> EntryMetrics {
    let (q3, hbond_f1, aligned_residues) = match dssp {
        Some(dssp) => compare_dssp(prediction, dssp),
        None => (None, None, 0),
    };

    EntryMetrics {
        id: id.to_string(),
        q3,
        hbond_f1,
        sasa: sasa.and_then(|s| compare_sasa(prediction, s)),
        aligned_residues,
        plip: plip.map(|g| compare_interactions(prediction, g)),
        ring: ring.map(|g| compare_interactions(prediction, g)),
    }
}

/// Per-type precision/recall/F1 of the prediction against an interaction golden.
///
/// Only the kinds the reference reports are scored. A kind the reference reports
/// but found none of still scores, against whatever pixelfold found; a kind the
/// reference does not report at all is left out.
///
/// A ligand-scoped reference (PLIP) covers only interactions that touch a ligand,
/// so pixelfold's predicted set is restricted to its own ligand-touching edges
/// first. Without that, every protein-protein interaction the reference never
/// looked at would count against precision.
fn compare_interactions(prediction: &Prediction, golden: &InteractionGolden) -> Vec<KindMetric> {
    let mut reference: HashMap<&str, HashSet<(ResidueKey, ResidueKey)>> = HashMap::new();
    for edge in &golden.interactions {
        reference
            .entry(edge.kind.as_str())
            .or_default()
            .insert(ordered_pair(edge.a.clone(), edge.b.clone()));
    }

    let touches_ligand = |pair: &(ResidueKey, ResidueKey)| {
        prediction.ligands.contains(&pair.0) || prediction.ligands.contains(&pair.1)
    };

    let empty = HashSet::new();
    golden
        .kinds
        .iter()
        .map(|kind| {
            let all = prediction.interactions.get(kind).unwrap_or(&empty);
            let predicted: HashSet<(ResidueKey, ResidueKey)> = if golden.ligand_scoped {
                all.iter().filter(|p| touches_ligand(p)).cloned().collect()
            } else {
                all.clone()
            };
            let found = reference.get(kind.as_str()).unwrap_or(&empty);

            KindMetric {
                kind: kind.clone(),
                prf1: set_prf1(&predicted, found),
                predicted: predicted.len(),
                reference: found.len(),
            }
        })
        .collect()
}

fn compare_dssp(prediction: &Prediction, dssp: &DsspGolden) -> (Option<f64>, Option<PrF1>, usize) {
    let predicted: HashMap<&ResidueKey, Ss3> = prediction.ss.iter().map(|(k, s)| (k, *s)).collect();

    let mut pred_ss = Vec::new();
    let mut ref_ss = Vec::new();
    for residue in &dssp.residues {
        if let Some(&p) = predicted.get(&residue.key) {
            pred_ss.push(p);
            ref_ss.push(reduce_dssp_letter(&residue.ss));
        }
    }
    let aligned = pred_ss.len();
    let q3 = q3_agreement(&pred_ss, &ref_ss);

    let hbond_f1 = dssp.hbonds.as_ref().map(|edges| {
        let reference: HashSet<(ResidueKey, ResidueKey)> = edges
            .iter()
            .map(|e| (e.donor.clone(), e.acceptor.clone()))
            .collect();
        set_prf1(&prediction.hbonds, &reference)
    });

    (q3, hbond_f1, aligned)
}

fn compare_sasa(prediction: &Prediction, sasa: &SasaGolden) -> Option<ErrorStats> {
    let predicted: HashMap<&ResidueKey, f64> =
        prediction.sasa.iter().map(|(k, v)| (k, *v)).collect();

    let mut pred_values = Vec::new();
    let mut ref_values = Vec::new();
    for residue in &sasa.residues {
        if let Some(&p) = predicted.get(&residue.key) {
            pred_values.push(p);
            ref_values.push(residue.sasa);
        }
    }

    error_stats(&pred_values, &ref_values)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::golden::{DsspResidue, InteractionEdge, SasaResidue};

    fn key(seq: i64) -> ResidueKey {
        ResidueKey {
            chain: "A".to_string(),
            seq,
            icode: None,
        }
    }

    fn prediction() -> Prediction {
        Prediction {
            ss: vec![(key(1), Ss3::Helix)],
            hbonds: HashSet::new(),
            sasa: vec![],
            interactions: HashMap::new(),
            ligands: HashSet::new(),
        }
    }

    // A prediction compared against a golden built from its own output must
    // score perfectly; this locks the alignment and metric plumbing.
    #[test]
    fn self_comparison_is_perfect() {
        let prediction = Prediction {
            ss: vec![
                (key(1), Ss3::Helix),
                (key(2), Ss3::Helix),
                (key(3), Ss3::Coil),
            ],
            hbonds: [(key(1), key(3))].into_iter().collect(),
            sasa: vec![(key(1), 10.0), (key(2), 20.0), (key(3), 30.0)],
            interactions: HashMap::new(),
            ligands: HashSet::new(),
        };

        let dssp = DsspGolden {
            residues: vec![
                DsspResidue {
                    key: key(1),
                    ss: "H".into(),
                },
                DsspResidue {
                    key: key(2),
                    ss: "H".into(),
                },
                DsspResidue {
                    key: key(3),
                    ss: "C".into(),
                },
            ],
            hbonds: Some(vec![crate::golden::HBondEdge {
                donor: key(1),
                acceptor: key(3),
            }]),
        };
        let sasa = SasaGolden {
            residues: vec![
                SasaResidue {
                    key: key(1),
                    sasa: 10.0,
                },
                SasaResidue {
                    key: key(2),
                    sasa: 20.0,
                },
                SasaResidue {
                    key: key(3),
                    sasa: 30.0,
                },
            ],
        };

        let m = compare("self", &prediction, Some(&dssp), Some(&sasa), None, None);
        assert_eq!(m.q3, Some(1.0));
        assert_eq!(m.aligned_residues, 3);
        assert_eq!(m.hbond_f1.unwrap().f1, 1.0);
        let sasa = m.sasa.unwrap();
        assert_eq!(sasa.mae, 0.0);
        assert!((sasa.pearson - 1.0).abs() < 1e-9);
        assert!(m.plip.is_none() && m.ring.is_none());
    }

    #[test]
    fn residues_missing_from_the_prediction_are_skipped() {
        let dssp = DsspGolden {
            residues: vec![
                DsspResidue {
                    key: key(1),
                    ss: "H".into(),
                },
                DsspResidue {
                    key: key(2),
                    ss: "E".into(),
                }, // not in the prediction
            ],
            hbonds: None,
        };
        let m = compare("x", &prediction(), Some(&dssp), None, None, None);
        assert_eq!(m.aligned_residues, 1); // only residue 1 aligned
        assert_eq!(m.q3, Some(1.0));
        assert!(m.hbond_f1.is_none()); // reference omitted H-bonds
    }

    fn edge(kind: &str, a: i64, b: i64) -> InteractionEdge {
        InteractionEdge {
            kind: kind.to_string(),
            a: key(a),
            b: key(b),
        }
    }

    /// A prediction and a golden that agree on one type and disagree on another,
    /// with the two directions of one edge treated as the same edge.
    #[test]
    fn interactions_score_per_type_and_ignore_direction() {
        let mut interactions = HashMap::new();
        interactions.insert(
            "hydrogen-bond".to_string(),
            [ordered_pair(key(1), key(5)), ordered_pair(key(2), key(6))]
                .into_iter()
                .collect(),
        );
        interactions.insert(
            "salt-bridge".to_string(),
            [ordered_pair(key(1), key(9))].into_iter().collect(),
        );
        let prediction = Prediction {
            ss: vec![],
            hbonds: HashSet::new(),
            sasa: vec![],
            interactions,
            ligands: HashSet::new(),
        };

        let golden = InteractionGolden {
            kinds: vec![
                "hydrogen-bond".to_string(),
                "salt-bridge".to_string(),
                "pi-stacking".to_string(),
            ],
            ligand_scoped: false,
            interactions: vec![
                // Same two hydrogen bonds, one written in the reverse direction.
                edge("hydrogen-bond", 5, 1),
                edge("hydrogen-bond", 2, 6),
                // A salt bridge the reference has and pixelfold missed, plus one
                // they share.
                edge("salt-bridge", 1, 9),
                edge("salt-bridge", 3, 8),
            ],
        };

        let metrics = compare_interactions(&prediction, &golden);
        let by_kind: HashMap<&str, &KindMetric> =
            metrics.iter().map(|m| (m.kind.as_str(), m)).collect();

        // Hydrogen bonds match exactly despite the reversed direction.
        assert_eq!(by_kind["hydrogen-bond"].prf1.f1, 1.0);

        // Salt bridge: pixelfold found 1 of the reference's 2, all correct.
        let salt = &by_kind["salt-bridge"].prf1;
        assert!((salt.precision - 1.0).abs() < 1e-9);
        assert!((salt.recall - 0.5).abs() < 1e-9);

        // Pi-stacking is a reported kind that both found none of: a vacuous 1.0.
        assert_eq!(by_kind["pi-stacking"].prf1.f1, 1.0);
        assert_eq!(by_kind["pi-stacking"].predicted, 0);
    }

    /// A kind the reference does not report is not scored, so pixelfold finding
    /// some of it is not counted against precision.
    #[test]
    fn a_kind_the_reference_omits_is_not_scored() {
        let mut interactions = HashMap::new();
        interactions.insert(
            "disulfide".to_string(),
            [ordered_pair(key(1), key(2))].into_iter().collect(),
        );
        let prediction = Prediction {
            ss: vec![],
            hbonds: HashSet::new(),
            sasa: vec![],
            interactions,
            ligands: HashSet::new(),
        };

        // PLIP does not detect disulfides, so its golden does not list the kind.
        let golden = InteractionGolden {
            kinds: vec!["hydrogen-bond".to_string()],
            ligand_scoped: false,
            interactions: vec![],
        };

        let metrics = compare_interactions(&prediction, &golden);
        assert!(metrics.iter().all(|m| m.kind != "disulfide"));
        assert_eq!(metrics.len(), 1);
    }

    /// A ligand-scoped reference (PLIP) is compared only against pixelfold's
    /// ligand-touching edges, so the protein-protein interactions it never looked
    /// at are not charged against precision.
    #[test]
    fn a_ligand_scoped_reference_ignores_protein_protein_edges() {
        // Residue 9 is the ligand. Pixelfold found a ligand hydrogen bond (1-9)
        // and a protein-protein one (1-2); PLIP reports only the ligand's.
        let mut interactions = HashMap::new();
        interactions.insert(
            "hydrogen-bond".to_string(),
            [ordered_pair(key(1), key(9)), ordered_pair(key(1), key(2))]
                .into_iter()
                .collect(),
        );
        let prediction = Prediction {
            ss: vec![],
            hbonds: HashSet::new(),
            sasa: vec![],
            interactions,
            ligands: [key(9)].into_iter().collect(),
        };

        let golden = InteractionGolden {
            kinds: vec!["hydrogen-bond".to_string()],
            ligand_scoped: true,
            interactions: vec![edge("hydrogen-bond", 1, 9)],
        };

        let metric = &compare_interactions(&prediction, &golden)[0];
        // Only the ligand edge is scored: 1 predicted, matching 1 reference.
        assert_eq!(metric.predicted, 1);
        assert_eq!(
            metric.prf1.f1, 1.0,
            "the protein-protein edge is out of scope"
        );

        // Without the scope, the protein-protein edge would drop precision.
        let unscoped = InteractionGolden {
            ligand_scoped: false,
            ..golden
        };
        assert_eq!(compare_interactions(&prediction, &unscoped)[0].predicted, 2);
    }
}
