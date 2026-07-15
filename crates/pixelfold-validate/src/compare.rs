//! Compare a prediction against golden references, aligning by residue.

use std::collections::{HashMap, HashSet};

use crate::analysis::{Prediction, reduce_dssp_letter};
use crate::golden::{DsspGolden, ResidueKey, SasaGolden};
use crate::metrics::{ErrorStats, PrF1, Ss3, error_stats, q3_agreement, set_prf1};

/// Metrics for one benchmark entry; each field is `None` when its golden is
/// absent (or, for H-bonds, when the reference does not report them).
pub struct EntryMetrics {
    pub id: String,
    pub q3: Option<f64>,
    pub hbond_f1: Option<PrF1>,
    pub sasa: Option<ErrorStats>,
    pub aligned_residues: usize,
}

/// Compare a prediction against whichever golden references are available.
pub fn compare(
    id: &str,
    prediction: &Prediction,
    dssp: Option<&DsspGolden>,
    sasa: Option<&SasaGolden>,
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
    }
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
    use crate::golden::{DsspResidue, SasaResidue};

    fn key(seq: i64) -> ResidueKey {
        ResidueKey {
            chain: "A".to_string(),
            seq,
            icode: None,
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

        let m = compare("self", &prediction, Some(&dssp), Some(&sasa));
        assert_eq!(m.q3, Some(1.0));
        assert_eq!(m.aligned_residues, 3);
        assert_eq!(m.hbond_f1.unwrap().f1, 1.0);
        let sasa = m.sasa.unwrap();
        assert_eq!(sasa.mae, 0.0);
        assert!((sasa.pearson - 1.0).abs() < 1e-9);
    }

    #[test]
    fn residues_missing_from_the_prediction_are_skipped() {
        let prediction = Prediction {
            ss: vec![(key(1), Ss3::Helix)],
            hbonds: HashSet::new(),
            sasa: vec![],
        };
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
        let m = compare("x", &prediction, Some(&dssp), None);
        assert_eq!(m.aligned_residues, 1); // only residue 1 aligned
        assert_eq!(m.q3, Some(1.0));
        assert!(m.hbond_f1.is_none()); // reference omitted H-bonds
    }
}
