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
    /// Every interaction reference, in the order reports list them.
    pub const ALL: [Tool; 2] = [Tool::Plip, Tool::Ring];

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
/// Pixelfold's predicted set is first restricted to the part of the structure
/// the reference looked at, or every interaction it never examined would count
/// against precision. A ligand-scoped reference (PLIP) covers only edges that
/// touch a ligand it profiled, so a HETATM group it chose not to profile is out
/// too; a reference without water (RING) covers only edges that touch none; a
/// reference with a sequence gap (RING) never pairs amino acids of one chain
/// closer than that gap.
fn compare_interactions(prediction: &Prediction, golden: &InteractionGolden) -> Vec<KindMetric> {
    let mut reference: HashMap<&str, HashSet<(ResidueKey, ResidueKey)>> = HashMap::new();
    for edge in &golden.interactions {
        reference
            .entry(edge.kind.as_str())
            .or_default()
            .insert(ordered_pair(edge.a.clone(), edge.b.clone()));
    }

    let ligands: HashSet<&ResidueKey> = golden.ligands.iter().collect();
    let touches_ligand =
        |pair: &(ResidueKey, ResidueKey)| ligands.contains(&pair.0) || ligands.contains(&pair.1);
    let in_scope = |pair: &(ResidueKey, ResidueKey)| {
        (!golden.ligand_scoped || touches_ligand(pair))
            && (!golden.water_excluded || !touches(&prediction.waters, pair))
            && apart_in_sequence(
                &prediction.chain_positions,
                pair,
                golden.sequence_separation,
            )
    };

    let empty = HashSet::new();
    golden
        .kinds
        .iter()
        .map(|kind| {
            let all = prediction.interactions.get(kind).unwrap_or(&empty);
            let predicted: HashSet<(ResidueKey, ResidueKey)> =
                all.iter().filter(|p| in_scope(p)).cloned().collect();
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

fn touches(residues: &HashSet<ResidueKey>, pair: &(ResidueKey, ResidueKey)) -> bool {
    residues.contains(&pair.0) || residues.contains(&pair.1)
}

/// False for two standard amino acids of one chain fewer than `separation`
/// positions apart. Positions count residues, not author numbers, so insertion
/// codes widen the distance; ligands, modified residues and nucleotides have no
/// position and pair at any distance, as do residues of different chains.
fn apart_in_sequence(
    positions: &HashMap<ResidueKey, usize>,
    pair: &(ResidueKey, ResidueKey),
    separation: u32,
) -> bool {
    match (positions.get(&pair.0), positions.get(&pair.1)) {
        (Some(a), Some(b)) if pair.0.chain == pair.1.chain => a.abs_diff(*b) >= separation as usize,
        _ => true,
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
            waters: HashSet::new(),
            chain_positions: HashMap::new(),
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
            waters: HashSet::new(),
            chain_positions: HashMap::new(),
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
            waters: HashSet::new(),
            chain_positions: HashMap::new(),
        };

        let golden = InteractionGolden {
            kinds: vec![
                "hydrogen-bond".to_string(),
                "salt-bridge".to_string(),
                "pi-stacking".to_string(),
            ],
            ligand_scoped: false,
            ligands: vec![],
            water_excluded: false,
            sequence_separation: 0,
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
            waters: HashSet::new(),
            chain_positions: HashMap::new(),
        };

        // PLIP does not detect disulfides, so its golden does not list the kind.
        let golden = InteractionGolden {
            kinds: vec!["hydrogen-bond".to_string()],
            ligand_scoped: false,
            ligands: vec![],
            water_excluded: false,
            sequence_separation: 0,
            interactions: vec![],
        };

        let metrics = compare_interactions(&prediction, &golden);
        assert!(metrics.iter().all(|m| m.kind != "disulfide"));
        assert_eq!(metrics.len(), 1);
    }

    /// A ligand-scoped reference (PLIP) is compared only against pixelfold's
    /// edges to a ligand it profiled, so neither the protein-protein interactions
    /// it never looked at nor a HETATM group it chose to skip (a crystallisation
    /// additive, a modified residue) is charged against precision.
    #[test]
    fn a_ligand_scoped_reference_scores_only_the_ligands_it_profiled() {
        // Residue 9 is the ligand PLIP profiled and 7 a glycerol it skipped.
        // Pixelfold found a hydrogen bond to each, and a protein-protein one.
        let mut interactions = HashMap::new();
        interactions.insert(
            "hydrogen-bond".to_string(),
            [
                ordered_pair(key(1), key(9)),
                ordered_pair(key(1), key(7)),
                ordered_pair(key(1), key(2)),
            ]
            .into_iter()
            .collect(),
        );
        let prediction = Prediction {
            ss: vec![],
            hbonds: HashSet::new(),
            sasa: vec![],
            interactions,
            waters: HashSet::new(),
            chain_positions: HashMap::new(),
        };

        let golden = InteractionGolden {
            kinds: vec!["hydrogen-bond".to_string()],
            ligand_scoped: true,
            ligands: vec![key(9)],
            water_excluded: false,
            sequence_separation: 0,
            interactions: vec![edge("hydrogen-bond", 1, 9)],
        };

        let metric = &compare_interactions(&prediction, &golden)[0];
        assert_eq!(metric.predicted, 1);
        assert_eq!(
            metric.prf1.f1, 1.0,
            "only the profiled ligand's edge is in scope"
        );

        let unscoped = InteractionGolden {
            ligand_scoped: false,
            ..golden
        };
        assert_eq!(compare_interactions(&prediction, &unscoped)[0].predicted, 3);
    }

    fn chain_key(chain: &str, seq: i64) -> ResidueKey {
        ResidueKey {
            chain: chain.to_string(),
            seq,
            icode: None,
        }
    }

    /// A reference that leaves water out (RING by default) is not charged with
    /// pixelfold's edges to a water, which it never saw.
    #[test]
    fn a_reference_without_water_ignores_water_edges() {
        // Residue 9 is a water coordinating the metal at residue 1.
        let mut interactions = HashMap::new();
        interactions.insert(
            "metal-coordination".to_string(),
            [ordered_pair(key(1), key(9)), ordered_pair(key(1), key(2))]
                .into_iter()
                .collect(),
        );
        let prediction = Prediction {
            waters: [key(9)].into_iter().collect(),
            interactions,
            ..prediction()
        };

        let golden = InteractionGolden {
            kinds: vec!["metal-coordination".to_string()],
            ligand_scoped: false,
            ligands: vec![],
            water_excluded: true,
            sequence_separation: 0,
            interactions: vec![edge("metal-coordination", 1, 2)],
        };

        let metric = &compare_interactions(&prediction, &golden)[0];
        assert_eq!(metric.predicted, 1);
        assert_eq!(metric.prf1.f1, 1.0, "the water edge is out of scope");

        let with_water = InteractionGolden {
            water_excluded: false,
            ..golden
        };
        assert_eq!(
            compare_interactions(&prediction, &with_water)[0].predicted,
            2
        );
    }

    /// A reference that never pairs residues of one chain closer than its
    /// sequence gap (RING's 3) is compared only against pixelfold's pairs at
    /// least that far apart; pairs across chains are always in scope.
    #[test]
    fn a_sequence_gap_drops_close_pairs_within_a_chain() {
        let mut interactions = HashMap::new();
        interactions.insert(
            "hydrogen-bond".to_string(),
            [
                ordered_pair(chain_key("A", 1), chain_key("A", 3)), // 2 apart: out
                ordered_pair(chain_key("A", 1), chain_key("A", 4)), // 3 apart: in
                ordered_pair(chain_key("A", 5), chain_key("A", 3)), // reversed, 2 apart: out
                ordered_pair(chain_key("A", 1), chain_key("B", 2)), // across chains: in
                ordered_pair(chain_key("A", 20), chain_key("A", 21)), // 3 positions apart: in
                ordered_pair(chain_key("A", 99), chain_key("A", 100)), // a ligand: in
            ]
            .into_iter()
            .collect(),
        );
        // Positions count residues, so two inserted residues (20A, 20B) put 20
        // and 21 three apart; residue 100 is a ligand and has no position.
        let chain_positions = [
            (chain_key("A", 1), 0),
            (chain_key("A", 3), 2),
            (chain_key("A", 4), 3),
            (chain_key("A", 5), 4),
            (chain_key("A", 20), 19),
            (chain_key("A", 21), 22),
            (chain_key("A", 99), 98),
            (chain_key("B", 2), 1),
        ]
        .into_iter()
        .collect();
        let prediction = Prediction {
            interactions,
            chain_positions,
            ..prediction()
        };

        let golden = InteractionGolden {
            kinds: vec!["hydrogen-bond".to_string()],
            ligand_scoped: false,
            ligands: vec![],
            water_excluded: false,
            sequence_separation: 3,
            interactions: vec![],
        };
        assert_eq!(compare_interactions(&prediction, &golden)[0].predicted, 4);

        let no_gap = InteractionGolden {
            sequence_separation: 0,
            ..golden
        };
        assert_eq!(compare_interactions(&prediction, &no_gap)[0].predicted, 6);
    }
}
