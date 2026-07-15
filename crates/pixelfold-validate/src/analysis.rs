//! Extract pixelfold's comparable output from a loaded structure.

use std::collections::{HashMap, HashSet};

use pixelfold_core::sasa::SurfaceCalculator;
use pixelfold_core::structure::Atom;
use pixelfold_core::{Protein, SecondaryStructure};

use crate::golden::ResidueKey;
use crate::metrics::Ss3;

/// Pixelfold's per-residue prediction for one structure.
pub struct Prediction {
    /// Secondary structure per residue, in chain then sequence order.
    pub ss: Vec<(ResidueKey, Ss3)>,
    /// Backbone hydrogen bonds as `(donor, acceptor)` residue pairs.
    pub hbonds: HashSet<(ResidueKey, ResidueKey)>,
    /// Solvent-accessible surface area per residue (A^2), residue order.
    pub sasa: Vec<(ResidueKey, f64)>,
}

/// Run the comparable analyses (DSSP is already assigned on load; SASA is
/// recomputed per atom) and collect them per residue.
pub fn predict(protein: &Protein) -> Prediction {
    let mut ss = Vec::new();
    let mut order: Vec<ResidueKey> = Vec::new();
    let mut seen: HashSet<ResidueKey> = HashSet::new();
    for atom in &protein.atoms {
        let key = residue_key(atom);
        if seen.insert(key.clone()) {
            order.push(key.clone());
            ss.push((key, reduce(atom.secondary_structure)));
        }
    }

    let mut hbonds = HashSet::new();
    for hb in &protein.hbonds {
        if let (Some(donor), Some(acceptor)) = (
            protein.atoms.get(hb.donor_atom_idx),
            protein.atoms.get(hb.acceptor_atom_idx),
        ) {
            hbonds.insert((residue_key(donor), residue_key(acceptor)));
        }
    }

    let atom_sasa = SurfaceCalculator::default().calculate_atom_sasa(&protein.atoms);
    let mut sasa_by_residue: HashMap<ResidueKey, f64> = HashMap::new();
    for (atom, &area) in protein.atoms.iter().zip(atom_sasa.iter()) {
        *sasa_by_residue.entry(residue_key(atom)).or_insert(0.0) += area as f64;
    }

    let sasa = order
        .into_iter()
        .map(|key| {
            let area = sasa_by_residue.get(&key).copied().unwrap_or(0.0);
            (key, area)
        })
        .collect();

    Prediction { ss, hbonds, sasa }
}

/// The DSSP single-letter code reduced to three states, for golden references.
pub fn reduce_dssp_letter(code: &str) -> Ss3 {
    match code.chars().next() {
        Some('H' | 'G' | 'I') => Ss3::Helix,
        Some('E' | 'B') => Ss3::Strand,
        _ => Ss3::Coil,
    }
}

fn reduce(ss: SecondaryStructure) -> Ss3 {
    match ss {
        SecondaryStructure::Helix => Ss3::Helix,
        SecondaryStructure::Sheet => Ss3::Strand,
        // Q3 has no turn state; DSSP T and S both fold into coil.
        SecondaryStructure::Turn | SecondaryStructure::Coil => Ss3::Coil,
    }
}

fn residue_key(atom: &Atom) -> ResidueKey {
    ResidueKey {
        chain: atom.chain_id.as_str().to_string(),
        seq: atom.residue_seq as i64,
        icode: atom.insertion_code.map(|c| c.to_string()),
    }
}
