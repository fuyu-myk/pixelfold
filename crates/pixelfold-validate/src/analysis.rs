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

    // SASA is compared against FreeSASA run on the protein alone, so match its
    // atom set: drop HETATM (waters, ligands, ions) and hydrogens.
    let sasa_atoms: Vec<Atom> = protein
        .atoms
        .iter()
        .filter(|atom| !atom.is_hetatm && !atom.element.as_str().eq_ignore_ascii_case("H"))
        .cloned()
        .collect();
    let atom_sasa = SurfaceCalculator::default().calculate_atom_sasa(&sasa_atoms);
    let mut sasa_by_residue: HashMap<ResidueKey, f64> = HashMap::new();
    for (atom, &area) in sasa_atoms.iter().zip(atom_sasa.iter()) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use glam::Vec3;
    use pixelfold_core::{FixedStr, Protein, SecondaryStructure};

    fn atom(element: &str, is_hetatm: bool, residue_name: &str, seq: u32, pos: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(element),
            element: FixedStr::new(element),
            residue_name: FixedStr::new(residue_name),
            residue_seq: seq,
            chain_id: FixedStr::new("A"),
            is_hetatm,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: pos,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    fn protein(atoms: Vec<Atom>) -> Protein {
        Protein {
            atoms,
            title: String::new(),
            surface_points: Vec::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        }
    }

    fn residue_sasa(prediction: &Prediction, seq: i64) -> Option<f64> {
        prediction
            .sasa
            .iter()
            .find(|(key, _)| key.seq == seq)
            .map(|(_, area)| *area)
    }

    // SASA must match FreeSASA's protein-only atom set: an overlapping water
    // (HETATM) and a hydrogen must neither shield nor contribute to a protein
    // residue's surface. Regression guard for the atom-set filter in `predict`.
    #[test]
    fn sasa_excludes_hetatm_and_hydrogen() {
        let carbon = atom("C", false, "ALA", 1, Vec3::ZERO);
        let baseline = predict(&protein(vec![carbon.clone()]));

        // Overlapping water + hydrogen would fully bury the carbon if counted.
        let water = atom("O", true, "HOH", 2, Vec3::ZERO);
        let hydrogen = atom("H", false, "ALA", 1, Vec3::ZERO);
        let with_solvent = predict(&protein(vec![carbon, water, hydrogen]));

        assert_eq!(residue_sasa(&baseline, 1), residue_sasa(&with_solvent, 1));
        assert!(residue_sasa(&baseline, 1).unwrap() > 0.0);
        // The water residue contributes no surface of its own.
        assert_eq!(residue_sasa(&with_solvent, 2), Some(0.0));
    }
}
