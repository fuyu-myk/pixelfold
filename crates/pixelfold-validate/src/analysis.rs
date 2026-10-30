//! Extract pixelfold's comparable output from a loaded structure.

use std::collections::{HashMap, HashSet};

use pixelfold_core::interactions::detect;
use pixelfold_core::sasa::SurfaceCalculator;
use pixelfold_core::structure::Atom;
use pixelfold_core::{Protein, SecondaryStructure};

use crate::golden::ResidueKey;
use crate::metrics::Ss3;

/// A non-covalent interaction reduced to what is compared against a reference:
/// an unordered residue pair, keyed by interaction type.
pub type InteractionEdges = HashMap<String, HashSet<(ResidueKey, ResidueKey)>>;

/// Pixelfold's per-residue prediction for one structure.
pub struct Prediction {
    /// Secondary structure per residue, in chain then sequence order.
    pub ss: Vec<(ResidueKey, Ss3)>,
    /// Backbone hydrogen bonds as `(donor, acceptor)` residue pairs.
    pub hbonds: HashSet<(ResidueKey, ResidueKey)>,
    /// Solvent-accessible surface area per residue (A^2), residue order.
    pub sasa: Vec<(ResidueKey, f64)>,
    /// Non-covalent interactions as unordered residue pairs, per type.
    pub interactions: InteractionEdges,
    /// The water residues. A reference that leaves water out (RING) is compared
    /// only against edges touching none of these.
    pub waters: HashSet<ResidueKey>,
    /// Each standard amino acid's position in its chain, counting every residue
    /// of the chain in file order. A reference with a sequence gap (RING) pairs
    /// two of these only when their positions are at least the gap apart.
    pub chain_positions: HashMap<ResidueKey, usize>,
}

/// Order a residue pair canonically, so an edge and its reverse are one key.
pub fn ordered_pair(a: ResidueKey, b: ResidueKey) -> (ResidueKey, ResidueKey) {
    if a <= b { (a, b) } else { (b, a) }
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

    // SASA is compared against FreeSASA run on the polymer alone, so match its
    // atom set: drop waters, ligands, ions and hydrogens.
    let sasa_atoms: Vec<Atom> = protein
        .atoms
        .iter()
        .filter(|atom| atom.is_polymer() && !atom.element.as_str().eq_ignore_ascii_case("H"))
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

    let waters = protein
        .atoms
        .iter()
        .filter(|atom| is_water(atom.residue_name.as_str()))
        .map(residue_key)
        .collect();

    Prediction {
        ss,
        hbonds,
        sasa,
        interactions: interaction_edges(protein),
        waters,
        chain_positions: chain_positions(protein),
    }
}

/// True when a residue name denotes water.
fn is_water(residue_name: &str) -> bool {
    matches!(residue_name, "HOH" | "WAT" | "H2O" | "DOD" | "TIP")
}

/// True for the 20 standard amino acids, the only residues RING's sequence gap
/// applies to.
fn is_standard_amino_acid(residue_name: &str) -> bool {
    matches!(
        residue_name,
        "ALA"
            | "ARG"
            | "ASN"
            | "ASP"
            | "CYS"
            | "GLN"
            | "GLU"
            | "GLY"
            | "HIS"
            | "ILE"
            | "LEU"
            | "LYS"
            | "MET"
            | "PHE"
            | "PRO"
            | "SER"
            | "THR"
            | "TRP"
            | "TYR"
            | "VAL"
    )
}

/// Number every residue of each chain in file order, keeping the positions of
/// the standard amino acids.
fn chain_positions(protein: &Protein) -> HashMap<ResidueKey, usize> {
    let mut next_in_chain: HashMap<String, usize> = HashMap::new();
    let mut seen: HashSet<ResidueKey> = HashSet::new();
    let mut positions = HashMap::new();

    for atom in &protein.atoms {
        let key = residue_key(atom);
        if !seen.insert(key.clone()) {
            continue;
        }
        let next = next_in_chain.entry(key.chain.clone()).or_insert(0);
        if is_standard_amino_acid(atom.residue_name.as_str()) {
            positions.insert(key, *next);
        }
        *next += 1;
    }

    positions
}

/// Every detected interaction reduced to an unordered residue pair, grouped by
/// type. Interactions of one type between the same residue pair collapse, since
/// the reference tools report a residue pair once per type rather than once per
/// atom pair.
fn interaction_edges(protein: &Protein) -> InteractionEdges {
    let mut edges: InteractionEdges = HashMap::new();

    for interaction in detect(protein) {
        let (Some(&a), Some(&b)) = (interaction.atoms_a.first(), interaction.atoms_b.first())
        else {
            continue;
        };
        let (key_a, key_b) = (
            residue_key(&protein.atoms[a]),
            residue_key(&protein.atoms[b]),
        );
        if key_a == key_b {
            continue; // an interaction within one residue is not an edge
        }

        edges
            .entry(interaction.kind.label().to_string())
            .or_default()
            .insert(ordered_pair(key_a, key_b));
    }

    edges
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
        seq: i64::from(atom.residue_seq),
        icode: atom.insertion_code.map(|c| c.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::Vec3;
    use pixelfold_core::{FixedStr, Molecule, Protein, SecondaryStructure};

    fn atom(element: &str, is_hetatm: bool, residue_name: &str, seq: i32, pos: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(element),
            element: FixedStr::new(element),
            residue_name: FixedStr::new(residue_name),
            residue_seq: seq,
            chain_id: FixedStr::new("A"),
            is_hetatm,
            molecule: pixelfold_core::Molecule::Unknown,
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

    // A modified residue such as selenomethionine is a HETATM record of the
    // polymer, so FreeSASA is run with it; a ligand stays out.
    #[test]
    fn sasa_counts_modified_residues_but_not_ligands() {
        let molecule = |molecule, atom: Atom| Atom { molecule, ..atom };
        let carbon = molecule(Molecule::Protein, atom("C", false, "ALA", 1, Vec3::ZERO));
        let selenium = molecule(
            Molecule::Protein,
            atom("SE", true, "MSE", 2, Vec3::new(3.0, 0.0, 0.0)),
        );
        let glycerol = molecule(
            Molecule::NonPolymer,
            atom("C", true, "GOL", 3, Vec3::new(-1.0, 0.0, 0.0)),
        );

        let alone = predict(&protein(vec![carbon.clone()]));
        let with_mse = predict(&protein(vec![carbon.clone(), selenium.clone()]));
        let with_ligand = predict(&protein(vec![carbon, selenium, glycerol]));

        let area =
            |prediction: &Prediction, seq| residue_sasa(prediction, seq).expect("residue present");
        assert!(area(&with_mse, 2) > 0.0);
        assert!(area(&with_mse, 1) < area(&alone, 1));
        assert_eq!(area(&with_ligand, 1), area(&with_mse, 1));
        assert_eq!(area(&with_ligand, 3), 0.0);
    }
}
