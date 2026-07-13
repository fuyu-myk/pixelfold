use glam::Vec3;

use crate::parser::HBond;
use crate::sasa::SurfacePoint;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DisplayMode {
    AllAtoms,
    Backbone,
}

#[derive(Clone)]
pub struct Atom {
    pub serial: u32,
    pub name: String,
    /// Element symbol (e.g. "C", "Ca", "Fe"), parsed from the structure's
    /// element column so a calcium ion is distinct from a C-alpha carbon.
    pub element: String,
    pub residue_name: String,
    pub residue_seq: u32,
    pub chain_id: String,
    pub is_hetatm: bool,
    pub altloc: Option<char>,
    pub occupancy: f32,
    pub insertion_code: Option<char>,
    pub model_number: usize,
    pub position: Vec3,
    pub b_factor: f32,
    pub secondary_structure: SecondaryStructure,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SecondaryStructure {
    Helix,
    Sheet,
    Turn,
    Coil,
}

impl SecondaryStructure {
    pub fn color_rgb(&self) -> (u8, u8, u8) {
        match self {
            SecondaryStructure::Helix => (255, 100, 100), // Red
            SecondaryStructure::Sheet => (100, 100, 255), // Blue
            SecondaryStructure::Turn => (255, 255, 100),  // Yellow
            SecondaryStructure::Coil => (100, 255, 100),  // Green
        }
    }

    pub fn color_256(&self) -> u8 {
        match self {
            SecondaryStructure::Helix => 196, // Red
            SecondaryStructure::Sheet => 33,  // Blue
            SecondaryStructure::Turn => 226,  // Yellow
            SecondaryStructure::Coil => 46,   // Green
        }
    }
}

pub struct Protein {
    pub atoms: Vec<Atom>,
    pub title: String,
    pub surface_points: Vec<SurfacePoint>,
    pub hbonds: Vec<HBond>,
    /// Set when the file declares a biological assembly that differs from the
    /// deposited asymmetric unit (the coordinates shown). `None` when the
    /// deposited coordinates already are the biological unit, or when no
    /// assembly is declared.
    pub assembly: Option<BiologicalAssembly>,
}

/// A biological assembly whose symmetry expansion differs from the deposited
/// asymmetric unit. Present only when the primary assembly applies more than the
/// identity operator, so the deposited coordinates are a partial view of the
/// functional unit. Assembly generation itself is not yet performed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BiologicalAssembly {
    /// Number of symmetry operators the primary assembly applies. Greater than
    /// one means the deposited coordinates are only part of the biological unit.
    pub operator_count: usize,
    /// Human-readable oligomeric state stated in the file, if any (mmCIF
    /// `oligomeric_details` or PDB `REMARK 350` biological-unit line, for
    /// example "tetrameric"). `None` when the file does not state one.
    pub oligomer: Option<String>,
}

pub struct SecondaryStructureAssignment {
    hbond_distance_threshold: f32, // in Å
    hbond_energy_threshold: f32,   // in kcal/mol
}

impl Default for SecondaryStructureAssignment {
    fn default() -> Self {
        Self {
            hbond_distance_threshold: 3.5,
            hbond_energy_threshold: -0.5,
        }
    }
}

impl SecondaryStructureAssignment {
    pub fn new(distance_threshold: f32, energy_threshold: f32) -> Self {
        Self {
            hbond_distance_threshold: distance_threshold,
            hbond_energy_threshold: energy_threshold,
        }
    }

    /// Assign secondary structure to residues based on hydrogen bonding patterns
    pub fn assign(&self, residues: &[Vec<Atom>]) -> Vec<SecondaryStructure> {
        let mut assignments = vec![SecondaryStructure::Coil; residues.len()];

        self.detect_helices(residues, &mut assignments);
        self.detect_sheets(residues, &mut assignments);
        self.detect_turns(residues, &mut assignments);

        assignments
    }

    /// Detect alpha-helices based on the heuristic of repeated i to i + 4 hydrogen bonds
    ///
    /// Alpha-helices have the following H-bond pattern:
    ///
    /// Residue i: N-H bonded to C=O of residue i + 4
    ///
    /// Typical i to i + 4 Cα distance: ~5.4 Å
    fn detect_helices(&self, residues: &[Vec<Atom>], assignments: &mut [SecondaryStructure]) {
        for i in 0..residues.len().saturating_sub(4) {
            let helix_window = 4;
            let mut helix_count = 0;

            for offset in 0..helix_window {
                if i + offset + 4 < residues.len()
                    && self.can_hbond(&residues[i + offset], &residues[i + offset + 4])
                {
                    helix_count += 1;
                }
            }

            if helix_count >= 3 {
                for offset in 0..=helix_window {
                    if i + offset < assignments.len() {
                        assignments[i + offset] = SecondaryStructure::Helix;
                    }
                }
            }
        }
    }

    /// Detect beta-sheets based on the heuristic of inter-strand (distant) hydrogen bonds
    ///
    /// Beta-sheets have residues H-bonded to non-adjacent residues.
    /// If residue i H-bonds with i ± 2 or i ± 3 (and not i + 4), it's likely a sheet.
    fn detect_sheets(&self, residues: &[Vec<Atom>], assignments: &mut [SecondaryStructure]) {
        for i in 0..residues.len() {
            if assignments[i] == SecondaryStructure::Helix {
                continue;
            }

            let mut sheet_hbond_count = 0;

            for j in 2..=3 {
                if i + j < residues.len() && self.can_hbond(&residues[i], &residues[i + j]) {
                    sheet_hbond_count += 1;
                }

                if i >= j && self.can_hbond(&residues[i - j], &residues[i]) {
                    sheet_hbond_count += 1;
                }
            }

            if sheet_hbond_count >= 2 {
                assignments[i] = SecondaryStructure::Sheet;
            }
        }
    }

    /// Detect turns based on regions of high curvature between secondary structures
    ///
    /// Following the simple heuristic of transitions between helices/sheets and coils
    fn detect_turns(&self, residues: &[Vec<Atom>], assignments: &mut [SecondaryStructure]) {
        for i in 1..assignments.len().saturating_sub(1) {
            let prev = assignments[i - 1];
            let curr = assignments[i];
            let next = assignments[i + 1];

            if curr == SecondaryStructure::Coil
                && (prev != SecondaryStructure::Coil || next != SecondaryStructure::Coil)
            {
                if let (Some(ca_prev), Some(ca_curr), Some(ca_next)) = (
                    self.get_ca(&residues[i - 1]),
                    self.get_ca(&residues[i]),
                    self.get_ca(&residues[i + 1]),
                ) {
                    let angle = self.compute_angle(ca_prev, ca_curr, ca_next);

                    if angle > 70.0 {
                        assignments[i] = SecondaryStructure::Turn;
                    }
                }
            }
        }
    }

    /// Simple heuristic to determine if two residues can form a hydrogen bond
    /// Given an O from the donor residue and an N-H from the acceptor residue,
    /// check if O...N distance is within the threshold
    fn can_hbond(&self, donor_residue: &[Atom], acceptor_residue: &[Atom]) -> bool {
        let o_donor = self.get_atom_by_name(donor_residue, "O");
        let n_acceptor = self.get_atom_by_name(acceptor_residue, "N");

        if let (Some(o), Some(n)) = (o_donor, n_acceptor) {
            let distance = (o.position - n.position).length();
            (2.4..=3.5).contains(&distance)
        } else {
            false
        }
    }

    /// Get coordinates of the Cα atom of a residue
    fn get_ca(&self, residue: &[Atom]) -> Option<Vec3> {
        self.get_atom_by_name(residue, "CA")
            .map(|atom| atom.position)
    }

    /// Get an atom by name from a residue
    fn get_atom_by_name<'a>(&self, residue: &'a [Atom], name: &str) -> Option<&'a Atom> {
        residue.iter().find(|atom| atom.name == name)
    }

    /// Compute the angle between three 3D points (a-b-c)
    fn compute_angle(&self, a: Vec3, b: Vec3, c: Vec3) -> f32 {
        let ab = (a - b).normalize();
        let cb = (c - b).normalize();
        let dot_product = ab.dot(cb).clamp(-1.0, 1.0);

        dot_product.acos().to_degrees()
    }
}

/// Calculate B-factor percentiles for normalization
pub fn calculate_bfactor_range(protein: &Protein) -> (f32, f32) {
    let mut b_factors: Vec<f32> = protein
        .atoms
        .iter()
        .map(|a| a.b_factor)
        .filter(|b| b.is_finite())
        .collect();
    if b_factors.is_empty() {
        return (0.0, 100.0);
    }
    b_factors.sort_by(f32::total_cmp);

    // Avoiding outliers
    let idx_min = (b_factors.len() as f32 * 0.02) as usize;
    let idx_max = (b_factors.len() as f32 * 0.98) as usize;

    let b_min = b_factors[idx_min.min(b_factors.len() - 1)];
    let b_max = b_factors[idx_max.min(b_factors.len() - 1)];

    (b_min, b_max)
}

/// Get indices of C-alpha atoms in the protein
pub fn get_calpha_indices(protein: &Protein) -> Vec<usize> {
    protein
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.name == "CA")
        .map(|(idx, _)| idx)
        .collect()
}

/// Get pairs of C-alpha indices that should be connected
/// Only connects sequential C-alphas in the same chain within distance threshold
pub fn get_calpha_connections(protein: &Protein, ca_indices: &[usize]) -> Vec<(usize, usize)> {
    let mut connections = Vec::new();
    let distance_threshold = 4.2; // Angstroms (typical CA-CA distance is ~3.8Å)

    // Sort C-alphas by chain first, then by residue sequence
    let mut sorted_cas: Vec<(String, u32, usize)> = ca_indices
        .iter()
        .map(|&idx| {
            let atom = &protein.atoms[idx];
            (atom.chain_id.clone(), atom.residue_seq, idx)
        })
        .collect();

    sorted_cas.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    for i in 0..sorted_cas.len().saturating_sub(1) {
        let (chain1, res1, idx1) = &sorted_cas[i];
        let (chain2, res2, idx2) = &sorted_cas[i + 1];

        // Only connect if same chain
        if chain1 == chain2 {
            // Exactly sequential (res2 = res1 + 1)
            if *res2 == *res1 + 1 {
                let atom1 = &protein.atoms[*idx1];
                let atom2 = &protein.atoms[*idx2];
                let distance = (atom2.position - atom1.position).length();

                if distance >= 2.5 && distance <= distance_threshold {
                    connections.push((*idx1, *idx2));
                }
            }
        }
    }

    connections
}

/// Policy for resolving atoms modelled in multiple alternate locations.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum AltlocPolicy {
    /// Keep the highest-occupancy conformer per atom; ties prefer altloc 'A'.
    #[default]
    Occupancy,
    /// Keep only altloc 'A' and atoms with no altloc.
    A,
    /// Keep only altloc 'B' and atoms with no altloc.
    B,
    /// Keep every conformer.
    All,
}

/// Apply an alternate-location policy. Atoms with no altloc are always kept.
pub fn filter_altlocs(atoms: Vec<Atom>, policy: AltlocPolicy) -> Vec<Atom> {
    match policy {
        AltlocPolicy::All => atoms,
        AltlocPolicy::A => atoms
            .into_iter()
            .filter(|a| a.altloc.is_none() || a.altloc == Some('A'))
            .collect(),
        AltlocPolicy::B => atoms
            .into_iter()
            .filter(|a| a.altloc.is_none() || a.altloc == Some('B'))
            .collect(),
        AltlocPolicy::Occupancy => {
            // Keep the best conformer per (chain, residue, insertion, atom name).
            let mut best: std::collections::HashMap<(String, u32, Option<char>, String), usize> =
                std::collections::HashMap::new();
            let mut kept: Vec<Atom> = Vec::new();

            for atom in atoms {
                if atom.altloc.is_none() {
                    kept.push(atom);
                    continue;
                }

                let key = (
                    atom.chain_id.clone(),
                    atom.residue_seq,
                    atom.insertion_code,
                    atom.name.clone(),
                );
                match best.get(&key).copied() {
                    None => {
                        best.insert(key, kept.len());
                        kept.push(atom);
                    }
                    Some(idx) => {
                        let existing = &kept[idx];
                        // Higher occupancy wins; on a tie the earlier altloc ('A' < 'B') wins.
                        let replace = atom.occupancy > existing.occupancy
                            || (atom.occupancy == existing.occupancy
                                && atom.altloc < existing.altloc);
                        if replace {
                            kept[idx] = atom;
                        }
                    }
                }
            }

            kept
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom_with_bfactor(b: f32) -> Atom {
        Atom {
            serial: 0,
            name: "CA".to_string(),
            element: "C".to_string(),
            residue_name: "ALA".to_string(),
            residue_seq: 1,
            chain_id: "A".to_string(),
            is_hetatm: false,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: Vec3::ZERO,
            b_factor: b,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    #[test]
    fn bfactor_range_does_not_panic_on_nan() {
        // Real PDB files contain NaN coordinates and B-factors; the depth/value
        // sorts must not panic on them.
        let protein = Protein {
            atoms: vec![
                atom_with_bfactor(10.0),
                atom_with_bfactor(f32::NAN),
                atom_with_bfactor(50.0),
            ],
            title: String::new(),
            surface_points: Vec::new(),
            hbonds: Vec::new(),
            assembly: None,
        };

        let (b_min, b_max) = calculate_bfactor_range(&protein);
        assert!(b_min.is_finite() && b_max.is_finite());
    }

    fn altloc_atom(name: &str, altloc: Option<char>, occupancy: f32) -> Atom {
        Atom {
            serial: 0,
            name: name.to_string(),
            element: "C".to_string(),
            residue_name: "SER".to_string(),
            residue_seq: 1,
            chain_id: "A".to_string(),
            is_hetatm: false,
            altloc,
            occupancy,
            insertion_code: None,
            model_number: 1,
            position: Vec3::ZERO,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    #[test]
    fn altloc_occupancy_keeps_highest_and_breaks_ties_with_a() {
        // Same atom, two conformers: B has higher occupancy so it wins.
        let kept = filter_altlocs(
            vec![
                altloc_atom("OG", Some('A'), 0.4),
                altloc_atom("OG", Some('B'), 0.6),
            ],
            AltlocPolicy::Occupancy,
        );
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].altloc, Some('B'));

        // On an occupancy tie, altloc 'A' wins even if listed second.
        let tie = filter_altlocs(
            vec![
                altloc_atom("OG", Some('B'), 0.5),
                altloc_atom("OG", Some('A'), 0.5),
            ],
            AltlocPolicy::Occupancy,
        );
        assert_eq!(tie.len(), 1);
        assert_eq!(tie[0].altloc, Some('A'));
    }

    #[test]
    fn altloc_a_keeps_a_and_none_drops_b() {
        let kept = filter_altlocs(
            vec![
                altloc_atom("OG", Some('A'), 0.5),
                altloc_atom("OG", Some('B'), 0.5),
                altloc_atom("N", None, 1.0),
            ],
            AltlocPolicy::A,
        );
        assert_eq!(kept.len(), 2);
        assert!(kept.iter().all(|a| a.altloc != Some('B')));
    }

    #[test]
    fn altloc_all_keeps_every_conformer() {
        let kept = filter_altlocs(
            vec![
                altloc_atom("OG", Some('A'), 0.5),
                altloc_atom("OG", Some('B'), 0.5),
            ],
            AltlocPolicy::All,
        );
        assert_eq!(kept.len(), 2);
    }
}
