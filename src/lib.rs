use glam::Vec3;

pub mod parser;
pub mod renderer;

#[derive(Clone)]
pub struct Atom {
    pub serial: u32,
    pub name: String,
    pub residue_name: String,
    pub residue_seq: u32,
    pub chain_id: String,
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
            SecondaryStructure::Helix => (255, 100, 100),   // Red
            SecondaryStructure::Sheet => (100, 100, 255),   // Blue
            SecondaryStructure::Turn => (255, 255, 100),    // Yellow
            SecondaryStructure::Coil => (100, 255, 100),    // Green
        }
    }

    pub fn color_256(&self) -> u8 {
        match self {
            SecondaryStructure::Helix => 196,   // Red
            SecondaryStructure::Sheet => 33,    // Blue
            SecondaryStructure::Turn => 226,    // Yellow
            SecondaryStructure::Coil => 46,     // Green
        }
    }
}

pub struct Protein {
    pub atoms: Vec<Atom>,
    pub title: String,
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
                if i + offset + 4 < residues.len() {
                    if self.can_hbond(&residues[i + offset], &residues[i + offset + 4]) {
                        helix_count += 1;
                    }
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
                if i + j < residues.len() {
                    if self.can_hbond(&residues[i], &residues[i + j]) {
                        sheet_hbond_count += 1;
                    }
                }

                if i >= j {
                    if self.can_hbond(&residues[i - j], &residues[i]) {
                        sheet_hbond_count += 1;
                    }
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

            if curr == SecondaryStructure::Coil && (prev != SecondaryStructure::Coil || next != SecondaryStructure::Coil) {
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
        self.get_atom_by_name(residue, "CA").map(|atom| atom.position)
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