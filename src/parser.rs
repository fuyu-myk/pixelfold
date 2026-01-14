use anyhow::Result;
use glam::Vec3;
use pdbtbx::*;
use std::{
    path::Path,
    collections::HashMap
};

use crate::{Atom, Protein, SecondaryStructure, visualization::surface::SurfaceCalculator};

/// Load a protein structure from a PDB or mmCIF file
pub fn load_protein<P: AsRef<Path>>(path: P) -> Result<Protein> {
    load_protein_with_options(path, false)
}

/// Load a protein structure with additional options
pub fn load_protein_with_options<P: AsRef<Path>>(path: P, skip_surface: bool) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path.to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;
    
    let (pdb, errors) = pdbtbx::open(path_str)
        .map_err(|e| anyhow::anyhow!("Failed to open file: {:?}", e))?;
    
    if !errors.is_empty() {
        eprintln!("Warning: {} errors while parsing file", errors.len());
    }

    let title = pdb.identifier.as_deref().unwrap_or("Unknown").to_string();

    // Use atoms_with_hierarchy to get residue information
    let mut atoms: Vec<Atom> = pdb
        .atoms_with_hierarchy()
        .map(|hierarchy| {
            let atom = hierarchy.atom();
            let residue = hierarchy.residue();
            let conformer = hierarchy.conformer();
            let chain = hierarchy.chain();
            
            Atom {
                serial: atom.serial_number() as u32,
                name: atom.name().to_string(),
                residue_name: conformer.name().to_string(),
                residue_seq: residue.serial_number() as u32,
                chain_id: chain.id().to_string(),
                position: Vec3::new(atom.x() as f32, atom.y() as f32, atom.z() as f32),
                b_factor: atom.b_factor() as f32,
                secondary_structure: SecondaryStructure::Coil, // Default to Coil
            }
        })
        .collect();

    let hbonds = assign_secondary_structures(&mut atoms);

    // Compute solvent-accessible surface (unless skipped for performance)
    let surface_points = if skip_surface {
        Vec::new()
    } else {
        let surface_calculator = SurfaceCalculator::default();
        surface_calculator.calculate_surface(&atoms)
    };

    Ok(Protein { atoms, title, surface_points, hbonds })
}

/// Load only backbone atoms (CA, C, N, O) for simplified rendering
pub fn load_protein_backbone<P: AsRef<Path>>(path: P) -> Result<Protein> {
    load_protein_backbone_with_options(path, false)
}

/// Load backbone atoms with additional options
pub fn load_protein_backbone_with_options<P: AsRef<Path>>(path: P, skip_surface: bool) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path.to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;
    
    let (pdb, errors) = pdbtbx::open(path_str)
        .map_err(|e| anyhow::anyhow!("Failed to open file: {:?}", e))?;
    
    if !errors.is_empty() {
        eprintln!("Warning: {} errors while parsing file", errors.len());
    }

    let title = pdb.identifier.as_deref().unwrap_or("Unknown").to_string();

    let mut atoms: Vec<Atom> = pdb
        .atoms_with_hierarchy()
        .filter(|hierarchy| {
            let name = hierarchy.atom().name();
            name == "CA" || name == "C" || name == "N" || name == "O"
        })
        .map(|hierarchy| {
            let atom = hierarchy.atom();
            let residue = hierarchy.residue();
            let conformer = hierarchy.conformer();
            let chain = hierarchy.chain();
            
            Atom {
                serial: atom.serial_number() as u32,
                name: atom.name().to_string(),
                residue_name: conformer.name().to_string(),
                residue_seq: residue.serial_number() as u32,
                chain_id: chain.id().to_string(),
                position: Vec3::new(atom.x() as f32, atom.y() as f32, atom.z() as f32),
                b_factor: atom.b_factor() as f32,
                secondary_structure: SecondaryStructure::Coil, // Default to Coil
            }
        })
        .collect();

    let hbonds = assign_secondary_structures(&mut atoms);

    // Compute solvent-accessible surface (unless skipped for performance)
    let surface_points = if skip_surface {
        Vec::new()
    } else {
        let surface_calculator = SurfaceCalculator::default();
        surface_calculator.calculate_surface(&atoms)
    };

    Ok(Protein { atoms, title, surface_points, hbonds })
}

/// Load only CA (alpha carbon) atoms for minimal rendering of large proteins
pub fn load_protein_ca_only<P: AsRef<Path>>(path: P) -> Result<Protein> {
    load_protein_ca_only_with_options(path, false)
}

/// Load CA atoms with additional options
pub fn load_protein_ca_only_with_options<P: AsRef<Path>>(path: P, skip_surface: bool) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path.to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;
    
    let (pdb, errors) = pdbtbx::open(path_str)
        .map_err(|e| anyhow::anyhow!("Failed to open file: {:?}", e))?;
    
    if !errors.is_empty() {
        eprintln!("Warning: {} errors while parsing file", errors.len());
    }

    let title = pdb.identifier.as_deref().unwrap_or("Unknown").to_string();

    let mut atoms: Vec<Atom> = pdb
        .atoms_with_hierarchy()
        .filter(|hierarchy| hierarchy.atom().name() == "CA")
        .map(|hierarchy| {
            let atom = hierarchy.atom();
            let residue = hierarchy.residue();
            let conformer = hierarchy.conformer();
            let chain = hierarchy.chain();
            
            Atom {
                serial: atom.serial_number() as u32,
                name: atom.name().to_string(),
                residue_name: conformer.name().to_string(),
                residue_seq: residue.serial_number() as u32,
                chain_id: chain.id().to_string(),
                position: Vec3::new(atom.x() as f32, atom.y() as f32, atom.z() as f32),
                b_factor: atom.b_factor() as f32,
                secondary_structure: SecondaryStructure::Coil, // Default to Coil
            }
        })
        .collect();

    let hbonds = assign_secondary_structures(&mut atoms);

    // Compute solvent-accessible surface (unless skipped for performance)
    let surface_points = if skip_surface {
        Vec::new()
    } else {
        let surface_calculator = SurfaceCalculator::default();
        surface_calculator.calculate_surface(&atoms)
    };

    Ok(Protein { atoms, title, surface_points, hbonds })
}

/// DSSP hydrogen bond for secondary structure assignment and visualization
#[derive(Clone, Copy, Debug)]
pub struct HBond {
    pub donor_residue: usize,        // Residue index (NH donor)
    pub acceptor_residue: usize,     // Residue index (CO acceptor)
    pub donor_atom_idx: usize,       // Atom index of N in donor residue
    pub acceptor_atom_idx: usize,    // Atom index of O in acceptor residue
    pub energy: f32,                 // kcal/mol
}

/// DSSP-based secondary structure assigner
/// 
/// Reference: Kabsch, W., & Sander, C. (1983). Dictionary of protein secondary structure: 
/// Pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22(12), 2577–2637.
struct DSSPAssigner {
    energy_threshold: f32,  // kcal/mol, typically -0.5
}

impl Default for DSSPAssigner {
    fn default() -> Self {
        Self {
            energy_threshold: -0.5,
        }
    }
}

impl DSSPAssigner {
    /// Assign secondary structure using DSSP algorithm and return H-bonds
    fn assign(&self, residues: &[Vec<Atom>]) -> (Vec<SecondaryStructure>, Vec<HBond>) {
        let n = residues.len();
        let mut ss = vec![SecondaryStructure::Coil; n];

        let hbonds = self.find_hbonds(residues);
        self.assign_from_hbonds(&hbonds, &mut ss, n);

        (ss, hbonds)
    }

    /// Find all backbone H-bonds in the protein
    /// 
    /// DSSP definition:
    /// E = 0.084 * { 1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN } * 332 kcal/mol
    /// 
    /// where:
    /// - r_ON = distance from carbonyl O to backbone N
    /// - r_CH = distance from carbonyl C to backbone H
    /// - r_OH = distance from carbonyl O to backbone H
    /// - r_CN = distance from carbonyl C to backbone N
    fn find_hbonds(&self, residues: &[Vec<Atom>]) -> Vec<HBond> {
        let mut hbonds = Vec::new();

        for (i, res_i) in residues.iter().enumerate() {
            for (j, res_j) in residues.iter().enumerate() {
                if i == j || (i as i32 - j as i32).abs() == 1 {
                    continue; // Skip same and adjacent residues
                }

                // Get backbone atoms
                let Some((c_i, o_i)) = self.get_co_atoms(res_i) else {
                    continue;
                };
                let Some((n_j, h_j)) = self.get_nh_atoms(res_j) else {
                    continue;
                };

                let energy = self.calculate_hbond_energy(c_i, o_i, n_j, h_j);

                if energy < self.energy_threshold {
                    hbonds.push(HBond {
                        donor_residue: j,
                        acceptor_residue: i,
                        donor_atom_idx: 0,  // Mapped later
                        acceptor_atom_idx: 0,  // Mapped later
                        energy,
                    });
                }
            }
        }

        hbonds
    }

    /// Calculate DSSP hydrogen bond energy
    /// 
    /// E = 0.084 * { 1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN } * 332 kcal/mol
    fn calculate_hbond_energy(
        &self,
        c: Vec3,
        o: Vec3,
        n: Vec3,
        h: Vec3,
    ) -> f32 {
        const K1: f32 = 0.084;
        const K2: f32 = 332.0;

        let r_on = (o - n).length();
        let r_ch = (c - h).length();
        let r_oh = (o - h).length();
        let r_cn = (c - n).length();

        // Avoid division by zero
        if r_on < 0.1 || r_ch < 0.1 || r_oh < 0.1 || r_cn < 0.1 {
            return 0.0;
        }

        K1 * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn) * K2
    }

    /// Get C and O atoms from a residue's carbonyl group
    fn get_co_atoms(&self, residue: &[Atom]) -> Option<(Vec3, Vec3)> {
        let c = residue.iter().find(|a| a.name == "C")?.position;
        let o = residue.iter().find(|a| a.name == "O")?.position;

        Some((c, o))
    }

    /// Get N and H atoms from a residue's amide group
    /// 
    /// Note: PDB files typically don't include H atoms unless explicitly modeled.
    /// H position is computed 1.0 Å from N, opposite the direction of C-N bond.
    fn get_nh_atoms(&self, residue: &[Atom]) -> Option<(Vec3, Vec3)> {
        let n = residue.iter().find(|a| a.name == "N")?.position;
        let ca = residue.iter().find(|a| a.name == "CA")?.position;

        // H position: N + normalized(N - CA) * 1.0
        let h_direction = (n - ca).normalize();
        let h = n + h_direction * 1.0;

        Some((n, h))
    }

    /// Assign secondary structure based on H-bond patterns
    fn assign_from_hbonds(
        &self,
        hbonds: &[HBond],
        ss: &mut [SecondaryStructure],
        n: usize,
    ) {
        // For each residue, track which residues donate H-bonds to it
        let mut acceptor_map: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut donor_map: Vec<Vec<usize>> = vec![Vec::new(); n];

        for hbond in hbonds {
            acceptor_map[hbond.acceptor_residue].push(hbond.donor_residue);
            donor_map[hbond.donor_residue].push(hbond.acceptor_residue);
        }

        // Detect helices (i -> i+3, i+4, i+5 patterns)
        for i in 0..n {
            for &acceptor in &donor_map[i] {
                let offset = (acceptor as i32 - i as i32).abs();

                // α-helix: i to i+4 H-bonds (most common)
                if offset == 4 {
                    ss[i] = SecondaryStructure::Helix;
                    ss[acceptor] = SecondaryStructure::Helix;
                }

                // Mark residues in between for continuous helices
                if offset == 4 && i + 4 < n {
                    for ss in ss.iter_mut().take(i + 4).skip(i + 1) {
                        if *ss == SecondaryStructure::Coil {
                            *ss = SecondaryStructure::Helix;
                        }
                    }
                }

                // 3₁₀ helix: i to i+3 H-bonds (mapped to Helix for visualization)
                if offset == 3 && ss[i] == SecondaryStructure::Coil {
                    ss[i] = SecondaryStructure::Helix;
                    ss[acceptor] = SecondaryStructure::Helix;
                }

                // π-helix: i to i+5 H-bonds (mapped to Helix for visualization)
                if offset == 5 && ss[i] == SecondaryStructure::Coil {
                    ss[i] = SecondaryStructure::Helix;
                    ss[acceptor] = SecondaryStructure::Helix;
                }
            }
        }

        // Detect β-sheets (parallel and antiparallel)
        for i in 0..n {
            if ss[i] == SecondaryStructure::Helix {
                continue; // Helices take precedence
            }

            // Look for non-local H-bonds (sheet pattern)
            for &donor in &acceptor_map[i] {
                let offset = (donor as i32 - i as i32).abs();

                // Sheet H-bonds are typically non-local (offset > 4)
                if offset > 4 {
                    ss[i] = SecondaryStructure::Sheet;
                    ss[donor] = SecondaryStructure::Sheet;
                }
            }

            for &acceptor in &donor_map[i] {
                let offset = (acceptor as i32 - i as i32).abs();

                if offset > 4 {
                    ss[i] = SecondaryStructure::Sheet;
                    ss[acceptor] = SecondaryStructure::Sheet;
                }
            }
        }

        // Detect turns in non-helix, non-sheet regions
        for i in 1..n.saturating_sub(1) {
            if ss[i] == SecondaryStructure::Coil {
                let prev = ss.get(i.saturating_sub(1)).copied().unwrap_or(SecondaryStructure::Coil);
                let next = ss.get(i + 1).copied().unwrap_or(SecondaryStructure::Coil);

                // Turn = transition between secondary structures
                if prev != SecondaryStructure::Coil || next != SecondaryStructure::Coil {
                    ss[i] = SecondaryStructure::Turn;
                }
            }
        }
    }
}

/// Group atoms by residue and assign secondary structures using DSSP
/// 
/// Returns H-bonds with atom indices mapped from residue indices
fn assign_secondary_structures(atoms: &mut [Atom]) -> Vec<HBond> {
    let mut residues: HashMap<u32, Vec<Atom>> = HashMap::new();
    for atom in atoms.iter() {
        residues.entry(atom.residue_seq)
            .or_default()
            .push((*atom).clone());
    }

    let mut residue_numbers: Vec<u32> = residues.keys().copied().collect();
    residue_numbers.sort_unstable();

    let residue_vec: Vec<Vec<Atom>> = residue_numbers
        .iter()
        .map(|&num| residues[&num].clone())
        .collect();

    let assigner = DSSPAssigner::default();
    let (assignments, mut hbonds) = assigner.assign(&residue_vec);

    let mut ss_map: HashMap<u32, SecondaryStructure> = HashMap::new();
    for (i, &residue_num) in residue_numbers.iter().enumerate() {
        ss_map.insert(residue_num, assignments[i]);
    }
    
    for atom in atoms.iter_mut() {
        if let Some(&ss) = ss_map.get(&atom.residue_seq) {
            atom.secondary_structure = ss;
        }
    }

    // Map residue indices to atom indices for H-bonds
    for hbond in hbonds.iter_mut() {
        let donor_res_num = residue_numbers[hbond.donor_residue];
        let acceptor_res_num = residue_numbers[hbond.acceptor_residue];

        // N atom in donor residue
        if let Some(donor_n_idx) = atoms.iter().position(|a| {
            a.residue_seq == donor_res_num && a.name == "N"
        }) {
            hbond.donor_atom_idx = donor_n_idx;
        }

        // O atom in acceptor residue
        if let Some(acceptor_o_idx) = atoms.iter().position(|a| {
            a.residue_seq == acceptor_res_num && a.name == "O"
        }) {
            hbond.acceptor_atom_idx = acceptor_o_idx;
        }
    }

    hbonds
}
