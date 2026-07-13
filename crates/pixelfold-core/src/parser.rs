use anyhow::Result;
use glam::Vec3;
use pdbtbx::*;
use std::{
    collections::{BTreeSet, HashMap},
    path::Path,
};

use crate::{
    fixed_str::FixedStr,
    sasa::SurfaceCalculator,
    structure::{AltlocPolicy, Atom, Protein, SecondaryStructure, filter_altlocs},
};

/// Build `Atom` from a pdbtbx hierarchy item.
fn build_atom(hierarchy: &impl ContainsAtomConformerResidueChainModel) -> Atom {
    let atom = hierarchy.atom();
    let residue = hierarchy.residue();
    let conformer = hierarchy.conformer();
    let chain = hierarchy.chain();
    let element = atom
        .element()
        .map(|e| e.symbol().to_string())
        .unwrap_or_else(|| infer_element(atom.name()));
    let altloc = conformer
        .alternative_location()
        .and_then(|s| s.chars().next());
    let insertion_code = residue.insertion_code().and_then(|s| s.chars().next());
    let position = Vec3::new(atom.x() as f32, atom.y() as f32, atom.z() as f32);

    Atom {
        serial: atom.serial_number() as u32,
        name: FixedStr::new(atom.name()),
        element: FixedStr::new(&element),
        residue_name: FixedStr::new(conformer.name()),
        residue_seq: residue.serial_number() as u32,
        chain_id: FixedStr::new(chain.id()),
        is_hetatm: atom.hetero(),
        altloc,
        occupancy: atom.occupancy() as f32,
        insertion_code,
        model_number: hierarchy.model().serial_number(),
        position,
        b_factor: atom.b_factor() as f32,
        secondary_structure: SecondaryStructure::Coil,
    }
}

/// Best-effort element inference from an atom name when the element column is
/// absent. Uses the leading letter, which is correct for standard protein
/// atoms but cannot disambiguate metals (for example a calcium named "CA")
/// without the column; those are expected to carry an explicit element symbol.
fn infer_element(name: &str) -> String {
    name.chars()
        .find(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase().to_string())
        .unwrap_or_default()
}

/// Format a one-line warning for identifiers that overflowed their fixed field
/// capacity, or `None` when nothing was truncated.
fn truncation_warning(truncated: &BTreeSet<(&'static str, String)>) -> Option<String> {
    if truncated.is_empty() {
        return None;
    }
    let list: Vec<String> = truncated
        .iter()
        .map(|(field, value)| format!("{field} '{value}'"))
        .collect();
    Some(format!(
        "{} identifier(s) exceeded the fixed field capacity and were truncated: {}",
        truncated.len(),
        list.join(", ")
    ))
}

/// Load a protein structure from a PDB or mmCIF file
pub fn load_protein<P: AsRef<Path>>(path: P) -> Result<Protein> {
    load_protein_with_options(path, false, AltlocPolicy::default())
}

/// Load a protein structure with additional options
pub fn load_protein_with_options<P: AsRef<Path>>(
    path: P,
    skip_surface: bool,
    altloc: AltlocPolicy,
) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;

    let (pdb, errors) = ReadOptions::default()
        .set_level(StrictnessLevel::Loose)
        .set_only_first_model(true)
        .read(path_str)
        .map_err(|e| anyhow::anyhow!("Failed to open file: {:?}", e))?;

    if !errors.is_empty() {
        eprintln!("Warning: {} errors while parsing file", errors.len());
    }

    let title = pdb.identifier.as_deref().unwrap_or("Unknown").to_string();

    // Build atoms, noting any identifier that overflowed its fixed-capacity
    // field.
    let mut truncated: BTreeSet<(&'static str, String)> = BTreeSet::new();
    let atoms: Vec<Atom> = pdb
        .atoms_with_hierarchy()
        .map(|hierarchy| {
            let atom = build_atom(&hierarchy);
            if atom.name.as_str() != hierarchy.atom().name() {
                truncated.insert(("atom name", hierarchy.atom().name().to_string()));
            }
            if atom.residue_name.as_str() != hierarchy.conformer().name() {
                truncated.insert(("residue name", hierarchy.conformer().name().to_string()));
            }
            if atom.chain_id.as_str() != hierarchy.chain().id() {
                truncated.insert(("chain id", hierarchy.chain().id().to_string()));
            }

            atom
        })
        .collect();

    if let Some(message) = truncation_warning(&truncated) {
        eprintln!("Warning: {message}");
    }

    let mut atoms = filter_altlocs(atoms, altloc);

    let hbonds = assign_secondary_structures(&mut atoms);

    // Compute solvent-accessible surface (unless skipped for performance)
    let surface_points = if skip_surface {
        Vec::new()
    } else {
        let surface_calculator = SurfaceCalculator::default();
        surface_calculator.calculate_surface(&atoms)
    };

    // A structure whose biological unit differs from the deposited asymmetric
    // unit is only detectable from the raw file; a gzip or read error just
    // means no warning.
    let assembly = std::fs::read_to_string(path)
        .ok()
        .and_then(|text| crate::assembly::detect_partial_assembly(&text));

    Ok(Protein {
        atoms,
        title,
        surface_points,
        hbonds,
        assembly,
    })
}

/// Load only backbone atoms (CA, C, N, O) for simplified rendering
pub fn load_protein_backbone<P: AsRef<Path>>(path: P) -> Result<Protein> {
    load_protein_backbone_with_options(path, false)
}

/// Load backbone atoms with additional options
pub fn load_protein_backbone_with_options<P: AsRef<Path>>(
    path: P,
    skip_surface: bool,
) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;

    let (pdb, errors) = ReadOptions::default()
        .set_level(StrictnessLevel::Loose)
        .set_only_first_model(true)
        .read(path_str)
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
        .map(|hierarchy| build_atom(&hierarchy))
        .collect();

    let hbonds = assign_secondary_structures(&mut atoms);

    // Compute solvent-accessible surface (unless skipped for performance)
    let surface_points = if skip_surface {
        Vec::new()
    } else {
        let surface_calculator = SurfaceCalculator::default();
        surface_calculator.calculate_surface(&atoms)
    };

    Ok(Protein {
        atoms,
        title,
        surface_points,
        hbonds,
        assembly: None,
    })
}

/// Load only CA (alpha carbon) atoms for minimal rendering of large proteins
pub fn load_protein_ca_only<P: AsRef<Path>>(path: P) -> Result<Protein> {
    load_protein_ca_only_with_options(path, false)
}

/// Load CA atoms with additional options
pub fn load_protein_ca_only_with_options<P: AsRef<Path>>(
    path: P,
    skip_surface: bool,
) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;

    let (pdb, errors) = ReadOptions::default()
        .set_level(StrictnessLevel::Loose)
        .set_only_first_model(true)
        .read(path_str)
        .map_err(|e| anyhow::anyhow!("Failed to open file: {:?}", e))?;

    if !errors.is_empty() {
        eprintln!("Warning: {} errors while parsing file", errors.len());
    }

    let title = pdb.identifier.as_deref().unwrap_or("Unknown").to_string();

    let mut atoms: Vec<Atom> = pdb
        .atoms_with_hierarchy()
        .filter(|hierarchy| hierarchy.atom().name() == "CA")
        .map(|hierarchy| build_atom(&hierarchy))
        .collect();

    let hbonds = assign_secondary_structures(&mut atoms);

    // Compute solvent-accessible surface (unless skipped for performance)
    let surface_points = if skip_surface {
        Vec::new()
    } else {
        let surface_calculator = SurfaceCalculator::default();
        surface_calculator.calculate_surface(&atoms)
    };

    Ok(Protein {
        atoms,
        title,
        surface_points,
        hbonds,
        assembly: None,
    })
}

/// DSSP hydrogen bond for secondary structure assignment and visualization
#[derive(Clone, Copy, Debug)]
pub struct HBond {
    pub donor_residue: usize,     // Residue index (NH donor)
    pub acceptor_residue: usize,  // Residue index (CO acceptor)
    pub donor_atom_idx: usize,    // Atom index of N in donor residue
    pub acceptor_atom_idx: usize, // Atom index of O in acceptor residue
    pub energy: f32,              // kcal/mol
}

/// DSSP-based secondary structure assigner
///
/// Reference: Kabsch, W., & Sander, C. (1983). Dictionary of protein secondary structure:
/// Pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22(12), 2577–2637.
struct DSSPAssigner {
    energy_threshold: f32, // kcal/mol, typically -0.5
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
                        donor_atom_idx: 0,    // Mapped later
                        acceptor_atom_idx: 0, // Mapped later
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
    fn calculate_hbond_energy(&self, c: Vec3, o: Vec3, n: Vec3, h: Vec3) -> f32 {
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
    fn assign_from_hbonds(&self, hbonds: &[HBond], ss: &mut [SecondaryStructure], n: usize) {
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
                let prev = ss
                    .get(i.saturating_sub(1))
                    .copied()
                    .unwrap_or(SecondaryStructure::Coil);
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
        residues
            .entry(atom.residue_seq)
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
        if let Some(donor_n_idx) = atoms
            .iter()
            .position(|a| a.residue_seq == donor_res_num && a.name == "N")
        {
            hbond.donor_atom_idx = donor_n_idx;
        }

        // O atom in acceptor residue
        if let Some(acceptor_o_idx) = atoms
            .iter()
            .position(|a| a.residue_seq == acceptor_res_num && a.name == "O")
        {
            hbond.acceptor_atom_idx = acceptor_o_idx;
        }
    }

    hbonds
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    #[test]
    fn build_atom_reads_element_and_hetatm() {
        // Minimal mmCIF: an ATOM C-alpha carbon (name "CA", element "C") and a
        // HETATM calcium ion (name "CA", element "Ca"). The element column, not
        // the name, must distinguish them.
        let cif = "\
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA ALA A 1 0.0 0.0 0.0 1.0 0.0
HETATM 2 Ca CA CA B 2 5.0 0.0 0.0 1.0 0.0
";
        let (pdb, _errors) = ReadOptions::default()
            .set_level(StrictnessLevel::Loose)
            .set_format(Format::Mmcif)
            .read_raw(BufReader::new(Cursor::new(cif)))
            .expect("minimal mmCIF should parse");

        let atoms: Vec<Atom> = pdb.atoms_with_hierarchy().map(|h| build_atom(&h)).collect();
        assert_eq!(atoms.len(), 2);

        let carbon = atoms
            .iter()
            .find(|a| !a.is_hetatm)
            .expect("the ATOM record");
        assert_eq!(carbon.name, "CA"); // C-alpha, named CA
        assert_eq!(carbon.element.as_str().to_uppercase(), "C"); // but its element is carbon

        let calcium = atoms
            .iter()
            .find(|a| a.is_hetatm)
            .expect("the HETATM record");
        assert_eq!(calcium.element.as_str().to_uppercase(), "CA"); // calcium, not carbon
    }

    #[test]
    fn truncation_warning_lists_offenders() {
        let mut set = BTreeSet::new();
        assert_eq!(truncation_warning(&set), None);

        set.insert(("chain id", "AAAAA".to_string()));
        let message = truncation_warning(&set).expect("non-empty set warns");
        assert!(message.contains("chain id 'AAAAA'"), "{message}");
        assert!(message.contains("1 identifier"), "{message}");
    }

    #[test]
    fn loader_truncates_overlong_chain_id_without_crashing() {
        // A 5-character chain id exceeds the 4-byte field and must truncate
        // rather than panic; the loader still returns both atoms.
        let cif = "\
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA ALA AAAAA AAAAA 1 0.0 0.0 0.0 1.0 0.0
ATOM 2 C CA ALA AAAAA AAAAA 2 3.8 0.0 0.0 1.0 0.0
";
        let path = std::env::temp_dir().join("pf_truncation_test.cif");
        std::fs::write(&path, cif).unwrap();
        let protein = load_protein_with_options(&path, true, AltlocPolicy::All).unwrap();

        assert_eq!(protein.atoms.len(), 2);
        assert_eq!(protein.atoms[0].chain_id.as_str(), "AAAA"); // truncated from "AAAAA"
        let _ = std::fs::remove_file(&path);
    }
}
