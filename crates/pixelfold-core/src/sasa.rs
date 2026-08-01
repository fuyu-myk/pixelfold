use crate::structure::Atom;
use nalgebra::Point3;
use rust_sasa::Atom as SasaAtom;

/// Shrake-Rupley algorithm for solvent-accessible surface calculation
pub struct SurfaceCalculator {
    probe_radius: f32,      // Solvent probe radius (typically 1.4Å for water)
    points_per_atom: usize, // Number of test points per atom sphere
}

impl Default for SurfaceCalculator {
    fn default() -> Self {
        Self {
            probe_radius: 1.4,
            points_per_atom: 100,
        }
    }
}

impl SurfaceCalculator {
    /// Create a new surface calculator with custom parameters
    pub fn new(probe_radius: f32, points_per_atom: usize) -> Self {
        Self {
            probe_radius,
            points_per_atom,
        }
    }

    /// Per-atom solvent-accessible surface area (square Angstroms), aligned to
    /// `atoms`. The raw Shrake-Rupley output, reported by `pixelfold sasa`, used
    /// for validation against reference tools, and summed per residue for the
    /// network pane's burial colouring.
    pub fn calculate_atom_sasa(&self, atoms: &[Atom]) -> Vec<f32> {
        if atoms.is_empty() {
            return Vec::new();
        }

        let sasa_atoms = to_sasa_atoms(atoms);
        rust_sasa::calculate_sasa_internal(
            &sasa_atoms,
            self.probe_radius,
            self.points_per_atom,
            true, // parallel
        )
    }
}

/// Convert atoms to the rust-sasa input representation.
fn to_sasa_atoms(atoms: &[Atom]) -> Vec<SasaAtom> {
    atoms
        .iter()
        .enumerate()
        .map(|(idx, atom)| SasaAtom {
            position: Point3::new(atom.position.x, atom.position.y, atom.position.z),
            radius: get_vdw_radius(atom.element.as_str()),
            id: idx,
            parent_id: Some(atom.residue_seq as isize),
        })
        .collect()
}

/// Van der Waals radius (Angstroms) for an element symbol.
///
/// Values from Bondi (1964) J. Phys. Chem. 68: 441-451. Keyed on the element,
/// not the atom name, so a calcium ion ("Ca") and a C-alpha carbon ("C") get
/// their correct, different radii.
pub fn get_vdw_radius(element: &str) -> f32 {
    match element.to_uppercase().as_str() {
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "S" => 1.80,
        "P" => 1.80,
        "H" => 1.20,
        "F" => 1.47,
        "CL" => 1.75,
        "BR" => 1.85,
        "I" => 1.98,
        "FE" => 1.80,
        "ZN" => 1.39,
        "MG" => 1.73,
        "CA" => 1.97, // calcium
        _ => 1.70,    // default to carbon
    }
}

/// Get Kyte-Doolittle hydrophobicity scale value for amino acid residue,
/// Kyte & Doolittle (1982) J. Mol. Biol. 157: 105-132
///
/// Positive values = hydrophobic, Negative values = hydrophilic
///
/// Range: -4.5 (most hydrophilic) to +4.5 (most hydrophobic)
pub fn get_hydrophobicity(residue_name: &str) -> f32 {
    match residue_name.to_uppercase().as_str() {
        // Most hydrophobic
        "ILE" | "I" => 4.5,
        "VAL" | "V" => 4.2,
        "LEU" | "L" => 3.8,
        "PHE" | "F" => 2.8,
        "CYS" | "C" => 2.5,
        "MET" | "M" => 1.9,
        "ALA" | "A" => 1.8,

        // Moderately hydrophobic
        "GLY" | "G" => -0.4,
        "THR" | "T" => -0.7,
        "SER" | "S" => -0.8,
        "TRP" | "W" => -0.9,
        "TYR" | "Y" => -1.3,
        "PRO" | "P" => -1.6,
        "HIS" | "H" => -3.2,

        // Hydrophilic
        "GLU" | "E" => -3.5,
        "GLN" | "Q" => -3.5,
        "ASP" | "D" => -3.5,
        "ASN" | "N" => -3.5,
        "LYS" | "K" => -3.9,

        // Most hydrophilic
        "ARG" | "R" => -4.5,

        // Unknown/non-standard
        _ => 0.0,
    }
}

/// Convert hydrophobicity value to RGB color
///
/// Uses a gradient: Blue (hydrophilic) -> White (neutral) -> Orange (hydrophobic)
pub fn hydrophobicity_to_color(hydrophobicity: f32) -> (u8, u8, u8) {
    // Normalization (Kyte-Doolittle range is -4.5 to +4.5)
    let min_hydro = -4.5;
    let max_hydro = 4.5;

    let t = ((hydrophobicity - min_hydro) / (max_hydro - min_hydro)).clamp(0.0, 1.0);

    // Blue (hydrophilic) -> White (neutral) -> Orange (hydrophobic)
    let (r, g, b) = if t < 0.5 {
        // Blue to White (hydrophilic to neutral)
        let local_t = t * 2.0; // Map [0, 0.5] to [0, 1]
        let r = (0.0 + 255.0 * local_t) as u8;
        let g = (100.0 + 155.0 * local_t) as u8;
        let b = 255;
        (r, g, b)
    } else {
        // White to Orange (neutral to hydrophobic)
        let local_t = (t - 0.5) * 2.0;
        let r = 255;
        let g = (255.0 - 100.0 * local_t) as u8;
        let b = (255.0 - 255.0 * local_t) as u8;
        (r, g, b)
    };

    (r, g, b)
}

/// Theoretical maximum solvent-accessible surface area (Å²) of each standard
/// residue, from Tien et al. 2013 (Gly-X-Gly tripeptides, all conformations
/// enumerated). Relative solvent accessibility is a residue's SASA over this.
///
/// Source: Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013),
/// "Maximum allowed solvent accessibilities of residues in proteins",
/// PLoS ONE 8(11): e80635, theoretical column.
pub fn max_asa(residue_name: &str) -> Option<f32> {
    let max = match residue_name {
        "ALA" => 129.0,
        "ARG" => 274.0,
        "ASN" => 195.0,
        "ASP" => 193.0,
        "CYS" => 167.0,
        "GLN" => 225.0,
        "GLU" => 223.0,
        "GLY" => 104.0,
        "HIS" => 224.0,
        "ILE" => 197.0,
        "LEU" => 201.0,
        "LYS" => 236.0,
        "MET" => 224.0,
        "PHE" => 240.0,
        "PRO" => 159.0,
        "SER" => 155.0,
        "THR" => 172.0,
        "TRP" => 285.0,
        "TYR" => 263.0,
        "VAL" => 174.0,
        _ => return None,
    };
    Some(max)
}

/// Relative solvent accessibility: absolute SASA over the residue's theoretical
/// maximum. `None` for residues with no reference (ligands, modified or
/// non-standard residues). Not clamped: a residue more exposed than the
/// reference tripeptide can read above 1.
pub fn relative_sasa(absolute: f32, residue_name: &str) -> Option<f32> {
    max_asa(residue_name).map(|max| absolute / max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_asa_covers_the_twenty_standard_residues() {
        for residue in [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        ] {
            assert!(
                max_asa(residue).is_some(),
                "{residue} missing from the table"
            );
        }
        assert_eq!(max_asa("TRP"), Some(285.0));
        assert_eq!(max_asa("GLY"), Some(104.0));
        assert_eq!(max_asa("STI"), None, "a ligand has no reference maximum");
    }

    #[test]
    fn relative_sasa_normalises_against_the_maximum() {
        // A fully exposed alanine (max 129 Å²) is 1.0; half that is 0.5.
        assert!((relative_sasa(129.0, "ALA").unwrap() - 1.0).abs() < 1e-6);
        assert!((relative_sasa(64.5, "ALA").unwrap() - 0.5).abs() < 1e-6);
        assert!(relative_sasa(50.0, "STI").is_none());
    }

    #[test]
    fn test_vdw_radii() {
        assert_eq!(get_vdw_radius("C"), 1.70);
        assert_eq!(get_vdw_radius("N"), 1.55);
        assert_eq!(get_vdw_radius("O"), 1.52);
        assert_eq!(get_vdw_radius("S"), 1.80);
    }

    #[test]
    fn vdw_distinguishes_calcium_from_carbon() {
        // The element column, not the atom name, drives the radius: a calcium
        // ion (element "Ca") and a C-alpha carbon (element "C") differ.
        assert_eq!(get_vdw_radius("Ca"), 1.97);
        assert_eq!(get_vdw_radius("C"), 1.70);
        assert_ne!(get_vdw_radius("Ca"), get_vdw_radius("C"));
    }

    #[test]
    fn test_hydrophobicity_scale() {
        // Most hydrophobic
        assert_eq!(get_hydrophobicity("ILE"), 4.5);
        assert_eq!(get_hydrophobicity("I"), 4.5);

        // Most hydrophilic
        assert_eq!(get_hydrophobicity("ARG"), -4.5);
        assert_eq!(get_hydrophobicity("R"), -4.5);

        // Neutral-ish
        assert_eq!(get_hydrophobicity("GLY"), -0.4);

        // Unknown
        assert_eq!(get_hydrophobicity("XYZ"), 0.0);
    }

    #[test]
    fn test_hydrophobicity_colors() {
        // Hydrophilic (blue-ish)
        let (r, g, b) = hydrophobicity_to_color(-4.5);
        assert!(b > r && b > g, "Hydrophilic should be blue-ish");

        // Hydrophobic (orange-ish)
        let (r, g, b) = hydrophobicity_to_color(4.5);
        assert!(r == 255 && g > b, "Hydrophobic should be orange-ish");

        // Neutral (white-ish)
        let (r, g, b) = hydrophobicity_to_color(0.0);
        assert!(r > 200 && g > 200 && b > 200, "Neutral should be white-ish");
    }
}
