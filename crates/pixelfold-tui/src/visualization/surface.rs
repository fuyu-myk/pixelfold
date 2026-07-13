use glam::Vec3;
use nalgebra::Point3;
use rust_sasa::Atom as SasaAtom;
use crate::Atom;


/// A point on the solvent-accessible surface
#[derive(Clone, Debug)]
pub struct SurfacePoint {
    pub position: Vec3,
    pub atom_idx: usize,        // Parent atom index
    pub hydrophobicity: f32,    // Kyte-Doolittle scale value
}

/// Shrake-Rupley algorithm for solvent-accessible surface calculation
pub struct SurfaceCalculator {
    probe_radius: f32,          // Solvent probe radius (typically 1.4Å for water)
    points_per_atom: usize,     // Number of test points per atom sphere
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

    /// Calculate solvent-accessible surface points for a protein
    pub fn calculate_surface(&self, atoms: &[Atom]) -> Vec<SurfacePoint> {
        if atoms.is_empty() {
            return Vec::new();
        }

        // Convert to rust-sasa format
        let sasa_atoms: Vec<SasaAtom> = atoms
            .iter()
            .enumerate()
            .map(|(idx, atom)| {
                let vdw_radius = get_vdw_radius(&atom.name);
                SasaAtom {
                    position: Point3::new(atom.position.x, atom.position.y, atom.position.z),
                    radius: vdw_radius,
                    id: idx,
                    parent_id: Some(atom.residue_seq as isize),
                }
            })
            .collect();

        let sasa_values = rust_sasa::calculate_sasa_internal(
            &sasa_atoms,
            self.probe_radius,
            self.points_per_atom,
            true  // Parallel computation for better performance
        );
        
        // Generate surface points based on SASA values
        let mut surface_points = Vec::new();
        let sphere_points = self.generate_fibonacci_sphere(self.points_per_atom);
        
        for (atom_idx, atom) in atoms.iter().enumerate() {
            let sasa_value = sasa_values.get(atom_idx).copied().unwrap_or(0.0);
            
            // Only generate surface points for atoms with accessible surface
            if sasa_value > 0.1 {  // Threshold (avoids atoms with tiny SASA)
                let vdw_radius = get_vdw_radius(&atom.name);
                let surface_radius = vdw_radius + self.probe_radius;
                let hydrophobicity = get_hydrophobicity(&atom.residue_name);
                
                // Calculate how many points to generate based on SASA ratio
                // Full sphere: 4π(r + probe)²
                let full_surface = 4.0 * std::f32::consts::PI * surface_radius * surface_radius;
                let accessibility_ratio = (sasa_value / full_surface).min(1.0);
                let num_points = (self.points_per_atom as f32 * accessibility_ratio).max(1.0) as usize;
                
                for &unit_point in sphere_points.iter().take(num_points) {
                    let surface_pos = atom.position + unit_point * surface_radius;
                    surface_points.push(SurfacePoint {
                        position: surface_pos,
                        atom_idx,
                        hydrophobicity,
                    });
                }
            }
        }

        surface_points
    }

    /// Generate uniformly distributed points on a unit sphere using Fibonacci spiral
    fn generate_fibonacci_sphere(&self, num_points: usize) -> Vec<Vec3> {
        let mut points = Vec::with_capacity(num_points);
        let golden_ratio = (1.0 + 5.0_f32.sqrt()) / 2.0;
        let angle_increment = std::f32::consts::PI * 2.0 * golden_ratio;

        for i in 0..num_points {
            let t = i as f32 / num_points as f32;
            let inclination = (1.0 - 2.0 * t).acos();
            let azimuth = angle_increment * i as f32;

            let x = inclination.sin() * azimuth.cos();
            let y = inclination.sin() * azimuth.sin();
            let z = inclination.cos();

            points.push(Vec3::new(x, y, z));
        }

        points
    }
}

/// Get van der Waals radius for an atom based on element
/// 
/// Values are obtained from Bondi (1964) J. Phys. Chem. 68: 441-451
pub fn get_vdw_radius(atom_name: &str) -> f32 {
    // Extract element from atom name (first 1-2 characters, trimmed)
    let element = atom_name.trim_start_matches(|c: char| c.is_numeric())
        .chars()
        .take(2)
        .collect::<String>()
        .trim()
        .to_uppercase();

    match element.as_str() {
        "CA" if atom_name.starts_with("CA") && atom_name.len() > 2 => 1.97,  // Calcium ion
        "C" | "CA" | "CB" | "CG" | "CD" | "CE" | "CZ" => 1.70,  // Carbon
        "N" | "NE" | "NH" | "NZ" | "ND" => 1.55,                // Nitrogen
        "O" | "OD" | "OE" | "OG" | "OH" | "OXT" => 1.52,        // Oxygen
        "S" | "SD" | "SG" => 1.80,                              // Sulfur
        "P" => 1.80,                                            // Phosphorus
        "H" => 1.20,                                            // Hydrogen
        "F" => 1.47,                                            // Fluorine
        "CL" => 1.75,                                           // Chlorine
        "BR" => 1.85,                                           // Bromine
        "I" => 1.98,                                            // Iodine
        "FE" => 1.80,                                           // Iron
        "ZN" => 1.39,                                           // Zinc
        "MG" => 1.73,                                           // Magnesium
        _ => 1.70,                                              // Default to carbon
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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fibonacci_sphere_generation() {
        let calc = SurfaceCalculator::default();
        let points = calc.generate_fibonacci_sphere(100);
        
        assert_eq!(points.len(), 100);
        
        // All points should be approximately unit length
        for point in points {
            let length = point.length();
            assert!((length - 1.0).abs() < 0.01, "Point length {} not close to 1.0", length);
        }
    }

    #[test]
    fn test_vdw_radii() {
        assert_eq!(get_vdw_radius("C"), 1.70);
        assert_eq!(get_vdw_radius("CA"), 1.70); // C-alpha carbon
        assert_eq!(get_vdw_radius("N"), 1.55);
        assert_eq!(get_vdw_radius("O"), 1.52);
        assert_eq!(get_vdw_radius("S"), 1.80);
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
