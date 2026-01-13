use std::collections::HashSet;

use anyhow::Result;
use glam::Vec3;

pub mod parser;
pub mod renderer;
pub mod surface;
pub mod network;

pub use surface::SurfacePoint;

use crate::{parser::HBond, renderer::Camera};


pub struct App {
    pub protein: Option<Protein>,
    pub camera: Camera,
    pub redraw_needed: bool,
    pub projected_atom_cache: Option<Vec<renderer::ProjectedAtom>>,
    pub ca_indices: Vec<usize>,
    pub backbone_connections: Vec<(usize, usize)>,
    pub residue_colors: std::collections::HashMap<u32, SecondaryStructure>,
    pub inspect_mode: bool,
    pub residue_highlight: bool,
    pub selected_atom_idx: Option<usize>,
    pub candidate_atoms: Vec<(usize, f32)>, // (atom_idx, distance_along_ray)
    pub candidate_selection_idx: usize, // Index into candidate_atoms
    pub last_canvas_width: f32,
    pub last_canvas_height: f32,
    pub highlighted_atom_indices: HashSet<usize>, // Precomputed set of atoms to highlight
    pub residue_highlight_distance_threshold: f32, // Screen-space distance in pixels
    pub display_mode: DisplayMode,
    pub show_connections: bool,
    pub use_bfactor_colors: bool,
    pub show_surface: bool,
    pub surface_point_density: usize, // Points per atom for surface calculation
    pub show_hydrogen_bonds: bool,
    pub show_hbond_network: bool,
    pub hbond_energy_threshold: f32, // kcal/mol, default -0.5
    pub hbond_graph: Option<network::HBondGraph>,
    pub network_analysis: Option<network::NetworkAnalysis>,
}

impl App {
    pub fn new() -> Self {
        Self {
            protein: None,
            camera: Camera::new(),
            redraw_needed: true,
            projected_atom_cache: None,
            ca_indices: Vec::new(),
            backbone_connections: Vec::new(),
            residue_colors: std::collections::HashMap::new(),
            inspect_mode: false,
            residue_highlight: false,
            selected_atom_idx: None,
            candidate_atoms: Vec::new(),
            candidate_selection_idx: 0,
            last_canvas_width: 0.0,
            last_canvas_height: 0.0,
            highlighted_atom_indices: HashSet::new(),
            residue_highlight_distance_threshold: 50.0, // 50 pixels default
            display_mode: DisplayMode::AllAtoms,
            show_connections: true,
            use_bfactor_colors: false,
            show_surface: false,
            surface_point_density: 100, // Default 100 points per atom
            show_hydrogen_bonds: false,
            show_hbond_network: false,
            hbond_energy_threshold: -0.5, // Default DSSP threshold
            hbond_graph: None,
            network_analysis: None,
        }
    }

    pub fn load_protein(&mut self, path: &str, width: f32, height: f32, skip_surface: bool) -> Result<()> {
        self.protein = Some(parser::load_protein_with_options(path, skip_surface)?);
        
        // Auto-frame the protein when loaded
        if let Some(ref protein) = self.protein {
            renderer::auto_frame_protein(protein, &mut self.camera, width, height);
            
            self.compute_static_geometry();
        }
        
        Ok(())
    }
    
    /// Pre-computes static geometry that doesn't change during rotation
    pub fn compute_static_geometry(&mut self) {
        if let Some(ref protein) = self.protein {
            // Compute CA indices
            self.ca_indices = protein.atoms.iter()
                .enumerate()
                .filter(|(_, atom)| atom.name == "CA")
                .map(|(idx, _)| idx)
                .collect();
            
            // Compute backbone connections
            self.backbone_connections = get_calpha_connections(protein, &self.ca_indices);
            
            // Compute residue colors
            self.residue_colors.clear();
            for atom in &protein.atoms {
                self.residue_colors.entry(atom.residue_seq)
                    .or_insert(atom.secondary_structure);
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DisplayMode {
    AllAtoms,
    Backbone,
}

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
    pub surface_points: Vec<SurfacePoint>,
    pub hbonds: Vec<HBond>,
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
                if i + offset + 4 < residues.len() &&
                   self.can_hbond(&residues[i + offset], &residues[i + offset + 4]) {
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

/// Map B-factor value to RGB color using blue –> cyan –> yellow –> red gradient
/// Normalizes B-factors using percentile-based scaling
pub fn bfactor_to_color(b_factor: f32, b_min: f32, b_max: f32) -> (u8, u8, u8) {
    if b_max <= b_min {
        return (128, 128, 128); // Gray for invalid range
    }
    
    let t = ((b_factor - b_min) / (b_max - b_min)).clamp(0.0, 1.0);
    
    // Blue –> Cyan –> Yellow –> Red gradient
    // Blue (0, 100, 255) at t = 0
    // Cyan (0, 255, 255) at t = 0.33
    // Yellow (255, 255, 0) at t = 0.67
    // Red (255, 0, 0) at t = 1.0
    
    let (r, g, b) = if t < 0.33 {
        // Blue to Cyan
        let local_t = t / 0.33;
        let r = 0.0;
        let g = 100.0 + (255.0 - 100.0) * local_t;
        let b = 255.0;
        (r, g, b)
    } else if t < 0.67 {
        // Cyan to Yellow
        let local_t = (t - 0.33) / 0.34;
        let r = 255.0 * local_t;
        let g = 255.0;
        let b = 255.0 * (1.0 - local_t);
        (r, g, b)
    } else {
        // Yellow to Red
        let local_t = (t - 0.67) / 0.33;
        let r = 255.0;
        let g = 255.0 * (1.0 - local_t);
        let b = 0.0;
        (r, g, b)
    };
    
    (r as u8, g as u8, b as u8)
}

/// Calculate B-factor percentiles for normalization
pub fn calculate_bfactor_range(protein: &Protein) -> (f32, f32) {
    if protein.atoms.is_empty() {
        return (0.0, 100.0);
    }
    
    let mut b_factors: Vec<f32> = protein.atoms.iter().map(|a| a.b_factor).collect();
    b_factors.sort_by(|a, b| a.partial_cmp(b).unwrap());
    
    // Avoiding outliers
    let idx_min = (b_factors.len() as f32 * 0.02) as usize;
    let idx_max = (b_factors.len() as f32 * 0.98) as usize;
    
    let b_min = b_factors[idx_min.min(b_factors.len() - 1)];
    let b_max = b_factors[idx_max.min(b_factors.len() - 1)];
    
    (b_min, b_max)
}

/// Get indices of C-alpha atoms in the protein
pub fn get_calpha_indices(protein: &Protein) -> Vec<usize> {
    protein.atoms.iter()
        .enumerate()
        .filter(|(_, atom)| atom.name == "CA")
        .map(|(idx, _)| idx)
        .collect()
}

/// Get pairs of C-alpha indices that should be connected
/// Only connects sequential C-alphas in the same chain within distance threshold
pub fn get_calpha_connections(protein: &Protein, ca_indices: &[usize]) -> Vec<(usize, usize)> {
    let mut connections = Vec::new();
    let distance_threshold = 4.2; // Angstroms - typical CA-CA distance is ~3.8Å
    
    // Sort C-alphas by chain first, then by residue sequence
    let mut sorted_cas: Vec<(String, u32, usize)> = ca_indices.iter()
        .map(|&idx| {
            let atom = &protein.atoms[idx];
            (atom.chain_id.clone(), atom.residue_seq, idx)
        })
        .collect();
    
    sorted_cas.sort_by(|a, b| {
        a.0.cmp(&b.0).then(a.1.cmp(&b.1))
    });
    
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

/// Draws a line between two points using Bresenham-like algorithm
pub fn draw_line(x0: f32, y0: f32, x1: f32, y1: f32) -> Vec<(f64, f64)> {
    let mut points = Vec::new();
    
    let dx = (x1 - x0).abs();
    let dy = (y1 - y0).abs();
    
    let sx = if x0 < x1 { 1.0 } else { -1.0 };
    let sy = if y0 < y1 { 1.0 } else { -1.0 };
    
    let mut err = dx - dy;
    let mut x = x0;
    let mut y = y0;
    
    // Starting point
    points.push((x as f64, y as f64));
    
    // Draw line using Bresenham's algorithm
    let max_steps = (dx.max(dy) as usize).max(1) + 2; // Safety limit
    
    for _ in 0..max_steps {
        if (x - x1).abs() < 0.5 && (y - y1).abs() < 0.5 {
            // Add endpoint if not already close to last point
            if let Some(&last) = points.last() {
                let dist = ((last.0 - x1 as f64).powi(2) + (last.1 - y1 as f64).powi(2)).sqrt();
                if dist > 1.0 {
                    points.push((x1 as f64, y1 as f64));
                }
            }
            break;
        }
        
        let e2 = 2.0 * err;
        
        // Move in x direction
        if e2 > -dy {
            err -= dy;
            x += sx;
        }
        
        // Move in y direction
        if e2 < dx {
            err += dx;
            y += sy;
        }
        
        points.push((x as f64, y as f64));
    }
    
    points
}

/// Map H-bond energy to color (cyan → yellow → orange)
/// 
/// Weak bonds (~ -0.5 kcal/mol): cyan; 
/// Medium bonds (~ -2.0 kcal/mol): yellow; 
/// Strong bonds (~ -5.0 kcal/mol): orange
pub fn hbond_energy_to_color(energy: f32) -> (u8, u8, u8) {
    // Normalization
    let t = ((energy.abs() - 0.5) / 4.5).clamp(0.0, 1.0);
    
    if t < 0.5 {
        // Cyan (0, 255, 255) to Yellow (255, 255, 0)
        let local_t = t / 0.5;
        let r = (255.0 * local_t) as u8;
        let g = 255;
        let b = (255.0 * (1.0 - local_t)) as u8;
        (r, g, b)
    } else {
        // Yellow (255, 255, 0) to Orange (255, 165, 0)
        let local_t = (t - 0.5) / 0.5;
        let r = 255;
        let g = (255.0 - 90.0 * local_t) as u8;
        let b = 0;
        (r, g, b)
    }
}

/// Draw dashed line by sampling every nth point
pub fn draw_dashed_line(x0: f32, y0: f32, x1: f32, y1: f32, dash_spacing: usize) -> Vec<(f64, f64)> {
    let full_line = draw_line(x0, y0, x1, y1);
    full_line
        .into_iter()
        .enumerate()
        .filter(|(i, _)| i % dash_spacing == 0)
        .map(|(_, p)| p)
        .collect()
}

/// Update the set of highlighted atom indices based on screen-space proximity
/// to the selected atom. Only highlights atoms in the same residue that are
/// within the distance threshold in screen-space pixels.
pub fn update_highlighted_atoms(app: &mut App, width: f32, height: f32) {
    app.highlighted_atom_indices.clear();
    
    let protein = match &app.protein {
        Some(p) => p,
        None => return,
    };
    
    let selected_idx = match app.selected_atom_idx {
        Some(idx) => idx,
        None => return,
    };
    
    if width == 0.0 || height == 0.0 {
        return;
    }
    
    let selected_atom = &protein.atoms[selected_idx];
    let selected_residue_seq = selected_atom.residue_seq;
    let selected_chain_id = &selected_atom.chain_id;
    let selected_position = selected_atom.position;
    
    // Project all atoms to screen space
    if app.camera.cached_view_matrix.is_none() {
        app.camera.get_view_matrix();
    }

    let projected = renderer::project_protein(protein, &app.camera, width, height);
    
    // Get selected atom's screen position
    let selected_screen = &projected[selected_idx];
    let selected_screen_x = selected_screen.x;
    let selected_screen_y = selected_screen.y;
    
    // Find all atoms in the same residue (same chain and residue_seq) within distance thresholds
    for (idx, atom) in protein.atoms.iter().enumerate() {
        // Residue matches both chain_id and residue_seq
        if atom.chain_id == *selected_chain_id && atom.residue_seq == selected_residue_seq {
            // 3D distance check 
            let distance_3d = (atom.position - selected_position).length();
            if distance_3d <= 15.0 { // 15Å threshold
                let proj = &projected[idx];
                let dx = proj.x - selected_screen_x;
                let dy = proj.y - selected_screen_y;
                let screen_distance = (dx * dx + dy * dy).sqrt();
                
                if screen_distance <= app.residue_highlight_distance_threshold {
                    app.highlighted_atom_indices.insert(idx);
                }
            }
        }
    }
}

/// Pick atoms near the click position using screen-space distance
/// Returns a sorted list of (atom_idx, distance_in_pixels) pairs
pub fn pick_atoms_along_ray(
    protein: &Protein,
    camera: &mut Camera,
    click_x: f32,
    click_y: f32,
    width: f32,
    height: f32,
) -> Vec<(usize, f32)> {
    // Project all atoms to screen space
    if camera.cached_view_matrix.is_none() {
        camera.get_view_matrix();
    }

    let projected = renderer::project_protein(protein, camera, width, height);
    
    let click_radius = 10.0; // Base radius in pixels
    
    let mut candidates: Vec<(usize, f32, f32)> = Vec::new(); // (idx, screen_distance, depth)
    
    // Find atoms within click radius in screen space
    for (i, proj_atom) in projected.iter().enumerate() {
        let dx = proj_atom.x - click_x;
        let dy = proj_atom.y - click_y;
        let screen_distance = (dx * dx + dy * dy).sqrt();
        
        if screen_distance <= click_radius {
            candidates.push((i, screen_distance, proj_atom.depth));
        }
    }
    
    // Sort by screen distance first, then by depth (closer to camera wins ties)
    candidates.sort_by(|a, b| {
        a.1.partial_cmp(&b.1)
            .unwrap()
            .then(a.2.partial_cmp(&b.2).unwrap())
    });
    
    // Return (atom_idx, screen_distance) pairs
    candidates.into_iter()
        .map(|(idx, dist, _depth)| (idx, dist))
        .collect()
}
