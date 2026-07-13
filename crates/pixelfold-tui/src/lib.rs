use std::collections::HashSet;

use anyhow::Result;
use pixelfold_core::get_calpha_connections;

pub mod search;
pub mod visualization;

pub use pixelfold_core::sasa::SurfacePoint;
pub use pixelfold_core::{Atom, DisplayMode, Protein, SecondaryStructure, get_calpha_indices};

use crate::visualization::renderer::Camera;

pub struct App {
    pub protein: Option<Protein>,
    pub camera: Camera,
    pub redraw_needed: bool,
    pub projected_atom_cache: Option<Vec<visualization::renderer::ProjectedAtom>>,
    pub ca_indices: Vec<usize>,
    pub backbone_connections: Vec<(usize, usize)>,
    pub residue_colors: std::collections::HashMap<u32, SecondaryStructure>,
    pub inspect_mode: bool,
    pub residue_highlight: bool,
    pub selected_atom_idx: Option<usize>,
    pub candidate_atoms: Vec<(usize, f32)>, // (atom_idx, distance_along_ray)
    pub candidate_selection_idx: usize,     // Index into candidate_atoms
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
    pub hbond_graph: Option<pixelfold_core::rin::HBondGraph>,
    pub network_analysis: Option<pixelfold_core::rin::NetworkAnalysis>,
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

    pub fn load_protein(
        &mut self,
        path: &str,
        width: f32,
        height: f32,
        skip_surface: bool,
    ) -> Result<()> {
        self.protein = Some(pixelfold_core::parser::load_protein_with_options(
            path,
            skip_surface,
        )?);

        // Auto-frame the protein when loaded
        if let Some(ref protein) = self.protein {
            visualization::renderer::auto_frame_protein(protein, &mut self.camera, width, height);

            self.compute_static_geometry();
        }

        Ok(())
    }

    /// Pre-computes static geometry that doesn't change during rotation
    pub fn compute_static_geometry(&mut self) {
        if let Some(ref protein) = self.protein {
            // Compute CA indices
            self.ca_indices = protein
                .atoms
                .iter()
                .enumerate()
                .filter(|(_, atom)| atom.name == "CA")
                .map(|(idx, _)| idx)
                .collect();

            // Compute backbone connections
            self.backbone_connections = get_calpha_connections(protein, &self.ca_indices);

            // Compute residue colors
            self.residue_colors.clear();
            for atom in &protein.atoms {
                self.residue_colors
                    .entry(atom.residue_seq)
                    .or_insert(atom.secondary_structure);
            }
        }
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
pub fn draw_dashed_line(
    x0: f32,
    y0: f32,
    x1: f32,
    y1: f32,
    dash_spacing: usize,
) -> Vec<(f64, f64)> {
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

    let projected = visualization::renderer::project_protein(protein, &app.camera, width, height);

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
            if distance_3d <= 15.0 {
                // 15Å threshold
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

    let projected = visualization::renderer::project_protein(protein, camera, width, height);

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
    candidates
        .into_iter()
        .map(|(idx, dist, _depth)| (idx, dist))
        .collect()
}
