use std::collections::HashSet;

use anyhow::Result;
use pixelfold_core::{DisplayMode, Protein, SecondaryStructure, get_calpha_connections};
use pixelfold_render::renderer::{self, Camera};

pub mod search;

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
            renderer::auto_frame_protein(protein, &mut self.camera, width, height);

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
    candidates
        .into_iter()
        .map(|(idx, dist, _depth)| (idx, dist))
        .collect()
}
