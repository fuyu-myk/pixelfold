use anyhow::Result;
use crossterm::event::{self, KeyCode, KeyModifiers};
use pixelfold::{
    parser,
    renderer::{self, Camera},
    surface,
    network,
    Protein,
    SecondaryStructure,
};
use ratatui::prelude::*;
use ratatui::widgets::canvas::{Canvas, Points};
use std::collections::HashSet;
use std::env;


#[derive(Clone, Copy, Debug, PartialEq)]
enum DisplayMode {
    AllAtoms,
    Backbone,
}

struct App {
    protein: Option<Protein>,
    camera: Camera,
    inspect_mode: bool,
    residue_highlight: bool,
    selected_atom_idx: Option<usize>,
    candidate_atoms: Vec<(usize, f32)>, // (atom_idx, distance_along_ray)
    candidate_selection_idx: usize, // Index into candidate_atoms
    last_canvas_width: f32,
    last_canvas_height: f32,
    highlighted_atom_indices: HashSet<usize>, // Precomputed set of atoms to highlight
    residue_highlight_distance_threshold: f32, // Screen-space distance in pixels
    display_mode: DisplayMode,
    show_connections: bool,
    use_bfactor_colors: bool,
    show_surface: bool,
    surface_point_density: usize, // Points per atom for surface calculation
    show_hydrogen_bonds: bool,
    show_hbond_network: bool,
    hbond_energy_threshold: f32, // kcal/mol, default -0.5
    hbond_graph: Option<network::HBondGraph>,
    network_analysis: Option<network::NetworkAnalysis>,
}

impl App {
    fn new() -> Self {
        Self {
            protein: None,
            camera: Camera::new(),
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

    fn load_protein(&mut self, path: &str, width: f32, height: f32, skip_surface: bool) -> Result<()> {
        self.protein = Some(parser::load_protein_with_options(path, skip_surface)?);
        
        // Auto-frame the protein when loaded
        if let Some(ref protein) = self.protein {
            renderer::auto_frame_protein(protein, &mut self.camera, width, height);
        }
        
        Ok(())
    }
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    
    // Parse command-line arguments
    let mut protein_path: Option<String> = None;
    let mut skip_surface = false;
    
    for arg in args.iter().skip(1) {
        if arg == "--no-surface" {
            skip_surface = true;
        } else if arg == "--help" || arg == "-h" {
            println!("PixelFold - Terminal-based 3D protein structure viewer");
            println!();
            println!("Usage: pixelfold <path/to/protein.pdb|protein.cif> [OPTIONS]");
            println!();
            println!("Options:");
            println!("  --no-surface    Skip surface calculation for faster loading (large proteins)");
            println!("  --help, -h      Show this help message");
            println!();
            println!("Controls:");
            println!("  WASDZX    Rotate structure");
            println!("  +/-       Zoom in/out");
            println!("  Arrows    Pan view (or adjust settings when mode active)");
            println!("  1         Show all atoms");
            println!("  2         Show backbone atoms");
            println!("  V         Toggle surface visualization");
            println!("    ↑↓      Adjust surface density when surface visible");
            println!("  C         Toggle backbone connections");
            println!("  B         Toggle B-factor coloring");
            println!("  H         Toggle hydrogen bond display");
            println!("    ↑↓      Adjust H-bond energy threshold when H-bonds visible");
            println!("  N         Toggle H-bond network analysis overlay");
            println!("  F         Frame the protein in view");
            println!("  I         Inspect mode");
            println!("    ↑↓      Cycle through nearby atoms");
            println!("  R         Toggle residue highlighting (in inspect mode)");
            println!("  Q         Quit");
            return Ok(());
        } else if !arg.starts_with('-') {
            protein_path = Some(arg.clone());
        }
    }
    
    let mut app = App::new();
    let mut terminal = ratatui::init();
    
    crossterm::execute!(
        std::io::stdout(),
        crossterm::event::EnableMouseCapture
    )?;
    
    if let Some(path) = protein_path {
        // Initial terminal size
        let size = terminal.size()?;
        let width = size.width as f32 * 2.0;  // Braille canvas width
        let height = size.height as f32 * 4.0; // Braille canvas height
        
        match app.load_protein(&path, width, height, skip_surface) {
            Ok(_) => {},
            Err(e) => eprintln!("Failed to load protein: {}", e),
        }
    }

    let result = run_app(&mut terminal, &mut app);
    
    crossterm::execute!(
        std::io::stdout(),
        crossterm::event::DisableMouseCapture
    )?;
    
    ratatui::restore();
    result
}

fn run_app(terminal: &mut ratatui::Terminal<impl ratatui::backend::Backend>, app: &mut App) -> Result<()> {
    loop {
        let size = terminal.size()?;
        let canvas_width = size.width as f32 * 2.0;
        let canvas_height = size.height as f32 * 4.0;
        
        terminal.draw(|frame| ui(frame, app))?;

        if event::poll(std::time::Duration::from_millis(16))? {
            match event::read()? {
                event::Event::Key(key) => {
                    handle_input(app, key.code, key.modifiers, canvas_width, canvas_height)?;

                    if matches!(key.code, KeyCode::Char('q')) {
                        break;
                    }
                }
                event::Event::Mouse(mouse) => {
                    let area = Rect {
                        x: 0,
                        y: 0,
                        width: size.width,
                        height: size.height,
                    };
                    handle_mouse(app, mouse, area)?;
                }
                _ => {}
            }
        }
    }
    Ok(())
}

/// Map B-factor value to RGB color using blue –> cyan –> yellow –> red gradient
/// Normalizes B-factors using percentile-based scaling
fn bfactor_to_color(b_factor: f32, b_min: f32, b_max: f32) -> (u8, u8, u8) {
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
fn calculate_bfactor_range(protein: &Protein) -> (f32, f32) {
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
fn get_calpha_indices(protein: &Protein) -> Vec<usize> {
    protein.atoms.iter()
        .enumerate()
        .filter(|(_, atom)| atom.name == "CA")
        .map(|(idx, _)| idx)
        .collect()
}

/// Get pairs of C-alpha indices that should be connected
/// Only connects sequential C-alphas in the same chain within distance threshold
fn get_calpha_connections(protein: &Protein, ca_indices: &[usize]) -> Vec<(usize, usize)> {
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
fn draw_line(x0: f32, y0: f32, x1: f32, y1: f32) -> Vec<(f64, f64)> {
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
fn hbond_energy_to_color(energy: f32) -> (u8, u8, u8) {
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
fn draw_dashed_line(x0: f32, y0: f32, x1: f32, y1: f32, dash_spacing: usize) -> Vec<(f64, f64)> {
    let full_line = draw_line(x0, y0, x1, y1);
    full_line
        .into_iter()
        .enumerate()
        .filter(|(i, _)| i % dash_spacing == 0)
        .map(|(_, p)| p)
        .collect()
}

fn handle_input(app: &mut App, key: KeyCode, _modifiers: KeyModifiers, width: f32, height: f32) -> Result<()> {
    let rotation_speed = 0.1;
    let zoom_speed = 0.1;
    let pan_speed = 5.0;

    match key {
        KeyCode::Char('q') => {}
        
        // Rotation controls (WASD)
        KeyCode::Char('w') => app.camera.rotate(rotation_speed, 0.0, 0.0),
        KeyCode::Char('s') => app.camera.rotate(-rotation_speed, 0.0, 0.0),
        KeyCode::Char('a') => app.camera.rotate(0.0, rotation_speed, 0.0),
        KeyCode::Char('d') => app.camera.rotate(0.0, -rotation_speed, 0.0),
        KeyCode::Char('z') => app.camera.rotate(0.0, 0.0, rotation_speed),
        KeyCode::Char('x') => app.camera.rotate(0.0, 0.0, -rotation_speed),
        
        // Zoom controls
        KeyCode::Char('+') | KeyCode::Char('=') => app.camera.adjust_zoom(zoom_speed),
        KeyCode::Char('-') | KeyCode::Char('_') => app.camera.adjust_zoom(-zoom_speed),
        
        // Pan controls (Arrow keys)
        // Cycle through candidates in inspect mode (up and down)
        // Adjust surface density (up and down)
        // Adjust H-bond energy threshold (up and down when H-bonds visible)
        KeyCode::Up => {
            if app.show_hydrogen_bonds && !app.show_surface && !app.inspect_mode {
                // More negative = stronger bonds only
                app.hbond_energy_threshold = (app.hbond_energy_threshold - 0.5).max(-10.0);
                
                // Recompute network analysis if network mode is active
                if app.show_hbond_network {
                    if let Some(ref graph) = app.hbond_graph {
                        let filtered = graph.filter_by_energy(app.hbond_energy_threshold);
                        app.network_analysis = Some(filtered.analyze());
                    }
                }
            } else if app.show_surface {
                app.surface_point_density = (app.surface_point_density + 25).min(500);
                if let Some(ref mut protein) = app.protein {
                    let surface_calculator = surface::SurfaceCalculator::new(1.4, app.surface_point_density);
                    protein.surface_points = surface_calculator.calculate_surface(&protein.atoms);
                }
            } else if app.inspect_mode && !app.candidate_atoms.is_empty() {
                if app.candidate_selection_idx == 0 {
                    app.candidate_selection_idx = app.candidate_atoms.len() - 1;
                } else {
                    app.candidate_selection_idx -= 1;
                }
                app.selected_atom_idx = Some(app.candidate_atoms[app.candidate_selection_idx].0);
                update_highlighted_atoms(app, width, height);
            } else {
                app.camera.pan_camera(0.0, -pan_speed);
            }
        }
        KeyCode::Down => {
            if app.show_hydrogen_bonds && !app.show_surface && !app.inspect_mode {
                // Less negative = more bonds shown
                app.hbond_energy_threshold = (app.hbond_energy_threshold + 0.5).min(-0.1);
                
                // Recompute network analysis if network mode is active
                if app.show_hbond_network {
                    if let Some(ref graph) = app.hbond_graph {
                        let filtered = graph.filter_by_energy(app.hbond_energy_threshold);
                        app.network_analysis = Some(filtered.analyze());
                    }
                }
            } else if app.show_surface {
                app.surface_point_density = (app.surface_point_density.saturating_sub(25)).max(100);
                if let Some(ref mut protein) = app.protein {
                    let surface_calculator = surface::SurfaceCalculator::new(1.4, app.surface_point_density);
                    protein.surface_points = surface_calculator.calculate_surface(&protein.atoms);
                }
            } else if app.inspect_mode && !app.candidate_atoms.is_empty() {
                app.candidate_selection_idx = (app.candidate_selection_idx + 1) % app.candidate_atoms.len();
                app.selected_atom_idx = Some(app.candidate_atoms[app.candidate_selection_idx].0);
                update_highlighted_atoms(app, width, height);
            } else {
                app.camera.pan_camera(0.0, pan_speed);
            }
        }
        KeyCode::Left => app.camera.pan_camera(pan_speed, 0.0),
        KeyCode::Right => app.camera.pan_camera(-pan_speed, 0.0),

        // Inspect mode
        KeyCode::Char('i') => {
            app.inspect_mode = !app.inspect_mode;

            if !app.inspect_mode {
                app.selected_atom_idx = None;
                app.candidate_atoms.clear();
                app.candidate_selection_idx = 0;
                app.residue_highlight = false;
            }
        }
        KeyCode::Char('r') => {
            if app.inspect_mode {
                app.residue_highlight = !app.residue_highlight;
            }
        }
        
        // Reset view
        KeyCode::Char('f') => {
            if let Some(ref protein) = app.protein {
                renderer::auto_frame_protein(protein, &mut app.camera, width, height);
            }
        }
        
        // Display mode controls
        KeyCode::Char('1') => {
            app.display_mode = DisplayMode::AllAtoms;
        }
        KeyCode::Char('2') => {
            app.display_mode = DisplayMode::Backbone;
        }
        
        // Toggle backbone connections
        KeyCode::Char('c') => {
            app.show_connections = !app.show_connections;
        }
        
        // Toggle B-factor coloring
        KeyCode::Char('b') => {
            app.use_bfactor_colors = !app.use_bfactor_colors;
        }
        
        // Toggle surface visualization
        KeyCode::Char('v') => {
            app.show_surface = !app.show_surface;
            
            // Compute surface on-demand if not already computed
            if app.show_surface {
                if let Some(ref mut protein) = app.protein {
                    if protein.surface_points.is_empty() {
                        let surface_calculator = surface::SurfaceCalculator::new(1.4, app.surface_point_density);
                        protein.surface_points = surface_calculator.calculate_surface(&protein.atoms);
                    }
                }
            }
        }

        // Toggle Hydrogen bond display
        KeyCode::Char('h') => {
            app.show_hydrogen_bonds = !app.show_hydrogen_bonds;

            // Build H-bond graph on first activation
            if app.show_hydrogen_bonds && app.hbond_graph.is_none() {
                if let Some(ref protein) = app.protein {
                    app.hbond_graph = Some(network::HBondGraph::build(protein));
                }
            }
        }
        
        // Toggle H-bond network visualization overlay
        KeyCode::Char('n') => {
            app.show_hbond_network = !app.show_hbond_network;
            
            // Compute network analysis on first activation
            if app.show_hbond_network {
                if app.hbond_graph.is_none() {
                    if let Some(ref protein) = app.protein {
                        app.hbond_graph = Some(network::HBondGraph::build(protein));
                    }
                }
                
                if let Some(ref graph) = app.hbond_graph {
                    let filtered = graph.filter_by_energy(app.hbond_energy_threshold);
                    app.network_analysis = Some(filtered.analyze());
                }
            }
        }
        
        _ => {}
    }

    Ok(())
}

fn handle_mouse(app: &mut App, mouse: event::MouseEvent, _terminal_size: Rect) -> Result<()> {
    use event::MouseEventKind;
    
    if !app.inspect_mode {
        return Ok(());
    }
    
    match mouse.kind {
        MouseEventKind::Down(event::MouseButton::Left) => {
            // Use the canvas dimensions from the last render
            let canvas_width = app.last_canvas_width;
            let canvas_height = app.last_canvas_height;
            
            if canvas_width == 0.0 || canvas_height == 0.0 {
                return Ok(()); // Not yet rendered
            }
            
            // Convert terminal coordinates to canvas coordinates
            let click_x = mouse.column as f32 * 2.0;
            let click_y = canvas_height - (mouse.row as f32 * 4.0);
            
            if let Some(ref protein) = app.protein {
                let candidates = pick_atoms_along_ray(
                    protein,
                    &app.camera,
                    click_x,
                    click_y,
                    canvas_width,
                    canvas_height,
                );
                
                if !candidates.is_empty() {
                    app.candidate_atoms = candidates.into_iter().take(5).collect(); // Top 5 candidates
                    app.candidate_selection_idx = 0;
                    app.selected_atom_idx = Some(app.candidate_atoms[0].0);
                    update_highlighted_atoms(app, canvas_width, canvas_height);
                } else {
                    app.selected_atom_idx = None;
                    app.candidate_atoms.clear();
                    app.candidate_selection_idx = 0;
                    app.highlighted_atom_indices.clear();
                }
            }
        }
        _ => {}
    }
    
    Ok(())
}

/// Update the set of highlighted atom indices based on screen-space proximity
/// to the selected atom. Only highlights atoms in the same residue that are
/// within the distance threshold in screen-space pixels.
fn update_highlighted_atoms(app: &mut App, width: f32, height: f32) {
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
fn pick_atoms_along_ray(
    protein: &Protein,
    camera: &Camera,
    click_x: f32,
    click_y: f32,
    width: f32,
    height: f32,
) -> Vec<(usize, f32)> {
    // Project all atoms to screen space
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

fn ui(frame: &mut Frame, app: &mut App) {
    let area = frame.area();
    
    if let Some(ref protein) = app.protein {
        // Split layout if atom is selected (left: main view, right: info panel)
        let main_area = if app.selected_atom_idx.is_some() {
            let chunks = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([
                    Constraint::Percentage(75),
                    Constraint::Percentage(25),
                ])
                .split(area);
            chunks[0]
        } else {
            area
        };
        
        // Store canvas dimensions for click detection
        app.last_canvas_width = main_area.width as f32 * 2.0;
        app.last_canvas_height = main_area.height as f32 * 4.0;
        
        // For highlighting residues
        let selected_residue_seq = app.selected_atom_idx
            .map(|idx| protein.atoms[idx].residue_seq);
        
        let canvas = Canvas::default()
            .marker(ratatui::symbols::Marker::Braille)
            .x_bounds([0.0, main_area.width as f64 * 2.0])
            .y_bounds([0.0, main_area.height as f64 * 4.0]) // Braille: 2 * 4 pixels
            .paint(|ctx| {
                // Project all atoms to 2D
                let width = main_area.width as f32 * 2.0;
                let height = main_area.height as f32 * 4.0;
                
                // Determine which atoms to display based on display mode
                let display_indices: Vec<usize> = match app.display_mode {
                    DisplayMode::AllAtoms => (0..protein.atoms.len()).collect(),
                    DisplayMode::Backbone => get_calpha_indices(protein),
                };
                
                // Calculate B-factor range if using B-factor coloring
                let (b_min, b_max) = if app.use_bfactor_colors {
                    calculate_bfactor_range(protein)
                } else {
                    (0.0, 100.0)
                };
                
                let mut residue_colors: std::collections::HashMap<u32, SecondaryStructure> = 
                    std::collections::HashMap::new();
                
                // Assign colors based on original atom order
                for atom in &protein.atoms {
                    residue_colors.entry(atom.residue_seq)
                        .or_insert(atom.secondary_structure);
                }
                
                let all_projected: Vec<renderer::ProjectedAtom> = 
                    renderer::project_protein(protein, &app.camera, width, height);
                
                // Filter and prepare atoms for rendering
                let mut projected_with_idx: Vec<(usize, &renderer::ProjectedAtom)> = 
                    display_indices.iter()
                        .map(|&idx| (idx, &all_projected[idx]))
                        .collect();
                
                // Sort by depth (painter's algorithm)
                projected_with_idx.sort_by(|a, b| a.1.depth.partial_cmp(&b.1.depth).unwrap());
                
                let selected_atom_idx = app.selected_atom_idx;
                let residue_highlight_active = app.residue_highlight && selected_residue_seq.is_some();
                
                // Draw backbone connections first (if enabled and in backbone mode)
                if app.show_connections && !app.residue_highlight && app.display_mode == DisplayMode::Backbone {
                    let ca_indices = get_calpha_indices(protein);
                    let connections = get_calpha_connections(protein, &ca_indices);
                    
                    for (idx1, idx2) in connections {
                        let proj1 = &all_projected[idx1];
                        let proj2 = &all_projected[idx2];
                        
                        // Only draw if both atoms are on screen
                        if proj1.x >= 0.0 && proj1.x <= width && 
                           proj1.y >= 0.0 && proj1.y <= height &&
                           proj2.x >= 0.0 && proj2.x <= width && 
                           proj2.y >= 0.0 && proj2.y <= height {
                            
                            let line_points = draw_line(proj1.x, proj1.y, proj2.x, proj2.y);
                            
                            // Color the line based on coloring mode (av of two atoms)
                            let color = if app.use_bfactor_colors {
                                let b1 = protein.atoms[idx1].b_factor;
                                let b2 = protein.atoms[idx2].b_factor;
                                let avg_b = (b1 + b2) / 2.0;
                                let (r, g, b) = bfactor_to_color(avg_b, b_min, b_max);
                                Color::Rgb(r, g, b)
                            } else {
                                // Use secondary structure color (slightly dimmed for lines)
                                let residue_seq = protein.atoms[idx1].residue_seq;
                                if let Some(&ss) = residue_colors.get(&residue_seq) {
                                    let (r, g, b) = ss.color_rgb();
                                    Color::Rgb(
                                        (r as f32 * 0.8) as u8,
                                        (g as f32 * 0.8) as u8,
                                        (b as f32 * 0.8) as u8,
                                    )
                                } else {
                                    Color::Gray
                                }
                            };
                            
                            ctx.draw(&Points {
                                coords: &line_points,
                                color,
                            });
                        }
                    }
                }
                
                // Draw hydrogen bonds (if enabled)
                if app.show_hydrogen_bonds {
                    // Filter H-bonds by energy threshold
                    let visible_hbonds: Vec<&pixelfold::parser::HBond> = protein.hbonds
                        .iter()
                        .filter(|hb| hb.energy < app.hbond_energy_threshold)
                        .collect();
                    
                    for hbond in visible_hbonds {
                        let donor_idx = hbond.donor_atom_idx;
                        let acceptor_idx = hbond.acceptor_atom_idx;
                        
                        // Skip if indices are out of bounds
                        if donor_idx >= all_projected.len() || acceptor_idx >= all_projected.len() {
                            continue;
                        }
                        
                        let proj_donor = &all_projected[donor_idx];
                        let proj_acceptor = &all_projected[acceptor_idx];
                        
                        // Only draw if both atoms are on screen
                        if proj_donor.x >= 0.0 && proj_donor.x <= width && 
                           proj_donor.y >= 0.0 && proj_donor.y <= height &&
                           proj_acceptor.x >= 0.0 && proj_acceptor.x <= width && 
                           proj_acceptor.y >= 0.0 && proj_acceptor.y <= height {
                            
                            let line_points = draw_dashed_line(
                                proj_donor.x, 
                                proj_donor.y, 
                                proj_acceptor.x, 
                                proj_acceptor.y,
                                3  // Dash spacing
                            );
                            
                            let (r, g, b) = hbond_energy_to_color(hbond.energy);
                            
                            ctx.draw(&Points {
                                coords: &line_points,
                                color: Color::Rgb(r, g, b),
                            });
                        }
                    }
                }
                
                // Surface rendering - when enabled, surface replaces atom visualization
                if app.show_surface {
                    let projected_surface: Vec<renderer::ProjectedSurfacePoint> = 
                        renderer::project_surface(protein, &app.camera, width, height);
                    
                    // Sort by depth
                    let mut surface_with_depth: Vec<&renderer::ProjectedSurfacePoint> = 
                        projected_surface.iter().collect();
                    surface_with_depth.sort_by(|a, b| a.depth.partial_cmp(&b.depth).unwrap());
                    
                    for proj_surface in surface_with_depth {
                        if proj_surface.x >= 0.0 && proj_surface.x <= width && 
                           proj_surface.y >= 0.0 && proj_surface.y <= height {
                            let point = vec![(proj_surface.x as f64, proj_surface.y as f64)];
                            let (r, g, b) = surface::hydrophobicity_to_color(proj_surface.hydrophobicity);
                            
                            ctx.draw(&Points {
                                coords: &point,
                                color: Color::Rgb(r, g, b),
                            });
                        }
                    }
                }
                
                // Pass 1: draw non-selected atoms (skip if surface hydrophobicity is shown)
                if !app.show_surface {
                    for (original_idx, proj_atom) in projected_with_idx.iter() {
                        if proj_atom.x >= 0.0 && proj_atom.x <= width && 
                        proj_atom.y >= 0.0 && proj_atom.y <= height {
                            let is_selected_atom = selected_atom_idx == Some(*original_idx);
                            let is_highlighted = app.highlighted_atom_indices.contains(original_idx);
                            
                            // Skip selected atom (drawn in Pass 3)
                            // Skip highlighted atoms only if highlighting is active (drawn in Pass 2)
                            let skip_for_pass2 = residue_highlight_active && is_highlighted;
                            
                            if !is_selected_atom && !skip_for_pass2 {
                                let point = vec![(proj_atom.x as f64, proj_atom.y as f64)];
                                
                                let color = if residue_highlight_active {
                                    // Dark grey for non-residue atoms when highlighting
                                    Color::Rgb(75, 75, 75)
                                } else if app.use_bfactor_colors {
                                    // B-factor coloring
                                    let b_factor = protein.atoms[*original_idx].b_factor;
                                    let (r, g, b) = bfactor_to_color(b_factor, b_min, b_max);
                                    Color::Rgb(r, g, b)
                                } else {
                                    // Secondary structure coloring
                                    let residue_seq = protein.atoms[*original_idx].residue_seq;
                                    if let Some(&ss) = residue_colors.get(&residue_seq) {
                                        let (r, g, b) = ss.color_rgb();
                                        Color::Rgb(r, g, b)
                                    } else {
                                        Color::White
                                    }
                                };
                                
                                ctx.draw(&Points {
                                    coords: &point,
                                    color,
                                });
                            }
                        }
                    }
                
                    // Pass 2: draw atoms in the highlighted set (but not the selected atom itself)
                    if residue_highlight_active {
                        for (original_idx, proj_atom) in projected_with_idx.iter() {
                            if proj_atom.x >= 0.0 && proj_atom.x <= width && 
                               proj_atom.y >= 0.0 && proj_atom.y <= height {
                                let is_selected_atom = selected_atom_idx == Some(*original_idx);
                                let is_highlighted = app.highlighted_atom_indices.contains(original_idx);
                                
                                if !is_selected_atom && is_highlighted {
                                    let point = vec![(proj_atom.x as f64, proj_atom.y as f64)];
                                    
                                    // White for residue atoms
                                    ctx.draw(&Points {
                                        coords: &point,
                                        color: Color::White,
                                    });
                                }
                            }
                        }
                    }
                    
                    // Pass 3: draw selected atom on top in cyan
                    if let Some(selected_idx) = selected_atom_idx {
                        let proj_atom = &all_projected[selected_idx];
                        if proj_atom.x >= 0.0 && proj_atom.x <= width && 
                           proj_atom.y >= 0.0 && proj_atom.y <= height {
                            ctx.draw(&Points {
                                coords: &[
                                    (proj_atom.x as f64, proj_atom.y as f64),       // center
                                    (proj_atom.x as f64 - 1.0, proj_atom.y as f64), // left
                                    (proj_atom.x as f64 + 1.0, proj_atom.y as f64), // right
                                    (proj_atom.x as f64, proj_atom.y as f64 - 1.0), // down
                                    (proj_atom.x as f64, proj_atom.y as f64 + 1.0), // up
                                ],
                                color: Color::Cyan,
                            });
                        }
                    }
                }
            });

        frame.render_widget(canvas, main_area);
        
        // Render info overlay
        let display_mode_text = match app.display_mode {
            DisplayMode::AllAtoms => "All Atoms",
            DisplayMode::Backbone => "C-Alpha",
        };
        let connections_text = if app.show_connections { "Connected" } else { "Dots" };
        let color_mode_text = if app.use_bfactor_colors { "B-factor" } else { "Sec. Struct." };
        let mode_text = if app.inspect_mode { " | [INSPECT]" } else { "" };
        let highlight_text = if app.residue_highlight && app.selected_atom_idx.is_some() {
            " | [RESIDUE]"
        } else {
            ""
        };
        let surface_text = if app.show_surface {
            format!(" | [SURFACE: {}pts]", app.surface_point_density)
        } else {
            String::new()
        };
        
        let hbond_text = if app.show_hydrogen_bonds {
            let visible_count = protein.hbonds.iter()
                .filter(|hb| hb.energy < app.hbond_energy_threshold)
                .count();
            format!(" | [H-BONDS: {} @ {:.1} kcal/mol]", visible_count, app.hbond_energy_threshold)
        } else {
            String::new()
        };
        
        let network_text = if app.show_hbond_network {
            if let Some(ref analysis) = app.network_analysis {
                let component_count = analysis.connected_components.len();
                let largest_component = analysis.connected_components
                    .iter()
                    .map(|c| c.len())
                    .max()
                    .unwrap_or(0);
                format!(" | [NETWORK: {} components, largest = {}]", component_count, largest_component)
            } else {
                String::from(" | [NETWORK]")
            }
        } else {
            String::new()
        };
        
        let atom_count = match app.display_mode {
            DisplayMode::AllAtoms => protein.atoms.len(),
            DisplayMode::Backbone => get_calpha_indices(protein).len(),
        };
        
        let info_text = format!(
            " {} | {} atoms | [{}] [{}] [{}] | Zoom: {:.1}x{}{}{}{}{} ",
            protein.title,
            atom_count,
            display_mode_text,
            connections_text,
            color_mode_text,
            app.camera.zoom,
            mode_text,
            highlight_text,
            surface_text,
            hbond_text,
            network_text
        );
        
        let info_area = Rect {
            x: 0,
            y: 0,
            width: main_area.width,
            height: 1,
        };
        
        frame.render_widget(
            ratatui::widgets::Paragraph::new(info_text)
                .alignment(Alignment::Center)
                .style(Style::default().fg(Color::White)),
            info_area,
        );
        
        // Render atom info panel if an atom is selected
        if let Some(atom_idx) = app.selected_atom_idx {
            let info_chunks = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([
                    Constraint::Percentage(75),
                    Constraint::Percentage(25),
                ])
                .split(area);
            
            let panel_area = info_chunks[1];
            let atom = &protein.atoms[atom_idx];
            
            let ss_name = match atom.secondary_structure {
                SecondaryStructure::Helix => "α-Helix",
                SecondaryStructure::Sheet => "β-Sheet",
                SecondaryStructure::Turn => "Turn",
                SecondaryStructure::Coil => "Coil",
            };
            
            let (r, g, b) = atom.secondary_structure.color_rgb();
            let ss_color = Color::Rgb(r, g, b);
            
            let mut info_lines = vec![
                Line::from(vec![
                    Span::styled("┌─── Atom Info ───┐", Style::default().fg(Color::Yellow).bold()),
                ]),
                Line::from(""),
                Line::from(vec![
                    Span::styled("  Atom:    ", Style::default().fg(Color::Gray)),
                    Span::styled(atom.name.to_string(), Style::default().fg(Color::Yellow).bold()),
                ]),
                Line::from(vec![
                    Span::styled("  Serial:  ", Style::default().fg(Color::Gray)),
                    Span::styled(format!("{}", atom.serial), Style::default().fg(Color::White)),
                ]),
                Line::from(""),
                Line::from(vec![
                    Span::styled("  Residue: ", Style::default().fg(Color::Gray)),
                    Span::styled(atom.residue_name.to_string(), Style::default().fg(Color::White).bold()),
                ]),
                Line::from(vec![
                    Span::styled("  Number:  ", Style::default().fg(Color::Gray)),
                    Span::styled(format!("{}", atom.residue_seq), Style::default().fg(Color::White)),
                ]),
                Line::from(vec![
                    Span::styled("  Chain:   ", Style::default().fg(Color::Gray)),
                    Span::styled(atom.chain_id.to_string(), Style::default().fg(Color::Cyan).bold()),
                ]),
                Line::from(""),
                Line::from(vec![
                    Span::styled("  Structure:", Style::default().fg(Color::Gray)),
                ]),
                Line::from(vec![
                    Span::styled(format!("  {}", ss_name), Style::default().fg(ss_color).bold()),
                ]),
                Line::from(""),
                Line::from(vec![
                    Span::styled("  Position:", Style::default().fg(Color::Gray)),
                ]),
                Line::from(vec![
                    Span::styled(format!("    x: {:.2}", atom.position.x), Style::default().fg(Color::White)),
                ]),
                Line::from(vec![
                    Span::styled(format!("    y: {:.2}", atom.position.y), Style::default().fg(Color::White)),
                ]),
                Line::from(vec![
                    Span::styled(format!("    z: {:.2}", atom.position.z), Style::default().fg(Color::White)),
                ]),
                Line::from(""),
                Line::from(vec![
                    Span::styled("  B-factor: ", Style::default().fg(Color::Gray)),
                    Span::styled(format!("{:.2}", atom.b_factor), Style::default().fg(Color::White)),
                ]),
                Line::from(""),
                Line::from(vec![
                    Span::styled(
                        if app.residue_highlight {
                            "  Press 'r' to hide residue"
                        } else {
                            "  Press 'r' to highlight residue"
                        },
                        Style::default().fg(Color::DarkGray).italic()
                    ),
                ]),
                Line::from(""),
            ];
            
            // Candidates section if there are multiple candidates
            if !app.candidate_atoms.is_empty() {
                info_lines.push(Line::from(vec![
                    Span::styled("├─ Nearby Atoms ──┤", Style::default().fg(Color::Cyan)),
                ]));
                info_lines.push(Line::from(""));
                info_lines.push(Line::from(vec![
                    Span::styled("  Use ↑↓ to cycle", Style::default().fg(Color::DarkGray).italic()),
                ]));
                info_lines.push(Line::from(""));
                
                for (i, &(candidate_idx, distance)) in app.candidate_atoms.iter().enumerate() {
                    let candidate = &protein.atoms[candidate_idx];
                    let is_current = i == app.candidate_selection_idx;
                    let prefix = if is_current { "→ " } else { "  " };
                    let color = if is_current { Color::Yellow } else { Color::White };
                    
                    info_lines.push(Line::from(vec![
                        Span::styled(prefix, Style::default().fg(color)),
                        Span::styled(candidate.name.to_string(), Style::default().fg(color).bold()),
                        Span::styled(format!(" ({})", candidate.residue_name), Style::default().fg(Color::Gray)),
                    ]));
                    info_lines.push(Line::from(vec![
                        Span::styled(format!("    d={:.1}px", distance), Style::default().fg(Color::DarkGray)),
                    ]));
                }
                info_lines.push(Line::from(""));
            }
            
            info_lines.push(Line::from(vec![
                Span::styled("└─────────────────┘", Style::default().fg(Color::Yellow)),
            ]));
            
            let info_paragraph = ratatui::widgets::Paragraph::new(info_lines)
                .style(Style::default().bg(Color::Black))
                .block(
                    ratatui::widgets::Block::default()
                        .borders(ratatui::widgets::Borders::ALL)
                        .border_style(Style::default().fg(Color::Yellow))
                );
            
            frame.render_widget(info_paragraph, panel_area);
        }
    } else {
        // No protein loaded - show help text
        let help_text = [
            "PixelFold - 3D Protein Viewer",
            "",
            "Usage: pixelfold <protein.pdb>",
            "",
            "Supports PDB and mmCIF formats",
            "",
            "Press 'q' to quit",
        ];
        
        let paragraph = ratatui::widgets::Paragraph::new(help_text.join("\n"))
            .alignment(Alignment::Center)
            .style(Style::default().fg(Color::Gray));
        
        frame.render_widget(paragraph, area);
    }
}
