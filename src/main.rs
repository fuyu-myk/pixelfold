use anyhow::Result;
use crossterm::event::{self, KeyCode, KeyModifiers};
use pixelfold::visualization::{network, renderer, surface};
use pixelfold::{
    App, DisplayMode, SecondaryStructure,
};
use ratatui::prelude::*;
use ratatui::widgets::canvas::{Canvas, Points};
use std::env;


fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    
    // Parse command-line arguments
    let mut protein_name: Option<String> = None;
    let mut skip_surface = false;
    
    for arg in args.iter().skip(1) {
        if arg == "--no-surface" {
            skip_surface = true;
        } else if arg == "--help" || arg == "-h" {
            println!("PixelFold - Terminal-based 3D protein structure viewer");
            println!();
            println!("--- Searching and fetching structures ---");
            println!("Usage: pixelfold [OPTIONS] <optional_protein_name>");
            println!();
            println!("Options (compulsory):");
            println!("  --fetch, -f     Opens a search interface to fetch protein structures by name");
            println!();
            println!("--- Visualization ---");
            println!("Usage: pixelfold <path/to/protein.pdb|protein.cif> [OPTIONS]");
            println!("       pixelfold <protein> [OPTIONS]");
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
            protein_name = Some(arg.clone());
        }
    }
    
    let mut app = App::new();
    let mut terminal = ratatui::init();
    
    crossterm::execute!(
        std::io::stdout(),
        crossterm::event::EnableMouseCapture
    )?;
    
    let result = if args.iter().any(|arg| arg == "--fetch" || arg == "-f") {
        pixelfold::search::fetch_structures(&mut terminal, protein_name)
    } else {
        if let Some(name) = protein_name {
            // Initial terminal size
            let size = terminal.size()?;
            let width = size.width as f32 * 2.0;   // Braille canvas width
            let height = size.height as f32 * 4.0; // Braille canvas height

            let project_dir = env!("CARGO_MANIFEST_DIR");
            
            let path = if name.starts_with("data") {
                format!("{}/{}", project_dir, name)
            } else {
                if name.ends_with(".cif") {
                    format!("{}/data/{}", project_dir, name)
                } else if name.ends_with(".pdb") {
                    format!("{}/data/{}", project_dir, name)
                } else {
                    format!("{}/data/{}.cif", project_dir, name)
                }
            };
            
            match app.load_protein(&path, width, height, skip_surface) {
                Ok(_) => {},
                Err(e) => eprintln!("Failed to load protein: {}", e),
            }
        }

        run_app(&mut terminal, &mut app)
    };
    
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
        
        if app.redraw_needed {
            terminal.draw(|frame| ui(frame, app))?;
            app.redraw_needed = false;
        }

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

fn handle_input(app: &mut App, key: KeyCode, _modifiers: KeyModifiers, width: f32, height: f32) -> Result<()> {
    let rotation_speed = 0.1;
    let zoom_speed = 0.1;
    let pan_speed = 5.0;

    match key {
        KeyCode::Char('q') => {}
        
        // Rotation controls (WASD)
        KeyCode::Char('w') => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.rotate(rotation_speed, 0.0, 0.0);
        }
        KeyCode::Char('s') => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.rotate(-rotation_speed, 0.0, 0.0);
        }
        KeyCode::Char('a') => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.rotate(0.0, rotation_speed, 0.0);
        }
        KeyCode::Char('d') => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.rotate(0.0, -rotation_speed, 0.0);
        }
        KeyCode::Char('z') => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.rotate(0.0, 0.0, rotation_speed);
        }
        KeyCode::Char('x') => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.rotate(0.0, 0.0, -rotation_speed);
        }
        
        // Zoom controls
        KeyCode::Char('+') | KeyCode::Char('=') => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.adjust_zoom(zoom_speed);
        }
        KeyCode::Char('-') | KeyCode::Char('_') => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.adjust_zoom(-zoom_speed);
        }
        
        // Pan controls (Arrow keys)
        // Cycle through candidates in inspect mode (up and down)
        // Adjust surface density (up and down)
        // Adjust H-bond energy threshold (up and down when H-bonds visible)
        KeyCode::Up => {
            if app.show_hydrogen_bonds && !app.show_surface && !app.inspect_mode {
                // More negative = stronger bonds only
                app.hbond_energy_threshold = (app.hbond_energy_threshold - 0.5).max(-10.0);
                app.redraw_needed = true;
                
                // Recompute network analysis if network mode is active
                if app.show_hbond_network {
                    if let Some(ref graph) = app.hbond_graph {
                        let filtered = graph.filter_by_energy(app.hbond_energy_threshold);
                        app.network_analysis = Some(filtered.analyze());
                    }
                }
            } else if app.show_surface {
                app.surface_point_density = (app.surface_point_density + 25).min(500);
                app.redraw_needed = true;
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
                app.redraw_needed = true;
                pixelfold::update_highlighted_atoms(app, width, height);
            } else {
                app.redraw_needed = true;
                app.projected_atom_cache = None;
                app.camera.pan_camera(0.0, -pan_speed);
            }
        }
        KeyCode::Down => {
            if app.show_hydrogen_bonds && !app.show_surface && !app.inspect_mode {
                // Less negative = more bonds shown
                app.hbond_energy_threshold = (app.hbond_energy_threshold + 0.5).min(-0.1);
                app.redraw_needed = true;
                
                // Recompute network analysis if network mode is active
                if app.show_hbond_network {
                    if let Some(ref graph) = app.hbond_graph {
                        let filtered = graph.filter_by_energy(app.hbond_energy_threshold);
                        app.network_analysis = Some(filtered.analyze());
                    }
                }
            } else if app.show_surface {
                app.surface_point_density = (app.surface_point_density.saturating_sub(25)).max(100);
                app.redraw_needed = true;
                if let Some(ref mut protein) = app.protein {
                    let surface_calculator = surface::SurfaceCalculator::new(1.4, app.surface_point_density);
                    protein.surface_points = surface_calculator.calculate_surface(&protein.atoms);
                }
            } else if app.inspect_mode && !app.candidate_atoms.is_empty() {
                app.candidate_selection_idx = (app.candidate_selection_idx + 1) % app.candidate_atoms.len();
                app.selected_atom_idx = Some(app.candidate_atoms[app.candidate_selection_idx].0);
                app.redraw_needed = true;
                pixelfold::update_highlighted_atoms(app, width, height);
            } else {
                app.redraw_needed = true;
                app.projected_atom_cache = None;
                app.camera.pan_camera(0.0, pan_speed);
            }
        }
        KeyCode::Left => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.pan_camera(pan_speed, 0.0);
        }
        KeyCode::Right => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.pan_camera(-pan_speed, 0.0);
        }

        // Inspect mode
        KeyCode::Char('i') => {
            app.inspect_mode = !app.inspect_mode;
            app.redraw_needed = true;

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
                app.redraw_needed = true;
            }
        }
        
        // Reset view
        KeyCode::Char('f') => {
            if let Some(ref protein) = app.protein {
                renderer::auto_frame_protein(protein, &mut app.camera, width, height);
                app.redraw_needed = true;
                app.projected_atom_cache = None;
            }
        }
        
        // Display mode controls
        KeyCode::Char('1') => {
            app.display_mode = DisplayMode::AllAtoms;
            app.redraw_needed = true;
        }
        KeyCode::Char('2') => {
            app.display_mode = DisplayMode::Backbone;
            app.redraw_needed = true;
        }
        
        // Toggle backbone connections
        KeyCode::Char('c') => {
            app.show_connections = !app.show_connections;
            app.redraw_needed = true;
        }
        
        // Toggle B-factor coloring
        KeyCode::Char('b') => {
            app.use_bfactor_colors = !app.use_bfactor_colors;
            app.redraw_needed = true;
        }
        
        // Toggle surface visualization
        KeyCode::Char('v') => {
            app.show_surface = !app.show_surface;
            app.redraw_needed = true;
            
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
            app.redraw_needed = true;

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
            app.redraw_needed = true;
            
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
            app.redraw_needed = true;
            app.projected_atom_cache = None;

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
                let candidates = pixelfold::pick_atoms_along_ray(
                    protein,
                    &mut app.camera,
                    click_x,
                    click_y,
                    canvas_width,
                    canvas_height,
                );
                
                if !candidates.is_empty() {
                    app.candidate_atoms = candidates.into_iter().take(5).collect(); // Top 5 candidates
                    app.candidate_selection_idx = 0;
                    app.selected_atom_idx = Some(app.candidate_atoms[0].0);
                    pixelfold::update_highlighted_atoms(app, canvas_width, canvas_height);
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

fn ui(frame: &mut Frame, app: &mut App) {
    let area = frame.area();
    
    if let Some(ref protein) = app.protein {
        // Info bar (1 line)
        let vertical_chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(1),  // Top info bar
                Constraint::Min(0),     // Content area (canvas + optional inspect panel)
            ])
            .split(area);
        
        let info_bar_area = vertical_chunks[0];
        let content_area = vertical_chunks[1];
        
        // Canvas and optional inspect panel (left: main view, right: info panel)
        let main_area = if app.selected_atom_idx.is_some() {
            let chunks = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([
                    Constraint::Percentage(75),
                    Constraint::Percentage(25),
                ])
                .split(content_area);
            chunks[0]
        } else {
            content_area
        };
        
        // Store canvas dimensions for click detection
        app.last_canvas_width = main_area.width as f32 * 2.0;
        app.last_canvas_height = main_area.height as f32 * 4.0;

        if app.camera.cached_view_matrix.is_none() {
            app.camera.get_view_matrix();
        }
        
        let width = main_area.width as f32 * 2.0;
        let height = main_area.height as f32 * 4.0;
        
        if app.projected_atom_cache.is_none() {
            app.projected_atom_cache = Some(renderer::project_protein(protein, &app.camera, width, height));
        }
        
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
                    DisplayMode::Backbone => app.ca_indices.clone(),
                };
                
                // Calculate B-factor range if using B-factor coloring
                let (b_min, b_max) = if app.use_bfactor_colors {
                    pixelfold::calculate_bfactor_range(protein)
                } else {
                    (0.0, 100.0)
                };
                
                let all_projected: Vec<renderer::ProjectedAtom> = match &app.projected_atom_cache {
                    Some(cache) => cache.clone(),
                    None => renderer::project_protein(protein, &app.camera, width, height),
                };
                
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
                    let connections = &app.backbone_connections;
                    
                    for &(idx1, idx2) in connections {
                        let proj1 = &all_projected[idx1];
                        let proj2 = &all_projected[idx2];
                        
                        // Only draw if both atoms are on screen
                        if proj1.x >= 0.0 && proj1.x <= width && 
                           proj1.y >= 0.0 && proj1.y <= height &&
                           proj2.x >= 0.0 && proj2.x <= width && 
                           proj2.y >= 0.0 && proj2.y <= height {
                            
                            let line_points = pixelfold::draw_line(proj1.x, proj1.y, proj2.x, proj2.y);
                            
                            // Color the line based on coloring mode (av of two atoms)
                            let color = if app.use_bfactor_colors {
                                let b1 = protein.atoms[idx1].b_factor;
                                let b2 = protein.atoms[idx2].b_factor;
                                let avg_b = (b1 + b2) / 2.0;
                                let (r, g, b) = pixelfold::bfactor_to_color(avg_b, b_min, b_max);
                                Color::Rgb(r, g, b)
                            } else {
                                // Use secondary structure color (slightly dimmed for lines)
                                let residue_seq = protein.atoms[idx1].residue_seq;
                                if let Some(&ss) = app.residue_colors.get(&residue_seq) {
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
                            
                            let line_points = pixelfold::draw_dashed_line(
                                proj_donor.x, 
                                proj_donor.y, 
                                proj_acceptor.x, 
                                proj_acceptor.y,
                                3  // Dash spacing
                            );
                            
                            let (r, g, b) = pixelfold::hbond_energy_to_color(hbond.energy);
                            
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
                                    let (r, g, b) = pixelfold::bfactor_to_color(b_factor, b_min, b_max);
                                    Color::Rgb(r, g, b)
                                } else {
                                    // Secondary structure coloring
                                    let residue_seq = protein.atoms[*original_idx].residue_seq;
                                    if let Some(&ss) = app.residue_colors.get(&residue_seq) {
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
            DisplayMode::Backbone => pixelfold::get_calpha_indices(protein).len(),
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
        
        // Render protein info at the top bar
        frame.render_widget(
            ratatui::widgets::Paragraph::new(info_text)
                .alignment(Alignment::Center)
                .style(Style::default().fg(Color::White)),
            info_bar_area,
        );
        
        // Render atom info panel if an atom is selected
        if let Some(atom_idx) = app.selected_atom_idx {
            let info_chunks = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([
                    Constraint::Percentage(75),
                    Constraint::Percentage(25),
                ])
                .split(content_area);
            
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
