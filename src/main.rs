use anyhow::Result;
use crossterm::event::{self, KeyCode, KeyModifiers};
use pixelfold::{
    parser,
    renderer::{self, Camera},
    Protein,
    SecondaryStructure,
};
use ratatui::prelude::*;
use ratatui::widgets::canvas::{Canvas, Points};
use std::collections::HashSet;
use std::env;


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
        }
    }

    fn load_protein(&mut self, path: &str, width: f32, height: f32) -> Result<()> {
        self.protein = Some(parser::load_protein(path)?);
        
        // Auto-frame the protein when loaded
        if let Some(ref protein) = self.protein {
            renderer::auto_frame_protein(protein, &mut self.camera, width, height);
        }
        
        Ok(())
    }
}

fn main() -> Result<()> {
    let mut app = App::new();
    let mut terminal = ratatui::init();
    
    crossterm::execute!(
        std::io::stdout(),
        crossterm::event::EnableMouseCapture
    )?;
    
    let args: Vec<String> = env::args().collect();
    if args.len() > 1 {
        // Initial terminal size
        let size = terminal.size()?;
        let width = size.width as f32 * 2.0;  // Braille canvas width
        let height = size.height as f32 * 4.0; // Braille canvas height
        
        match app.load_protein(&args[1], width, height) {
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
        
        // Pan controls (Arrow keys) - or cycle through candidates in inspect mode
        KeyCode::Up => {
            if app.inspect_mode && !app.candidate_atoms.is_empty() {
                if app.candidate_selection_idx == 0 {
                    app.candidate_selection_idx = app.candidate_atoms.len() - 1;
                } else {
                    app.candidate_selection_idx -= 1;
                }
                app.selected_atom_idx = Some(app.candidate_atoms[app.candidate_selection_idx].0);
                update_highlighted_atoms(app, width, height);
            } else {
                app.camera.pan_camera(0.0, pan_speed);
            }
        }
        KeyCode::Down => {
            if app.inspect_mode && !app.candidate_atoms.is_empty() {
                app.candidate_selection_idx = (app.candidate_selection_idx + 1) % app.candidate_atoms.len();
                app.selected_atom_idx = Some(app.candidate_atoms[app.candidate_selection_idx].0);
                update_highlighted_atoms(app, width, height);
            } else {
                app.camera.pan_camera(0.0, -pan_speed);
            }
        }
        KeyCode::Left => app.camera.pan_camera(-pan_speed, 0.0),
        KeyCode::Right => app.camera.pan_camera(pan_speed, 0.0),

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
                
                // Store canvas dimensions for click detection (done in closure, so need mutable access)
                // We'll store it after the paint closure
                
                let mut residue_colors: std::collections::HashMap<u32, SecondaryStructure> = 
                    std::collections::HashMap::new();
                
                // Assign colors based on original atom order
                for atom in &protein.atoms {
                    residue_colors.entry(atom.residue_seq)
                        .or_insert(atom.secondary_structure);
                }
                
                // Project and sort for rendering
                let mut projected_with_idx: Vec<(usize, renderer::ProjectedAtom)> = 
                    renderer::project_protein(protein, &app.camera, width, height)
                        .into_iter()
                        .enumerate()
                        .collect();
                
                // Sort by depth (painter's algorithm)
                projected_with_idx.sort_by(|a, b| a.1.depth.partial_cmp(&b.1.depth).unwrap());
                
                let selected_atom_idx = app.selected_atom_idx;
                let residue_highlight_active = app.residue_highlight && selected_residue_seq.is_some();
                
                // Pass 1: draw non-selected atoms
                for (original_idx, proj_atom) in projected_with_idx.iter() {
                    if proj_atom.x >= 0.0 && proj_atom.x <= width && 
                       proj_atom.y >= 0.0 && proj_atom.y <= height {
                        let is_selected_atom = selected_atom_idx == Some(*original_idx);
                        let is_highlighted = app.highlighted_atom_indices.contains(original_idx);
                        let residue_seq = protein.atoms[*original_idx].residue_seq;
                        
                        // Skip selected atom (drawn in Pass 3)
                        // Skip highlighted atoms only if highlighting is active (drawn in Pass 2)
                        let skip_for_pass2 = residue_highlight_active && is_highlighted;
                        
                        if !is_selected_atom && !skip_for_pass2 {
                            let point = vec![(proj_atom.x as f64, proj_atom.y as f64)];
                            
                            if residue_highlight_active {
                                // Dark grey for non-residue atoms when highlighting
                                ctx.draw(&Points {
                                    coords: &point,
                                    color: Color::Rgb(75, 75, 75),
                                });
                            } else {
                                if let Some(&ss) = residue_colors.get(&residue_seq) {
                                    let (r, g, b) = ss.color_rgb();
                                    ctx.draw(&Points {
                                        coords: &point,
                                        color: Color::Rgb(r, g, b),
                                    });
                                }
                            }
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
                    if let Some((_, proj_atom)) = projected_with_idx.iter()
                        .find(|(idx, _)| *idx == selected_idx) {
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
        let mode_text = if app.inspect_mode { " | [INSPECT MODE]" } else { "" };
        let highlight_text = if app.residue_highlight && app.selected_atom_idx.is_some() {
            " | [RESIDUE HIGHLIGHTED]"
        } else {
            ""
        };
        let info_text = format!(
            " {} | Atoms: {} | Zoom: {:.1}x{}{} | Controls: WASDZX = rotate, +/- = zoom, arrows = pan, f = frame, i = inspect, q = quit ",
            protein.title,
            protein.atoms.len(),
            app.camera.zoom,
            mode_text,
            highlight_text
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
                    Span::styled(format!("{}", atom.name), Style::default().fg(Color::Yellow).bold()),
                ]),
                Line::from(vec![
                    Span::styled("  Serial:  ", Style::default().fg(Color::Gray)),
                    Span::styled(format!("{}", atom.serial), Style::default().fg(Color::White)),
                ]),
                Line::from(""),
                Line::from(vec![
                    Span::styled("  Residue: ", Style::default().fg(Color::Gray)),
                    Span::styled(format!("{}", atom.residue_name), Style::default().fg(Color::White).bold()),
                ]),
                Line::from(vec![
                    Span::styled("  Number:  ", Style::default().fg(Color::Gray)),
                    Span::styled(format!("{}", atom.residue_seq), Style::default().fg(Color::White)),
                ]),
                Line::from(vec![
                    Span::styled("  Chain:   ", Style::default().fg(Color::Gray)),
                    Span::styled(format!("{}", atom.chain_id), Style::default().fg(Color::Cyan).bold()),
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
                        Span::styled(format!("{}", candidate.name), Style::default().fg(color).bold()),
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
        let help_text = vec![
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
