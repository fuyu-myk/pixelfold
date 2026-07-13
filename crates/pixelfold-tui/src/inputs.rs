use anyhow::Result;
use crossterm::event::{self, KeyCode, KeyModifiers};

use pixelfold_core::{DisplayMode, rin, sasa};
use pixelfold_render::renderer;

use crate::App;

pub(crate) fn handle_input(app: &mut App, key: KeyCode, _modifiers: KeyModifiers) -> Result<()> {
    let rotation_speed = 0.1;
    let zoom_speed = 0.1;
    let pan_speed = 5.0;

    // Input, rendering, and mouse picking share one coordinate system.
    let width = app.last_canvas_width;
    let height = app.last_canvas_height;

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
                    let surface_calculator =
                        sasa::SurfaceCalculator::new(1.4, app.surface_point_density);
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
                crate::update_highlighted_atoms(app, width, height);
            } else {
                app.redraw_needed = true;
                app.projected_atom_cache = None;
                app.camera.pan_camera(0.0, pan_speed);
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
                    let surface_calculator =
                        sasa::SurfaceCalculator::new(1.4, app.surface_point_density);
                    protein.surface_points = surface_calculator.calculate_surface(&protein.atoms);
                }
            } else if app.inspect_mode && !app.candidate_atoms.is_empty() {
                app.candidate_selection_idx =
                    (app.candidate_selection_idx + 1) % app.candidate_atoms.len();
                app.selected_atom_idx = Some(app.candidate_atoms[app.candidate_selection_idx].0);
                app.redraw_needed = true;
                crate::update_highlighted_atoms(app, width, height);
            } else {
                app.redraw_needed = true;
                app.projected_atom_cache = None;
                app.camera.pan_camera(0.0, -pan_speed);
            }
        }
        KeyCode::Left => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.pan_camera(-pan_speed, 0.0);
        }
        KeyCode::Right => {
            app.redraw_needed = true;
            app.projected_atom_cache = None;
            app.camera.pan_camera(pan_speed, 0.0);
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
                app.projected_atom_cache = None;
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
                        let surface_calculator =
                            sasa::SurfaceCalculator::new(1.4, app.surface_point_density);
                        protein.surface_points =
                            surface_calculator.calculate_surface(&protein.atoms);
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
                    app.hbond_graph = Some(rin::HBondGraph::build(protein));
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
                        app.hbond_graph = Some(rin::HBondGraph::build(protein));
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

pub(crate) fn handle_mouse(app: &mut App, mouse: event::MouseEvent) -> Result<()> {
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
                let candidates = crate::pick_atoms_along_ray(
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
                    crate::update_highlighted_atoms(app, canvas_width, canvas_height);
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
