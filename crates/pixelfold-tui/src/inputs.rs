use anyhow::Result;
use crossterm::event::{self, KeyCode, KeyModifiers};

use pixelfold_core::DisplayMode;
use pixelfold_render::renderer;

use crate::App;

const SLAB_DEFAULT: (f32, f32) = (0.3, 0.7);
const SLAB_STEP: f32 = 0.05;
/// Smallest half-width of the slab, so narrowing never inverts the band.
const SLAB_MIN_HALF: f32 = 0.025;

pub(crate) fn handle_input(app: &mut App, key: KeyCode, _modifiers: KeyModifiers) -> Result<()> {
    let rotation_speed = 0.1;
    let zoom_speed = 0.1;
    let pan_speed = 5.0;

    match key {
        KeyCode::Char('q') => {}

        // Rotation (WASD + z/x for roll)
        KeyCode::Char('w') => rotate(app, rotation_speed, 0.0, 0.0),
        KeyCode::Char('s') => rotate(app, -rotation_speed, 0.0, 0.0),
        KeyCode::Char('a') => rotate(app, 0.0, rotation_speed, 0.0),
        KeyCode::Char('d') => rotate(app, 0.0, -rotation_speed, 0.0),
        KeyCode::Char('z') => rotate(app, 0.0, 0.0, rotation_speed),
        KeyCode::Char('x') => rotate(app, 0.0, 0.0, -rotation_speed),

        // Zoom
        KeyCode::Char('+') | KeyCode::Char('=') => {
            app.camera.adjust_zoom(zoom_speed);
            app.redraw_needed = true;
        }
        KeyCode::Char('-') | KeyCode::Char('_') => {
            app.camera.adjust_zoom(-zoom_speed);
            app.redraw_needed = true;
        }

        // Up/Down adjust the hydrogen-bond energy threshold when bonds are shown,
        // otherwise they pan vertically.
        KeyCode::Up => {
            if app.show_hydrogen_bonds && !app.inspect_mode {
                app.hbond_energy_threshold = (app.hbond_energy_threshold - 0.5).max(-10.0);
            } else {
                app.camera.pan_camera(0.0, pan_speed);
            }
            app.redraw_needed = true;
        }
        KeyCode::Down => {
            if app.show_hydrogen_bonds && !app.inspect_mode {
                app.hbond_energy_threshold = (app.hbond_energy_threshold + 0.5).min(-0.1);
            } else {
                app.camera.pan_camera(0.0, -pan_speed);
            }
            app.redraw_needed = true;
        }
        KeyCode::Left => {
            app.camera.pan_camera(-pan_speed, 0.0);
            app.redraw_needed = true;
        }
        KeyCode::Right => {
            app.camera.pan_camera(pan_speed, 0.0);
            app.redraw_needed = true;
        }

        // Inspect mode (click to pick atoms)
        KeyCode::Char('i') => {
            app.inspect_mode = !app.inspect_mode;
            app.redraw_needed = true;
            if !app.inspect_mode {
                app.selected_atom_idx = None;
                app.residue_highlight = false;
                app.highlighted_atom_indices.clear();
            }
        }
        KeyCode::Char('r') => {
            if app.inspect_mode && app.selected_atom_idx.is_some() {
                app.residue_highlight = !app.residue_highlight;
                app.redraw_needed = true;
            }
        }

        // Reset the view to fit the structure.
        KeyCode::Char('f') => {
            let dims = app
                .last_frame
                .as_ref()
                .map(|frame| (frame.fb.width() as f32, frame.fb.height() as f32));
            if let (Some((width, height)), Some(protein)) = (dims, app.protein.as_ref()) {
                renderer::auto_frame_protein(protein, &mut app.camera, width, height);
                app.redraw_needed = true;
            }
        }

        // Display mode
        KeyCode::Char('1') => {
            app.display_mode = DisplayMode::AllAtoms;
            app.redraw_needed = true;
        }
        KeyCode::Char('2') => {
            app.display_mode = DisplayMode::Backbone;
            app.redraw_needed = true;
        }

        KeyCode::Char('c') => {
            app.show_connections = !app.show_connections;
            app.redraw_needed = true;
        }
        KeyCode::Char('b') => {
            app.use_bfactor_colors = !app.use_bfactor_colors;
            app.redraw_needed = true;
        }
        KeyCode::Char('h') => {
            app.show_hydrogen_bonds = !app.show_hydrogen_bonds;
            app.redraw_needed = true;
        }

        // Toggle the linked network pane. On first open the network is built on a
        // worker thread, so a large structure's detection, layout, and accessibility
        // pass never freeze the event loop; the pane shows a placeholder until the
        // result arrives (see `run_app`).
        KeyCode::Char('g') => {
            app.show_network = !app.show_network;
            if app.show_network
                && app.network_view.is_none()
                && app.network_build.is_none()
                && let Some(protein) = app.protein.clone()
            {
                let (tx, rx) = std::sync::mpsc::channel();
                std::thread::spawn(move || {
                    let view = crate::network::NetworkView::build(&protein);
                    tx.send(view).ok();
                });

                app.network_build = Some(rx);
            }
            app.redraw_needed = true;
        }

        // Switch the network layout between force-directed and spatial.
        KeyCode::Char('l') => {
            if app.show_network
                && let Some(view) = app.network_view.as_mut()
            {
                view.toggle_mode();
                app.redraw_needed = true;
            }
        }

        // Switch the network's node colouring between secondary structure and burial.
        KeyCode::Char('e') => {
            if app.show_network
                && let Some(view) = app.network_view.as_mut()
            {
                view.toggle_color();
                app.redraw_needed = true;
            }
        }

        // Cycle the network's interaction-kind filter (all kinds, then each in turn).
        KeyCode::Char('t') => {
            if app.show_network
                && let Some(view) = app.network_view.as_mut()
            {
                view.cycle_kind_filter();
                app.redraw_needed = true;
            }
        }

        // Cycle the network's chain filter (all, intra-chain, inter-chain interface).
        KeyCode::Char('y') => {
            if app.show_network
                && let Some(view) = app.network_view.as_mut()
            {
                view.cycle_chain_filter();
                app.redraw_needed = true;
            }
        }

        // Slab (clipping planes): k toggles, [ ] move the band in depth, , . resize it.
        KeyCode::Char('k') => {
            app.slab = match app.slab {
                None => Some(SLAB_DEFAULT),
                Some(_) => None,
            };
            app.redraw_needed = true;
        }
        KeyCode::Char('[') => adjust_slab(app, |s| shift_slab(s, -SLAB_STEP)),
        KeyCode::Char(']') => adjust_slab(app, |s| shift_slab(s, SLAB_STEP)),
        KeyCode::Char(',') => adjust_slab(app, |s| resize_slab(s, -SLAB_STEP)),
        KeyCode::Char('.') => adjust_slab(app, |s| resize_slab(s, SLAB_STEP)),

        _ => {}
    }

    Ok(())
}

fn rotate(app: &mut App, pitch: f32, yaw: f32, roll: f32) {
    app.camera.rotate(pitch, yaw, roll);
    app.redraw_needed = true;
}

fn adjust_slab(app: &mut App, f: impl Fn((f32, f32)) -> (f32, f32)) {
    if let Some(slab) = app.slab {
        app.slab = Some(f(slab));
        app.redraw_needed = true;
    }
}

/// Move the band nearer or farther, keeping its width and staying in `0..1`.
fn shift_slab((far, near): (f32, f32), delta: f32) -> (f32, f32) {
    let width = near - far;
    let new_far = (far + delta).clamp(0.0, 1.0 - width);
    (new_far, new_far + width)
}

/// Widen (delta > 0) or narrow (delta < 0) the band about its centre.
fn resize_slab((far, near): (f32, f32), delta: f32) -> (f32, f32) {
    let center = (far + near) / 2.0;
    let half = ((near - far) / 2.0 + delta).clamp(SLAB_MIN_HALF, 0.5);
    ((center - half).max(0.0), (center + half).min(1.0))
}

pub(crate) fn handle_mouse(app: &mut App, mouse: event::MouseEvent) -> Result<()> {
    use event::MouseEventKind;

    if !matches!(mouse.kind, MouseEventKind::Down(event::MouseButton::Left)) {
        return Ok(());
    }

    // A click inside the network pane selects a residue there, which the
    // structure view then highlights.
    if let Some(area) = app.network_area
        && within(area, mouse.column, mouse.row)
    {
        if let Some(node) = app
            .network_view
            .as_ref()
            .and_then(|view| view.pick(mouse.column, mouse.row))
        {
            let atom = app.network_view.as_ref().unwrap().atom_of_node(node);
            app.selected_atom_idx = Some(atom);
            crate::update_highlighted_atoms(app);
            app.redraw_needed = true;
        }

        return Ok(());
    }

    // Otherwise a click in the 3D view picks an atom, in inspect mode only.
    if !app.inspect_mode {
        return Ok(());
    }
    match crate::pick_at(app, mouse.column, mouse.row) {
        Some(idx) => {
            app.selected_atom_idx = Some(idx);
            crate::update_highlighted_atoms(app);
        }
        None => {
            app.selected_atom_idx = None;
            app.highlighted_atom_indices.clear();
        }
    }
    app.redraw_needed = true;

    Ok(())
}

/// Whether a terminal cell falls inside `area`.
fn within(area: ratatui::layout::Rect, column: u16, row: u16) -> bool {
    column >= area.x && column < area.x + area.width && row >= area.y && row < area.y + area.height
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: (f32, f32), b: (f32, f32)) {
        assert!(
            (a.0 - b.0).abs() < 1e-4 && (a.1 - b.1).abs() < 1e-4,
            "{a:?} != {b:?}"
        );
    }

    #[test]
    fn shift_keeps_width_and_stays_in_range() {
        close(shift_slab((0.3, 0.7), 0.05), (0.35, 0.75));
        // Clamped at the near edge without inverting the band.
        close(shift_slab((0.3, 0.7), 1.0), (0.6, 1.0));
        close(shift_slab((0.3, 0.7), -1.0), (0.0, 0.4));
    }

    #[test]
    fn resize_never_inverts_the_band() {
        let (far, near) = resize_slab((0.45, 0.55), -1.0);
        assert!(
            near - far >= 2.0 * SLAB_MIN_HALF - f32::EPSILON,
            "band collapsed"
        );
        assert!(far >= 0.0 && near <= 1.0);
    }
}
