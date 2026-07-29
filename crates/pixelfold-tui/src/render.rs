//! The interactive view's render path: current camera and toggles to a fresh
//! RGBA framebuffer. Atoms are z-buffered sphere impostors; the backbone and
//! hydrogen-bond overlays are depth-tested lines on top. The atom-id buffer the
//! rasteriser fills is what mouse picking reads back.

use pixelfold_core::{DisplayMode, Protein, calculate_bfactor_range};
use pixelfold_render::{
    Coloring, Framebuffer, NO_ID, RenderOptions, Scene, bfactor_to_color, draw_segment,
    hbond_energy_to_color, rasterize,
};

use crate::App;

pub(crate) const BACKGROUND: [u8; 4] = [13, 13, 18, 255];
const FOG_STRENGTH: f32 = 0.45;
const SELECTED_COLOR: [u8; 3] = [0, 255, 255];
const HIGHLIGHT_COLOR: [u8; 3] = [255, 255, 255];
/// Backbone lines are dimmed slightly below their residue colour so the atoms
/// read as the foreground.
const CONNECTION_DIM: f32 = 0.8;
/// Spacing between the dots of a dashed hydrogen-bond line.
const HBOND_DASH: usize = 3;

/// Rasterise the current view into a fresh framebuffer sized `width` x `height`.
pub(crate) fn render_scene(app: &App, width: u32, height: u32) -> Framebuffer {
    let mut fb = Framebuffer::new(width, height, BACKGROUND);
    let Some(protein) = app.protein.as_ref() else {
        return fb;
    };

    let indices: Vec<usize> = match app.display_mode {
        DisplayMode::AllAtoms => (0..protein.atoms.len()).collect(),
        DisplayMode::Backbone => app.ca_indices.clone(),
    };

    let coloring = if app.use_bfactor_colors {
        Coloring::Bfactor
    } else {
        Coloring::SecondaryStructure
    };
    let mut scene = Scene::from_protein(protein, coloring, Some(&indices));
    recolor_selection(app, &mut scene);

    let opts = RenderOptions {
        radius_scale: app.radius_scale,
        depth_cue: if app.show_depth_cue {
            FOG_STRENGTH
        } else {
            0.0
        },
        slab: app.slab,
        ..RenderOptions::default()
    };
    rasterize(&scene, &app.camera, &mut fb, &opts);

    draw_overlays(app, protein, &mut fb);
    fb
}

/// Tint the selected atom, and its whole residue when residue highlight is on.
/// Colours are overridden in place on the scene's parallel arrays, keyed by
/// source atom index.
fn recolor_selection(app: &App, scene: &mut Scene) {
    let Some(selected) = app.selected_atom_idx else {
        return;
    };
    let highlight = app.residue_highlight;

    for (slot, &atom_idx) in scene.atom_idx.iter().enumerate() {
        let atom_idx = atom_idx as usize;
        if atom_idx == selected {
            scene.colors[slot] = SELECTED_COLOR;
        } else if highlight && app.highlighted_atom_indices.contains(&atom_idx) {
            scene.colors[slot] = HIGHLIGHT_COLOR;
        }
    }
}

fn draw_overlays(app: &App, protein: &Protein, fb: &mut Framebuffer) {
    let highlight = app.residue_highlight && app.selected_atom_idx.is_some();

    if app.show_connections && !highlight && app.display_mode == DisplayMode::Backbone {
        let (b_min, b_max) = if app.use_bfactor_colors {
            calculate_bfactor_range(protein)
        } else {
            (0.0, 100.0)
        };

        for &(i, j) in &app.backbone_connections {
            let color = if app.use_bfactor_colors {
                let avg = (protein.atoms[i].b_factor + protein.atoms[j].b_factor) / 2.0;
                let (r, g, b) = bfactor_to_color(avg, b_min, b_max);
                [r, g, b]
            } else {
                match app.residue_colors.get(&protein.atoms[i].residue_seq) {
                    Some(ss) => {
                        let (r, g, b) = ss.color_rgb();
                        [
                            (r as f32 * CONNECTION_DIM) as u8,
                            (g as f32 * CONNECTION_DIM) as u8,
                            (b as f32 * CONNECTION_DIM) as u8,
                        ]
                    }
                    None => [128, 128, 128],
                }
            };
            draw_segment(
                fb,
                &app.camera,
                protein.atoms[i].position,
                protein.atoms[j].position,
                color,
                NO_ID,
                0,
            );
        }
    }

    if app.show_hydrogen_bonds {
        for hbond in protein
            .hbonds
            .iter()
            .filter(|hb| hb.energy < app.hbond_energy_threshold)
        {
            let (donor, acceptor) = (hbond.donor_atom_idx, hbond.acceptor_atom_idx);
            if donor >= protein.atoms.len() || acceptor >= protein.atoms.len() {
                continue;
            }
            let (r, g, b) = hbond_energy_to_color(hbond.energy);
            draw_segment(
                fb,
                &app.camera,
                protein.atoms[donor].position,
                protein.atoms[acceptor].position,
                [r, g, b],
                NO_ID,
                HBOND_DASH,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::Vec3;
    use pixelfold_core::{Atom, FixedStr, SecondaryStructure};
    use pixelfold_render::renderer::auto_frame_protein;

    fn atom(name: &str, element: &str, pos: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(element),
            residue_name: FixedStr::new("ALA"),
            residue_seq: 1,
            chain_id: FixedStr::new("A"),
            is_hetatm: false,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: pos,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    fn protein(atoms: Vec<Atom>) -> Protein {
        Protein {
            atoms,
            title: String::new(),
            surface_points: Vec::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        }
    }

    fn framed_app(protein: Protein, w: u32, h: u32) -> App {
        let mut app = App::new();
        let mut camera = pixelfold_render::renderer::Camera::new();
        auto_frame_protein(&protein, &mut camera, w as f32, h as f32);
        app.camera = camera;
        app.protein = Some(protein);
        app.compute_static_geometry();
        app.show_depth_cue = false;
        app
    }

    #[test]
    fn backbone_mode_draws_only_calpha_atoms() {
        // A CA (index 0) and a side-chain CB (index 1) in one residue.
        let p = protein(vec![
            atom("CA", "C", Vec3::new(0.0, 0.0, 0.0)),
            atom("CB", "C", Vec3::new(1.5, 0.0, 0.0)),
        ]);
        let mut app = framed_app(p, 64, 64);
        app.display_mode = DisplayMode::Backbone;

        let fb = render_scene(&app, 64, 64);
        assert!(fb.ids().contains(&0), "the C-alpha must be drawn");
        assert!(
            !fb.ids().contains(&1),
            "backbone mode must not draw the side-chain atom"
        );
    }

    #[test]
    fn selected_atom_is_recoloured_cyan() {
        let p = protein(vec![atom("CA", "C", Vec3::ZERO)]);
        let mut app = framed_app(p, 48, 48);
        app.selected_atom_idx = Some(0);

        let fb = render_scene(&app, 48, 48);

        let mut covered = 0;
        for y in 0..fb.height() {
            for x in 0..fb.width() {
                if fb.id_at(x, y) == 0 {
                    covered += 1;
                    let [r, g, b, _] = fb.pixel(x, y);
                    // Cyan survives Lambert shading as red == 0 and green == blue.
                    assert_eq!(r, 0, "selected atom kept a red channel at ({x},{y})");
                    assert_eq!(g, b, "selected atom is not cyan at ({x},{y})");
                }
            }
        }
        assert!(covered > 0, "the selected atom drew no pixels");
    }
}
