//! Z-buffered sphere-impostor rasteriser.
//!
//! Each atom is a screen-space disc; the front hemisphere's depth is solved per
//! pixel and depth-tested, so occlusion is exact and per-pixel with no sort and
//! no painter's-algorithm bleeding. The surface normal falls out of the same
//! solve, giving Lambert shading for free, and an optional depth cue fades
//! distant atoms toward the background.

use std::f32::consts::FRAC_1_SQRT_2;

use glam::Vec3;

use crate::framebuffer::Framebuffer;
use crate::renderer::Camera;
use crate::scene::Scene;

/// Knobs the rasteriser exposes; [`Default`] is a sensible space-filling look.
#[derive(Clone, Copy, Debug)]
pub struct RenderOptions {
    /// Multiplies every van der Waals radius. 1.0 is full space-filling.
    pub radius_scale: f32,
    /// Fraction of a distant atom's colour replaced by the background. 0 is off.
    pub depth_cue: f32,
    /// Baseline lighting so faces pointing away from the light are not black.
    pub ambient: f32,
    /// Light direction in screen space (x right, y down, z toward the viewer).
    /// Normalised on use; the default is a headlight, so sphere centres are the
    /// brightest points.
    pub light: Vec3,
}

impl Default for RenderOptions {
    fn default() -> Self {
        Self {
            radius_scale: 1.0,
            depth_cue: 0.45,
            ambient: 0.15,
            light: Vec3::Z,
        }
    }
}

/// Rasterise every atom in `scene` into `fb` from `camera`'s viewpoint.
pub fn rasterize(scene: &Scene, camera: &Camera, fb: &mut Framebuffer, opts: &RenderOptions) {
    let width = fb.width() as f32;
    let height = fb.height() as f32;
    let zoom = camera.zoom;
    let light = opts.light.normalize_or_zero();
    let bg = fb.background();
    // Skip the depth-range pass entirely when the cue is off.
    let (zmin, fog_span, fog_on) = if opts.depth_cue > 0.0 {
        let (zmin, zmax) = depth_range(scene, camera, width, height);
        let span = zmax - zmin;
        (zmin, span, span > f32::EPSILON)
    } else {
        (0.0, 0.0, false)
    };

    for i in 0..scene.len() {
        let r_ang = scene.radii[i] * opts.radius_scale;
        let (sx, sy_up, zc) = camera.project_point_with_depth(scene.positions[i], width, height);
        // Screen y is lower-left; the framebuffer is top-left.
        let cx = sx;
        let cy = (height - 1.0) - sy_up;
        if !(cx.is_finite() && cy.is_finite() && zc.is_finite() && r_ang.is_finite())
            || r_ang <= 0.0
        {
            continue;
        }

        let r_px = r_ang * zoom;
        let color = scene.colors[i];
        let id = scene.atom_idx[i];

        // When the centre lands on integer pixel coordinates the nearest pixel
        // centre is 1/sqrt(2) away, so a disc smaller than that can cover no
        // pixel. Below the threshold, draw the centre directly so tiny atoms
        // never vanish.
        if r_px < FRAC_1_SQRT_2 {
            if let (Some(px), Some(py)) = (in_bounds(cx, fb.width()), in_bounds(cy, fb.height())) {
                let shaded = shade(color, 1.0, opts.ambient);
                let final_color = fog(shaded, bg, zc, zmin, fog_span, fog_on, opts.depth_cue);
                fb.test_and_set(px, py, zc, final_color, id);
            }
            continue;
        }

        let r_px2 = r_px * r_px;
        let x0 = (cx - r_px).floor().max(0.0) as i64;
        let x1 = (cx + r_px).ceil().min(width - 1.0) as i64;
        let y0 = (cy - r_px).floor().max(0.0) as i64;
        let y1 = (cy + r_px).ceil().min(height - 1.0) as i64;
        if x1 < x0 || y1 < y0 {
            continue;
        }

        for py in y0..=y1 {
            let dy = py as f32 + 0.5 - cy;
            for px in x0..=x1 {
                let dx = px as f32 + 0.5 - cx;
                let d2 = dx * dx + dy * dy;
                if d2 > r_px2 {
                    continue;
                }

                let dz_px = (r_px2 - d2).max(0.0).sqrt();
                let z = zc + dz_px / zoom;

                let n = Vec3::new(dx, dy, dz_px) / r_px;
                let intensity = n.dot(light).max(0.0);
                let shaded = shade(color, intensity, opts.ambient);
                let final_color = fog(shaded, bg, z, zmin, fog_span, fog_on, opts.depth_cue);

                fb.test_and_set(px as u32, py as u32, z, final_color, id);
            }
        }
    }
}

/// The near and far bounds of the scene's projected atom centres, used to scale
/// the depth cue. Non-finite projections are skipped.
fn depth_range(scene: &Scene, camera: &Camera, width: f32, height: f32) -> (f32, f32) {
    let mut zmin = f32::INFINITY;
    let mut zmax = f32::NEG_INFINITY;
    for &pos in &scene.positions {
        let (_, _, z) = camera.project_point_with_depth(pos, width, height);
        if z.is_finite() {
            zmin = zmin.min(z);
            zmax = zmax.max(z);
        }
    }
    if zmin.is_finite() {
        (zmin, zmax)
    } else {
        (0.0, 0.0)
    }
}

/// A whole-pixel coordinate inside `[0, extent)`, or `None` when off-screen.
fn in_bounds(coord: f32, extent: u32) -> Option<u32> {
    let rounded = coord.round();
    if rounded >= 0.0 && rounded < extent as f32 {
        Some(rounded as u32)
    } else {
        None
    }
}

fn shade(color: [u8; 3], intensity: f32, ambient: f32) -> [u8; 3] {
    let scale = (ambient + (1.0 - ambient) * intensity).clamp(0.0, 1.0);
    [
        (color[0] as f32 * scale) as u8,
        (color[1] as f32 * scale) as u8,
        (color[2] as f32 * scale) as u8,
    ]
}

#[allow(clippy::too_many_arguments)]
fn fog(
    color: [u8; 3],
    bg: [u8; 4],
    z: f32,
    zmin: f32,
    fog_span: f32,
    fog_on: bool,
    strength: f32,
) -> [u8; 4] {
    if !fog_on {
        return [color[0], color[1], color[2], 255];
    }

    // Farther atoms (smaller z) take more background.
    let t = (1.0 - (z - zmin) / fog_span).clamp(0.0, 1.0) * strength;
    [
        lerp(color[0], bg[0], t),
        lerp(color[1], bg[1], t),
        lerp(color[2], bg[2], t),
        255,
    ]
}

fn lerp(a: u8, b: u8, t: f32) -> u8 {
    (a as f32 * (1.0 - t) + b as f32 * t).round() as u8
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::framebuffer::NO_ID;
    use crate::scene::Coloring;
    use pixelfold_core::{Atom, FixedStr, Protein, SecondaryStructure};

    const BG: [u8; 4] = [0, 0, 0, 255];

    fn lit_opts() -> RenderOptions {
        // Fog off and no ambient floor keeps the pixel-level assertions exact.
        RenderOptions {
            radius_scale: 1.0,
            depth_cue: 0.0,
            ambient: 0.0,
            light: Vec3::Z,
        }
    }

    fn atom_at(element: &str, pos: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new("X"),
            element: FixedStr::new(element),
            residue_name: FixedStr::new("LIG"),
            residue_seq: 1,
            chain_id: FixedStr::new("A"),
            is_hetatm: true,
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

    /// A camera looking straight down -z with a chosen zoom, centred so a point
    /// at the origin lands at the middle of a `size`x`size` framebuffer.
    fn centred_camera(zoom: f32) -> Camera {
        let mut cam = Camera::new();
        cam.center = Vec3::ZERO;
        cam.zoom = zoom;
        cam
    }

    fn luminance(p: [u8; 4]) -> u32 {
        p[0] as u32 + p[1] as u32 + p[2] as u32
    }

    #[test]
    fn single_sphere_is_a_shaded_disc_brightest_at_the_centre() {
        let size = 41u32;
        let scene = Scene::from_protein(
            &protein(vec![atom_at("O", Vec3::ZERO)]),
            Coloring::Element,
            None,
        );
        let cam = centred_camera(8.0); // 1.52 A * 8 ~ 12 px radius
        let mut fb = Framebuffer::new(size, size, BG);
        rasterize(&scene, &cam, &mut fb, &lit_opts());

        // The projected centre, in framebuffer pixels (with the y-flip), is the
        // nearest point of the sphere and must be the brightest.
        let (sx, sy_up, _) = cam.project_point_with_depth(Vec3::ZERO, size as f32, size as f32);
        let cx = sx.floor() as u32;
        let cy = ((size as f32 - 1.0) - sy_up).floor() as u32;
        let centre = fb.pixel(cx, cy);
        assert!(luminance(centre) > 0, "centre must be painted");
        assert_ne!(fb.id_at(cx, cy), NO_ID);

        // A corner is outside the disc and keeps the background.
        assert_eq!(fb.pixel(0, 0), BG);
        assert_eq!(fb.id_at(0, 0), NO_ID);

        // No covered pixel is brighter than the centre, and no depth is NaN.
        for y in 0..size {
            for x in 0..size {
                assert!(
                    luminance(fb.pixel(x, y)) <= luminance(centre),
                    "pixel ({x},{y}) brighter than the centre"
                );
                if fb.id_at(x, y) != NO_ID {
                    assert!(fb.depth_at(x, y).is_finite(), "NaN depth at ({x},{y})");
                }
            }
        }

        // The disc dims from centre to rim: a mid-radius pixel is dimmer.
        let mid = fb.pixel(cx + 8, cy);
        assert!(luminance(mid) < luminance(centre));
        assert!(luminance(mid) > 0);
    }

    #[test]
    fn nearer_sphere_wins_every_contested_pixel() {
        // Two oxygens on the same screen spot, one 100 A nearer the viewer.
        let far = atom_at("O", Vec3::new(0.0, 0.0, -50.0));
        let near = atom_at("N", Vec3::new(0.0, 0.0, 50.0));
        let cam = centred_camera(10.0);
        let size = 61u32;

        let scene = Scene::from_protein(
            &protein(vec![far.clone(), near.clone()]),
            Coloring::Element,
            None,
        );
        let mut fb = Framebuffer::new(size, size, BG);
        rasterize(&scene, &cam, &mut fb, &lit_opts());

        let near_id = 1u32; // second atom
        let mut contested = 0;
        for y in 0..size {
            for x in 0..size {
                let id = fb.id_at(x, y);
                if id != NO_ID {
                    assert_eq!(id, near_id, "far atom leaked through at ({x},{y})");
                    assert!(fb.depth_at(x, y).is_finite());
                    contested += 1;
                }
            }
        }
        assert!(contested > 0, "the near sphere must cover pixels");
    }

    #[test]
    fn occlusion_is_independent_of_draw_order() {
        let far = atom_at("O", Vec3::new(0.0, 0.0, -50.0));
        let near = atom_at("N", Vec3::new(0.5, 0.0, 50.0));
        let cam = centred_camera(10.0);
        let size = 61u32;

        let forward = Scene::from_protein(
            &protein(vec![far.clone(), near.clone()]),
            Coloring::Element,
            None,
        );
        let reversed = Scene::from_protein(&protein(vec![near, far]), Coloring::Element, None);

        let mut a = Framebuffer::new(size, size, BG);
        let mut b = Framebuffer::new(size, size, BG);
        rasterize(&forward, &cam, &mut a, &lit_opts());
        rasterize(&reversed, &cam, &mut b, &lit_opts());

        // The winning fragment is the same regardless of which atom drew first;
        // only the atom ids swap.
        for y in 0..size {
            for x in 0..size {
                assert_eq!(a.pixel(x, y), b.pixel(x, y), "colour differs at ({x},{y})");
                assert_eq!(
                    a.depth_at(x, y),
                    b.depth_at(x, y),
                    "depth differs at ({x},{y})"
                );
            }
        }
    }

    #[test]
    fn empty_scene_leaves_the_background() {
        let scene = Scene::from_protein(&protein(vec![]), Coloring::Element, None);
        let cam = centred_camera(10.0);
        let mut fb = Framebuffer::new(8, 8, BG);
        rasterize(&scene, &cam, &mut fb, &RenderOptions::default());
        assert!(fb.color().iter().all(|&p| p == BG));
        assert!(fb.ids().iter().all(|&id| id == NO_ID));
    }

    #[test]
    fn non_finite_positions_are_skipped_without_panicking() {
        let scene = Scene::from_protein(
            &protein(vec![atom_at("C", Vec3::new(f32::NAN, 0.0, 0.0))]),
            Coloring::Element,
            None,
        );
        let cam = centred_camera(10.0);
        let mut fb = Framebuffer::new(16, 16, BG);
        rasterize(&scene, &cam, &mut fb, &RenderOptions::default());
        assert!(fb.color().iter().all(|&p| p == BG));
    }

    #[test]
    fn small_radius_atoms_never_vanish() {
        // An even framebuffer projects an origin atom to integer pixel
        // coordinates, the worst case for the disc coverage test: the nearest
        // pixel centre is 1/sqrt(2) away.
        let size = 40u32;
        let scene = Scene::from_protein(
            &protein(vec![atom_at("C", Vec3::ZERO)]),
            Coloring::Element,
            None,
        );
        for &zoom in &[0.2_f32, 0.35, 0.41] {
            let cam = centred_camera(zoom); // r_px = 1.70 * zoom, i.e. 0.34, 0.60, 0.70
            let mut fb = Framebuffer::new(size, size, BG);
            rasterize(&scene, &cam, &mut fb, &lit_opts());

            let painted: Vec<u32> = (0..size * size)
                .filter(|&i| fb.ids()[i as usize] != NO_ID)
                .collect();
            assert!(
                !painted.is_empty(),
                "atom vanished at r_px = {}",
                1.70 * zoom
            );
            for &i in &painted {
                let (x, y) = (i % size, i / size);
                assert_eq!(fb.id_at(x, y), 0);
                assert!(fb.depth_at(x, y).is_finite());
            }
        }
    }

    #[test]
    fn depth_cue_fades_distant_atoms_toward_the_background() {
        // With the cue on, the far atom is pulled harder toward the background,
        // so its brightest pixel is dimmer (background here is black).
        let near = atom_at("O", Vec3::new(-2.0, 0.0, 40.0));
        let far = atom_at("O", Vec3::new(2.0, 0.0, -40.0));
        let cam = centred_camera(10.0);
        let size = 81u32;
        let scene = Scene::from_protein(&protein(vec![near, far]), Coloring::Element, None);

        let mut fb = Framebuffer::new(size, size, BG);
        let opts = RenderOptions {
            radius_scale: 1.0,
            depth_cue: 0.6,
            ambient: 0.15,
            light: Vec3::Z,
        };
        rasterize(&scene, &cam, &mut fb, &opts);

        let brightest = |xs: std::ops::Range<u32>| {
            let mut best = 0;
            for y in 0..size {
                for x in xs.clone() {
                    best = best.max(luminance(fb.pixel(x, y)));
                }
            }
            best
        };
        let near_side = brightest(0..size / 2);
        let far_side = brightest(size / 2 + 1..size);
        assert!(
            near_side > far_side,
            "near atom ({near_side}) should outshine the fogged far atom ({far_side})"
        );
    }
}
