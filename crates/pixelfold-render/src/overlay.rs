//! Depth-aware line overlays drawn into the framebuffer after the spheres.
//!
//! Bonds and backbone connections are thin lines that must sit in the same
//! z-buffer as the atoms, so a line is hidden wherever an atom is in front of
//! it. Each pixel's depth is interpolated between the projected endpoints and
//! goes through the same [`Framebuffer::test_and_set`] the rasteriser uses.

use glam::Vec3;

use crate::draw::{draw_dashed_line, draw_line};
use crate::framebuffer::Framebuffer;
use crate::renderer::Camera;

/// Draw a world-space segment from `a` to `b` into `fb`, depth-tested against
/// whatever the rasteriser already wrote. `dash` greater than zero draws a
/// dashed line with that spacing; zero draws it solid. Pixels are tagged with
/// `id` (pass [`crate::framebuffer::NO_ID`] for a non-pickable overlay).
pub fn draw_segment(
    fb: &mut Framebuffer,
    camera: &Camera,
    a: Vec3,
    b: Vec3,
    color: [u8; 3],
    id: u32,
    dash: usize,
) {
    let width = fb.width() as f32;
    let height = fb.height() as f32;
    let (ax, ay_up, az) = camera.project_point_with_depth(a, width, height);
    let (bx, by_up, bz) = camera.project_point_with_depth(b, width, height);
    if !(ax.is_finite()
        && ay_up.is_finite()
        && az.is_finite()
        && bx.is_finite()
        && by_up.is_finite()
        && bz.is_finite())
    {
        return;
    }

    // The camera's screen y points up; the framebuffer's origin is top-left,
    // the same flip the rasteriser applies.
    let ay = (height - 1.0) - ay_up;
    let by = (height - 1.0) - by_up;

    let points = if dash > 0 {
        draw_dashed_line(ax, ay, bx, by, dash)
    } else {
        draw_line(ax, ay, bx, by)
    };
    if points.is_empty() {
        return;
    }

    let rgba = [color[0], color[1], color[2], 255];
    let last = (points.len() - 1).max(1) as f32;
    for (i, (px, py)) in points.into_iter().enumerate() {
        let t = i as f32 / last;
        let z = az + (bz - az) * t;
        let x = px.round();
        let y = py.round();
        if x < 0.0 || y < 0.0 || x >= width as f64 || y >= height as f64 {
            continue;
        }

        fb.test_and_set(x as u32, y as u32, z, rgba, id);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const BG: [u8; 4] = [0, 0, 0, 255];

    fn centred_camera(zoom: f32) -> Camera {
        let mut cam = Camera::new();
        cam.center = Vec3::ZERO;
        cam.zoom = zoom;
        cam
    }

    #[test]
    fn draws_pixels_along_a_segment() {
        let mut fb = Framebuffer::new(32, 32, BG);
        let cam = centred_camera(2.0);
        draw_segment(
            &mut fb,
            &cam,
            Vec3::new(-5.0, 0.0, 0.0),
            Vec3::new(5.0, 0.0, 0.0),
            [0, 255, 0],
            7,
            0,
        );
        let painted = fb.ids().iter().filter(|&&id| id == 7).count();
        assert!(painted > 0, "the segment drew no pixels");
    }

    #[test]
    fn nearer_fragments_occlude_the_segment() {
        let mut fb = Framebuffer::new(32, 32, BG);
        // A wall in front of everything (large depth is nearer).
        for y in 0..32 {
            for x in 0..32 {
                fb.test_and_set(x, y, 1000.0, [255, 0, 0, 255], 1);
            }
        }
        let cam = centred_camera(2.0);
        // The segment sits at view depth 0, well behind the wall.
        draw_segment(
            &mut fb,
            &cam,
            Vec3::new(-5.0, 0.0, 0.0),
            Vec3::new(5.0, 0.0, 0.0),
            [0, 255, 0],
            7,
            0,
        );
        assert!(
            !fb.ids().contains(&7),
            "an occluded segment must not overwrite nearer fragments"
        );
    }
}
