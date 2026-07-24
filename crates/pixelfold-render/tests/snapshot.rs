//! Golden-image regression test for the rasteriser.
//!
//! A fixed synthetic scene is rendered at a fixed camera and compared against a
//! committed PNG. This needs no structure file, so it runs in CI. It guards the
//! whole pipeline (projection, occlusion, shading, depth cue, PNG encode) the
//! way the per-pixel unit tests guard the individual invariants.
//!
//! Regenerate the golden after an intended visual change with:
//!   PIXELFOLD_UPDATE_GOLDEN=1 cargo test -p pixelfold-render --test snapshot

use std::path::PathBuf;

use glam::Vec3;

use pixelfold_render::renderer::Camera;
use pixelfold_render::{Framebuffer, RenderOptions, Scene, decode_rgba, rasterize, to_png_bytes};

const SIZE: u32 = 128;
const BACKGROUND: [u8; 4] = [13, 13, 18, 255];

/// A pixel counts as changed when any channel moves more than this. The
/// rasteriser has no anti-aliasing, so a floating-point rounding difference at a
/// sphere rim flips a whole pixel; the tolerated-count check below absorbs a
/// handful of those while still catching real regressions.
const CHANNEL_TOLERANCE: i32 = 16;
const MAX_CHANGED_FRACTION: f64 = 0.01;

fn golden_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/synthetic.png")
}

/// A small deterministic scene: a central carbon flanked by oxygens and
/// nitrogens, a sulfur nearest to the viewer (it must occlude the carbon) and a
/// carbon farthest away (it takes the most fog).
fn synthetic_scene() -> Scene {
    let atoms: &[(&str, Vec3, [u8; 3], f32)] = &[
        ("C", Vec3::new(0.0, 0.0, 0.0), [144, 144, 144], 1.70),
        ("O", Vec3::new(1.6, 0.0, 0.4), [255, 13, 13], 1.52),
        ("N", Vec3::new(-1.6, 0.0, 0.4), [48, 80, 248], 1.55),
        ("O", Vec3::new(0.0, 1.6, -0.4), [255, 13, 13], 1.52),
        ("N", Vec3::new(0.0, -1.6, -0.4), [48, 80, 248], 1.55),
        ("S", Vec3::new(0.0, 0.0, 2.2), [255, 255, 48], 1.80),
        ("C", Vec3::new(0.0, 0.0, -2.4), [144, 144, 144], 1.70),
    ];

    Scene {
        positions: atoms.iter().map(|a| a.1).collect(),
        radii: atoms.iter().map(|a| a.3).collect(),
        colors: atoms.iter().map(|a| a.2).collect(),
        atom_idx: (0..atoms.len() as u32).collect(),
    }
}

fn render_synthetic() -> Framebuffer {
    let scene = synthetic_scene();
    let mut camera = Camera::new();
    camera.center = Vec3::ZERO;
    camera.zoom = 14.0;

    let mut fb = Framebuffer::new(SIZE, SIZE, BACKGROUND);
    rasterize(&scene, &camera, &mut fb, &RenderOptions::default());
    fb
}

#[test]
fn synthetic_scene_matches_the_golden_png() {
    let fb = render_synthetic();
    let rendered = to_png_bytes(&fb).expect("encode");

    if std::env::var_os("PIXELFOLD_UPDATE_GOLDEN").is_some() {
        std::fs::create_dir_all(golden_path().parent().unwrap()).expect("create golden dir");
        std::fs::write(golden_path(), &rendered).expect("write golden");
        return;
    }

    let golden_bytes = std::fs::read(golden_path()).unwrap_or_else(|_| {
        panic!(
            "golden missing; regenerate with PIXELFOLD_UPDATE_GOLDEN=1 cargo test -p \
             pixelfold-render --test snapshot"
        )
    });

    let (gw, gh, golden) = decode_rgba(&golden_bytes).expect("decode golden");
    assert_eq!((gw, gh), (SIZE, SIZE), "golden dimensions changed");

    let current = fb.color();
    assert_eq!(current.len(), golden.len(), "pixel count differs");

    let changed = current
        .iter()
        .zip(&golden)
        .filter(|(a, b)| (0..4).any(|c| (a[c] as i32 - b[c] as i32).abs() > CHANNEL_TOLERANCE))
        .count();

    let fraction = changed as f64 / current.len() as f64;
    assert!(
        fraction <= MAX_CHANGED_FRACTION,
        "{changed} of {} pixels changed ({:.2}%); render differs from the golden",
        current.len(),
        fraction * 100.0
    );

    let painted: Vec<&[u8; 4]> = current
        .iter()
        .zip(&golden)
        .filter(|(_, g)| **g != BACKGROUND)
        .map(|(c, _)| c)
        .collect();
    assert!(
        painted.len() > 500,
        "only {} pixels painted; scene looks empty",
        painted.len()
    );

    let mut total = 0u64;
    for (c, g) in current.iter().zip(&golden) {
        if *g != BACKGROUND {
            total += (0..4)
                .map(|i| (c[i] as i32 - g[i] as i32).unsigned_abs() as u64)
                .sum::<u64>();
        }
    }
    let mean = total as f64 / (painted.len() * 4) as f64;
    assert!(
        mean < 3.0,
        "mean per-channel drift over painted pixels is {mean:.2}; shading regressed"
    );
}
