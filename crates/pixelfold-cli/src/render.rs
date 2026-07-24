//! The `render` subcommand: a structure to a PNG, headless.

use std::path::Path;

use anyhow::{Context, Result};
use pixelfold_core::{Protein, get_calpha_indices};
use pixelfold_render::renderer::{Camera, auto_frame_points};
use pixelfold_render::{Coloring, Framebuffer, RenderOptions, Scene, rasterize, write_png};

/// Largest image dimension accepted, bounding the framebuffer allocation.
const MAX_DIMENSION: u32 = 8192;

const BACKGROUND: [u8; 4] = [13, 13, 18, 255];
const FOG_STRENGTH: f32 = 0.45;

pub struct RenderParams {
    pub width: u32,
    pub height: u32,
    pub coloring: Coloring,
    pub backbone: bool,
    pub radius_scale: f32,
    pub depth_cue: bool,
}

impl RenderParams {
    fn validate(&self) -> Result<()> {
        if self.width == 0 || self.height == 0 {
            anyhow::bail!("image dimensions must be non-zero");
        }
        if self.width > MAX_DIMENSION || self.height > MAX_DIMENSION {
            anyhow::bail!("image dimensions must not exceed {MAX_DIMENSION}");
        }
        if !(self.radius_scale.is_finite() && self.radius_scale > 0.0) {
            anyhow::bail!("radius scale must be a positive number");
        }

        Ok(())
    }
}

/// Render `indices` of `protein` (all atoms when `None`) to a PNG at `out`.
pub fn render(
    protein: &Protein,
    indices: Option<&[usize]>,
    params: &RenderParams,
    out: &Path,
) -> Result<()> {
    params.validate()?;

    let shown = shown_indices(protein, indices, params.backbone);
    if shown.is_empty() {
        anyhow::bail!("nothing to render (the selection matched no displayable atoms)");
    }

    let scene = Scene::from_protein(protein, params.coloring, Some(&shown));

    let mut camera = Camera::new();
    auto_frame_points(
        &scene.positions,
        &mut camera,
        params.width as f32,
        params.height as f32,
    );

    let mut fb = Framebuffer::new(params.width, params.height, BACKGROUND);
    let opts = RenderOptions {
        radius_scale: params.radius_scale,
        depth_cue: if params.depth_cue { FOG_STRENGTH } else { 0.0 },
        ..RenderOptions::default()
    };
    rasterize(&scene, &camera, &mut fb, &opts);

    write_png(&fb, out).with_context(|| format!("could not write {}", out.display()))?;
    eprintln!(
        "wrote {} ({}x{}, {} atoms)",
        out.display(),
        params.width,
        params.height,
        shown.len()
    );

    Ok(())
}

/// The atom indices to draw: the C-alpha trace or every atom, intersected with a
/// selection when one is given.
fn shown_indices(protein: &Protein, selection: Option<&[usize]>, backbone: bool) -> Vec<usize> {
    let base: Vec<usize> = if backbone {
        get_calpha_indices(protein)
    } else {
        (0..protein.atoms.len()).collect()
    };

    match selection {
        None => base,
        Some(selection) => {
            let keep: std::collections::HashSet<usize> = selection.iter().copied().collect();
            base.into_iter().filter(|i| keep.contains(i)).collect()
        }
    }
}
