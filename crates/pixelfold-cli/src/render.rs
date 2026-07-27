//! The `render` subcommand: a structure to a PNG file or straight to the terminal.

use std::io::{self, IsTerminal, Write};
use std::path::Path;

use anyhow::{Context, Result};
use pixelfold_core::{Protein, get_calpha_indices};
use pixelfold_render::renderer::{Camera, auto_frame_points};
use pixelfold_render::sink::{self, Protocol};
use pixelfold_render::{
    Coloring, Framebuffer, RenderOptions, Scene, rasterize, to_png_bytes, write_png,
};

/// Largest image dimension accepted, bounding the framebuffer allocation.
const MAX_DIMENSION: u32 = 8192;

const BACKGROUND: [u8; 4] = [13, 13, 18, 255];
const FOG_STRENGTH: f32 = 0.45;

/// Defaults for a file, where there is no terminal to fit.
const DEFAULT_WIDTH: u32 = 1200;
const DEFAULT_HEIGHT: u32 = 900;

/// Approximate character-cell size in pixels, used to fit graphics-protocol
/// output to the terminal. Fonts vary; this only needs to be close.
const CELL_W: u32 = 8;
const CELL_H: u32 = 16;

/// The terminal protocol to use when no output file is given.
#[derive(Clone, Copy, Debug)]
pub enum TerminalTarget {
    Auto,
    Kitty,
    Iterm2,
    HalfBlock,
}

pub struct RenderParams {
    /// `None` fits the terminal (or the file default).
    pub width: Option<u32>,
    pub height: Option<u32>,
    pub coloring: Coloring,
    pub backbone: bool,
    pub radius_scale: f32,
    pub depth_cue: bool,
}

/// Where the render goes, resolved once so the framebuffer can be sized for it.
enum Destination<'a> {
    File(&'a Path),
    /// A graphics or text escape written to the terminal.
    Terminal(Protocol),
    /// A raw PNG to a non-terminal stdout, so `render X > out.png` works.
    PipePng,
}

/// Render `indices` of `protein` (all atoms when `None`). Writes a PNG to
/// `output` when given, otherwise emits to the terminal (or a piped PNG).
pub fn render(
    protein: &Protein,
    indices: Option<&[usize]>,
    params: &RenderParams,
    output: Option<&Path>,
    target: TerminalTarget,
) -> Result<()> {
    let shown = shown_indices(protein, indices, params.backbone);
    if shown.is_empty() {
        anyhow::bail!("nothing to render (the selection matched no displayable atoms)");
    }

    let dest = destination(output, target);
    let (width, height) = effective_dims(&dest, params);
    validate(width, height, params.radius_scale)?;

    let scene = Scene::from_protein(protein, params.coloring, Some(&shown));
    let mut camera = Camera::new();
    auto_frame_points(&scene.positions, &mut camera, width as f32, height as f32);

    let mut fb = Framebuffer::new(width, height, BACKGROUND);
    let opts = RenderOptions {
        radius_scale: params.radius_scale,
        depth_cue: if params.depth_cue { FOG_STRENGTH } else { 0.0 },
        ..RenderOptions::default()
    };
    rasterize(&scene, &camera, &mut fb, &opts);

    match dest {
        Destination::File(path) => {
            write_png(&fb, path).with_context(|| format!("could not write {}", path.display()))?;
            eprintln!(
                "wrote {} ({width}x{height}, {} atoms)",
                path.display(),
                shown.len()
            );

            Ok(())
        }
        Destination::PipePng => {
            let mut stdout = io::stdout().lock();
            let bytes = to_png_bytes(&fb).context("encoding the image as PNG")?;
            write_all(&mut stdout, &bytes)
        }
        Destination::Terminal(protocol) => {
            let mut stdout = io::stdout().lock();
            let bytes =
                sink::to_terminal(&fb, protocol).context("encoding the image for the terminal")?;

            write_all(&mut stdout, &bytes)?;
            write_all(&mut stdout, b"\n")?;
            flush(&mut stdout)
        }
    }
}

/// Resolve where output goes. An explicit protocol is always honored (so its
/// escape can be captured to a file); `Auto` detects the terminal when stdout is
/// one, and otherwise pipes a PNG.
fn destination(output: Option<&Path>, target: TerminalTarget) -> Destination<'_> {
    if let Some(path) = output {
        return Destination::File(path);
    }

    match target {
        TerminalTarget::Kitty => Destination::Terminal(Protocol::Kitty),
        TerminalTarget::Iterm2 => Destination::Terminal(Protocol::Iterm2),
        TerminalTarget::HalfBlock => Destination::Terminal(Protocol::HalfBlock),
        TerminalTarget::Auto if io::stdout().is_terminal() => {
            Destination::Terminal(sink::detect(|k| std::env::var(k).ok()))
        }
        TerminalTarget::Auto => Destination::PipePng,
    }
}

/// The render dimensions: an explicit `--width`/`--height` wins per axis;
/// otherwise a file uses the fixed defaults and a terminal is fitted to its size.
fn effective_dims(dest: &Destination, params: &RenderParams) -> (u32, u32) {
    match dest {
        Destination::File(_) | Destination::PipePng => (
            params.width.unwrap_or(DEFAULT_WIDTH),
            params.height.unwrap_or(DEFAULT_HEIGHT),
        ),
        Destination::Terminal(protocol) => {
            let (cols, rows) = crossterm::terminal::size().unwrap_or((80, 24));
            let (fit_w, fit_h) = fit_dims(*protocol, cols, rows);
            (
                params.width.unwrap_or(fit_w),
                params.height.unwrap_or(fit_h),
            )
        }
    }
}

/// Dimensions that fill the terminal's text area. Half blocks pack two pixels
/// per row and one per column; the graphics protocols display a PNG, so the cell
/// size is approximated. One row is left for the prompt, and the result is
/// clamped so auto-fit never trips the dimension limit.
fn fit_dims(protocol: Protocol, cols: u16, rows: u16) -> (u32, u32) {
    let cols = cols.max(1) as u32;
    let rows = rows.saturating_sub(1).max(1) as u32;
    let (w, h) = match protocol {
        Protocol::HalfBlock => (cols, rows * 2),
        Protocol::Kitty | Protocol::Iterm2 => (cols * CELL_W, rows * CELL_H),
    };

    (w.min(MAX_DIMENSION), h.min(MAX_DIMENSION))
}

fn validate(width: u32, height: u32, radius_scale: f32) -> Result<()> {
    if width == 0 || height == 0 {
        anyhow::bail!("image dimensions must be non-zero");
    }
    if width > MAX_DIMENSION || height > MAX_DIMENSION {
        anyhow::bail!("image dimensions must not exceed {MAX_DIMENSION}");
    }
    if !(radius_scale.is_finite() && radius_scale > 0.0) {
        anyhow::bail!("radius scale must be a positive number");
    }

    Ok(())
}

/// Write everything, treating a closed pipe as success.
fn write_all(out: &mut impl Write, bytes: &[u8]) -> Result<()> {
    match out.write_all(bytes) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == io::ErrorKind::BrokenPipe => Ok(()),
        Err(error) => Err(error).context("writing to standard output"),
    }
}

fn flush(out: &mut impl Write) -> Result<()> {
    match out.flush() {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == io::ErrorKind::BrokenPipe => Ok(()),
        Err(error) => Err(error).context("flushing standard output"),
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn halfblock_fit_fills_columns_and_double_rows() {
        // One pixel per column, two per row, leaving a row for the prompt.
        assert_eq!(fit_dims(Protocol::HalfBlock, 80, 24), (80, 46));
    }

    #[test]
    fn graphics_fit_approximates_cell_pixels() {
        assert_eq!(
            fit_dims(Protocol::Kitty, 80, 24),
            (80 * CELL_W, 23 * CELL_H)
        );
    }

    #[test]
    fn fit_never_exceeds_the_dimension_limit() {
        let (w, h) = fit_dims(Protocol::Kitty, u16::MAX, u16::MAX);
        assert!(w <= MAX_DIMENSION && h <= MAX_DIMENSION);
    }
}
