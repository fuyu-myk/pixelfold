use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{Context, Result};
use pixelfold_core::{
    AltlocPolicy, DisplayMode, Protein, SecondaryStructure, get_calpha_connections,
};
use pixelfold_render::Framebuffer;
use pixelfold_render::renderer::{self, Camera};
use ratatui::layout::Rect;
use ratatui_image::picker::Picker;

pub mod search;

mod app;
mod inputs;
mod iterm2;
mod network;
mod render;
mod ui;

/// The framebuffer of the last drawn frame, kept so a mouse click can be mapped
/// straight to the atom under it through the id buffer. `area` is where the image
/// was placed; the framebuffer may be a lower resolution that the terminal scaled
/// to fill it, so picking maps proportionally rather than by a fixed cell size.
pub struct RenderedFrame {
    pub fb: Framebuffer,
    pub area: Rect,
}

pub struct App {
    pub protein: Option<Protein>,
    pub camera: Camera,
    pub redraw_needed: bool,
    pub ca_indices: Vec<usize>,
    pub backbone_connections: Vec<(usize, usize)>,
    pub residue_colors: HashMap<u32, SecondaryStructure>,
    pub inspect_mode: bool,
    pub residue_highlight: bool,
    pub selected_atom_idx: Option<usize>,
    /// Every atom sharing the selected atom's chain and residue number.
    pub highlighted_atom_indices: HashSet<usize>,
    pub display_mode: DisplayMode,
    pub show_connections: bool,
    pub use_bfactor_colors: bool,
    pub show_hydrogen_bonds: bool,
    pub hbond_energy_threshold: f32,
    /// Multiplies van der Waals radii; 1.0 is full space-filling.
    pub radius_scale: f32,
    pub show_depth_cue: bool,
    /// Depth-band clip `(far, near)` in `0..1`, or `None` for the whole depth.
    pub slab: Option<(f32, f32)>,
    pub last_frame: Option<RenderedFrame>,
    /// The linked network pane: shown beside the structure, built on first open.
    pub show_network: bool,
    pub(crate) network_view: Option<network::NetworkView>,
    /// A network build running on a worker thread, delivering its result here.
    /// Set while the pane is computing to prevent blocking of the event loop.
    pub(crate) network_build: Option<std::sync::mpsc::Receiver<network::NetworkView>>,
    /// Where the network pane last drew, for routing clicks to it.
    pub network_area: Option<Rect>,
}

impl Default for App {
    fn default() -> Self {
        Self::new()
    }
}

impl App {
    pub fn new() -> Self {
        Self {
            protein: None,
            camera: Camera::new(),
            redraw_needed: true,
            ca_indices: Vec::new(),
            backbone_connections: Vec::new(),
            residue_colors: HashMap::new(),
            inspect_mode: false,
            residue_highlight: false,
            selected_atom_idx: None,
            highlighted_atom_indices: HashSet::new(),
            display_mode: DisplayMode::AllAtoms,
            show_connections: true,
            use_bfactor_colors: false,
            show_hydrogen_bonds: false,
            hbond_energy_threshold: -0.5,
            radius_scale: 1.0,
            show_depth_cue: true,
            slab: None,
            last_frame: None,
            show_network: false,
            network_view: None,
            network_build: None,
            network_area: None,
        }
    }

    pub fn load_protein(
        &mut self,
        path: &Path,
        width: f32,
        height: f32,
        skip_surface: bool,
        altloc: AltlocPolicy,
    ) -> Result<()> {
        self.protein = Some(pixelfold_core::parser::load_protein_with_options(
            path,
            skip_surface,
            altloc,
        )?);

        if let Some(ref protein) = self.protein {
            renderer::auto_frame_protein(protein, &mut self.camera, width, height);
            self.compute_static_geometry();
        }

        Ok(())
    }

    /// Pre-compute geometry that does not change as the camera moves.
    pub fn compute_static_geometry(&mut self) {
        if let Some(ref protein) = self.protein {
            self.ca_indices = protein
                .atoms
                .iter()
                .enumerate()
                .filter(|(_, atom)| atom.name == "CA")
                .map(|(idx, _)| idx)
                .collect();

            self.backbone_connections = get_calpha_connections(protein, &self.ca_indices);

            self.residue_colors.clear();
            for atom in &protein.atoms {
                self.residue_colors
                    .entry(atom.residue_seq)
                    .or_insert(atom.secondary_structure);
            }
        }
    }
}

/// Load a structure and run the interactive viewer.
pub fn view(path: &Path, skip_surface: bool, altloc: AltlocPolicy) -> Result<()> {
    let mut terminal = ratatui::init();

    // Query the terminal for its graphics protocol and font.
    let picker = Picker::from_query_stdio().unwrap_or_else(|_| Picker::halfblocks());

    let result = run_view(&mut terminal, &picker, path, skip_surface, altloc);

    let _ = crossterm::execute!(std::io::stdout(), crossterm::event::DisableMouseCapture);
    ratatui::restore();

    result
}

fn run_view(
    terminal: &mut ratatui::DefaultTerminal,
    picker: &Picker,
    path: &Path,
    skip_surface: bool,
    altloc: AltlocPolicy,
) -> Result<()> {
    crossterm::execute!(std::io::stdout(), crossterm::event::EnableMouseCapture)?;

    let (cols, rows) =
        crossterm::terminal::size().context("failed to query terminal size (not a terminal?)")?;
    let (font_w, font_h) = picker.font_size();
    let width = (cols as u32 * font_w as u32).max(1) as f32;
    let height = (rows as u32 * font_h as u32).max(1) as f32;

    let mut app = App::new();
    app.load_protein(path, width, height, skip_surface, altloc)?;

    crate::app::run_app(terminal, &mut app, picker)
}

/// Run the RCSB search-and-fetch interface, caching downloads to `cache_dir`.
pub fn search(query: Option<String>, cache_dir: std::path::PathBuf) -> Result<()> {
    with_terminal(|terminal| crate::search::fetch_structures(terminal, query, cache_dir))
}

/// Set up the terminal (alternate screen + mouse capture), run `f`, then restore.
fn with_terminal<F>(f: F) -> Result<()>
where
    F: FnOnce(&mut ratatui::DefaultTerminal) -> Result<()>,
{
    let mut terminal = ratatui::init();
    crossterm::execute!(std::io::stdout(), crossterm::event::EnableMouseCapture)?;

    let result = f(&mut terminal);
    let _ = crossterm::execute!(std::io::stdout(), crossterm::event::DisableMouseCapture);

    ratatui::restore();
    result
}

/// Recompute the highlighted set as every atom in the selected atom's residue
/// (same chain and residue number). Cleared when nothing is selected.
pub fn update_highlighted_atoms(app: &mut App) {
    app.highlighted_atom_indices.clear();

    let (Some(protein), Some(selected)) = (&app.protein, app.selected_atom_idx) else {
        return;
    };

    let selected_atom = &protein.atoms[selected];
    let chain = selected_atom.chain_id;
    let residue_seq = selected_atom.residue_seq;

    for (idx, atom) in protein.atoms.iter().enumerate() {
        if atom.chain_id == chain && atom.residue_seq == residue_seq {
            app.highlighted_atom_indices.insert(idx);
        }
    }
}

/// The atom under a terminal cell, via the id buffer of the last drawn frame.
/// Returns `None` for a click outside the image or on background. The click's
/// position within the image area is mapped proportionally onto the framebuffer,
/// so it is correct whatever resolution that frame was drawn at.
pub fn pick_at(app: &App, column: u16, row: u16) -> Option<usize> {
    let frame = app.last_frame.as_ref()?;
    let area = frame.area;
    if column < area.x || row < area.y || area.width == 0 || area.height == 0 {
        return None;
    }

    let local_col = (column - area.x) as u32;
    let local_row = (row - area.y) as u32;
    if local_col >= area.width as u32 || local_row >= area.height as u32 {
        return None;
    }

    let fx = (local_col as f32 + 0.5) / area.width as f32;
    let fy = (local_row as f32 + 0.5) / area.height as f32;
    let px = ((fx * frame.fb.width() as f32) as u32).min(frame.fb.width().saturating_sub(1));
    let py = ((fy * frame.fb.height() as f32) as u32).min(frame.fb.height().saturating_sub(1));

    match frame.fb.id_at(px, py) {
        pixelfold_render::NO_ID => None,
        id => Some(id as usize),
    }
}
