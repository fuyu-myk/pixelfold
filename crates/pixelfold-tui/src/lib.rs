use std::collections::HashSet;
use std::path::Path;

use anyhow::{Context, Result};
use pixelfold_core::parser::ParsedStructure;
use pixelfold_core::{AltlocPolicy, DisplayMode, Protein, get_calpha_connections};
use pixelfold_render::Framebuffer;
use pixelfold_render::renderer::{self, Camera};
use ratatui::layout::Rect;
use ratatui_image::picker::{Picker, ProtocolType};

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
    pub inspect_mode: bool,
    pub residue_highlight: bool,
    pub selected_atom_idx: Option<usize>,
    /// Every atom sharing the selected atom's chain, residue number, and
    /// insertion code.
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
    /// The parsed file the shown model came from, kept so stepping between
    /// models never re-reads it.
    pub(crate) source: Option<ParsedStructure>,
    pub(crate) altloc: AltlocPolicy,
    /// The file's model numbers, and the position of the shown one among them.
    pub(crate) models: Vec<usize>,
    pub(crate) model_position: usize,
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
            source: None,
            altloc: AltlocPolicy::default(),
            models: Vec::new(),
            model_position: 0,
        }
    }

    /// Load model `model` of the file at `path`, or its first model when `None`,
    /// and frame it.
    pub fn load_protein(
        &mut self,
        path: &Path,
        width: f32,
        height: f32,
        altloc: AltlocPolicy,
        model: Option<usize>,
    ) -> Result<()> {
        let source = ParsedStructure::read(path)?;
        let protein = source.protein(model, altloc)?;
        self.models = source.model_numbers();
        self.model_position = model
            .and_then(|number| self.models.iter().position(|&m| m == number))
            .unwrap_or(0);
        self.source = Some(source);
        self.altloc = altloc;

        renderer::auto_frame_protein(&protein, &mut self.camera, width, height);
        self.protein = Some(protein);
        self.compute_static_geometry();

        Ok(())
    }

    /// Show the next model of a multi-model file (`delta` 1) or the previous one
    /// (-1), wrapping at either end. The camera stays put; the selection and the
    /// network pane, which index the old model's atoms, start over.
    pub fn step_model(&mut self, delta: isize) -> Result<()> {
        let count = self.models.len();
        let Some(source) = self.source.as_ref().filter(|_| count > 1) else {
            return Ok(());
        };

        let position = (self.model_position as isize + delta).rem_euclid(count as isize) as usize;
        self.protein = Some(source.protein(Some(self.models[position]), self.altloc)?);
        self.model_position = position;
        self.compute_static_geometry();

        self.selected_atom_idx = None;
        self.residue_highlight = false;
        self.highlighted_atom_indices.clear();
        self.network_view = None;
        self.network_build = None;
        if self.show_network {
            self.start_network_build();
        }
        self.redraw_needed = true;

        Ok(())
    }

    /// Build the network pane's view on a worker thread, so a large structure's
    /// detection, layout, and accessibility pass never freeze the event loop;
    /// the pane shows a placeholder until the result arrives (see `run_app`).
    pub(crate) fn start_network_build(&mut self) {
        let Some(protein) = self.protein.clone() else {
            return;
        };

        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let view = network::NetworkView::build(&protein);
            tx.send(view).ok();
        });
        self.network_build = Some(rx);
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
        }
    }
}

/// Which terminal graphics protocol the viewer draws the structure with.
///
/// Detection queries the terminal, which can be wrong. For instance, an emulator
/// that claims a protocol it does not paint leaves the structure pane blank,
/// with the rest of the interface intact. `HalfBlock` renders are forced on any
/// truecolor terminal, to overcome this limitation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ViewProtocol {
    #[default]
    Auto,
    Kitty,
    Iterm2,
    HalfBlock,
}

/// Load a structure and run the interactive viewer.
pub fn view(
    path: &Path,
    altloc: AltlocPolicy,
    model: Option<usize>,
    protocol: ViewProtocol,
) -> Result<()> {
    let mut terminal = ratatui::init();

    let picker = build_picker(protocol);

    let result = run_view(&mut terminal, &picker, path, altloc, model);

    let _ = crossterm::execute!(std::io::stdout(), crossterm::event::DisableMouseCapture);
    ratatui::restore();

    result
}

/// The picker the structure pane draws through. `Auto` queries the terminal and
/// falls back to half blocks when the query fails; the explicit modes keep the
/// queried font size but override what the terminal claimed to support.
fn build_picker(protocol: ViewProtocol) -> Picker {
    let queried = Picker::from_query_stdio();
    let forced = match protocol {
        ViewProtocol::Auto => return queried.unwrap_or_else(|_| Picker::halfblocks()),
        ViewProtocol::Kitty => ProtocolType::Kitty,
        ViewProtocol::Iterm2 => ProtocolType::Iterm2,
        ViewProtocol::HalfBlock => ProtocolType::Halfblocks,
    };

    let mut picker = queried.unwrap_or_else(|_| Picker::halfblocks());
    picker.set_protocol_type(forced);

    picker
}

fn run_view(
    terminal: &mut ratatui::DefaultTerminal,
    picker: &Picker,
    path: &Path,
    altloc: AltlocPolicy,
    model: Option<usize>,
) -> Result<()> {
    crossterm::execute!(std::io::stdout(), crossterm::event::EnableMouseCapture)?;

    let (cols, rows) =
        crossterm::terminal::size().context("failed to query terminal size (not a terminal?)")?;
    let (font_w, font_h) = picker.font_size();
    let width = (cols as u32 * font_w as u32).max(1) as f32;
    let height = (rows as u32 * font_h as u32).max(1) as f32;

    let mut app = App::new();
    app.load_protein(path, width, height, altloc, model)?;

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
/// (same chain, residue number, and insertion code). Cleared when nothing is
/// selected.
pub fn update_highlighted_atoms(app: &mut App) {
    app.highlighted_atom_indices.clear();

    let (Some(protein), Some(selected)) = (&app.protein, app.selected_atom_idx) else {
        return;
    };

    let selected_atom = &protein.atoms[selected];
    let chain = selected_atom.chain_id;
    let residue_seq = selected_atom.residue_seq;
    let insertion_code = selected_atom.insertion_code;

    for (idx, atom) in protein.atoms.iter().enumerate() {
        if atom.chain_id == chain
            && atom.residue_seq == residue_seq
            && atom.insertion_code == insertion_code
        {
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

#[cfg(test)]
mod tests {
    use super::*;
    use glam::Vec3;
    use pixelfold_core::{Atom, FixedStr, SecondaryStructure};

    fn atom(name: &str, insertion_code: Option<char>) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new("C"),
            residue_name: FixedStr::new("ALA"),
            residue_seq: 1,
            chain_id: FixedStr::new("L"),
            is_hetatm: false,
            molecule: pixelfold_core::Molecule::Unknown,
            altloc: None,
            occupancy: 1.0,
            insertion_code,
            model_number: 1,
            position: Vec3::ZERO,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    #[test]
    fn highlighting_a_residue_leaves_its_insertion_code_neighbours_alone() {
        let mut app = App::new();
        app.protein = Some(Protein {
            atoms: vec![
                atom("N", Some('A')),
                atom("CA", Some('A')),
                atom("N", None),
                atom("CA", None),
            ],
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        });
        app.selected_atom_idx = Some(1);

        update_highlighted_atoms(&mut app);

        assert_eq!(app.highlighted_atom_indices, HashSet::from([0, 1]));
    }

    #[test]
    fn stepping_models_wraps_and_starts_the_selection_over() {
        let cif = "\
data_NMR
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA GLY A 1 0.0 0.0 0.0 1
ATOM 2 C CA GLY A 2 3.8 0.0 0.0 1
ATOM 3 C CA GLY A 1 0.5 0.0 0.0 2
ATOM 4 C CA GLY A 2 4.3 0.0 0.0 2
";
        let path = std::env::temp_dir().join("pf_tui_two_models.cif");
        std::fs::write(&path, cif).unwrap();
        let mut app = App::new();
        app.load_protein(&path, 64.0, 64.0, AltlocPolicy::Occupancy, None)
            .unwrap();
        let _ = std::fs::remove_file(&path);
        let shown = |app: &App| app.protein.as_ref().unwrap().atoms[0].model_number;

        app.selected_atom_idx = Some(1);
        app.step_model(1).unwrap();
        assert_eq!(shown(&app), 2);
        assert_eq!(app.selected_atom_idx, None);

        app.step_model(1).unwrap();
        assert_eq!(shown(&app), 1);
        app.step_model(-1).unwrap();
        assert_eq!(shown(&app), 2);
    }
}
