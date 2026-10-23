use std::collections::HashSet;
use std::path::Path;

use anyhow::{Context, Result};
use pixelfold_core::parser::{AssemblyChoice, Load, ParsedStructure};
use pixelfold_core::{Protein, get_calpha_connections};
use pixelfold_render::Framebuffer;
use pixelfold_render::renderer::{self, Camera};
use ratatui::layout::Rect;
use ratatui_image::picker::{Picker, ProtocolType};

mod search;

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
/// Which atoms the viewer draws.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum DisplayMode {
    AllAtoms,
    Backbone,
}

pub(crate) struct RenderedFrame {
    pub(crate) fb: Framebuffer,
    pub(crate) area: Rect,
}

pub(crate) struct App {
    pub(crate) protein: Option<Protein>,
    pub(crate) camera: Camera,
    pub(crate) redraw_needed: bool,
    pub(crate) ca_indices: Vec<usize>,
    pub(crate) backbone_connections: Vec<(usize, usize)>,
    pub(crate) inspect_mode: bool,
    pub(crate) residue_highlight: bool,
    pub(crate) selected_atom_idx: Option<usize>,
    /// Every atom sharing the selected atom's chain, residue number, and
    /// insertion code.
    pub(crate) highlighted_atom_indices: HashSet<usize>,
    pub(crate) display_mode: DisplayMode,
    pub(crate) show_connections: bool,
    pub(crate) use_bfactor_colors: bool,
    pub(crate) show_hydrogen_bonds: bool,
    pub(crate) hbond_energy_threshold: f32,
    /// Multiplies van der Waals radii; 1.0 is full space-filling.
    pub(crate) radius_scale: f32,
    pub(crate) show_depth_cue: bool,
    /// Depth-band clip `(far, near)` in `0..1`, or `None` for the whole depth.
    pub(crate) slab: Option<(f32, f32)>,
    pub(crate) last_frame: Option<RenderedFrame>,
    /// The linked network pane: shown beside the structure, built on first open.
    pub(crate) show_network: bool,
    pub(crate) network_view: Option<network::NetworkView>,
    /// A network build running on a worker thread, delivering its result here.
    /// Set while the pane is computing to prevent blocking of the event loop.
    pub(crate) network_build: Option<std::sync::mpsc::Receiver<network::NetworkView>>,
    /// Where the network pane last drew, for routing clicks to it.
    pub(crate) network_area: Option<Rect>,
    /// The parsed file the shown model came from, kept so stepping between
    /// models never re-reads it.
    pub(crate) source: Option<ParsedStructure>,
    /// The model, assembly, and altloc policy of what is shown.
    pub(crate) load: Load,
    /// The file's model numbers, and the position of the shown one among them.
    pub(crate) models: Vec<usize>,
    pub(crate) model_position: usize,
    /// The ids of the biological assemblies the file declares.
    pub(crate) assemblies: Vec<String>,
}

impl Default for App {
    fn default() -> Self {
        Self::new()
    }
}

impl App {
    pub(crate) fn new() -> Self {
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
            load: Load::default(),
            models: Vec::new(),
            model_position: 0,
            assemblies: Vec::new(),
        }
    }

    /// Load what `load` names from the file at `path` and frame it.
    pub(crate) fn load_protein(
        &mut self,
        path: &Path,
        width: f32,
        height: f32,
        load: Load,
    ) -> Result<()> {
        let source = ParsedStructure::read(path)?;
        let protein = source.protein(&load)?;
        self.models = source.model_numbers();
        self.model_position = load
            .model
            .and_then(|number| self.models.iter().position(|&m| m == number))
            .unwrap_or(0);
        self.assemblies = source.assembly_ids();
        self.source = Some(source);
        self.load = load;

        renderer::auto_frame_protein(&protein, &mut self.camera, width, height);
        self.protein = Some(protein);
        self.compute_static_geometry();

        Ok(())
    }

    /// Show the next model of a multi-model file (`delta` 1) or the previous one
    /// (-1), wrapping at either end. The camera stays put.
    pub(crate) fn step_model(&mut self, delta: isize) -> Result<()> {
        let count = self.models.len();
        if count < 2 {
            return Ok(());
        }

        let position = (self.model_position as isize + delta).rem_euclid(count as isize) as usize;
        self.show(Load {
            model: Some(self.models[position]),
            ..self.load.clone()
        })?;
        self.model_position = position;

        Ok(())
    }

    /// Switch between the deposited coordinates and the primary biological
    /// assembly, framing the result in `frame` (width, height) when given, since
    /// the extent changes. Does nothing for a file that declares no assembly.
    pub(crate) fn toggle_assembly(&mut self, frame: Option<(f32, f32)>) -> Result<()> {
        if self.assemblies.is_empty() {
            return Ok(());
        }

        let assembly = match self.load.assembly {
            Some(_) => None,
            None => Some(AssemblyChoice::Primary),
        };
        self.show(Load {
            assembly,
            ..self.load.clone()
        })?;
        if let (Some((width, height)), Some(protein)) = (frame, self.protein.as_ref()) {
            renderer::auto_frame_protein(protein, &mut self.camera, width, height);
        }

        Ok(())
    }

    /// The assembly shown, by id, or `None` for the deposited coordinates.
    pub(crate) fn shown_assembly(&self) -> Option<&str> {
        match self.load.assembly.as_ref()? {
            AssemblyChoice::Id(id) => Some(id.as_str()),
            AssemblyChoice::Primary => self
                .assemblies
                .iter()
                .find(|id| *id == "1")
                .or_else(|| self.assemblies.first())
                .map(String::as_str),
        }
    }

    /// Build what `load` names from the kept file in place of what is shown.
    /// The selection and the network pane, which index the old atoms, start
    /// over.
    fn show(&mut self, load: Load) -> Result<()> {
        let Some(source) = self.source.as_ref() else {
            return Ok(());
        };

        self.protein = Some(source.protein(&load)?);
        self.load = load;
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
    pub(crate) fn compute_static_geometry(&mut self) {
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
pub fn view(path: &Path, load: Load, protocol: ViewProtocol) -> Result<()> {
    let mut terminal = ratatui::init();

    let picker = build_picker(protocol);

    let result = run_view(&mut terminal, &picker, path, load);

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
    load: Load,
) -> Result<()> {
    crossterm::execute!(std::io::stdout(), crossterm::event::EnableMouseCapture)?;

    let (cols, rows) =
        crossterm::terminal::size().context("failed to query terminal size (not a terminal?)")?;
    let (font_w, font_h) = picker.font_size();
    let width = (cols as u32 * font_w as u32).max(1) as f32;
    let height = (rows as u32 * font_h as u32).max(1) as f32;

    let mut app = App::new();
    app.load_protein(path, width, height, load)?;

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
pub(crate) fn update_highlighted_atoms(app: &mut App) {
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
pub(crate) fn pick_at(app: &App, column: u16, row: u16) -> Option<usize> {
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
        app.load_protein(&path, 64.0, 64.0, Load::default())
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

    #[test]
    fn the_assembly_toggles_on_and_off_and_starts_the_selection_over() {
        let cif = "\
data_ASM
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
_atom_site.auth_seq_id
_atom_site.auth_asym_id
ATOM 1 C CA GLY A 1 1.0 0.0 0.0 1 A
ATOM 2 C CA GLY A 2 4.8 0.0 0.0 2 A
#
_pdbx_struct_oper_list.id 1
_pdbx_struct_oper_list.type 'identity operation'
_pdbx_struct_oper_list.matrix[1][1] 1
_pdbx_struct_oper_list.matrix[1][2] 0
_pdbx_struct_oper_list.matrix[1][3] 0
_pdbx_struct_oper_list.vector[1] 0
_pdbx_struct_oper_list.matrix[2][1] 0
_pdbx_struct_oper_list.matrix[2][2] 1
_pdbx_struct_oper_list.matrix[2][3] 0
_pdbx_struct_oper_list.vector[2] 0
_pdbx_struct_oper_list.matrix[3][1] 0
_pdbx_struct_oper_list.matrix[3][2] 0
_pdbx_struct_oper_list.matrix[3][3] 1
_pdbx_struct_oper_list.vector[3] 0
#
_pdbx_struct_assembly_gen.assembly_id 1
_pdbx_struct_assembly_gen.oper_expression 1,1
_pdbx_struct_assembly_gen.asym_id_list A
#
";
        let path = std::env::temp_dir().join("pf_tui_assembly.cif");
        std::fs::write(&path, cif).unwrap();
        let mut app = App::new();
        app.load_protein(&path, 64.0, 64.0, Load::default())
            .unwrap();
        let _ = std::fs::remove_file(&path);
        let chains = |app: &App| -> Vec<String> {
            let protein = app.protein.as_ref().unwrap();
            protein
                .atoms
                .iter()
                .map(|a| a.chain_id.to_string())
                .collect()
        };
        assert_eq!(app.shown_assembly(), None);

        app.selected_atom_idx = Some(0);
        app.toggle_assembly(Some((64.0, 64.0))).unwrap();
        assert_eq!(chains(&app), ["A", "A", "A-2", "A-2"]);
        assert_eq!(app.shown_assembly(), Some("1"));
        assert_eq!(app.selected_atom_idx, None);

        app.toggle_assembly(None).unwrap();
        assert_eq!(chains(&app), ["A", "A"]);
    }

    #[test]
    fn a_file_without_assemblies_ignores_the_toggle() {
        let mut app = App::new();
        let path = std::env::temp_dir().join("pf_tui_no_assembly.cif");
        std::fs::write(
            &path,
            "data_X\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
             _atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
             _atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
             ATOM 1 C CA GLY A 1 0.0 0.0 0.0\n",
        )
        .unwrap();
        app.load_protein(&path, 64.0, 64.0, Load::default())
            .unwrap();
        let _ = std::fs::remove_file(&path);

        app.toggle_assembly(None).unwrap();

        assert_eq!(app.protein.as_ref().unwrap().atoms.len(), 1);
        assert_eq!(app.shown_assembly(), None);
    }
}
