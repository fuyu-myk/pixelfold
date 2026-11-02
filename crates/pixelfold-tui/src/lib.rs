use std::collections::HashSet;
use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use pixelfold_core::parser::{AssemblyChoice, Load, ParsedStructure};
use pixelfold_core::{Protein, get_calpha_connections};
use pixelfold_render::Framebuffer;
use pixelfold_render::renderer::{self, Camera};
use ratatui::layout::Rect;
use ratatui_image::picker::{Picker, ProtocolType};

mod search;

mod app;
mod framing;
mod inputs;
mod iterm2;
mod network;
mod pose_table;
mod poses;
mod render;
mod ui;

pub use pose_table::Pose;

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

/// Why a network build that delivered nothing failed: its worker stopped.
pub(crate) const BUILD_STOPPED: &str = "the network build stopped";

/// The result a network build's worker delivers: the view, or why none.
pub(crate) type NetworkBuild = std::sync::mpsc::Receiver<Result<network::NetworkView, String>>;

pub(crate) struct RenderedFrame {
    pub(crate) fb: Framebuffer,
    pub(crate) area: Rect,
    /// Set once other atoms replace those drawn, whose ids the frame holds.
    pub(crate) stale: bool,
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
    pub(crate) network_build: Option<NetworkBuild>,
    /// Whether the build in flight only lays a network already built over the
    /// shown model, rather than detecting its interactions.
    pub(crate) network_remap: bool,
    /// Why the last network build failed, shown in the pane in its place.
    pub(crate) network_error: Option<String>,
    /// Where the network pane last drew, for routing clicks to it.
    pub(crate) network_area: Option<Rect>,
    /// Whose network the pane shows, and the ensemble's frequency floor.
    pub(crate) network_scope: network::NetworkScope,
    /// The parsed file the shown model came from, kept so stepping between
    /// models never re-reads it, and shared with an ensemble network's build.
    pub(crate) source: Option<Arc<ParsedStructure>>,
    /// What reading the shown file skipped or cut short, for the header: the
    /// viewer owns the terminal, so the parser prints nothing.
    pub(crate) warnings: Vec<String>,
    /// The model, assembly, and altloc policy of what is shown.
    pub(crate) load: Load,
    /// The file's model numbers, and the position of the shown one among them.
    pub(crate) models: Vec<usize>,
    pub(crate) model_position: usize,
    /// The ids of the biological assemblies the file declares.
    pub(crate) assemblies: Vec<String>,
    /// The poses being triaged, when the viewer shows poses, not a file.
    pub(crate) poses: Option<pose_table::PoseTable>,
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
            network_remap: false,
            network_error: None,
            network_area: None,
            network_scope: network::NetworkScope::default(),
            source: None,
            warnings: Vec::new(),
            load: Load::default(),
            models: Vec::new(),
            model_position: 0,
            assemblies: Vec::new(),
            poses: None,
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
        self.warnings = source.warnings().to_vec();
        self.source = Some(Arc::new(source));
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
        // A pose's model is the pose table's to change.
        if count < 2 || self.poses.is_some() {
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
        // A pose's atoms are those it was scored on.
        if self.assemblies.is_empty() || self.poses.is_some() {
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
    /// An ensemble's network, unchanged by which model is shown, is only laid
    /// over the new model's atoms.
    fn show(&mut self, load: Load) -> Result<()> {
        let Some(source) = self.source.as_ref() else {
            return Ok(());
        };

        let same_ensemble = Load {
            model: None,
            ..load.clone()
        } == Load {
            model: None,
            ..self.load.clone()
        };
        let protein = source.protein(&load)?;
        self.load = load;
        self.replace_protein(protein, self.network_scope.is_ensemble() && same_ensemble);

        Ok(())
    }

    /// Show `protein` in place of what is shown. The selection and the network
    /// pane, which index the old atoms, start over; when `relay`, the network
    /// is only laid over the new atoms.
    fn replace_protein(&mut self, protein: Protein, relay: bool) {
        self.protein = Some(protein);
        self.compute_static_geometry();
        if let Some(frame) = self.last_frame.as_mut() {
            frame.stale = true;
        }

        self.selected_atom_idx = None;
        self.residue_highlight = false;
        self.highlighted_atom_indices.clear();
        match (relay, self.network_view.take(), self.network_build.take()) {
            (true, Some(view), _) => {
                self.network_remap = true;
                self.start_network_remap(move || Ok(view));
            }
            // A build or remap still in flight is waited for, not repeated.
            (true, None, Some(pending)) => self.start_network_remap(move || {
                pending
                    .recv()
                    .map_err(|_| BUILD_STOPPED.to_string())
                    .and_then(|built| built)
            }),
            _ if self.show_network => self.start_network_build(),
            _ => {}
        }
        self.redraw_needed = true;
    }

    /// Step the network pane through the shown model's network, the
    /// ensemble's, and the ensemble's with ever fewer rare edges, then back,
    /// opening the pane; a closed pane showing an ensemble only opens. Does
    /// nothing for a file of one model, or for poses, which are no ensemble.
    pub(crate) fn cycle_network_scope(&mut self) {
        if self.models.len() < 2 || self.poses.is_some() {
            return;
        }

        self.redraw_needed = true;
        if !self.show_network {
            self.show_network = true;
            if self.network_scope.is_ensemble() {
                if self.network_view.is_none() && self.network_build.is_none() {
                    self.start_network_build();
                }
                return;
            }
        }

        let next = self.network_scope.next();
        let rebuild = next.is_ensemble() != self.network_scope.is_ensemble();
        self.network_scope = next;
        if rebuild {
            self.network_view = None;
            self.start_network_build();
        } else if let Some(view) = self.network_view.as_mut() {
            view.set_min_frequency(next.min_frequency());
        }
    }

    /// Show a finished network build: the view, under the floor the scope has
    /// now, or why it failed.
    pub(crate) fn install_network(&mut self, built: Result<network::NetworkView, String>) {
        match built {
            Ok(mut view) => {
                view.set_min_frequency(self.network_scope.min_frequency());
                self.network_view = Some(view);
            }
            Err(error) => self.network_error = Some(error),
        }
        self.network_build = None;
        self.network_remap = false;
        self.redraw_needed = true;
    }

    /// Build the network pane's view on a worker thread, so a large structure's
    /// detection, layout, and accessibility pass never freeze the event loop;
    /// the pane shows a placeholder until the result arrives (see `run_app`).
    pub(crate) fn start_network_build(&mut self) {
        let Some(protein) = self.protein.clone() else {
            return;
        };
        let ensemble = self
            .source
            .clone()
            .filter(|_| self.network_scope.is_ensemble())
            .map(|source| (source, self.load.clone(), self.models.clone()));

        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let built = match ensemble {
                Some((source, load, models)) => {
                    network::NetworkView::build_ensemble(&protein, &source, &load, &models)
                        .map_err(|error| format!("{error:#}"))
                }
                None => Ok(network::NetworkView::build(&protein)),
            };
            tx.send(built).ok();
        });
        self.network_build = Some(rx);
        self.network_remap = false;
        self.network_error = None;
    }

    /// Lay the view `built` gives over the shown model on a worker thread, as
    /// [`App::start_network_build`] builds one.
    fn start_network_remap(
        &mut self,
        built: impl FnOnce() -> Result<network::NetworkView, String> + Send + 'static,
    ) {
        let Some(protein) = self.protein.clone() else {
            return;
        };

        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            tx.send(built().map(|view| view.remapped(&protein))).ok();
        });
        self.network_build = Some(rx);
        self.network_error = None;
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
    run_viewer(protocol, |app, width, height| {
        app.load_protein(path, width, height, load)
    })
}

/// Run the interactive viewer on ligand poses ranked best first, listed in a
/// table under the structure; the pose picked is drawn on the reference's
/// pocket. SD file poses are laid into `receptor`. `load` names the assembly
/// and altloc policy of every pose's model.
pub fn view_poses(
    poses: Vec<Pose>,
    receptor: Option<Protein>,
    load: Load,
    protocol: ViewProtocol,
) -> Result<()> {
    run_viewer(protocol, |app, _, _| app.load_poses(poses, receptor, load))
}

/// Run the viewer on what `open` shows, framed in the terminal's pixels.
fn run_viewer(
    protocol: ViewProtocol,
    open: impl FnOnce(&mut App, f32, f32) -> Result<()>,
) -> Result<()> {
    require_terminal(
        "the viewer",
        "for a report, run interactions, rin, poses, ss or sasa, or render -o FILE",
    )?;
    let mut terminal = ratatui::init();

    let picker = build_picker(protocol);

    let result = run_view(&mut terminal, &picker, open);

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
    open: impl FnOnce(&mut App, f32, f32) -> Result<()>,
) -> Result<()> {
    crossterm::execute!(std::io::stdout(), crossterm::event::EnableMouseCapture)?;

    let (cols, rows) =
        crossterm::terminal::size().context("failed to query terminal size (not a terminal?)")?;
    let (font_w, font_h) = picker.font_size();
    let width = (cols as u32 * font_w as u32).max(1) as f32;
    let height = (rows as u32 * font_h as u32).max(1) as f32;

    let mut app = App::new();
    open(&mut app, width, height)?;

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
    require_terminal(
        "the search interface",
        "to download an entry, give its 4-character PDB id to any command",
    )?;
    let mut terminal = ratatui::init();
    crossterm::execute!(std::io::stdout(), crossterm::event::EnableMouseCapture)?;

    let result = f(&mut terminal);
    let _ = crossterm::execute!(std::io::stdout(), crossterm::event::DisableMouseCapture);

    ratatui::restore();
    result
}

/// Fail unless standard input and output are a terminal, which `what` draws on
/// and reads keys from: ratatui panics without one, as in an agent's shell or
/// with output piped. `instead` says what to run there.
fn require_terminal(what: &str, instead: &str) -> Result<()> {
    use std::io::IsTerminal;

    if std::io::stdin().is_terminal() && std::io::stdout().is_terminal() {
        return Ok(());
    }

    anyhow::bail!("{what} needs a terminal, and standard input or output is not one; {instead}")
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
    let frame = app.last_frame.as_ref().filter(|frame| !frame.stale)?;
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
        id => Some(id as usize).filter(|&atom| {
            app.protein
                .as_ref()
                .is_some_and(|protein| atom < protein.atoms.len())
        }),
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

    /// Install the network build in flight, once it finishes.
    fn finish_network_build(app: &mut App) -> &network::NetworkView {
        let built = app
            .network_build
            .take()
            .expect("a network build in flight")
            .recv()
            .expect("the build finishes");
        app.install_network(built);
        app.network_view.as_ref().expect("the build succeeds")
    }

    /// An app showing a file of three models, Cys1 and Cys2 bridged in the
    /// first two.
    fn three_models(name: &str) -> App {
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
ATOM 1 S SG CYS A 1 0.0 0.0 0.0 1
ATOM 2 S SG CYS A 2 2.05 0.0 0.0 1
ATOM 3 S SG CYS A 1 0.0 0.0 0.0 2
ATOM 4 S SG CYS A 2 2.05 0.0 0.0 2
ATOM 5 S SG CYS A 1 0.0 0.0 0.0 3
ATOM 6 S SG CYS A 2 6.0 0.0 0.0 3
";
        let path = std::env::temp_dir().join(format!("pf_tui_{name}.cif"));
        std::fs::write(&path, cif).unwrap();
        let mut app = App::new();
        app.load_protein(&path, 64.0, 64.0, Load::default())
            .unwrap();
        let _ = std::fs::remove_file(&path);
        app
    }

    /// The network's model count and its one edge's frequency.
    fn bridged(view: &network::NetworkView) -> (usize, f32) {
        let network = view.network();
        (network.models, network.frequency(&network.edges[0]))
    }

    fn shown(app: &App) -> usize {
        app.protein.as_ref().unwrap().atoms[0].model_number
    }

    #[test]
    fn the_ensemble_network_spans_every_model_and_outlives_a_model_step() {
        let mut app = three_models("ensemble");

        app.cycle_network_scope();
        assert!(app.show_network, "the scope opens the pane");
        assert_eq!(bridged(finish_network_build(&mut app)), (3, 2.0 / 3.0));

        // The next model keeps the network and its floor, laid over its atoms.
        app.cycle_network_scope();
        app.step_model(1).unwrap();
        assert!(app.network_remap);
        let view = finish_network_build(&mut app);
        assert_eq!(bridged(view), (3, 2.0 / 3.0));
        assert_eq!(view.min_frequency(), 0.5);

        // The last floor, then back to the shown model's own network.
        app.cycle_network_scope();
        assert!(app.network_build.is_none(), "a floor needs no build");
        app.cycle_network_scope();
        assert_eq!(bridged(finish_network_build(&mut app)), (1, 1.0));

        let mut lone = App::new();
        lone.cycle_network_scope();
        assert!(!lone.show_network, "a single model has no ensemble");
    }

    // Keys pressed before a build lands: a floor is not lost and a model step
    // waits for the build in flight rather than starting the ensemble again.
    #[test]
    fn keys_pressed_while_the_ensemble_builds_are_not_lost() {
        let mut app = three_models("racing");

        app.cycle_network_scope();
        app.cycle_network_scope();
        app.step_model(1).unwrap();
        app.step_model(1).unwrap();
        assert!(!app.network_remap, "still waiting on the first build");
        let view = finish_network_build(&mut app);

        assert_eq!(bridged(view), (3, 2.0 / 3.0));
        assert_eq!(view.min_frequency(), 0.5);
        assert_eq!(shown(&app), 3);
        let atoms = app.protein.as_ref().unwrap().atoms.len();
        let view = app.network_view.as_ref().unwrap();
        assert!((0..2).all(|node| view.atom_of_node(node).is_some_and(|atom| atom < atoms)));
    }

    #[test]
    fn reopening_the_pane_on_an_ensemble_keeps_its_floor() {
        let mut app = three_models("reopen");
        app.cycle_network_scope();
        finish_network_build(&mut app);

        app.show_network = false;
        app.cycle_network_scope();

        assert!(app.show_network);
        assert_eq!(app.network_view.as_ref().unwrap().min_frequency(), 0.0);
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
