//! Triaging poses in the viewer: the pose picked from the table is drawn on
//! the reference's pocket, its ligand tinted, the network pane following.

use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use pixelfold_core::Protein;
use pixelfold_core::parser::{Load, ParsedStructure};

use crate::framing::frame_pose;
use crate::pose_table::{Pose, PoseTable};
use crate::{App, DisplayMode};

/// The frame a pose is drawn in: the reference's, or its file's own.
#[derive(PartialEq)]
enum Placement<'a> {
    Reference,
    Own(&'a Path),
}

fn placement(pose: &Pose) -> Placement<'_> {
    match pose.superposition {
        Some(_) => Placement::Reference,
        None => Placement::Own(&pose.path),
    }
}

impl App {
    /// Triage `poses`, given best first: the table lists them and the best is
    /// drawn, framed once the pane is. SD file poses are laid into
    /// `receptor`. `load` names the assembly and altloc policy of every pose's
    /// model.
    pub(crate) fn load_poses(
        &mut self,
        poses: Vec<Pose>,
        receptor: Option<Protein>,
        load: Load,
    ) -> Result<()> {
        anyhow::ensure!(!poses.is_empty(), "there are no poses to show");
        self.poses = Some(PoseTable::new(poses, receptor));
        self.load = load;
        // Space-filling atoms bury a ligand in its pocket; the trace does not.
        self.display_mode = DisplayMode::Backbone;

        self.show_pose(0)
    }

    /// Draw the pose listed `delta` rows from the one drawn, wrapping, past
    /// any that cannot be drawn; the table names the last of those.
    pub(crate) fn step_pose(&mut self, delta: isize) {
        let Some(table) = self.poses.as_ref() else {
            return;
        };
        let (count, mut rank) = (table.poses.len(), table.shown.unwrap_or(0));

        let mut failure = None;
        for _ in 0..count {
            rank = self
                .poses
                .as_ref()
                .map_or(rank, |table| table.stepped_from(rank, delta));
            match self.show_pose(rank) {
                Ok(()) => break,
                Err(error) => failure = Some(self.failure(rank, &error)),
            }
        }
        if let Some(table) = self.poses.as_mut() {
            table.error = failure;
        }
        self.redraw_needed = true;
    }

    /// Draw the pose ranked `rank`, laid on the reference's pocket when it can
    /// be. The camera stays put when the pose is drawn in the frame of the one
    /// it replaces, and is framed at the next drawing otherwise. Nothing
    /// changes when it cannot be drawn.
    fn show_pose(&mut self, rank: usize) -> Result<()> {
        let Some(table) = self.poses.as_ref() else {
            return Ok(());
        };
        if table.shown == Some(rank) {
            return Ok(());
        }
        let pose = &table.poses[rank];
        let replaced = table.shown_pose();
        let reread = replaced.is_none_or(|shown| shown.path != pose.path);
        let reframe = replaced.is_none_or(|shown| placement(shown) != placement(pose));
        let (path, model, superposition) = (pose.path.clone(), pose.model, pose.superposition);
        let beyond = pose.ligand.last().copied();
        let load = Load {
            model: Some(model),
            ..self.load.clone()
        };

        let (mut protein, source) = match (&pose.docked, &table.receptor) {
            (Some(ligand), Some(receptor)) => {
                let mut protein = receptor.clone();
                protein.atoms.extend(ligand.iter().cloned());
                (protein, None)
            }
            (Some(_), None) => anyhow::bail!("no receptor to lay {} into", path.display()),
            (None, _) => {
                let source = match self.source.clone().filter(|_| !reread) {
                    Some(source) => source,
                    None => Arc::new(
                        ParsedStructure::read(&path)
                            .with_context(|| format!("failed to load {}", path.display()))?,
                    ),
                };
                let protein = source
                    .protein(&load)
                    .with_context(|| format!("{}, model {model}", path.display()))?;
                (protein, Some(source))
            }
        };
        if beyond.is_some_and(|atom| atom >= protein.atoms.len()) {
            anyhow::bail!("{} has changed since it was scored", path.display());
        }
        if let Some(superposition) = superposition {
            for atom in &mut protein.atoms {
                atom.position = superposition.apply(atom.position);
            }
        }

        self.models = source
            .as_ref()
            .map_or_else(|| vec![model], |source| source.model_numbers());
        self.model_position = self.models.iter().position(|&m| m == model).unwrap_or(0);
        self.source = source;
        self.load = load;
        if let Some(table) = self.poses.as_mut() {
            table.shown = Some(rank);
            table.reframe |= reframe;
        }
        self.replace_protein(protein, false);

        Ok(())
    }

    /// Why the pose ranked `rank` could not be drawn, naming it.
    fn failure(&self, rank: usize, error: &anyhow::Error) -> String {
        let named = self.poses.as_ref().map(|table| &table.poses[rank]);
        match named {
            Some(pose) => format!(
                "could not draw {} model {}: {error:#}",
                pose.label, pose.model
            ),
            None => format!("{error:#}"),
        }
    }

    /// Frame the pose drawn in `width` by `height` pixels, if drawing it in a
    /// new frame, or `f`, asked for it: the pane's size is known only when it
    /// is drawn.
    pub(crate) fn frame_pose_if_asked(&mut self, width: f32, height: f32) {
        let (Some(table), Some(protein)) = (self.poses.as_mut(), self.protein.as_ref()) else {
            return;
        };
        if std::mem::take(&mut table.reframe) {
            frame_pose(&mut self.camera, protein, table.ligand(), width, height);
        }
    }

    /// Frame the pose drawn afresh at the next drawing.
    pub(crate) fn reframe_pose(&mut self) {
        if let Some(table) = self.poses.as_mut() {
            table.reframe = true;
            self.redraw_needed = true;
        }
    }

    /// Sort the pose table by the next column in turn.
    pub(crate) fn cycle_pose_sort(&mut self) {
        if let Some(table) = self.poses.as_mut() {
            table.cycle_sort();
            self.redraw_needed = true;
        }
    }

    /// Show or hide the pose table.
    pub(crate) fn toggle_pose_table(&mut self) {
        if let Some(table) = self.poses.as_mut() {
            table.visible = !table.visible;
            self.redraw_needed = true;
        }
    }

    /// Draw the pose listed at a terminal cell, as the table was last drawn.
    /// Returns whether the cell falls in the table, which takes the click
    /// wherever it lands there.
    pub(crate) fn pick_pose(&mut self, column: u16, row: u16) -> bool {
        let Some(table) = self
            .poses
            .as_ref()
            .filter(|table| table.covers(column, row))
        else {
            return false;
        };
        if let Some(rank) = table.pose_at(column, row) {
            let failure = self
                .show_pose(rank)
                .err()
                .map(|error| self.failure(rank, &error));
            if let Some(table) = self.poses.as_mut() {
                table.error = failure;
            }
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;
    use crate::pose_table::PoseSort;
    use glam::Vec3;
    use pixelfold_core::rmsd::Superposition;
    use pixelfold_render::renderer::Camera;
    use ratatui::Terminal;
    use ratatui::backend::TestBackend;

    /// A file of one model per shift, each three glycine Cαs and, behind
    /// them as first seen, a one-atom ligand, moved that far along x.
    fn poses_file(name: &str, shifts: &[f32]) -> PathBuf {
        let mut cif = String::from(
            "data_POSES\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
             _atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
             _atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n\
             _atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\n",
        );
        let atoms = [
            ("ATOM", "CA", "GLY", "A", 1, Vec3::ZERO),
            ("ATOM", "CA", "GLY", "A", 2, Vec3::new(3.8, 0.0, 0.0)),
            ("ATOM", "CA", "GLY", "A", 3, Vec3::new(3.8, 3.8, 0.0)),
            ("HETATM", "C1", "LIG", "B", 1, Vec3::new(2.0, 2.0, -2.0)),
        ];
        let mut id = 0;
        for (model, shift) in shifts.iter().enumerate() {
            for (group, name, residue, chain, seq, position) in atoms {
                id += 1;
                let at = position + Vec3::X * *shift;
                cif += &format!(
                    "{group} {id} C {name} {residue} {chain} {seq} {} {} {} {}\n",
                    at.x,
                    at.y,
                    at.z,
                    model + 1
                );
            }
        }
        let path = std::env::temp_dir().join(format!("pf_tui_poses_{name}.cif"));
        std::fs::write(&path, cif).unwrap();
        path
    }

    /// Model `model` of `path` as a pose, its ligand the fourth atom.
    fn pose(path: &Path, model: usize, superposition: Option<Superposition>) -> Pose {
        Pose {
            label: path.file_name().unwrap().to_string_lossy().into_owned(),
            path: path.to_path_buf(),
            model,
            copy: None,
            similarity: 1.0,
            recovery: 1.0,
            rmsd: Some(0.0),
            clashes: 0,
            inverted: Some(Vec::new()),
            lost: Vec::new(),
            ligand: vec![3],
            superposition,
            docked: None,
        }
    }

    /// What lays a model moved `shift` along x back on the first.
    fn shifted_back(shift: f32) -> Superposition {
        let points = [Vec3::ZERO, Vec3::X, Vec3::Y, Vec3::Z];
        let moved: Vec<Vec3> = points.iter().map(|&p| p + Vec3::X * shift).collect();
        Superposition::of(&moved, &points)
    }

    /// An app triaging `poses`, drawn once in 64 by 64 pixels.
    fn app_of(poses: Vec<Pose>) -> App {
        let mut app = App::new();
        app.load_poses(poses, None, Load::default()).unwrap();
        app.frame_pose_if_asked(64.0, 64.0);
        app
    }

    fn first_atom(app: &App) -> Vec3 {
        app.protein.as_ref().unwrap().atoms[0].position
    }

    fn shown(app: &App) -> usize {
        app.poses.as_ref().unwrap().shown.unwrap()
    }

    #[test]
    fn each_pose_is_drawn_on_the_reference_pocket_its_ligand_tinted() {
        let path = poses_file("pocket", &[0.0, 10.0]);
        let mut app = app_of(vec![
            pose(&path, 1, Some(shifted_back(0.0))),
            pose(&path, 2, Some(shifted_back(10.0))),
        ]);
        assert_eq!(shown(&app), 0, "the best pose is shown first");
        assert_eq!(
            app.display_mode,
            DisplayMode::Backbone,
            "poses open on the trace, which hides no ligand"
        );

        app.step_pose(1);
        let _ = std::fs::remove_file(&path);

        assert_eq!(shown(&app), 1);
        assert_eq!(app.protein.as_ref().unwrap().atoms[0].model_number, 2);
        assert!(first_atom(&app).length() < 1e-4, "{}", first_atom(&app));
        let fb = crate::render::render_scene(&app, 64, 64);
        let ligand: Vec<[u8; 4]> = (0..fb.height())
            .flat_map(|y| (0..fb.width()).map(move |x| (x, y)))
            .filter(|&(x, y)| fb.id_at(x, y) == 3)
            .map(|(x, y)| fb.pixel(x, y))
            .collect();
        assert!(!ligand.is_empty(), "the ligand drew no pixels");
        // Magenta survives shading as no green and red equal to blue.
        assert!(ligand.iter().all(|&[r, g, b, _]| g == 0 && r == b && r > 0));
    }

    #[test]
    fn the_poses_step_in_the_order_the_table_is_sorted_in() {
        let path = poses_file("order", &[0.0, 0.0, 0.0]);
        let measured = |model, rmsd| Pose {
            rmsd,
            ..pose(&path, model, None)
        };
        let mut app = app_of(vec![
            measured(1, Some(2.0)),
            measured(2, None),
            measured(3, Some(0.5)),
        ]);

        app.cycle_pose_sort();
        let table = app.poses.as_ref().unwrap();
        assert_eq!(table.sort, PoseSort::Rmsd);
        assert_eq!(table.order, [2, 0, 1], "the unmeasured last");
        assert_eq!(shown(&app), 0, "a sort keeps the pose shown");

        app.step_pose(1);
        assert_eq!(shown(&app), 1);
        app.step_pose(1);
        assert_eq!(shown(&app), 2, "the last wraps to the first");
        app.step_pose(-1);
        assert_eq!(shown(&app), 1);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn the_camera_stays_put_between_poses_on_the_reference_pocket() {
        let path = poses_file("camera", &[0.0, 10.0, 30.0]);
        let mut app = app_of(vec![
            pose(&path, 1, Some(shifted_back(0.0))),
            pose(&path, 2, Some(shifted_back(10.0))),
            pose(&path, 3, None),
        ]);
        app.camera.zoom = 7.0;

        app.step_pose(1);
        app.frame_pose_if_asked(64.0, 64.0);
        assert_eq!(app.camera.zoom, 7.0, "a pose on the pocket keeps the view");

        app.step_pose(1);
        app.frame_pose_if_asked(64.0, 64.0);
        let _ = std::fs::remove_file(&path);
        assert_ne!(app.camera.zoom, 7.0, "a pose in its own frame is framed");
        assert!(
            (first_atom(&app).x - 30.0).abs() < 1e-4,
            "where its file puts it"
        );
    }

    #[test]
    fn a_pose_framed_afresh_turns_its_ligand_to_the_viewer() {
        let path = poses_file("facing", &[0.0]);
        let app = app_of(vec![pose(&path, 1, None)]);
        let _ = std::fs::remove_file(&path);

        let depth = |atom: &pixelfold_core::Atom| {
            app.camera
                .project_point_with_depth(atom.position, 64.0, 64.0)
                .2
        };
        let atoms = &app.protein.as_ref().unwrap().atoms;
        let ligand = app.camera.orientation * atoms[3].position;
        let blocked = atoms[..3].iter().any(|atom| {
            let at = app.camera.orientation * atom.position;
            at.z > ligand.z && (at - ligand).truncate().length() < 2.0 * 1.7
        });
        assert!(!blocked, "an atom stands between the ligand and the viewer");
        assert!(depth(&atoms[3]).is_finite());
    }

    #[test]
    fn the_first_pose_waits_for_the_pane_to_be_framed_and_f_frames_it_again() {
        let path = poses_file("first-frame", &[0.0]);
        let mut app = App::new();
        app.load_poses(vec![pose(&path, 1, None)], None, Load::default())
            .unwrap();
        let _ = std::fs::remove_file(&path);
        assert_eq!(app.camera.zoom, Camera::new().zoom, "no pane to fit yet");

        app.frame_pose_if_asked(64.0, 32.0);
        let framed = (app.camera.zoom, app.camera.orientation);
        app.frame_pose_if_asked(10.0, 10.0);
        assert_eq!((app.camera.zoom, app.camera.orientation), framed, "once");

        app.camera.rotate(0.4, 0.9, 0.0);
        crate::inputs::handle_input(
            &mut app,
            crossterm::event::KeyCode::Char('f'),
            crossterm::event::KeyModifiers::NONE,
        )
        .unwrap();
        app.frame_pose_if_asked(64.0, 32.0);
        assert_eq!(
            (app.camera.zoom, app.camera.orientation),
            framed,
            "f looks at the ligand again"
        );
    }

    #[test]
    fn a_click_after_a_pose_change_picks_nothing_until_it_is_drawn() {
        let first = poses_file("stale-first", &[0.0]);
        let second = poses_file("stale-second", &[0.0]);
        let mut app = app_of(vec![pose(&first, 1, None), pose(&second, 1, None)]);
        let area = ratatui::layout::Rect::new(0, 0, 64, 32);
        app.last_frame = Some(crate::RenderedFrame {
            fb: crate::render::render_scene(&app, 64, 64),
            area,
            stale: false,
        });
        let cell = (0..area.height)
            .flat_map(|row| (0..area.width).map(move |column| (column, row)))
            .find(|&(column, row)| crate::pick_at(&app, column, row).is_some())
            .expect("an atom under some cell");

        app.step_pose(1);
        let _ = std::fs::remove_file(&first);
        let _ = std::fs::remove_file(&second);

        assert_eq!(crate::pick_at(&app, cell.0, cell.1), None);
    }

    #[test]
    fn a_pose_that_cannot_be_drawn_is_stepped_past_and_named() {
        let path = poses_file("broken", &[0.0, 10.0]);
        let missing = std::env::temp_dir().join("pf_tui_poses_missing.cif");
        let changed = Pose {
            ligand: vec![99],
            ..pose(&path, 1, None)
        };
        let mut app = app_of(vec![
            pose(&path, 1, None),
            pose(&missing, 1, None),
            changed,
            pose(&path, 2, None),
        ]);

        app.step_pose(1);
        let table = app.poses.as_ref().unwrap();
        assert_eq!(table.shown, Some(3), "past both");
        let error = table.error.clone().unwrap_or_default();
        assert!(error.contains("changed since"), "{error}");

        app.step_pose(-1);
        assert_eq!(shown(&app), 0, "back past both");
        let error = app
            .poses
            .as_ref()
            .unwrap()
            .error
            .clone()
            .unwrap_or_default();
        assert!(error.contains("pf_tui_poses_missing.cif"), "{error}");

        app.step_pose(0);
        let _ = std::fs::remove_file(&path);
        assert_eq!(
            app.poses.as_ref().unwrap().error,
            None,
            "drawn, said nothing"
        );
    }

    #[test]
    fn the_pose_drawn_is_not_drawn_again() {
        let path = poses_file("again", &[0.0]);
        let mut app = app_of(vec![pose(&path, 1, None)]);
        app.selected_atom_idx = Some(1);

        app.step_pose(1);
        let _ = std::fs::remove_file(&path);

        assert_eq!(app.selected_atom_idx, Some(1), "the selection stays");
    }

    #[test]
    fn the_file_keys_leave_a_pose_alone() {
        let path = poses_file("file-keys", &[0.0, 10.0]);
        let mut app = app_of(vec![pose(&path, 1, Some(shifted_back(0.0)))]);
        app.assemblies = vec!["1".to_string()];

        app.step_model(1).unwrap();
        app.toggle_assembly(None).unwrap();
        let _ = std::fs::remove_file(&path);

        assert_eq!(app.protein.as_ref().unwrap().atoms[0].model_number, 1);
        assert_eq!(app.load.assembly, None);
    }

    #[test]
    fn a_pose_in_another_file_is_read_from_that_file() {
        let first = poses_file("first", &[0.0]);
        let second = poses_file("second", &[50.0, 60.0]);
        let mut app = app_of(vec![pose(&first, 1, None), pose(&second, 2, None)]);

        app.step_pose(1);
        assert!((first_atom(&app).x - 60.0).abs() < 1e-4);
        assert_eq!(app.models, [1, 2]);

        app.step_pose(1);
        let _ = std::fs::remove_file(&first);
        let _ = std::fs::remove_file(&second);
        assert!(first_atom(&app).length() < 1e-4, "and back");
        assert_eq!(app.models, [1]);
    }

    #[test]
    fn the_network_pane_follows_the_pose_and_spans_no_ensemble() {
        let path = poses_file("network", &[0.0, 10.0]);
        let mut app = app_of(vec![pose(&path, 1, None), pose(&path, 2, None)]);
        app.show_network = true;

        app.step_pose(1);
        let _ = std::fs::remove_file(&path);
        assert!(
            app.network_build.is_some(),
            "the pose shown gets its network"
        );

        app.cycle_network_scope();
        assert!(!app.network_scope.is_ensemble(), "poses are no ensemble");
    }

    #[test]
    fn the_table_lists_the_poses_and_a_click_on_a_row_shows_that_pose() {
        let path = poses_file("click", &[0.0, 0.0, 0.0]);
        let mut app = app_of(vec![
            pose(&path, 1, None),
            Pose {
                rmsd: Some(1.234),
                lost: vec!["HIS A/57 hydrogen-bond".to_string()],
                ..pose(&path, 2, None)
            },
            pose(&path, 3, None),
        ]);
        let mut terminal = Terminal::new(TestBackend::new(100, 8)).unwrap();
        let mut draw = |app: &mut App| {
            terminal
                .draw(|frame| app.poses.as_mut().unwrap().render(frame, frame.area()))
                .unwrap();
            let buffer = terminal.backend().buffer();
            buffer
                .content
                .iter()
                .map(|cell| cell.symbol())
                .collect::<String>()
        };

        let text = draw(&mut app);
        assert!(text.contains("sorted by similarity"), "{text}");
        // Below the top border and the header, the second row.
        assert!(app.pick_pose(10, 3));
        assert_eq!(shown(&app), 1);
        assert!(!app.pick_pose(10, 20), "a click elsewhere");

        let text = draw(&mut app);
        assert!(text.contains("1.23"), "{text}");
        assert!(text.contains("HIS A/57 hydrogen-bond"), "{text}");

        // A sort not drawn yet leaves a click on the rows as they were drawn.
        app.cycle_pose_sort();
        assert!(app.pick_pose(10, 4));
        let _ = std::fs::remove_file(&path);
        assert_eq!(shown(&app), 2, "the third row drawn");
    }

    #[test]
    fn a_docked_pose_is_drawn_in_the_receptor_and_read_from_no_file() {
        // The receptor is the poses file's model without its ligand; each
        // docked ligand is a copy of that ligand moved along x, and the SD
        // file they came from no longer exists.
        let path = poses_file("receptor", &[0.0]);
        let model = ParsedStructure::read(&path)
            .unwrap()
            .protein(&Load::default())
            .unwrap();
        let _ = std::fs::remove_file(&path);
        let receptor = model.without(|index, _| index == 3);
        let gone = std::env::temp_dir().join("pf_tui_poses_docked_gone.sdf");
        let docked = |number: usize, shift: f32| {
            let mut ligand = model.atoms[3].clone();
            ligand.position.x += shift;
            Pose {
                docked: Some(vec![ligand]),
                ..pose(&gone, number, Some(shifted_back(0.0)))
            }
        };

        let mut app = App::new();
        app.load_poses(
            vec![docked(1, 0.0), docked(2, 5.0)],
            Some(receptor),
            Load::default(),
        )
        .unwrap();
        app.step_pose(1);

        let table = app.poses.as_ref().unwrap();
        assert_eq!(table.error, None);
        assert_eq!(shown(&app), 1);
        let atoms = &app.protein.as_ref().unwrap().atoms;
        assert_eq!(atoms.len(), 4);
        assert_eq!(atoms[3].position, Vec3::new(7.0, 2.0, -2.0));
    }
}
