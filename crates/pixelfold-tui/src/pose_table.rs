//! The pose table: ligand poses scored against a reference, listed best
//! first or sorted by another score, and the rows a click lands on.

use std::path::{Path, PathBuf};

use pixelfold_core::rmsd::Superposition;
use pixelfold_core::{Atom, Protein};
use ratatui::Frame;
use ratatui::layout::{Constraint, Position, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::Line;
use ratatui::widgets::{Block, Borders, Row, Table, TableState};

/// The most rows the table shows before it scrolls.
const MAX_ROWS: usize = 8;

/// One pose of a ligand scored against a reference, as `pixelfold poses`
/// ranks it, with what drawing it needs.
#[derive(Clone, Debug)]
pub struct Pose {
    /// The pose file as it was named.
    pub label: String,
    pub path: PathBuf,
    pub model: usize,
    /// The copy of the ligand scored, when the pose holds several.
    pub copy: Option<String>,
    pub similarity: f32,
    pub recovery: f32,
    /// The ligand's RMSD in angstroms on the reference's pocket, when measured.
    pub rmsd: Option<f32>,
    pub clashes: usize,
    /// The reference stereocentres the pose inverts, or `None` when its atoms
    /// could not be paired with the reference's.
    pub inverted: Option<Vec<String>>,
    /// The checks of its own geometry the ligand fails, each with how many
    /// fail it, or `None` when its atoms could not be paired.
    pub geometry: Option<Vec<String>>,
    /// The reference contacts the pose lacks.
    pub lost: Vec<String>,
    /// The scored copy's atoms, as indices into the model loaded.
    pub ligand: Vec<usize>,
    /// What lays the pose on the reference's pocket, or `None` when too few
    /// pocket atoms have a counterpart.
    pub superposition: Option<Superposition>,
    /// An SD file pose's ligand, laid into the receptor the viewer is given;
    /// `None` for a pose drawn from its file.
    pub docked: Option<Vec<Atom>>,
}

/// What the pose table is sorted by, best first; ties keep the rank.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum PoseSort {
    #[default]
    Rank,
    Rmsd,
    Clashes,
    Inverted,
}

impl PoseSort {
    fn next(self) -> Self {
        match self {
            PoseSort::Rank => PoseSort::Rmsd,
            PoseSort::Rmsd => PoseSort::Clashes,
            PoseSort::Clashes => PoseSort::Inverted,
            PoseSort::Inverted => PoseSort::Rank,
        }
    }

    /// The column the sort orders by.
    fn column(self) -> &'static str {
        match self {
            PoseSort::Rank => "similarity",
            PoseSort::Rmsd => "RMSD",
            PoseSort::Clashes => "clashes",
            PoseSort::Inverted => "inverted",
        }
    }

    /// The `poses`, by rank, in the order this sort lists them.
    fn order(self, poses: &[Pose]) -> Vec<usize> {
        let mut order: Vec<usize> = (0..poses.len()).collect();
        match self {
            PoseSort::Rank => {}
            PoseSort::Rmsd => order.sort_by(|&a, &b| {
                let rmsd = |pose: usize| poses[pose].rmsd.unwrap_or(f32::INFINITY);
                rmsd(a).total_cmp(&rmsd(b))
            }),
            PoseSort::Clashes => order.sort_by_key(|&pose| poses[pose].clashes),
            PoseSort::Inverted => order
                .sort_by_key(|&pose| poses[pose].inverted.as_ref().map_or(usize::MAX, Vec::len)),
        }

        order
    }
}

/// `label` cut to `width` columns from the left, so a path keeps its file.
fn elided(label: &str, width: usize) -> String {
    let count = label.chars().count();
    if count <= width {
        return label.to_string();
    }

    let tail: String = label.chars().skip(count + 1 - width.max(1)).collect();
    format!("\u{2026}{tail}")
}

/// The poses, in rank order, as the table lists them under its sort.
pub(crate) struct PoseTable {
    pub(crate) poses: Vec<Pose>,
    /// What SD file poses are laid into.
    pub(crate) receptor: Option<Protein>,
    pub(crate) order: Vec<usize>,
    pub(crate) sort: PoseSort,
    /// The pose drawn, by rank.
    pub(crate) shown: Option<usize>,
    pub(crate) visible: bool,
    /// Whether the pose drawn waits to be framed, which needs the pane's size.
    pub(crate) reframe: bool,
    /// Why the last pose asked for could not be drawn.
    pub(crate) error: Option<String>,
    state: TableState,
    /// Where the table last drew, and the poses its rows listed, for routing
    /// clicks to it.
    area: Option<Rect>,
    drawn: Vec<usize>,
}

impl PoseTable {
    /// The table of `poses`, given best first, listed by rank, SD file poses
    /// laid into `receptor`.
    pub(crate) fn new(mut poses: Vec<Pose>, receptor: Option<Protein>) -> Self {
        for pose in &mut poses {
            pose.ligand.sort_unstable();
        }

        Self {
            order: PoseSort::Rank.order(&poses),
            poses,
            receptor,
            sort: PoseSort::Rank,
            shown: None,
            visible: true,
            reframe: false,
            error: None,
            state: TableState::default(),
            area: None,
            drawn: Vec::new(),
        }
    }

    pub(crate) fn shown_pose(&self) -> Option<&Pose> {
        self.shown.map(|rank| &self.poses[rank])
    }

    /// The atoms of the ligand of the pose drawn, in order.
    pub(crate) fn ligand(&self) -> &[usize] {
        self.shown_pose().map_or(&[], |pose| &pose.ligand)
    }

    /// Whether `atom` of the model drawn belongs to the pose's ligand.
    pub(crate) fn is_ligand(&self, atom: usize) -> bool {
        self.shown_pose()
            .is_some_and(|pose| pose.ligand.binary_search(&atom).is_ok())
    }

    /// The pose drawn, named with its rank among them all.
    pub(crate) fn caption(&self) -> Option<String> {
        let rank = self.shown?;
        let pose = &self.poses[rank];
        let name = Path::new(&pose.label)
            .file_name()
            .map_or(pose.label.as_str(), |name| {
                name.to_str().unwrap_or(&pose.label)
            });
        Some(format!(
            "{} model {} (#{} of {})",
            name,
            pose.model,
            rank + 1,
            self.poses.len()
        ))
    }

    /// Sort the table by the next column in turn, keeping the pose drawn.
    pub(crate) fn cycle_sort(&mut self) {
        self.sort = self.sort.next();
        self.order = self.sort.order(&self.poses);
    }

    /// Where the pose drawn is listed.
    fn listed_at(&self) -> usize {
        self.shown
            .and_then(|shown| self.order.iter().position(|&rank| rank == shown))
            .unwrap_or(0)
    }

    /// The pose listed `delta` rows from the pose ranked `rank`, wrapping.
    pub(crate) fn stepped_from(&self, rank: usize, delta: isize) -> usize {
        let at = self
            .order
            .iter()
            .position(|&listed| listed == rank)
            .unwrap_or(0);
        let count = self.order.len() as isize;
        self.order[(at as isize + delta).rem_euclid(count) as usize]
    }

    /// The lines the table takes: its rows, up to a few, a header and borders.
    pub(crate) fn height(&self) -> u16 {
        (self.poses.len().min(MAX_ROWS) + 3) as u16
    }

    /// Draw the table in `area`, the pose drawn marked, and what it lacks along
    /// the bottom border.
    pub(crate) fn render(&mut self, frame: &mut Frame, area: Rect) {
        let copies = self.poses.iter().any(|pose| pose.copy.is_some());
        let mut columns = vec![
            ("#", Constraint::Length(4)),
            ("pose", Constraint::Fill(1)),
            ("model", Constraint::Length(5)),
        ];
        if copies {
            columns.push(("copy", Constraint::Length(9)));
        }
        columns.extend([
            ("similarity", Constraint::Length(11)),
            ("recovery", Constraint::Length(9)),
            ("RMSD", Constraint::Length(5)),
            ("clashes", Constraint::Length(8)),
            ("inverted", Constraint::Length(9)),
        ]);
        // The pose column fills what the others and the gaps between them leave.
        let fixed: u16 = columns
            .iter()
            .map(|&(_, width)| match width {
                Constraint::Length(length) => length + 1,
                _ => 1,
            })
            .sum();
        let label_width = usize::from(area.width.saturating_sub(fixed + 1));
        let header = Row::new(columns.iter().map(|&(name, _)| {
            if name == self.sort.column() {
                format!("{name}\u{25be}")
            } else {
                name.to_string()
            }
        }))
        .style(
            Style::default()
                .fg(Color::Gray)
                .add_modifier(Modifier::BOLD),
        );

        let rows = self.order.iter().map(|&rank| {
            let pose = &self.poses[rank];
            let mut cells = vec![
                (rank + 1).to_string(),
                elided(&pose.label, label_width),
                pose.model.to_string(),
            ];
            if copies {
                cells.push(pose.copy.clone().unwrap_or_default());
            }
            cells.extend([
                format!("{:.3}", pose.similarity),
                format!("{:.3}", pose.recovery),
                pose.rmsd
                    .map_or("-".to_string(), |rmsd| format!("{rmsd:.2}")),
                pose.clashes.to_string(),
                pose.inverted
                    .as_ref()
                    .map_or("-".to_string(), |centres| centres.len().to_string()),
            ]);
            Row::new(cells)
        });

        let title = format!(
            " Poses ({}) \u{b7} sorted by {} ",
            self.poses.len(),
            self.sort.column()
        );
        let bottom = match &self.error {
            Some(error) => Line::from(format!(" {error} ")).style(Style::default().fg(Color::Red)),
            None => Line::from(self.details(usize::from(area.width.saturating_sub(4)))),
        };
        let block = Block::default()
            .borders(Borders::ALL)
            .border_style(Style::default().fg(Color::DarkGray))
            .title(title)
            .title_bottom(bottom);
        let table = Table::new(rows, columns.iter().map(|&(_, width)| width))
            .header(header)
            .block(block)
            .row_highlight_style(
                Style::default()
                    .bg(Color::Rgb(45, 45, 70))
                    .add_modifier(Modifier::BOLD),
            );

        self.state.select(self.shown.map(|_| self.listed_at()));
        frame.render_stateful_widget(table, area, &mut self.state);
        self.area = Some(area);
        self.drawn.clone_from(&self.order);
    }

    /// What the pose drawn lacks, in about `width` columns: the reference's
    /// frame, when it is drawn in its own, the centres it inverts, the checks
    /// of its geometry it fails and the contacts it loses.
    fn details(&self, width: usize) -> String {
        let Some(pose) = self.shown_pose() else {
            return String::new();
        };
        let mut parts = Vec::new();
        if pose.superposition.is_none() {
            parts.push("in its own frame: too few pocket atoms matched".to_string());
        }
        if let Some(inverted) = pose.inverted.as_ref().filter(|centres| !centres.is_empty()) {
            parts.push(format!("inverts {}", inverted.join(", ")));
        }
        if let Some(failed) = pose.geometry.as_ref().filter(|failed| !failed.is_empty()) {
            parts.push(format!("strays in {}", failed.join(", ")));
        }
        if pose.lost.is_empty() {
            parts.push("keeps every reference contact".to_string());
        } else {
            let used: usize = parts.iter().map(|part| part.chars().count() + 3).sum();
            let room = width.saturating_sub(used + "lost: ".len() + 2);
            parts.push(format!("lost: {}", crate::ui::fitted(&pose.lost, room)));
        }

        format!(" {} ", parts.join(" \u{b7} "))
    }

    /// Whether a terminal cell falls in the table as last drawn.
    pub(crate) fn covers(&self, column: u16, row: u16) -> bool {
        self.visible
            && self
                .area
                .is_some_and(|area| area.contains(Position::new(column, row)))
    }

    /// The pose listed at a terminal cell of the table as last drawn.
    pub(crate) fn pose_at(&self, column: u16, row: u16) -> Option<usize> {
        let area = self.area?;
        // Rows start below the top border and the header.
        let first = area.y + 2;
        let inside = column > area.x && column + 1 < area.x + area.width;
        if !inside || row < first || row + 1 >= area.y + area.height {
            return None;
        }

        self.drawn
            .get(self.state.offset() + usize::from(row - first))
            .copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A pose scored with `clashes` and `inverted` centres, otherwise perfect.
    fn scored(clashes: usize, inverted: Option<&[&str]>) -> Pose {
        Pose {
            label: "pose.cif".to_string(),
            path: PathBuf::from("pose.cif"),
            model: 1,
            copy: None,
            similarity: 1.0,
            recovery: 1.0,
            rmsd: Some(0.0),
            clashes,
            inverted: inverted.map(|names| names.iter().map(|name| name.to_string()).collect()),
            geometry: None,
            lost: Vec::new(),
            ligand: Vec::new(),
            superposition: None,
            docked: None,
        }
    }

    #[test]
    fn a_label_too_long_for_its_column_keeps_its_end() {
        assert_eq!(elided("runs/seed-1/model_0.cif", 12), "\u{2026}model_0.cif");
        assert_eq!(elided("model_0.cif", 12), "model_0.cif");
    }

    #[test]
    fn the_sort_steps_through_clashes_and_inversions_back_to_the_rank() {
        let mut table = PoseTable::new(
            vec![
                scored(3, Some(&[])),
                scored(0, None),
                scored(3, Some(&["LIG B/1 C1"])),
                scored(1, Some(&[])),
            ],
            None,
        );

        table.cycle_sort();
        table.cycle_sort();
        assert_eq!(table.sort, PoseSort::Clashes);
        assert_eq!(table.order, [1, 3, 0, 2], "ties keep the rank");
        table.cycle_sort();
        assert_eq!(table.sort, PoseSort::Inverted);
        assert_eq!(table.order, [0, 3, 2, 1], "the unpaired last");
        table.cycle_sort();
        assert_eq!(table.order, [0, 1, 2, 3]);
    }
}
