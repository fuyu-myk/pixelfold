//! Drawing the pane on a braille canvas, and mapping a click back to the node
//! drawn nearest it.

use ratatui::prelude::*;
use ratatui::widgets::canvas::{Canvas, Circle, Line as CanvasLine, Points};
use ratatui::widgets::{Block, Borders};

use pixelfold_render::renderer::Camera;

use super::palette::{
    burial_color, edge_color, edge_intensity, frequency_intensity, kind_short, legend, node_color,
    scale_color,
};
use super::{ChainFilter, EdgeVis, LayoutMode, NetworkView, NodeColorMode, edge_passes};

/// A click within this many cells of a node selects it.
const PICK_RADIUS_CELLS: i32 = 3;

const SELECTED_COLOR: Color = Color::Rgb(0, 255, 255);
/// Radius of a hub's halo at the highest centrality; leaves get none.
const HUB_MAX_RADIUS: f64 = 0.035;
/// Normalised centrality below which a node draws as a bare point, no halo.
const HUB_THRESHOLD: f32 = 0.15;

impl NetworkView {
    /// How brightly `edge` draws: by the fraction of the models it is present
    /// in for an ensemble, by its energy against the strongest on screen,
    /// `max_energy`, for one model.
    pub(super) fn brightness(&self, edge: &EdgeVis, max_energy: f32) -> f32 {
        if self.is_ensemble() {
            frequency_intensity(edge.frequency)
        } else {
            edge_intensity(edge.energy, max_energy)
        }
    }

    /// The node positions to draw: the fixed force-directed layout, or each
    /// residue projected through `camera` so the graph tracks the structure.
    fn display_positions(&self, camera: &Camera) -> Vec<(f32, f32)> {
        match self.mode {
            LayoutMode::Force => self.positions.clone(),
            LayoutMode::Spatial => {
                let raw = self
                    .node_pos3
                    .iter()
                    .map(|&pos| {
                        let (x, y, _) = camera.project_point_with_depth(pos, 1.0, 1.0);
                        (x, y)
                    })
                    .collect();
                normalize_positions(raw)
            }
        }
    }

    /// The node nearest a click, within the pick radius.
    pub fn pick(&self, column: u16, row: u16) -> Option<usize> {
        let mut best: Option<(usize, i32)> = None;
        for (i, &(c, r)) in self.node_cell.iter().enumerate() {
            let dc = c as i32 - column as i32;
            let dr = r as i32 - row as i32;
            let d2 = dc * dc + dr * dr;
            if d2 <= PICK_RADIUS_CELLS * PICK_RADIUS_CELLS && best.is_none_or(|(_, b)| d2 < b) {
                best = Some((i, d2));
            }
        }

        best.map(|(i, _)| i)
    }

    /// Render into `area`, emphasising `selected` when it is a node in this view.
    /// `camera` is used to place the nodes in spatial mode.
    pub fn render(
        &mut self,
        frame: &mut Frame,
        area: Rect,
        selected: Option<usize>,
        camera: &Camera,
    ) {
        let display = self.display_positions(camera);

        // The drawable region inside the one-cell border. Node cell positions are
        // recorded against it so a click maps back to the nearest node.
        let inner_w = area.width.saturating_sub(2).max(1);
        let inner_h = area.height.saturating_sub(2).max(1);
        let (ox, oy) = (area.x + 1, area.y + 1);
        self.node_cell = display
            .iter()
            .map(|&(x, y)| {
                let col = ox + (x.clamp(0.0, 1.0) * (inner_w - 1) as f32).round() as u16;
                // Canvas y runs bottom-up; screen rows run top-down.
                let row = oy + ((1.0 - y.clamp(0.0, 1.0)) * (inner_h - 1) as f32).round() as u16;
                (col, row)
            })
            .collect();

        // The bracket tags the layout, colour scheme, and any active filter.
        let mut tags = vec![
            match self.mode {
                LayoutMode::Force => "force",
                LayoutMode::Spatial => "spatial",
            },
            match self.color_mode {
                NodeColorMode::Structure => "struct",
                NodeColorMode::Burial => "burial",
            },
        ];
        if let Some(kind) = self.kind_filter {
            tags.push(kind_short(kind));
        }
        match self.chain_filter {
            ChainFilter::All => {}
            ChainFilter::IntraOnly => tags.push("intra"),
            ChainFilter::InterOnly => tags.push("inter"),
        }
        let floor = format!("≥{:.0}%", self.min_frequency * 100.0);
        if self.is_ensemble() {
            tags.push("ensemble");
            if self.min_frequency > 0.0 {
                tags.push(&floor);
            }
        }
        let tag = tags.join(" · ");

        let visible_edges = self
            .edges
            .iter()
            .filter(|e| edge_passes(e, self.kind_filter, self.chain_filter, self.min_frequency))
            .count();
        let models = match self.network.models {
            1 => String::new(),
            models => format!(", {models} models"),
        };
        let title = match selected {
            Some(node) => format!(
                " Interaction Network [{tag}] — {} {} ",
                self.network.nodes[node].id, self.network.nodes[node].resn
            ),
            None => format!(
                " Interaction Network [{tag}] — {} residues, {visible_edges} edges{models} ",
                self.network.nodes.len(),
            ),
        };

        let legend = legend(&self.edges);

        let network = &self.network;
        let positions = &display;
        let edges = &self.edges;
        let hub = &self.hub;
        let rsa = &self.rsa;
        let color_mode = self.color_mode;
        let kind_filter = self.kind_filter;
        let chain_filter = self.chain_filter;
        let min_frequency = self.min_frequency;
        // Brightest visible edge sets the top of the (log) intensity scale, so a
        // filtered view rescales to the edges it actually shows.
        let max_energy = edges
            .iter()
            .filter(|e| edge_passes(e, kind_filter, chain_filter, min_frequency))
            .map(|e| e.energy)
            .fold(0.0f32, f32::max);
        let brightness: Vec<f32> = edges
            .iter()
            .map(|e| self.brightness(e, max_energy))
            .collect();
        // Dimmest first: a residue pair joined by two kinds draws both along
        // one line, and the canvas keeps the last colour written to a cell.
        let mut order: Vec<usize> = (0..edges.len()).collect();
        order.sort_by(|&a, &b| brightness[a].total_cmp(&brightness[b]));

        let canvas = Canvas::default()
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .border_style(Style::default().fg(Color::DarkGray))
                    .title(title)
                    .title_bottom(legend),
            )
            .marker(symbols::Marker::Braille)
            .x_bounds([0.0, 1.0])
            .y_bounds([0.0, 1.0])
            .paint(move |ctx| {
                for &i in &order {
                    let (e, bright) = (&edges[i], brightness[i]);
                    if !edge_passes(e, kind_filter, chain_filter, min_frequency) {
                        continue;
                    }
                    let (ux, uy) = positions[e.u];
                    let (vx, vy) = positions[e.v];
                    let color = scale_color(edge_color(e.kind), bright);
                    ctx.draw(&CanvasLine::new(
                        ux as f64, uy as f64, vx as f64, vy as f64, color,
                    ));
                }
                ctx.layer();

                for (i, &(x, y)) in positions.iter().enumerate() {
                    let coords = [(x as f64, y as f64)];
                    let color = match color_mode {
                        NodeColorMode::Structure => node_color(network.nodes[i].ss),
                        NodeColorMode::Burial => burial_color(rsa[i]),
                    };
                    ctx.draw(&Points {
                        coords: &coords,
                        color,
                    });
                    // Betweenness hubs get a halo whose radius grows with rank.
                    if hub[i] > HUB_THRESHOLD {
                        ctx.draw(&Circle {
                            x: x as f64,
                            y: y as f64,
                            radius: hub[i] as f64 * HUB_MAX_RADIUS,
                            color,
                        });
                    }
                }

                if let Some(node) = selected {
                    let (x, y) = positions[node];
                    ctx.draw(&Circle {
                        x: x as f64,
                        y: y as f64,
                        radius: 0.03,
                        color: SELECTED_COLOR,
                    });
                }
            });

        frame.render_widget(canvas, area);
    }
}

/// Rescale positions to fill the unit square, preserving aspect. Used to fit the
/// spatial projection into the same `[0, 1]` frame the force layout uses.
pub(super) fn normalize_positions(mut raw: Vec<(f32, f32)>) -> Vec<(f32, f32)> {
    if raw.is_empty() {
        return raw;
    }

    let (mut min_x, mut min_y) = (f32::INFINITY, f32::INFINITY);
    let (mut max_x, mut max_y) = (f32::NEG_INFINITY, f32::NEG_INFINITY);
    for &(x, y) in &raw {
        min_x = min_x.min(x);
        min_y = min_y.min(y);
        max_x = max_x.max(x);
        max_y = max_y.max(y);
    }

    let extent = (max_x - min_x).max(max_y - min_y).max(1e-4);
    let (center_x, center_y) = ((min_x + max_x) * 0.5, (min_y + max_y) * 0.5);
    for p in &mut raw {
        p.0 = 0.5 + (p.0 - center_x) / extent;
        p.1 = 0.5 + (p.1 - center_y) / extent;
    }

    raw
}
