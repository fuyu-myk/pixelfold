//! The linked residue-interaction-network pane.
//!
//! A force-directed 2D view of the residue interaction network, drawn on a
//! braille canvas beside the 3D structure. Clicking a node selects that residue,
//! which the structure pane highlights; picking a residue in the structure
//! highlights its node here. One `selected_atom_idx` drives both panes, so the
//! two views stay in sync.
//!
//! For a multi-model file the pane can show the network of the whole ensemble
//! instead: every edge any model has, brightened by the fraction of models it
//! is present in, laid over the atoms of whichever model the 3D view shows.

use std::collections::HashMap;

use glam::Vec3;
use ratatui::prelude::*;
use ratatui::widgets::canvas::{Canvas, Circle, Line as CanvasLine, Points};
use ratatui::widgets::{Block, Borders};

use anyhow::Result;
use pixelfold_core::interactions::InteractionKind;
use pixelfold_core::parser::{Load, ParsedStructure};
use pixelfold_core::rin::{self, Ensemble, Network};
use pixelfold_core::sasa::{SurfaceCalculator, relative_sasa};
use pixelfold_core::{Protein, interactions};
use pixelfold_render::renderer::Camera;

/// How the network's nodes are placed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LayoutMode {
    /// Force-directed: positions from the graph's topology, independent of the
    /// structure's orientation.
    #[default]
    Force,
    /// Spatial: each residue at its 3D position projected through the camera, so
    /// the graph rotates with the structure.
    Spatial,
}

/// How the network's nodes are coloured.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NodeColorMode {
    /// By DSSP secondary structure, the same palette the 3D view uses.
    #[default]
    Structure,
    /// By burial, or relative solvent accessibility.
    Burial,
}

/// Fruchterman-Reingold steps run once when the pane is opened.
const LAYOUT_ITERATIONS: usize = 300;
/// A click within this many cells of a node selects it.
const PICK_RADIUS_CELLS: i32 = 3;

const SELECTED_COLOR: Color = Color::Rgb(0, 255, 255);
/// Radius of a hub's halo at the highest centrality; leaves get none.
const HUB_MAX_RADIUS: f64 = 0.035;
/// Normalised centrality below which a node draws as a bare point, no halo.
const HUB_THRESHOLD: f32 = 0.15;
/// The fractions of an ensemble's models an edge must be present in to draw,
/// stepped through in turn: every edge, then those in half the models, then
/// those in nearly all.
const FREQUENCY_FLOORS: [f32; 3] = [0.0, 0.5, 0.9];

/// What the pane shows: the network of the model the 3D view shows, or of
/// every model folded into one with the edges in fewer than a floor of them
/// hidden. Kept by the app rather than the view, so a key pressed while a view
/// is still being built never loses it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NetworkScope {
    #[default]
    Model,
    Ensemble {
        /// The index of the floor in [`FREQUENCY_FLOORS`].
        floor: usize,
    },
}

impl NetworkScope {
    /// The next scope: the ensemble after the model, each floor in turn, and
    /// the model again after the last.
    pub fn next(self) -> Self {
        match self {
            NetworkScope::Model => NetworkScope::Ensemble { floor: 0 },
            NetworkScope::Ensemble { floor } if floor + 1 < FREQUENCY_FLOORS.len() => {
                NetworkScope::Ensemble { floor: floor + 1 }
            }
            NetworkScope::Ensemble { .. } => NetworkScope::Model,
        }
    }

    pub fn is_ensemble(self) -> bool {
        matches!(self, NetworkScope::Ensemble { .. })
    }

    /// The least fraction of the models an edge drawn is present in.
    pub fn min_frequency(self) -> f32 {
        match self {
            NetworkScope::Model => 0.0,
            NetworkScope::Ensemble { floor } => FREQUENCY_FLOORS[floor],
        }
    }
}

/// One edge as the pane draws it: its two node indices, kind (which sets the
/// hue), summed energy (which sets the brightness of one model's network),
/// whether it crosses chains (for the chain filter), and the fraction of the
/// models it is present in (the brightness of an ensemble's network).
struct EdgeVis {
    u: usize,
    v: usize,
    kind: InteractionKind,
    energy: f32,
    interchain: bool,
    frequency: f32,
}

/// Which edges the chain filter keeps.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum ChainFilter {
    /// Every edge.
    #[default]
    All,
    /// Only edges within a single chain.
    IntraOnly,
    /// Only edges between different chains, the inter-chain interface.
    InterOnly,
}

pub struct NetworkView {
    network: Network,
    mode: LayoutMode,
    color_mode: NodeColorMode,
    /// Force-directed layout position of each node in the unit square, y up.
    positions: Vec<(f32, f32)>,
    /// The 3D position of each node's representative atom, projected through the
    /// camera in spatial mode.
    node_pos3: Vec<Vec3>,
    /// Each edge for drawing: node indices, kind (hue), and energy (brightness).
    edges: Vec<EdgeVis>,
    /// The interaction kinds actually present, in display order, for cycling the
    /// kind filter over only the kinds this network has.
    present_kinds: Vec<InteractionKind>,
    /// Show only this interaction kind's edges, or every kind when `None`.
    kind_filter: Option<InteractionKind>,
    /// Which edges the chain filter keeps.
    chain_filter: ChainFilter,
    /// The least fraction of an ensemble's models an edge drawn is present in.
    min_frequency: f32,
    /// How many of the models each pair of nodes (the lower index first)
    /// interacts in, by any kind.
    touches: HashMap<(usize, usize), usize>,
    /// Normalised betweenness centrality per node (`sqrt` scaled to `0..1`), for
    /// hub sizing.
    hub: Vec<f32>,
    /// Relative solvent accessibility per node, `None` for residues with no
    /// reference maximum (ligands, non-standard). Drives the burial colour mode.
    rsa: Vec<Option<f32>>,
    /// A representative atom for each node, for selecting into the 3D view;
    /// `None` for a residue of an ensemble's network the shown model lacks.
    repr_atom: Vec<Option<usize>>,
    /// The node each interacting atom belongs to, for reading the selection back.
    node_of_atom: HashMap<usize, usize>,
    /// Where each node last drew, in terminal cells, for click picking.
    node_cell: Vec<(u16, u16)>,
}

impl NetworkView {
    /// Detect the interactions, build and lay out the network, and map residues
    /// to atoms both ways.
    pub fn build(protein: &Protein) -> Self {
        let detected = interactions::detect(protein);
        let network = rin::build(protein, &detected);
        Self::from_network(protein, network)
    }

    /// Detect the interactions of every one of `models` of `source`, built as
    /// `load` says, and lay the network of them all over `protein`, the model
    /// shown.
    pub fn build_ensemble(
        protein: &Protein,
        source: &ParsedStructure,
        load: &Load,
        models: &[usize],
    ) -> Result<Self> {
        let mut ensemble = Ensemble::default();
        for &model in models {
            let built = source.protein(&Load {
                model: Some(model),
                ..load.clone()
            })?;
            ensemble.add(&built, &interactions::detect(&built));
        }

        Ok(Self::from_ensemble(protein, &ensemble))
    }

    /// Lay `ensemble`'s network over `protein`, knowing how many models each
    /// residue pair interacts in by any kind, which its edges, one per kind,
    /// cannot say.
    fn from_ensemble(protein: &Protein, ensemble: &Ensemble) -> Self {
        let mut view = Self::from_network(protein, ensemble.network(0.0));
        let index: HashMap<&str, usize> = view
            .network
            .nodes
            .iter()
            .enumerate()
            .map(|(i, node)| (node.id.as_str(), i))
            .collect();
        view.touches = ensemble
            .contacts()
            .map(|(source, target, models)| ((index[source], index[target]), models))
            .collect();

        view
    }

    /// Lay out an already-built network and map its residues to `protein`'s
    /// atoms both ways. Split from [`build`] so the mapping can be tested without
    /// running the interaction engine.
    fn from_network(protein: &Protein, network: Network) -> Self {
        let positions: Vec<(f32, f32)> = rin::force_directed(&network, LAYOUT_ITERATIONS)
            .positions
            .iter()
            .map(|p| (p.x, p.y))
            .collect();

        let edges: Vec<EdgeVis> = {
            let index: HashMap<&str, usize> = network
                .nodes
                .iter()
                .enumerate()
                .map(|(i, node)| (node.id.as_str(), i))
                .collect();
            network
                .edges
                .iter()
                .filter_map(|edge| {
                    Some(EdgeVis {
                        u: *index.get(edge.source.as_str())?,
                        v: *index.get(edge.target.as_str())?,
                        kind: edge.kind,
                        energy: edge.energy,
                        interchain: edge.interchain,
                        frequency: network.frequency(edge),
                    })
                })
                .collect()
        };

        // One model's pairs interact in it or not; an ensemble's are counted by
        // `from_ensemble`.
        let touches = edges.iter().map(|e| ((e.u, e.v), 1)).collect();

        // The kinds actually present, in the legend's fixed order.
        let present: std::collections::BTreeSet<InteractionKind> =
            edges.iter().map(|e| e.kind).collect();
        let present_kinds: Vec<InteractionKind> = KIND_ORDER
            .into_iter()
            .filter(|k| present.contains(k))
            .collect();

        // Size hubs by betweenness centrality; the square root lifts the many
        // mid-ranked residues out of the shadow of the few very high ones.
        let betweenness = rin::analyze(&network).betweenness;
        let max = betweenness.iter().copied().fold(0.0f64, f64::max);
        let hub: Vec<f32> = betweenness
            .iter()
            .map(|&b| {
                if max > 0.0 {
                    (b / max).sqrt() as f32
                } else {
                    0.0
                }
            })
            .collect();

        let atoms = AtomMap::of(protein, &network);

        Self {
            network,
            mode: LayoutMode::default(),
            color_mode: NodeColorMode::default(),
            positions,
            edges,
            present_kinds,
            kind_filter: None,
            chain_filter: ChainFilter::default(),
            min_frequency: FREQUENCY_FLOORS[0],
            touches,
            hub,
            rsa: atoms.rsa,
            repr_atom: atoms.repr_atom,
            node_of_atom: atoms.node_of_atom,
            node_pos3: atoms.node_pos3,
            node_cell: Vec::new(),
        }
    }

    /// This view laid over `protein`, another model of the same file: the
    /// network, its layout and every setting kept, the residues mapped to the
    /// new model's atoms and placed and buried as it has them.
    pub fn remapped(self, protein: &Protein) -> Self {
        let atoms = AtomMap::of(protein, &self.network);

        Self {
            rsa: atoms.rsa,
            repr_atom: atoms.repr_atom,
            node_of_atom: atoms.node_of_atom,
            node_pos3: atoms.node_pos3,
            ..self
        }
    }

    pub fn network(&self) -> &Network {
        &self.network
    }

    /// Whether the network is an ensemble's, drawn by how often each edge is
    /// present rather than by its energy.
    fn is_ensemble(&self) -> bool {
        self.network.models > 1
    }

    /// Hide the edges in fewer than `fraction` of the models.
    pub fn set_min_frequency(&mut self, fraction: f32) {
        self.min_frequency = fraction;
    }

    #[cfg(test)]
    pub fn min_frequency(&self) -> f32 {
        self.min_frequency
    }

    /// The residues `node` interacts with and how many of the models it does
    /// in, by any kind of interaction, the most often first.
    pub fn partners(&self, node: usize) -> Vec<(&str, usize)> {
        let mut partners: Vec<(usize, usize)> = self
            .touches
            .iter()
            .filter_map(|(&(u, v), &models)| match (u == node, v == node) {
                (true, _) => Some((v, models)),
                (_, true) => Some((u, models)),
                _ => None,
            })
            .collect();
        partners.sort_by_key(|&(partner, models)| (std::cmp::Reverse(models), partner));

        partners
            .into_iter()
            .map(|(partner, models)| (self.network.nodes[partner].id.as_str(), models))
            .collect()
    }

    /// How brightly `edge` draws: by the fraction of the models it is present
    /// in for an ensemble, by its energy against the strongest on screen,
    /// `max_energy`, for one model.
    fn brightness(&self, edge: &EdgeVis, max_energy: f32) -> f32 {
        if self.is_ensemble() {
            frequency_intensity(edge.frequency)
        } else {
            edge_intensity(edge.energy, max_energy)
        }
    }

    /// Switch between the force-directed and spatial layouts.
    pub fn toggle_mode(&mut self) {
        self.mode = match self.mode {
            LayoutMode::Force => LayoutMode::Spatial,
            LayoutMode::Spatial => LayoutMode::Force,
        };
    }

    /// Switch node colouring between secondary structure and burial.
    pub fn toggle_color(&mut self) {
        self.color_mode = match self.color_mode {
            NodeColorMode::Structure => NodeColorMode::Burial,
            NodeColorMode::Burial => NodeColorMode::Structure,
        };
    }

    /// Step the kind filter: every kind, then each present kind in turn, and
    /// back to every kind.
    pub fn cycle_kind_filter(&mut self) {
        self.kind_filter = match self.kind_filter {
            None => self.present_kinds.first().copied(),
            Some(current) => match self.present_kinds.iter().position(|&k| k == current) {
                Some(i) => self.present_kinds.get(i + 1).copied(),
                None => None,
            },
        };
    }

    /// Step the chain filter: all edges, then intra-chain only, then the
    /// inter-chain interface only, and back.
    pub fn cycle_chain_filter(&mut self) {
        self.chain_filter = match self.chain_filter {
            ChainFilter::All => ChainFilter::IntraOnly,
            ChainFilter::IntraOnly => ChainFilter::InterOnly,
            ChainFilter::InterOnly => ChainFilter::All,
        };
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

    /// A representative atom of a node, to select into the 3D view.
    pub fn atom_of_node(&self, node: usize) -> Option<usize> {
        self.repr_atom[node]
    }

    /// The node an atom belongs to, if the atom is in the network.
    pub fn node_of_atom(&self, atom: usize) -> Option<usize> {
        self.node_of_atom.get(&atom).copied()
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

/// A network's residues mapped to the atoms of one model, and what that model
/// says of each: where it sits and how buried it is.
struct AtomMap {
    /// A representative atom for each node, `None` for a residue of an
    /// ensemble's network the model lacks.
    repr_atom: Vec<Option<usize>>,
    /// The node each atom of a residue in the network belongs to.
    node_of_atom: HashMap<usize, usize>,
    /// Each node's representative atom's position; a residue the model lacks
    /// sits at the middle of those it has.
    node_pos3: Vec<Vec3>,
    /// Each node's relative solvent accessibility, `None` for residues with no
    /// reference maximum (ligands, non-standard) or none in the model.
    rsa: Vec<Option<f32>>,
}

impl AtomMap {
    fn of(protein: &Protein, network: &Network) -> Self {
        let key_index: HashMap<(String, i32, Option<char>), usize> = network
            .nodes
            .iter()
            .enumerate()
            .map(|(i, node)| ((node.chain.clone(), node.resi, node.icode), i))
            .collect();

        let mut repr_atom = vec![None; network.nodes.len()];
        let mut node_of_atom = HashMap::new();
        for (atom_idx, atom) in protein.atoms.iter().enumerate() {
            let key = (
                atom.chain_id.as_str().to_string(),
                atom.residue_seq,
                atom.insertion_code,
            );
            if let Some(&node) = key_index.get(&key) {
                node_of_atom.insert(atom_idx, node);
                repr_atom[node].get_or_insert(atom_idx);
            }
        }

        let present: Vec<Vec3> = repr_atom
            .iter()
            .flatten()
            .map(|&atom| protein.atoms[atom].position)
            .collect();
        let middle = present.iter().sum::<Vec3>() / present.len().max(1) as f32;
        let node_pos3 = repr_atom
            .iter()
            .map(|atom| atom.map_or(middle, |atom| protein.atoms[atom].position))
            .collect();

        let atom_sasa = SurfaceCalculator::default().calculate_atom_sasa(&protein.atoms);
        let mut node_sasa = vec![0.0f32; network.nodes.len()];
        for (&atom_idx, &node) in &node_of_atom {
            node_sasa[node] += atom_sasa[atom_idx];
        }
        let rsa: Vec<Option<f32>> = network
            .nodes
            .iter()
            .enumerate()
            .map(|(i, node)| repr_atom[i].and_then(|_| relative_sasa(node_sasa[i], &node.resn)))
            .collect();

        Self {
            repr_atom,
            node_of_atom,
            node_pos3,
            rsa,
        }
    }
}

/// Rescale positions to fill the unit square, preserving aspect. Used to fit the
/// spatial projection into the same `[0, 1]` frame the force layout uses.
fn normalize_positions(mut raw: Vec<(f32, f32)>) -> Vec<(f32, f32)> {
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

/// The secondary-structure palette, keyed by the DSSP letter carried on a node.
fn node_color(ss: char) -> Color {
    match ss {
        'H' | 'G' | 'I' => Color::Rgb(255, 100, 100),
        'E' | 'B' => Color::Rgb(100, 100, 255),
        'T' | 'S' => Color::Rgb(255, 200, 100),
        _ => Color::Gray,
    }
}

/// The burial palette, keyed by a node's relative solvent accessibility. Buried
/// residues (low RSA, the structural core) read warm; solvent-facing residues
/// (high RSA) read cool. Residues with no reference maximum (ligands,
/// non-standard) are drawn neutral. RSA above 1 is clamped: a residue more
/// exposed than the reference tripeptide still reads fully solvent-facing.
fn burial_color(rsa: Option<f32>) -> Color {
    match rsa {
        None => Color::DarkGray,
        Some(r) => {
            let exposed = r.clamp(0.0, 1.0);
            let lerp = |buried: f32, solvent: f32| (buried + (solvent - buried) * exposed) as u8;
            Color::Rgb(lerp(230.0, 70.0), lerp(140.0, 150.0), lerp(55.0, 235.0))
        }
    }
}

/// Whether an edge passes the kind, chain and frequency filters and should be
/// drawn.
fn edge_passes(
    edge: &EdgeVis,
    kind_filter: Option<InteractionKind>,
    chain_filter: ChainFilter,
    min_frequency: f32,
) -> bool {
    let kind_ok = kind_filter.is_none_or(|k| k == edge.kind);
    let chain_ok = match chain_filter {
        ChainFilter::All => true,
        ChainFilter::IntraOnly => !edge.interchain,
        ChainFilter::InterOnly => edge.interchain,
    };

    kind_ok && chain_ok && edge.frequency >= min_frequency
}

/// Perceptual brightness for an edge of a given summed energy, relative to the
/// strongest edge on screen. The scale is logarithmic so a covalent disulfide
/// (~250 kJ/mol) does not blow out the many weaker non-covalent edges, and it
/// never falls to zero, so even the faintest hydrophobic contact still draws.
fn edge_intensity(energy: f32, max_energy: f32) -> f32 {
    const FLOOR: f32 = 0.4;
    let denom = (1.0 + max_energy.max(0.0)).ln();
    let norm = if denom > 0.0 {
        (1.0 + energy.max(0.0)).ln() / denom
    } else {
        1.0
    };

    FLOOR + (1.0 - FLOOR) * norm.clamp(0.0, 1.0)
}

/// Brightness for an edge of an ensemble present in `frequency` of its models,
/// rising linearly from the same floor as [`edge_intensity`], so an edge in a
/// single model still draws.
fn frequency_intensity(frequency: f32) -> f32 {
    const FLOOR: f32 = 0.4;

    FLOOR + (1.0 - FLOOR) * frequency.clamp(0.0, 1.0)
}

/// Dim an RGB colour toward black by `t` in `0..1`; non-RGB colours pass through.
fn scale_color(color: Color, t: f32) -> Color {
    match color {
        Color::Rgb(r, g, b) => {
            let t = t.clamp(0.0, 1.0);
            Color::Rgb(
                (r as f32 * t) as u8,
                (g as f32 * t) as u8,
                (b as f32 * t) as u8,
            )
        }
        other => other,
    }
}

/// The edge palette, one colour per interaction kind.
/// Exhaustive on purpose: a new interaction type must choose a colour.
fn edge_color(kind: InteractionKind) -> Color {
    match kind {
        InteractionKind::HydrogenBond => Color::Rgb(90, 160, 230),
        InteractionKind::SaltBridge => Color::Rgb(240, 150, 60),
        InteractionKind::Hydrophobic => Color::Rgb(130, 140, 90),
        InteractionKind::PiStacking => Color::Rgb(180, 120, 240),
        InteractionKind::PiCation => Color::Rgb(230, 100, 190),
        InteractionKind::PiHydrogenBond => Color::Rgb(120, 220, 90),
        InteractionKind::HalogenBond => Color::Rgb(80, 210, 190),
        InteractionKind::WaterBridge => Color::Rgb(110, 150, 200),
        InteractionKind::MetalCoordination => Color::Rgb(210, 190, 80),
        InteractionKind::Disulfide => Color::Rgb(230, 220, 90),
    }
}

/// A short legend label for each interaction kind.
fn kind_short(kind: InteractionKind) -> &'static str {
    match kind {
        InteractionKind::HydrogenBond => "H-bond",
        InteractionKind::SaltBridge => "salt",
        InteractionKind::Hydrophobic => "phobic",
        InteractionKind::PiStacking => "π-stack",
        InteractionKind::PiCation => "π-cation",
        InteractionKind::PiHydrogenBond => "π-H",
        InteractionKind::HalogenBond => "halogen",
        InteractionKind::WaterBridge => "water",
        InteractionKind::MetalCoordination => "metal",
        InteractionKind::Disulfide => "S-S",
    }
}

/// The interaction kinds in a fixed display order.
const KIND_ORDER: [InteractionKind; 10] = [
    InteractionKind::HydrogenBond,
    InteractionKind::SaltBridge,
    InteractionKind::Hydrophobic,
    InteractionKind::PiStacking,
    InteractionKind::PiCation,
    InteractionKind::PiHydrogenBond,
    InteractionKind::HalogenBond,
    InteractionKind::WaterBridge,
    InteractionKind::MetalCoordination,
    InteractionKind::Disulfide,
];

/// A one-line legend of the interaction kinds present, each a coloured swatch and
/// its label, in the fixed order.
fn legend(edges: &[EdgeVis]) -> Line<'static> {
    let present: std::collections::BTreeSet<InteractionKind> =
        edges.iter().map(|e| e.kind).collect();

    let mut spans = Vec::new();
    for kind in KIND_ORDER.into_iter().filter(|k| present.contains(k)) {
        spans.push(Span::styled(" ● ", Style::default().fg(edge_color(kind))));
        spans.push(Span::styled(
            kind_short(kind),
            Style::default().fg(Color::Gray),
        ));
    }

    Line::from(spans)
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::Vec3;
    use pixelfold_core::interactions::InteractionKind;
    use pixelfold_core::{Atom, FixedStr, Interaction, SecondaryStructure};

    fn atom(chain: &str, seq: i32, name: &str, residue: &str) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(&name[..1]),
            residue_name: FixedStr::new(residue),
            residue_seq: seq,
            chain_id: FixedStr::new(chain),
            is_hetatm: false,
            molecule: pixelfold_core::Molecule::Unknown,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: Vec3::ZERO,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    fn hbond(a: usize, b: usize) -> Interaction {
        Interaction {
            kind: InteractionKind::HydrogenBond,
            atoms_a: vec![a],
            atoms_b: vec![b],
            distance: 3.0,
            angle: None,
            bridge: None,
        }
    }

    #[test]
    fn normalize_fits_positions_into_the_unit_square() {
        assert!(normalize_positions(vec![]).is_empty());
        assert_eq!(normalize_positions(vec![(5.0, 5.0)]), vec![(0.5, 0.5)]);

        // Aspect-preserving: the larger (y) axis spans the frame, both centred.
        let out = normalize_positions(vec![(2.0, 2.0), (4.0, 6.0)]);
        for &(x, y) in &out {
            assert!(
                (0.0..=1.0).contains(&x) && (0.0..=1.0).contains(&y),
                "{x},{y}"
            );
        }
        assert!((out[0].1 - 0.0).abs() < 1e-6 && (out[1].1 - 1.0).abs() < 1e-6);
    }

    #[test]
    fn the_residue_atom_mapping_is_consistent_both_ways() {
        // Residue A/1 (atoms 0,1) hydrogen-bonds A/2 (atoms 2,3); atom 4 is an
        // isolated residue that never interacts.
        let protein = Protein {
            atoms: vec![
                atom("A", 1, "N", "LYS"),
                atom("A", 1, "NZ", "LYS"),
                atom("A", 2, "O", "GLU"),
                atom("A", 2, "OE1", "GLU"),
                atom("B", 5, "CB", "LEU"),
            ],
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        };
        let network = rin::build(&protein, &[hbond(0, 2)]);
        let view = NetworkView::from_network(&protein, network);

        // Both atoms of A/1 map to A/1's node, and that node's representative
        // atom belongs to A/1.
        let node = view.node_of_atom(0).expect("atom 0 is in the network");
        assert_eq!(
            view.node_of_atom(1),
            Some(node),
            "both A/1 atoms share a node"
        );
        assert_ne!(view.node_of_atom(2), Some(node), "A/2 is a different node");
        let repr = view.atom_of_node(node);
        assert!(
            repr == Some(0) || repr == Some(1),
            "the node's atom is one of A/1's"
        );

        // The non-interacting residue is not in the network.
        assert_eq!(view.node_of_atom(4), None);

        // One laid-out position per node.
        assert_eq!(view.positions.len(), view.network.nodes.len());
        // And one burial value per node.
        assert_eq!(view.rsa.len(), view.network.nodes.len());
    }

    #[test]
    fn toggle_color_switches_the_node_scheme() {
        let protein = Protein {
            atoms: vec![atom("A", 1, "N", "LYS"), atom("A", 2, "O", "GLU")],
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        };
        let network = rin::build(&protein, &[hbond(0, 1)]);
        let mut view = NetworkView::from_network(&protein, network);

        assert_eq!(view.color_mode, NodeColorMode::Structure);
        view.toggle_color();
        assert_eq!(view.color_mode, NodeColorMode::Burial);
        view.toggle_color();
        assert_eq!(view.color_mode, NodeColorMode::Structure);
    }

    #[test]
    fn burial_color_separates_buried_from_exposed() {
        // Buried reads warm (red over blue), exposed reads cool (blue over red).
        let Color::Rgb(br, _, bb) = burial_color(Some(0.0)) else {
            panic!("expected an rgb colour");
        };
        assert!(br > bb, "a buried residue should read warm");

        let Color::Rgb(er, _, eb) = burial_color(Some(1.0)) else {
            panic!("expected an rgb colour");
        };
        assert!(eb > er, "an exposed residue should read cool");

        // Above the reference maximum still reads fully exposed, not past it.
        assert_eq!(burial_color(Some(2.0)), burial_color(Some(1.0)));
        // A residue with no reference maximum is neutral.
        assert_eq!(burial_color(None), Color::DarkGray);
    }

    #[test]
    fn edge_intensity_rises_with_energy_but_compresses_the_disulfide() {
        // Disulfide (251) sets the top; a hydrogen bond (20) and hydrophobic
        // contact (0.3) sit below it, brighter with energy but never dark.
        let max = 251.0;
        let weak = edge_intensity(0.3, max);
        let hbond = edge_intensity(20.0, max);
        let disulfide = edge_intensity(251.0, max);

        assert!(weak < hbond && hbond < disulfide, "intensity tracks energy");
        assert!(weak >= 0.4, "the weakest edge still draws above the floor");
        assert!((disulfide - 1.0).abs() < 1e-6, "the strongest edge is full");
        // Log compression: the 12x energy gap is nowhere near a 12x brightness
        // gap, so one covalent bond does not blow out the scale.
        assert!(
            disulfide / hbond < 2.0,
            "disulfide does not dominate visually"
        );
    }

    #[test]
    fn scale_color_dims_toward_black() {
        assert_eq!(
            scale_color(Color::Rgb(200, 100, 50), 0.5),
            Color::Rgb(100, 50, 25)
        );
        assert_eq!(
            scale_color(Color::Rgb(200, 100, 50), 1.0),
            Color::Rgb(200, 100, 50)
        );
        // A non-rgb colour is left untouched.
        assert_eq!(scale_color(Color::Gray, 0.5), Color::Gray);
    }

    #[test]
    fn edge_passes_honours_kind_and_chain_filters() {
        let intra_hbond = EdgeVis {
            u: 0,
            v: 1,
            kind: InteractionKind::HydrogenBond,
            energy: 20.0,
            interchain: false,
            frequency: 1.0,
        };
        let inter_salt = EdgeVis {
            u: 0,
            v: 2,
            kind: InteractionKind::SaltBridge,
            energy: 17.0,
            interchain: true,
            frequency: 1.0,
        };

        // No filter: everything passes.
        assert!(edge_passes(&intra_hbond, None, ChainFilter::All, 0.0));
        assert!(edge_passes(&inter_salt, None, ChainFilter::All, 0.0));

        // Kind filter keeps only the chosen kind.
        let only_hbond = Some(InteractionKind::HydrogenBond);
        assert!(edge_passes(&intra_hbond, only_hbond, ChainFilter::All, 0.0));
        assert!(!edge_passes(&inter_salt, only_hbond, ChainFilter::All, 0.0));

        // Chain filter splits intra from inter.
        assert!(edge_passes(&intra_hbond, None, ChainFilter::IntraOnly, 0.0));
        assert!(!edge_passes(
            &intra_hbond,
            None,
            ChainFilter::InterOnly,
            0.0
        ));
        assert!(edge_passes(&inter_salt, None, ChainFilter::InterOnly, 0.0));
        assert!(!edge_passes(&inter_salt, None, ChainFilter::IntraOnly, 0.0));
    }

    #[test]
    fn cycle_chain_filter_rotates_all_intra_inter() {
        let protein = Protein {
            atoms: vec![atom("A", 1, "N", "LYS"), atom("A", 2, "O", "GLU")],
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        };
        let network = rin::build(&protein, &[hbond(0, 1)]);
        let mut view = NetworkView::from_network(&protein, network);

        assert_eq!(view.chain_filter, ChainFilter::All);
        view.cycle_chain_filter();
        assert_eq!(view.chain_filter, ChainFilter::IntraOnly);
        view.cycle_chain_filter();
        assert_eq!(view.chain_filter, ChainFilter::InterOnly);
        view.cycle_chain_filter();
        assert_eq!(view.chain_filter, ChainFilter::All);
    }

    #[test]
    fn cycle_kind_filter_walks_present_kinds_then_clears() {
        // One residue pair joined by both a hydrogen bond and a salt bridge, so
        // two kinds are present.
        let protein = Protein {
            atoms: vec![
                atom("A", 1, "N", "LYS"),
                atom("A", 1, "NZ", "LYS"),
                atom("A", 2, "O", "GLU"),
                atom("A", 2, "OE1", "GLU"),
            ],
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        };
        let salt = Interaction {
            kind: InteractionKind::SaltBridge,
            atoms_a: vec![1],
            atoms_b: vec![3],
            distance: 3.5,
            angle: None,
            bridge: None,
        };
        let network = rin::build(&protein, &[hbond(0, 2), salt]);
        let mut view = NetworkView::from_network(&protein, network);

        assert_eq!(
            view.present_kinds,
            vec![InteractionKind::HydrogenBond, InteractionKind::SaltBridge]
        );
        assert_eq!(view.kind_filter, None);
        view.cycle_kind_filter();
        assert_eq!(view.kind_filter, Some(InteractionKind::HydrogenBond));
        view.cycle_kind_filter();
        assert_eq!(view.kind_filter, Some(InteractionKind::SaltBridge));
        view.cycle_kind_filter();
        assert_eq!(
            view.kind_filter, None,
            "past the last kind clears the filter"
        );
    }

    /// Three residues, A/1 hydrogen-bonded to A/2 and in hydrophobic contact
    /// with B/5, at `x` along the x axis.
    fn residues(x: f32) -> Protein {
        let at = |mut atom: Atom| {
            atom.position.x += x;
            atom
        };
        Protein {
            atoms: vec![
                at(atom("A", 1, "N", "LYS")),
                at(atom("A", 2, "O", "GLU")),
                at(atom("B", 5, "CB", "LEU")),
            ],
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        }
    }

    fn contact(a: usize, b: usize) -> Interaction {
        Interaction {
            kind: InteractionKind::Hydrophobic,
            ..hbond(a, b)
        }
    }

    /// Four models: the hydrogen bond in every one, the contact in one, and in
    /// the last A/1 and A/2 in hydrophobic contact with no hydrogen bond.
    fn ensemble() -> rin::Ensemble {
        let mut ensemble = rin::Ensemble::default();
        ensemble.add(&residues(0.0), &[hbond(0, 1), contact(0, 2)]);
        for _ in 0..2 {
            ensemble.add(&residues(0.0), &[hbond(0, 1)]);
        }
        ensemble.add(&residues(0.0), &[contact(0, 1)]);
        ensemble
    }

    #[test]
    fn an_ensemble_brightens_frequent_edges_and_a_floor_hides_rare_ones() {
        let mut view = NetworkView::from_ensemble(&residues(0.0), &ensemble());
        let edge = |kind, v| {
            view.edges
                .iter()
                .find(|e| e.kind == kind && e.v == v)
                .expect("an edge of the kind")
        };
        let (hbond, contact) = (
            edge(InteractionKind::HydrogenBond, 1),
            edge(InteractionKind::Hydrophobic, 2),
        );
        assert_eq!((hbond.frequency, contact.frequency), (0.75, 0.25));
        assert!(view.brightness(hbond, 0.0) > view.brightness(contact, 0.0));
        let shown = |view: &NetworkView| {
            view.edges
                .iter()
                .filter(|e| edge_passes(e, None, ChainFilter::All, view.min_frequency))
                .count()
        };
        assert_eq!(shown(&view), 3);

        view.set_min_frequency(0.5);
        assert_eq!(shown(&view), 1);
    }

    #[test]
    fn the_scope_steps_through_the_floors_and_back_to_the_model() {
        let mut scope = NetworkScope::default();
        let mut floors = Vec::new();
        for _ in 0..4 {
            scope = scope.next();
            floors.push(scope.is_ensemble().then(|| scope.min_frequency()));
        }

        assert_eq!(floors, [Some(0.0), Some(0.5), Some(0.9), None]);
    }

    // A/1 and A/2 hydrogen-bond in three models and touch only hydrophobically
    // in the fourth: they interact in every model, though no edge says so.
    #[test]
    fn a_residue_s_partners_count_the_models_they_touch_it_in_by_any_kind() {
        let view = NetworkView::from_ensemble(&residues(0.0), &ensemble());
        let a1 = view.node_of_atom(0).expect("A/1 is in the network");

        assert_eq!(view.partners(a1), [("A/2", 4), ("B/5", 1)]);
    }

    #[test]
    fn a_network_laid_over_another_model_keeps_its_layout_and_settings() {
        let mut view = NetworkView::from_ensemble(&residues(0.0), &ensemble());
        view.toggle_mode();
        view.set_min_frequency(0.5);
        let positions = view.positions.clone();

        let moved = view.remapped(&residues(10.0));

        assert_eq!(moved.positions, positions);
        assert_eq!(moved.mode, LayoutMode::Spatial);
        assert_eq!(moved.min_frequency, 0.5);
        assert_eq!(moved.partners(0), [("A/2", 4), ("B/5", 1)]);
        assert_eq!(moved.node_pos3[0].x, 10.0, "the new model's coordinates");
        assert_eq!(moved.node_of_atom(2), Some(2));
    }

    #[test]
    fn a_residue_the_shown_model_lacks_maps_to_no_atom() {
        let view = NetworkView::from_ensemble(&residues(0.0), &ensemble());
        let without_b5 = Protein {
            atoms: residues(0.0).atoms[..2].to_vec(),
            ..residues(0.0)
        };

        let moved = view.remapped(&without_b5);

        assert_eq!(moved.atom_of_node(2), None);
        assert_eq!(moved.rsa[2], None);
        assert_eq!(moved.atom_of_node(1), Some(1));
    }
}
