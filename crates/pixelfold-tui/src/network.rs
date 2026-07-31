//! The linked residue-interaction-network pane.
//!
//! A force-directed 2D view of the residue interaction network, drawn on a
//! braille canvas beside the 3D structure. Clicking a node selects that residue,
//! which the structure pane highlights; picking a residue in the structure
//! highlights its node here. One `selected_atom_idx` drives both panes, so the
//! two views stay in sync.

use std::collections::HashMap;

use ratatui::prelude::*;
use ratatui::widgets::canvas::{Canvas, Circle, Line as CanvasLine, Points};
use ratatui::widgets::{Block, Borders};

use pixelfold_core::interactions::InteractionKind;
use pixelfold_core::rin::{self, Network};
use pixelfold_core::{Protein, interactions};

/// Fruchterman-Reingold steps run once when the pane is opened.
const LAYOUT_ITERATIONS: usize = 300;
/// A click within this many cells of a node selects it.
const PICK_RADIUS_CELLS: i32 = 3;

const SELECTED_COLOR: Color = Color::Rgb(0, 255, 255);
/// Radius of a hub's halo at the highest centrality; leaves get none.
const HUB_MAX_RADIUS: f64 = 0.035;
/// Normalised centrality below which a node draws as a bare point, no halo.
const HUB_THRESHOLD: f32 = 0.15;

pub struct NetworkView {
    network: Network,
    /// Layout position of each node in the unit square, y up.
    positions: Vec<(f32, f32)>,
    /// Each edge as its two node indices and its interaction kind, for colouring.
    edges: Vec<(usize, usize, InteractionKind)>,
    /// Normalised betweenness centrality per node (`sqrt` scaled to `0..1`), for
    /// hub sizing.
    hub: Vec<f32>,
    /// A representative atom for each node, for selecting into the 3D view.
    repr_atom: Vec<usize>,
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

    /// Lay out an already-built network and map its residues to `protein`'s
    /// atoms both ways. Split from [`build`] so the mapping can be tested without
    /// running the interaction engine.
    fn from_network(protein: &Protein, network: Network) -> Self {
        let positions: Vec<(f32, f32)> = rin::force_directed(&network, LAYOUT_ITERATIONS)
            .positions
            .iter()
            .map(|p| (p.x, p.y))
            .collect();

        let edges = {
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
                    Some((
                        *index.get(edge.source.as_str())?,
                        *index.get(edge.target.as_str())?,
                        edge.kind,
                    ))
                })
                .collect()
        };

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

        let key_index: HashMap<(String, u32, Option<char>), usize> = network
            .nodes
            .iter()
            .enumerate()
            .map(|(i, node)| ((node.chain.clone(), node.resi, node.icode), i))
            .collect();

        let mut repr_atom = vec![0usize; network.nodes.len()];
        let mut assigned = vec![false; network.nodes.len()];
        let mut node_of_atom = HashMap::new();
        for (atom_idx, atom) in protein.atoms.iter().enumerate() {
            let key = (
                atom.chain_id.as_str().to_string(),
                atom.residue_seq,
                atom.insertion_code,
            );
            if let Some(&node) = key_index.get(&key) {
                node_of_atom.insert(atom_idx, node);
                if !assigned[node] {
                    repr_atom[node] = atom_idx;
                    assigned[node] = true;
                }
            }
        }

        Self {
            network,
            positions,
            edges,
            hub,
            repr_atom,
            node_of_atom,
            node_cell: Vec::new(),
        }
    }

    /// A representative atom of a node, to select into the 3D view.
    pub fn atom_of_node(&self, node: usize) -> usize {
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
    pub fn render(&mut self, frame: &mut Frame, area: Rect, selected: Option<usize>) {
        // The drawable region inside the one-cell border. Node cell positions are
        // recorded against it so a click maps back to the nearest node.
        let inner_w = area.width.saturating_sub(2).max(1);
        let inner_h = area.height.saturating_sub(2).max(1);
        let (ox, oy) = (area.x + 1, area.y + 1);
        self.node_cell = self
            .positions
            .iter()
            .map(|&(x, y)| {
                let col = ox + (x.clamp(0.0, 1.0) * (inner_w - 1) as f32).round() as u16;
                // Canvas y runs bottom-up; screen rows run top-down.
                let row = oy + ((1.0 - y.clamp(0.0, 1.0)) * (inner_h - 1) as f32).round() as u16;
                (col, row)
            })
            .collect();

        let title = match selected {
            Some(node) => format!(
                " Interaction Network — {} {} ",
                self.network.nodes[node].id, self.network.nodes[node].resn
            ),
            None => format!(
                " Interaction Network — {} residues, {} edges ",
                self.network.nodes.len(),
                self.network.edges.len()
            ),
        };

        let legend = legend(&self.edges);

        let network = &self.network;
        let positions = &self.positions;
        let edges = &self.edges;
        let hub = &self.hub;

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
                for &(u, v, kind) in edges {
                    let (ux, uy) = positions[u];
                    let (vx, vy) = positions[v];
                    ctx.draw(&CanvasLine::new(
                        ux as f64,
                        uy as f64,
                        vx as f64,
                        vy as f64,
                        edge_color(kind),
                    ));
                }
                ctx.layer();

                for (i, &(x, y)) in positions.iter().enumerate() {
                    let coords = [(x as f64, y as f64)];
                    let color = node_color(network.nodes[i].ss);
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

/// The secondary-structure palette, keyed by the DSSP letter carried on a node.
fn node_color(ss: char) -> Color {
    match ss {
        'H' | 'G' | 'I' => Color::Rgb(255, 100, 100),
        'E' | 'B' => Color::Rgb(100, 100, 255),
        'T' | 'S' => Color::Rgb(255, 200, 100),
        _ => Color::Gray,
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
        InteractionKind::HalogenBond => "halogen",
        InteractionKind::WaterBridge => "water",
        InteractionKind::MetalCoordination => "metal",
        InteractionKind::Disulfide => "S-S",
    }
}

/// The interaction kinds in a fixed display order.
const KIND_ORDER: [InteractionKind; 9] = [
    InteractionKind::HydrogenBond,
    InteractionKind::SaltBridge,
    InteractionKind::Hydrophobic,
    InteractionKind::PiStacking,
    InteractionKind::PiCation,
    InteractionKind::HalogenBond,
    InteractionKind::WaterBridge,
    InteractionKind::MetalCoordination,
    InteractionKind::Disulfide,
];

/// A one-line legend of the interaction kinds present, each a coloured swatch and
/// its label, in the fixed order.
fn legend(edges: &[(usize, usize, InteractionKind)]) -> Line<'static> {
    let present: std::collections::BTreeSet<InteractionKind> =
        edges.iter().map(|&(_, _, kind)| kind).collect();

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

    fn atom(chain: &str, seq: u32, name: &str, residue: &str) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(&name[..1]),
            residue_name: FixedStr::new(residue),
            residue_seq: seq,
            chain_id: FixedStr::new(chain),
            is_hetatm: false,
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
            surface_points: Vec::new(),
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
        assert!(repr == 0 || repr == 1, "the node's atom is one of A/1's");

        // The non-interacting residue is not in the network.
        assert_eq!(view.node_of_atom(4), None);

        // One laid-out position per node.
        assert_eq!(view.positions.len(), view.network.nodes.len());
    }
}
