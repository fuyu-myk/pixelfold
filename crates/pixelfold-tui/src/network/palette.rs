//! The pane's colours: nodes by secondary structure or burial, edges by
//! interaction kind, brightened by energy or by how often they are present.

use ratatui::prelude::*;

use pixelfold_core::interactions::InteractionKind;

use super::EdgeVis;

/// The secondary-structure palette, keyed by the DSSP letter carried on a node.
pub(super) fn node_color(ss: char) -> Color {
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
pub(super) fn burial_color(rsa: Option<f32>) -> Color {
    match rsa {
        None => Color::DarkGray,
        Some(r) => {
            let exposed = r.clamp(0.0, 1.0);
            let lerp = |buried: f32, solvent: f32| (buried + (solvent - buried) * exposed) as u8;
            Color::Rgb(lerp(230.0, 70.0), lerp(140.0, 150.0), lerp(55.0, 235.0))
        }
    }
}

/// Perceptual brightness for an edge of a given summed energy, relative to the
/// strongest edge on screen. The scale is logarithmic so a covalent disulfide
/// (~250 kJ/mol) does not blow out the many weaker non-covalent edges, and it
/// never falls to zero, so even the faintest hydrophobic contact still draws.
pub(super) fn edge_intensity(energy: f32, max_energy: f32) -> f32 {
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
pub(super) fn frequency_intensity(frequency: f32) -> f32 {
    const FLOOR: f32 = 0.4;

    FLOOR + (1.0 - FLOOR) * frequency.clamp(0.0, 1.0)
}

/// Dim an RGB colour toward black by `t` in `0..1`; non-RGB colours pass through.
pub(super) fn scale_color(color: Color, t: f32) -> Color {
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
pub(super) fn edge_color(kind: InteractionKind) -> Color {
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
pub(super) fn kind_short(kind: InteractionKind) -> &'static str {
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
pub(super) const KIND_ORDER: [InteractionKind; 10] = [
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
pub(super) fn legend(edges: &[EdgeVis]) -> Line<'static> {
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
