use image::{DynamicImage, RgbaImage};
use ratatui::prelude::*;
use ratatui::widgets::{Block, Borders, Paragraph};
use ratatui_image::picker::{Picker, ProtocolType};
use ratatui_image::{Resize, StatefulImage};

use pixelfold_core::{DisplayMode, SecondaryStructure};
use pixelfold_render::Framebuffer;

use crate::{App, RenderedFrame, render};

/// Safety cap on a framebuffer dimension.
const MAX_FB_DIMENSION: u32 = 8192;

/// Opaque backing colour, matching the framebuffer's own background.
const BG: Color = Color::Rgb(
    render::BACKGROUND[0],
    render::BACKGROUND[1],
    render::BACKGROUND[2],
);

/// Rows reserved under the network graph for the selected residue's summary.
const SUMMARY_HEIGHT: u16 = 7;

pub(crate) fn ui(frame: &mut Frame, app: &mut App, picker: &Picker) {
    let area = frame.area();

    frame.render_widget(Block::default().style(Style::default().bg(BG)), area);

    if app.protein.is_none() {
        render_help(frame, area);
        return;
    }

    let assembly_note = assembly_note(app);

    let header_height = if assembly_note.is_some() { 2 } else { 1 };
    let vertical_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(header_height),
            Constraint::Min(0),
            Constraint::Length(1),
        ])
        .split(area);

    let (info_bar_area, assembly_note_area) = if assembly_note.is_some() {
        let header_chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([Constraint::Length(1), Constraint::Length(1)])
            .split(vertical_chunks[0]);
        (header_chunks[0], Some(header_chunks[1]))
    } else {
        (vertical_chunks[0], None)
    };
    let content_area = vertical_chunks[1];
    let footer_area = vertical_chunks[2];

    // The right slot holds the network pane (wider) when it is open, otherwise
    // the atom inspect panel when a residue is selected.
    let (main_area, right_area) = if app.show_network {
        let chunks = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(55), Constraint::Percentage(45)])
            .split(content_area);
        (chunks[0], Some(chunks[1]))
    } else if app.selected_atom_idx.is_some() {
        let chunks = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(75), Constraint::Percentage(25)])
            .split(content_area);
        (chunks[0], Some(chunks[1]))
    } else {
        (content_area, None)
    };

    let (font_w, font_h) = picker.font_size();
    let fb_w = (main_area.width as u32 * font_w as u32).clamp(1, MAX_FB_DIMENSION);
    let fb_h = (main_area.height as u32 * font_h as u32).clamp(1, MAX_FB_DIMENSION);

    // The scene is rasterised into an owned framebuffer, so the immutable borrow
    // of `app` ends before the frame is stashed back into it for picking.
    let fb = render::render_scene(app, fb_w, fb_h);

    if main_area.width > 0 && main_area.height > 0 {
        let use_inplace =
            picker.protocol_type() == ProtocolType::Iterm2 && std::env::var_os("TMUX").is_none();

        if use_inplace {
            crate::iterm2::place(frame.buffer_mut(), &fb, main_area, BG);
        } else if let Some(dyn_img) = framebuffer_to_image(&fb) {
            let mut protocol = picker.new_resize_protocol(dyn_img);
            frame.render_stateful_widget(
                StatefulImage::default().resize(Resize::Fit(None)),
                main_area,
                &mut protocol,
            );
        }
    }

    app.last_frame = Some(RenderedFrame {
        fb,
        area: main_area,
        font: (font_w, font_h),
    });

    render_info_bar(frame, app, info_bar_area);

    if let (Some(note_area), Some(note)) = (assembly_note_area, assembly_note) {
        frame.render_widget(
            Paragraph::new(note)
                .alignment(Alignment::Center)
                .style(Style::default().fg(Color::Yellow).bold()),
            note_area,
        );
    }

    if app.show_network {
        let selected_node = app.selected_atom_idx.and_then(|atom| {
            app.network_view
                .as_ref()
                .and_then(|view| view.node_of_atom(atom))
        });
        let (net_area, summary_area) = match right_area {
            Some(right) if app.selected_atom_idx.is_some() => {
                let chunks = Layout::default()
                    .direction(Direction::Vertical)
                    .constraints([Constraint::Min(0), Constraint::Length(SUMMARY_HEIGHT)])
                    .split(right);
                (Some(chunks[0]), Some(chunks[1]))
            }
            other => (other, None),
        };
        if let (Some(area), Some(view)) = (net_area, app.network_view.as_mut()) {
            view.render(frame, area, selected_node, &app.camera);
        }
        app.network_area = net_area;
        if let (Some(area), Some(atom_idx)) = (summary_area, app.selected_atom_idx) {
            render_residue_summary(frame, app, atom_idx, area);
        }
    } else {
        app.network_area = None;
        if let (Some(right), Some(atom_idx)) = (right_area, app.selected_atom_idx) {
            render_inspect_panel(frame, app, atom_idx, right);
        }
    }

    frame.render_widget(
        Paragraph::new(
            "WASD rotate  ·  +/- zoom  ·  arrows pan  ·  i inspect  ·  g network  ·  l layout  ·  1/2 atoms/Cα  ·  c b h  ·  k [ ] , . slab  ·  f reset  ·  q quit",
        )
        .alignment(Alignment::Center)
        .style(Style::default().fg(Color::DarkGray)),
        footer_area,
    );
}

/// A header line for a structure whose biological assembly differs from the
/// deposited coordinates, so the partial view is not mistaken for complete.
fn assembly_note(app: &App) -> Option<String> {
    let assembly = app.protein.as_ref()?.assembly.as_ref()?;
    Some(match &assembly.oligomer {
        Some(oligomer) => format!(
            " Showing asymmetric unit; biological assembly is {} ({} symmetry operators, not yet generated) ",
            oligomer, assembly.operator_count
        ),
        None => format!(
            " Showing asymmetric unit; biological assembly applies {} symmetry operators (not yet generated) ",
            assembly.operator_count
        ),
    })
}

fn framebuffer_to_image(fb: &Framebuffer) -> Option<DynamicImage> {
    let buf: Vec<u8> = fb.color().iter().flatten().copied().collect();
    RgbaImage::from_raw(fb.width(), fb.height(), buf).map(DynamicImage::ImageRgba8)
}

fn render_info_bar(frame: &mut Frame, app: &App, area: Rect) {
    let Some(protein) = app.protein.as_ref() else {
        return;
    };

    let display_mode_text = match app.display_mode {
        DisplayMode::AllAtoms => "All Atoms",
        DisplayMode::Backbone => "C-Alpha",
    };
    let connections_text = if app.show_connections {
        "Connected"
    } else {
        "Dots"
    };
    let color_mode_text = if app.use_bfactor_colors {
        "B-factor"
    } else {
        "Sec. Struct."
    };
    let mode_text = if app.inspect_mode { " | [INSPECT]" } else { "" };
    let highlight_text = if app.residue_highlight && app.selected_atom_idx.is_some() {
        " | [RESIDUE]"
    } else {
        ""
    };
    let hbond_text = if app.show_hydrogen_bonds {
        let visible = protein
            .hbonds
            .iter()
            .filter(|hb| hb.energy < app.hbond_energy_threshold)
            .count();
        format!(
            " | [H-BONDS: {} @ {:.1} kcal/mol]",
            visible, app.hbond_energy_threshold
        )
    } else {
        String::new()
    };
    let slab_text = match app.slab {
        Some((far, near)) => format!(" | [SLAB {far:.2}-{near:.2}]"),
        None => String::new(),
    };

    let atom_count = match app.display_mode {
        DisplayMode::AllAtoms => protein.atoms.len(),
        DisplayMode::Backbone => app.ca_indices.len(),
    };

    let info_text = format!(
        " {} | {} atoms | [{}] [{}] [{}] | Zoom: {:.1}x{}{}{}{} ",
        protein.title,
        atom_count,
        display_mode_text,
        connections_text,
        color_mode_text,
        app.camera.zoom,
        mode_text,
        highlight_text,
        hbond_text,
        slab_text
    );

    frame.render_widget(
        Paragraph::new(info_text)
            .alignment(Alignment::Center)
            .style(Style::default().fg(Color::White)),
        area,
    );
}

fn render_inspect_panel(frame: &mut Frame, app: &App, atom_idx: usize, area: Rect) {
    let Some(protein) = app.protein.as_ref() else {
        return;
    };
    let atom = &protein.atoms[atom_idx];

    let ss_name = ss_name(atom.secondary_structure);
    let (r, g, b) = atom.secondary_structure.color_rgb();
    let ss_color = Color::Rgb(r, g, b);

    let label = Style::default().fg(Color::Gray);
    let value = Style::default().fg(Color::White);

    let info_lines = vec![
        Line::from(vec![Span::styled(
            "┌─── Atom Info ───┐",
            Style::default().fg(Color::Yellow).bold(),
        )]),
        Line::from(""),
        Line::from(vec![
            Span::styled("  Atom:    ", label),
            Span::styled(
                atom.name.to_string(),
                Style::default().fg(Color::Yellow).bold(),
            ),
        ]),
        Line::from(vec![
            Span::styled("  Serial:  ", label),
            Span::styled(format!("{}", atom.serial), value),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("  Residue: ", label),
            Span::styled(
                atom.residue_name.to_string(),
                Style::default().fg(Color::White).bold(),
            ),
        ]),
        Line::from(vec![
            Span::styled("  Number:  ", label),
            Span::styled(format!("{}", atom.residue_seq), value),
        ]),
        Line::from(vec![
            Span::styled("  Chain:   ", label),
            Span::styled(
                atom.chain_id.to_string(),
                Style::default().fg(Color::Cyan).bold(),
            ),
        ]),
        Line::from(""),
        Line::from(vec![Span::styled("  Structure:", label)]),
        Line::from(vec![Span::styled(
            format!("  {ss_name}"),
            Style::default().fg(ss_color).bold(),
        )]),
        Line::from(""),
        Line::from(vec![Span::styled("  Position:", label)]),
        Line::from(vec![Span::styled(
            format!("    x: {:.2}", atom.position.x),
            value,
        )]),
        Line::from(vec![Span::styled(
            format!("    y: {:.2}", atom.position.y),
            value,
        )]),
        Line::from(vec![Span::styled(
            format!("    z: {:.2}", atom.position.z),
            value,
        )]),
        Line::from(""),
        Line::from(vec![
            Span::styled("  B-factor: ", label),
            Span::styled(format!("{:.2}", atom.b_factor), value),
        ]),
        Line::from(""),
        Line::from(vec![Span::styled(
            if app.residue_highlight {
                "  Press 'r' to hide residue"
            } else {
                "  Press 'r' to highlight residue"
            },
            Style::default().fg(Color::DarkGray).italic(),
        )]),
        Line::from(""),
        Line::from(vec![Span::styled(
            "└─────────────────┘",
            Style::default().fg(Color::Yellow),
        )]),
    ];

    frame.render_widget(
        Paragraph::new(info_lines)
            .style(Style::default().bg(BG))
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .border_style(Style::default().fg(Color::Yellow)),
            ),
        area,
    );
}

fn ss_name(ss: SecondaryStructure) -> &'static str {
    match ss {
        SecondaryStructure::Helix => "α-Helix",
        SecondaryStructure::Sheet => "β-Sheet",
        SecondaryStructure::Turn => "Turn",
        SecondaryStructure::Coil => "Coil",
    }
}

/// The compact residue block shown under the network graph: the identity and
/// secondary structure of the selected residue, at the residue scale the network
/// is read at.
fn render_residue_summary(frame: &mut Frame, app: &App, atom_idx: usize, area: Rect) {
    let Some(protein) = app.protein.as_ref() else {
        return;
    };

    let atom = &protein.atoms[atom_idx];
    let (r, g, b) = atom.secondary_structure.color_rgb();

    let label = Style::default().fg(Color::Gray);
    let lines = vec![
        Line::from(vec![
            Span::styled(
                format!("{}/{}", atom.chain_id, atom.residue_seq),
                Style::default().fg(Color::Cyan).bold(),
            ),
            Span::styled(
                format!("  {}", atom.residue_name),
                Style::default().fg(Color::White).bold(),
            ),
        ]),
        Line::from(Span::styled(
            ss_name(atom.secondary_structure),
            Style::default().fg(Color::Rgb(r, g, b)).bold(),
        )),
        Line::from(vec![
            Span::styled("B-factor ", label),
            Span::styled(
                format!("{:.2}", atom.b_factor),
                Style::default().fg(Color::White),
            ),
        ]),
    ];

    frame.render_widget(
        Paragraph::new(lines).style(Style::default().bg(BG)).block(
            Block::default()
                .borders(Borders::ALL)
                .border_style(Style::default().fg(Color::DarkGray))
                .title(" Residue "),
        ),
        area,
    );
}

fn render_help(frame: &mut Frame, area: Rect) {
    let help_text = [
        "PixelFold - 3D Protein Viewer",
        "",
        "Usage: pixelfold <protein.pdb>",
        "",
        "Supports PDB and mmCIF formats",
        "",
        "Press 'q' to quit",
    ];

    frame.render_widget(
        Paragraph::new(help_text.join("\n"))
            .alignment(Alignment::Center)
            .style(Style::default().fg(Color::Gray)),
        area,
    );
}
