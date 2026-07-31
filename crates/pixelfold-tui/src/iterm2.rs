//! Opaque, in-place iTerm2 image placement.
//!
//! ratatui-image's iTerm2 path erases the image region before repainting it,
//! which it needs for transparent images. On each redraw that erase shows a
//! one-frame gap between the old image and the new one (a flash). The framebuffer
//! is fully opaque, so the new image can be drawn straight over the old with no
//! erase and no gap. This writes a bespoke OSC 1337 escape into the ratatui
//! buffer the way ratatui-image does, minus the erase, and sized to the exact
//! cell box so it never overflows its area.

use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::style::Color;

use pixelfold_render::Framebuffer;
use pixelfold_render::sink::base64;

/// The inline-image escape for `fb`, scaled to exactly `cols` x `rows` cells,
/// opaque, and leaving the cursor in place so it composes inside a ratatui frame.
pub(crate) fn escape(fb: &Framebuffer, cols: u16, rows: u16) -> Option<String> {
    let png = pixelfold_render::to_png_bytes_fast(fb).ok()?;
    let payload = base64::encode(&png);

    let mut out = String::with_capacity(payload.len() + 96);
    out.push_str("\x1b]1337;File=inline=1;size=");
    out.push_str(&png.len().to_string());
    out.push_str(";width=");
    out.push_str(&cols.to_string());
    out.push_str(";height=");
    out.push_str(&rows.to_string());
    out.push_str(";preserveAspectRatio=0;doNotMoveCursor=1:");
    out.push_str(&payload);
    out.push('\x07');

    Some(out)
}

/// Place `fb` into `area` of `buf` as one opaque image. The escape lives in the
/// top-left cell; the rest of the area is marked skipped so ratatui does not
/// paint over the image.
pub(crate) fn place(buf: &mut Buffer, fb: &Framebuffer, area: Rect, bg: Color) {
    if area.width == 0 || area.height == 0 {
        return;
    }
    let Some(escape) = escape(fb, area.width, area.height) else {
        return;
    };

    if let Some(cell) = buf.cell_mut((area.x, area.y)) {
        cell.set_symbol(&escape);
        cell.set_bg(bg);
    }
    for y in area.top()..area.bottom() {
        for x in area.left()..area.right() {
            if x == area.x && y == area.y {
                continue;
            }

            if let Some(cell) = buf.cell_mut((x, y)) {
                cell.set_skip(true);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escape_is_opaque_in_place_and_cell_sized() {
        let fb = Framebuffer::new(4, 4, [13, 13, 18, 255]);
        let seq = escape(&fb, 20, 10).expect("encode");

        assert!(seq.starts_with("\x1b]1337;File=inline=1;size="));
        assert!(seq.contains(";width=20;height=10;"));
        assert!(seq.contains("preserveAspectRatio=0"));
        assert!(seq.contains("doNotMoveCursor=1"));
        assert!(seq.ends_with('\u{7}'), "missing BEL terminator");
        // No erase sequences (the whole point).
        assert!(!seq.contains("\x1b[") || !seq.contains("X\x1b["));
    }
}
