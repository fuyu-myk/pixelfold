//! HalfBlock Unicode fallback: two vertical pixels per character cell.
//!
//! The upper half block carries a truecolor foreground (the top pixel) over a
//! truecolor background (the bottom pixel), so each pixel keeps its own colour.
//! This is the floor of the fallback ladder for terminals with no image
//! protocol.

use crate::framebuffer::Framebuffer;

const UPPER_HALF: &str = "\u{2580}"; // the upper half block, over the lower half

/// Append a `u8` in decimal without allocating or a formatter.
fn push_u8(out: &mut String, n: u8) {
    if n >= 100 {
        out.push((b'0' + n / 100) as char);
    }
    if n >= 10 {
        out.push((b'0' + (n / 10) % 10) as char);
    }
    out.push((b'0' + n % 10) as char);
}

fn push_color(out: &mut String, lead: &str, c: [u8; 4]) {
    out.push_str(lead);
    push_u8(out, c[0]);
    out.push(';');
    push_u8(out, c[1]);
    out.push(';');
    push_u8(out, c[2]);
    out.push('m');
}

/// Encode a framebuffer as ANSI truecolor half-block text. Each line ends with a
/// colour reset and a newline, so the result prints directly to a truecolor
/// terminal. An odd final pixel row pairs with the background.
pub fn encode(fb: &Framebuffer) -> String {
    let (w, h) = (fb.width(), fb.height());
    let mut out = String::with_capacity(w as usize * h as usize / 2 * 40);

    let mut y = 0;
    while y < h {
        for x in 0..w {
            let top = fb.pixel(x, y);
            let bottom = if y + 1 < h {
                fb.pixel(x, y + 1)
            } else {
                fb.background()
            };

            push_color(&mut out, "\x1b[38;2;", top);
            push_color(&mut out, "\x1b[48;2;", bottom);
            out.push_str(UPPER_HALF);
        }

        out.push_str("\x1b[0m\n");
        y += 2;
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn push_u8_covers_the_range() {
        let mut s = String::new();
        for (n, want) in [(0u8, "0"), (5, "5"), (42, "42"), (100, "100"), (255, "255")] {
            s.clear();
            push_u8(&mut s, n);
            assert_eq!(s, want);
        }
    }

    #[test]
    fn two_rows_collapse_to_one_line_with_fg_top_and_bg_bottom() {
        let mut fb = Framebuffer::new(1, 2, [0, 0, 0, 255]);
        fb.test_and_set(0, 0, 0.0, [10, 20, 30, 255], 0); // top
        fb.test_and_set(0, 1, 0.0, [40, 50, 60, 255], 0); // bottom

        let out = encode(&fb);
        assert_eq!(
            out,
            "\x1b[38;2;10;20;30m\x1b[48;2;40;50;60m\u{2580}\x1b[0m\n"
        );
    }

    #[test]
    fn odd_height_pairs_the_last_row_with_the_background() {
        let bg = [7, 8, 9, 255];
        let mut fb = Framebuffer::new(1, 3, bg);
        fb.test_and_set(0, 2, 0.0, [100, 100, 100, 255], 0);

        let out = encode(&fb);
        // Two text lines: rows 0-1, then row 2 over the background.
        assert_eq!(out.matches('\n').count(), 2);
        assert!(
            out.contains("\x1b[48;2;7;8;9m"),
            "last row must sit on the background"
        );
    }
}
