//! iTerm2 inline-image protocol (also understood by WezTerm and others).
//!
//! The framebuffer is transmitted as a base64 PNG inside an OSC 1337 escape and
//! displayed at the cursor.

use crate::framebuffer::Framebuffer;
use crate::sink::{base64, png};

/// Encode a framebuffer as an iTerm2 inline-image escape ready to write to a
/// terminal.
pub fn encode(fb: &Framebuffer) -> Result<Vec<u8>, png::EncodingError> {
    let image = png::to_png_bytes(fb)?;
    let payload = base64::encode(&image);

    let mut out = Vec::with_capacity(payload.len() + 64);
    out.extend_from_slice(b"\x1b]1337;File=inline=1;size=");
    out.extend_from_slice(image.len().to_string().as_bytes());
    out.extend_from_slice(b";width=auto;height=auto:");
    out.extend_from_slice(payload.as_bytes());
    out.push(0x07); // BEL terminator

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wraps_a_base64_png_in_the_osc_1337_escape() {
        let fb = Framebuffer::new(3, 2, [0, 0, 0, 255]);
        let out = encode(&fb).expect("encode");

        let image = png::to_png_bytes(&fb).unwrap();
        let header = format!(
            "\x1b]1337;File=inline=1;size={};width=auto;height=auto:",
            image.len()
        );
        assert!(out.starts_with(header.as_bytes()), "OSC header wrong");
        assert_eq!(*out.last().unwrap(), 0x07, "missing BEL terminator");

        let body = &out[header.len()..out.len() - 1];
        assert_eq!(body, base64::encode(&image).as_bytes());
    }
}
