//! Kitty graphics protocol encoder.
//!
//! The framebuffer is transmitted as a PNG, so the payload is compressed, which
//! is what makes full-resolution rendering viable over SSH. The base64 payload
//! is split into 4096-byte chunks the way the protocol requires, and the image
//! is displayed at the cursor.

use crate::framebuffer::Framebuffer;
use crate::sink::{base64, png};

/// Maximum base64 payload bytes per escape chunk, fixed by the protocol.
const CHUNK: usize = 4096;

/// Encode a framebuffer as a Kitty escape sequence that transmits and displays
/// it. The bytes are ready to write to a terminal.
pub fn encode(fb: &Framebuffer) -> Result<Vec<u8>, png::EncodingError> {
    let payload = base64::encode(&png::to_png_bytes(fb)?);
    Ok(escape(payload.as_bytes()))
}

/// Wrap a base64 payload in the chunked APC escapes. The first chunk carries the
/// control data (`a=T` transmit-and-display, `f=100` PNG); every chunk carries
/// `m=1` until the last, which carries `m=0`.
fn escape(payload: &[u8]) -> Vec<u8> {
    if payload.is_empty() {
        return Vec::new();
    }

    let mut out = Vec::with_capacity(payload.len() + payload.len() / CHUNK * 16 + 32);
    let mut chunks = payload.chunks(CHUNK).peekable();
    let mut first = true;
    while let Some(chunk) = chunks.next() {
        let last = chunks.peek().is_none();
        out.extend_from_slice(b"\x1b_G");
        if first {
            out.extend_from_slice(b"a=T,f=100,");
            first = false;
        }

        out.extend_from_slice(if last { b"m=0" } else { b"m=1" });
        out.push(b';');
        out.extend_from_slice(chunk);
        out.extend_from_slice(b"\x1b\\");
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn count(haystack: &[u8], needle: &[u8]) -> usize {
        haystack
            .windows(needle.len())
            .filter(|w| *w == needle)
            .count()
    }

    #[test]
    fn single_chunk_is_a_complete_transmit_and_display_escape() {
        let fb = Framebuffer::new(2, 2, [0, 0, 0, 255]);
        let out = encode(&fb).expect("encode");

        assert!(
            out.starts_with(b"\x1b_Ga=T,f=100,m=0;"),
            "control data wrong"
        );
        assert!(out.ends_with(b"\x1b\\"), "missing ST terminator");
        assert_eq!(count(&out, b"\x1b_G"), 1, "small image should be one chunk");

        // The payload between the control data and the terminator is the base64
        // PNG.
        let expected = base64::encode(&png::to_png_bytes(&fb).unwrap());
        let start = b"\x1b_Ga=T,f=100,m=0;".len();
        let body = &out[start..out.len() - 2];
        assert_eq!(body, expected.as_bytes());
    }

    /// Split the escape stream into its `(control, payload)` chunks.
    fn parse_chunks(out: &[u8]) -> Vec<(Vec<u8>, Vec<u8>)> {
        let find = |h: &[u8], n: &[u8]| h.windows(n.len()).position(|w| w == n);
        let mut chunks = Vec::new();
        let mut rest = out;
        while let Some(g) = find(rest, b"\x1b_G") {
            rest = &rest[g + 3..];
            let semi = find(rest, b";").expect("control data ends with ;");
            let control = rest[..semi].to_vec();
            rest = &rest[semi + 1..];
            let st = find(rest, b"\x1b\\").expect("chunk ends with ST");
            chunks.push((control, rest[..st].to_vec()));
            rest = &rest[st + 2..];
        }
        chunks
    }

    #[test]
    fn a_large_image_is_split_into_chunks_bounded_by_the_protocol() {
        // High-entropy pixels keep the PNG from compressing below one 4096-byte
        // chunk, so the chunking path is exercised.
        let mut fb = Framebuffer::new(256, 256, [0, 0, 0, 255]);
        for y in 0..256u32 {
            for x in 0..256u32 {
                let h = x
                    .wrapping_mul(2654435761)
                    .wrapping_add(y.wrapping_mul(2246822519));
                fb.test_and_set(
                    x,
                    y,
                    0.0,
                    [h as u8, (h >> 8) as u8, (h >> 16) as u8, 255],
                    0,
                );
            }
        }

        let expected = base64::encode(&png::to_png_bytes(&fb).unwrap());
        assert!(expected.len() > CHUNK, "test scene too small to chunk");

        let out = encode(&fb).expect("encode");
        let chunks = parse_chunks(&out);
        assert!(chunks.len() > 1, "expected multiple chunks");
        assert!(out.ends_with(b"\x1b\\"), "sequence not terminated");

        let last = chunks.len() - 1;
        let mut reassembled = Vec::new();
        for (i, (control, payload)) in chunks.iter().enumerate() {
            // No chunk exceeds the protocol's payload bound.
            assert!(
                payload.len() <= CHUNK,
                "chunk {i} over the {CHUNK}-byte bound"
            );
            // Only the first chunk carries the image control data.
            if i == 0 {
                assert_eq!(control, b"a=T,f=100,m=1", "first chunk control data wrong");
            } else if i == last {
                assert_eq!(control, b"m=0", "final chunk must end the transmission");
            } else {
                assert_eq!(control, b"m=1", "continuation chunk control data wrong");
            }
            reassembled.extend_from_slice(payload);
        }
        // The chunks concatenate back to exactly the base64 PNG, so no byte is
        // dropped or duplicated at a boundary.
        assert_eq!(reassembled, expected.as_bytes());
    }

    #[test]
    fn empty_framebuffer_payload_is_never_produced() {
        // A zero-area framebuffer still yields a valid (tiny) PNG, so the escape
        // is non-empty; the empty-payload guard is for the internal helper.
        assert!(escape(b"").is_empty());
    }
}
