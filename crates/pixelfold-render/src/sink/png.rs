//! Framebuffer to PNG. The headless output of `pixelfold render`, and the basis
//! of the renderer's snapshot tests.
//!
//! A later terminal sink reuses [`to_png_bytes`] as the compressed Kitty payload.

use std::io::{self, Write};
use std::path::Path;

use crate::framebuffer::Framebuffer;

/// Re-exported so sibling sink modules can name the PNG error without reaching
/// past the local `png` module to the external crate of the same name.
pub use png::EncodingError;

/// Encode a framebuffer as an RGBA8 PNG to `writer`, at the default compression.
pub fn encode<W: Write>(fb: &Framebuffer, writer: W) -> Result<(), png::EncodingError> {
    encode_with(fb, writer, false)
}

fn encode_with<W: Write>(
    fb: &Framebuffer,
    writer: W,
    fast: bool,
) -> Result<(), png::EncodingError> {
    let mut encoder = png::Encoder::new(writer, fb.width(), fb.height());
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);

    if fast {
        encoder.set_compression(png::Compression::Fast);
    }

    let mut writer = encoder.write_header()?;
    writer.write_image_data(fb.color().as_flattened())?;

    Ok(())
}

/// Encode a framebuffer as a PNG byte buffer at the default compression.
pub fn to_png_bytes(fb: &Framebuffer) -> Result<Vec<u8>, png::EncodingError> {
    let mut bytes = Vec::new();
    encode(fb, &mut bytes)?;
    Ok(bytes)
}

/// Encode a framebuffer as a PNG byte buffer at the fastest compression.
pub fn to_png_bytes_fast(fb: &Framebuffer) -> Result<Vec<u8>, png::EncodingError> {
    let mut bytes = Vec::new();
    encode_with(fb, &mut bytes, true)?;
    Ok(bytes)
}

/// Write a framebuffer to `path` as a PNG.
pub fn write_png(fb: &Framebuffer, path: &Path) -> Result<(), png::EncodingError> {
    let file = std::fs::File::create(path)?;
    encode(fb, std::io::BufWriter::new(file))
}

/// Decode an RGBA8 PNG into width, height, and row-major pixels. Only the RGBA8
/// PNGs this module writes are supported; anything else is an error.
pub fn decode_rgba(bytes: &[u8]) -> Result<(u32, u32, Vec<[u8; 4]>), png::DecodingError> {
    let decoder = png::Decoder::new(io::Cursor::new(bytes));
    let mut reader = decoder.read_info()?;
    let size = reader.output_buffer_size().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "PNG output buffer size overflow",
        )
    })?;
    let mut buf = vec![0u8; size];
    let info = reader.next_frame(&mut buf)?;

    if info.color_type != png::ColorType::Rgba || info.bit_depth != png::BitDepth::Eight {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "expected an RGBA8 PNG").into());
    }

    buf.truncate(info.buffer_size());
    let pixels = buf
        .chunks_exact(4)
        .map(|c| [c[0], c[1], c[2], c[3]])
        .collect();

    Ok((info.width, info.height, pixels))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn png_roundtrips_the_framebuffer() {
        let mut fb = Framebuffer::new(5, 3, [10, 20, 30, 255]);
        fb.test_and_set(0, 0, 0.0, [200, 100, 50, 255], 1);
        fb.test_and_set(4, 2, 0.0, [1, 2, 3, 255], 2);

        let bytes = to_png_bytes(&fb).expect("encode");
        let (w, h, pixels) = decode_rgba(&bytes).expect("decode");

        assert_eq!((w, h), (5, 3));
        assert_eq!(pixels.len(), 15);
        assert_eq!(pixels, fb.color().to_vec());
    }

    #[test]
    fn decode_rejects_non_png_input() {
        assert!(decode_rgba(b"not a png").is_err());
    }
}
