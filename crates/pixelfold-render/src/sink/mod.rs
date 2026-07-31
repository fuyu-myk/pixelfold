//! Framebuffer sinks: PNG for headless output and snapshot tests, and terminal
//! graphics protocols for showing pixels in a terminal.

pub mod base64;
pub mod detect;
pub mod halfblock;
pub mod iterm2;
pub mod kitty;
pub mod png;

pub use detect::{Protocol, detect, is_ssh};
pub use png::{decode_rgba, encode, to_png_bytes, to_png_bytes_fast, write_png};

use crate::framebuffer::Framebuffer;

/// Encode a framebuffer for a terminal protocol, ready to write to the terminal.
pub fn to_terminal(fb: &Framebuffer, protocol: Protocol) -> Result<Vec<u8>, png::EncodingError> {
    match protocol {
        Protocol::Kitty => kitty::encode(fb),
        Protocol::Iterm2 => iterm2::encode(fb),
        Protocol::HalfBlock => Ok(halfblock::encode(fb).into_bytes()),
    }
}
