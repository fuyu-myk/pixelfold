//! Pixelfold rendering: scene to RGBA framebuffer to terminal and PNG sinks.

#![forbid(unsafe_code)]

pub mod draw;
pub mod framebuffer;
pub mod overlay;
pub mod raster;
pub mod renderer;
pub mod scene;
pub mod sink;

pub use draw::{
    bfactor_to_color, draw_dashed_line, draw_line, element_to_color, hbond_energy_to_color,
};
pub use framebuffer::{Framebuffer, NO_ID};
pub use overlay::draw_segment;
pub use raster::{RenderOptions, rasterize};
pub use scene::{Coloring, Scene};
pub use sink::{decode_rgba, to_png_bytes, write_png};
