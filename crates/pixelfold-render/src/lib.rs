//! Pixelfold rendering: scene to RGBA framebuffer to terminal and PNG sinks.

pub mod draw;
pub mod renderer;

pub use draw::{bfactor_to_color, draw_dashed_line, draw_line, hbond_energy_to_color};
