use anyhow::Result;
use crossterm::event::{self, KeyCode};

use crate::App;
use crate::inputs::{handle_input, handle_mouse};
use crate::ui::ui;

pub(crate) fn run_app(
    terminal: &mut ratatui::Terminal<impl ratatui::backend::Backend>,
    app: &mut App,
) -> Result<()> {
    loop {
        if app.redraw_needed {
            terminal.draw(|frame| ui(frame, app))?;
            app.redraw_needed = false;
        }

        if event::poll(std::time::Duration::from_millis(16))? {
            match event::read()? {
                event::Event::Key(key) => {
                    handle_input(app, key.code, key.modifiers)?;

                    if matches!(key.code, KeyCode::Char('q')) {
                        break;
                    }
                }
                event::Event::Mouse(mouse) => {
                    handle_mouse(app, mouse)?;
                }
                event::Event::Resize(_, _) => {
                    // Dimensions changed: drop the projection cache (keyed to the old
                    // size) and repaint.
                    app.projected_atom_cache = None;
                    app.redraw_needed = true;
                }
                _ => {}
            }
        }
    }

    Ok(())
}
