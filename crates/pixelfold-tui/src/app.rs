use anyhow::Result;
use crossterm::event::{self, KeyCode};
use ratatui::layout::Rect;

use crate::App;
use crate::inputs::{handle_input, handle_mouse};
use crate::ui::ui;

pub(crate) fn run_app(
    terminal: &mut ratatui::Terminal<impl ratatui::backend::Backend>,
    app: &mut App,
) -> Result<()> {
    loop {
        let size = terminal.size()?;
        let canvas_width = size.width as f32 * 2.0;
        let canvas_height = size.height as f32 * 4.0;

        if app.redraw_needed {
            terminal.draw(|frame| ui(frame, app))?;
            app.redraw_needed = false;
        }

        if event::poll(std::time::Duration::from_millis(16))? {
            match event::read()? {
                event::Event::Key(key) => {
                    handle_input(app, key.code, key.modifiers, canvas_width, canvas_height)?;

                    if matches!(key.code, KeyCode::Char('q')) {
                        break;
                    }
                }
                event::Event::Mouse(mouse) => {
                    let area = Rect {
                        x: 0,
                        y: 0,
                        width: size.width,
                        height: size.height,
                    };
                    handle_mouse(app, mouse, area)?;
                }
                _ => {}
            }
        }
    }

    Ok(())
}
