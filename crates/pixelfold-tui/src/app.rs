use anyhow::Result;
use crossterm::event::{self, KeyCode};
use ratatui_image::picker::Picker;

use crate::App;
use crate::inputs::{handle_input, handle_mouse};
use crate::ui::ui;

pub(crate) fn run_app(
    terminal: &mut ratatui::Terminal<impl ratatui::backend::Backend>,
    app: &mut App,
    picker: &Picker,
) -> Result<()> {
    loop {
        // Install a finished network build.
        if let Some(rx) = &app.network_build {
            match rx.try_recv() {
                Ok(view) => {
                    app.network_view = Some(view);
                    app.network_build = None;
                    app.redraw_needed = true;
                }
                Err(std::sync::mpsc::TryRecvError::Empty) => {}
                Err(std::sync::mpsc::TryRecvError::Disconnected) => {
                    // The worker died before sending; drop the build so a later
                    // toggle can start a new one.
                    app.network_build = None;
                }
            }
        }

        if app.redraw_needed {
            terminal.draw(|frame| ui(frame, app, picker))?;
            app.redraw_needed = false;
        }

        if !event::poll(std::time::Duration::from_millis(16))? {
            continue;
        }

        // Drain every buffered event before drawing again.
        let mut quit = false;
        loop {
            match event::read()? {
                event::Event::Key(key) => {
                    handle_input(app, key.code, key.modifiers)?;
                    if matches!(key.code, KeyCode::Char('q')) {
                        quit = true;
                        break;
                    }
                }
                event::Event::Mouse(mouse) => {
                    handle_mouse(app, mouse)?;
                }
                event::Event::Resize(_, _) => {
                    app.redraw_needed = true;
                }
                _ => {}
            }

            if !event::poll(std::time::Duration::ZERO)? {
                break;
            }
        }

        if quit {
            break;
        }
    }

    Ok(())
}
