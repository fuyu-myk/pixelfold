use anyhow::Result;
use crossterm::event::{self, KeyCode, KeyModifiers};
use ratatui::prelude::*;
use ratatui::widgets::{Block, List, ListItem, ListState, Paragraph};
use ratatui::{Frame, layout::Rect};
use std::sync::mpsc::{self, Receiver, TryRecvError};

use crate::search::types::{FoundProtein, PageState};
use crate::search::client::RCSBClient;
use crate::search::download::DownloadManager;

pub mod client;
pub mod download;
pub mod types;


enum SearchMessage {
    Success(Vec<FoundProtein>),
    Error(String),
}

enum DownloadMessage {
    Progress(String, usize),  // (current_file, completed_count)
    Complete(usize),          // total_downloaded
    Error(String),
}

enum AppState {
    Input(SearchInput),
    Searching(SearchInput, Receiver<SearchMessage>),
    Results(SearchInput, SearchSelection),
    Downloading(SearchInput, SearchSelection, DownloadProgress, Receiver<DownloadMessage>),
    Done(String),
}

struct DownloadProgress {
    total: usize,
    completed: usize,
    current: String,
}

impl DownloadProgress {
    fn new(total: usize) -> Self {
        Self {
            total,
            completed: 0,
            current: String::new(),
        }
    }
}

pub fn fetch_structures(terminal: &mut ratatui::Terminal<impl ratatui::backend::Backend>, protein_name: Option<String>) -> Result<()> {
    let mut state = AppState::Input(SearchInput::new());
    
    if let Some(name) = protein_name {
        if let AppState::Input(ref mut input) = state {
            input.query = name;
        }
    }

    let mut redraw_needed = true;
    
    loop {
        match &mut state {
            AppState::Searching(input, rx) => {
                match rx.try_recv() {
                    Ok(SearchMessage::Success(proteins)) => {
                        let selection = SearchSelection::new(proteins);
                        let mut result_input = input.clone();

                        result_input.is_loading = false;
                        state = AppState::Results(result_input, selection);
                        redraw_needed = true;
                    }
                    Ok(SearchMessage::Error(err)) => {
                        let mut err_input = input.clone();

                        err_input.is_loading = false;
                        err_input.err_msg = Some(err);
                        state = AppState::Input(err_input);
                        redraw_needed = true;
                    }
                    Err(TryRecvError::Empty) => {
                        redraw_needed = true;
                    }
                    Err(TryRecvError::Disconnected) => {
                        let mut err_input = input.clone();

                        err_input.is_loading = false;
                        err_input.err_msg = Some("Search thread disconnected".to_string());
                        state = AppState::Input(err_input);
                        redraw_needed = true;
                    }
                }
            }
            AppState::Downloading(_, _, progress, rx) => {
                match rx.try_recv() {
                    Ok(DownloadMessage::Progress(current, completed)) => {
                        progress.current = current;
                        progress.completed = completed;
                        redraw_needed = true;
                    }
                    Ok(DownloadMessage::Complete(total)) => {
                        let message = format!(
                            "Successfully downloaded {} structure(s) to data/\nPress Esc to exit.",
                            total
                        );

                        state = AppState::Done(message);
                        redraw_needed = true;
                    }
                    Ok(DownloadMessage::Error(err)) => {
                        let message = format!(
                            "Download failed: {}\nPress Esc to exit.",
                            err
                        );

                        state = AppState::Done(message);
                        redraw_needed = true;
                    }
                    Err(TryRecvError::Empty) => {
                        redraw_needed = true;
                    }
                    Err(TryRecvError::Disconnected) => {
                        let message = "Download thread disconnected\nPress any key to exit.".to_string();

                        state = AppState::Done(message);
                        redraw_needed = true;
                    }
                }
            }
            _ => {}
        }
        
        if redraw_needed {
            terminal.draw(|frame| ui(frame, &mut state))?;
            redraw_needed = false;
        }

        if event::poll(std::time::Duration::from_millis(16))? {
            match event::read()? {
                event::Event::Key(key) => {
                    let should_quit = handle_input(&mut state, key.code, key.modifiers)?;
                    redraw_needed = true;

                    if should_quit {
                        break;
                    }
                }
                event::Event::Mouse(mouse) => {
                    handle_mouse(&mut state, mouse)?;
                    redraw_needed = true;
                }
                _ => {}
            }
        }
    }

    Ok(())
}

fn ui(frame: &mut Frame, state: &mut AppState) {
    let area = frame.area();
    
    match state {
        AppState::Input(input) => {
            input.render(frame, area);
        }
        AppState::Searching(input, _) => {
            input.render(frame, area);
        }
        AppState::Results(input, selection) => {
            let chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Length(5),
                    Constraint::Min(0),
                    Constraint::Length(2),
                ])
                .split(area);
            
            input.render(frame, chunks[0]);
            selection.render(frame, chunks[1]);

            let selected_count = selection.get_selected_pdb_ids().len();
            if selected_count > 0 {
                let summary = Paragraph::new(format!(
                    "{} structure(s) selected for download.",
                    selected_count
                ))
                .style(Style::default().fg(Color::Green).bold());

                let summary_area = Rect {
                    x: area.x + 1,
                    y: area.y + area.height.saturating_sub(2),
                    width: area.width.saturating_sub(2),
                    height: 1,
                };
                frame.render_widget(summary, summary_area);
            }
        }
        AppState::Downloading(input, selection, progress, _) => {
            let chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Length(5),
                    Constraint::Min(0),
                    Constraint::Length(2),
                ])
                .split(area);
            
            input.render(frame, chunks[0]);
            selection.render(frame, chunks[1]);
            
            let progress_text = format!(
                "Downloading: {} ({}/{})",
                progress.current,
                progress.completed,
                progress.total
            );
            let progress_widget = Paragraph::new(progress_text)
                .style(Style::default().fg(Color::Cyan))
                .block(Block::bordered().title("Download Progress"));
            
            frame.render_widget(progress_widget, chunks[2]);
        }
        AppState::Done(message) => {
            let paragraph = Paragraph::new(message.as_str())
                .style(Style::default().fg(Color::Green))
                .block(
                    Block::bordered()
                        .title("Complete")
                        .title_alignment(Alignment::Center)
                );
            
            frame.render_widget(paragraph, area);
        }
    }
}

fn handle_input(state: &mut AppState, key: KeyCode, modifiers: KeyModifiers) -> Result<bool> {
    match key {
        KeyCode::Esc => {
            if matches!(state, AppState::Downloading(..)) {
                return Ok(false);
            }
            return Ok(true);
        }
        KeyCode::Enter => {
            match state {
                AppState::Input(input) => {
                    if !input.query.is_empty() {
                        let query = input.query.clone();
                        let mut new_input = input.clone();
                        new_input.is_loading = true;
                        
                        let (tx, rx) = mpsc::channel();
                        *state = AppState::Searching(new_input, rx);
                        
                        std::thread::spawn(move || {
                            let rt = tokio::runtime::Runtime::new().unwrap();
                            let client = RCSBClient::new();
                            
                            match rt.block_on(client.search(&query)) {
                                Ok(results) => {
                                    if results.is_empty() {
                                        let _ = tx.send(SearchMessage::Error("No results found".to_string()));
                                    } else {
                                        let proteins: Vec<FoundProtein> = results
                                            .into_iter()
                                            .map(FoundProtein::new)
                                            .collect();
                                        let _ = tx.send(SearchMessage::Success(proteins));
                                    }
                                }
                                Err(e) => {
                                    let _ = tx.send(SearchMessage::Error(format!("Search failed: {}", e)));
                                }
                            }
                        });
                    }
                }
                AppState::Results(input, selection) => {
                    let selected_ids = selection.get_selected_pdb_ids();
                    if !selected_ids.is_empty() {
                        let progress = DownloadProgress::new(selected_ids.len());
                        let input_copy = input.clone();
                        let selection_copy = selection.clone();
                        
                        let (tx, rx) = mpsc::channel();
                        *state = AppState::Downloading(input_copy, selection_copy, progress, rx);
                        
                        let output_dir = std::path::PathBuf::from("data");
                        
                        std::thread::spawn(move || {
                            let rt = tokio::runtime::Runtime::new().unwrap();
                            let manager = DownloadManager::new(output_dir, 3);
                            
                            let mut completed = 0;
                            let total = selected_ids.len();
                            
                            match rt.block_on(manager.download_multiple(
                                selected_ids.clone(),
                                |pdb_id, success| {
                                    if success {
                                        completed += 1;
                                        let _ = tx.send(DownloadMessage::Progress(pdb_id.clone(), completed));
                                    }
                                }
                            )) {
                                Ok(_) => {
                                    let _ = tx.send(DownloadMessage::Complete(total));
                                }
                                Err(e) => {
                                    let _ = tx.send(DownloadMessage::Error(format!("{}", e)));
                                }
                            }
                        });
                    }
                }
                AppState::Done(_) => {
                    return Ok(true);
                }
                _ => {}
            }
        }
        _ => {
            match state {
                AppState::Input(input) => {
                    input.handle_key(key);
                }
                AppState::Results(input, selection) => {
                    if modifiers.contains(KeyModifiers::CONTROL) && key == KeyCode::Char('c') {
                        return Ok(true);
                    }

                    if key == KeyCode::F(1) {
                        input.clear();
                        selection.selected_index = None;

                        for protein in &mut selection.results {
                            protein.selected = false;
                        }

                        *state = AppState::Input(input.clone());

                        return Ok(false);
                    }
                    selection.handle_key(key);
                }
                _ => {}
            }
        }
    }
    
    Ok(false)
}

fn handle_mouse(state: &mut AppState, mouse: event::MouseEvent) -> Result<()> {
    if let AppState::Results(_, selection) = state {
        selection.handle_mouse(mouse);
    }
    Ok(())
}

#[derive(Clone)]
struct SearchInput {
    query: String,
    is_loading: bool,
    err_msg: Option<String>,
}

impl SearchInput {
    fn new() -> Self {
        Self {
            query: String::new(),
            is_loading: false,
            err_msg: None,
        }
    }

    fn handle_key(&mut self, key: KeyCode) {
        match key {
            KeyCode::Char(c) => {
                self.query.push(c);
            }
            KeyCode::F(1) => self.clear(),
            KeyCode::Backspace => {
                self.query.pop();
            }
            _ => {}
        }
        self.err_msg = None;
    }

    fn clear(&mut self) {
        self.query.clear();
        self.err_msg = None;
    }

    fn render(&self, f: &mut Frame, area: Rect) {
        let content = if self.is_loading {
            format!("Searching for '{}'...", self.query)
        } else if self.query.is_empty() {
            "Type to search for protein structures, then press Enter".to_string()
        } else {
            format!("Search: {} (Press Enter to search, F1 to clear)", self.query)
        };

        let color = if self.err_msg.is_some() {
            ratatui::style::Color::Red
        } else if self.is_loading {
            ratatui::style::Color::Yellow
        } else {
            ratatui::style::Color::White
        };

        let paragraph = Paragraph::new(content)
            .style(Style::default().fg(color))
            .alignment(Alignment::Center)
            .block(
                Block::bordered()
                    .title("Protein Structure Search")
                    .title_alignment(Alignment::Center)
                    .border_style(Style::default().fg(color)),
            );

        f.render_widget(paragraph, area);

        if let Some(err) = &self.err_msg {
            let err_para = Paragraph::new(format!("Error: {}", err))
                .style(Style::default().fg(ratatui::style::Color::Red));
            let err_area = Rect {
                x: area.x + 1,
                y: area.y + area.height.saturating_sub(3),
                width: area.width.saturating_sub(2),
                height: 2,
            };

            f.render_widget(err_para, err_area);
        }
    }
}

#[derive(Clone)]
struct SearchSelection {
    results: Vec<FoundProtein>,
    pagination: PageState,
    selected_index: Option<usize>,
    last_list_area: Rect,
}

impl SearchSelection {
    fn new(results: Vec<FoundProtein>) -> Self {
        let total = results.len();

        Self {
            results,
            pagination: PageState::new(total),
            selected_index: Some(0),
            last_list_area: Rect::default(),
        }
    }

    fn current_page_indices(&self, page_size: usize) -> Vec<usize> {
        let start = self.pagination.current_page * page_size;
        let end = (start + page_size).min(self.results.len());

        (start..end).collect()
    }

    fn update_page_size(&mut self, page_size: usize) {
        let old_page_size = self.pagination.items_per_page;
        if old_page_size != page_size {
            let current_item = self.selected_index;
            self.pagination.items_per_page = page_size;
            self.pagination.total_pages = (self.results.len() + page_size - 1) / page_size;

            if let Some(idx) = current_item {
                self.pagination.current_page = idx / page_size;
            }
        }
    }

    fn get_display_index(&self, page_size: usize) -> Option<usize> {
        self.selected_index.map(|idx| idx % page_size)
    }

    fn handle_key(&mut self, key: KeyCode) {
        match key {
            KeyCode::Left => {
                self.selected_index = None;
                self.pagination.previous_page();
            }
            KeyCode::Right => {
                self.selected_index = None;
                self.pagination.next_page();
            }
            _ => {}
        }
    }

    fn handle_mouse(&mut self, event: event::MouseEvent) {
        if event.kind == event::MouseEventKind::Down(event::MouseButton::Left) {
            if let Some(idx) = self.get_clicked_index(event.column, event.row) {
                if idx < self.results.len() {
                    self.selected_index = Some(idx);
                    self.toggle_selection(idx);
                }
            }
        }
    }

    fn get_clicked_index(&self, _x: u16, y: u16) -> Option<usize> {
        let area = self.last_list_area;
        
        if y < area.y || y >= area.y + area.height {
            return None;
        }
        
        let relative_y = y - area.y;
        
        if relative_y < 1 || relative_y >= area.height.saturating_sub(1) {
            return None;
        }
        
        let available_height = area.height.saturating_sub(2) as usize;
        let lines_per_item = (available_height / 10).max(1);
        let display_idx = ((relative_y - 1) as usize) / lines_per_item;
        
        if display_idx >= 10 {
            return None;
        }
        
        let page_start = self.pagination.current_page * self.pagination.items_per_page;
        let clicked_idx = page_start + display_idx;
        
        if clicked_idx < self.results.len() {
            Some(clicked_idx)
        } else {
            None
        }
    }

    fn toggle_selection(&mut self, idx: usize) {
        if idx < self.results.len() {
            self.results[idx].selected = !self.results[idx].selected;
        }
    }

    fn get_selected_pdb_ids(&self) -> Vec<String> {
        self.results
            .iter()
            .filter(|protein| protein.selected)
            .map(|protein| protein.pdb_id.clone())
            .collect()
    }

    fn render(&mut self, f: &mut Frame, area: Rect) {
        self.last_list_area = area;
        
        let page_size = 10;
        self.update_page_size(page_size);
        
        let available_height = area.height.saturating_sub(2) as usize;
        let lines_per_item = (available_height / page_size).max(1);
        
        let indices = self.current_page_indices(page_size);

        let items: Vec<ListItem> = indices
            .iter()
            .map(|&global_idx| {
                let result = &self.results[global_idx];
                let is_selected = result.selected;

                let text = result.fmt_for_display();

                let style = if is_selected {
                    Style::default()
                        .bg(Color::Green)
                        .fg(Color::Black)
                } else {
                    Style::default()
                };

                let mut lines: Vec<Line> = vec![Line::from(text)];

                for _ in 1..lines_per_item {
                    lines.push(Line::from(""));
                }
                
                ListItem::new(lines).style(style)
            })
            .collect();

        let mut list_state = ListState::default();
        list_state.select(self.get_display_index(page_size));

        let list = List::new(items).block(
            Block::bordered()
                .title(format!(
                    "Search results (Page {}/{}) | Click: Select/Deselect, ← →: Navigate Pages | Enter: Download selected | Esc: Quit",
                    self.pagination.current_page + 1,
                    self.pagination.total_pages,
                ))
                .title_alignment(Alignment::Center),
        );
        f.render_stateful_widget(list, area, &mut list_state);
    }
}
