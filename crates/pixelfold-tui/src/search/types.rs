use serde::Deserialize;

use pixelfold_fetch::types::SearchResult;

#[derive(Debug, Deserialize, Clone)]
pub struct FoundProtein {
    pub pdb_id: String,
    pub title: String,
    pub date: String,
    pub selected: bool,
}

impl FoundProtein {
    pub fn new(result: SearchResult) -> Self {
        Self {
            pdb_id: result.identifier,
            title: String::new(),
            date: String::new(),
            selected: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct PageState {
    pub current_page: usize,
    pub total_pages: usize,
    pub items_per_page: usize,
}

impl PageState {
    pub fn new(total: usize) -> Self {
        let total_pages = total.div_ceil(10);

        Self {
            current_page: 0,
            total_pages: total_pages.max(1),
            items_per_page: 10,
        }
    }

    pub fn next_page(&mut self) {
        if self.current_page < self.total_pages - 1 {
            self.current_page += 1;
        } else {
            self.current_page = 0;
        }
    }

    pub fn previous_page(&mut self) {
        if self.current_page > 0 {
            self.current_page -= 1;
        } else {
            self.current_page = self.total_pages - 1;
        }
    }
}
