use serde::{Deserialize, Serialize};


/// Response structure from RCSB search API
#[derive(Debug, Deserialize)]
pub struct SearchResponse {
    #[serde(default)]
    pub result_set: Vec<SearchResult>,
    pub total_count: usize,
}

/// Individual search result from RCSB
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SearchResult {
    pub identifier: String,
    pub score: f64,
}

#[derive(Debug, Deserialize, Clone)]
pub struct FoundProtein {
    pub pdb_id: String,
    pub title: String,
    pub resolution: Option<f64>,
    pub date: String,
    pub selected: bool,
}

impl FoundProtein {
    pub fn new(result: SearchResult) -> Self {
        Self {
            pdb_id: result.identifier,
            title: String::new(),
            resolution: None,
            date: String::new(),
            selected: false,
        }
    }

    pub fn fmt_for_display(&self) -> String {
        if self.title.is_empty() {
            format!("{}", self.pdb_id)
        } else {
            let res = self.resolution
                .map(|r| format!(" - {:.2}Å", r))
                .unwrap_or_else(|| " - N/A".to_string());

            format!(
                "{}: {}{} ({})",
                self.pdb_id,
                self.title,
                res,
                self.date
            )
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
        let total_pages = (total + 9) / 10;

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

    pub fn get_range(&self) -> (usize, usize) {
        let start = self.current_page * self.items_per_page;
        let end = ((self.current_page + 1) * self.items_per_page)
            .min(self.current_page * self.items_per_page + self.items_per_page);

        (start, end)
    }
}
