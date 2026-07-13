use serde::{Deserialize, Serialize};

/// Response structure from RCSB search API
#[derive(Debug, Deserialize)]
pub struct SearchResponse {
    #[serde(default)]
    pub result_set: Vec<SearchResult>,
    pub total_count: usize,
}

/// Detailed entry data fetched from RCSB
#[derive(Debug, Deserialize)]
pub struct SearchData {
    #[serde(rename = "struct")]
    pub struct_info: StructInfo,
    pub rcsb_accession_info: AccessionInfo,
}

#[derive(Debug, Deserialize)]
pub struct StructInfo {
    pub title: String,
}

#[derive(Debug, Deserialize)]
pub struct AccessionInfo {
    pub revision_date: String,
}

/// Individual search result from RCSB
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SearchResult {
    pub identifier: String,
    pub score: f64,
}
