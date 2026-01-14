use anyhow::Result;
use reqwest::Client;

use crate::search::types::{SearchResponse, SearchResult};


#[derive(Debug, Clone)]
pub struct RCSBClient {
    client: reqwest::Client,
    url: String,
}

impl RCSBClient {
    pub fn new() -> Self {
        Self {
            client: Client::new(),
            url: "https://search.rcsb.org/rcsbsearch/v2/query".to_string(),
        }
    }

    pub async fn search(&self, query: &str) -> Result<Vec<SearchResult>> {
        let query = build_query(query);
        let res = self.client
            .post(&self.url)
            .json(&query)
            .timeout(std::time::Duration::from_secs(30))
            .send()
            .await?;

        if !res.status().is_success() {
            anyhow::bail!(
                "RCSB API error: {} - {}",
                res.status(),
                res.text().await?
            );
        }

        let data: SearchResponse = res.json().await?;
        Ok(data.result_set)
    }

    pub async fn download_cif(&self, pdb_id: &str) -> Result<bytes::Bytes> {
        let url = format!("https://files.wwpdb.org/download/{}.cif.gz", pdb_id.to_uppercase());

        let res = self.client
            .get(&url)
            .timeout(std::time::Duration::from_secs(60))
            .send()
            .await?;

        if !res.status().is_success() {
            anyhow::bail!(
                "Failed to download CIF for {}: {}",
                pdb_id,
                res.status(),
            );
        }

        Ok(res.bytes().await?)
    }
}

impl Default for RCSBClient {
    fn default() -> Self {
        Self::new()
    }
}

fn build_query(query: &str) -> serde_json::Value {
    serde_json::json!({
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": query
            }
        },
        "return_type": "entry",
        "request_options": {
            "sort": [{
                "sort_by": "score",
                "direction": "desc"
            }],
            "paginate": {
                "start": 0,
                "rows": 100
            }
        }
    })
}
