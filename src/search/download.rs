use std::{fs, io::{Read, Write}, path::PathBuf};

use anyhow::Result;
use flate2::read::GzDecoder;
use futures::stream::StreamExt;

use crate::search::client::RCSBClient;


pub struct DownloadManager {
    rscb_client: RCSBClient,
    output_dir: PathBuf,
    concurrency: usize,
}

impl DownloadManager {
    pub fn new(output_dir: PathBuf, concurrency: usize) -> Self {
        Self {
            rscb_client: RCSBClient::new(),
            output_dir,
            concurrency,
        }
    }

    pub async fn download_multiple(
        &self,
        pdb_ids: Vec<String>,
        mut progress_callback: impl FnMut(String, bool) + Send
    ) -> Result<()> {
        fs::create_dir_all(&self.output_dir)?;

        let results = tokio_stream::iter(pdb_ids)
            .map(|id| {
                let client = self.rscb_client.clone();
                let output_dir = self.output_dir.clone();

                async move {
                    Self::download_single(&client, &id, &output_dir).await
                }
            })
            .buffer_unordered(self.concurrency)
            .collect::<Vec<Result<String>>>()
            .await;

        for result in results {
            match result {
                Ok(pdb_id) => progress_callback(pdb_id, true),
                Err(e) => {
                    eprintln!("Download error: {}", e);
                }
            }
        }

        Ok(())
    }

    pub async fn download_single(
        client: &RCSBClient,
        pdb_id: &str,
        output_dir: &PathBuf,
    ) -> Result<String> {
        let bytes = client.download_cif(pdb_id).await?;
        let decompressed = Self::decompress(&bytes)?;

        let file_pth = output_dir.join(format!("{}.cif", pdb_id.to_uppercase()));
        let mut file = fs::File::create(&file_pth)?;
        file.write_all(&decompressed)?;

        Ok(pdb_id.to_uppercase())
    }

    /// Decompress gzipped data to raw bytes
    fn decompress(data: &[u8]) -> Result<Vec<u8>> {
        let mut decoder = GzDecoder::new(data);
        let mut decompressed = Vec::new();
        decoder.read_to_end(&mut decompressed)?;

        Ok(decompressed)
    }
}

impl Clone for DownloadManager {
    fn clone(&self) -> Self {
        Self {
            rscb_client: self.rscb_client.clone(),
            output_dir: self.output_dir.clone(),
            concurrency: self.concurrency,
        }
    }
}
