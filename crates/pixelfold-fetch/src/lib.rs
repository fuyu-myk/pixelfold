//! Pixelfold fetch: RCSB structure search and download.

use std::io::Read;
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use flate2::read::GzDecoder;

pub mod client;
pub mod download;
pub mod types;

/// Download a structure's mmCIF from RCSB and write it to `dest_dir/<ID>.cif`.
///
/// Blocking; returns the path to the written file.
pub fn fetch_cif(pdb_id: &str, dest_dir: &Path) -> Result<PathBuf> {
    let id = pdb_id.to_uppercase();
    let url = format!("https://files.wwpdb.org/download/{id}.cif.gz");

    let response = reqwest::blocking::Client::new()
        .get(&url)
        .timeout(std::time::Duration::from_secs(60))
        .send()?;
    if !response.status().is_success() {
        bail!(
            "failed to download {id} from RCSB (HTTP {})",
            response.status()
        );
    }

    let gz = response.bytes()?;
    let mut cif = Vec::new();
    GzDecoder::new(&gz[..]).read_to_end(&mut cif)?;

    std::fs::create_dir_all(dest_dir)?;
    let dest = dest_dir.join(format!("{id}.cif"));
    std::fs::write(&dest, cif)?;

    Ok(dest)
}
