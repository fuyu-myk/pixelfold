//! Pixelfold fetch: RCSB and AlphaFold structure search and download.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::Duration;

use anyhow::{Context, Result, anyhow, bail};
use flate2::read::GzDecoder;
use serde::Deserialize;

pub mod client;
pub mod download;
pub mod types;

/// Download a structure's mmCIF and write it to `dest_dir/<ID>.cif`.
///
/// Ids prefixed `AF-` (for example `AF-P01542`) are AlphaFold DB UniProt
/// accessions, fetched from AlphaFold DB; every other id is a PDB id fetched
/// from the wwPDB. Blocking; returns the path to the written file.
pub fn fetch_cif(pdb_id: &str, dest_dir: &Path) -> Result<PathBuf> {
    let id = pdb_id.to_uppercase();
    let cif = match id.strip_prefix("AF-") {
        Some(accession) => fetch_alphafold_cif(accession)?,
        None => fetch_pdb_cif(&id)?,
    };

    std::fs::create_dir_all(dest_dir)?;
    let dest = dest_dir.join(format!("{id}.cif"));
    std::fs::write(&dest, cif)?;

    Ok(dest)
}

/// Download the Chemical Component Dictionary entry of component `code` from
/// RCSB and write it to `dest_dir/<CODE>.cif`, or `None` when the dictionary
/// has no such component.
///
/// A code is one to five letters or digits; anything else (a co-folding
/// tool's `LIG_A`, a path) is refused before any request. Blocking; returns
/// the path to the written file.
pub fn fetch_ccd(code: &str, dest_dir: &Path) -> Result<Option<PathBuf>> {
    let code = code.to_uppercase();
    if !is_component_code(&code) {
        bail!("{code} is not a chemical component code");
    }
    let url = format!("https://files.rcsb.org/ligands/download/{code}.cif");
    let response = http_client()?.get(&url).send()?;
    if response.status() == reqwest::StatusCode::NOT_FOUND {
        return Ok(None);
    }
    if !response.status().is_success() {
        bail!(
            "failed to download the dictionary entry {code} from RCSB (HTTP {})",
            response.status()
        );
    }

    let text = response.bytes()?;
    std::fs::create_dir_all(dest_dir)?;
    let dest = dest_dir.join(format!("{code}.cif"));
    // Written aside and renamed into place, so an interrupted download never
    // leaves a cached entry short of its bonds.
    let partial = dest_dir.join(format!("{code}.cif.part"));
    std::fs::write(&partial, text)?;
    std::fs::rename(&partial, &dest)?;

    Ok(Some(dest))
}

/// Whether `code` can name a chemical component: one to five ASCII letters
/// or digits.
pub fn is_component_code(code: &str) -> bool {
    (1..=5).contains(&code.len()) && code.chars().all(|c| c.is_ascii_alphanumeric())
}

/// A blocking HTTP client with a descriptive User-Agent and request timeout.
/// Some hosts (for example the AlphaFold DB API) reject requests that send no
/// User-Agent with `403 Forbidden`.
fn http_client() -> Result<reqwest::blocking::Client> {
    reqwest::blocking::Client::builder()
        .user_agent(concat!("pixelfold/", env!("CARGO_PKG_VERSION")))
        .timeout(Duration::from_secs(60))
        .build()
        .context("building HTTP client")
}

/// Fetch a gzipped mmCIF for a PDB id from the wwPDB and decompress it.
fn fetch_pdb_cif(id: &str) -> Result<Vec<u8>> {
    let url = format!("https://files.wwpdb.org/download/{id}.cif.gz");
    let response = http_client()?.get(&url).send()?;
    if !response.status().is_success() {
        bail!(
            "failed to download {id} from wwPDB (HTTP {})",
            response.status()
        );
    }

    let gz = response.bytes()?;
    let mut cif = Vec::new();
    GzDecoder::new(&gz[..]).read_to_end(&mut cif)?;

    Ok(cif)
}

/// Fetch an AlphaFold DB ModelCIF for a UniProt accession. The file URL and
/// its version is resolved through the AlphaFold DB API. The model file
/// itself is not compressed.
fn fetch_alphafold_cif(accession: &str) -> Result<Vec<u8>> {
    let client = http_client()?;
    let api = format!("https://alphafold.ebi.ac.uk/api/prediction/{accession}");
    let meta = client.get(&api).send()?;
    if !meta.status().is_success() {
        bail!(
            "failed to resolve AlphaFold prediction for {accession} (HTTP {})",
            meta.status()
        );
    }

    let predictions: Vec<AlphaFoldPrediction> = serde_json::from_slice(&meta.bytes()?)
        .with_context(|| format!("parsing AlphaFold API response for {accession}"))?;
    let cif_url = predictions
        .into_iter()
        .next()
        .ok_or_else(|| anyhow!("no AlphaFold prediction for {accession}"))?
        .cif_url;

    let response = client.get(&cif_url).send()?;
    if !response.status().is_success() {
        bail!(
            "failed to download AlphaFold model {cif_url} (HTTP {})",
            response.status()
        );
    }

    Ok(response.bytes()?.to_vec())
}

/// One prediction record from the AlphaFold DB API.
#[derive(Deserialize)]
struct AlphaFoldPrediction {
    #[serde(rename = "cifUrl")]
    cif_url: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_component_code_is_one_to_five_letters_or_digits() {
        for code in ["A", "STI", "0G6", "A1AAT"] {
            assert!(is_component_code(code), "{code}");
        }
        for code in ["", "LIG_A", "ABCDEF", "../x", "S I"] {
            assert!(!is_component_code(code), "{code}");
        }
    }
}
