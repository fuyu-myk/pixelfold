use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

/// The outcome of resolving a structure argument.
#[derive(Debug, PartialEq)]
pub enum Resolution {
    /// A local file to open directly.
    Found(PathBuf),
    /// A PDB id that is not cached; the caller should fetch it to `dest`.
    Fetch { id: String, dest: PathBuf },
}

/// Resolve a structure argument.
///
/// Resolution order:
/// 1. an existing path (absolute or relative) is opened as given
/// 2. a 4-character PDB id resolves to `<cache_dir>/<ID>.cif`, either found in
///    the cache or requested for fetch
/// 3. otherwise an error explaining what was expected
pub fn resolve(input: &str, cache_dir: &Path) -> Result<Resolution> {
    let candidate = Path::new(input);
    if candidate.exists() {
        return Ok(Resolution::Found(candidate.to_path_buf()));
    }

    if is_pdb_id(input) {
        let id = input.to_uppercase();
        let dest = cache_dir.join(format!("{id}.cif"));
        if dest.exists() {
            return Ok(Resolution::Found(dest));
        }

        return Ok(Resolution::Fetch { id, dest });
    }

    bail!("'{input}' is not an existing file and not a 4-character PDB id")
}

/// A PDB id is four characters: a leading digit 1-9 followed by three
/// alphanumeric characters (for example 4HHB or 1CRN).
pub fn is_pdb_id(s: &str) -> bool {
    let bytes = s.as_bytes();
    bytes.len() == 4
        && bytes[0].is_ascii_digit()
        && bytes[0] != b'0'
        && bytes[1..].iter().all(u8::is_ascii_alphanumeric)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn recognizes_valid_pdb_ids() {
        assert!(is_pdb_id("4HHB"));
        assert!(is_pdb_id("1crn"));
        assert!(is_pdb_id("6ZzR"));
    }

    #[test]
    fn rejects_non_pdb_ids() {
        assert!(!is_pdb_id("prot")); // no leading digit
        assert!(!is_pdb_id("0ABC")); // leading zero
        assert!(!is_pdb_id("12345")); // too long
        assert!(!is_pdb_id("1CR")); // too short
        assert!(!is_pdb_id("1-AB")); // non-alphanumeric
    }

    #[test]
    fn existing_path_used_as_given() {
        let dir = tempfile::tempdir().unwrap();
        let file = dir.path().join("prot.cif");
        fs::write(&file, "data").unwrap();
        let cache = tempfile::tempdir().unwrap();

        let resolved = resolve(file.to_str().unwrap(), cache.path()).unwrap();
        assert_eq!(resolved, Resolution::Found(file));
    }

    #[test]
    fn pdb_id_resolves_to_cached_file() {
        let cache = tempfile::tempdir().unwrap();
        let cached = cache.path().join("4HHB.cif");
        fs::write(&cached, "data").unwrap();

        let resolved = resolve("4hhb", cache.path()).unwrap();
        assert_eq!(resolved, Resolution::Found(cached));
    }

    #[test]
    fn uncached_pdb_id_requests_fetch() {
        let cache = tempfile::tempdir().unwrap();
        let resolved = resolve("4hhb", cache.path()).unwrap();
        assert_eq!(
            resolved,
            Resolution::Fetch {
                id: "4HHB".to_string(),
                dest: cache.path().join("4HHB.cif"),
            }
        );
    }

    #[test]
    fn unknown_input_errors() {
        let cache = tempfile::tempdir().unwrap();
        let err = resolve("not_a_file", cache.path()).unwrap_err().to_string();
        assert!(err.contains("not an existing file"), "got: {err}");
    }
}
