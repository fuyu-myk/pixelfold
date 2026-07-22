//! The headless commands over real structures, through the built binary.
//!
//! The structures come from the benchmark cache; it skips rather than fails
//! if the structures are not available.

use std::path::{Path, PathBuf};
use std::process::Command;

fn cached(id: &str) -> Option<PathBuf> {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(format!("../../benchmarks/cache/{id}.cif"))
        .canonicalize()
        .ok()
}

struct Run {
    stdout: String,
    stderr: String,
    ok: bool,
}

impl Run {
    fn rows(&self) -> Vec<&str> {
        self.stdout.lines().skip(1).collect()
    }
}

/// Run the binary against a cached structure, or nothing when it is absent.
fn pixelfold(id: &str, args: &[&str]) -> Option<Run> {
    let structure = cached(id)?;
    let output = Command::new(env!("CARGO_BIN_EXE_pixelfold"))
        .args(args)
        .arg(structure)
        .output()
        .expect("the binary runs");

    Some(Run {
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        ok: output.status.success(),
    })
}

macro_rules! run {
    ($id:expr, $args:expr) => {
        match pixelfold($id, $args) {
            Some(run) => run,
            None => return,
        }
    };
}

#[test]
fn secondary_structure_reports_one_row_per_residue() {
    let run = run!("1CRN", &["ss", "--format", "tsv"]);
    assert!(run.ok, "command failed: {}", run.stderr);

    assert_eq!(
        run.stdout.lines().next(),
        Some("chain\tresi\ticode\tresn\tss")
    );
    // Crambin is 46 residues, plus whatever solvent the file carries.
    assert!(run.rows().len() >= 46, "expected a row per residue");
    assert!(run.rows().iter().any(|row| row.starts_with("A\t46\t")));
    assert!(
        run.rows().iter().all(|line| line.split('\t').count() == 5),
        "every row has every column"
    );
}

#[test]
fn a_selection_narrows_the_report() {
    let run = run!("1CRN", &["sasa", "--format", "tsv", "--select", "resi 1-5"]);
    assert!(run.ok, "command failed: {}", run.stderr);

    assert_eq!(run.rows().len(), 5);
    assert!(run.rows()[0].starts_with("A\t1\t"));
}

#[test]
fn interactions_are_json_parseable_and_typed() {
    let run = run!(
        "1CRN",
        &["interactions", "--format", "json", "--type", "disulfide"]
    );
    assert!(run.ok, "command failed: {}", run.stderr);

    let parsed: serde_json::Value = serde_json::from_str(&run.stdout).expect("valid json");
    let found = parsed.as_array().expect("an array");
    // Crambin has three disulfide bridges and the filter admits nothing else.
    assert_eq!(found.len(), 3);
    assert!(found.iter().all(|row| row["kind"] == "disulfide"));
    assert!(found.iter().all(|row| row["atoms_a"] == "SG"));
}

#[test]
fn a_selection_matching_nothing_is_an_error_rather_than_an_empty_report() {
    let run = run!("1CRN", &["sasa", "--select", "chain ZZZ"]);

    assert!(!run.ok, "expected a failure, got: {}", run.stdout);
    assert!(
        run.stderr.contains("matched no atoms"),
        "unhelpful error: {}",
        run.stderr
    );
}

#[test]
fn a_malformed_selection_reports_the_query() {
    let run = run!("1CRN", &["ss", "--select", "resi and"]);

    assert!(!run.ok);
    assert!(
        run.stderr.contains("selection"),
        "unhelpful error: {}",
        run.stderr
    );
}

/// Ordered waters sit on the surface being measured, so including them occludes it.
#[test]
fn surface_area_is_of_the_polymer_alone() {
    let run = run!("1UBQ", &["sasa", "--format", "json"]);
    assert!(run.ok, "command failed: {}", run.stderr);

    let parsed: serde_json::Value = serde_json::from_str(&run.stdout).expect("valid json");
    let rows = parsed.as_array().expect("an array");

    assert!(
        rows.iter().all(|row| row["resn"] != "HOH"),
        "solvent is not part of the surface"
    );
    // Ubiquitin is 76 residues; the file also carries 58 waters.
    assert_eq!(rows.len(), 76);

    // Glu18 is exposed. With the waters counted as structure it reads 41.
    let glu18 = rows
        .iter()
        .find(|row| row["resi"] == 18)
        .expect("residue 18");
    let area = glu18["sasa"].as_f64().expect("a number");
    assert!(
        (95.0..110.0).contains(&area),
        "expected the exposed area, got {area}"
    );
}

/// A chymotrypsin-numbered structure distinguishes residues only by insertion
/// code. Dropping it emits byte-identical rows for different residues.
#[test]
fn insertion_codes_keep_residues_distinct() {
    let run = run!("1PPB", &["ss", "--format", "tsv"]);
    assert!(run.ok, "command failed: {}", run.stderr);

    let mut rows = run.rows();
    let total = rows.len();
    rows.sort_unstable();
    rows.dedup();
    assert_eq!(total, rows.len(), "duplicate rows: residues are not keyed");

    assert!(
        run.rows().iter().any(|row| {
            let mut cells = row.split('\t');
            cells.nth(2).is_some_and(|icode| !icode.is_empty())
        }),
        "no insertion code survived into the report"
    );
}

/// `pixelfold ... | head` is a normal thing to do, in every format.
#[test]
fn a_closed_pipe_is_not_a_failure_in_any_format() {
    let Some(structure) = cached("1CRN") else {
        return;
    };

    for format in ["table", "tsv", "json"] {
        let status = Command::new("sh")
            .arg("-c")
            .arg(format!(
                "set -o pipefail; {} interactions {} --format {format} | head -1 > /dev/null",
                env!("CARGO_BIN_EXE_pixelfold"),
                structure.display(),
            ))
            .status()
            .expect("the shell runs");

        assert!(
            status.success(),
            "{format} reported a closed pipe as an error"
        );
    }
}

/// The residue interaction network is built from the whole interaction engine,
/// aggregated to the residue level, and exports as node-link JSON.
#[test]
fn the_network_exports_as_node_link_json() {
    let run = run!("1CRN", &["rin", "--format", "json", "--type", "disulfide"]);
    assert!(run.ok, "command failed: {}", run.stderr);

    let parsed: serde_json::Value = serde_json::from_str(&run.stdout).expect("valid json");
    let edges = parsed["edges"].as_array().expect("edges");
    // Crambin's three disulfides join three residue pairs.
    assert_eq!(edges.len(), 3);
    assert!(edges.iter().all(|e| e["kind"] == "disulfide"));

    let nodes = parsed["nodes"].as_array().expect("nodes");
    // Six cysteines, each an endpoint of one bridge.
    assert_eq!(nodes.len(), 6);
    assert!(nodes.iter().all(|n| n["resn"] == "CYS" && n["degree"] == 1));
}

/// The GraphML export is well-formed XML that an importer can read.
#[test]
fn the_network_exports_as_valid_graphml() {
    let run = run!("1CRN", &["rin", "--format", "graphml"]);
    assert!(run.ok, "command failed: {}", run.stderr);

    assert!(
        run.stdout
            .starts_with(r#"<?xml version="1.0" encoding="UTF-8"?>"#)
    );
    // Every opened element is closed.
    for tag in ["node", "edge"] {
        assert_eq!(
            run.stdout.matches(&format!("<{tag} ")).count(),
            run.stdout.matches(&format!("</{tag}>")).count(),
            "unbalanced {tag} tags",
        );
    }
}

/// `-o` writes the network to a file rather than standard output.
#[test]
fn the_network_can_be_written_to_a_file() {
    let Some(structure) = cached("1CRN") else {
        return;
    };
    let dir = std::env::temp_dir().join(format!("pixelfold-rin-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temp dir");
    let out = dir.join("net.json");

    let status = Command::new(env!("CARGO_BIN_EXE_pixelfold"))
        .args(["rin", "--format", "json", "-o"])
        .arg(&out)
        .arg(&structure)
        .status()
        .expect("the binary runs");
    assert!(status.success());

    let text = std::fs::read_to_string(&out).expect("the file was written");
    let parsed: serde_json::Value = serde_json::from_str(&text).expect("valid json");
    assert!(parsed["nodes"].is_array());

    std::fs::remove_dir_all(&dir).ok();
}

/// Reading a file needs no cache, so an unusable cache directory must not stop
/// it. Only a download touches the cache.
#[test]
fn a_local_file_does_not_need_a_usable_cache_directory() {
    let run = run!(
        "1CRN",
        &["ss", "--format", "tsv", "--cache-dir", "/nonexistentroot/x"]
    );

    assert!(run.ok, "a local file needed the cache: {}", run.stderr);
    assert!(!run.rows().is_empty());
    assert!(
        !Path::new("/nonexistentroot").exists(),
        "the cache directory was created for a purely local read"
    );
}
