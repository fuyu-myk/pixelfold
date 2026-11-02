//! The headless commands over a real structure, through the built binary.
//!
//! The fixture is committed beside these tests rather than taken from the
//! benchmark cache, which is gitignored: a cache-sourced fixture is absent on
//! every clean checkout, and these tests used to skip themselves when it was,
//! so CI passed them without running the binary once.

use std::collections::HashMap;
use std::ffi::OsStr;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// A structure committed under `tests/fixtures`.
fn fixture(id: &str) -> PathBuf {
    fixture_file(&format!("{id}.cif"))
}

/// The legacy PDB twin of a committed fixture.
fn pdb_fixture(id: &str) -> PathBuf {
    fixture_file(&format!("{id}.pdb"))
}

fn fixture_file(name: &str) -> PathBuf {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name);

    path.canonicalize()
        .unwrap_or_else(|e| panic!("missing test fixture {}: {e}", path.display()))
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

/// Run the binary against a committed fixture.
fn pixelfold(id: &str, args: &[&str]) -> Run {
    run_on(&fixture(id), args)
}

/// Run the binary against a structure file.
fn run_on(structure: &Path, args: &[&str]) -> Run {
    run_args(args.iter().map(OsStr::new).chain([structure.as_os_str()]))
}

/// Run the binary with exactly these arguments.
fn run_args<'a>(args: impl IntoIterator<Item = &'a OsStr>) -> Run {
    let output = Command::new(env!("CARGO_BIN_EXE_pixelfold"))
        .args(args)
        .output()
        .expect("the binary runs");

    Run {
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        ok: output.status.success(),
    }
}

macro_rules! run {
    ($id:expr, $args:expr) => {
        pixelfold($id, $args)
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
///
/// The pipe is built here rather than by a shell.
#[test]
fn a_closed_pipe_is_not_a_failure_in_any_format() {
    let structure = fixture("1PPB");

    for format in ["table", "tsv", "json"] {
        let mut child = Command::new(env!("CARGO_BIN_EXE_pixelfold"))
            .args(["interactions", "--format", format])
            .arg(&structure)
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .expect("the binary runs");

        // Emulate what `| head -1` does.
        let stdout = child.stdout.take().expect("stdout was piped");
        let mut reader = BufReader::new(stdout);
        let mut first = String::new();
        reader
            .read_line(&mut first)
            .expect("the first line arrives");
        drop(reader);

        let status = child.wait().expect("the child exits");
        assert!(
            status.success(),
            "{format} reported a closed pipe as an error (exit {status})"
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

/// `--analyze` reports the network's structure instead of exporting the graph.
#[test]
fn the_analysis_report_summarises_the_network_structure() {
    // Crambin's three disulfides are three isolated residue pairs: three
    // components, no residue bridges another, and (leaves never disconnect a
    // graph) no cut residues.
    let run = run!("1CRN", &["rin", "--type", "disulfide", "--analyze"]);
    assert!(run.ok, "command failed: {}", run.stderr);

    assert!(
        run.stdout.contains("6 residues, 3 edges, 3 components"),
        "summary line missing: {}",
        run.stdout
    );
    assert!(run.stdout.contains("Top hubs by betweenness centrality:"));
    assert!(run.stdout.contains("Articulation points (0 cut residues):"));
}

/// `--mvsj` emits a MolViewSpec scene tree that Mol* can load.
#[test]
fn mvsj_export_is_a_valid_molviewspec_scene() {
    let run = run!("1CRN", &["rin", "--mvsj"]);
    assert!(run.ok, "command failed: {}", run.stderr);

    let scene: serde_json::Value = serde_json::from_str(&run.stdout).expect("valid json");
    assert_eq!(scene["metadata"]["version"], "1");

    let download = &scene["root"]["children"][0];
    assert_eq!(download["kind"], "download");
    assert!(
        download["params"]["url"]
            .as_str()
            .unwrap()
            .starts_with("file://"),
        "a local structure becomes a file URL"
    );

    let structure = &download["children"][0]["children"][0];
    assert_eq!(structure["params"]["type"], "model");

    let components = structure["children"].as_array().expect("components");
    assert_eq!(components[0]["params"]["selector"], "polymer");
    assert!(
        components[1]["params"]["selector"].is_array(),
        "the network's residues are a selector array"
    );
}

// Mol* draws a file's first model unless told otherwise, whatever model the
// network was built over.
#[test]
fn a_scene_shows_the_models_its_network_spans() {
    let dir = tempfile::tempdir().unwrap();
    let path = with_a_second_model("1CRN", dir.path(), 40, 60.0);
    let shown = |args: &[&str]| -> Vec<serde_json::Value> {
        let run = run_on(&path, &[&["rin", "--mvsj"], args].concat());
        assert!(run.ok, "{}", run.stderr);
        let scene: serde_json::Value = serde_json::from_str(&run.stdout).unwrap();
        let parse = &scene["root"]["children"][0]["children"][0];
        parse["children"]
            .as_array()
            .unwrap()
            .iter()
            .map(|structure| structure["params"]["model_index"].clone())
            .collect()
    };

    assert_eq!(shown(&[]), [0]);
    assert_eq!(shown(&["--model", "2"]), [1]);
    assert_eq!(shown(&["--model", "all"]), [0, 1]);
}

// An agent's shell gives no terminal: the interactive commands must say so
// and point at the headless ones, not panic.
#[test]
fn the_interactive_commands_refuse_without_a_terminal() {
    let fixture = fixture("1PPB");
    let path = fixture.to_str().unwrap();
    for args in [
        vec![path],
        vec!["view", path],
        vec![
            "poses",
            "--reference",
            path,
            "--ligand",
            "0G6",
            path,
            "--view",
        ],
        vec!["fetch", "kinase"],
    ] {
        let run = run_args(args.iter().map(OsStr::new));

        assert!(!run.ok, "{args:?}");
        assert!(
            run.stderr.contains("needs a terminal"),
            "{args:?}: {}",
            run.stderr
        );
        assert!(!run.stderr.contains("panicked"), "{args:?}: {}", run.stderr);
    }
}

// The structural report and a route, for a script or an agent: the text
// stays the default, JSON when asked for.
#[test]
fn the_network_report_and_a_route_come_as_json_when_asked() {
    let report = run!("1CRN", &["rin", "--analyze", "--format", "json"]);
    let report: serde_json::Value = serde_json::from_str(&report.stdout).expect("json");
    let residues = report["residues"].as_array().unwrap();
    assert_eq!(report["models"], 1);
    assert_eq!(
        residues.len(),
        report["components"]
            .as_array()
            .unwrap()
            .iter()
            .map(|c| c.as_array().unwrap().len())
            .sum::<usize>()
    );
    let between = |at: usize| residues[at]["betweenness"].as_f64().unwrap();
    assert!(between(0) >= between(1), "the hubs first");
    for field in ["id", "resn", "ss", "degree", "betweenness", "articulation"] {
        assert!(
            residues[0].get(field).is_some(),
            "no {field}: {}",
            residues[0]
        );
    }

    let text = run!("1CRN", &["rin", "--analyze"]);
    assert!(
        text.stdout.starts_with("Residue interaction network"),
        "{}",
        text.stdout
    );

    let route = run!(
        "1CRN",
        &["rin", "--path", "A/3", "A/40", "--format", "json"]
    );
    let route: serde_json::Value = serde_json::from_str(&route.stdout).expect("json");
    assert_eq!(route["from"], "A/3");
    assert_eq!(route["to"], "A/40");
    let path = route["path"].as_array().unwrap();
    assert_eq!(path.first().unwrap()["id"], "A/3");
    assert_eq!(path.last().unwrap()["id"], "A/40");
    assert_eq!(route["steps"], path.len() - 1);

    let graphml = run!("1CRN", &["rin", "--analyze", "--format", "graphml"]);
    assert!(!graphml.ok);
    assert!(graphml.stderr.contains("--analyze"), "{}", graphml.stderr);
}

/// `--path` reports the shortest chain of interactions between two residues.
#[test]
fn shortest_path_reports_a_route_between_residues() {
    // Path a residue that is certainly in the network (the first node) to itself:
    // a zero-length route that always exists.
    let json = run!("1CRN", &["rin", "--format", "json"]);
    let parsed: serde_json::Value = serde_json::from_str(&json.stdout).expect("valid json");
    let id = parsed["nodes"][0]["id"]
        .as_str()
        .expect("a node id")
        .to_owned();

    let run = run!("1CRN", &["rin", "--path", id.as_str(), id.as_str()]);
    assert!(run.ok, "command failed: {}", run.stderr);
    assert!(run.stdout.contains("Shortest path"), "{}", run.stdout);
    assert!(
        run.stdout.contains("1 residues"),
        "a self-path is one residue: {}",
        run.stdout
    );
}

/// `--path` errors when a named residue is not in the network.
#[test]
fn shortest_path_errors_for_an_absent_residue() {
    let run = run!("1CRN", &["rin", "--path", "Z/999", "A/1"]);
    assert!(!run.ok, "should fail for a residue not in the network");
    assert!(run.stderr.contains("not in the network"), "{}", run.stderr);
}

/// `sasa` reports relative accessibility beside the absolute area.
#[test]
fn sasa_reports_relative_accessibility() {
    let run = run!("1CRN", &["sasa", "--format", "tsv"]);
    assert!(run.ok, "command failed: {}", run.stderr);

    assert_eq!(
        run.stdout.lines().next(),
        Some("chain\tresi\ticode\tresn\tsasa\trsa"),
    );
    // Crambin is all standard residues, so every row carries a numeric RSA.
    for row in run.stdout.lines().skip(1) {
        let rsa = row.split('\t').nth(5).expect("an rsa column");
        assert!(rsa.parse::<f32>().is_ok(), "rsa not numeric in {row:?}");
    }
}

/// `-o` writes the network to a file rather than standard output.
#[test]
fn the_network_can_be_written_to_a_file() {
    let structure = fixture("1CRN");
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

/// Parse (width, height) from a PNG's IHDR without decoding the image.
fn png_dimensions(bytes: &[u8]) -> Option<(u32, u32)> {
    const SIGNATURE: &[u8] = &[137, 80, 78, 71, 13, 10, 26, 10];
    if bytes.len() < 24 || &bytes[..8] != SIGNATURE || &bytes[12..16] != b"IHDR" {
        return None;
    }
    let width = u32::from_be_bytes(bytes[16..20].try_into().ok()?);
    let height = u32::from_be_bytes(bytes[20..24].try_into().ok()?);

    Some((width, height))
}

#[test]
fn render_writes_a_png_of_the_requested_size() {
    let structure = fixture("1CRN");

    let dir = std::env::temp_dir().join(format!("pixelfold-render-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temp dir");
    let out = dir.join("crambin.png");

    let output = Command::new(env!("CARGO_BIN_EXE_pixelfold"))
        .args(["render", "--width", "200", "--height", "150", "-o"])
        .arg(&out)
        .arg(&structure)
        .output()
        .expect("the binary runs");
    assert!(
        output.status.success(),
        "render failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bytes = std::fs::read(&out).expect("the PNG was written");
    assert_eq!(
        png_dimensions(&bytes),
        Some((200, 150)),
        "output is not a 200x150 PNG"
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn render_rejects_an_empty_selection() {
    let structure = fixture("1CRN");

    let dir = std::env::temp_dir().join(format!("pixelfold-render-empty-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temp dir");
    let out = dir.join("nothing.png");

    let output = Command::new(env!("CARGO_BIN_EXE_pixelfold"))
        .args(["render", "--select", "resn STI", "-o"])
        .arg(&out)
        .arg(&structure)
        .output()
        .expect("the binary runs");

    assert!(
        !output.status.success(),
        "a selection of nothing should fail"
    );
    assert!(
        !out.exists(),
        "no PNG should be written for an empty selection"
    );

    std::fs::remove_dir_all(&dir).ok();
}

/// Run render with the given extra args and return raw stdout bytes.
fn render_stdout(args: &[&str]) -> Option<Vec<u8>> {
    let structure = fixture("1CRN");
    let output = Command::new(env!("CARGO_BIN_EXE_pixelfold"))
        .arg("render")
        .args(args)
        .arg(&structure)
        .output()
        .expect("the binary runs");
    assert!(
        output.status.success(),
        "render failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    Some(output.stdout)
}

/// Without a file and with a piped (non-terminal) stdout, render writes a PNG so
/// it composes in a pipeline.
#[test]
fn render_without_output_pipes_a_png() {
    let Some(out) = render_stdout(&["--width", "80", "--height", "60"]) else {
        return;
    };
    assert_eq!(
        png_dimensions(&out),
        Some((80, 60)),
        "piped output is not a PNG"
    );
}

/// An explicit protocol is honored even when piped, so its escape can be
/// captured. Kitty output is the graphics-protocol APC escape, not a PNG.
#[test]
fn render_protocol_kitty_emits_the_graphics_escape() {
    let Some(out) = render_stdout(&["--protocol", "kitty", "--width", "40", "--height", "30"])
    else {
        return;
    };
    assert!(out.starts_with(b"\x1b_Ga=T,f=100,"), "not a Kitty escape");
    assert!(out.ends_with(b"\x1b\\\n"), "Kitty escape not terminated");
}

/// A depth-band slab renders through the same path and yields a valid PNG.
#[test]
fn render_slab_produces_a_valid_png() {
    let Some(out) = render_stdout(&["--slab", "0.4", "0.6", "--width", "60", "--height", "45"])
    else {
        return;
    };
    assert_eq!(
        png_dimensions(&out),
        Some((60, 45)),
        "sliced output is not a PNG"
    );
}

#[test]
fn render_slab_rejects_an_inverted_range() {
    let structure = fixture("1CRN");
    let dir = std::env::temp_dir().join(format!("pixelfold-slab-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temp dir");
    let out = dir.join("bad.png");

    let output = Command::new(env!("CARGO_BIN_EXE_pixelfold"))
        .args(["render", "--slab", "0.8", "0.2", "-o"])
        .arg(&out)
        .arg(&structure)
        .output()
        .expect("the binary runs");

    assert!(!output.status.success(), "FAR >= NEAR should fail");
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("FAR"),
        "error should explain the slab range"
    );
    std::fs::remove_dir_all(&dir).ok();
}

/// Half-block output is ANSI truecolor text ending each row with a reset.
#[test]
fn render_protocol_halfblock_emits_ansi_text() {
    let Some(out) = render_stdout(&[
        "--protocol",
        "half-block",
        "--width",
        "20",
        "--height",
        "10",
    ]) else {
        return;
    };
    assert!(out.starts_with(b"\x1b[38;2;"), "not truecolor foreground");
    // The upper half block, U+2580, is E2 96 80 in UTF-8.
    assert!(
        out.windows(3).any(|w| w == [0xe2, 0x96, 0x80]),
        "no half-block glyphs emitted"
    );
    assert!(out.windows(4).any(|w| w == b"\x1b[0m"), "no row reset");
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

const TWO_MODELS: &str = "\
data_NMR
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA GLY A 1 0.0 0.0 0.0 1
ATOM 2 C CA GLY A 2 3.8 0.0 0.0 1
ATOM 3 C CA GLY A 1 0.5 0.0 0.0 2
ATOM 4 C CA GLY A 2 4.3 0.0 0.0 2
";

#[test]
fn model_chooses_which_model_of_a_file_is_analysed() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("nmr.cif");
    std::fs::write(&path, TWO_MODELS).unwrap();

    let second = run_on(
        &path,
        &[
            "ss", "--format", "tsv", "--model", "2", "--select", "model 2",
        ],
    );
    assert!(second.ok, "command failed: {}", second.stderr);
    assert_eq!(second.rows().len(), 2);

    let absent = run_on(&path, &["ss", "--model", "5"]);
    assert!(!absent.ok);
    assert!(
        absent
            .stderr
            .contains("there is no model 5; the file has models 1-2"),
        "{}",
        absent.stderr
    );
}

fn two_model_file(dir: &tempfile::TempDir) -> PathBuf {
    let path = dir.path().join("nmr.cif");
    std::fs::write(&path, TWO_MODELS).unwrap();
    path
}

#[test]
fn model_all_reports_every_model_with_its_number() {
    let dir = tempfile::tempdir().unwrap();
    let path = two_model_file(&dir);

    let tsv = run_on(&path, &["ss", "--format", "tsv", "--model", "all"]);
    assert!(tsv.ok, "command failed: {}", tsv.stderr);
    assert!(tsv.stdout.starts_with("model\tchain\t"), "{}", tsv.stdout);
    let models: Vec<&str> = tsv
        .rows()
        .iter()
        .map(|row| row.split('\t').next().unwrap())
        .collect();
    assert_eq!(models, ["1", "1", "2", "2"]);

    let json = run_on(&path, &["sasa", "--format", "json", "--model", "all"]);
    assert!(json.ok, "command failed: {}", json.stderr);
    let records: serde_json::Value = serde_json::from_str(&json.stdout).unwrap();
    assert_eq!(records[2]["model"], 2);
    assert_eq!(records[2]["resi"], 1);
}

// A selection may match in some models and not others; only a selection that
// matches nowhere is an error.
#[test]
fn model_all_errors_only_when_no_model_matches() {
    let dir = tempfile::tempdir().unwrap();
    let path = two_model_file(&dir);

    let second = run_on(
        &path,
        &[
            "ss", "--format", "tsv", "--model", "all", "--select", "model 2",
        ],
    );
    assert!(second.ok, "command failed: {}", second.stderr);
    assert_eq!(second.rows().len(), 2);

    let nowhere = run_on(&path, &["ss", "--model", "all", "--select", "resn TRP"]);
    assert!(!nowhere.ok);
    assert!(
        nowhere.stderr.contains("matched no atoms in any model"),
        "{}",
        nowhere.stderr
    );
}

#[test]
fn single_model_commands_refuse_every_model() {
    let dir = tempfile::tempdir().unwrap();
    let path = two_model_file(&dir);
    let png = dir.path().join("out.png");

    let run = run_on(
        &path,
        &["render", "--model", "all", "-o", png.to_str().unwrap()],
    );
    assert!(!run.ok, "render should refuse --model all");
    assert!(run.stderr.contains("one model at a time"), "{}", run.stderr);
}

/// A copy of fixture `id` in `dir` with a second model after the first, the
/// same but for residue `moved` of chain A, `shift` A along x.
fn with_a_second_model(id: &str, dir: &Path, moved: i32, shift: f32) -> PathBuf {
    let text = std::fs::read_to_string(fixture(id)).unwrap();
    let column = atom_site_columns(&text);
    let is_atom = |line: &str| line.starts_with("ATOM") || line.starts_with("HETATM");

    let second: Vec<String> = text
        .lines()
        .filter(|line| is_atom(line))
        .map(|line| {
            let mut fields: Vec<String> = line.split_whitespace().map(String::from).collect();
            let serial: usize = fields[column["id"]].parse().unwrap();
            fields[column["id"]] = (serial + 10_000).to_string();
            fields[column["pdbx_PDB_model_num"]] = "2".to_string();
            if fields[column["auth_asym_id"]] == "A"
                && fields[column["auth_seq_id"]] == moved.to_string()
            {
                let x: f32 = fields[column["Cartn_x"]].parse().unwrap();
                fields[column["Cartn_x"]] = format!("{:.3}", x + shift);
            }
            fields.join(" ")
        })
        .collect();
    let after = text
        .lines()
        .enumerate()
        .filter(|&(_, line)| is_atom(line))
        .last()
        .unwrap()
        .0
        + 1;
    let lines: Vec<String> = text
        .lines()
        .take(after)
        .map(String::from)
        .chain(second)
        .chain(text.lines().skip(after).map(String::from))
        .collect();

    let path = dir.join(format!("{id}-two-models.cif"));
    std::fs::write(&path, lines.join("\n")).unwrap();
    path
}

// Crambin with Cys40 carried off in its second model: the 3-40 disulfide is
// present in one model of two, the 4-32 disulfide in both.
#[test]
fn the_network_of_every_model_weighs_each_edge_by_the_models_it_is_in() {
    let dir = tempfile::tempdir().unwrap();
    let path = with_a_second_model("1CRN", dir.path(), 40, 60.0);

    let run = run_on(
        &path,
        &[
            "rin",
            "--model",
            "all",
            "--type",
            "disulfide",
            "--format",
            "json",
        ],
    );
    assert!(run.ok, "command failed: {}", run.stderr);
    let network: serde_json::Value = serde_json::from_str(&run.stdout).unwrap();
    assert_eq!(network["models"], 2);
    let frequency = |source: &str, target: &str| {
        network["edges"]
            .as_array()
            .unwrap()
            .iter()
            .find(|edge| edge["source"] == source && edge["target"] == target)
            .map(|edge| (edge["frequency"].as_f64().unwrap(), edge["models"].clone()))
    };
    assert_eq!(frequency("A/3", "A/40"), Some((0.5, 1.into())));
    assert_eq!(frequency("A/4", "A/32"), Some((1.0, 2.into())));

    let frequent = run_on(
        &path,
        &[
            "rin",
            "--model",
            "all",
            "--min-frequency",
            "1",
            "--format",
            "tsv",
        ],
    );
    assert!(frequent.ok, "command failed: {}", frequent.stderr);
    assert!(frequent.stdout.starts_with("source\ttarget\t"));
    assert!(
        frequent
            .stdout
            .lines()
            .next()
            .unwrap()
            .ends_with("\tmodels\tfrequency")
    );
    assert!(!frequent.rows().is_empty());
    assert!(
        frequent
            .rows()
            .iter()
            .all(|row| !row.contains("A/40\t") && row.ends_with("\t2\t1.000")),
        "{}",
        frequent.stdout
    );
}

#[test]
fn a_frequency_threshold_needs_every_model_and_a_fraction() {
    let dir = tempfile::tempdir().unwrap();
    let path = two_model_file(&dir);

    let one = run_on(&path, &["rin", "--min-frequency", "0.5"]);
    assert!(!one.ok);
    assert!(one.stderr.contains("--model all"), "{}", one.stderr);

    let beyond = run_on(&path, &["rin", "--model", "all", "--min-frequency", "1.5"]);
    assert!(!beyond.ok);
    assert!(
        beyond.stderr.contains("between 0 and 1"),
        "{}",
        beyond.stderr
    );
}

/// Fixture `id` and a copy with residue `moved` of chain A carried `shift` A
/// along x, as separate files.
fn two_samples(id: &str, dir: &Path, moved: i32, shift: f32) -> [PathBuf; 2] {
    let copy = edited_fixture(id, dir, |fields, column| {
        if fields[column["auth_asym_id"]] == "A"
            && fields[column["auth_seq_id"]] == moved.to_string()
        {
            let x: f32 = fields[column["Cartn_x"]].parse().unwrap();
            fields[column["Cartn_x"]] = format!("{:.3}", x + shift);
        }
    });

    [fixture(id), copy]
}

fn network_of(run: &Run) -> serde_json::Value {
    assert!(run.ok, "command failed: {}", run.stderr);
    serde_json::from_str(&run.stdout).unwrap()
}

/// The frequency and model count of the edge from `source` to `target`.
fn edge_weight(network: &serde_json::Value, source: &str, target: &str) -> Option<(f64, u64)> {
    network["edges"]
        .as_array()
        .unwrap()
        .iter()
        .find(|edge| edge["source"] == source && edge["target"] == target)
        .map(|edge| {
            (
                edge["frequency"].as_f64().unwrap(),
                edge["models"].as_u64().unwrap(),
            )
        })
}

// Crambin and a copy with Cys40 carried off, as a co-folding model writes two
// samples: the 3-40 disulfide is in one file of two, the 4-32 disulfide in both.
#[test]
fn the_network_of_several_files_weighs_each_edge_by_the_files_it_is_in() {
    let dir = tempfile::tempdir().unwrap();
    let [crambin, moved] = two_samples("1CRN", dir.path(), 40, 60.0);
    let rin = |extra: &[&str], files: &[&Path]| {
        let args = ["rin", "--type", "disulfide", "--format", "json"];
        network_of(&run_args(
            args.iter()
                .chain(extra)
                .map(OsStr::new)
                .chain(files.iter().map(|file| file.as_os_str())),
        ))
    };

    let both = rin(&[], &[&crambin, &moved]);
    assert_eq!(both["models"], 2);
    assert_eq!(edge_weight(&both, "A/3", "A/40"), Some((0.5, 1)));
    assert_eq!(edge_weight(&both, "A/4", "A/32"), Some((1.0, 2)));

    let frequent = rin(&["--min-frequency", "1"], &[&crambin, &moved]);
    assert_eq!(edge_weight(&frequent, "A/3", "A/40"), None);
    assert_eq!(edge_weight(&frequent, "A/4", "A/32"), Some((1.0, 2)));

    let ensemble = with_a_second_model("1CRN", dir.path(), 40, 60.0);
    let every = rin(&["--model", "all"], &[&ensemble, &crambin]);
    assert_eq!(
        every["models"], 3,
        "both models of the first file, then the second"
    );
    assert_eq!(
        edge_weight(&every, "A/3", "A/40").map(|(_, models)| models),
        Some(2)
    );
}

#[test]
fn files_that_name_a_residue_differently_are_refused() {
    let dir = tempfile::tempdir().unwrap();
    let mutant = edited_fixture("1CRN", dir.path(), |fields, column| {
        if fields[column["auth_seq_id"]] == "1" {
            fields[column["label_comp_id"]] = "ALA".to_string();
            fields[column["auth_comp_id"]] = "ALA".to_string();
        }
    });

    let run = run_args([
        OsStr::new("rin"),
        fixture("1CRN").as_os_str(),
        mutant.as_os_str(),
    ]);

    assert!(
        !run.ok,
        "a renumbered or different sequence must not fold in"
    );
    for named in ["1CRN-edited.cif", "A/1", "THR", "ALA"] {
        assert!(run.stderr.contains(named), "no {named}: {}", run.stderr);
    }
}

/// Crambin with Ile25 also modelled as a leucine, conformer B at 0.4
/// occupancy, as 1EJG has it.
fn microheterogeneous(dir: &Path) -> PathBuf {
    let text = std::fs::read_to_string(fixture("1CRN")).unwrap();
    let column = atom_site_columns(&text);
    let mut lines = Vec::new();
    let mut leucine = Vec::new();
    for line in text.lines() {
        let fields: Vec<String> = line.split_whitespace().map(String::from).collect();
        let is_atom = line.starts_with("ATOM");
        if !(is_atom && fields[column["auth_seq_id"]] == "25") {
            if is_atom {
                lines.append(&mut leucine);
            }
            lines.push(line.to_string());
            continue;
        }
        let mut isoleucine = fields.clone();
        isoleucine[column["label_alt_id"]] = "A".to_string();
        isoleucine[column["occupancy"]] = "0.60".to_string();
        lines.push(isoleucine.join(" "));
        let mut other = fields;
        let serial: usize = other[column["id"]].parse().unwrap();
        other[column["id"]] = (serial + 10_000).to_string();
        other[column["label_alt_id"]] = "B".to_string();
        other[column["occupancy"]] = "0.40".to_string();
        for name in ["label_comp_id", "auth_comp_id"] {
            other[column[name]] = "LEU".to_string();
        }
        for name in ["label_atom_id", "auth_atom_id"] {
            let renamed = match other[column[name]].as_str() {
                "CG1" => "CG",
                "CG2" => "CD2",
                kept => kept,
            };
            other[column[name]] = renamed.to_string();
        }
        leucine.push(other.join(" "));
    }

    let path = dir.join("micro.cif");
    std::fs::write(&path, lines.join("\n")).unwrap();
    path
}

#[test]
fn a_residue_holding_two_names_folds_in_beside_itself() {
    let dir = tempfile::tempdir().unwrap();
    let micro = microheterogeneous(dir.path());

    let run = run_args([
        OsStr::new("rin"),
        OsStr::new("--format"),
        OsStr::new("json"),
        micro.as_os_str(),
        micro.as_os_str(),
    ]);

    assert_eq!(network_of(&run)["models"], 2);
}

#[test]
fn each_file_must_hold_the_model_asked_for() {
    let dir = tempfile::tempdir().unwrap();
    let ensemble = with_a_second_model("1CRN", dir.path(), 40, 60.0);

    let run = run_args([
        OsStr::new("rin"),
        OsStr::new("--model"),
        OsStr::new("2"),
        ensemble.as_os_str(),
        fixture("1CRN").as_os_str(),
    ]);

    assert!(!run.ok);
    assert!(run.stderr.contains("1CRN.cif"), "{}", run.stderr);
    assert!(run.stderr.contains("there is no model 2"), "{}", run.stderr);
}

#[test]
fn a_model_that_fails_to_load_under_model_all_is_named() {
    let dir = tempfile::tempdir().unwrap();
    let path = two_model_file(&dir);

    let run = run_on(&path, &["ss", "--model", "all", "--assembly=9"]);

    assert!(!run.ok);
    assert!(run.stderr.contains("model 1"), "{}", run.stderr);
}

#[test]
fn model_takes_a_number_or_all() {
    let dir = tempfile::tempdir().unwrap();
    let path = two_model_file(&dir);

    let run = run_on(&path, &["ss", "--model", "first"]);
    assert!(!run.ok);
    assert!(
        run.stderr.contains("a model number or `all`"),
        "{}",
        run.stderr
    );
}

// The PDB and mmCIF readers must agree wherever a legacy PDB file carries the
// same information. Thrombin's light chain runs 1H, 1G, ... 1A, 1, the order a
// reader once re-sorted. A PDB file carries no chemical component definitions,
// so thrombin's inhibitor stays untyped there and its interactions differ by
// design.
#[test]
fn pdb_and_mmcif_files_of_one_entry_report_alike() {
    let cases: [(&str, &[&str]); 2] = [
        ("1UBQ", &["ss", "sasa", "interactions"]),
        ("1PPB", &["ss", "sasa"]),
    ];
    for (id, reports) in cases {
        for report in reports {
            let sorted = |path: &Path| {
                let run = run_on(path, &[report, "--format", "tsv"]);
                assert!(run.ok, "{id} {report} failed: {}", run.stderr);
                let mut rows: Vec<String> = run.stdout.lines().map(str::to_owned).collect();
                rows.sort();
                rows
            };
            assert_eq!(
                sorted(&pdb_fixture(id)),
                sorted(&fixture(id)),
                "{id} {report}"
            );
        }
    }
}

const ASSEMBLY: &str = "\
data_ASM
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
ATOM 1 C CA GLY A 1 1.0 0.0 0.0 1 A
ATOM 2 C CA GLY A 2 4.8 0.0 0.0 2 A
ATOM 3 C CA GLY C 1 1.0 20.0 0.0 1 C
#
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 0
2 'crystal symmetry operation' -1 0 0 10 0 -1 0 0 0 0 1 0
#
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1,2 A
2 1 C
#
";

fn assembly_file(dir: &tempfile::TempDir) -> PathBuf {
    let path = dir.path().join("assembly.cif");
    std::fs::write(&path, ASSEMBLY).unwrap();
    path
}

fn chains(run: &Run) -> Vec<String> {
    run.rows()
        .iter()
        .map(|row| row.split('\t').next().unwrap().to_string())
        .collect()
}

#[test]
fn assembly_reports_each_copy_under_its_own_chain() {
    let dir = tempfile::tempdir().unwrap();
    let path = assembly_file(&dir);

    let deposited = run_on(&path, &["ss", "--format", "tsv"]);
    assert!(deposited.ok, "command failed: {}", deposited.stderr);
    assert_eq!(chains(&deposited), ["A", "A", "C"]);

    let primary = run_on(&path, &["ss", "--format", "tsv", "--assembly"]);
    assert!(primary.ok, "command failed: {}", primary.stderr);
    assert_eq!(chains(&primary), ["A", "A", "A-2", "A-2"]);

    let second = run_on(&path, &["ss", "--format", "tsv", "--assembly=2"]);
    assert!(second.ok, "command failed: {}", second.stderr);
    assert_eq!(chains(&second), ["C"]);

    let copy = run_on(
        &path,
        &[
            "ss",
            "--format",
            "tsv",
            "--assembly",
            "--select",
            "chain A-2",
        ],
    );
    assert!(copy.ok, "command failed: {}", copy.stderr);
    assert_eq!(copy.rows().len(), 2);
}

#[test]
fn an_unknown_assembly_and_an_assembly_scene_are_refused() {
    let dir = tempfile::tempdir().unwrap();
    let path = assembly_file(&dir);

    let unknown = run_on(&path, &["ss", "--assembly=9"]);
    assert!(!unknown.ok);
    assert!(
        unknown
            .stderr
            .contains("there is no assembly 9; the file has assemblies 1, 2"),
        "{}",
        unknown.stderr
    );

    let scene = run_on(&path, &["rin", "--mvsj", "--assembly"]);
    assert!(!scene.ok);
    assert!(scene.stderr.contains("--mvsj"), "{}", scene.stderr);
}

/// A copy of fixture `id` in `dir`, with `edit` applied to the fields of every
/// atom row; `edit` gets the column indices by name.
fn edited_fixture(
    id: &str,
    dir: &Path,
    mut edit: impl FnMut(&mut [String], &HashMap<&str, usize>),
) -> PathBuf {
    let text = std::fs::read_to_string(fixture(id)).unwrap();
    let columns = atom_site_columns(&text);

    let edited: Vec<String> = text
        .lines()
        .map(|line| {
            if !(line.starts_with("ATOM") || line.starts_with("HETATM")) {
                return line.to_string();
            }
            let mut fields: Vec<String> = line.split_whitespace().map(String::from).collect();
            edit(&mut fields, &columns);
            fields.join(" ")
        })
        .collect();

    let path = dir.join(format!("{id}-edited.cif"));
    std::fs::write(&path, edited.join("\n")).unwrap();
    path
}

/// The index of each `_atom_site` column of an mmCIF file, by name.
fn atom_site_columns(text: &str) -> HashMap<&str, usize> {
    text.lines()
        .filter_map(|line| line.strip_prefix("_atom_site."))
        .enumerate()
        .map(|(index, name)| (name.trim(), index))
        .collect()
}

/// A copy of fixture `id` in `dir` holding a second copy of the `ligand`
/// component, as chain Y moved `shift` A along x and written before the first.
fn with_another_copy(id: &str, dir: &Path, ligand: &str, shift: f32) -> PathBuf {
    let text = std::fs::read_to_string(fixture(id)).unwrap();
    let column = atom_site_columns(&text);
    let is_ligand = |line: &str| {
        line.starts_with("HETATM")
            && line.split_whitespace().nth(column["auth_comp_id"]) == Some(ligand)
    };

    let copy: Vec<String> = text
        .lines()
        .filter(|line| is_ligand(line))
        .map(|line| {
            let mut fields: Vec<String> = line.split_whitespace().map(String::from).collect();
            let serial: usize = fields[column["id"]].parse().unwrap();
            let x: f32 = fields[column["Cartn_x"]].parse().unwrap();
            fields[column["id"]] = (serial + 10_000).to_string();
            fields[column["Cartn_x"]] = format!("{:.3}", x + shift);
            fields[column["label_asym_id"]] = "Y".to_string();
            fields[column["auth_asym_id"]] = "Y".to_string();
            fields[column["auth_seq_id"]] = "1".to_string();
            fields.join(" ")
        })
        .collect();
    let first = text.lines().position(is_ligand).unwrap();
    let lines: Vec<String> = text
        .lines()
        .take(first)
        .map(String::from)
        .chain(copy)
        .chain(text.lines().skip(first).map(String::from))
        .collect();

    let path = dir.join(format!("{id}-two-copies.cif"));
    std::fs::write(&path, lines.join("\n")).unwrap();
    path
}

/// A copy of fixture `id` in `dir` holding a second copy of the whole complex,
/// written after the first and moved 60 A along x: its chains renamed by
/// `chains` and its `ligand` component as chain Y. The first copy's ligand is
/// moved `nudge` A along x.
fn with_another_complex(
    id: &str,
    dir: &Path,
    chains: &[(&str, &str)],
    ligand: &str,
    nudge: f32,
) -> PathBuf {
    let text = std::fs::read_to_string(fixture(id)).unwrap();
    let column = atom_site_columns(&text);
    let renamed: HashMap<&str, &str> = chains.iter().copied().collect();
    let is_atom = |line: &str| line.starts_with("ATOM") || line.starts_with("HETATM");
    let moved = |fields: &mut [String], shift: f32| {
        let x: f32 = fields[column["Cartn_x"]].parse().unwrap();
        fields[column["Cartn_x"]] = format!("{:.3}", x + shift);
    };
    let last = text
        .lines()
        .enumerate()
        .filter(|&(_, line)| is_atom(line))
        .map(|(index, _)| index)
        .last()
        .unwrap();

    let mut lines: Vec<String> = Vec::new();
    let mut second: Vec<String> = Vec::new();
    for (index, line) in text.lines().enumerate() {
        if !is_atom(line) {
            lines.push(line.to_string());
            continue;
        }
        let mut fields: Vec<String> = line.split_whitespace().map(String::from).collect();
        let is_ligand = fields[column["auth_comp_id"]] == ligand;
        let mut copy = fields.clone();
        if is_ligand {
            moved(&mut fields, nudge);
        }
        lines.push(fields.join(" "));

        let serial: usize = copy[column["id"]].parse().unwrap();
        let chain = copy[column["auth_asym_id"]].clone();
        let name = if is_ligand {
            "Y"
        } else {
            renamed[chain.as_str()]
        };
        copy[column["id"]] = (serial + 10_000).to_string();
        copy[column["auth_asym_id"]] = name.to_string();
        copy[column["label_asym_id"]] = format!("{}2", copy[column["label_asym_id"]]);
        moved(&mut copy, 60.0);
        second.push(copy.join(" "));
        if index == last {
            lines.append(&mut second);
        }
    }

    let path = dir.join(format!("{id}-two-complexes.cif"));
    std::fs::write(&path, lines.join("\n")).unwrap();
    path
}

/// Run `poses` with 1PPB as the reference and its PPACK inhibitor (0G6) as the
/// ligand, written as JSON.
fn poses_against(reference: &Path, poses: &[&Path]) -> Run {
    poses_with(&["--ligand", "0G6"], reference, poses)
}

/// Run `poses` with the ligand picked by `ligand`, written as JSON.
fn poses_with(ligand: &[&str], reference: &Path, poses: &[&Path]) -> Run {
    let head = ["poses", "--format", "json"]
        .into_iter()
        .chain(ligand.iter().copied())
        .chain(["--reference"])
        .map(OsStr::new);
    run_args(
        head.chain([reference.as_os_str()])
            .chain(poses.iter().map(|pose| pose.as_os_str())),
    )
}

/// A copy of fixture `id` as co-folding output writes it: each polymer chain
/// renamed by `chains` and numbered from 1 without insertion codes, and the
/// `ligand` component a chain of its own, moved `shift` A along x.
fn cofolded(id: &str, dir: &Path, chains: &[(&str, &str)], ligand: &str, shift: f32) -> PathBuf {
    let renamed: HashMap<&str, &str> = chains.iter().copied().collect();
    let mut numbers: HashMap<(String, String, String), usize> = HashMap::new();
    let mut counts: HashMap<String, usize> = HashMap::new();
    let own = dir.join(format!("cofolded-{shift}"));
    std::fs::create_dir_all(&own).unwrap();
    edited_fixture(id, &own, |fields, column| {
        let chain = fields[column["auth_asym_id"]].clone();
        if fields[column["auth_comp_id"]] == ligand {
            fields[column["auth_asym_id"]] = "Z".to_string();
            fields[column["auth_seq_id"]] = "1".to_string();
            let x: f32 = fields[column["Cartn_x"]].parse().unwrap();
            fields[column["Cartn_x"]] = format!("{:.3}", x + shift);
        } else if let Some(name) = renamed.get(chain.as_str()) {
            let residue = (
                chain.clone(),
                fields[column["auth_seq_id"]].clone(),
                fields[column["pdbx_PDB_ins_code"]].clone(),
            );
            let number = *numbers.entry(residue).or_insert_with(|| {
                let count = counts.entry(chain).or_insert(0);
                *count += 1;
                *count
            });
            fields[column["auth_asym_id"]] = name.to_string();
            fields[column["auth_seq_id"]] = number.to_string();
            fields[column["pdbx_PDB_ins_code"]] = "?".to_string();
        }
    })
}

fn records(run: &Run) -> Vec<serde_json::Value> {
    assert!(run.ok, "command failed: {}", run.stderr);
    let parsed: serde_json::Value = serde_json::from_str(&run.stdout).expect("valid json");

    parsed.as_array().expect("an array").clone()
}

#[test]
fn the_pose_viewer_takes_no_output_format_and_its_scores_stay_in_it() {
    let reference = fixture("1PPB");
    let poses = |flags: &[&str]| {
        let head = ["poses", "--ligand", "0G6"]
            .iter()
            .chain(flags)
            .chain(&["--reference"])
            .map(OsStr::new);
        run_args(head.chain([reference.as_os_str(), reference.as_os_str()]))
    };

    let formatted = poses(&["--view", "--format", "tsv"]);
    assert!(!formatted.ok);
    assert!(
        formatted.stderr.contains("cannot be used with"),
        "{}",
        formatted.stderr
    );
    let protocol = poses(&["--protocol", "half-block"]);
    assert!(!protocol.ok);
    assert!(protocol.stderr.contains("--view"), "{}", protocol.stderr);

    // What drawing a pose needs is not written out.
    let found = records(&poses_against(&reference, &[&reference]));
    assert!(found[0].get("ligand").is_none(), "{}", found[0]);
    assert!(found[0].get("superposition").is_none(), "{}", found[0]);
}

#[test]
fn poses_rank_by_the_reference_contacts_they_keep() {
    let dir = tempfile::tempdir().unwrap();
    // The inhibitor slid 2 A along x, into the protein.
    let displaced = edited_fixture("1PPB", dir.path(), |fields, column| {
        if fields[column["auth_comp_id"]] == "0G6" {
            let x: f32 = fields[column["Cartn_x"]].parse().unwrap();
            fields[column["Cartn_x"]] = format!("{:.3}", x + 2.0);
        }
    });

    let found = records(&poses_against(
        &fixture("1PPB"),
        &[&displaced, &fixture("1PPB")],
    ));
    assert_eq!(found.len(), 2);

    // The identical pose ranks first, whatever the order given: every contact
    // kept, and its covalent bonds to Ser195 and His57 are not clashes.
    let same = &found[0];
    assert!(
        same["pose"].as_str().unwrap().ends_with("1PPB.cif"),
        "{same}"
    );
    assert_eq!(same["model"], 1);
    assert_eq!(same["similarity"], 1.0);
    assert_eq!(same["recovery"], 1.0);
    assert_eq!(same["clashes"], 0);
    assert_eq!(same["rmsd"], 0.0);
    assert_eq!(same["lost"], serde_json::json!([]));
    assert!(same["contacts"].as_u64().unwrap() > 0);
    assert_eq!(same["contacts"], same["shared"]);

    let moved = &found[1];
    assert!(moved["recovery"].as_f64().unwrap() < 1.0, "{moved}");
    assert!(moved["clashes"].as_u64().unwrap() > 0, "{moved}");
    assert!(!moved["lost"].as_array().unwrap().is_empty(), "{moved}");
    // Every ligand atom moved 2 A and the pocket not at all.
    assert_eq!(moved["rmsd"], 2.0, "{moved}");
}

// The cyclic urea's hydroxyl oxygen O4 carried 2 A off along x: its bond, and
// the angles across it, stray from the reference's.
#[test]
fn a_pose_whose_ligand_strays_from_the_reference_s_geometry_says_so() {
    let dir = tempfile::tempdir().unwrap();
    let strayed = edited_fixture("1HVR", dir.path(), |fields, column| {
        if fields[column["auth_comp_id"]] == "XK2" && fields[column["auth_atom_id"]] == "O4" {
            let x: f32 = fields[column["Cartn_x"]].parse().unwrap();
            fields[column["Cartn_x"]] = format!("{:.3}", x + 2.0);
        }
    });

    let found = records(&poses_with(
        &["--ligand", "XK2"],
        &fixture("1HVR"),
        &[&fixture("1HVR"), &strayed],
    ));
    let geometry = |name: &str| {
        found
            .iter()
            .find(|record| record["pose"].as_str().unwrap().ends_with(name))
            .map(|record| record["geometry"].clone())
            .unwrap()
    };

    assert_eq!(
        geometry("1HVR.cif"),
        serde_json::json!({"bonds": 0, "angles": 0, "internal_clashes": 0, "unflat": 0})
    );
    let strayed = geometry("1HVR-edited.cif");
    assert!(strayed["bonds"].as_u64().unwrap() >= 1, "{strayed}");
    assert!(strayed["angles"].as_u64().unwrap() >= 1, "{strayed}");

    let table = run_args(
        ["poses", "--format", "tsv", "--ligand", "XK2", "--reference"]
            .map(OsStr::new)
            .into_iter()
            .chain([
                fixture("1HVR").as_os_str(),
                dir.path().join("1HVR-edited.cif").as_os_str(),
            ]),
    );
    assert!(table.ok, "{}", table.stderr);
    let header: Vec<&str> = table.stdout.lines().next().unwrap().split('\t').collect();
    let column = header.iter().position(|&name| name == "geometry").unwrap();
    let cell = table.rows()[0].split('\t').nth(column).unwrap().to_string();
    assert!(cell.starts_with("bonds "), "{cell}");
}

#[test]
fn a_pose_that_inverts_a_stereocentre_of_the_ligand_names_it() {
    // The phenylalanine's amine and side chain trade places on its alpha
    // carbon, which turns it from D to L.
    let dir = tempfile::tempdir().unwrap();
    let text = std::fs::read_to_string(fixture("1PPB")).unwrap();
    let column = atom_site_columns(&text);
    let position = |name: &str| -> Vec<String> {
        let line = text
            .lines()
            .find(|line| {
                let fields: Vec<&str> = line.split_whitespace().collect();
                line.starts_with("HETATM")
                    && fields[column["auth_comp_id"]] == "0G6"
                    && fields[column["auth_atom_id"]] == name
            })
            .unwrap();
        let fields: Vec<&str> = line.split_whitespace().collect();
        ["Cartn_x", "Cartn_y", "Cartn_z"]
            .map(|axis| fields[column[axis]].to_string())
            .to_vec()
    };
    let (amine, side_chain) = (position("N"), position("CB"));
    let inverted = edited_fixture("1PPB", dir.path(), |fields, column| {
        if fields[column["auth_comp_id"]] != "0G6" {
            return;
        }
        let to = match fields[column["auth_atom_id"]].as_str() {
            "N" => &side_chain,
            "CB" => &amine,
            _ => return,
        };
        for (axis, value) in ["Cartn_x", "Cartn_y", "Cartn_z"].into_iter().zip(to) {
            fields[column[axis]] = value.clone();
        }
    });

    let found = records(&poses_against(
        &fixture("1PPB"),
        &[&inverted, &fixture("1PPB")],
    ));
    let by_file = |suffix: &str| {
        found
            .iter()
            .find(|record| record["pose"].as_str().unwrap().ends_with(suffix))
            .unwrap()
    };
    assert_eq!(by_file("1PPB.cif")["inverted"], serde_json::json!([]));
    assert_eq!(
        by_file("1PPB-edited.cif")["inverted"],
        serde_json::json!(["0G6 H/1 CA"])
    );
}

/// A copy of fixture `id` in `dir`, the whole complex turned a quarter about z
/// and lifted 30 A, as a co-folding model sits in a frame of its own.
fn turned_fixture(id: &str, dir: &Path) -> PathBuf {
    edited_fixture(id, dir, |fields, column| {
        let [x, y, z] = ["Cartn_x", "Cartn_y", "Cartn_z"]
            .map(|axis| fields[column[axis]].parse::<f32>().unwrap());
        fields[column["Cartn_x"]] = format!("{:.3}", -y);
        fields[column["Cartn_y"]] = format!("{:.3}", x);
        fields[column["Cartn_z"]] = format!("{:.3}", z + 30.0);
    })
}

#[test]
fn a_pose_in_another_frame_is_superposed_on_the_pocket() {
    let dir = tempfile::tempdir().unwrap();
    let turned = turned_fixture("1PPB", dir.path());

    let found = records(&poses_against(&fixture("1PPB"), &[&turned]));
    assert_eq!(found[0]["rmsd"], 0.0, "{}", found[0]);
    assert_eq!(found[0]["similarity"], 1.0, "{}", found[0]);
}

#[test]
fn a_pose_without_ligand_chemistry_borrows_the_reference_definition() {
    // The legacy PDB file carries no component definitions; the mmCIF does.
    let borrowed = poses_against(&fixture("1PPB"), &[&pdb_fixture("1PPB")]);
    let found = records(&borrowed);
    assert_eq!(found[0]["similarity"], 1.0, "{}", found[0]);
    assert!(!borrowed.stderr.contains("warning"), "{}", borrowed.stderr);

    // Neither file defines the inhibitor when the reference is the PDB file, so
    // it makes no contacts at all and there is nothing to rank against.
    let undefined = poses_against(&pdb_fixture("1PPB"), &[&pdb_fixture("1PPB")]);
    assert!(
        !undefined.ok,
        "expected a failure, got: {}",
        undefined.stdout
    );
    assert!(
        undefined.stderr.contains("warning: ligand 0G6 H/1")
            && undefined.stderr.contains("makes no contacts"),
        "{}",
        undefined.stderr
    );
}

/// 1PPB's chains as co-folding names them: the light chain first, as A.
const THROMBIN_CHAINS: [(&str, &str); 2] = [("L", "A"), ("H", "B")];

#[test]
fn a_pose_numbered_as_co_folding_output_is_aligned_with_the_reference() {
    let dir = tempfile::tempdir().unwrap();
    let same = cofolded("1PPB", dir.path(), &THROMBIN_CHAINS, "0G6", 0.0);

    let found = records(&poses_against(&fixture("1PPB"), &[&same]));
    assert_eq!(found[0]["similarity"], 1.0, "{}", found[0]);
    assert_eq!(found[0]["clashes"], 0, "{}", found[0]);

    // A displaced pose scores the same in either numbering, and its lost
    // contacts are named in the reference's, insertion codes and all.
    let displaced = edited_fixture("1PPB", dir.path(), |fields, column| {
        if fields[column["auth_comp_id"]] == "0G6" {
            let x: f32 = fields[column["Cartn_x"]].parse().unwrap();
            fields[column["Cartn_x"]] = format!("{:.3}", x + 2.0);
        }
    });
    let renumbered = cofolded("1PPB", dir.path(), &THROMBIN_CHAINS, "0G6", 2.0);
    let found = records(&poses_against(&fixture("1PPB"), &[&displaced, &renumbered]));
    for field in [
        "contacts",
        "shared",
        "recovery",
        "similarity",
        "rmsd",
        "clashes",
        "lost",
    ] {
        assert_eq!(found[0][field], found[1][field], "{field}");
    }
    let lost = found[1]["lost"].to_string();
    assert!(lost.contains(" H/60A "), "{lost}");
}

#[test]
fn a_homodimer_pose_with_its_chains_named_the_other_way_round_keeps_every_contact() {
    // HIV-1 protease binds its inhibitor between two identical chains, so only
    // the contacts tell which pose chain stands for which reference chain.
    let dir = tempfile::tempdir().unwrap();
    let swapped = cofolded("1HVR", dir.path(), &[("A", "B"), ("B", "A")], "XK2", 0.0);

    let found = records(&poses_with(
        &["--ligand", "XK2"],
        &fixture("1HVR"),
        &[&swapped],
    ));
    assert_eq!(found[0]["similarity"], 1.0, "{}", found[0]);
    assert_eq!(found[0]["clashes"], 0, "{}", found[0]);
    assert_eq!(found[0]["rmsd"], 0.0, "{}", found[0]);
}

#[test]
fn a_pose_of_another_protein_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let unrelated = edited_fixture("1PPB", dir.path(), |fields, column| {
        if fields[column["group_PDB"]] == "ATOM" {
            fields[column["auth_comp_id"]] = "ALA".to_string();
            fields[column["label_comp_id"]] = "ALA".to_string();
        }
    });

    let run = poses_against(&fixture("1PPB"), &[&unrelated]);
    assert!(!run.ok, "expected a failure, got: {}", run.stdout);
    assert!(
        run.stderr.contains("aligns by sequence"),
        "unhelpful error: {}",
        run.stderr
    );
}

#[test]
fn the_reference_ligand_can_be_picked_apart_from_the_poses() {
    // A selection naming the reference's chain finds no ligand in a pose that
    // puts it in a chain of its own; --reference-select picks it apart.
    let dir = tempfile::tempdir().unwrap();
    let pose = cofolded("1PPB", dir.path(), &THROMBIN_CHAINS, "0G6", 0.0);
    let reference_only = [
        "--ligand",
        "0G6",
        "--reference-select",
        "chain H and resn 0G6",
    ];

    let found = records(&poses_with(&reference_only, &fixture("1PPB"), &[&pose]));
    assert_eq!(found[0]["similarity"], 1.0, "{}", found[0]);
}

#[test]
fn a_pose_holding_several_copies_of_the_ligand_is_scored_by_the_one_in_place() {
    // A second inhibitor 40 A away, clear of the protein, written first.
    let dir = tempfile::tempdir().unwrap();
    let pose = with_another_copy("1HVR", dir.path(), "XK2", 40.0);

    let run = poses_with(&["--ligand", "XK2"], &fixture("1HVR"), &[&pose]);
    let found = records(&run);
    assert_eq!(found.len(), 1, "{found:?}");
    assert_eq!(found[0]["copy"], "A/263", "{}", found[0]);
    assert_eq!(found[0]["similarity"], 1.0, "{}", found[0]);
    assert_eq!(found[0]["rmsd"], 0.0, "{}", found[0]);
    assert_eq!(found[0]["clashes"], 0, "{}", found[0]);
    assert!(!run.stderr.contains("warning"), "{}", run.stderr);
}

#[test]
fn each_copy_of_the_ligand_is_paired_with_the_chains_around_it() {
    // A second complex 60 A away, as chains C and D, written after the first,
    // whose inhibitor is nudged 1 A: a copy matches exactly only when paired
    // through C and D.
    let dir = tempfile::tempdir().unwrap();
    let pose = with_another_complex("1HVR", dir.path(), &[("A", "C"), ("B", "D")], "XK2", 1.0);

    let found = records(&poses_with(
        &["--ligand", "XK2"],
        &fixture("1HVR"),
        &[&pose],
    ));
    assert_eq!(found[0]["copy"], "Y/263", "{}", found[0]);
    assert_eq!(found[0]["similarity"], 1.0, "{}", found[0]);
    assert_eq!(found[0]["rmsd"], 0.0, "{}", found[0]);
    assert_eq!(found[0]["clashes"], 0, "{}", found[0]);
}

#[test]
fn copies_named_otherwise_than_the_reference_ligand_are_still_told_apart() {
    // Output from a SMILES input names the ligand as it likes.
    let dir = tempfile::tempdir().unwrap();
    let pose = with_another_copy("1HVR", dir.path(), "XK2", 40.0);
    let text = std::fs::read_to_string(&pose).unwrap();
    let renamed: Vec<String> = text
        .lines()
        .map(|line| match line.starts_with("HETATM") {
            true => line.replace(" XK2 ", " LIG "),
            false => line.to_string(),
        })
        .collect();
    std::fs::write(&pose, renamed.join("\n")).unwrap();

    let apart = ["--select", "resn LIG", "--reference-select", "resn XK2"];
    let found = records(&poses_with(&apart, &fixture("1HVR"), &[&pose]));
    assert_eq!(found[0]["copy"], "A/263", "{}", found[0]);
    assert_eq!(found[0]["rmsd"], 0.0, "{}", found[0]);
}

#[test]
fn a_ligand_in_pieces_that_are_not_all_copies_is_scored_as_one_with_a_warning() {
    let dir = tempfile::tempdir().unwrap();
    let pose = with_another_copy("1HVR", dir.path(), "XK2", 40.0);
    let mixed = [
        "--select",
        "resn XK2 or (chain A and resi 25)",
        "--reference-select",
        "resn XK2",
    ];

    let run = poses_with(&mixed, &fixture("1HVR"), &[&pose]);
    let found = records(&run);
    assert_eq!(found[0]["copy"], serde_json::Value::Null, "{}", found[0]);
    assert!(
        run.stderr.contains("is 3 separate pieces"),
        "{}",
        run.stderr
    );
}

#[test]
fn a_reference_ligand_in_several_pieces_is_compared_whole() {
    let dir = tempfile::tempdir().unwrap();
    let both = with_another_copy("1HVR", dir.path(), "XK2", 40.0);

    let found = records(&poses_with(&["--ligand", "XK2"], &both, &[&both]));
    assert_eq!(found[0]["copy"], serde_json::Value::Null, "{}", found[0]);
    assert_eq!(found[0]["similarity"], 1.0, "{}", found[0]);
    assert_eq!(found[0]["rmsd"], 0.0, "{}", found[0]);
}

#[test]
fn a_ligand_in_a_chain_of_its_own_keeps_its_covalent_attachment() {
    // Co-folding output writes the ligand as its own chain, numbered from 1.
    let dir = tempfile::tempdir().unwrap();
    let rechained = edited_fixture("1PPB", dir.path(), |fields, column| {
        if fields[column["auth_comp_id"]] == "0G6" {
            fields[column["auth_asym_id"]] = "X".to_string();
            fields[column["auth_seq_id"]] = "1".to_string();
        }
    });

    let found = records(&poses_against(&fixture("1PPB"), &[&rechained]));
    assert_eq!(found[0]["similarity"], 1.0, "{}", found[0]);
    assert_eq!(found[0]["clashes"], 0, "{}", found[0]);
}

/// An SD file in `dir` with a record per offset in `shifts`: the heavy atoms
/// of the `ligand` component of the mmCIF file `source`, listed last first as
/// a docking program orders them by its own lights, moved along x by the
/// offset, with the bonds their distances imply.
fn docked_sdf(source: &Path, ligand: &str, dir: &Path, shifts: &[f32]) -> PathBuf {
    let text = std::fs::read_to_string(source).unwrap();
    let column = atom_site_columns(&text);
    let mut atoms: Vec<(String, [f32; 3])> = text
        .lines()
        .filter(|line| line.starts_with("HETATM"))
        .map(|line| line.split_whitespace().collect::<Vec<&str>>())
        .filter(|fields| {
            fields[column["auth_comp_id"]] == ligand
                && fields[column["type_symbol"]] != "H"
                && matches!(fields[column["label_alt_id"]], "." | "A")
        })
        .map(|fields| {
            let position = ["Cartn_x", "Cartn_y", "Cartn_z"]
                .map(|axis| fields[column[axis]].parse::<f32>().unwrap());
            (fields[column["type_symbol"]].to_string(), position)
        })
        .collect();
    atoms.reverse();
    let mut bonds = Vec::new();
    for (a, (_, first)) in atoms.iter().enumerate() {
        for (b, (_, second)) in atoms.iter().enumerate().skip(a + 1) {
            let squared: f32 = (0..3)
                .map(|axis| (first[axis] - second[axis]).powi(2))
                .sum();
            if squared < 1.9 * 1.9 {
                bonds.push((a + 1, b + 1));
            }
        }
    }

    let mut out = String::new();
    for shift in shifts {
        out.push_str(&format!(
            "{ligand}\n  docked\n\n{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
            atoms.len(),
            bonds.len()
        ));
        for (element, [x, y, z]) in &atoms {
            out.push_str(&format!(
                "{:>10.4}{y:>10.4}{z:>10.4} {element:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
                x + shift
            ));
        }
        for (a, b) in &bonds {
            out.push_str(&format!("{a:>3}{b:>3}  1  0\n"));
        }
        out.push_str("M  END\n$$$$\n");
    }

    let path = dir.join(format!("{ligand}-docked.sdf"));
    std::fs::write(&path, out).unwrap();
    path
}

#[test]
fn a_docked_copy_of_the_reference_ligand_scores_as_the_reference_itself_does() {
    // The receptor is the reference with its ligand left out, so each record
    // sits where the reference's ligand did, the second moved 2 A.
    let dir = tempfile::tempdir().unwrap();
    let docked = docked_sdf(&fixture("1PPB"), "0G6", dir.path(), &[0.0, 2.0]);

    let found = records(&poses_against(
        &fixture("1PPB"),
        &[&fixture("1PPB"), &docked],
    ));
    let itself = found
        .iter()
        .find(|record| record["pose"].as_str().unwrap().ends_with(".cif"))
        .unwrap();
    let model = |number: u64| {
        found
            .iter()
            .find(|record| {
                record["pose"].as_str().unwrap().ends_with(".sdf") && record["model"] == number
            })
            .unwrap()
    };
    for key in [
        "contacts",
        "similarity",
        "rmsd",
        "clashes",
        "inverted",
        "lost",
    ] {
        assert_eq!(model(1)[key], itself[key], "{key}: {}", model(1));
    }
    assert_eq!(model(1)["rmsd"], 0.0, "{}", model(1));
    assert!(model(1).get("docked").is_none(), "{}", model(1));
    let moved = model(2)["rmsd"].as_f64().unwrap();
    assert!((1.5..=2.0).contains(&moved), "{}", model(2));
}

#[test]
fn a_docked_pose_in_a_receptor_of_its_own_is_laid_on_the_reference_pocket() {
    // The receptor, turned into a frame of its own, still holds its copy of
    // the ligand, which is left out.
    // Its definition of the ligand calls every nitrogen a carbon; the docked
    // ligand is named as the reference's, so the reference's definition holds.
    let dir = tempfile::tempdir().unwrap();
    let turned = misdefined(&turned_fixture("1PPB", dir.path()), "0G6");
    let docked = docked_sdf(&turned, "0G6", dir.path(), &[0.0]);

    let run = poses_with(
        &["--ligand", "0G6", "--receptor", turned.to_str().unwrap()],
        &fixture("1PPB"),
        &[&docked],
    );
    let found = records(&run);
    assert_eq!(found[0]["rmsd"], 0.0, "{}", found[0]);
    assert_eq!(found[0]["similarity"], 1.0, "{}", found[0]);
    assert!(run.stderr.contains("0G6"), "{}", run.stderr);
}

#[test]
fn a_receptor_needs_an_sd_file_pose_and_a_model_its_record() {
    let dir = tempfile::tempdir().unwrap();
    let reference = fixture("1PPB");

    let stray = poses_with(
        &["--ligand", "0G6", "--receptor", reference.to_str().unwrap()],
        &reference,
        &[&reference],
    );
    assert!(
        !stray.ok && stray.stderr.contains("--receptor"),
        "{}",
        stray.stderr
    );

    let docked = docked_sdf(&reference, "0G6", dir.path(), &[0.0, 2.0]);
    let beyond = poses_with(&["--ligand", "0G6", "--model", "3"], &reference, &[&docked]);
    assert!(
        !beyond.ok && beyond.stderr.contains("2 records"),
        "{}",
        beyond.stderr
    );
}

/// A copy of the mmCIF file `path` beside it whose definition of `component`
/// gives each of its nitrogens as carbon.
fn misdefined(path: &Path, component: &str) -> PathBuf {
    let text = std::fs::read_to_string(path).unwrap();
    let mut in_atoms = false;
    let lines: Vec<String> = text
        .lines()
        .map(|line| {
            if line.starts_with("_chem_comp_atom.") {
                in_atoms = true;
            } else if in_atoms && (line.starts_with('#') || line.starts_with("loop_")) {
                in_atoms = false;
            }
            let mut fields: Vec<&str> = line.split_whitespace().collect();
            if in_atoms && fields.first() == Some(&component) && fields.get(2) == Some(&"N") {
                fields[2] = "C";
                return fields.join(" ");
            }
            line.to_string()
        })
        .collect();

    let copy = path.with_file_name(format!(
        "{}-misdefined.cif",
        path.file_stem().unwrap().to_string_lossy()
    ));
    std::fs::write(&copy, lines.join("\n")).unwrap();
    copy
}

/// A cache directory holding the committed Chemical Component Dictionary
/// entries, so a run fetches nothing.
fn seeded_cache() -> tempfile::TempDir {
    let cache = tempfile::tempdir().unwrap();
    let entries = cache.path().join("ccd");
    std::fs::create_dir_all(&entries).unwrap();
    for entry in std::fs::read_dir(fixture_file("ccd")).unwrap() {
        let path = entry.unwrap().path();
        std::fs::copy(&path, entries.join(path.file_name().unwrap())).unwrap();
    }
    cache
}

/// Run `poses` without a reference, as JSON, against the dictionary entries
/// in `cache`.
fn checked(cache: &Path, args: &[&str], poses: &[&Path]) -> Run {
    let head = ["poses", "--format", "json", "--cache-dir"]
        .into_iter()
        .map(OsStr::new)
        .chain([cache.as_os_str()])
        .chain(args.iter().map(OsStr::new));
    run_args(head.chain(poses.iter().map(|pose| pose.as_os_str())))
}

/// A ligand atom's name, element and position.
type LigandAtom = (String, String, [f32; 3]);

/// The heavy atoms of `ligand` in fixture `id`, and the bonds its definition
/// gives between them, each with its SD bond type.
fn ligand_of(id: &str, ligand: &str) -> (Vec<LigandAtom>, Vec<(usize, usize, u8)>) {
    let text = std::fs::read_to_string(fixture(id)).unwrap();
    let column = atom_site_columns(&text);
    let atoms: Vec<LigandAtom> = text
        .lines()
        .filter(|line| line.starts_with("HETATM"))
        .map(|line| line.split_whitespace().collect::<Vec<&str>>())
        .filter(|fields| {
            fields[column["auth_comp_id"]] == ligand
                && fields[column["type_symbol"]] != "H"
                && matches!(fields[column["label_alt_id"]], "." | "A")
        })
        .map(|fields| {
            let position = ["Cartn_x", "Cartn_y", "Cartn_z"]
                .map(|axis| fields[column[axis]].parse::<f32>().unwrap());
            (
                fields[column["auth_atom_id"]].to_string(),
                fields[column["type_symbol"]].to_string(),
                position,
            )
        })
        .collect();
    let index: HashMap<&str, usize> = atoms
        .iter()
        .enumerate()
        .map(|(at, (name, _, _))| (name.as_str(), at))
        .collect();
    let bonds = text
        .lines()
        .map(|line| line.split_whitespace().collect::<Vec<&str>>())
        .filter(|fields| fields.len() >= 4 && fields[0] == ligand)
        .filter_map(|fields| {
            let (a, b) = (*index.get(fields[1])?, *index.get(fields[2])?);
            let order = match fields[3] {
                "DOUB" => 2,
                "TRIP" => 3,
                "AROM" => 4,
                _ => 1,
            };
            Some((a, b, order))
        })
        .collect();

    (atoms, bonds)
}

/// The centre of the heavy atoms given.
fn centre(atoms: &[LigandAtom]) -> [f32; 3] {
    let mut sum = [0.0; 3];
    for (_, _, position) in atoms {
        for axis in 0..3 {
            sum[axis] += position[axis] / atoms.len() as f32;
        }
    }
    sum
}

/// An SD file of `ligand` from fixture `id`, one record as the crystal holds
/// it, then, when `mirrored`, one reflected through the plane across its
/// centre normal to x, which inverts every stereocentre.
fn ligand_sdf(id: &str, ligand: &str, dir: &Path, mirrored: bool) -> PathBuf {
    let (atoms, bonds) = ligand_of(id, ligand);
    let middle = centre(&atoms)[0];
    let mut out = String::new();
    for reflect in [false, true].into_iter().take(1 + usize::from(mirrored)) {
        out.push_str(&format!(
            "{ligand}\n  docked\n\n{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
            atoms.len(),
            bonds.len()
        ));
        for (_, element, [x, y, z]) in &atoms {
            let x = if reflect { 2.0 * middle - x } else { *x };
            out.push_str(&format!(
                "{x:>10.4}{y:>10.4}{z:>10.4} {element:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n"
            ));
        }
        for (a, b, order) in &bonds {
            out.push_str(&format!("{:>3}{:>3}{order:>3}  0\n", a + 1, b + 1));
        }
        out.push_str("M  END\n$$$$\n");
    }

    let path = dir.join(format!("{ligand}-docked.sdf"));
    std::fs::write(&path, out).unwrap();
    path
}

/// A copy of fixture `id` named `name` in `dir`, its `ligand` reflected
/// through the plane across its centre normal to x.
fn with_mirrored_ligand(id: &str, ligand: &str, dir: &Path) -> PathBuf {
    let middle = centre(&ligand_of(id, ligand).0)[0];
    let own = dir.join("mirrored");
    std::fs::create_dir_all(&own).unwrap();
    edited_fixture(id, &own, |fields, column| {
        if fields[column["auth_comp_id"]] == ligand {
            let x: f32 = fields[column["Cartn_x"]].parse().unwrap();
            fields[column["Cartn_x"]] = format!("{:.3}", 2.0 * middle - x);
        }
    })
}

fn names_of(value: &serde_json::Value) -> Vec<String> {
    let mut names: Vec<String> = value
        .as_array()
        .expect("a list")
        .iter()
        .map(|name| name.as_str().unwrap().to_string())
        .collect();
    names.sort();
    names
}

#[test]
fn a_crystal_ligand_checked_without_a_reference_passes_every_check() {
    let cache = seeded_cache();
    let run = checked(cache.path(), &["--ligand", "XK2"], &[&fixture("1HVR")]);
    let records = records(&run);

    assert_eq!(records.len(), 1);
    let pose = &records[0];
    assert_eq!(pose["failed"], serde_json::json!([]), "{pose}");
    assert_eq!(pose["unchecked"], serde_json::json!([]));
    assert_eq!(pose["stereo"], serde_json::json!([]));
    assert_eq!(pose["clashes"], 0);
    assert_eq!(pose["geometry"]["bonds"], 0);
    assert_eq!(pose["geometry"]["unflat"], 0);
    assert!(pose["contacts"].as_u64().unwrap() > 5, "{pose}");
    assert!(pose["confidence"].is_null());
    for field in ["rmsd", "recovery", "similarity", "lost", "copy"] {
        assert!(pose.get(field).is_none(), "{field} needs a reference");
    }
}

#[test]
fn a_mirrored_ligand_fails_its_stereochemistry_and_ranks_after_the_crystal() {
    let cache = seeded_cache();
    let dir = tempfile::tempdir().unwrap();
    let mirrored = with_mirrored_ligand("1HVR", "XK2", dir.path());
    let crystal = fixture("1HVR");
    let run = checked(cache.path(), &["--ligand", "XK2"], &[&mirrored, &crystal]);
    let records = records(&run);

    assert_eq!(records[0]["pose"], crystal.to_str().unwrap());
    let flipped = &records[1];
    assert!(
        names_of(&flipped["failed"]).contains(&"stereo".to_string()),
        "{flipped}"
    );
    assert_eq!(
        names_of(&flipped["stereo"]),
        [
            "XK2 A/263 C3",
            "XK2 A/263 C4",
            "XK2 A/263 C5",
            "XK2 A/263 C6"
        ]
    );
}

#[test]
fn a_ligand_named_as_a_dictionary_component_it_is_not_has_its_geometry_unchecked() {
    // XK2 renamed LIG, as co-folding output names a ligand given as SMILES,
    // so the file defines no LIG: the dictionary's LIG is another molecule,
    // so nothing is borrowed.
    let cache = seeded_cache();
    let dir = tempfile::tempdir().unwrap();
    let renamed = edited_fixture("1HVR", dir.path(), |fields, column| {
        for name in ["label_comp_id", "auth_comp_id"] {
            if fields[column[name]] == "XK2" {
                fields[column[name]] = "LIG".to_string();
            }
        }
    });

    let run = checked(cache.path(), &["--ligand", "LIG"], &[&renamed]);
    let pose = &records(&run)[0];
    assert_eq!(names_of(&pose["unchecked"]), ["geometry", "stereo"]);
    assert!(pose["geometry"].is_null() && pose["stereo"].is_null());
    assert!(run.stderr.contains("LIG"), "{}", run.stderr);
}

#[test]
fn docked_poses_without_a_reference_are_checked_on_their_own_bond_table() {
    let cache = seeded_cache();
    let dir = tempfile::tempdir().unwrap();
    let docked = ligand_sdf("1HVR", "XK2", dir.path(), true);
    let receptor = fixture("1HVR");

    // Named after the dictionary entry, so their stereochemistry is checked;
    // the receptor's own XK2 is left out, as the poses take its place.
    let run = checked(
        cache.path(),
        &["--ligand", "XK2", "--receptor", receptor.to_str().unwrap()],
        &[&docked],
    );
    let records = records(&run);
    assert_eq!(records.len(), 2);
    assert_eq!(records[0]["model"], 1);
    assert_eq!(
        records[0]["failed"],
        serde_json::json!([]),
        "{}",
        records[0]
    );
    assert_eq!(records[1]["model"], 2);
    assert_eq!(names_of(&records[1]["stereo"]).len(), 4);
    assert!(run.stderr.contains("left out"), "{}", run.stderr);

    // Unnamed, the records are checked on their bonds alone.
    let run = checked(
        cache.path(),
        &[
            "--select",
            "ligand",
            "--receptor",
            receptor.to_str().unwrap(),
        ],
        &[&docked],
    );
    let unnamed = self::records(&run);
    assert!(unnamed.iter().all(|pose| pose["stereo"].is_null()));
    assert!(unnamed.iter().all(|pose| !pose["geometry"].is_null()));
    assert!(run.stderr.contains("--ligand"), "{}", run.stderr);
}

#[test]
fn co_folding_confidence_is_read_beside_each_model_and_breaks_ties() {
    let cache = seeded_cache();
    let dir = tempfile::tempdir().unwrap();
    let crystal = std::fs::read(fixture("1HVR")).unwrap();
    // AlphaFold Server's names, and Boltz's.
    let alphafold = dir.path().join("fold_job_model_0.cif");
    std::fs::write(&alphafold, &crystal).unwrap();
    std::fs::write(
        dir.path().join("fold_job_summary_confidences_0.json"),
        r#"{"ranking_score": 0.42, "chain_iptm": [0.31, 0.5], "iptm": 0.4, "ptm": 0.8}"#,
    )
    .unwrap();
    let boltz = dir.path().join("job_model_1.cif");
    std::fs::write(&boltz, &crystal).unwrap();
    std::fs::write(
        dir.path().join("confidence_job_model_1.json"),
        r#"{"confidence_score": 0.9, "ligand_iptm": 0.77, "ptm": 0.8}"#,
    )
    .unwrap();

    let run = checked(cache.path(), &["--ligand", "XK2"], &[&alphafold, &boltz]);
    let records = records(&run);
    assert_eq!(records[0]["pose"], boltz.to_str().unwrap());
    assert_eq!(records[0]["confidence"]["source"], "boltz");
    assert_eq!(records[0]["confidence"]["score"], 0.9);
    assert_eq!(records[0]["confidence"]["ligand_iptm"], 0.77);
    assert_eq!(records[1]["confidence"]["source"], "alphafold3");
    assert_eq!(records[1]["confidence"]["score"], 0.42);
    // AlphaFold 3 lists its chains' ipTM in the model's chain order, and the
    // ligand sits in chain A, the first.
    assert_eq!(records[1]["confidence"]["ligand_iptm"], 0.31);
}

#[test]
fn a_ligand_code_in_lower_case_is_the_same_component() {
    let cache = seeded_cache();
    let dir = tempfile::tempdir().unwrap();
    let docked = ligand_sdf("1HVR", "XK2", dir.path(), true);
    let receptor = fixture("1HVR");

    let run = checked(
        cache.path(),
        &["--ligand", "xk2", "--receptor", receptor.to_str().unwrap()],
        &[&docked],
    );
    let records = records(&run);
    assert!(run.stderr.contains("left out"), "{}", run.stderr);
    assert_eq!(records[0]["clashes"], 0);
    assert_eq!(names_of(&records[1]["stereo"]).len(), 4);
}

#[test]
fn options_that_compare_with_a_reference_need_one() {
    let cache = seeded_cache();
    let pose = fixture("1HVR");
    for args in [
        &["--ligand", "XK2", "--reference-select", "resn XK2"][..],
        &["--ligand", "XK2", "--view"][..],
    ] {
        let run = checked(cache.path(), args, &[&pose]);
        assert!(!run.ok, "{args:?} ran");
    }

    let dir = tempfile::tempdir().unwrap();
    let docked = ligand_sdf("1HVR", "XK2", dir.path(), false);
    let run = checked(cache.path(), &["--ligand", "XK2"], &[&docked]);
    assert!(!run.ok);
    assert!(run.stderr.contains("--receptor"), "{}", run.stderr);
}
