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

    for args in [
        vec!["rin", "--model", "all"],
        vec!["render", "--model", "all", "-o", png.to_str().unwrap()],
    ] {
        let run = run_on(&path, &args);
        assert!(!run.ok, "{args:?} should refuse --model all");
        assert!(run.stderr.contains("one model at a time"), "{}", run.stderr);
    }
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
    edit: impl Fn(&mut [String], &HashMap<&str, usize>),
) -> PathBuf {
    let text = std::fs::read_to_string(fixture(id)).unwrap();
    let columns: HashMap<&str, usize> = text
        .lines()
        .filter_map(|line| line.strip_prefix("_atom_site."))
        .enumerate()
        .map(|(index, name)| (name.trim(), index))
        .collect();

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

/// Run `poses` with 1PPB as the reference and its PPACK inhibitor (0G6) as the
/// ligand, written as JSON.
fn poses_against(reference: &Path, poses: &[&Path]) -> Run {
    let head = [
        "poses",
        "--ligand",
        "0G6",
        "--format",
        "json",
        "--reference",
    ]
    .map(OsStr::new);
    run_args(
        head.into_iter()
            .chain([reference.as_os_str()])
            .chain(poses.iter().map(|pose| pose.as_os_str())),
    )
}

fn records(run: &Run) -> Vec<serde_json::Value> {
    assert!(run.ok, "command failed: {}", run.stderr);
    let parsed: serde_json::Value = serde_json::from_str(&run.stdout).expect("valid json");

    parsed.as_array().expect("an array").clone()
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
    assert_eq!(same["lost"], serde_json::json!([]));
    assert!(same["contacts"].as_u64().unwrap() > 0);
    assert_eq!(same["contacts"], same["shared"]);

    let moved = &found[1];
    assert!(moved["recovery"].as_f64().unwrap() < 1.0, "{moved}");
    assert!(moved["clashes"].as_u64().unwrap() > 0, "{moved}");
    assert!(!moved["lost"].as_array().unwrap().is_empty(), "{moved}");
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

#[test]
fn a_pose_numbered_unlike_the_reference_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let shifted = edited_fixture("1PPB", dir.path(), |fields, column| {
        if fields[column["group_PDB"]] == "ATOM" && fields[column["auth_asym_id"]] == "H" {
            let number: i32 = fields[column["auth_seq_id"]].parse().unwrap();
            fields[column["auth_seq_id"]] = (number + 1).to_string();
        }
    });

    let run = poses_against(&fixture("1PPB"), &[&shifted]);
    assert!(!run.ok, "expected a failure, got: {}", run.stdout);
    assert!(
        run.stderr
            .contains("in the reference; a pose must number its residues as the reference does"),
        "unhelpful error: {}",
        run.stderr
    );
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
