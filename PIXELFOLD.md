# PIXELFOLD.md

Pixelfold loads `PIXELFOLD.md` from the workspace root as agent memory (similar to AGENTS.md / CLAUDE.md). This file is also the project's living architecture doc.

**Read it before making changes, and update it whenever you change architecture, add a crate/module, or ship a user-facing feature.** Keeping this file accurate is part of "done" for any change.

## Project

**Pixelfold**: open-source, terminal-based protein structure visualization and analysis tool, written in Rust and Ratatui. It features a braille-based rendering engine for high-resolution visualization of protein structures, hydrogen bond networks, and solvent-accessible surfaces.

- Run/build: `cargo build --release`, and `cargo run --release`
- Rust checks: `cargo clippy --workspace --all-targets -- -D warnings` and `cargo test --workspace --all-targets`.

## Quality bar

Production-grade or it does not ship. Every change is judged against all of these, not just "it works":

- **Correctness**: edge cases, failure modes, concurrent access. No "works for now".
- **Performance**: speed is the product. For every change ask: how much time it takes, whether it adds unnecessary computations or allocations, whether it triggers extra re-renders or wasted work.
- **UI/UX (TUI)**: polished, professional, premium. Every state and detail considered.
- **Architecture**: new or changed logic lives in pure, dependency-light functions (functional core). Keeps it testable without a later rewrite.

Verify before claiming done: `cargo clippy`, `cargo test`. A change to a core subsystem (main TUI render, analysis calculations) needs a test that locks the invariant.

## Conventions

- **Comments**: default to none, the code should explain itself. If genuinely needed, 1-2 lines on *why*, never *what*. No AI-generic filler.
  - Avoid the use of "we" in comments, be explicit and direct: "this runs on startup" instead of "we run this on startup".
- **No em-dash** anywhere: code, comments, commits, docs.
- **No emojis** anywhere.
- **References to plans**: references to phases in the implementation plan/roadmap should be avoided in code/comments/docs.
- **Immutability**: every state transition returns a new value; never mutate in place. (`Tab` is `Clone + PartialEq` but not `Eq` because it carries an `f64` zoom.)

## Architecture

Living map of where code lives and where it is going. This doubles as the progress tracker: keep the status markers current as work lands. Status: `[x]` done, `[~]` partial, `[ ]` planned.

### Target: Cargo workspace

The Cargo workspace exists. `pixelfold-core` holds the extracted analysis code (`structure`, `parser`, `sasa`, `rin`); `pixelfold-render` holds camera/projection (`renderer`) and the color/line helpers (`draw`); `pixelfold-fetch` holds the RCSB client + download + API types; `pixelfold-tui` is now a library holding the UI, `App` state, screen-space picking, and the app driver (`run`); `pixelfold-cli` is the `pixelfold` binary that calls `pixelfold_tui::run()`; `validate` is still a stub awaiting its code. The direction is that the library has lasting value independent of the terminal UI. Rules:

- `core` carries zero terminal or rendering dependencies; it is usable as a plain library.
- `core` uses fixed-size byte fields in hot structs, not `String` (atom names `[u8; 4]`, residue names `[u8; 3]`, chain ids `[u8; 2]`).
- Everything the TUI does, the CLI can do headless. A TUI-only feature is in the wrong crate.

Crates:

- `[ ]` `crates/pixelfold-core/` : the library (pure, dependency-light)
  - `[~]` `structure/` : Atom/Residue/Chain/Model + mmCIF/PDB ingest. (today: `core/src/structure.rs` holds `Atom`/`Protein`; ingest in `core/src/parser.rs`; `Atom` carries element/is_hetatm/altloc/occupancy/insertion_code/model_number; hot string fields are fixed-size byte arrays via `FixedStr<N>` (name 4, residue 6, chain 4, element 2) in `core/src/fixed_str.rs`; the parser warns once per load if any identifier overflows and is truncated; no first-class Residue/Chain; entity_id not exposed by pdbtbx; single-model only)
  - `[~]` `spatial/` : uniform grid / cell list, cell size = max interaction cutoff (~6 A). O(N) build, O(1) neighbour queries; powers interactions, `within` selections, SASA neighbours, clash detection. (today: `core/src/spatial.rs` `Grid` (`build` + `within` / `for_each_within`); consumed by DSSP's H-bond search; interactions/selections/SASA still to adopt it)
  - `[x]` `select/` : selection language, `nom` parser to atom-index bitset. (today: `core/src/select/`; `mod.rs` the `Selection` AST + `select(query, protein)` entry + `AtomSet` (a `fixedbitset` over atom indices with union/intersect/complement); `parser.rs` a recursive-descent `nom` parser (or/and/not precedence, parens, `byres`/`within` binding to a primary, reserved keywords); `eval.rs` evaluates to an `AtomSet`, using the spatial grid for `within` and residue keys for `byres`. Predicates: chain/resi/resn/name/element/model/ss/bfactor/plddt/protein/nucleic/ligand/water/metal/hetatm/all/none)
  - `[~]` `dssp/` : secondary structure + Kabsch-Sander backbone H-bond energies. (today: `core/src/dssp.rs`; chain-aware residue model, amide hydrogen inferred correctly (opposite the preceding carbonyl, skipping proline and chain starts/breaks), Kabsch-Sander energies, and DSSP pattern rules (consecutive n-turns for alpha/3-10/pi helices, parallel/antiparallel bridges for sheets, turns) reduced to helix/sheet/turn/coil. Reproduces ubiquitin's fold; isolated single bridges are still called sheet, and bend (S) is not computed. The H-bond search uses the spatial grid (O(N)); the bridge scan is still O(N^2))
  - `[~]` `sasa/` : Shrake-Rupley, wraps `rust-sasa`. (today: `core/src/sasa.rs`)
  - `[~]` `interactions/` : all interaction types (H-bond, salt bridge, hydrophobic, pi-stacking, pi-cation, halogen, water bridge, metal, disulfide) with angular terms. (today: `core/src/interactions/`; `params.rs` holds every threshold taken from PLIP's `config.py` with the PLIP variable named and the RING version and difference recorded per value; `topology.rs` tabulates the heavy-atom bonds of the 20 standard residues, which keeps atom roles derived rather than transcribed; `classify.rs` holds the atom-role rules: residue tables at the pH 7.4 protonation model for the standard residues, and for everything else the chemistry the file states about its own components, so a ligand's apolar carbons and H-bond acceptors are typed by OpenBabel's own rules (nitro, sulfone, ester and aryl-ether oxygens; nitrogen by degree and hybridisation decoded from bond orders and aromatic flags); `hydrogens.rs` infers polar hydrogens, placing them once where the heavy atoms fix them and against the acceptor under test where a rotation is free, from the residue tables for a standard residue and from the component definition otherwise. A definition states the free molecule, so two corrections make it the polymer's: every bond to a neighbouring component accounts for one of the free component's hydrogens, and the pH 7.4 protonation model already applied to Asp and Glu extends to any carboxylate or phosphate. The converse correction is deliberately not made, a neighbour absent by chemistry being indistinguishable from one absent through disorder, which costs the hydroxyl at a chain's 5' end; `aromatic.rs` perceives rings, tabulated for Phe/Tyr/His/Trp and found as five- and six-membered cycles over the atoms the file flags aromatic for everything else, each reduced to a centre and a Newell-fitted plane normal; `connectivity.rs` builds the bond graph, joining the residue table, the file's own component definitions, and every linkage between components: the peptide bond by name, and the rest (phosphodiester, glycosidic, covalent ligand) by covalent radii from the coordinates, since neither source describes a component's bonds to its neighbours. Water and metals are excluded from that perception, a coordination contact and a solvent atom modelled over a side chain both falling in the same distance range without being bonds; `detectors.rs` implements hydrogen bond, disulfide, metal coordination, salt bridge, pi-stacking, pi-cation, hydrophobic contact, water bridge, and halogen bond over the spatial grid, with `detect` building the bond graph once and sharing it. `Interaction` carries an optional `bridge` atom, the mediating water of a water bridge. A ligand, a nucleotide, and a modified residue all contribute apolar carbons, H-bond donors and acceptors, and rings; what a component still lacks is a formal charge, so it takes no part in salt bridges or pi-cation)
  - `[~]` `rin/` : petgraph residue interaction network + centrality + motifs. (today: `core/src/rin/`; `hbond.rs` the backbone H-bond graph the DSSP pass feeds, used by the viewer overlay, with degree centrality, connected components, and motif detection; `network.rs` the residue interaction network proper, built over the whole interaction engine: one node per residue, one edge per residue pair per interaction kind, the atom-level interactions aggregated to the residue scale, undirected, with per-node degree. Betweenness and RSA node attributes still to come with the network view)
- `[~]` `crates/pixelfold-render/` : scene to RGBA framebuffer to sinks
  - `[~]` `raster/` : (today: `render/src/renderer.rs` orthographic projection + camera, plus `render/src/draw.rs` Bresenham line + b-factor/hbond/hydrophobicity color mapping; z-buffer, sphere impostors, depth cue, slab, ID buffer are planned)
  - `[ ]` `sink/` : kitty, sixel, iterm2, unicode, png
- `[~]` `crates/pixelfold-fetch/` : RCSB search + structure download (addition beyond roadmap layout; used by both the CLI `fetch` subcommand and the TUI search mode). (today: `client.rs` async reqwest wrapper, `download.rs` concurrent gzip download for the search TUI, `types.rs` API response types, and `fetch_cif` a blocking single-structure download for the CLI resolver)
- `[~]` `crates/pixelfold-cli/` : headless subcommands (the product): view, rin, interactions, sasa, ss, render, fetch, validate. (today: `clap` subcommands `view`/`interactions`/`ss`/`sasa`/`rin`/`fetch`, each over a path/PDB-id resolver that auto-fetches uncached ids into the XDG cache. Every analysis takes `--select` (the selection language); the flat reports take `--format table|tsv|json` and `interactions` additionally `--type` per kind and `--ligand CODE`; `rin` takes `--type`, `--format json|graphml|tsv` (node-link JSON, GraphML for Cytoscape/Gephi, or a flat edge list), and `-o` to write a file. A bare `pixelfold <structure>` still opens the viewer. `render` still to come, with the renderer)
- `[~]` `crates/pixelfold-tui/` : ratatui front-end (the demo), now a library. (today: `args`/`app`/`inputs`/`ui` split the entry + event loop + input + render; `lib.rs` `App` state; `search/` search mode)
- `[~]` `crates/pixelfold-validate/` : validation harness, precision/recall vs DSSP, FreeSASA, PLIP, RING. (today: compares pixelfold's DSSP (Q3, H-bond edge F1) and SASA (MAE/median/Pearson) against committed golden files and prints a markdown report; PLIP/RING interaction validation waits on the interaction engine; golden files themselves come from the dockerised reference tools)
- `[~]` `benchmarks/` : pinned PDB manifest + golden reference files. (today: `manifest.toml` with a confident starter set per stratum (expand toward the target counts); `golden/` for `<ID>.dssp.json` / `<ID>.freesasa.json`; `tools/` a pinned DSSP 4 + FreeSASA Docker image + `generate_golden.py` that produces the golden set; `thresholds.toml` the `--check` regression gates (all off until baselines exist); structures cache under `benchmarks/cache/`, gitignored. A CI job runs the harness once golden files land)

Two binaries now exist (`pixelfold` from cli, `pixelfold-validate`); `default-members = ["crates/pixelfold-cli"]` keeps bare `cargo run` pointed at `pixelfold`.

### Current source map

`pixelfold-core/src/`:

- `structure.rs` : domain types (`Atom`, `Protein`, `BiologicalAssembly`, `DisplayMode`, `SecondaryStructure`) + pure helpers (b-factor range, C-alpha indices/connections).
- `fixed_str.rs` : `FixedStr<N>`, an inline fixed-capacity ASCII string (Copy, zero-alloc, cheap HashMap key) for hot `Atom` fields; truncates content longer than `N` bytes.
- `spatial.rs` : `Grid`, a uniform-grid spatial index (O(N) build, O(1) `within`/`for_each_within` neighbour queries); cell size = max interaction cutoff.
- `interactions/` : the interaction engine. `mod.rs` (`Interaction`, `InteractionKind`, `detect`), `params.rs` (thresholds from PLIP source, RING differences noted), `topology.rs` (standard-residue heavy-atom bonds), `connectivity.rs` (the whole structure's bond graph, residue table plus file component definitions plus the peptide bond), `classify.rs` (residue atom-role tables, plus the component-driven rules for everything else), `hydrogens.rs` (polar hydrogen inference, fixed, rotatable, and undetermined between sites), `aromatic.rs` (ring perception, centre and normal), `detectors.rs` (per-type detection over the grid).

  The engine's hydrogen bonds are not DSSP's. `dssp.rs` scores backbone pairs by Kabsch-Sander energy to assign secondary structure; the engine applies PLIP's geometry to backbone and side chains alike, and the two disagree on purpose. On ubiquitin the engine reproduces 83% of DSSP's backbone bonds, and every one it drops is explained: DSSP has no angular term at all, so it admits bonds the geometry rejects, and it lets one amide hydrogen serve two acceptors where the engine allows one bond per hydrogen.
- `parser.rs` : pdbtbx PDB/mmCIF loading (`build_atom`) + altloc policy (`filter_altlocs`) + fixed-capacity truncation warning; calls `dssp::assign` after loading.
- `dssp.rs` : DSSP secondary-structure assignment: chain-aware residue model, amide-hydrogen inference, Kabsch-Sander H-bond energies, and the `HBond` type consumed by the renderer and RIN.
- `components.rs` : `Dictionary` of chemical component definitions, read from the `_chem_comp_atom` / `_chem_comp_bond` categories the structure file carries about its own residues and ligands. Supplies element, aromatic flag, bond order, and the hydrogen count of each heavy atom, which together are what the atom-role rules need. No chemical component dictionary needs downloading; a file without those categories yields an empty dictionary and ligands stay untyped. Note the definitions describe the free component, so a standard residue's hydrogen count is the free amino acid's and not the polymer's: the residue tables stay authoritative for those.
- `mmcif.rs` : the looped and single-record mmCIF readers shared by `assembly.rs` and `components.rs`.
- `assembly.rs` : pure `detect_partial_assembly` reading `_pdbx_struct_assembly` / `REMARK 350` from the raw file (pdbtbx does not parse them), flagging when the biological unit needs symmetry expansion beyond the deposited coordinates. Detection only; generation is later work.
- `sasa.rs` : Shrake-Rupley SASA (wraps `rust-sasa`): surface points for rendering + `calculate_atom_sasa` (per-atom area, for validation) + vdW radius / hydrophobicity tables.
- `rin/` : residue interaction networks. `hbond.rs` the backbone H-bond graph (centrality, components, motifs) the viewer overlay uses; `network.rs` the interaction network over the whole engine (residue nodes, typed edges aggregated to the residue scale).
- `select/` : the selection language. `mod.rs` (`Selection` AST, `select` entry, `AtomSet` bitset), `parser.rs` (`nom` recursive descent), `eval.rs` (evaluation to an `AtomSet`, grid-backed `within`).

`pixelfold-render/src/`:

- `renderer.rs` : `Camera` (quaternion orientation, screen-space pan, fit-relative zoom, orthographic projection; recomputed per frame, no cached matrix), `ProjectedAtom`/`ProjectedSurfacePoint`, projection + auto-frame
- `draw.rs` : Bresenham line/dashed-line + b-factor/hbond/hydrophobicity color mapping

`pixelfold-fetch/src/`:

- `lib.rs` : `fetch_cif(id, dir)` blocking single-structure download for the CLI resolver
- `client.rs` : `RCSBClient` async reqwest wrapper (search, entry data, CIF download)
- `download.rs` : `DownloadManager` concurrent gzip download to a caller-supplied dir (search TUI)
- `types.rs` : RCSB API response types (`SearchResponse`, `SearchData`, `SearchResult`, ...)

`pixelfold-cli/src/`:

- `main.rs` : `clap` subcommand definitions + dispatch; XDG cache dir (via `dirs`); the bare-structure viewer path kept for compatibility
- `load.rs` : resolve, fetch on miss, load without the surface (nothing headless renders), and evaluate `--select`; a selection matching nothing is an error rather than an empty report
- `analyse.rs` : what each subcommand computes, as records; an interaction is reported when *either* side is selected, which is what asking for a ligand's interactions means. Distances and angles are rounded once here so every format agrees and none claims more precision than the coordinates carry. Every record carries the insertion code, without which a chymotrypsin-numbered structure emits byte-identical rows for different residues. SASA is computed over the polymer alone (no solvent, no hydrogens), the same atom set `pixelfold-validate` compares against FreeSASA: counting the ordered waters occludes the surface being measured and triples the error against the reference
- `report.rs` : the `Row` trait and the three writers (aligned table, TSV, JSON). An empty result still prints its columns, so nothing-found reads differently from wrong-command; the table drops a column empty in every row, the machine formats keep every column so their shape is stable
- `graph.rs` : the residue interaction network writers: node-link JSON, GraphML (with XML escaping), and a flat edge-list TSV. Core stays serde-free, so the JSON is built from local structs referencing the core `Network`
- `resolve.rs` : structure resolver returning a `Resolution` (existing path / cached id / fetch request / error), unit-tested
- `tests/headless.rs` : the built binary run over a cached structure, covering each format and the selection error paths

`pixelfold-tui/src/`:

- `lib.rs` : mutable `App` state, screen-space picking/highlight, and the public entry points `view(path, skip_surface)` / `search(query, cache_dir)` over a `with_terminal` helper
- `app.rs` : `run_app()` event loop
- `inputs.rs` : `handle_input` (keyboard) + `handle_mouse` (click picking)
- `ui.rs` : `ui()` per-frame render (canvas, info bar, inspect panel, plus a highlighted header line when the loaded coordinates are a partial biological assembly)
- `search/mod.rs` : search-mode TUI (state machine, widgets, thread/runtime orchestration over `pixelfold-fetch`)
- `search/types.rs` : `FoundProtein` UI row model + `PageState` pagination

`pixelfold-validate/src/`:

- `main.rs` : CLI (`--manifest`/`--golden-dir`/`--cache-dir`/`--offline`); loads the manifest, resolves + loads each structure, runs the analyses, compares against golden, prints the report
- `metrics.rs` : pure metrics (Q3 agreement, edge precision/recall/F1, SASA MAE/median/Pearson)
- `golden.rs` : serde schema + loader for the golden JSON files (per-residue SS + H-bonds; per-residue SASA)
- `analysis.rs` : extract pixelfold's per-residue prediction (SS, H-bond edges, per-residue SASA) from a loaded `Protein`
- `compare.rs` : align a prediction to golden by residue and score it
- `report.rs` : aggregate per-entry metrics into a markdown table (per-entry breakdown + summary)
- `check.rs` : `--check` regression gate comparing aggregate metrics to `benchmarks/thresholds.toml`

`benchmarks/`:

- `manifest.toml` : stratified PDB entry list (confident starter set)
- `golden/` : reference outputs, `<ID>.dssp.json` (mkdssp) and `<ID>.freesasa.json` (FreeSASA); PLIP/RING added in the interaction work
- `cache/` : downloaded structures (gitignored)

Structure resolution now goes through the CLI: an existing path is used as given, a 4-character PDB id is looked up in (and downloads land in) the XDG cache dir (`dirs::cache_dir()/pixelfold`, override with `--cache-dir`). The old `env!("CARGO_MANIFEST_DIR")` lookup is gone. `crates/pixelfold-tui/data/` remains only as optional local sample structures (gitignored), reachable by explicit path.

### Load-bearing decisions

- One uniform grid built once at load; every neighbour query goes through it.
- Everything renders through one RGBA framebuffer, so real pixels, headless PNG export, and snapshot tests all come from one path.
- Analysis lives in pure, testable functions (functional core), independent of the render and TUI crates.

The roadmap, milestone sequence, and the full backlog of known correctness/performance debt live in `ROADMAP.md`; the status markers above track which target modules exist yet.
