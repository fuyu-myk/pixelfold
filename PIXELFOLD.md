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
  - `[x]` `select/` : selection language, `nom` parser to atom-index bitset. (today: `core/src/select/`; `mod.rs` the `Selection` AST + `select(query, protein)` entry + `AtomSet` (a `fixedbitset` over atom indices with union/intersect/complement); `parser.rs` a recursive-descent `nom` parser (or/and/not precedence, parens, `byres`/`within` binding to a primary, reserved keywords); `eval.rs` evaluates to an `AtomSet`, using the spatial grid for `within` and residue keys for `byres`. Predicates: chain/resi/resn/name/element/model/ss/bfactor/plddt/protein/nucleic/ligand/water/metal/hetatm/all/none. Ligand rings and full ligand typing still wait on the CCD)
  - `[~]` `dssp/` : secondary structure + Kabsch-Sander backbone H-bond energies. (today: `core/src/dssp.rs`; chain-aware residue model, amide hydrogen inferred correctly (opposite the preceding carbonyl, skipping proline and chain starts/breaks), Kabsch-Sander energies, and DSSP pattern rules (consecutive n-turns for alpha/3-10/pi helices, parallel/antiparallel bridges for sheets, turns) reduced to helix/sheet/turn/coil. Reproduces ubiquitin's fold; isolated single bridges are still called sheet, and bend (S) is not computed. The H-bond search uses the spatial grid (O(N)); the bridge scan is still O(N^2))
  - `[~]` `sasa/` : Shrake-Rupley, wraps `rust-sasa`. (today: `core/src/sasa.rs`)
  - `[~]` `interactions/` : all interaction types (H-bond, salt bridge, hydrophobic, pi-stacking, pi-cation, halogen, water bridge, metal, disulfide) with angular terms. (today: `core/src/interactions/`; `params.rs` holds every threshold taken from PLIP's `config.py` with the PLIP variable named and the RING version and difference recorded per value; `topology.rs` tabulates the heavy-atom bonds of the 20 standard residues, which keeps atom roles derived rather than transcribed; `classify.rs` holds the atom-role rules: residue tables at the pH 7.4 protonation model for the standard residues, and for everything else the chemistry the file states about its own components, so a ligand's apolar carbons and H-bond acceptors are typed by OpenBabel's own rules (nitro, sulfone, ester and aryl-ether oxygens; nitrogen by degree and hybridisation decoded from bond orders and aromatic flags); `hydrogens.rs` infers polar hydrogens, placing them once where the heavy atoms fix them and against the acceptor under test where a rotation is free; `aromatic.rs` perceives the aromatic rings of Phe/Tyr/His/Trp and reduces each to a centre and a Newell-fitted plane normal; `connectivity.rs` builds the bond graph, joining the residue table, the file's own component definitions, and the peptide bond; `detectors.rs` implements hydrogen bond, disulfide, metal coordination, salt bridge, pi-stacking, pi-cation, hydrophobic contact, water bridge, and halogen bond over the spatial grid. `Interaction` carries an optional `bridge` atom, the mediating water of a water bridge. A ligand now contributes apolar carbons and H-bond acceptors; it does not yet donate a hydrogen bond nor carry a perceived ring)
  - `[~]` `rin/` : petgraph residue interaction network + centrality + motifs. (today: H-bond only in `core/src/rin.rs`)
- `[~]` `crates/pixelfold-render/` : scene to RGBA framebuffer to sinks
  - `[~]` `raster/` : (today: `render/src/renderer.rs` orthographic projection + camera, plus `render/src/draw.rs` Bresenham line + b-factor/hbond/hydrophobicity color mapping; z-buffer, sphere impostors, depth cue, slab, ID buffer are planned)
  - `[ ]` `sink/` : kitty, sixel, iterm2, unicode, png
- `[~]` `crates/pixelfold-fetch/` : RCSB search + structure download (addition beyond roadmap layout; used by both the CLI `fetch` subcommand and the TUI search mode). (today: `client.rs` async reqwest wrapper, `download.rs` concurrent gzip download for the search TUI, `types.rs` API response types, and `fetch_cif` a blocking single-structure download for the CLI resolver)
- `[~]` `crates/pixelfold-cli/` : headless subcommands (the product): view, rin, interactions, sasa, ss, render, fetch, validate. (today: `clap` arg parsing, a path/PDB-id resolver that auto-fetches uncached ids into the XDG cache, then dispatches to `pixelfold_tui::view`/`search`. Named subcommands still to come)
- `[~]` `crates/pixelfold-tui/` : ratatui front-end (the demo), now a library. (today: `args`/`app`/`inputs`/`ui` split the entry + event loop + input + render; `lib.rs` `App` state; `search/` search mode)
- `[~]` `crates/pixelfold-validate/` : validation harness, precision/recall vs DSSP, FreeSASA, PLIP, RING. (today: compares pixelfold's DSSP (Q3, H-bond edge F1) and SASA (MAE/median/Pearson) against committed golden files and prints a markdown report; PLIP/RING interaction validation waits on the interaction engine; golden files themselves come from the dockerised reference tools)
- `[~]` `benchmarks/` : pinned PDB manifest + golden reference files. (today: `manifest.toml` with a confident starter set per stratum (expand toward the target counts); `golden/` for `<ID>.dssp.json` / `<ID>.freesasa.json`; `tools/` a pinned DSSP 4 + FreeSASA Docker image + `generate_golden.py` that produces the golden set; `thresholds.toml` the `--check` regression gates (all off until baselines exist); structures cache under `benchmarks/cache/`, gitignored. A CI job runs the harness once golden files land)

Two binaries now exist (`pixelfold` from cli, `pixelfold-validate`); `default-members = ["crates/pixelfold-cli"]` keeps bare `cargo run` pointed at `pixelfold`.

### Current source map

`pixelfold-core/src/`:

- `structure.rs` : domain types (`Atom`, `Protein`, `BiologicalAssembly`, `DisplayMode`, `SecondaryStructure`) + pure helpers (b-factor range, C-alpha indices/connections).
- `fixed_str.rs` : `FixedStr<N>`, an inline fixed-capacity ASCII string (Copy, zero-alloc, cheap HashMap key) for hot `Atom` fields; truncates content longer than `N` bytes.
- `spatial.rs` : `Grid`, a uniform-grid spatial index (O(N) build, O(1) `within`/`for_each_within` neighbour queries); cell size = max interaction cutoff.
- `interactions/` : the interaction engine. `mod.rs` (`Interaction`, `InteractionKind`, `detect`), `params.rs` (thresholds from PLIP source, RING differences noted), `topology.rs` (standard-residue heavy-atom bonds), `connectivity.rs` (the whole structure's bond graph, residue table plus file component definitions plus the peptide bond), `classify.rs` (residue atom-role tables), `hydrogens.rs` (polar hydrogen inference, fixed and rotatable), `aromatic.rs` (ring perception, centre and normal), `detectors.rs` (per-type detection over the grid).

  The engine's hydrogen bonds are not DSSP's. `dssp.rs` scores backbone pairs by Kabsch-Sander energy to assign secondary structure; the engine applies PLIP's geometry to backbone and side chains alike, and the two disagree on purpose. On ubiquitin the engine reproduces 83% of DSSP's backbone bonds, and every one it drops is explained: DSSP has no angular term at all, so it admits bonds the geometry rejects, and it lets one amide hydrogen serve two acceptors where the engine allows one bond per hydrogen.
- `parser.rs` : pdbtbx PDB/mmCIF loading (`build_atom`) + altloc policy (`filter_altlocs`) + fixed-capacity truncation warning; calls `dssp::assign` after loading.
- `dssp.rs` : DSSP secondary-structure assignment: chain-aware residue model, amide-hydrogen inference, Kabsch-Sander H-bond energies, and the `HBond` type consumed by the renderer and RIN.
- `components.rs` : `Dictionary` of chemical component definitions, read from the `_chem_comp_atom` / `_chem_comp_bond` categories the structure file carries about its own residues and ligands. Supplies element, aromatic flag, bond order, and the hydrogen count of each heavy atom, which together are what the atom-role rules need. No chemical component dictionary needs downloading; a file without those categories yields an empty dictionary and ligands stay untyped. Note the definitions describe the free component, so a standard residue's hydrogen count is the free amino acid's and not the polymer's: the residue tables stay authoritative for those.
- `mmcif.rs` : the looped and single-record mmCIF readers shared by `assembly.rs` and `components.rs`.
- `assembly.rs` : pure `detect_partial_assembly` reading `_pdbx_struct_assembly` / `REMARK 350` from the raw file (pdbtbx does not parse them), flagging when the biological unit needs symmetry expansion beyond the deposited coordinates. Detection only; generation is later work.
- `sasa.rs` : Shrake-Rupley SASA (wraps `rust-sasa`): surface points for rendering + `calculate_atom_sasa` (per-atom area, for validation) + vdW radius / hydrophobicity tables.
- `rin.rs` : petgraph H-bond residue interaction network + centrality/components/motifs.
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

- `main.rs` : `clap` parser (`--no-surface`, `--cache-dir`, `--altloc`, `--fetch`), XDG cache dir (via `dirs`), fetch-on-miss, dispatch to `pixelfold_tui::view`/`search`
- `resolve.rs` : structure resolver returning a `Resolution` (existing path / cached id / fetch request / error), unit-tested

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
