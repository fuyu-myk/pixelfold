# Changelog

All notable changes to PixelFold are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Interaction agreement against RING 4.0**, the whole-structure reference
  beside PLIP's ligand sites. Its golden files cover all 205 benchmark entries
  and come from the public RING 4.0 web server
  (`benchmarks/tools/generate_ring.py`), since RING's standalone binary may not
  be passed on. Seven RING types map onto pixelfold's; van der Waals contacts,
  pi-hydrogen bonds and generic contacts have no counterpart and are left
  unscored. Pixelfold is scored only where RING looks: edges to water, and
  amino acids of one chain fewer than three positions apart, are left out.
  Per-type F1 runs from 0.32 (pi-stacking, where RING's wider window finds
  seven times as many) to 1.00 (disulfide), with a mean of 0.729, and
  `--check` now fails below 0.65 (`min_ring_f1`).

### Fixed

- **`pixelfold sasa` left out modified residues.** Surface atoms were chosen by
  record type, so a selenomethionine or phosphoserine inside a chain, written as
  a HETATM record, was neither reported nor allowed to shield the residues
  beside it, which read too exposed: 3G5S's 14 MSE were missing, and `--select
  'resn MSE'` printed an empty table. The surface is now every heavy atom of a
  polymer entity (`Atom::is_polymer`); a file without entity records keeps the
  record-type rule. Selenium also had no radius of its own and took carbon's
  1.70 A; it now takes Bondi's 1.90 A. The FreeSASA golden files are
  regenerated on the same atom set, the 51 entries holding modified residues
  change, and 8A71's modified DNA is scored for the first time: SASA mean
  absolute error against FreeSASA moves from 4.93 to 4.77 A^2 over 203 entries.
- **The validation harness charged pixelfold for ligands PLIP never profiled.**
  It decided what a ligand was by the HETATM record, while PLIP leaves out
  modified residues of a chain and sets aside some crystallisation additives,
  so pixelfold's contacts to those counted as false positives. PLIP golden
  files now list the residues PLIP profiled, and only edges touching one of
  them are scored. Hydrogen-bond precision moves from 0.58 to 0.64 and
  water-bridge precision from 0.50 to 0.53; no recall changes. A golden file
  that fails to parse, or a ligand-scoped one whose interactions name a ligand
  it does not list, now stops the harness instead of skipping its entry, which
  had left the gates it feeds silently unchecked.

## [2.0.0] (2026-10-22)

The structure-IO milestone. PixelFold now reads structures through
`mmcif-model` and `mmcif-cif` instead of pdbtbx, and comparing the two readers
across the benchmark corpus turned up a string of numbering and conformer bugs,
fixed here. Multi-model files and biological assemblies become first-class:
pick a model or report on all of them, and build an assembly with each copy of
a chain named. Loading a large structure is about 30 times faster.

This is a major release because the library API changes. `Atom::residue_seq`
is signed and `Atom::chain_id` is an eight-byte `ChainId`; loading goes through
`parser::Load`; and `pixelfold-tui` and `pixelfold-core` now expose only what
the other crates use (see Removed). The command line keeps every option it had;
single-model output keeps its shape, with residues listed in label-chain order.

### Added

- **Biological assemblies**, built on request with `--assembly` (the primary
  assembly) or `--assembly=ID`, on every subcommand and in the viewer, where `u`
  switches between the deposited coordinates and the primary assembly. Copies
  of a chain are named `A`, `A-2`, `A-3`, with a ligand or water following its
  protein copy, so `select chain A-2` picks one copy and the residue network
  shows the interfaces between copies (6ZPA's hexamer, 7AEW's two copies of
  its asymmetric unit). Where the deposited coordinates hold more than one
  biological unit, the assembly keeps only its own chains. An unknown id is an
  error naming those the file declares; a legacy PDB file builds the assemblies
  its `REMARK 350` declares; and `rin --mvsj` refuses `--assembly`, since Mol*
  names copies by operator. Library users name what to build with
  `parser::Load` (model, assembly, altloc policy), which
  `ParsedStructure::protein`, `App::load_protein` and `pixelfold_tui::view`
  now take.
- **`--model N|all`**, a global option choosing which model of a multi-model
  file (an NMR ensemble, say) every subcommand and the viewer load; the first
  model remains the default, and an absent one is an error naming the models
  the file has. `select model N` now has something to match. `--model all`
  makes `ss`, `sasa` and `interactions` report every model, with a leading
  `model` column in the table and TSV and a `model` field in the JSON; output
  for a single model keeps its shape. A selection is then an error only when it
  matches nothing in any model. `rin`, `render` and the viewer work on one model
  at a time and refuse `all`.
- **Model stepping in the viewer**: `m` and `M` show the next and previous
  model, wrapping at either end, without re-reading the file. The camera stays
  put, the info bar names the model, and the selection and network pane start
  over for the new model's atoms. Library users read a file once with
  `parser::ParsedStructure` and build any model from it; `App::load_protein`
  and `pixelfold_tui::view` take the model to open.

### Changed

- **Chain ids hold eight bytes**, as `ChainId`, so a four-character author id
  keeps room for a copy number (`AAA-2`); atom names stay four.
- **Loading a large structure is about 30 times faster.** Secondary-structure
  assignment tried every pair of residues for a beta bridge; it now tries only
  the pairs within a step of some hydrogen bond's two ends, the only pairs a
  bridge pattern can match. On 1AON (58,870 atoms) the pass fell from 1.7 s to
  11 ms and the whole load from about 1.8 s to 60 ms, with every assignment
  across the benchmark corpus unchanged.
- **Element symbols are matched without allocating.** Radii, colours, metal
  tests and ligand atom typing upper-cased each symbol into a new string, on
  every atom and every frame; building a frame's scene for 1AON fell from
  2.4 ms to 0.6 ms.

- **Structures load through `mmcif-model` and `mmcif-cif`** in place of pdbtbx.
  One parse now serves the coordinates, the component definitions and the
  assembly records, so the hand-written CIF reader and the space-group retry are
  gone. Reports list residues in label-chain order (polymer chains, then
  ligands, then waters, as an mmCIF file writes them) rather than grouped by
  author chain, and an atom group lists its atoms in file order. Gzipped files
  and RCSB's numbered assembly files (`.pdb1`) still load, BinaryCIF (`.bcif`)
  loads too, and a file whose assembly or entity records cannot be read loads
  its coordinates with a warning instead of failing; a file whose extension
  names no structure format is refused. A missing B-factor is missing, shown as
  `n/a` and coloured grey, rather than read as 0.

### Removed

- **Most of `pixelfold-tui`'s public API.** It exposes `view`, `search` and
  `ViewProtocol`; `App`, its fields and methods, `RenderedFrame`, the picking
  helpers, and the `search` module are internal, and `DisplayMode` moves here
  from `pixelfold-core`.
- **`pixelfold-core` internals from the public API**: the `assembly`, `dssp`
  and `spatial` modules (`HBond` stays, re-exported), the interaction
  submodules (`aromatic`, `classify`, `connectivity`, `hydrogens`, `params`,
  `topology`) and `Charge`, the selection syntax tree (`Selection`, `Cmp`,
  `ResiRange`, `AtomSet::of`), `filter_altlocs`, and
  `Dictionary::from_block`, which took an `mmcif-cif` type. No `mmcif-*` type
  is part of the public API, so a future mmcif upgrade is not a breaking change.
- **Unused code**: `rin::HBondGraph` with its motifs and analysis,
  `load_protein_backbone`, `load_protein_ca_only`,
  `SecondaryStructure::color_256`, `params::GRID_CELL_SIZE`, and
  `Grid::within`.
- `App::residue_colors` in `pixelfold-tui`. It mapped a residue number to a
  secondary structure with no chain, and every atom already carries its own
  residue's secondary structure.

### Fixed

- **Backbone connections took another chain's colour.** In backbone mode a
  connection was coloured by the first residue anywhere in the structure with
  the same number, so in a multi-chain structure chain B's residue 10 showed
  chain A's residue 10 secondary structure. Each connection now takes its own
  residue's.
- **Modified residues were selected as ligands.** `protein` matched only the 20
  standard amino acid names and `ligand` any non-water, non-metal HETATM, so a
  selenomethionine or phosphotyrosine inside a chain was a ligand and not
  protein (3G5S's 14 MSE, 9I24's PTR, 8OW6's HYP). `protein`, `nucleic`,
  `ligand` and `water` now follow the file's entities, recorded on each atom as
  `Atom::molecule`; a free amino acid declared as a non-polymer stays a ligand.
  A file without entity records keeps the old name and record-type rules.
- **Hydrophobic contacts depended on atom order.** The reduction that keeps
  the closest contact per atom and residue ran from whichever side came first
  in the file, and the two directions can keep different numbers of contacts,
  so a residue network's summed energies followed atom order. It now runs from
  the ligand side, as PLIP's does, then from the earlier residue: in 1IEP,
  imatinib's contacts with Lys 271 count 2 instead of 1. Hydrophobic contacts,
  pi-stacking and disulfides list their sides in that order too. No residue
  pair gains or loses an interaction across the benchmark corpus, and the
  validation report is unchanged.
- **Residues differing only by insertion code were merged** in two places.
  Selecting an atom highlighted every residue with its chain and number, so
  clicking thrombin's `1A` lit up `1C`, `1B`, `1A` and `1` together; and
  `HBondGraph` made one node of them all, pooling their atoms and hydrogen
  bonds. Both now key a residue by chain, number and insertion code;
  `ResidueNode` carries `insertion_code` and `HBondGraph::get_node` takes it.
- **The default altloc policy kept the minor conformer** wherever a side chain
  was split and the backbone was not. pdbtbx copied the shared atoms into every
  conformer at a fraction of their occupancy, so the backbone vote always tied
  and the first label won: 76 residues across 18 of the 206 benchmark entries
  kept the less occupied side chain (1GCI SER 3 kept A at 0.47 over B at 0.53).
  Shared atoms now stay unlabelled at their own occupancy, `--altloc b` no
  longer drops the backbone of a residue modelled only in A, and `--altloc all`
  no longer reports shared atoms twice.
- **Negative residue numbers wrapped** to 4294967295 and up, so those residues
  printed wrong, could not be selected, and never matched a reference in
  `pixelfold-validate`. Residue numbers are signed, and `resi` takes negative
  numbers and ranges (`resi -5--1`).
- **Unknown atoms (UNX) were typed as uranium**, the element inferred from the
  name `UNK` because the `X` element could not be read. They keep element `X`.
- **A numeric-looking component id was rewritten** (`4E2` read as `400`), so the
  ligand printed under the wrong name. Component ids are kept verbatim.
- **A legacy PDB file's insertion-code residues were re-sorted**, so a chain
  numbered `1C, 1B, 1A, 1` (thrombin's light chain) was read as `1, 1A, 1B, 1C`,
  and secondary structure and the bond graph saw the wrong neighbours. The PDB
  and mmCIF files of one entry now give the same result (3PTB, 5JZY, 5LCE,
  6T3Q, 6T4A did not).
- **Two kinds of partial biological assembly went unflagged**: an asymmetric
  unit copied whole across separate generation records (7AEW) and a chain group
  placed by a crystal operator (3BOG). An assembly is now partial when it places
  a chain more than once or places its chains by different operators.

## [1.0.0] (2026-08-02)

The rendering and network milestone. The viewer stops drawing an approximation of
the structure and draws real pixels, through the same framebuffer path the new
headless `render` subcommand writes, and the residue interaction network becomes
a first-class linked view with its own analytics, energies, and filters. The
interaction geometry is untouched by this release, so the agreement figures
published for 0.6.0 stand as measured.

### Added

- **Z-buffered sphere-impostor rasteriser** (`render/src/raster.rs`). Each atom is
  a screen-space disc whose front-hemisphere depth is solved per pixel and depth
  tested, so there is no sort and no painter's-algorithm bleeding, and the surface
  normal falls out of the same solve for Lambert shading. It writes an RGBA8
  framebuffer carrying parallel depth and atom-id buffers, with depth-cue fog over
  distance. Everything renders through that one buffer: the interactive viewer,
  PNG export, and the tests.
- **Terminal graphics sinks** (`render/src/sink/`): the Kitty graphics protocol (a
  base64 PNG payload chunked at 4096 bytes), iTerm2's OSC 1337 inline image, and a
  truecolor upper-half-block fallback packing two vertical pixels per cell.
  `detect.rs` reads terminal capability and SSH from the environment and
  `to_terminal` dispatches on what it finds. The payload is a compressed PNG, so a
  remote session does not pay for raw pixels. Sixel is reachable in the
  interactive viewer through `ratatui-image` but not yet through `to_terminal`.
- **`render` subcommand**, space-filling spheres with no terminal UI involved.
  With `-o` it writes a PNG (1200x900 unless told otherwise); without one it emits
  a terminal graphics escape when standard output is a terminal, sized to the cell
  grid so half-block output cannot overflow it, or pipes a PNG when redirected. It
  takes `--width`/`--height`, `--color element|bfactor|ss`, `--backbone`,
  `--select`, `--radius-scale`, `--no-depth-cue`, `--slab FAR NEAR`, and
  `--protocol auto|kitty|iterm2|half-block`, and it errors when the selection
  leaves nothing to draw rather than writing an empty image.
- **Slab (clipping planes)**, a depth band that slices a dense structure open:
  `--slab FAR NEAR` headless, and in the viewer `k` to toggle, `[` and `]` to move
  the band through depth, `,` and `.` to narrow or widen it.
- **ID-buffer mouse picking.** The rasteriser records the source atom index at
  every pixel it wins, so a click reads the atom under it directly instead of
  hunting for the nearest projected centre. The last frame's buffer and image area
  are kept for exactly this.
- **Depth-aware overlay lines** (`render/src/overlay.rs`). Backbone connections and
  hydrogen bonds interpolate their endpoint depths against the rasteriser's own
  z-buffer, so a bond that passes behind an atom is hidden by it.
- **Golden-image regression test** over a fixed synthetic scene at a fixed camera,
  which needs no structure file and so runs in CI. It guards projection,
  occlusion, shading, depth cue, and the PNG encode together, with a per-channel
  tolerance and a changed-pixel budget so a rounding difference at a sphere rim
  does not fail the run.
- **Linked network pane** (`g`), the residue interaction network drawn beside the
  3D view. Selection is shared both ways, so clicking a node highlights that
  residue in the structure and picking an atom locates its node. Two layouts (`l`):
  a deterministic Fruchterman-Reingold force-directed layout
  (`core/src/rin/layout.rs`), which reads topology independent of orientation, and
  a spatial layout that projects every residue through the 3D camera, so the graph
  turns with the structure. Edges are coloured by interaction kind with a legend on
  the border; nodes are haloed by betweenness centrality so hubs stand out.
- **Network analytics** (`core/src/rin/analysis.rs`): connected components,
  betweenness centrality by Brandes, articulation points, and a breadth-first
  shortest path between two residues. `pixelfold rin --analyze` prints them as a
  structural report (components, hub residues by betweenness, cut residues)
  instead of exporting the graph, and `--path FROM TO` prints the shortest chain of
  interactions between two residues, or says they lie in different components.
- **Per-edge interaction energies** (`core/src/interactions/energy.rs`), one
  literature-sourced magnitude in kJ/mol per interaction kind with a primary
  citation for each, on the intrinsic single-contact basis RING's own weights use
  rather than gas-phase extremes or net folding free energy. Every network edge
  sums the contacts folded into it (weight times count) and carries that sum into
  node-link JSON, GraphML, and the TSV edge list, and the pane brightens edges by
  it on a log scale. The four kinds RING scores land inside these independently
  sourced ranges. Two values carry caveats stated at the table: the hydrophobic
  weight is per atom-atom contact, the scale detection works at, so an edge's many
  contacts sum to the per-residue-pair magnitude; and the disulfide is the covalent
  S-S bond energy, an order of magnitude above every non-covalent contact, so it is
  compressed at the display layer rather than distorted in the table.
- **Relative solvent accessibility.** `relative_sasa` divides per-residue area by
  Tien 2013's theoretical maximum. `pixelfold sasa` reports it as an `rsa` column
  beside the absolute area, and in the network pane `e` toggles node colouring
  between secondary structure and burial (per-node RSA, summed once on open from
  per-atom SASA), so the structural core and the solvent-facing surface separate at
  a glance.
- **Network filters**, both reflected in the pane title with a live count of
  visible edges: `t` steps the interaction-kind filter (all kinds, then each kind
  present in turn) and `y` the chain filter (all, intra-chain, or the inter-chain
  interface only). The energy intensity scale rescales to whatever remains visible,
  so a filtered view still spans its own range.
- **MolViewSpec export**: `pixelfold rin --mvsj` writes a `.mvsj` scene for Mol\*,
  the structure with the network's residues highlighted. A 4-character id resolves
  to its canonical RCSB mmCIF URL, which any Mol\* can fetch; a local path becomes
  a `file://` URL for a desktop or locally served viewer.
- The `rin` TSV edge list prints aligned fixed-width columns when it lands on a
  terminal and genuine tab-separated values when written to a file or a pipe, so
  reading it by eye and piping it into a script no longer trade against each other.
- `rust-toolchain.toml` and `rustfmt.toml` pin the toolchain and the 2024 style
  edition, so local formatting and CI cannot diverge by version.

### Changed

- The viewer draws the shared RGBA framebuffer, shown through `ratatui-image`
  (iTerm2, Kitty, or Sixel, with a half-block fallback), instead of the braille
  canvas it rendered to before. `pixelfold-tui` is split accordingly into `app`
  (event loop), `inputs` (keyboard and mouse), `render` (framebuffer build), and
  `ui` (per-frame layout), and the viewer's framebuffer is the same one `render`
  writes to a PNG.
- The network is built on a worker thread the first time the pane opens.
  Detection, layout, and the accessibility pass no longer freeze the event loop on
  a large structure; the pane shows a placeholder until the result arrives. The
  loop polls for it inside the 16 ms wait it already performs, so frame timing is
  unchanged.
- The interactive iTerm2 draw encodes each frame at the fastest compression
  (`to_png_bytes_fast`), roughly halving the per-frame cost. For a local terminal a
  larger payload is cheaper than the time spent shrinking it. Files written by
  `render` keep the default compression, where the trade runs the other way.
- The event loop drains every buffered event before drawing again, so a burst of
  key repeats advances the view once rather than once per event.
- The framebuffer is capped in its largest dimension and scaled to fill the image
  area, rather than always allocated at the full cell grid times font size.

### Removed

- Solvent-surface point generation: `SurfacePoint`,
  `SurfaceCalculator::calculate_surface`, and the Fibonacci-sphere distribution
  behind it. Nothing had drawn the points since the braille viewer, yet every load
  paid to build them. Per-atom Shrake-Rupley areas (`calculate_atom_sasa`) are
  untouched: they remain what `pixelfold sasa` reports, what the validation harness
  compares against FreeSASA, and what the network pane sums per residue for burial
  colouring.
- The `--no-surface` flag, on `pixelfold view` and on the bare-structure viewer
  path. It skipped a calculation that no longer happens.

### Fixed

- Mouse picking maps a click proportionally onto the framebuffer instead of
  through a fixed cell-to-pixel size. The framebuffer can be a lower resolution
  than its image area and the terminal scales it to fit, so the fixed mapping
  selected the wrong atom, or none, whenever the two disagreed.
- The headless CLI tests took their structures from the gitignored benchmark cache
  and skipped themselves when it was absent, so on a clean checkout, CI included,
  all 24 reported as passing without running the binary once. The fixtures
  (crambin, ubiquitin, thrombin) are now committed under `tests/fixtures` and a
  missing one is a hard failure.

## [0.6.0] (2026-07-24)

The interaction-engine milestone. PixelFold now computes non-covalent
interactions and a residue interaction network locally, and publishes their
agreement with reference implementations.

### Added

- **Non-covalent interaction engine** over the uniform spatial grid: hydrogen
  bonds, salt bridges, hydrophobic contacts, pi-stacking, pi-cation, halogen
  bonds, water bridges, metal coordination, and disulfides. Detection uses
  angular terms, not distance-only cutoffs, with thresholds derived from PLIP's
  source.
- **Ligand support**: the PDB Chemical Component Dictionary carried in a file is
  parsed for atom elements and covalent connectivity, and a bond graph derives
  atom typing for arbitrary components. Aromatic ring perception and polar
  hydrogen inference at pH 7.4.
- **Selection language** (a `nom` parser to an atom bitset): `chain`, `resi`,
  `resn`, `within X of <sel>`, `byres`, `protein`/`ligand`/`water`/`metal`, and
  boolean combinators, available on every subcommand through `--select`.
- **Headless CLI subcommands**: `interactions`, `rin`, `ss`, `sasa`, and
  `fetch`, alongside the `view` TUI. `interactions` writes table, TSV, or JSON;
  `rin` writes node-link JSON, GraphML (for Cytoscape and Gephi), or a TSV edge
  list.
- **Validation harness** (`pixelfold-validate`): per-entry comparison against
  DSSP 4, FreeSASA, and PLIP across a stratified 205-entry benchmark, a CI
  regression gate (`benchmarks/thresholds.toml`), and the README validation
  table. Reference goldens are produced by pinned Docker images under
  `benchmarks/tools/`.
- Cargo workspace split into `pixelfold-{core,render,fetch,tui,cli,validate}`;
  `core` carries no terminal or rendering dependencies.

### Validation

- Secondary structure (DSSP 4): Q3 97.1% over 194 entries; backbone H-bond edge
  F1 0.97 over 205.
- Solvent-accessible surface (FreeSASA): mean absolute error 4.93 A^2, median
  3.77, Pearson r 0.99 over 202 entries.
- Interactions (PLIP, 60 protein-ligand complexes, per-type F1): halogen-bond
  0.98, pi-stacking 0.92, hydrophobic 0.87, pi-cation 0.84, salt-bridge 0.83,
  metal-coordination 0.81, hydrogen-bond 0.65, water-bridge 0.59.

### Fixed

- Alternate locations resolve to one conformer per residue, chosen by backbone
  occupancy so a disordered sidechain never shifts the main chain. A residue is
  no longer left chimeric, and an atom the chosen conformer omits keeps its own
  copy rather than being dropped.
- Under `--altloc all`, interaction edges duplicated across conformers are
  collapsed to a single edge.
- Backbone hydrogen bonds are capped at the two most favourable per donor, as
  DSSP 4 records them.
- The PLIP golden generator takes each interaction's ligand from the
  interaction's own fields rather than the binding-site key, correcting metal
  residue attribution in the reference set.

### Known limitations

- RING 4.0, the whole-structure protein-protein reference, is a license-gated
  binary and is not generated here; the harness reads a golden of the same shape
  once one is available.
- Halogen-bond recall deliberately excludes organic fluorine (Auffinger et al.,
  2004); this is a documented divergence from PLIP, which counts it.
