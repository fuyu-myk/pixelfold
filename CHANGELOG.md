# Changelog

All notable changes to PixelFold are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
