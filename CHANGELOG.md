# Changelog

All notable changes to PixelFold are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
