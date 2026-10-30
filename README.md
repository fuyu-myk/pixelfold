# PixelFold

[![CI](https://github.com/fuyu-myk/pixelfold/actions/workflows/ci.yml/badge.svg)](https://github.com/fuyu-myk/pixelfold/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A protein structure viewer and residue interaction network engine that runs
entirely in your terminal.

![PixelFold demo](https://raw.githubusercontent.com/fuyu-myk/pixelfold/main/docs/demo.gif)

## Validation

The viewer is the demo. The interaction engine is the product, and it is measured
against existing trusted tools.

<!-- VALIDATION:START -->
Validated 205 benchmark entries against DSSP 4, FreeSASA, PLIP and RING.

| Metric | Reference | Entries | Result |
| --- | --- | --- | --- |
| Secondary structure (Q3) | DSSP 4 | 194 | 97.1% |
| Backbone H-bond edges (F1) | DSSP 4 | 205 | 0.97 |
| SASA mean abs error | FreeSASA | 203 | 4.77 A^2 |
| SASA median abs error | FreeSASA | 203 | 3.61 A^2 |
| SASA correlation (Pearson r) | FreeSASA | 203 | 0.99 |

Per-type interaction agreement (precision / recall / F1):

| Interaction type | Reference | Entries | Pred edges | Ref edges | Precision | Recall | F1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| hydrogen-bond | PLIP | 60 | 676 | 519 | 0.64 | 0.82 | 0.68 |
| salt-bridge | PLIP | 60 | 129 | 127 | 0.91 | 0.93 | 0.87 |
| hydrophobic | PLIP | 60 | 435 | 377 | 0.87 | 0.96 | 0.90 |
| pi-stacking | PLIP | 60 | 23 | 22 | 0.97 | 0.95 | 0.92 |
| pi-cation | PLIP | 60 | 6 | 17 | 0.97 | 0.88 | 0.84 |
| halogen-bond | PLIP | 60 | 4 | 5 | 1.00 | 0.98 | 0.98 |
| water-bridge | PLIP | 60 | 511 | 345 | 0.59 | 0.83 | 0.65 |
| metal-coordination | PLIP | 60 | 444 | 338 | 0.82 | 0.92 | 0.81 |
| hydrogen-bond | RING | 205 | 49374 | 43512 | 0.79 | 0.90 | 0.83 |
| salt-bridge | RING | 205 | 4251 | 1545 | 0.38 | 1.00 | 0.50 |
| pi-stacking | RING | 205 | 517 | 3756 | 1.00 | 0.25 | 0.32 |
| pi-cation | RING | 205 | 453 | 454 | 0.72 | 0.71 | 0.62 |
| halogen-bond | RING | 205 | 7 | 9 | 0.99 | 0.98 | 0.97 |
| metal-coordination | RING | 205 | 567 | 558 | 0.89 | 0.91 | 0.86 |
| disulfide | RING | 205 | 249 | 249 | 1.00 | 1.00 | 1.00 |
<!-- VALIDATION:END -->

These are agreement figures and never improvement claims. Where PixelFold and the
reference disagree, the reference is right by construction. The harness lives in
this repository and runs in CI, allowing the numbers above to be reproduced, rather than simply
asserted. See [Validation and methods](#validation-and-methods) for what the
weaker rows mean.

## Install

```bash
cargo install pixelfold
```

This needs Rust 1.88 or newer. Each [GitHub release](https://github.com/fuyu-myk/pixelfold/releases)
also carries prebuilt binaries for macOS, Linux and Windows, with shell and
PowerShell installers. From v2.1.0 the Linux binaries need glibc 2.28 or newer,
which covers RHEL/Rocky 8 and later, the usual cluster login node; earlier ones
need glibc 2.34.

## A 30 second demo

**Dihydrofolate reductase with its inhibitor methotrexate**. How does this drug
interact with the enzyme, and does it matter structurally?

```console
$ pixelfold interactions 4DFR --ligand MTX
kind           chain_a  resi_a  resn_a  atoms_a     chain_b  resi_b  resn_b  atoms_b  distance  angle
-------------  -------  ------  ------  ----------  -------  ------  ------  -------  --------  -----
hydrogen-bond  A            52  ARG     NH2         A           161  MTX     O            2.88  152.8
hydrogen-bond  A            52  ARG     NH2         A           161  MTX     OE2          3.43  138.2
hydrogen-bond  A            57  ARG     NH1         A           161  MTX     O1           2.57  155.1
hydrogen-bond  A            57  ARG     NH2         A           161  MTX     O2           2.72  170.8
... 53 rows
```

Next, analyze the network to find where the important residues are located.

```console
$ pixelfold rin 4DFR --analyze
Residue interaction network: 311 residues, 629 edges, 5 components
...
Top hubs by betweenness centrality:
  residue    resn  ss  degree      between
  A/23       ASN   C        5     13880.74
  A/24       LEU   C        5     12961.77
  B/162      MTX   C       17     10909.25
  A/161      MTX   C       15     10552.69
  B/22       TRP   C        7     10378.71
```

Both copies of the inhibitor rank among the most central residues in the
structure. Both are also articulation points, where the network will split if one is removed.
This local computation takes about a second, and provides insight to the fold.

Everything here is scriptable. `--format json|tsv|graphml` feeds a pipeline,
`--select` narrows any command with the same query language, and `pixelfold 4DFR`
opens the interactive viewer over the same data.

## Why it exists

Structures usually live where the compute is, like a cluster login node, a remote
workstation, or a container. To look at one normally, one must copy it to their machine,
open a desktop application, before any analysis can be made.

PixelFold runs in the shell that already has the file, with no need for X11, a display server,
a GPU, browser upload, and Python environment. It is a
single static binary, which drops into a Nextflow or Snakemake step as easily as
`grep`. The analysis is in text, so it diffs, greps, and pipes, while being
deterministic.

## Honest comparison

PixelFold is is in its infancy. It is not competitive on rendering, breadth, or
maturity with any established viewer, and it does not claim to be the reference on interaction
chemistry. Below is a compilation of how different tools compare.

| Tool | Where it beats PixelFold | Where PixelFold differs |
| --- | --- | --- |
| [PyMOL](https://pymol.org) | Ray-traced publication figures, cartoons, surfaces, electrostatics, trajectories, twenty years of plugins, an embedded Python API, and near-universal adoption | Terminal-native with no display server; MIT throughout; the network is the primary output rather than a side feature |
| [Mol\*](https://molstar.org) | GPU rendering at a quality and scale PixelFold cannot approach, million-atom and integrative models, MD streaming, zero install, and it powers RCSB, PDBe and AlphaFold DB | No browser and no upload, so embargoed coordinates stays local. PixelFold complements through its exportable `.mvsj` scenes |
| [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) | Cryo-EM density, structure building and refinement, VR, AlphaFold workflows, a bundle ecosystem, and NIH-funded maintenance | MIT versus a non-commercial licence, so it can be vendored freely; a small static binary rather than a GB-class install |
| [RING 4.0](https://ring.biocomputingup.it) | The peer-reviewed network reference the field accepts. Ligand coverage across the whole PDB chemical component dictionary, plus pi-hydrogen bonds, which PixelFold does not detect. RING-PyMOL already does ensemble and MD contact analysis | Readable MIT source you can fork or embed as a Rust library, offline with no web upload, and graph analytics in the same binary. Distribution differences, not accuracy claims |
| [PLIP](https://plip-tool.biotec.tu-dresden.de) | The yardstick for interaction detection. Its protonation-aware pipeline is more physically grounded, and it covers DNA, RNA and protein-protein interfaces | No Python or OpenBabel to install; profiles the whole structure as a graph rather than only ligand sites |
| [asciiMol](https://github.com/dewberryants/asciiMol) | Runs on genuinely any terminal, including plain ncurses, where PixelFold degrades to half-blocks and looks worse. `pip install` into an environment you already have. Native `.xyz` and SMILES | Protein-scale mmCIF and PDB with chains, conformers and a selection language, and the analysis engine behind it |

## Validation and methods

Reference implementations: DSSP 4 for secondary structure, FreeSASA for
solvent-accessible surface area, PLIP for the interactions of bound ligands, and
RING 4.0 for interactions across the whole structure, over a stratified benchmark
of 205 PDB entries, 60 of them protein-ligand complexes. Run it yourself with
`cargo run -p pixelfold-validate -- --check`.

Each interaction reference is scored only over the part of the structure it looks
at. Against PLIP, only PixelFold's edges that touch a ligand PLIP profiled count.
Against RING, edges to water are left out, as are pairs of amino acids fewer than
three residues apart in one chain, which RING never pairs at its defaults.

Interaction geometry follows PLIP's published thresholds, with the exact source
variable recorded next to every constant in
[`params.rs`](https://github.com/fuyu-myk/pixelfold/blob/main/crates/pixelfold-core/src/interactions/params.rs). Each known
difference from RING is also noted. Secondary structure is
Kabsch-Sander electrostatic hydrogen bonding with inferred amide hydrogens.

**The weak rows are a modelling choice, and they are the honest headline.**
Hydrogen bonds (F1 0.68, precision 0.64) and water bridges (F1 0.65, precision
0.59) are the two classes where PixelFold reports what the geometry allows while
PLIP commits to one reading of it. Of the 261 hydrogen-bonded residue pairs
PixelFold finds and PLIP does not, 94 come from a rotatable hydroxyl, which
PixelFold turns toward each acceptor and PLIP fixes wherever OpenBabel placed
its hydrogen; 77 also form a salt bridge, which PLIP then drops and RING keeps;
and the other 90 have nitrogen donors, many of them histidine tautomers, amide
flips and ligand nitrogens the coordinates leave open. A water bridges up to two
donors and two acceptors, one per lone pair and per hydrogen, where PLIP keeps
one donor per acceptor and drops a donor that already hydrogen-bonds the ligand.
Those two classes are roughly half of all protein-ligand interactions, so for
ligand-site work where a curated answer matters, use PLIP.

**Against RING, the weak rows come from following PLIP's definitions.** RING 4.0
draws three types differently, and where the two disagree RING is the reference.
Salt bridges (F1 0.50, precision 0.38) are a superset: PixelFold finds 1541 of
RING's 1545 and 2710 more, because PLIP pairs charge centres up to 5.5 A apart
where RING asks for about 4.0 A between atoms, and RING pairs only amino acids.
Pi-stacking (F1 0.32, recall 0.25) runs the other way: 515 of PixelFold's 517
stacks are in RING, but RING's wider window (about 6.5 A between ring centres,
against PLIP's 5.5 A with at most 2.0 A of lateral offset) finds seven times as
many. Pi-cation (F1 0.62) disagrees in both directions: PixelFold, like PLIP,
counts histidine as a cation and needs the cation within 2.0 A of the ring's
axis, while RING's cations here are only arginine and lysine, within about
5.0 A. Hydrogen bonds agree at F1 0.83 (precision 0.79): RING's shorter 3.9 A
cutoff and its added hydrogen-to-acceptor limit leave it fewer bonds. The
halogen-bond rows rest on five and nine reference edges, and an entry where
neither side finds one counts as agreement, so they say little.

Per-edge interaction energies are a literature-sourced scoring table with one
primary citation per value, in
[`energy.rs`](https://github.com/fuyu-myk/pixelfold/blob/main/crates/pixelfold-core/src/interactions/energy.rs). They are
representative magnitudes, not computed free energies.

## Architecture

A Cargo workspace. `pixelfold-core` has no terminal or rendering dependencies, so
the analysis is usable as a library.

- [`pixelfold-core`](https://github.com/fuyu-myk/pixelfold/blob/main/crates/pixelfold-core/src): parsing, DSSP, SASA, the
  interaction engine, the selection language, and the residue interaction network
  with its graph analytics
- [`pixelfold-render`](https://github.com/fuyu-myk/pixelfold/blob/main/crates/pixelfold-render/src): a z-buffered
  sphere-impostor rasteriser writing one RGBA framebuffer, and the sinks that put
  it on a terminal (Kitty, iTerm2, truecolor half-blocks) or in a PNG
- [`pixelfold-tui`](https://github.com/fuyu-myk/pixelfold/blob/main/crates/pixelfold-tui/src): the `ratatui` front end and the
  linked network pane
- [`pixelfold-cli`](https://github.com/fuyu-myk/pixelfold/blob/main/crates/pixelfold-cli/src): the `pixelfold` binary
- [`pixelfold-validate`](https://github.com/fuyu-myk/pixelfold/blob/main/crates/pixelfold-validate/src): the agreement harness

The CLI computes everything the viewer shows headlessly.

## Limitations

Being explicit about what this does not do:

- **The RING reference comes from its web server.** RING's standalone binary is
  licensed to one person and cannot be passed on, so the RING golden files were
  generated by the public RING 4.0 web server, which reports no version. The
  committed files keep the comparison reproducible offline, but regenerating
  them depends on that server, and its anonymous rate limit makes a full rerun
  take about a day and a half.
- **Hydrogen bonds and water bridges over-report** relative to PLIP, as above.
- **Salt bridges over-report and pi-stacking under-reports** relative to RING,
  as above, because both follow PLIP's geometry.
- **No cartoons or ribbons, no molecular surfaces, no ray tracing.** Space-filling
  spheres and bonds only.
- **No trajectories.** A multi-model file such as an NMR ensemble works model
  by model (`--model N`, or `--model all` for a report on every model), but a
  molecular-dynamics trajectory is not read.
- **No cryo-EM density, no electrostatics, no structure building or refinement.**
- **No pi-hydrogen bonds**, which RING 4.0 detects.
- **Rendering quality depends on your terminal.** Kitty and iTerm2 show real
  pixels; anywhere else falls back to half-blocks, which is coarse. Pass
  `--protocol half-block` if the structure pane comes up blank.
- **Not peer reviewed.** The validation harness is the evidence on offer.

## Usage

```bash
# Interactive viewer, from a file or a 4-character PDB id (fetched and cached)
pixelfold structure.cif
pixelfold 4HHB

# Non-covalent interactions
pixelfold interactions 4DFR --ligand MTX
pixelfold interactions 4DFR --type hydrogen-bond --format tsv

# The residue interaction network, for Cytoscape, Gephi, or a script
pixelfold rin 4DFR --format graphml -o network.graphml
pixelfold rin 4DFR --analyze
pixelfold rin 4DFR --path A/23 A/94
pixelfold rin 4DFR --mvsj -o scene.mvsj

# Secondary structure and solvent accessibility
pixelfold ss 1UBQ --format tsv
pixelfold sasa 1UBQ --format json

# Render to a PNG, or to the terminal
pixelfold render 4HHB -o hemoglobin.png --width 1600 --height 1200

# One model of an NMR ensemble, or every model with a model column
pixelfold ss 1D3Z --model 3
pixelfold sasa 1D3Z --model all --format tsv

# The biological assembly instead of the deposited coordinates; copies of a
# chain are named A, A-2, A-3
pixelfold rin 6ZPA --assembly --analyze
pixelfold interactions 7AEW --assembly --select 'chain AAA-2'

# Any command, narrowed by the selection language
pixelfold interactions 1STP --select 'byres within 5 of resn BTN'
```

### Viewer controls

| Key | Action |
| --- | --- |
| `WASD` | Rotate; `z`/`x` roll |
| `+`/`-`, arrows | Zoom, pan |
| `1`/`2` | All atoms, C-alpha backbone |
| `c` `b` `h` | Bonds, B-factor colouring, hydrogen bonds |
| `i`, click | Inspect mode, pick an atom; `r` highlights its residue |
| `g` | Linked interaction network pane |
| `l` `e` | Network layout (force or spatial), node colour (structure or burial) |
| `t` `y` | Filter network edges by interaction type, by chain |
| `k` `[` `]` `,` `.` | Slab: toggle, move, and resize a depth slice |
| `m`/`M` | Next, previous model of a multi-model file |
| `u` | Biological assembly on or off |
| `f`, `q` | Reset the view, quit |

## Citations

```bibtex
@article{Kabsch1983,
  author = {Kabsch, Wolfgang and Sander, Christian},
  title = {Dictionary of protein secondary structure: Pattern recognition of hydrogen-bonded and geometrical features},
  journal = {Biopolymers},
  volume = {22},
  number = {12},
  pages = {2577--2637},
  year = {1983},
  doi = {10.1002/bip.360221211}
}

@article{Kyte1982,
  author = {Kyte, Jack and Doolittle, Russell F.},
  title = {A simple method for displaying the hydropathic character of a protein},
  journal = {Journal of Molecular Biology},
  volume = {157},
  number = {1},
  pages = {105--132},
  year = {1982},
  doi = {10.1016/0022-2836(82)90515-0}
}

@article{Bondi1964,
  author = {Bondi, A.},
  title = {van der Waals Volumes and Radii},
  journal = {The Journal of Physical Chemistry},
  volume = {68},
  number = {3},
  pages = {441--451},
  year = {1964},
  doi = {10.1021/j100785a001}
}

@article{Shrake1973,
  author = {Shrake, A. and Rupley, J. A.},
  title = {Environment and exposure to solvent of protein atoms. Lysozyme and insulin},
  journal = {Journal of Molecular Biology},
  volume = {79},
  number = {2},
  pages = {351--371},
  year = {1973},
  doi = {10.1016/0022-2836(73)90011-9}
}
```
