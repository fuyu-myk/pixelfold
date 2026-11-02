---
name: pixelfold
description: Analyse protein structures and ligand poses with the pixelfold command line. Use it to rank co-folding or docking poses of a ligand against a crystal structure (contacts kept, pocket-aligned RMSD, clashes, inverted stereocentres, E/Z flips, PoseBusters-style geometry), to check poses PoseBusters-style when there is no crystal structure (clashes, geometry, stereochemistry against the Chemical Component Dictionary, AlphaFold 3 or Boltz confidence), to list the non-covalent interactions around a ligand or between chains, to build residue interaction networks of one structure, an NMR ensemble or a prediction's samples (hubs, cut residues, shortest routes), to report secondary structure or solvent accessibility per residue, or to render a structure to a PNG.
---

# pixelfold

`pixelfold` reads mmCIF, PDB and BinaryCIF files (gzipped or not) or fetches a
4-character PDB id from RCSB and caches it. Check it is installed with
`pixelfold --version`; if not, `cargo install pixelfold`, or take a binary from
https://github.com/fuyu-myk/pixelfold/releases.

## Rules

- Pass `--format json` to every report and parse standard output. Warnings and
  errors go to standard error; a failed command exits 1 with the reason.
  The fields are documented at
  https://github.com/fuyu-myk/pixelfold/blob/main/docs/json.md.
- Never run the interactive commands: `pixelfold FILE`, `pixelfold view`,
  `poses --view` or `pixelfold fetch`. They need a terminal and refuse without
  one. To get an entry, give its PDB id to any report.
- Read the warnings. They say when a ligand lacks chemical definitions, a
  selection covers several copies, or a pose gets no RMSD, and why.
- `pixelfold COMMAND --help` ends with worked examples.

## Selections

`--select` takes a query: `chain A and resi 10-40`, `resn STI`,
`byres within 5 of resn STI`, `protein`, `ligand`, `water`, `metal`,
`ss H`, `plddt < 70`, combined with `and`, `or`, `not` and parentheses.
`--ligand CODE` is shorthand for `--select 'resn CODE'`. The `ligand` keyword
leaves out water and metal ions; select a metal with `metal` or `resn ZN`.

## Ligand interactions

```bash
pixelfold interactions 4DFR --ligand MTX --format json
pixelfold interactions 4DFR --type hydrogen-bond --type salt-bridge --select 'chain A' --format json
```

Each record names both partners (`chain_a`, `resi_a`, `resn_a`, `atoms_a`, and
the same for `b`), the `kind`, the `distance` in angstroms and, for the kinds
judged by one, the `angle`. Kinds: hydrogen-bond, salt-bridge, hydrophobic,
pi-stacking, pi-cation, pi-hydrogen-bond, halogen-bond, water-bridge,
metal-coordination, disulfide.

## Ranking ligand poses

```bash
pixelfold poses --reference crystal.cif --ligand STI model_*.cif --format json
```

- Every model of every pose file is a pose; `--model N` keeps one.
- Residues are matched by aligning chain sequences, so a co-folding model may
  number and name chains as it likes. A pose file without the ligand's chemical
  definitions borrows the reference's by atom name, so keep the reference's
  ligand code and atom names where possible.
- When the reference holds several copies of the ligand, pick one:
  `--reference-select 'chain A and resn STI'`. Scoring a single-ligand pose
  against every copy caps its recovery.
- SD files of docked ligands (`.sdf`, `.sdf.gz`, `.mol`) are poses too, each
  record laid into the reference without its ligand, or into the protein given
  by `--receptor receptor.pdb`. A record's atoms take the reference ligand's
  names by bond graph, so the molecule must be the reference's.

Records come best first: most `similarity` (Tanimoto of residue-level
contacts), then `recovery`, then fewest `clashes`. Read with them:

- `rmsd`: heavy-atom RMSD in angstroms on the reference's pocket; 2 A is the
  usual success cutoff. `null` means no pocket or no atom pairing; a warning
  says which.
- `inverted`: reference stereocentres and E/Z double bonds the pose turns the
  other way. Non-empty means a different molecule, however well it sits.
- `geometry`: bonds, angles, internal clashes and non-flat rings or double
  bonds that stray from the reference ligand's own geometry. All zero is a
  sound pose.
- `lost`: the reference contacts the pose misses, e.g. `ASP A/381 salt-bridge`.

## Checking poses without a reference

```bash
pixelfold poses --ligand STI fold_job_model_*.cif --format json
pixelfold poses --receptor receptor.pdb --ligand STI docked.sdf --format json
```

- Leave out `--reference` when no crystal structure of the complex exists.
  Each pose is checked as PoseBusters checks a docked one; there is no RMSD and
  no contact recovery.
- Records come best first: fewest `failed` checks (an `unchecked` one counts
  as failed), then the highest `confidence.score`, read from AlphaFold 3's or Boltz's JSON beside each model
  file (keep the tools' file names), then the order given.
- `failed` names the checks failed: `clashes` (with the protein), `bonds`,
  `angles`, `internal_clashes`, `flatness`, `stereo`. Empty means the pose
  passes everything that could be checked; read `unchecked` too.
- `stereo` lists the centres and E/Z bonds turned the other way from the
  Chemical Component Dictionary entry of the ligand's code, fetched once and
  cached. Use the real component code (`--ligand STI`): a ligand named `LIG`
  from a SMILES input has no entry its atoms match, so its stereochemistry,
  and without definitions in the file its geometry, go `unchecked`.
- `.sdf` poses need `--receptor`; the receptor's own copy of the `--ligand`
  component is left out. Without `--ligand` an SD record is checked on its own
  bonds only.
- A covalent ligand's bond to the protein counts as clashes here.

## Residue interaction networks

```bash
pixelfold rin 4DFR --format json
pixelfold rin 4DFR --analyze --format json
pixelfold rin 4DFR --path A/23 A/94 --format json
pixelfold rin 1D3Z --model all --min-frequency 0.9 --format json
pixelfold rin sample_*.cif --min-frequency 0.5 --format json
```

- Residue ids are `chain/resi` with any insertion code: `A/60B`.
- `--analyze` ranks every residue by betweenness (hubs) and flags cut residues
  (`articulation`); `--path` gives the shortest route of interactions.
- `--model all` folds every model of a file into one network, and several files
  (a prediction's samples) fold the same way. Each edge carries `models` and
  `frequency`; `--min-frequency` keeps the frequent ones. The files must number
  residues alike: a residue named differently is refused.
- Other outputs: `--format graphml -o network.graphml` for Cytoscape or Gephi,
  `--format tsv`, and `--mvsj -o scene.mvsj`, a Mol* scene of the structure with
  the network's residues highlighted.

## Per-residue reports

```bash
pixelfold ss 1UBQ --format json
pixelfold sasa 1UBQ --format json
pixelfold ss 1D3Z --model all --format json
```

`ss` gives each residue's secondary structure (`H`, `E`, `T`, `C`, from
DSSP); `sasa` gives each polymer residue's accessible area in square
angstroms and `rsa`, relative to the residue's maximum. Under `--model all`,
each record starts with `model`.

## Figures

```bash
pixelfold render 4HHB -o hemoglobin.png --width 1600 --height 1200
pixelfold render 4DFR -o pocket.png --select 'byres within 6 of resn MTX' --color ss
```

Always pass `-o`: without it the image is drawn to the terminal.

## Other options

- `--model N|all`: which model of a multi-model file (the first by default).
- `--assembly`, `--assembly=ID`: the biological assembly instead of the
  deposited coordinates; copies of chain A are named `A-2`, `A-3`.
- `--altloc occupancy|a|b|all`: alternate conformers (the most occupied by
  default).
- A PDB-format file carries no chemical component definitions, so its ligand
  contacts are incomplete (a warning says so); prefer mmCIF.
