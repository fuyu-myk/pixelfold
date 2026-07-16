# Golden reference generation

The validation harness (`pixelfold-validate`) compares pixelfold's output against
committed golden files under `benchmarks/golden/`. Those golden files come from
the reference tools, run in a pinned Docker image so the results are
reproducible. This directory holds that image and the generator.

## What it produces

For each PDB id in `benchmarks/manifest.toml`:

- `benchmarks/golden/<ID>.dssp.json`: per-residue secondary structure and
  backbone hydrogen bonds, from **DSSP 4** (`mkdssp`).
- `benchmarks/golden/<ID>.freesasa.json`: per-residue solvent-accessible
  surface area, from **FreeSASA**.

Structures are cached (and reused) under `benchmarks/cache/` (gitignored).

Ids prefixed `AF-` (for example `AF-P01542`) are AlphaFold DB UniProt
accessions: the model is fetched from AlphaFold DB instead of the PDB, and its
golden files are named `AF-P01542.dssp.json` / `.freesasa.json` like any other
entry.

## Run it

From the repository root:

```bash
docker build -t pixelfold-golden benchmarks/tools
docker run --rm -v "$PWD/benchmarks:/work/benchmarks" pixelfold-golden \
    --manifest benchmarks/manifest.toml --golden-dir benchmarks/golden
```

To (re)generate a single entry, pass ids: `... pixelfold-golden 1UBQ 4HHB`.

Then commit the new/changed files under `benchmarks/golden/`, run the harness to
see the numbers, and update the README validation table (see the root README):

```bash
cargo run -p pixelfold-validate -- --report
```

## First-run checks

Residue keys are matched by author chain id and sequence number. On the first
run, spot-check that a couple of entries align (non-zero `Residues` in the
per-entry table, plausible Q3). Two common mismatch sources:

- **Author vs label identifiers** in mmCIF. Pixelfold reads author ids; the
  generator does too via Biopython, but confirm they agree for multi-chain
  entries.
- **Which atoms count toward SASA.** FreeSASA and pixelfold must include the
  same atom set (for example both including or both excluding waters) for the
  MAE to be meaningful; adjust the generator or the harness if they diverge.

If a version pin in the `Dockerfile` fails to resolve, bump it, rebuild, and note
the change with the regenerated golden files.
