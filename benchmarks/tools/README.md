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
  surface area of the polymer, from **FreeSASA**. The polymer is every residue
  of a polymer entity, so modified residues such as selenomethionine count even
  though they are HETATM records; waters, ligands and ions do not.
- `benchmarks/golden/<ID>.plip.json`: per-type non-covalent interactions of the
  structure's ligands, from **PLIP** (`{kinds, ligand_scoped, ligands,
  interactions}`).

Structures are cached (and reused) under `benchmarks/cache/` (gitignored).

There are two images, one per chemistry backend: `Dockerfile` for DSSP + FreeSASA
and `Dockerfile.plip` for PLIP. PLIP is the protein-ligand reference: it profiles
ligand binding sites, so a `.plip.json` covers only interactions that touch a
ligand, which is what `ligand_scoped` records, and `ligands` lists the residues
PLIP profiled. The harness restricts pixelfold to its edges touching one of those
before scoring, so a whole-structure engine is charged neither for the
protein-protein interactions PLIP never looks at nor for the HETATM groups PLIP
does not treat as ligands (modified residues of a chain, and crystallisation
additives it sets aside in some entries). The harness refuses a ligand-scoped
file with an interaction whose ligand is not listed. PLIP does not
detect disulfides, so that type is absent from `kinds` and left unscored.

RING 4.0 (the protein-protein reference, whole-structure) is not generated here:
its engine is a license-gated closed-source binary. The harness already reads a
`<ID>.ring.json` of the same shape once one is produced.

Ids prefixed `AF-` (for example `AF-P01542`) are AlphaFold DB UniProt
accessions: the model is fetched from AlphaFold DB instead of the PDB, and its
golden files are named `AF-P01542.dssp.json` / `.freesasa.json` like any other
entry.

## Run it

From the repository root:

```bash
# DSSP + FreeSASA
docker build -t pixelfold-golden benchmarks/tools
docker run --rm -v "$PWD/benchmarks:/work/benchmarks" pixelfold-golden \
    --manifest benchmarks/manifest.toml --golden-dir benchmarks/golden

# PLIP (reads the legacy PDB, cached separately from the mmCIF the harness uses)
docker build -f benchmarks/tools/Dockerfile.plip -t pixelfold-plip benchmarks/tools
docker run --rm -v "$PWD/benchmarks:/work/benchmarks" pixelfold-plip \
    --golden-dir benchmarks/golden --structures benchmarks/pdb-cache
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
  same atom set for the MAE to be meaningful. Both take the polymer entities'
  heavy atoms (falling back to ATOM records for a file without entities);
  adjust the generator and the harness together if that changes. FreeSASA's
  ProtOr classifier has no radii for most modified residues (cysteine sulfonic
  acid, D-amino acids, 8A71's modified DNA), so those atoms take the element
  radius FreeSASA's reader guesses for them; pixelfold assigns every radius by
  element.
- **PLIP residue alignment.** PLIP reads the legacy PDB, whose author chain and
  residue numbering match the mmCIF's. Confirm a known ligand's contacts line up
  with `pixelfold interactions <id> --ligand <code>`: on 1IEP, PLIP and pixelfold
  agree residue-for-residue on imatinib's hydrogen bonds, salt bridge, and stack.
  A structure too large for PDB format has no `.pdb` and is skipped; those are
  never protein-ligand entries.

If a version pin in a `Dockerfile` fails to resolve, bump it, rebuild, and note
the change with the regenerated golden files.
