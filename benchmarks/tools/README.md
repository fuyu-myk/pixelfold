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
- `benchmarks/golden/<ID>.ring.json`: per-type non-covalent interactions of the
  whole structure, from the **RING 4.0** web server (same shape without
  `ligands`, plus the scope fields `water_excluded` and `sequence_separation` and
  a `source` record).

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

RING 4.0 is the whole-structure reference, and it is not dockerised. Its
standalone binary is free only to academic users, on request, under a licence
issued to one person that forbids passing it on, so it cannot be built into an
image others pull. `generate_ring.py` instead submits each structure to the RING
4.0 web server (`https://ring.biocomputingup.it`, through its public job API at
`drmaatic.biocomputingup.it`, the same one RING's PyMOL plugin uses), which needs
no account. It runs on the host with the standard library alone (Python 3.11 or
newer) and sends one job at a time. The server reports no version, so each golden
records the server task (`ring-v4`), the parameters sent, and the date it ran.
Regenerating depends on that server staying up and keeping its behaviour, which
is the price of a reference anyone can run; the committed goldens keep the
comparison itself reproducible either way.

RING runs at its defaults over the whole structure, ligands included. Two of those
defaults narrow what it looks at, and the golden records both: it leaves water out
(`water_excluded`), and it never pairs two standard amino acids of one chain
fewer than three positions apart (`sequence_separation: 3`). Positions count
residues rather than author numbers, so an insertion code widens the distance,
and ligands, modified residues and nucleotides pair at any distance. The harness
drops pixelfold's edges outside that scope before scoring, as it does for PLIP's
ligand scope. Every server reply is kept under `benchmarks/ring-cache/`
(gitignored) and parsed in place of a new job, so changing the type mapping never
costs the server a rerun.

RING's types map onto pixelfold's as follows; the rest are left out of the golden
and unscored.

| RING | pixelfold | Note |
| --- | --- | --- |
| `HBOND` | `hydrogen-bond` | |
| `IONIC` | `salt-bridge` | |
| `PIPISTACK` | `pi-stacking` | |
| `PICATION` | `pi-cation` | |
| `HALOGEN` | `halogen-bond` | |
| `METAL_ION` | `metal-coordination` | |
| `SSBOND` | `disulfide` | |
| `VDW` | none | any two atoms whose van der Waals radii overlap; pixelfold's `hydrophobic` is two apolar carbons within 4.0 A, a different criterion |
| `PIHBOND` | none | pi-hydrogen bonds, which pixelfold does not detect |
| `IAC` | none | generic contacts, emitted only under RING's C-alpha, C-beta and lollipop node policies, not the closest-atom policy used here |
| none | `water-bridge` | RING has no water-mediated type |

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

# RING 4.0 (web server, no image; reads the mmCIF the harness caches)
python3 benchmarks/tools/generate_ring.py --golden-dir benchmarks/golden
```

`generate_ring.py` skips any entry whose server reply is already kept; pass
`--force` to resubmit. The server meters anonymous jobs: a fresh address gets
about twenty at once, then one every ten minutes as its credit refills, so a
whole-benchmark run from an empty `ring-cache` takes about a day and a half. It is
resumable, so stopping and rerunning picks up where it left off.

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
- **RING residue alignment.** RING reads the same mmCIF and names residues by
  author chain id (label chain ids only on request, which the generator never
  sends), so its node ids (`chain:seq:icode:name`, `_` for no insertion code)
  map straight onto the harness's residue keys.

If a version pin in a `Dockerfile` fails to resolve, bump it, rebuild, and note
the change with the regenerated golden files.
