# JSON output

Every command that reports data writes JSON with `--format json`, to standard
output or, for `rin`, to the file `-o` names. Warnings and errors go to
standard error, and a failed command exits with status 1, so standard output
holds only the result. Field names and their meanings hold within a major
version; a minor release may add fields. A field a record cannot fill is `null`.

## Shared fields

A residue is named by four fields, or, in `rin`, by one id:

| Field | Meaning |
|---|---|
| `chain` | The author chain id; in a biological assembly, a copy of chain A is `A-2`, `A-3`, ... |
| `resi` | The author residue number |
| `icode` | The insertion code, or `null` |
| `resn` | The component code (`ALA`, `HOH`, `STI`) |
| `id` | `chain/resi` with any insertion code appended: `A/60B` |
| `ss` | Secondary structure from DSSP, folded to four letters: `H` helix, `E` strand, `T` turn, `C` coil |

With `--model all`, every record of `interactions`, `ss` and `sasa` starts with
`model`, the model's number in the file.

## `interactions`

An array of interactions:

| Field | Meaning |
|---|---|
| `kind` | `hydrogen-bond`, `salt-bridge`, `hydrophobic`, `pi-stacking`, `pi-cation`, `pi-hydrogen-bond`, `halogen-bond`, `water-bridge`, `metal-coordination` or `disulfide` |
| `chain_a`, `resi_a`, `icode_a`, `resn_a` | The first partner's residue |
| `atoms_a` | Its atom names, joined by `+` for a group (a ring, a charged group) |
| `chain_b`, `resi_b`, `icode_b`, `resn_b`, `atoms_b` | The second partner, the same way |
| `distance` | In angstroms: atom to atom, or to a group's centre |
| `angle` | In degrees, for hydrogen bonds, halogen bonds, water bridges, pi-stacking and pi-hydrogen bonds, the angle each is judged by; else `null` |

## `ss`

An array of residues, each with the shared fields and `ss`.

## `sasa`

An array of polymer residues (modified residues included, ligands and water
left out):

| Field | Meaning |
|---|---|
| `sasa` | Solvent-accessible surface area, in square angstroms |
| `rsa` | The same relative to the residue's maximum (Tien et al., 2013), or `null` for a residue without one |

## `rin`

The network, as one object:

| Field | Meaning |
|---|---|
| `models` | How many models the network spans: 1, or more under `--model all` or with several files |
| `nodes` | Residues: the shared fields, `id`, `ss` (the most common letter across models) and `degree` |
| `edges` | One per residue pair per kind of interaction: `source` and `target` ids, `kind`, `count` (interactions, summed over models), `min_distance` (angstroms), `interchain`, `energy` (kJ/mol, summed), `models` (how many models hold it) and `frequency` (`models` over the network's, to three places) |

### `rin --analyze --format json`

| Field | Meaning |
|---|---|
| `models`, `edges` | As above |
| `components` | Each connected component's residue ids, the largest first |
| `residues` | Every residue, the highest betweenness first: `id`, `resn`, `ss`, `degree`, `betweenness` (the shortest paths between other residues through it) and `articulation` (whether removing it splits its component) |

### `rin --path FROM TO --format json`

| Field | Meaning |
|---|---|
| `from`, `to` | The two ids asked for |
| `steps` | Interactions along the shortest route, or `null` when the two are in different components |
| `path` | The residues along it in order, each `id`, `resn` and `ss`, or `null` |

## `poses`

An array of poses, the best first (most similar to the reference's contacts,
then the most recovered, then the fewest clashes):

| Field | Meaning |
|---|---|
| `pose`, `model` | The pose file as given, and the model (or SD record) number |
| `copy` | The ligand copy scored, by its first residue, when the pose holds several, else `null` |
| `contacts`, `shared` | The pose ligand's residue-level contacts, and how many the reference also has |
| `recovery` | `shared` over the reference's contacts |
| `similarity` | Tanimoto similarity of the two sets of contacts |
| `rmsd` | Heavy-atom RMSD in angstroms once the pose is superposed on the reference's pocket, or `null` |
| `clashes` | Ligand-protein heavy-atom pairs closer than 0.75 of their van der Waals radii summed |
| `inverted` | The reference's stereocentres (`STI A/201 C20`) and E/Z double bonds (`GPP A/101 C2=C3`) the pose turns the other way, or `null` when its atoms could not be paired with the reference's |
| `geometry` | What of the ligand's own geometry strays from the reference ligand's: `bonds`, `angles`, `internal_clashes`, `unflat` (aromatic rings and double bonds not flat) and `flat_rings` (saturated six-membered rings the reference holds puckered left flat), the last two `null` without the reference's bond orders; `null` when the atoms could not be paired |
| `lost` | The reference's contacts the pose lacks, each `RESN chain/resi kind` |

### `poses` without a reference

Without `--reference`, each pose is checked on its own, as PoseBusters checks a
docked pose, and the array comes best first: the fewest checks failed or
unchecked, then the fewest failed, then the highest `confidence.score`, then the
order given.

| Field | Meaning |
|---|---|
| `pose`, `model` | The pose file as given, and the model (or SD record) number |
| `failed` | The checks the pose fails: `clashes`, `metal_clashes`, `water_clashes`, `too_far`, `bonds`, `angles`, `internal_clashes`, `flatness`, `flat_rings`, `stereo` |
| `unchecked` | The checks that could not be made: `geometry` (no bond table for some ligand residue: no definition, dictionary entry or SMILES pairs with it), `stereo` (no dictionary entry its atoms match and no SMILES that pairs) |
| `contacts` | The ligand's residue-level contacts, a count; not ranked |
| `clashes` | Heavy-atom pairs of the ligand and the protein or an organic cofactor closer than 0.75 of their van der Waals radii summed, the covalent bonds the pose file records in `_struct_conn` (and pairs one bond from them) left out |
| `metal_clashes` | Pairs of the ligand and a metal or inorganic ion (sulfate, phosphate, an iron-sulfur cluster, or any chlorine or bromine atom, as PoseBusters takes them) closer than 0.75 of their covalent radii summed, which a coordination bond never is |
| `water_clashes` | Pairs of the ligand and a water oxygen closer than 0.75 of their van der Waals radii summed |
| `protein_distance` | The shortest distance in angstroms from a ligand heavy atom to a polymer one; over 5 fails `too_far`, as does `null` (no polymer) |
| `geometry` | The ligand's own geometry against the bounds RDKit gives its bonds: `bonds`, `angles`, `internal_clashes`, `unflat`, `flat_rings` (each a count); `null` when unchecked |
| `stereo` | The stereocentres (`STI A/201 C20`) and E/Z double bonds (`GPP A/101 C2=C3`) the pose turns the other way from those the dictionary labels or the SMILES string (`--smiles`, or the model's own) states; `null` when no residue could be checked |
| `confidence` | `null`, or the co-folding tool's confidence found beside the model file: `source` (`alphafold3` or `boltz`), `score` (`ranking_score` or `confidence_score`), `ligand_iptm` (the ligand chain's `chain_iptm`, or `ligand_iptm`; `null` when absent) |
