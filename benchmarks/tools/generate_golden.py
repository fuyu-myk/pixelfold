#!/usr/bin/env python3
"""Generate golden reference files for pixelfold-validate.

For every PDB id in the benchmark manifest, this fetches the mmCIF, runs DSSP 4
(``mkdssp``) and FreeSASA, and writes the reference outputs pixelfold-validate
consumes:

    <golden-dir>/<ID>.dssp.json      per-residue secondary structure + H-bonds
    <golden-dir>/<ID>.freesasa.json  per-residue solvent-accessible surface area

Run it inside the pinned Docker image (see Dockerfile) so the tool versions are
reproducible. Residue keys use the author chain id and sequence number to match
how pixelfold reads structures; verify the alignment on the first run.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import freesasa
import tomllib
from Bio.PDB import DSSP, MMCIFParser

RCSB_CIF = "https://files.rcsb.org/download/{id}.cif"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
HBOND_ENERGY_CUTOFF = -0.5  # kcal/mol, DSSP's definition


def structure_url(pdb_id: str) -> str:
    """Resolve the mmCIF download URL for a benchmark id.

    Ids prefixed ``AF-`` (for example ``AF-P01542``) are AlphaFold DB UniProt
    accessions; their file URL and version is resolved through the AlphaFold
    DB API. Every other id is a PDB id served by RCSB.
    """
    if pdb_id.upper().startswith("AF-"):
        accession = pdb_id.split("-", 1)[1]
        with urllib.request.urlopen(ALPHAFOLD_API.format(accession=accession)) as response:
            predictions = json.loads(response.read())
        if not predictions:
            raise ValueError(f"no AlphaFold prediction for {accession}")
        return predictions[0]["cifUrl"]
    return RCSB_CIF.format(id=pdb_id)


def fetch_structure(pdb_id: str, structures: Path) -> Path:
    path = structures / f"{pdb_id}.cif"
    if not path.exists():
        url = structure_url(pdb_id)
        print(f"  fetching {url}", file=sys.stderr)
        urllib.request.urlretrieve(url, path)
    return path


def residue_key(chain: str, resnum: int, icode: str) -> dict:
    icode = icode.strip()
    return {"chain": chain, "seq": resnum, "icode": icode or None}


def dssp_golden(pdb_id: str, cif: Path) -> dict:
    """Per-residue SS and backbone H-bonds from DSSP 4, via Biopython."""
    structure = MMCIFParser(QUIET=True).get_structure(pdb_id, str(cif))
    model = next(structure.get_models())
    dssp = DSSP(model, str(cif), dssp="mkdssp")

    # DSSP indexes residues 1..N; map that index to a residue key so H-bond
    # relative offsets can be resolved back to concrete residues.
    keys = list(dssp.keys())
    index_to_key: dict[int, dict] = {}
    residues = []
    for key in keys:
        chain = key[0]
        _hetflag, resnum, icode = key[1]
        record = dssp[key]
        dssp_index, _aa, ss = record[0], record[1], record[2]
        rkey = residue_key(chain, resnum, icode)
        index_to_key[dssp_index] = rkey
        # DSSP uses '-' for coil.
        residues.append({**rkey, "ss": ss if ss != "-" else "C"})

    hbonds = []
    for key in keys:
        record = dssp[key]
        donor_index = record[0]
        # NH-->O(1) and NH-->O(2): this residue's N-H donates to a partner C=O.
        for rel_idx_pos, energy_pos in ((6, 7), (10, 11)):
            rel, energy = record[rel_idx_pos], record[energy_pos]
            if rel == 0 or energy >= HBOND_ENERGY_CUTOFF:
                continue
            acceptor = index_to_key.get(donor_index + rel)
            if acceptor is not None:
                hbonds.append({"donor": index_to_key[donor_index], "acceptor": acceptor})

    return {"residues": residues, "hbonds": hbonds}


def freesasa_golden(pdb_id: str, cif: Path) -> dict:
    """Per-residue SASA (A^2) from FreeSASA.

    FreeSASA parses only PDB natively, so build its structure from the
    Biopython-parsed mmCIF. The defaults skip waters/HETATM and hydrogens and
    use the first model, matching how the DSSP path reads the structure.
    """
    freesasa.setVerbosity(freesasa.nowarnings)
    structure = MMCIFParser(QUIET=True).get_structure(pdb_id, str(cif))
    fs_structure = freesasa.structureFromBioPDB(structure)
    # Guard against 0-atom structures.
    if fs_structure.nAtoms() == 0:
        print("  no SASA-eligible atoms; writing empty SASA", file=sys.stderr)
        return {"residues": []}
    result = freesasa.calc(fs_structure)
    residues = []
    for chain, chain_residues in result.residueAreas().items():
        for resnum, area in chain_residues.items():
            # FreeSASA residue numbers may carry an insertion-code suffix.
            number, icode = resnum, ""
            if resnum and not resnum[-1].isdigit():
                number, icode = resnum[:-1], resnum[-1]
            residues.append({**residue_key(chain, int(number), icode), "sasa": area.total})
    return {"residues": residues}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("benchmarks/manifest.toml"))
    parser.add_argument("--golden-dir", type=Path, default=Path("benchmarks/golden"))
    parser.add_argument("--structures", type=Path, default=Path("benchmarks/cache"))
    parser.add_argument("ids", nargs="*", help="PDB ids to process (default: whole manifest)")
    args = parser.parse_args()

    args.golden_dir.mkdir(parents=True, exist_ok=True)
    args.structures.mkdir(parents=True, exist_ok=True)

    ids = args.ids
    if not ids:
        manifest = tomllib.loads(args.manifest.read_text())
        ids = [pid for stratum in manifest["strata"].values() for pid in stratum]

    failures = 0
    for pdb_id in (pid.upper() for pid in ids):
        print(f"{pdb_id}", file=sys.stderr)
        try:
            cif = fetch_structure(pdb_id, args.structures)
            (args.golden_dir / f"{pdb_id}.dssp.json").write_text(
                json.dumps(dssp_golden(pdb_id, cif), indent=1)
            )
            (args.golden_dir / f"{pdb_id}.freesasa.json").write_text(
                json.dumps(freesasa_golden(pdb_id, cif), indent=1)
            )
        except Exception as error:  # noqa: BLE001
            print(f"  FAILED: {error}", file=sys.stderr)
            failures += 1

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
