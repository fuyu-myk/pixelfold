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
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

import freesasa
import tomllib
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBExceptions import PDBConstructionException

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


def _optional(value: str) -> str | None:
    """Collapse mmCIF's absent markers (``.`` / ``?``) to ``None``."""
    return None if value in (".", "?", "") else value


def dssp_golden(pdb_id: str, cif: Path) -> dict:
    """Per-residue SS and backbone H-bonds from DSSP 4.

    Parses mkdssp's mmCIF output rather than its legacy fixed-column format:
    the legacy format cannot represent multi-character chain ids (e.g. 7P5R's
    ``AAA``), and MMCIF2Dict is a tokenizer, so it also tolerates altloc-only
    residues that trip Biopython's structure builder (e.g. 1EJG). This
    reproduces the Biopython/legacy output exactly on structures both handle.
    """
    with tempfile.NamedTemporaryFile(suffix=".cif") as out:
        subprocess.run(
            ["mkdssp", "--output-format=mmcif", str(cif), out.name],
            check=True,
            capture_output=True,
            text=True,
        )
        dssp = MMCIF2Dict(out.name)

    bridge = "_dssp_struct_bridge_pairs."
    if bridge + "auth_asym_id" not in dssp:
        # No protein residues (e.g. a pure nucleic-acid entry): DSSP emits no
        # struct categories at all.
        return {"residues": [], "hbonds": []}

    # Secondary structure is keyed by label ids ('.'/'-' is coil).
    summary = "_dssp_struct_summary."
    ss_by_label = {
        (asym, seq): ("C" if ss in (".", "-") else ss)
        for asym, seq, ss in zip(
            dssp[summary + "label_asym_id"],
            dssp[summary + "label_seq_id"],
            dssp[summary + "secondary_structure"],
        )
    }

    # One bridge-pairs row per residue carries its author id plus its H-bond
    # partners; acceptor_1/2 are the C=O that this residue's N-H donates to.
    residues, hbonds = [], []
    rows = zip(
        dssp[bridge + "auth_asym_id"],
        dssp[bridge + "auth_seq_id"],
        dssp[bridge + "pdbx_PDB_ins_code"],
        dssp[bridge + "label_asym_id"],
        dssp[bridge + "label_seq_id"],
    )
    for i, (asym, seq, icode, label_asym, label_seq) in enumerate(rows):
        donor = residue_key(asym, int(seq), _optional(icode) or "")
        residues.append({**donor, "ss": ss_by_label.get((label_asym, label_seq), "C")})
        for slot in (bridge + "acceptor_1", bridge + "acceptor_2"):
            partner = _optional(dssp[slot + "_auth_asym_id"][i])
            energy = _optional(dssp[slot + "_energy"][i])
            if partner is None or energy is None or float(energy) >= HBOND_ENERGY_CUTOFF:
                continue
            acceptor = residue_key(
                partner,
                int(dssp[slot + "_auth_seq_id"][i]),
                _optional(dssp[slot + "_pdbx_PDB_ins_code"][i]) or "",
            )
            hbonds.append({"donor": donor, "acceptor": acceptor})

    return {"residues": residues, "hbonds": hbonds}


def _freesasa_structure_from_dict(cif: Path) -> freesasa.Structure:
    """Build a FreeSASA structure straight from the mmCIF via MMCIF2Dict.

    Fallback for files Biopython's structure builder rejects (e.g. 1EJG, whose
    Ser22 has only altloc B/C and no primary conformer). Mirrors
    ``structureFromBioPDB``'s defaults, first model, skip HETATM/water and
    hydrogens, keep one conformer per atom, and yields identical SASA on
    inputs both paths accept.
    """
    atoms = MMCIF2Dict(str(cif))
    site = "_atom_site."
    structure = freesasa.Structure()
    seen: set[tuple[str, str, str]] = set()
    first_model = atoms[site + "pdbx_PDB_model_num"][0]
    rows = zip(
        atoms[site + "group_PDB"],
        atoms[site + "label_atom_id"],
        atoms[site + "type_symbol"],
        atoms[site + "auth_comp_id"],
        atoms[site + "auth_seq_id"],
        atoms[site + "auth_asym_id"],
        atoms[site + "pdbx_PDB_model_num"],
        atoms[site + "Cartn_x"],
        atoms[site + "Cartn_y"],
        atoms[site + "Cartn_z"],
    )
    for group, name, element, comp, seq, asym, model, x, y, z in rows:
        if group != "ATOM" or element == "H" or model != first_model:
            continue
        atom = name.strip('"')
        conformer = (asym, seq, atom)
        if conformer in seen:
            continue
        seen.add(conformer)
        structure.addAtom(f"{atom:<4.4}", f"{comp:>3.3}", seq, asym, float(x), float(y), float(z))
    return structure


def freesasa_golden(pdb_id: str, cif: Path) -> dict:
    """Per-residue SASA (A^2) from FreeSASA.

    FreeSASA parses only PDB natively, so build its structure from the
    Biopython-parsed mmCIF. The defaults skip waters/HETATM and hydrogens and
    use the first model, matching how the DSSP path reads the structure.
    """
    freesasa.setVerbosity(freesasa.nowarnings)
    try:
        structure = MMCIFParser(QUIET=True).get_structure(pdb_id, str(cif))
        fs_structure = freesasa.structureFromBioPDB(structure)
    except PDBConstructionException:
        fs_structure = _freesasa_structure_from_dict(cif)
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
