#!/usr/bin/env python3
"""Generate PLIP interaction golden files for pixelfold-validate.

For every id in the benchmark manifest, this runs PLIP over the ligand binding
sites and writes the per-type interaction edges pixelfold-validate compares
against:

    <golden-dir>/<ID>.plip.json   {"kinds": [...], "interactions": [{kind, a, b}]}

PLIP is a protein-ligand interaction profiler, so a golden file covers only the
interactions of a structure's ligands, not protein-protein contacts. Each edge
is an unordered residue pair: the binding-site residue on one side (PLIP's own
author numbering) and the ligand on the other. Residue keys use the author chain
id and sequence number to match how pixelfold reads structures; verify the
alignment against `pixelfold interactions <id>` on the first run.

Run inside the pinned Docker image (see Dockerfile.plip).
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import tomllib
from plip.structure.preparation import PDBComplex

# PLIP reads PDB, not mmCIF. The legacy PDB file carries the same author chain
# and residue numbering the mmCIF does, so PLIP's residues still align with
# pixelfold's. A structure too large for PDB format has no .pdb and is skipped;
# those are never protein-ligand entries, which is all PLIP scores.
RCSB_PDB = "https://files.rcsb.org/download/{id}.pdb"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

# The interaction-set attribute lists PLIP fills, mapped to pixelfold's type
# labels.
PLIP_LISTS: dict[str, tuple[str, ...]] = {
    "hydrogen-bond": ("hbonds_ldon", "hbonds_pdon"),
    "salt-bridge": ("saltbridge_lneg", "saltbridge_pneg"),
    "hydrophobic": ("hydrophobic_contacts",),
    "pi-stacking": ("pistacking",),
    "pi-cation": ("pication_laro", "pication_paro"),
    "halogen-bond": ("halogen_bonds",),
    "water-bridge": ("water_bridges",),
    "metal-coordination": ("metal_complexes",),
}


def structure_url(pdb_id: str) -> str:
    if pdb_id.upper().startswith("AF-"):
        accession = pdb_id.split("-", 1)[1]
        with urllib.request.urlopen(ALPHAFOLD_API.format(accession=accession)) as response:
            predictions = json.loads(response.read())
        if not predictions:
            raise ValueError(f"no AlphaFold prediction for {accession}")
        return predictions[0]["pdbUrl"]
    return RCSB_PDB.format(id=pdb_id)


def fetch_structure(pdb_id: str, structures: Path) -> Path:
    path = structures / f"{pdb_id}.pdb"
    if not path.exists():
        url = structure_url(pdb_id)
        print(f"  fetching {url}", file=sys.stderr)
        urllib.request.urlretrieve(url, path)
    return path


def residue_key(chain: str, seq: int, icode: str | None = None) -> dict:
    icode = (icode or "").strip()
    return {"chain": chain, "seq": int(seq), "icode": icode or None}


def edges_for(pdb_id: str, pdb: Path) -> list[dict]:
    """Every PLIP interaction as `{kind, a, b}`, an unordered residue pair."""
    mol = PDBComplex()
    mol.load_pdb(str(pdb))
    mol.analyze()

    seen: set[tuple[str, str, str]] = set()
    edges: list[dict] = []
    for site in mol.interaction_sets.values():
        for kind, lists in PLIP_LISTS.items():
            for attr in lists:
                for interaction in getattr(site, attr, []) or []:
                    # PLIP names both endpoints on every interaction: the binding-
                    # site residue in reschain/resnr, and the ligand (or, for a
                    # metal complex, the metal ion) in reschain_l/resnr_l. Take the
                    # ligand from those fields, not the binding-site key. Neither side
                    # exposes an insertion code, so it stays None.
                    partner = residue_key(interaction.reschain, interaction.resnr)
                    ligand = residue_key(interaction.reschain_l, interaction.resnr_l)
                    a, b = sorted(
                        (partner, ligand),
                        key=lambda r: (r["chain"], r["seq"], r["icode"] or ""),
                    )
                    fingerprint = (kind, json.dumps(a, sort_keys=True), json.dumps(b, sort_keys=True))
                    if fingerprint in seen:
                        continue
                    seen.add(fingerprint)
                    edges.append({"kind": kind, "a": a, "b": b})

    return edges


def plip_golden(pdb_id: str, pdb: Path) -> dict:
    return {
        "kinds": list(PLIP_LISTS),
        "ligand_scoped": True,
        "interactions": edges_for(pdb_id, pdb),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("benchmarks/manifest.toml"))
    parser.add_argument("--golden-dir", type=Path, default=Path("benchmarks/golden"))
    parser.add_argument("--structures", type=Path, default=Path("benchmarks/pdb-cache"))
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
        print(pdb_id, file=sys.stderr)
        try:
            pdb = fetch_structure(pdb_id, args.structures)
            golden = plip_golden(pdb_id, pdb)
            (args.golden_dir / f"{pdb_id}.plip.json").write_text(json.dumps(golden, indent=1))
            print(f"  {len(golden['interactions'])} interactions", file=sys.stderr)
        except Exception as error:  # noqa: BLE001
            print(f"  FAILED: {error}", file=sys.stderr)
            failures += 1

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
