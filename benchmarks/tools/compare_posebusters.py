#!/usr/bin/env python3
"""Measure how often pixelfold's reference-free pose checks agree with
PoseBusters', check by check, over the cases generate_posebusters.py wrote.

For each case it runs

    pixelfold poses --receptor CASE/receptor.pdb --ligand CODE CASE/poses.sdf \\
        --format json --cache-dir CASES

(the cases directory holds the dictionary entries under ccd/, so nothing is
fetched) and compares each pose's failed checks with PoseBusters' verdicts:

    bonds             Bond lengths
    angles            Bond angles
    internal_clashes  Internal steric clash
    flatness          Aromatic ring flatness and Double bond flatness
    flat_rings        Non-aromatic ring non-flatness
    clashes           Minimum distance to protein
    stereo (centres)  Tetrahedral chirality
    stereo (bonds)    Double bond stereochemistry

A pose a check could not be made of on either side (PoseBusters could not
sanitize it, or pixelfold had no dictionary entry) is left out of that check.
With --smiles each pose is checked against its ligand's canonical SMILES from
the dictionary entry (`--smiles`) instead of the entry itself (`--ligand`), as
a ligand without a code would be. With --entries it also reads each crystal
ligand from its entry's mmCIF file
(--cache), through the file's own component definitions, and counts how many
get their SD copy's verdicts. Runs on the host with the standard library
alone; pass the pixelfold binary with --pixelfold. Prints a markdown table,
and with --details each disagreement.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

# pixelfold's check, PoseBusters' columns (a pose passes when every one does),
# and how the row is labelled.
CHECKS = [
    ("bonds", ["bond_lengths"], "Bond lengths"),
    ("angles", ["bond_angles"], "Bond angles"),
    ("internal_clashes", ["internal_steric_clash"], "Internal steric clash"),
    ("flatness", ["aromatic_ring_flatness", "double_bond_flatness"], "Flatness (aromatic rings, double bonds)"),
    ("flat_rings", ["non-aromatic_ring_non-flatness"], "Non-aromatic ring non-flatness"),
    ("clashes", ["minimum_distance_to_protein"], "Minimum distance to protein"),
    ("centres", ["tetrahedral_chirality"], "Tetrahedral chirality"),
    ("double_bonds", ["double_bond_stereochemistry"], "Double bond stereochemistry"),
]

GEOMETRY = {"bonds", "angles", "internal_clashes", "flatness", "flat_rings"}


def pose_names(sdf: Path) -> list[str]:
    names, start = [], True
    for line in sdf.read_text().splitlines():
        if start:
            names.append(line.strip())
            start = False
        if line.strip() == "$$$$":
            start = True
    return names


def pixelfold_fails(record: dict, check: str) -> bool | None:
    """Whether pixelfold fails `check` for this record, or None when it
    could not make it."""
    if check in GEOMETRY:
        return None if record["geometry"] is None else check in record["failed"]
    if check == "clashes":
        return record["clashes"] > 0
    if record["stereo"] is None:
        return None
    labels = record["stereo"]
    if check == "centres":
        return any("=" not in label for label in labels)
    return any("=" in label for label in labels)


def canonical_smiles(ccd: Path) -> str | None:
    """The dictionary entry's canonical SMILES, OpenEye's where it gives one."""
    found = {}
    for line in ccd.read_text().splitlines():
        if "SMILES_CANONICAL" in line and '"' in line:
            program = "OpenEye" if "OpenEye" in line else "other"
            found.setdefault(program, line.rsplit('"', 2)[1])
    return found.get("OpenEye") or found.get("other")


def run_json(command: list[str]) -> list[dict] | None:
    run = subprocess.run(command, capture_output=True, text=True, check=False)
    return json.loads(run.stdout) if run.returncode == 0 else None


def entry_agreement(args) -> tuple[int, int]:
    """How many crystal ligands, read from their entry's mmCIF file with its own
    component definitions, get the verdicts their SD copy (record 1) gets."""
    verdicts = json.loads((args.cases / "verdicts.json").read_text())
    same = read = 0
    for case in sorted(verdicts):
        entry, code = case.split("_", 1)
        directory = args.cases / case
        common = ["--ligand", code, "--format", "json", "--cache-dir", str(args.cases)]
        docked = run_json([args.pixelfold, "poses", "--receptor", str(directory / "receptor.pdb"),
                           str(directory / "poses.sdf"), "--model", "1", *common])
        own = run_json([args.pixelfold, "poses", str(args.cache / f"{entry}.cif"), *common])
        if docked is None or own is None:
            print(f"{case}: pixelfold failed on the entry or its SD copy", file=sys.stderr)
            continue
        read += 1
        same += int(all(docked[0][field] == own[0][field] for field in ("geometry", "stereo")))
    return same, read


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cases", type=Path, default=Path("benchmarks/posebusters"))
    parser.add_argument("--pixelfold", default="target/release/pixelfold")
    parser.add_argument("--details", action="store_true")
    parser.add_argument("--entries", action="store_true", help="also read each crystal ligand from its entry")
    parser.add_argument(
        "--smiles", action="store_true",
        help="name each ligand by its dictionary entry's canonical SMILES (--smiles) rather than its code",
    )
    parser.add_argument("--cache", type=Path, default=Path("benchmarks/cache"))
    args = parser.parse_args()

    verdicts = json.loads((args.cases / "verdicts.json").read_text())
    counts = {check: {"poses": 0, "agree": 0, "both": 0, "posebusters": 0, "pixelfold": 0, "skipped": 0} for check, _, _ in CHECKS}
    disagreements = []
    for case, by_pose in sorted(verdicts.items()):
        directory = args.cases / case
        code = case.split("_", 1)[1]
        named = ["--ligand", code]
        if args.smiles:
            smiles = canonical_smiles(args.cases / "ccd" / f"{code}.cif")
            if smiles is None:
                print(f"{case}: no canonical SMILES", file=sys.stderr)
                continue
            named = ["--smiles", smiles]
        run = subprocess.run(
            [
                args.pixelfold, "poses", "--receptor", str(directory / "receptor.pdb"), *named,
                str(directory / "poses.sdf"), "--format", "json", "--cache-dir", str(args.cases),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if run.returncode != 0:
            print(f"{case}: pixelfold failed: {run.stderr.strip()}", file=sys.stderr)
            continue
        names = pose_names(directory / "poses.sdf")
        records = {names[record["model"] - 1]: record for record in json.loads(run.stdout)}
        for name, posebusters in by_pose.items():
            record = records[name]
            for check, columns, _ in CHECKS:
                outcomes = [posebusters.get(column) for column in columns]
                ours = pixelfold_fails(record, check)
                if any(outcome is None for outcome in outcomes) or ours is None or posebusters.get("sanitization") is False:
                    counts[check]["skipped"] += 1
                    continue
                theirs = not all(outcomes)
                tally = counts[check]
                tally["poses"] += 1
                if ours == theirs:
                    tally["agree"] += 1
                    tally["both"] += int(ours)
                elif theirs:
                    tally["posebusters"] += 1
                    disagreements.append((case, name, check, "PoseBusters fails, pixelfold passes"))
                else:
                    tally["pixelfold"] += 1
                    disagreements.append((case, name, check, "pixelfold fails, PoseBusters passes"))

    print(f"{len(verdicts)} ligands, {sum(len(poses) for poses in verdicts.values())} poses\n")
    print("| Check | Poses | Agreement | Both fail | PoseBusters only | pixelfold only |")
    print("|---|---:|---:|---:|---:|---:|")
    for check, _, label in CHECKS:
        tally = counts[check]
        share = 100 * tally["agree"] / tally["poses"] if tally["poses"] else 0.0
        print(
            f"| {label} | {tally['poses']} | {share:.1f}% | {tally['both']} | "
            f"{tally['posebusters']} | {tally['pixelfold']} |"
        )
    if args.entries:
        same, read = entry_agreement(args)
        print(f"\nRead from their entry files, {same} of {read} crystal ligands get the geometry and")
        print("stereochemistry verdicts of their SD copies.")
    if args.details:
        print()
        for case, name, check, how in disagreements:
            print(f"{case} {name} {check}: {how}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
