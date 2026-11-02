#!/usr/bin/env python3
"""Generate the PoseBusters verdicts pixelfold's reference-free pose checks are
compared against.

For every protein-ligand entry of the benchmark manifest, each drug-like ligand
(six or more heavy atoms, no metal, not a buffer or solvent) whose first copy
is complete becomes a case: the crystal ligand with the Chemical Component
Dictionary's bond orders and hydrogens added, and perturbed copies of it that
should fail one check or another, or sit at a check's edge:

    crystal        the crystal coordinates
    ideal          the dictionary's ideal conformer laid on the crystal
    mirrored       the crystal reflected (every stereocentre inverted)
    stretched-N    a terminal group pulled N A along its bond
    bent-N         a terminal group turned N degrees at its bond angle
    twisted        a rotatable bond turned to bring its ends closest
    turned         a rotatable bond turned halfway between
    puckered-N     an aromatic ring atom moved N A off the ring's plane
    flattened-N    a saturated six-membered ring pressed N of the way flat
    flipped        a labelled E/Z double bond turned 180 degrees
    shifted-N      the crystal moved N A towards the protein
    random-N       a random small rotation, shift and torsion

With --partners it takes instead the PARTNERS entries, whose ligands sit beside
metal ions, cofactors and waters, keeps those in the receptor, and moves the
crystal ligand for the distance checks:

    toward-K-N     the crystal moved N A towards the nearest atom of kind K
                   (ion, water, cofactor), as PoseBusters sorts them
    away-4.5       the crystal moved out from the protein's centre until no
    away-5.5       protein atom lies within 4.5 or 5.5 A of it
    away-far       10 A further out
    random-N       a random small rotation and shift

and writes, under --cases (benchmarks/posebusters by default, or
benchmarks/posebusters-partners with --partners):

    <ENTRY>_<CODE>/poses.sdf     every pose, its name the variant's
    <ENTRY>_<CODE>/receptor.pdb  the entry's polymer, no ligands or water (with
                                 --partners, everything but the ligand)
    ccd/<CODE>.cif               the dictionary entry (pixelfold's --cache-dir)
    verdicts.json                {case: {variant: {check: passed}}}

PoseBusters runs over the checks pixelfold makes, as its dock and redock
configurations set them: the dictionary's ideal conformer is the true molecule
whose InChI stereochemistry layers each pose's are compared with. The receptor
PDB is regenerated each run and not committed; the poses, entries and verdicts
are.

Run inside the pinned Docker image (see Dockerfile.posebusters).
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
import urllib.request
import zlib
from collections import deque
from pathlib import Path

import gemmi
import numpy as np
import tomllib
from posebusters import PoseBusters
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

RDLogger.DisableLog("rdApp.*")

RCSB_CIF = "https://files.wwpdb.org/download/{id}.cif.gz"
CCD_CIF = "https://files.rcsb.org/ligands/download/{code}.cif"
CCD_SDF = "https://files.rcsb.org/ligands/download/{code}_ideal.sdf"

# Crystallisation additives, buffers and solvents, which are not what a docking
# or co-folding pose is of.
ADDITIVES = {
    "GOL", "EDO", "PEG", "PGE", "PG4", "PG6", "1PE", "P6G", "MPD", "MRD", "SO4",
    "PO4", "ACT", "ACY", "DMS", "TRS", "MES", "EPE", "CIT", "FMT", "IMD", "BTB",
    "NO3", "SCN", "MLA", "IPA", "EOH", "MOH", "BOG", "LDA", "HEZ", "BU3",
}
METALS = {"MG", "FE", "ZN", "CU", "CO", "MN", "NI", "CA", "NA", "K", "SR", "LI", "MO", "W", "V"}

# Entries whose ligand sits beside metal ions, cofactors or waters, for the
# distance checks (--partners).
PARTNERS = [
    ("3HS4", "AZM"),  # zinc, glycerol, waters
    ("1A42", "BZU"),  # zinc, mercury
    ("1TLP", "RDF"),  # zinc, calcium
    ("1ATP", "ATP"),  # manganese
    ("2CPP", "CAM"),  # haem
    ("1RX2", "FOL"),  # NADP+, manganese
    ("1HSG", "MK1"),  # the flap water
]

# What PoseBusters takes for inorganic cofactors (posebusters/tools/protein.py).
ION_ELEMENTS = {"LI", "BE", "NA", "MG", "CL", "K", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "BR", "RB", "MO", "CD"}
ION_RESIDUES = {"FES", "MOS", "PO3", "PO4", "PPK", "SO3", "SO4", "VO4"}

# The checks pixelfold makes, as PoseBusters names its table's columns.
CHECKS = [
    "sanitization",
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "double_bond_flatness",
    "non-aromatic_ring_non-flatness",
    "minimum_distance_to_protein",
    "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors",
    "minimum_distance_to_waters",
    "protein-ligand_maximum_distance",
    "tetrahedral_chirality",
    "double_bond_stereochemistry",
]

DG = {
    "set15bounds": True,
    "scaleVDW": True,
    "doTriangleSmoothing": True,
    "useMacrocycle14config": False,
}

# The kinds of receptor atom PoseBusters' distance checks take one at a time.
PARTNER_TYPES = {"protein", "organic_cofactors", "inorganic_cofactors", "waters"}

LOAD = {"cleanup": False, "sanitize": False, "add_hs": False, "assign_stereo": False, "load_all": True}

CONFIG = {
    "modules": [
        {
            "name": "Loading",
            "function": "loading",
            "chosen_binary_test_output": ["mol_pred_loaded", "mol_true_loaded", "mol_cond_loaded"],
            "rename_outputs": {
                "mol_pred_loaded": "MOL_PRED loaded",
                "mol_true_loaded": "MOL_TRUE loaded",
                "mol_cond_loaded": "MOL_COND loaded",
            },
        },
        {
            "name": "Chemistry",
            "function": "rdkit_sanity",
            "chosen_binary_test_output": ["passes_rdkit_sanity_checks"],
            "rename_outputs": {"passes_rdkit_sanity_checks": "Sanitization"},
        },
        {
            "name": "Chemistry",
            "function": "identity",
            "parameters": {"inchi_options": "w"},
            "chosen_binary_test_output": ["stereo_dbond", "stereo_tetrahedral"],
            "rename_outputs": {
                "stereo_dbond": "Double bond stereochemistry",
                "stereo_tetrahedral": "Tetrahedral chirality",
            },
        },
        {
            "name": "Geometry",
            "function": "distance_geometry",
            "parameters": {
                "bound_matrix_params": DG,
                "threshold_bad_bond_length": 0.25,
                "threshold_bad_angle": 0.25,
                "threshold_clash": 0.3,
                "ignore_hydrogens": True,
                "sanitize": True,
            },
            "chosen_binary_test_output": [
                "bond_lengths_within_bounds",
                "bond_angles_within_bounds",
                "no_internal_clash",
            ],
            "rename_outputs": {
                "bond_lengths_within_bounds": "Bond lengths",
                "bond_angles_within_bounds": "Bond angles",
                "no_internal_clash": "Internal steric clash",
            },
        },
        {
            "name": "Ring flatness",
            "function": "flatness",
            "parameters": {
                "flat_systems": {
                    "aromatic_5_membered_rings_sp2": "[ar5^2]1[ar5^2][ar5^2][ar5^2][ar5^2]1",
                    "aromatic_6_membered_rings_sp2": "[ar6^2]1[ar6^2][ar6^2][ar6^2][ar6^2][ar6^2]1",
                },
                "threshold_flatness": 0.25,
            },
            "chosen_binary_test_output": ["flatness_passes"],
            "rename_outputs": {"flatness_passes": "Aromatic ring flatness"},
        },
        {
            "name": "Ring non-flatness",
            "function": "flatness",
            "parameters": {
                "check_nonflat": True,
                "flat_systems": {
                    "non-aromatic_6_membered_rings": "[C,O,S,N;R1]~1[C,O,S,N;R1][C,O,S,N;R1][C,O,S,N;R1][C,O,S,N;R1][C,O,S,N;R1]1",
                    "non-aromatic_6_membered_rings_db03_0": "[C;R1]~1[C;R1][C,O,S,N;R1]~[C,O,S,N;R1][C;R1][C;R1]1",
                    "non-aromatic_6_membered_rings_db03_1": "[C;R1]~1[C;R1][C;R1]~[C;R1][C,O,S,N;R1][C;R1]1",
                    "non-aromatic_6_membered_rings_db02_0": "[C;R1]~1[C;R1][C;R1][C,O,S,N;R1]~[C,O,S,N;R1][C;R1]1",
                    "non-aromatic_6_membered_rings_db02_1": "[C;R1]~1[C;R1][C,O,S,N;R1][C;R1]~[C;R1][C;R1]1",
                },
                "threshold_flatness": 0.05,
            },
            "chosen_binary_test_output": ["flatness_passes"],
            "rename_outputs": {"flatness_passes": "Non-aromatic ring non-flatness"},
        },
        {
            "name": "Double bond flatness",
            "function": "flatness",
            "parameters": {
                "flat_systems": {"trigonal_planar_double_bonds": "[C;X3;^2](*)(*)=[C;X3;^2](*)(*)"},
                "threshold_flatness": 0.25,
            },
            "chosen_binary_test_output": ["flatness_passes"],
            "rename_outputs": {"flatness_passes": "Double bond flatness"},
        },
        {
            "name": "Distance to protein",
            "function": "intermolecular_distance",
            "parameters": {
                "radius_type": "vdw",
                "radius_scale": 1.0,
                "clash_cutoff": 0.75,
                "ignore_types": ["hydrogens", "organic_cofactors", "inorganic_cofactors", "waters"],
                "max_distance": 5.0,
            },
            "chosen_binary_test_output": ["not_too_far_away", "no_clashes"],
            "rename_outputs": {
                "not_too_far_away": "Protein-ligand maximum distance",
                "no_clashes": "Minimum distance to protein",
            },
            "rename_suffix": "_protein",
        },
        *(
            {
                "name": f"Distance to {label}",
                "function": "intermolecular_distance",
                "parameters": {
                    "radius_type": radius,
                    "radius_scale": 1.0,
                    "clash_cutoff": 0.75,
                    "ignore_types": ["hydrogens", *sorted(PARTNER_TYPES - {kind})],
                    "max_distance": 5.0,
                },
                "chosen_binary_test_output": ["no_clashes"],
                "rename_outputs": {"no_clashes": f"Minimum distance to {label}"},
                "rename_suffix": f"_{kind}",
            }
            for kind, label, radius in (
                ("organic_cofactors", "organic cofactors", "vdw"),
                ("inorganic_cofactors", "inorganic cofactors", "covalent"),
                ("waters", "waters", "vdw"),
            )
        ),
    ],
    "loading": {
        "mol_pred": LOAD,
        "mol_true": LOAD,
        "mol_cond": {**LOAD, "proximityBonding": False},
    },
}


def fetch(url: str, dest: Path) -> Path:
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=60) as response:
            data = response.read()
        if url.endswith(".gz"):
            import gzip

            data = gzip.decompress(data)
        dest.write_bytes(data)
    return dest


def protein_ligand_ids(manifest: Path) -> list[str]:
    strata = tomllib.loads(manifest.read_text())["strata"]
    ids = []
    for name, entries in strata.items():
        if "ligand" in name:
            ids.extend(entries)
    return ids


def dictionary(path: Path) -> dict:
    """The heavy atoms of a dictionary entry in order, with whether each is a
    leaving atom, and its E or Z bonds by atom names."""
    block = gemmi.cif.read(str(path)).sole_block()
    atoms = block.find("_chem_comp_atom.", ["atom_id", "type_symbol", "pdbx_leaving_atom_flag"])
    heavy = [(row.str(0), row.str(1).upper(), row.str(2) == "Y") for row in atoms if row.str(1).upper() not in ("H", "D")]
    bonds = block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "pdbx_stereo_config"])
    stereo = [(row.str(0), row.str(1)) for row in bonds if row.str(2) in ("E", "Z")]
    missing = block.find_value("_chem_comp.pdbx_ideal_coordinates_missing_flag")
    return {"heavy": heavy, "stereo_bonds": stereo, "ideal_missing": missing == "Y"}


def first_complete_copy(structure: gemmi.Structure, code: str, names: list[str]):
    for chain in structure[0]:
        for residue in chain:
            if residue.name != code:
                continue
            atoms = {atom.name: atom for atom in residue if atom.element.name not in ("H", "D")}
            if set(atoms) == set(names):
                return chain.name, residue.seqid.num, atoms
    return None


def graph(mol: Chem.Mol) -> list[list[int]]:
    return [[n.GetIdx() for n in atom.GetNeighbors()] for atom in mol.GetAtoms()]


def side(neighbours: list[list[int]], start: int, block: int) -> set[int] | None:
    """The atoms reachable from `start` without passing `block`, or None when
    `block` is among them (the bond is in a ring)."""
    seen = {start}
    queue = deque([start])
    while queue:
        atom = queue.popleft()
        for other in neighbours[atom]:
            if atom == start and other == block:
                continue
            if other == block:
                return None
            if other not in seen:
                seen.add(other)
                queue.append(other)
    return seen


def heavy_count(mol: Chem.Mol, atoms) -> int:
    return sum(1 for a in atoms if mol.GetAtomWithIdx(a).GetAtomicNum() > 1)


def with_positions(mol: Chem.Mol, positions: np.ndarray, name: str) -> Chem.Mol:
    copy = Chem.Mol(mol)
    conformer = copy.GetConformer()
    for at, position in enumerate(positions):
        conformer.SetAtomPosition(at, position.tolist())
    copy.SetProp("_Name", name)
    return copy


def rotated(points: np.ndarray, origin: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    cos, sin = math.cos(angle), math.sin(angle)
    relative = points - origin
    return origin + relative * cos + np.cross(axis, relative) * sin + np.outer(relative @ axis, axis) * (1 - cos)


def terminal_bond(mol: Chem.Mol, neighbours):
    """A non-ring single bond between heavy atoms whose far side holds the
    fewest heavy atoms, with a heavy neighbour of the near end besides."""
    best = None
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bond.IsInRing() or bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        if mol.GetAtomWithIdx(a).GetAtomicNum() == 1 or mol.GetAtomWithIdx(b).GetAtomicNum() == 1:
            continue
        for near, far in ((a, b), (b, a)):
            group = side(neighbours, far, near)
            others = [n for n in neighbours[near] if n != far and mol.GetAtomWithIdx(n).GetAtomicNum() > 1]
            if group is None or not others:
                continue
            key = (heavy_count(mol, group), near, far)
            if best is None or key < best[0]:
                best = (key, near, far, group, others[0])
    return best


def rotatable_bond(mol: Chem.Mol, neighbours):
    """A non-ring single bond between heavy atoms each with another heavy
    neighbour, the far side the smaller."""
    best = None
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bond.IsInRing() or bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        heavy = lambda at, other: [n for n in neighbours[at] if n != other and mol.GetAtomWithIdx(n).GetAtomicNum() > 1]
        if not heavy(a, b) or not heavy(b, a):
            continue
        for near, far in ((a, b), (b, a)):
            group = side(neighbours, far, near)
            if group is None:
                continue
            size = heavy_count(mol, group)
            rest = mol.GetNumHeavyAtoms() - size
            if size < 2 or rest < 2:
                continue
            key = (abs(size - rest), near, far)
            if best is None or key < best[0]:
                best = (key, near, far, group)
    return best


def closest_contact(mol: Chem.Mol, positions: np.ndarray, distances: np.ndarray) -> float:
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    best = math.inf
    for i in heavy:
        for j in heavy:
            if j > i and distances[i][j] >= 3:
                best = min(best, float(np.linalg.norm(positions[i] - positions[j])))
    return best


def variants(mol: Chem.Mol, ideal: Chem.Mol, stereo_bonds, names: list[str], protein: np.ndarray, seed: int):
    """The crystal pose and its perturbed copies, each named."""
    positions = mol.GetConformer().GetPositions()
    neighbours = graph(mol)
    distances = Chem.GetDistanceMatrix(mol)
    out = [with_positions(mol, positions, "crystal")]

    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    aligned = Chem.Mol(ideal)
    rdMolAlign.AlignMol(aligned, mol, atomMap=[(i, i) for i in heavy])
    out.append(with_positions(mol, aligned.GetConformer().GetPositions(), "ideal"))

    centre = positions[heavy].mean(axis=0)
    mirrored = positions.copy()
    mirrored[:, 0] = 2 * centre[0] - mirrored[:, 0]
    out.append(with_positions(mol, mirrored, "mirrored"))

    terminal = terminal_bond(mol, neighbours)
    if terminal:
        _, near, far, group, other = terminal
        group = sorted(group)
        direction = positions[far] - positions[near]
        direction /= np.linalg.norm(direction)
        for amount in (0.25, 0.4, 0.6):
            moved = positions.copy()
            moved[group] += direction * amount
            out.append(with_positions(mol, moved, f"stretched-{amount}"))
        axis = np.cross(positions[other] - positions[near], positions[far] - positions[near])
        if np.linalg.norm(axis) > 1e-3:
            for degrees in (20, 40, 60):
                moved = positions.copy()
                moved[group] = rotated(positions[group], positions[near], axis, math.radians(degrees))
                out.append(with_positions(mol, moved, f"bent-{degrees}"))

    rotatable = rotatable_bond(mol, neighbours)
    if rotatable:
        _, near, far, group = rotatable
        group = sorted(group)
        axis = positions[far] - positions[near]
        trials = []
        for degrees in range(10, 360, 10):
            moved = positions.copy()
            moved[group] = rotated(positions[group], positions[far], axis, math.radians(degrees))
            trials.append((closest_contact(mol, moved, distances), degrees, moved))
        trials.sort(key=lambda trial: trial[0])
        out.append(with_positions(mol, trials[0][2], "twisted"))
        out.append(with_positions(mol, trials[len(trials) // 2][2], "turned"))

    ring = next((r for r in mol.GetRingInfo().AtomRings() if len(r) in (5, 6) and all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in r)), None)
    if ring:
        points = positions[list(ring)]
        _, _, vt = np.linalg.svd(points - points.mean(axis=0))
        normal = vt[-1]
        atom = ring[0]
        attached = [atom] + [n for n in neighbours[atom] if mol.GetAtomWithIdx(n).GetAtomicNum() == 1]
        for amount in (0.2, 0.35, 0.6):
            moved = positions.copy()
            moved[attached] += normal * amount
            out.append(with_positions(mol, moved, f"puckered-{amount}"))

    saturated = next(
        (
            r for r in mol.GetRingInfo().AtomRings()
            if len(r) == 6 and not any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in r)
        ),
        None,
    )
    if saturated:
        ring = list(saturated)
        points = positions[ring]
        middle = points.mean(axis=0)
        _, _, vt = np.linalg.svd(points - middle)
        normal = vt[-1]
        for share in (0.5, 1.0):
            moved = positions.copy()
            offsets = (points - middle) @ normal
            moved[ring] = points - share * np.outer(offsets, normal)
            out.append(with_positions(mol, moved, f"flattened-{share}"))

    index = {name: at for at, name in enumerate(names)}
    for first, second in stereo_bonds:
        a, b = index[first], index[second]
        group = side(neighbours, a, b)
        if group is None:
            continue
        group = sorted(group - {a})
        moved = positions.copy()
        moved[group] = rotated(positions[group], positions[a], positions[b] - positions[a], math.pi)
        out.append(with_positions(mol, moved, "flipped"))
        break

    nearest = protein[np.argmin(np.linalg.norm(protein - centre, axis=1))]
    toward = nearest - centre
    toward /= np.linalg.norm(toward)
    for amount in (0.5, 1.0, 1.5, 2.0):
        out.append(with_positions(mol, positions + toward * amount, f"shifted-{amount}"))

    rng = random.Random(seed)
    for number in range(1, 4):
        axis = np.array([rng.gauss(0, 1) for _ in range(3)])
        moved = rotated(positions, centre, axis, math.radians(rng.uniform(-20, 20)))
        moved = moved + np.array([rng.uniform(-1, 1) for _ in range(3)])
        if rotatable:
            _, near, far, group = rotatable
            group = sorted(group)
            moved[group] = rotated(moved[group], moved[far], moved[far] - moved[near], math.radians(rng.uniform(0, 360)))
        out.append(with_positions(mol, moved, f"random-{number}"))

    return out


def partner_atoms(receptor: gemmi.Structure) -> dict[str, np.ndarray]:
    """The receptor's heavy atom positions by the kind PoseBusters takes each
    for: protein, ion, water or cofactor."""
    kinds: dict[str, list] = {"protein": [], "ion": [], "water": [], "cofactor": []}
    for chain in receptor[0]:
        for residue in chain:
            for atom in residue:
                if residue.name == "HOH":
                    kind = "water"
                elif atom.element.name.upper() in ION_ELEMENTS or residue.name in ION_RESIDUES:
                    kind = "ion"
                elif residue.het_flag == "H":
                    kind = "cofactor"
                else:
                    kind = "protein"
                kinds[kind].append([atom.pos.x, atom.pos.y, atom.pos.z])
    return {kind: np.array(points).reshape(-1, 3) for kind, points in kinds.items()}


def nearest(points: np.ndarray, others: np.ndarray) -> tuple[float, int, int]:
    """The shortest distance between the two sets, and the pair that makes it."""
    gaps = np.linalg.norm(points[:, None, :] - others[None, :, :], axis=-1)
    at, other = np.unravel_index(gaps.argmin(), gaps.shape)
    return float(gaps[at, other]), int(at), int(other)


def partner_variants(mol: Chem.Mol, partners: dict[str, np.ndarray], seed: int):
    """The crystal pose moved towards the nearest ion, water and cofactor atom,
    out from the protein, and at random, each named."""
    positions = mol.GetConformer().GetPositions()
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    out = [with_positions(mol, positions, "crystal")]

    for kind in ("ion", "water", "cofactor"):
        points = partners[kind]
        if len(points) == 0:
            continue
        gap, at, other = nearest(positions[heavy], points)
        if gap > 6.0:
            continue
        toward = (points[other] - positions[heavy][at]) / gap
        for amount in (0.25, 0.5, 0.75, 1.0, 1.5):
            out.append(with_positions(mol, positions + toward * amount, f"toward-{kind}-{amount}"))

    protein = partners["protein"]
    centre = positions[heavy].mean(axis=0)
    outward = centre - protein.mean(axis=0)
    outward /= np.linalg.norm(outward)
    step = 0.0
    for mark in (4.5, 5.5):
        while nearest(positions[heavy] + outward * step, protein)[0] <= mark and step < 100:
            step += 0.25
        out.append(with_positions(mol, positions + outward * step, f"away-{mark}"))
    out.append(with_positions(mol, positions + outward * (step + 10), "away-far"))

    rng = random.Random(seed)
    for number in range(1, 4):
        axis = np.array([rng.gauss(0, 1) for _ in range(3)])
        moved = rotated(positions, centre, axis, math.radians(rng.uniform(-20, 20)))
        moved = moved + np.array([rng.uniform(-1, 1) for _ in range(3)])
        out.append(with_positions(mol, moved, f"random-{number}"))

    return out


def build_case(entry: str, code: str, args, ccd_dir: Path):
    cif = fetch(CCD_CIF.format(code=code), ccd_dir / f"{code}.cif")
    entry_info = dictionary(cif)
    if entry_info["ideal_missing"]:
        return None
    names = [name for name, _, _ in entry_info["heavy"]]
    sdf = fetch(CCD_SDF.format(code=code), args.cases / "ideal" / f"{code}.sdf")
    ideal = Chem.MolFromMolFile(str(sdf), removeHs=False)
    if ideal is None:
        return None
    template = Chem.RemoveHs(ideal)
    elements = [a.GetSymbol().upper() for a in template.GetAtoms()]
    if elements != [element for _, element, _ in entry_info["heavy"]]:
        return None

    path = fetch(RCSB_CIF.format(id=entry), args.cache / f"{entry}.cif")
    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    structure.remove_alternative_conformations()
    copy = first_complete_copy(structure, code, names)
    if copy is None:
        return None
    _, _, atoms = copy
    crystal = Chem.Mol(template)
    conformer = crystal.GetConformer()
    for at, name in enumerate(names):
        position = atoms[name].pos
        conformer.SetAtomPosition(at, (position.x, position.y, position.z))
    crystal = Chem.AddHs(crystal, addCoords=True)
    # The ideal conformer with its hydrogens added the same way, so its atoms
    # lie in the crystal pose's order.
    ideal_pose = Chem.AddHs(template, addCoords=True)

    receptor = structure.clone()
    if args.partners:
        for chain in receptor[0]:
            for at in reversed(range(len(chain))):
                if chain[at].name == code:
                    del chain[at]
    else:
        receptor.remove_ligands_and_waters()
    receptor.remove_hydrogens()
    receptor.remove_empty_chains()
    protein = np.array([[a.pos.x, a.pos.y, a.pos.z] for chain in receptor[0] for residue in chain for a in residue])

    case = args.cases / f"{entry}_{code}"
    case.mkdir(parents=True, exist_ok=True)
    receptor.write_pdb(str(case / "receptor.pdb"))
    seed = zlib.crc32(f"{entry}_{code}".encode())
    if args.partners:
        poses = partner_variants(crystal, partner_atoms(receptor), seed)
    else:
        poses = variants(crystal, ideal_pose, entry_info["stereo_bonds"], names, protein, seed)
    writer = Chem.SDWriter(str(case / "poses.sdf"))
    writer.SetKekulize(True)
    for pose in poses:
        writer.write(pose)
    writer.close()

    return case, sdf, [pose.GetProp("_Name") for pose in poses]


def bust(case: Path, ideal_sdf: Path, names: list[str]) -> dict:
    table = PoseBusters(config=CONFIG).bust(
        mol_pred=str(case / "poses.sdf"),
        mol_true=str(ideal_sdf),
        mol_cond=str(case / "receptor.pdb"),
    )
    verdicts = {}
    for (_, name, _), row in table.iterrows():
        verdicts[name] = {check: _passed(row.get(check)) for check in CHECKS}
    missing = set(names) - set(verdicts)
    if missing:
        raise RuntimeError(f"PoseBusters reported no row for {sorted(missing)}")
    return verdicts


def _passed(value) -> bool | None:
    """A check's outcome, or None when PoseBusters could not make it."""
    if value is None or value != value:
        return None
    return bool(value)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", type=Path, default=Path("benchmarks/manifest.toml"))
    parser.add_argument("--cache", type=Path, default=Path("benchmarks/cache"))
    parser.add_argument("--cases", type=Path)
    parser.add_argument("--only", nargs="*", help="entry ids to (re)generate")
    parser.add_argument("--partners", action="store_true", help="the PARTNERS entries, ions, cofactors and waters kept")
    args = parser.parse_args()
    if args.cases is None:
        args.cases = Path("benchmarks/posebusters-partners" if args.partners else "benchmarks/posebusters")

    ccd_dir = args.cases / "ccd"
    verdicts_path = args.cases / "verdicts.json"
    verdicts = json.loads(verdicts_path.read_text()) if verdicts_path.exists() else {}

    def record(built) -> None:
        case, ideal_sdf, names = built
        verdicts[case.name] = bust(case, ideal_sdf, names)
        print(f"{case.name}: {len(names)} poses", file=sys.stderr)
        verdicts_path.write_text(json.dumps(verdicts, indent=1, sort_keys=True) + "\n")

    if args.partners:
        for entry, code in PARTNERS:
            if args.only and entry not in args.only:
                continue
            built = build_case(entry, code, args, ccd_dir)
            if built is None:
                print(f"{entry} {code}: skipped (no complete copy, or no ideal model)", file=sys.stderr)
                continue
            record(built)
        return 0

    for entry in protein_ligand_ids(args.manifest):
        if args.only and entry not in args.only:
            continue
        structure = gemmi.read_structure(str(fetch(RCSB_CIF.format(id=entry), args.cache / f"{entry}.cif")))
        structure.setup_entities()
        codes = sorted(
            {
                r.name
                for chain in structure[0]
                for r in chain
                if r.het_flag == "H" and r.entity_type == gemmi.EntityType.NonPolymer
            }
            - ADDITIVES
        )
        for code in codes:
            residue = next(r for chain in structure[0] for r in chain if r.name == code)
            elements = {a.element.name.upper() for a in residue}
            if elements & METALS or sum(1 for a in residue if a.element.name not in ("H", "D")) < 6:
                continue
            try:
                built = build_case(entry, code, args, ccd_dir)
            except Exception as error:  # noqa: BLE001
                print(f"{entry} {code}: skipped ({error})", file=sys.stderr)
                continue
            if built is None:
                print(f"{entry} {code}: skipped (no complete copy, or no ideal model)", file=sys.stderr)
                continue
            record(built)

    return 0


if __name__ == "__main__":
    sys.exit(main())
