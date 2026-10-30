#!/usr/bin/env python3
"""Generate RING 4.0 interaction golden files for pixelfold-validate.

For every id in the benchmark manifest, this submits the mmCIF to the RING 4.0
web server (the public job API behind ring.biocomputingup.it), downloads the
edges file RING writes, and records the per-type interaction edges
pixelfold-validate compares against:

    <golden-dir>/<ID>.ring.json   {"kinds": [...], "interactions": [{kind, a, b}]}

RING's standalone binary is licensed to one person and may not be passed on, so
it cannot live in a pinned image the way DSSP, FreeSASA and PLIP do. The web
server needs no account, so it is the copy of RING anyone can run. It reports no
version, which is why each golden records the server task, the parameters sent,
and the date.

RING runs at its defaults over the whole structure, ligands included. Two of those
defaults narrow what it looks at, and the golden states both so the harness can
narrow pixelfold the same way: water is left out, and standard amino acids of one
chain fewer than three positions apart are never paired.

Each server reply is kept under --raw-dir, and a kept reply is parsed instead of
resubmitted, so a mapping change never costs the server another run. Jobs go one
at a time. Standard library only; needs Python 3.11 or newer (tomllib).
"""

from __future__ import annotations

import argparse
import datetime
import json
import sys
import time
import urllib.error
import urllib.request
import uuid
from pathlib import Path

import tomllib

RING_API = "https://drmaatic.biocomputingup.it"
RING_TASK = "ring-v4"
# Sent explicitly although both are RING's defaults, so a later change of server
# default cannot silently change what a regenerated golden means.
SEQUENCE_SEPARATION = 3
RING_PARAMS = {"seq_sep": str(SEQUENCE_SEPARATION), "policy": "closest"}

RCSB_CIF = "https://files.rcsb.org/download/{id}.cif"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

POLL_SECONDS = 3
# Short, because the server restarts its refill clock on the first look after a
# step, so the time between the step and that look is lost credit.
CREDIT_POLL_SECONDS = 5
PAUSE_SECONDS = 2
JOB_TIMEOUT_SECONDS = 1800
RETRIES = 5

JOB_SUCCEEDED = "job finished normally"
JOB_FAILED = (
    "job finished, but failed",
    "job has been rejected from the ws",
    "job has been deleted",
    "job was stopped by the user",
)

# RING interaction types mapped to pixelfold's type labels.
RING_KINDS: dict[str, str] = {
    "HBOND": "hydrogen-bond",
    "IONIC": "salt-bridge",
    "PIPISTACK": "pi-stacking",
    "PICATION": "pi-cation",
    "PIHBOND": "pi-hydrogen-bond",
    "HALOGEN": "halogen-bond",
    "METAL_ION": "metal-coordination",
    "SSBOND": "disulfide",
}

# RING types with no pixelfold counterpart, left out of the golden:
# - VDW is any two atoms whose van der Waals radii overlap; pixelfold's
#   hydrophobic contact is two apolar carbons within 4.0 A, a different criterion.
# - IAC, the generic inter-atomic contact, appears only under RING's CA, CB and
#   lollipop node policies, never under the closest-atom policy used here.
UNSCORED = frozenset({"VDW", "IAC"})


def structure_url(pdb_id: str) -> str:
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


def _open(request: urllib.request.Request | str) -> bytes:
    """Read a URL, retrying the failures a busy server recovers from."""
    last: Exception = RuntimeError("no attempt made")
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code < 500:
                raise
            last = error
        except urllib.error.URLError as error:
            last = error
        time.sleep(POLL_SECONDS * 2**attempt)
    raise last


def job_cost() -> int:
    """The credit one RING job draws from the server's anonymous budget."""
    return json.loads(_open(f"{RING_API}/task/{RING_TASK}/"))["required_tokens"]


def wait_for_credit(cost: int) -> None:
    """Sleep until the server's anonymous budget covers one job.

    Anonymous jobs draw on a credit bucket the server refills in whole steps. A
    refused job's Retry-After names only the next step, so this asks the server
    how long until the bucket holds enough. That estimate assumes a continuous
    refill and can reach zero a step early, so only the credit itself releases a
    job; the estimate just paces the next look.
    """
    while True:
        credit = json.loads(_open(f"{RING_API}/cpu-credit/?required={cost}"))
        if credit["wait_time"] == "infinite":
            raise RuntimeError(f"a RING job costs {cost} credits, more than the server ever grants")
        if credit["available_cpu_credit"] >= cost:
            return
        wait = max(credit["wait_time"], CREDIT_POLL_SECONDS)
        print(f"  waiting {wait} s for server credit", file=sys.stderr)
        time.sleep(wait + 1)


def _multipart(fields: dict[str, str], file_field: str, path: Path) -> tuple[bytes, str]:
    boundary = uuid.uuid4().hex
    parts = [
        f'--{boundary}\r\nContent-Disposition: form-data; name="{name}"\r\n\r\n{value}\r\n'.encode()
        for name, value in fields.items()
    ]
    parts.append(
        (
            f'--{boundary}\r\nContent-Disposition: form-data; name="{file_field}"; '
            f'filename="{path.name}"\r\nContent-Type: application/octet-stream\r\n\r\n'
        ).encode()
        + path.read_bytes()
        + b"\r\n"
    )
    parts.append(f"--{boundary}--\r\n".encode())
    return b"".join(parts), f"multipart/form-data; boundary={boundary}"


def submit(cif: Path, cost: int) -> str:
    """Submit one structure once the budget allows, returning the job id."""
    body, content_type = _multipart({"task": RING_TASK, **RING_PARAMS}, "input_file", cif)
    for _ in range(RETRIES):
        wait_for_credit(cost)
        request = urllib.request.Request(
            f"{RING_API}/job/", data=body, headers={"Content-Type": content_type}
        )
        try:
            return json.loads(_open(request))["uuid"]
        except urllib.error.HTTPError as error:
            # Another client on the same address can spend the credit first.
            if error.code != 429:
                raise
    raise RuntimeError(f"RING server still throttling after {RETRIES} credit waits")


def run_ring(cif: Path, cost: int) -> str:
    """Submit one structure to the RING web server and return its edges file."""
    job = submit(cif, cost)
    print(f"  job {job}", file=sys.stderr)

    deadline = time.monotonic() + JOB_TIMEOUT_SECONDS
    while True:
        status = json.loads(_open(f"{RING_API}/job/{job}/"))["status"]
        if status == JOB_SUCCEEDED:
            break
        if status in JOB_FAILED:
            raise RuntimeError(f"RING job {job}: {status}")
        if time.monotonic() > deadline:
            raise TimeoutError(f"RING job {job} still '{status}' after {JOB_TIMEOUT_SECONDS} s")
        time.sleep(POLL_SECONDS)

    files = json.loads(_open(f"{RING_API}/job/{job}/file/"))
    edges = [name for name in files if name.endswith("_ringEdges")]
    if not edges:
        raise RuntimeError(f"RING job {job} wrote no edges file: {files}")
    return _open(f"{RING_API}/job/{job}/file/{edges[0]}").decode()


def residue_key(node: str) -> dict:
    """A RING node id, `chain:seq:icode:name` with `_` for no insertion code."""
    chain, seq, icode, _name = node.split(":")
    return {"chain": chain, "seq": int(seq), "icode": None if icode == "_" else icode}


def edges_from(ring_edges: str) -> list[dict]:
    """Every scored RING interaction as `{kind, a, b}`, one per residue pair per type.

    RING writes one row per atom pair, so a residue pair can recur within a type;
    the harness compares residue pairs, so the repeats collapse here.
    """
    rows = [line.split("\t") for line in ring_edges.splitlines()[1:] if line.strip()]
    first_model = rows[0][11] if rows else None

    seen: set[tuple[str, str, str]] = set()
    edges: list[dict] = []
    for row in rows:
        node_a, interaction, node_b, model = row[0], row[1], row[2], row[11]
        if model != first_model:
            continue
        ring_type = interaction.split(":")[0]
        if ring_type in UNSCORED:
            continue
        if ring_type not in RING_KINDS:
            raise ValueError(f"unmapped RING interaction type {ring_type!r}")

        kind = RING_KINDS[ring_type]
        a, b = sorted(
            (residue_key(node_a), residue_key(node_b)),
            key=lambda r: (r["chain"], r["seq"], r["icode"] or ""),
        )
        if a == b:
            continue
        fingerprint = (kind, json.dumps(a, sort_keys=True), json.dumps(b, sort_keys=True))
        if fingerprint in seen:
            continue
        seen.add(fingerprint)
        edges.append({"kind": kind, "a": a, "b": b})

    return edges


def ring_golden(ring_edges: str, generated: str) -> dict:
    return {
        "source": {
            "tool": "RING 4.0 web server",
            "task": RING_TASK,
            "params": RING_PARAMS,
            "generated": generated,
        },
        "kinds": list(RING_KINDS.values()),
        "ligand_scoped": False,
        "water_excluded": True,
        "sequence_separation": SEQUENCE_SEPARATION,
        "interactions": edges_from(ring_edges),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("benchmarks/manifest.toml"))
    parser.add_argument("--golden-dir", type=Path, default=Path("benchmarks/golden"))
    parser.add_argument("--structures", type=Path, default=Path("benchmarks/cache"))
    parser.add_argument("--raw-dir", type=Path, default=Path("benchmarks/ring-cache"))
    parser.add_argument(
        "--force", action="store_true", help="resubmit even when a server reply is kept"
    )
    parser.add_argument("ids", nargs="*", help="PDB ids to process (default: whole manifest)")
    args = parser.parse_args()

    for directory in (args.golden_dir, args.structures, args.raw_dir):
        directory.mkdir(parents=True, exist_ok=True)

    ids = args.ids
    if not ids:
        manifest = tomllib.loads(args.manifest.read_text())
        ids = [pid for stratum in manifest["strata"].values() for pid in stratum]

    cost: int | None = None
    failures = 0
    for pdb_id in (pid.upper() for pid in ids):
        print(pdb_id, file=sys.stderr)
        raw = args.raw_dir / f"{pdb_id}_ringEdges"
        meta = args.raw_dir / f"{pdb_id}.generated"
        try:
            if args.force or not raw.exists():
                cost = cost if cost is not None else job_cost()
                cif = fetch_structure(pdb_id, args.structures)
                raw.write_text(run_ring(cif, cost))
                meta.write_text(datetime.date.today().isoformat())
                time.sleep(PAUSE_SECONDS)
            generated = meta.read_text().strip() if meta.exists() else "unknown"
            golden = ring_golden(raw.read_text(), generated)
            (args.golden_dir / f"{pdb_id}.ring.json").write_text(json.dumps(golden, indent=1))
            print(f"  {len(golden['interactions'])} interactions", file=sys.stderr)
        except urllib.error.HTTPError as error:
            if error.code in (401, 403):
                # A refusal (quota or permission) would refuse every later job too.
                print(f"  STOPPED: the RING server refused the job ({error})", file=sys.stderr)
                return 1
            print(f"  FAILED: {error}", file=sys.stderr)
            failures += 1
        except Exception as error:  # noqa: BLE001
            print(f"  FAILED: {error}", file=sys.stderr)
            failures += 1

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
