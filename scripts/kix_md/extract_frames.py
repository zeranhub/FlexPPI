"""Extract evenly spaced production frames from the downloaded KIX trajectories.

Writes 10 PDB snapshots per state (5, 15, ..., 95 ns of the 100 ns production)
into ``artifacts/external/kix_md_20260914/frames/<state>/`` and a manifest with
frame indices, times and SHA256 values.  No model is involved.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from flexppi.data.structures import parse_pdb  # noqa: E402

LOCAL_ROOT = ROOT / "artifacts/external/kix_md_20260914"
STATES = ("apo", "myb", "mll", "both")
FRAME_TIMES_NS = tuple(range(5, 100, 10))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def extract_state(state: str) -> dict:
    directory = LOCAL_ROOT / state
    topology = directory / "protein.pdb"
    segments = sorted((directory / "production").glob("segment_*.dcd"))
    if len(segments) != 100:
        raise ValueError("{}: expected 100 segments, found {}".format(state, len(segments)))
    reference = parse_pdb(str(topology))
    universe = mda.Universe(str(topology), [str(segment) for segment in segments])
    n_frames = universe.trajectory.n_frames
    if n_frames != 100000:
        raise ValueError("{}: expected 100000 frames, found {}".format(state, n_frames))
    targets = {time_ns * 1000 - 1: time_ns for time_ns in FRAME_TIMES_NS}
    output_directory = LOCAL_ROOT / "frames" / state
    output_directory.mkdir(parents=True, exist_ok=True)
    rows = []
    for timestep in universe.trajectory:
        time_ns = targets.get(timestep.frame)
        if time_ns is None:
            continue
        path = output_directory / "{}_t{:03d}ns.pdb".format(state, time_ns)
        universe.atoms.write(str(path))
        written = parse_pdb(str(path))
        if written.num_residues != reference.num_residues:
            raise ValueError("{} t{}: residue count differs".format(state, time_ns))
        for chain in reference.unique_chains():
            if written.chain_sequence(chain) != reference.chain_sequence(chain):
                raise ValueError("{} t{}: chain {} sequence differs".format(state, time_ns, chain))
        rows.append({"state": state, "time_ns": time_ns, "frame_index": int(timestep.frame),
                     "n_atoms": int(universe.atoms.n_atoms),
                     "n_residues": int(reference.num_residues),
                     "chains": reference.unique_chains(),
                     "pdb": str(path.relative_to(ROOT)),
                     "sha256": sha256(path)})
        print("{} t={} ns frame={} -> {}".format(state, time_ns, timestep.frame, path.name),
              flush=True)
    return {"state": state, "topology": str(topology), "segments": len(segments),
            "frames": rows}


def main() -> None:
    manifest = {"frame_times_ns": list(FRAME_TIMES_NS), "states": {}}
    for state in STATES:
        manifest["states"][state] = extract_state(state)
    output = LOCAL_ROOT / "frames" / "frames_manifest.json"
    output.write_text(json.dumps(manifest, indent=2))
    print(json.dumps({"wrote": str(output),
                      "n_frames": sum(len(v["frames"]) for v in manifest["states"].values())},
                     indent=2))


if __name__ == "__main__":
    main()
