"""Extract per-residue side-chain chi entropy (H_chi) from ITrace trajectories.

Uses exactly the FlexPPI label convention (``CHI_ATOMS`` dihedral quadruples,
3-bin discretisation, base-3 rotamer coding, Shannon entropy over frames) and
the same residue keys as the RMSF benchmark (chain / residue_number /
insertion_code / amino_acid_3), so the result can be merged onto
``common_residue_predictions.csv`` for a per-residue external H_chi target.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from flexppi.constants import CHI_ATOMS, THREE_TO_ONE


FIELDNAMES = [
    "pdb_id",
    "traj_id",
    "chain",
    "residue_number",
    "insertion_code",
    "amino_acid_3",
    "chi_entropy",
    "n_frames",
    "n_chi_angles",
]


def dihedral_degrees(points: np.ndarray) -> np.ndarray:
    """Signed dihedral angle in degrees for an (n, 4, 3) coordinate array."""

    b0 = -(points[:, 1] - points[:, 0])
    b1 = points[:, 2] - points[:, 1]
    b2 = points[:, 3] - points[:, 2]
    b1 = b1 / np.maximum(np.linalg.norm(b1, axis=-1, keepdims=True), 1e-8)
    v = b0 - np.sum(b0 * b1, axis=-1, keepdims=True) * b1
    w = b2 - np.sum(b2 * b1, axis=-1, keepdims=True) * b1
    x = np.sum(v * w, axis=-1)
    y = np.sum(np.cross(b1, v) * w, axis=-1)
    return np.degrees(np.arctan2(y, x))


def entropy(codes: np.ndarray) -> float:
    _, counts = np.unique(codes, return_counts=True)
    probabilities = counts.astype(np.float64) / max(int(counts.sum()), 1)
    return float(-np.sum(probabilities * np.log(np.maximum(probabilities, 1e-12))))


def residue_key(atom) -> Tuple[str, str, str, str]:
    chain = str(getattr(atom, "chainID", "") or "").strip()
    if not chain:
        chain = str(getattr(atom, "segid", "") or "").strip() or "_"
    insertion_code = str(getattr(atom, "icode", "") or "").strip()
    return (chain, str(getattr(atom, "resid", "")), insertion_code, str(atom.resname).upper())


def process_trajectory(
    pdb_id: str,
    traj_id: str,
    topology_path: str,
    trajectory_path: str,
) -> Tuple[Optional[List[Dict[str, object]]], Optional[str]]:
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise RuntimeError("MDAnalysis is required for trajectory labels") from exc
    universe = mda.Universe(topology_path, trajectory_path)
    chi_quads: List[List[int]] = []
    residue_to_chis: Dict[int, List[int]] = {}
    keys: Dict[int, Tuple[str, str, str, str]] = {}
    for residue_index, residue in enumerate(universe.residues):
        resname = str(residue.resname).upper()
        if resname not in CHI_ATOMS or resname not in THREE_TO_ONE:
            continue
        calpha = [atom for atom in residue.atoms if str(atom.name).upper() == "CA"]
        if not calpha:
            continue
        lookup = {str(atom.name).upper(): int(atom.index) for atom in residue.atoms}
        definitions = CHI_ATOMS[resname]
        local: List[int] = []
        for definition in definitions:
            if all(name in lookup for name in definition):
                local.append(len(chi_quads))
                chi_quads.append([lookup[name] for name in definition])
        if local and len(local) == len(definitions):
            residue_to_chis[residue_index] = local
            keys[residue_index] = residue_key(calpha[0])
    if not chi_quads:
        return [], None
    quad_array = np.asarray(chi_quads, dtype=np.int64)
    state_frames: List[np.ndarray] = []
    for timestep in universe.trajectory:
        angles = dihedral_degrees(timestep.positions[quad_array])
        bins = np.clip(np.floor((angles + 180.0) / 120.0).astype(np.int8), 0, 2)
        state_frames.append(bins)
    states = np.stack(state_frames, axis=0)
    n_frames = int(states.shape[0])
    rows: List[Dict[str, object]] = []
    for residue_index, chi_indices in residue_to_chis.items():
        codes = np.zeros(n_frames, dtype=np.int32)
        for power, chi_index in enumerate(chi_indices):
            codes += states[:, chi_index].astype(np.int32) * (3 ** power)
        chain, residue_number, insertion_code, amino_acid_3 = keys[residue_index]
        rows.append(
            {
                "pdb_id": pdb_id,
                "traj_id": traj_id,
                "chain": chain,
                "residue_number": residue_number,
                "insertion_code": insertion_code,
                "amino_acid_3": amino_acid_3,
                "chi_entropy": entropy(codes),
                "n_frames": n_frames,
                "n_chi_angles": len(chi_indices),
            }
        )
    return rows, None


def worker(args) -> Tuple[Optional[List[Dict[str, object]]], Optional[str]]:
    pdb_id, traj_id, topology_path, trajectory_path = args
    try:
        return process_trajectory(pdb_id, traj_id, topology_path, trajectory_path)
    except Exception as error:  # noqa: BLE001 - record any trajectory failure
        return None, f"{traj_id}: {type(error).__name__}: {error}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--itrace-root", default=r"E:\data\ITrace")
    parser.add_argument(
        "--output-dir",
        default=r"E:\Project\FlexPPI\output\itrace_hchi_benchmark",
    )
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    root = Path(args.itrace_root)
    output_dir = Path(args.output_dir)
    raw_dir = output_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((root / "tables" / "manifest.json").read_text(encoding="utf-8"))
    tasks: List[Tuple[str, str, str, str]] = []
    for row in manifest["trajectories"]:
        traj_id = str(row["traj_id"])
        pdb_id = str(row["pdb_id"]).lower()
        topology = root / "trajectories" / f"{traj_id}_topology.pdb"
        trajectory = root / "trajectories" / f"{traj_id}_md_processed.xtc"
        if not topology.exists() or not trajectory.exists():
            print(f"MISSING FILES for {traj_id}", flush=True)
            continue
        tasks.append((pdb_id, traj_id, str(topology), str(trajectory)))

    output_path = raw_dir / "chi_entropy_per_replica.csv"
    failures: List[Dict[str, str]] = []
    started = time.time()
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        with Pool(processes=args.workers) as pool:
            for index, (rows, failure) in enumerate(pool.imap_unordered(worker, tasks), 1):
                if failure is not None:
                    failures.append({"traj_id": failure.split(":")[0], "error": failure})
                elif rows:
                    writer.writerows(rows)
                if index % 50 == 0 or index == len(tasks):
                    print(
                        f"H_chi extraction {index}/{len(tasks)} "
                        f"(failures={len(failures)}, elapsed={time.time() - started:.0f}s)",
                        flush=True,
                    )
                handle.flush()
    report = {
        "trajectory_count": len(tasks),
        "failure_count": len(failures),
        "failures": failures,
        "elapsed_seconds": time.time() - started,
        "convention": (
            "CHI_ATOMS quadruples; 3-bin rotamer discretisation; base-3 code; "
            "Shannon entropy over all readable frames per replica"
        ),
    }
    (output_dir / "extraction_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
