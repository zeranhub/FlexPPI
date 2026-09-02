"""Recompute a small ThermoPCD apo/holo pilot with consistent chain fitting.

Three local pairs have an unambiguous production XTC on both sides after
checking the downloaded directory contents.  Each state is sampled to at most
5,000 frames, divided into five contiguous time blocks, and the target chain is
least-squares fitted independently inside every block before C-alpha RMSF is
computed.  This matches FlexPPI's label convention and avoids treating rigid
motion of an antigen relative to an antibody as internal flexibility.
"""

from __future__ import annotations

import argparse
import ctypes
import json
import math
import re
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from flexppi.data.labels import _chainwise_rmsf


PILOT_PAIRS = (
    # pdb, temperature, isolated-component seed, complex seed
    ("5OCX", 310.0, 1, 1),
    # The nominal seed-2 holo XTC ends at 37.56 ns; use two complete 100 ns
    # independent replicates instead.
    ("6BFS", 311.0, 1, 3),
    ("6SF6", 312.0, 2, 2),
    # Seed labels are local replicate identifiers, not coupled randomness
    # across separately initialized apo/holo ensembles.  These two systems
    # therefore use different clean local replicates at the same temperature.
    ("2XWT", 310.0, 1, 4),
    ("4G3Y", 310.0, 1, 4),
)
PDB_RE_TEMPLATE = r"(?i)(?<![A-Z0-9]){pdb}(?![A-Z0-9])"
TEMP_RE = re.compile(r"(?i)(298|309|310|311(?:\.5)?|312|313)\s*K\b")
SEED_FILE_RE = re.compile(r"(?i)\brs[-_ ]?(\d+)")


@dataclass
class State:
    pdb_id: str
    kind: str
    temperature_k: float
    seed: int
    directory: Path
    topology: Path
    trajectory: Path


def classify(text: str) -> str:
    lower = text.lower().replace("_", " ")
    return "holo" if re.search(r"\bcomplex\b|\bbound\b", lower) else "apo"


def windows_short_path(path: Path) -> str:
    """Return the ASCII 8.3 path needed by the XDR reader on Windows."""

    buffer = ctypes.create_unicode_buffer(32768)
    length = ctypes.windll.kernel32.GetShortPathNameW(str(path), buffer, len(buffer))
    if length == 0 or length >= len(buffer):
        return str(path)
    return buffer.value


def discover_state(root: Path, pdb_id: str, temperature: float, seed: int, kind: str) -> State:
    candidates: list[tuple[int, Path, Path, Path]] = []
    pdb_pattern = re.compile(PDB_RE_TEMPLATE.format(pdb=re.escape(pdb_id)))
    for directory in (root / "data").iterdir():
        if not directory.is_dir() or not pdb_pattern.search(directory.name):
            continue
        if classify(directory.name) != kind:
            continue
        temp_match = TEMP_RE.search(directory.name)
        if temp_match is None or float(temp_match.group(1)) != temperature:
            continue
        trajectories = [
            path
            for path in directory.iterdir()
            if path.is_file() and path.name.lower() in {"dynamic-nopbc.xtc", "system_reduced.xtc"}
        ]
        if not trajectories:
            continue
        directory_seed_match = re.search(r"(?i)random\s*seed\s*#\s*(\d+)", directory.name)
        directory_seed = int(directory_seed_match.group(1)) if directory_seed_match else None
        snapshots = [
            path
            for path in directory.iterdir()
            if path.is_file()
            and path.suffix.lower() == ".pdb"
            and f"{temperature:g}k" in path.name.lower()
            and "0ns" in path.name.lower()
            and (f"rs{seed}" in path.name.lower() or directory_seed == seed)
        ]
        if not snapshots:
            continue
        observed_seeds = {
            int(match.group(1))
            for path in directory.iterdir()
            if path.is_file()
            for match in [SEED_FILE_RE.search(path.name)]
            if match is not None
        }
        observed_temperatures = {
            float(match.group(1))
            for path in directory.iterdir()
            if path.is_file()
            for match in [TEMP_RE.search(path.name)]
            if match is not None
        }
        # A generic XTC inside a directory formed by merging two Dataverse
        # entries cannot be assigned to either seed or temperature safely.
        if (observed_seeds and observed_seeds != {seed}) or (
            not observed_seeds and directory_seed != seed
        ) or observed_temperatures != {temperature}:
            continue
        directory_seed_bonus = int(directory_seed == seed)
        score = directory_seed_bonus * 10 + int(trajectories[0].name.lower() == "dynamic-nopbc.xtc")
        candidates.append((score, directory, snapshots[0], trajectories[0]))
    if not candidates:
        raise FileNotFoundError(
            f"No unambiguous {kind} state for {pdb_id} {temperature:g} K seed {seed}"
        )
    _, directory, topology, trajectory = max(candidates, key=lambda item: item[0])
    return State(pdb_id, kind, temperature, seed, directory, topology, trajectory)


def ca_chains(universe: Any) -> dict[str, list[Any]]:
    chains: dict[str, list[Any]] = {}
    for residue in universe.residues:
        ca_atoms = [atom for atom in residue.atoms if str(atom.name).upper() == "CA"]
        if not ca_atoms:
            continue
        atom = ca_atoms[0]
        chain = str(getattr(atom, "chainID", "")).strip() or str(residue.segid).strip() or "_"
        chains.setdefault(chain, []).append(residue)
    return chains


def residue_names(residues: list[Any]) -> list[str]:
    return [str(residue.resname).upper() for residue in residues]


def match_target_chain(apo_universe: Any, holo_universe: Any) -> tuple[list[Any], list[Any], dict[str, Any]]:
    apo_chains = ca_chains(apo_universe)
    holo_chains = ca_chains(holo_universe)
    if not apo_chains or not holo_chains:
        raise ValueError("No C-alpha-containing chain found")
    # The curated pilot isolates one antigen/epitope chain in every apo state.
    apo_chain, apo_residues = max(apo_chains.items(), key=lambda item: len(item[1]))
    apo_names = residue_names(apo_residues)
    scored: list[tuple[float, int, str, list[Any]]] = []
    for chain, residues in holo_chains.items():
        names = residue_names(residues)
        identity = SequenceMatcher(a=apo_names, b=names, autojunk=False).ratio()
        scored.append((identity, -abs(len(names) - len(apo_names)), chain, residues))
    identity, _, holo_chain, holo_residues = max(scored, key=lambda item: (item[0], item[1]))
    if identity < 0.95 or len(apo_residues) != len(holo_residues):
        raise ValueError(
            f"Target-chain match failed: apo {apo_chain}/{len(apo_residues)}, "
            f"holo {holo_chain}/{len(holo_residues)}, identity={identity:.3f}"
        )
    return apo_residues, holo_residues, {
        "apo_chain": apo_chain,
        "holo_chain": holo_chain,
        "residue_count": len(apo_residues),
        "sequence_identity": identity,
        "apo_chain_count": len(apo_chains),
        "holo_chain_count": len(holo_chains),
    }


def ca_indices(residues: list[Any]) -> np.ndarray:
    return np.asarray(
        [next(int(atom.index) for atom in residue.atoms if str(atom.name).upper() == "CA") for residue in residues],
        dtype=np.int64,
    )


def load_ca_frames(
    universe: Any, indices: np.ndarray, max_frames: int, common_end_time_ps: float
) -> tuple[np.ndarray, dict[str, Any]]:
    total = len(universe.trajectory)
    dt = float(universe.trajectory.dt)
    eligible = min(total, int(math.floor(common_end_time_ps / dt + 1e-6)) + 1)
    stride = max(1, int(math.ceil(eligible / max_frames)))
    frames = np.stack(
        [
            np.asarray(ts.positions[indices], dtype=np.float32).copy()
            for ts in universe.trajectory[:eligible:stride]
        ],
        axis=0,
    )
    return frames, {
        "trajectory_frames": total,
        "trajectory_end_time_ps": float(universe.trajectory[-1].time),
        "eligible_frames_to_common_time": eligible,
        "frame_stride": stride,
        "sampled_frames": len(frames),
    }


def block_rmsf(frames: np.ndarray, block_count: int) -> np.ndarray:
    if len(frames) < block_count * 10:
        raise ValueError(f"Too few sampled frames: {len(frames)}")
    edges = np.linspace(0, len(frames), block_count + 1, dtype=int)
    chain_index = np.zeros(frames.shape[1], dtype=np.int64)
    return np.stack(
        [_chainwise_rmsf(frames[start:stop], chain_index) for start, stop in zip(edges[:-1], edges[1:])],
        axis=0,
    )


def interface_distance(holo_universe: Any, target_indices: np.ndarray) -> np.ndarray:
    holo_universe.trajectory[0]
    target = np.asarray(holo_universe.atoms[target_indices].positions, dtype=np.float64)
    target_set = set(int(value) for value in target_indices)
    partner_indices = [
        int(atom.index)
        for residues in ca_chains(holo_universe).values()
        for residue in residues
        for atom in residue.atoms
        if str(atom.name).upper() == "CA" and int(atom.index) not in target_set
    ]
    partner = np.asarray(holo_universe.atoms[partner_indices].positions, dtype=np.float64)
    return np.sqrt(((target[:, None, :] - partner[None, :, :]) ** 2).sum(axis=2)).min(axis=1)


def compare_source_labels(
    residue_rows: pd.DataFrame,
    source_path: Path,
    pdb_id: str,
    temperature: float,
    apo_seed: int,
    holo_seed: int,
) -> dict[str, Any]:
    if apo_seed != holo_seed:
        return {
            "source_label_comparison_available": False,
            "source_label_comparison_reason": "different independent replicate labels",
        }
    if not source_path.exists():
        return {}
    source = pd.read_csv(source_path)
    selected = source.loc[
        (source["pdb_id"] == pdb_id)
        & (source["temperature_k"] == temperature)
        & (source["seed"] == apo_seed)
    ].sort_values("residue_index")
    current = residue_rows.sort_values("target_index")
    if len(selected) != len(current):
        return {"source_label_residue_count": int(len(selected)), "source_label_comparison_available": False}
    x = selected["delta_rmsf_angstrom"].to_numpy(float)
    y = current["delta_rmsf_angstrom"].to_numpy(float)
    correlation = spearmanr(x, y).statistic if len(x) >= 3 else math.nan
    return {
        "source_label_residue_count": int(len(x)),
        "source_label_comparison_available": True,
        "source_vs_refit_delta_spearman": float(correlation),
        "source_vs_refit_delta_mae_angstrom": float(np.mean(np.abs(x - y))),
    }


def analyse_pair(
    root: Path,
    pdb_id: str,
    temperature: float,
    apo_seed: int,
    holo_seed: int,
    max_frames: int,
):
    import MDAnalysis as mda

    apo = discover_state(root, pdb_id, temperature, apo_seed, "apo")
    holo = discover_state(root, pdb_id, temperature, holo_seed, "holo")
    apo_universe = mda.Universe(windows_short_path(apo.topology), windows_short_path(apo.trajectory))
    holo_universe = mda.Universe(windows_short_path(holo.topology), windows_short_path(holo.trajectory))
    apo_residues, holo_residues, chain_qc = match_target_chain(apo_universe, holo_universe)
    apo_ca = ca_indices(apo_residues)
    holo_ca = ca_indices(holo_residues)
    apo_end_time = float(apo_universe.trajectory[-1].time)
    holo_end_time = float(holo_universe.trajectory[-1].time)
    common_end_time = min(apo_end_time, holo_end_time)
    apo_frames, apo_qc = load_ca_frames(apo_universe, apo_ca, max_frames, common_end_time)
    holo_frames, holo_qc = load_ca_frames(holo_universe, holo_ca, max_frames, common_end_time)
    apo_blocks = block_rmsf(apo_frames, 5)
    holo_blocks = block_rmsf(holo_frames, 5)
    delta_blocks = holo_blocks - apo_blocks
    delta_mean = delta_blocks.mean(axis=0)
    delta_sd = delta_blocks.std(axis=0, ddof=1)
    half_width = 2.7764451051977987 * delta_sd / math.sqrt(5)
    distances = interface_distance(holo_universe, holo_ca)
    residue_rows = pd.DataFrame(
        {
            "pdb_id": pdb_id,
            "temperature_k": temperature,
            "apo_seed": apo_seed,
            "holo_seed": holo_seed,
            "same_seed_label": apo_seed == holo_seed,
            "target_index": np.arange(len(delta_mean)),
            "residue_name": residue_names(holo_residues),
            "author_residue_id": [int(residue.resid) for residue in holo_residues],
            "static_distance_angstrom": distances,
            "interface_calpha8": distances <= 8.0,
            "apo_rmsf_angstrom": apo_blocks.mean(axis=0),
            "holo_rmsf_angstrom": holo_blocks.mean(axis=0),
            "delta_rmsf_angstrom": delta_mean,
            "delta_block_sd_angstrom": delta_sd,
            "ci95_low": delta_mean - half_width,
            "ci95_high": delta_mean + half_width,
            "stiffening": delta_mean < 0,
            "softening": delta_mean > 0,
            "robust_stiffening": delta_mean + half_width < 0,
            "robust_softening": delta_mean - half_width > 0,
        }
    )
    qc = {
        "pdb_id": pdb_id,
        "temperature_k": temperature,
        "apo_seed": apo_seed,
        "holo_seed": holo_seed,
        "same_seed_label": apo_seed == holo_seed,
        "apo_directory": str(apo.directory),
        "holo_directory": str(holo.directory),
        "apo_topology": str(apo.topology),
        "holo_topology": str(holo.topology),
        "apo_trajectory": str(apo.trajectory),
        "holo_trajectory": str(holo.trajectory),
        "common_end_time_ps": common_end_time,
        "common_duration_ns": common_end_time / 1000.0,
        "nominal_100ns_complete": common_end_time >= 99990.0,
        **chain_qc,
        **{f"apo_{key}": value for key, value in apo_qc.items()},
        **{f"holo_{key}": value for key, value in holo_qc.items()},
    }
    return residue_rows, qc


def region_summary(residue: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for pdb_id, system in residue.groupby("pdb_id", sort=False):
        for region, group in (
            ("all", system),
            ("interface_calpha8", system.loc[system["interface_calpha8"]]),
            ("noninterface_calpha8", system.loc[~system["interface_calpha8"]]),
        ):
            if group.empty:
                continue
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "region": region,
                    "residues": len(group),
                    "mean_delta_rmsf_angstrom": float(group["delta_rmsf_angstrom"].mean()),
                    "median_delta_rmsf_angstrom": float(group["delta_rmsf_angstrom"].median()),
                    "fraction_stiffening": float(group["stiffening"].mean()),
                    "fraction_softening": float(group["softening"].mean()),
                    "fraction_robust_stiffening": float(group["robust_stiffening"].mean()),
                    "fraction_robust_softening": float(group["robust_softening"].mean()),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(r"D:\data\ThermoPCD_trajectories"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/thermopcd_audit"))
    parser.add_argument("--max-frames", type=int, default=5000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    residue_parts: list[pd.DataFrame] = []
    qc_rows: list[dict[str, Any]] = []
    source_path = args.output_dir / "supplied_rmsf_pair_residues.csv"
    for pdb_id, temperature, apo_seed, holo_seed in PILOT_PAIRS:
        print(
            f"[{pdb_id}] {temperature:g} K apo-seed {apo_seed}, holo-seed {holo_seed}",
            flush=True,
        )
        residue, qc = analyse_pair(
            args.root, pdb_id, temperature, apo_seed, holo_seed, args.max_frames
        )
        qc.update(
            compare_source_labels(
                residue, source_path, pdb_id, temperature, apo_seed, holo_seed
            )
        )
        residue_parts.append(residue)
        qc_rows.append(qc)
    all_residue = pd.concat(residue_parts, ignore_index=True)
    summary = region_summary(all_residue)
    interface = summary.loc[summary["region"] == "interface_calpha8"]
    report = {
        "pilot_pairs": len(PILOT_PAIRS),
        "residues": int(len(all_residue)),
        "interface_residues_calpha8": int(all_residue["interface_calpha8"].sum()),
        "pdb_equal_interface_mean_delta_rmsf_angstrom": float(
            interface["mean_delta_rmsf_angstrom"].mean()
        ),
        "pdb_equal_interface_fraction_stiffening": float(interface["fraction_stiffening"].mean()),
        "pdb_equal_interface_fraction_robust_stiffening": float(
            interface["fraction_robust_stiffening"].mean()
        ),
        "pdb_equal_interface_fraction_robust_softening": float(
            interface["fraction_robust_softening"].mean()
        ),
        "systems_with_negative_interface_mean": int(
            (interface["mean_delta_rmsf_angstrom"] < 0).sum()
        ),
        "method": "five temporal blocks; target-chain-fitted C-alpha RMSF; t(4) 95% CI",
    }
    all_residue.to_csv(args.output_dir / "refit_pair_residue_response.csv", index=False)
    summary.to_csv(args.output_dir / "refit_pair_region_summary.csv", index=False)
    pd.DataFrame(qc_rows).to_csv(args.output_dir / "refit_pair_qc.csv", index=False)
    (args.output_dir / "refit_pair_summary.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
