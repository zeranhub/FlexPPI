"""Evaluate FlexPPI and chain-wise FastProtFlex on the ITrace MD dataset.

The primary target is C-alpha RMSF after fitting every protein chain
independently, matching the label definition used to train FlexPPI.  The
official ITrace whole-complex-fitted RMSF profile is retained as a sensitivity
analysis.  Three MD replicas are collapsed before model-level statistics, so
the independent statistical unit is a PDB complex rather than a trajectory.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import MDAnalysis as mda
import numpy as np
import pandas as pd
import scipy
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, wilcoxon
import torch

from flexppi.comparison import _bootstrap_mean_ci, _system_metrics
from flexppi.data.dataset import collate_states
from flexppi.data.structures import ProteinStructure, parse_pdb
from flexppi.inference import _structure_example, load_model
from flexppi.trainer import move_to_device


METRICS = (
    "spearman",
    "pearson",
    "z_mae",
    "z_rmse",
    "top20_auroc",
    "top20_average_precision",
    "top20_precision_at_k",
)
ADVANTAGE_DIRECTION = {
    "spearman": 1.0,
    "pearson": 1.0,
    "z_mae": -1.0,
    "z_rmse": -1.0,
    "top20_auroc": 1.0,
    "top20_average_precision": 1.0,
    "top20_precision_at_k": 1.0,
}
LABEL_COLUMNS = {
    "chain_fitted_internal_rmsf": "target_chainfit_rmsf_angstrom",
    "itrace_official_complex_fitted_rmsf": "target_official_rmsf_angstrom",
}
ROLE_FIELDS = {
    "mhc_alpha_chain": "mhc",
    "b2m_chain": "b2m",
    "peptide_chain": "peptide",
    "tcr_alpha_chain": "tcr_alpha",
    "tcr_beta_chain": "tcr_beta",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--itrace-root", default=r"E:\data\ITrace")
    parser.add_argument("--checkpoint", default="artifacts/pretrain/best.pt")
    parser.add_argument("--output-dir", default="output/itrace_benchmark")
    parser.add_argument(
        "--rscript",
        default=r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe",
    )
    parser.add_argument(
        "--fastprotflex-script",
        default="scripts/fastprotflex_batch.R",
    )
    parser.add_argument("--bootstrap-repetitions", type=int, default=2000)
    parser.add_argument("--api-workers", type=int, default=12)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--reuse", action="store_true")
    parser.add_argument("--skip-api", action="store_true")
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_dump(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8")


def load_itrace_records(root: Path, limit: Optional[int]) -> List[Dict[str, Any]]:
    manifest = json.loads((root / "tables" / "manifest.json").read_text(encoding="utf-8"))
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for row in manifest["trajectories"]:
        grouped.setdefault(str(row["pdb_id"]).lower(), []).append(dict(row))
    records = []
    for pdb_id in sorted(grouped):
        replicas = sorted(grouped[pdb_id], key=lambda row: row["traj_id"])
        run1 = next(
            (row for row in replicas if str(row["traj_id"]).endswith("_run1")),
            replicas[0],
        )
        topology_path = root / "trajectories" / f"{run1['traj_id']}_topology.pdb"
        records.append(
            {
                "pdb_id": pdb_id,
                "topology_path": str(topology_path.resolve()),
                "replicas": replicas,
                "peptide_seq": run1.get("peptide_seq"),
                "title": run1.get("title"),
                "resolution": run1.get("resolution"),
                "host_species": run1.get("host_species"),
                "antigen_category": run1.get("antigen_category"),
            }
        )
    return records[:limit] if limit is not None else records


def training_pdb_ids(paths: Sequence[Path]) -> set[str]:
    identifiers: set[str] = set()
    for path in paths:
        if not path.exists():
            continue
        for line in path.read_text(encoding="utf-8").splitlines():
            record = json.loads(line)
            for pdb_id in record.get("metadata", {}).get("pdb_ids", []):
                identifiers.add(str(pdb_id).upper())
            sample_id = str(record.get("sample_id", "")).upper()
            if len(sample_id) >= 4 and sample_id[0].isdigit():
                identifiers.add(sample_id[:4])
    return identifiers


def fetch_json(url: str, retries: int = 3) -> Dict[str, Any]:
    last_error: Optional[Exception] = None
    for attempt in range(retries):
        try:
            request = Request(url, headers={"User-Agent": "FlexPPI-ITrace-Benchmark/1.0"})
            with urlopen(request, timeout=45) as response:
                return json.loads(response.read().decode("utf-8"))
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as error:
            last_error = error
            time.sleep(1.0 + attempt)
    raise RuntimeError(f"Could not fetch {url}: {last_error}")


def fetch_itrace_annotation(
    trajectory_id: str,
) -> Tuple[str, Dict[str, str], List[Dict[str, Any]]]:
    payload = fetch_json(f"https://itrace.life/api/trajectories/{trajectory_id}/analysis")
    identity = payload.get("identity", {}).get("complex_identity", {})
    roles = {}
    for field, role in ROLE_FIELDS.items():
        chain = identity.get(field)
        if chain is not None and str(chain).strip():
            roles[str(chain).strip()] = role
    rows = []
    for residue in payload.get("rmsf_profile", {}).get("residues", []):
        rows.append(
            {
                "traj_id": trajectory_id,
                "chain": str(residue.get("chain", "")).strip() or "_",
                "residue_number": str(residue.get("resid", "")),
                "amino_acid_3": str(residue.get("resname", "")).upper(),
                "region": str(residue.get("region", "")),
                "official_rmsf_angstrom": float(residue["rmsf"]),
            }
        )
    return trajectory_id, roles, rows


def infer_roles_from_structure(
    structure: ProteinStructure,
    peptide_sequence: str,
) -> Dict[str, str]:
    chains = structure.unique_chains()
    roles: Dict[str, str] = {}
    peptide_matches = [
        chain for chain in chains if structure.chain_sequence(chain) == peptide_sequence
    ]
    if peptide_matches:
        roles[peptide_matches[0]] = "peptide"
    remaining = [chain for chain in chains if chain not in roles]
    lengths = {chain: len(structure.chain_sequence(chain)) for chain in remaining}
    if remaining:
        b2m = min(remaining, key=lambda chain: abs(lengths[chain] - 100))
        roles[b2m] = "b2m"
    remaining = [chain for chain in remaining if chain not in roles]
    if remaining:
        mhc = max(remaining, key=lambda chain: lengths[chain])
        roles[mhc] = "mhc"
    remaining = [chain for chain in remaining if chain not in roles]
    for index, chain in enumerate(sorted(remaining)):
        roles[chain] = "tcr_alpha" if index == 0 else "tcr_beta"
    return roles


def collect_api_data(
    records: Sequence[Mapping[str, Any]],
    official_path: Path,
    roles_path: Path,
    workers: int,
    skip_api: bool,
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, str]], List[Dict[str, str]]]:
    trajectory_ids = [
        str(replica["traj_id"]) for record in records for replica in record["replicas"]
    ]
    rows: List[Dict[str, Any]] = []
    roles_by_pdb: Dict[str, Dict[str, str]] = {}
    failures: List[Dict[str, str]] = []
    if not skip_api:
        with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
            futures = {
                executor.submit(fetch_itrace_annotation, trajectory_id): trajectory_id
                for trajectory_id in trajectory_ids
            }
            completed = 0
            for future in as_completed(futures):
                trajectory_id = futures[future]
                try:
                    returned_id, roles, profile_rows = future.result()
                    pdb_id = returned_id.split("_run", 1)[0].lower()
                    if roles:
                        roles_by_pdb.setdefault(pdb_id, roles)
                    rows.extend(profile_rows)
                except Exception as error:
                    failures.append({"traj_id": trajectory_id, "error": str(error)})
                completed += 1
                if completed % 50 == 0 or completed == len(futures):
                    print(
                        f"API annotations {completed}/{len(futures)} "
                        f"(failures={len(failures)})",
                        flush=True,
                    )
    if rows:
        official = pd.DataFrame(rows)
        official.to_csv(official_path, index=False)
    elif official_path.exists():
        official = pd.read_csv(official_path, dtype={"residue_number": str})
    else:
        official = pd.DataFrame(
            columns=[
                "traj_id",
                "chain",
                "residue_number",
                "amino_acid_3",
                "region",
                "official_rmsf_angstrom",
            ]
        )
    if roles_by_pdb:
        json_dump(roles_path, roles_by_pdb)
    elif roles_path.exists():
        roles_by_pdb = json.loads(roles_path.read_text(encoding="utf-8"))
    return official, roles_by_pdb, failures


def atom_group_keys(calpha: Any) -> List[Tuple[str, str, str, str]]:
    try:
        insertion_codes = [str(value).strip() for value in calpha.icodes]
    except (AttributeError, TypeError):
        insertion_codes = [""] * len(calpha)
    keys = []
    for index, atom in enumerate(calpha):
        chain = str(atom.chainID).strip()
        if not chain:
            chain = str(atom.segid).strip() or "_"
        keys.append(
            (
                chain,
                str(atom.resid),
                insertion_codes[index],
                str(atom.resname).upper(),
            )
        )
    return keys


def batched_kabsch(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    centered_mobile = mobile - mobile.mean(axis=1, keepdims=True)
    centered_reference = reference - reference.mean(axis=0, keepdims=True)
    covariance = np.einsum(
        "tni,nj->tij", centered_mobile, centered_reference, optimize=True
    )
    left, _, right_t = np.linalg.svd(covariance)
    rotation = np.einsum("tij,tjk->tik", left, right_t)
    negative = np.linalg.det(rotation) < 0
    if np.any(negative):
        left[negative, :, -1] *= -1
        rotation[negative] = np.einsum(
            "tij,tjk->tik", left[negative], right_t[negative]
        )
    return np.einsum("tni,tij->tnj", centered_mobile, rotation, optimize=True)


def chain_fitted_rmsf(coordinates: np.ndarray, chain_ids: Sequence[str]) -> np.ndarray:
    output = np.full(coordinates.shape[1], np.nan, dtype=np.float32)
    chain_array = np.asarray(chain_ids)
    for chain in list(dict.fromkeys(chain_ids)):
        indices = np.flatnonzero(chain_array == chain)
        if len(indices) < 3:
            continue
        aligned = batched_kabsch(coordinates[:, indices], coordinates[0, indices])
        mean = aligned.mean(axis=0)
        output[indices] = np.sqrt(
            np.mean(np.sum((aligned - mean) ** 2, axis=-1), axis=0)
        ).astype(np.float32)
    return output


def compute_replica_rmsf(topology_path: Path, trajectory_path: Path) -> List[Dict[str, Any]]:
    universe = mda.Universe(str(topology_path), str(trajectory_path))
    calpha = universe.select_atoms("protein and name CA")
    coordinates = np.empty(
        (len(universe.trajectory), len(calpha), 3),
        dtype=np.float32,
    )
    for frame_index, _ in enumerate(universe.trajectory):
        coordinates[frame_index] = calpha.positions
    keys = atom_group_keys(calpha)
    chain_ids = [key[0] for key in keys]
    rmsf = chain_fitted_rmsf(coordinates, chain_ids)
    return [
        {
            "chain": key[0],
            "residue_number": key[1],
            "insertion_code": key[2],
            "amino_acid_3": key[3],
            "chainfit_rmsf_angstrom": float(value),
            "frame_count": int(coordinates.shape[0]),
        }
        for key, value in zip(keys, rmsf)
        if np.isfinite(value)
    ]


def compute_chainfit_labels(
    root: Path,
    records: Sequence[Mapping[str, Any]],
    output_path: Path,
) -> Tuple[pd.DataFrame, List[Dict[str, str]]]:
    failures: List[Dict[str, str]] = []
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pdb_id",
        "traj_id",
        "chain",
        "residue_number",
        "insertion_code",
        "amino_acid_3",
        "chainfit_rmsf_angstrom",
        "frame_count",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        total = sum(len(record["replicas"]) for record in records)
        completed = 0
        for record in records:
            for replica in record["replicas"]:
                trajectory_id = str(replica["traj_id"])
                topology = root / "trajectories" / f"{trajectory_id}_topology.pdb"
                trajectory = root / "trajectories" / f"{trajectory_id}_md_processed.xtc"
                try:
                    rows = compute_replica_rmsf(topology, trajectory)
                    for row in rows:
                        writer.writerow(
                            {
                                "pdb_id": record["pdb_id"],
                                "traj_id": trajectory_id,
                                **row,
                            }
                        )
                except Exception as error:
                    failures.append({"traj_id": trajectory_id, "error": str(error)})
                completed += 1
                if completed % 10 == 0 or completed == total:
                    print(
                        f"Chain-fitted RMSF {completed}/{total} "
                        f"(failures={len(failures)})",
                        flush=True,
                    )
                handle.flush()
    return pd.read_csv(
        output_path,
        dtype={"residue_number": str, "insertion_code": str},
        low_memory=False,
    ), failures


def residue_consensus(
    chainfit: pd.DataFrame,
    official: pd.DataFrame,
    output_path: Path,
) -> pd.DataFrame:
    key = ["pdb_id", "chain", "residue_number", "insertion_code", "amino_acid_3"]
    chainfit = chainfit.copy()
    chainfit["pdb_id"] = chainfit["pdb_id"].str.lower()
    chainfit["insertion_code"] = chainfit["insertion_code"].fillna("").astype(str)
    grouped_chainfit = (
        chainfit.groupby(key, dropna=False)["chainfit_rmsf_angstrom"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(
            columns={
                "mean": "target_chainfit_rmsf_angstrom",
                "std": "chainfit_replica_std_angstrom",
                "count": "chainfit_replica_count",
            }
        )
    )
    if official.empty:
        consensus = grouped_chainfit
        consensus["target_official_rmsf_angstrom"] = np.nan
        consensus["official_replica_std_angstrom"] = np.nan
        consensus["official_replica_count"] = 0
        consensus["region"] = ""
    else:
        official = official.copy()
        official["pdb_id"] = official["traj_id"].str.split("_run").str[0].str.lower()
        official["insertion_code"] = ""
        official["residue_number"] = official["residue_number"].astype(str)
        official_key = [
            "pdb_id",
            "chain",
            "residue_number",
            "insertion_code",
            "amino_acid_3",
        ]
        values = (
            official.groupby(official_key, dropna=False)["official_rmsf_angstrom"]
            .agg(["mean", "std", "count"])
            .reset_index()
            .rename(
                columns={
                    "mean": "target_official_rmsf_angstrom",
                    "std": "official_replica_std_angstrom",
                    "count": "official_replica_count",
                }
            )
        )
        regions = (
            official.sort_values("traj_id")
            .groupby(official_key, dropna=False)["region"]
            .first()
            .reset_index()
        )
        values = values.merge(regions, on=official_key, how="left")
        consensus = grouped_chainfit.merge(values, on=key, how="left")
    consensus.to_csv(output_path, index=False)
    return consensus


def pairwise_replica_reliability(
    rows: pd.DataFrame,
    value_column: str,
    output_label: str,
) -> pd.DataFrame:
    records: List[Dict[str, Any]] = []
    for (pdb_id, chain), group in rows.groupby(["pdb_id", "chain"], sort=True):
        maps = {
            trajectory_id: {
                (
                    str(row.residue_number),
                    str(row.insertion_code or ""),
                    str(row.amino_acid_3),
                ): float(getattr(row, value_column))
                for row in subset.itertuples(index=False)
                if np.isfinite(float(getattr(row, value_column)))
            }
            for trajectory_id, subset in group.groupby("traj_id")
        }
        trajectory_ids = sorted(maps)
        correlations = []
        for left_index in range(len(trajectory_ids)):
            for right_index in range(left_index + 1, len(trajectory_ids)):
                left = maps[trajectory_ids[left_index]]
                right = maps[trajectory_ids[right_index]]
                common = sorted(set(left) & set(right))
                if len(common) < 4:
                    continue
                correlation = spearmanr(
                    [left[key] for key in common],
                    [right[key] for key in common],
                ).statistic
                if np.isfinite(correlation):
                    correlations.append(float(correlation))
        if correlations:
            records.append(
                {
                    "pdb_id": pdb_id,
                    "chain": chain,
                    "label_definition": output_label,
                    "pair_count": len(correlations),
                    "mean_pairwise_spearman": float(np.mean(correlations)),
                    "median_pairwise_spearman": float(np.median(correlations)),
                }
            )
    chain_rows = pd.DataFrame(records)
    if chain_rows.empty:
        return chain_rows
    complex_rows = (
        chain_rows.groupby(["pdb_id", "label_definition"], as_index=False)
        .agg(
            chain_count=("chain", "count"),
            mean_pairwise_spearman=("mean_pairwise_spearman", "mean"),
            median_pairwise_spearman=("median_pairwise_spearman", "median"),
        )
    )
    complex_rows["chain"] = "__complex_macro__"
    return pd.concat([chain_rows, complex_rows], ignore_index=True)


def write_fastprotflex_manifest(
    records: Sequence[Mapping[str, Any]],
    path: Path,
) -> Dict[str, Tuple[str, str]]:
    sample_map: Dict[str, Tuple[str, str]] = {}
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "pdb_path", "chain_filter"])
        for record in records:
            structure = parse_pdb(str(record["topology_path"]))
            for chain in structure.unique_chains():
                safe_chain = chain if chain != "_" else "blank"
                sample_id = f"{record['pdb_id']}__{safe_chain}"
                sample_map[sample_id] = (str(record["pdb_id"]), chain)
                writer.writerow([sample_id, record["topology_path"], chain])
    return sample_map


def run_fastprotflex(
    manifest_path: Path,
    prediction_path: Path,
    timing_path: Path,
    rscript: Path,
    script: Path,
) -> None:
    environment = os.environ.copy()
    environment.setdefault(
        "R_LIBS_USER",
        str(Path(os.environ.get("LOCALAPPDATA", "")) / "R" / "win-library" / "4.6"),
    )
    subprocess.run(
        [
            str(rscript),
            str(script),
            str(manifest_path),
            str(prediction_path),
            str(timing_path),
        ],
        check=True,
        env=environment,
    )


@torch.no_grad()
def run_flexppi(
    records: Sequence[Mapping[str, Any]],
    checkpoint: Path,
    prediction_path: Path,
    timing_path: Path,
) -> None:
    model, model_config, device = load_model(str(checkpoint), "auto")
    rmsf_index = model_config.metrics.index("rmsf")
    rows: List[Dict[str, Any]] = []
    timings: List[Dict[str, Any]] = []
    if records:
        first = records[0]
        structure = parse_pdb(str(first["topology_path"]))
        warmup = collate_states(
            [
                _structure_example(
                    str(first["topology_path"]),
                    structure.unique_chains(),
                    310.0,
                    200.0,
                    1.0,
                )
            ]
        )
        model.forward_state(move_to_device(warmup, device))
    for index, record in enumerate(records, start=1):
        started = time.perf_counter()
        try:
            structure = parse_pdb(str(record["topology_path"]))
            example = _structure_example(
                str(record["topology_path"]),
                structure.unique_chains(),
                310.0,
                200.0,
                1.0,
            )
            batch = move_to_device(collate_states([example]), device)
            output = model.forward_state(batch)
            prediction = output["mean"][0, :, rmsf_index].cpu().numpy()
            uncertainty = torch.exp(
                0.5 * output["log_variance"][0, :, rmsf_index]
            ).cpu().numpy()
            for residue_index in range(structure.num_residues):
                author_id = structure.author_residue_ids[residue_index]
                digits = "".join(
                    character
                    for character in author_id
                    if character.isdigit() or character == "-"
                )
                insertion = author_id[len(digits) :] if author_id.startswith(digits) else ""
                rows.append(
                    {
                        "pdb_id": record["pdb_id"],
                        "residue_index": residue_index,
                        "chain": structure.chain_ids[residue_index],
                        "residue_number": digits or author_id,
                        "insertion_code": insertion,
                        "amino_acid_1": structure.sequence[residue_index],
                        "amino_acid_3": structure.residue_names[residue_index],
                        "flexppi_rmsf_angstrom": float(prediction[residue_index]),
                        "flexppi_std_angstrom": float(uncertainty[residue_index]),
                    }
                )
            status, error = "ok", ""
        except Exception as exception:
            status, error = "error", str(exception)
        timings.append(
            {
                "pdb_id": record["pdb_id"],
                "residue_count": (
                    structure.num_residues if "structure" in locals() else np.nan
                ),
                "end_to_end_seconds": time.perf_counter() - started,
                "status": status,
                "error": error,
            }
        )
        if index % 10 == 0 or index == len(records):
            print(
                f"FlexPPI {index}/{len(records)} "
                f"(failures={sum(row['status'] != 'ok' for row in timings)})",
                flush=True,
            )
    pd.DataFrame(rows).to_csv(prediction_path, index=False)
    pd.DataFrame(timings).to_csv(timing_path, index=False)


def normalize_fastprotflex_predictions(
    prediction_path: Path,
    sample_map: Mapping[str, Tuple[str, str]],
) -> pd.DataFrame:
    frame = pd.read_csv(
        prediction_path,
        dtype={"residue_number": str, "insertion_code": str},
    )
    frame["pdb_id"] = frame["sample_id"].map(lambda value: sample_map[value][0])
    frame["requested_chain"] = frame["sample_id"].map(lambda value: sample_map[value][1])
    frame["insertion_code"] = frame["insertion_code"].fillna("")
    return frame


def interface_mask(
    structure: ProteinStructure,
    roles: Mapping[str, str],
    cutoff: float = 8.0,
) -> np.ndarray:
    mask = np.zeros(structure.num_residues, dtype=bool)
    present = structure.backbone_mask[:, 1]
    coordinates = structure.backbone[:, 1].astype(np.float64)
    pmhc = np.asarray(
        [
            present[index]
            and roles.get(structure.chain_ids[index]) in {"mhc", "b2m", "peptide"}
            for index in range(structure.num_residues)
        ]
    )
    tcr = np.asarray(
        [
            present[index]
            and roles.get(structure.chain_ids[index]) in {"tcr_alpha", "tcr_beta"}
            for index in range(structure.num_residues)
        ]
    )
    left = np.flatnonzero(pmhc)
    right = np.flatnonzero(tcr)
    if len(left) == 0 or len(right) == 0:
        return mask
    neighbours = cKDTree(coordinates[right]).query_ball_point(coordinates[left], cutoff)
    for left_position, local_right in enumerate(neighbours):
        if local_right:
            mask[left[left_position]] = True
            mask[right[np.asarray(local_right, dtype=int)]] = True
    return mask


def normalize_region(value: Any) -> str:
    text = str(value or "").strip().lower()
    replacements = {
        "α": "alpha",
        "β": "beta",
        "î±": "alpha",
        "î²": "beta",
    }
    for source, target in replacements.items():
        text = text.replace(source, target)
    if "cdr1" in text:
        return "cdr1"
    if "cdr2" in text:
        return "cdr2"
    if "cdr3" in text:
        return "cdr3"
    if "tcr-fr" in text or "tcr framework" in text:
        return "tcr_framework"
    if "peptide" in text:
        return "peptide"
    if "2m" in text:
        return "b2m"
    if "helix" in text or "mhc floor" in text:
        return "mhc"
    return "other"


def static_annotation_rows(
    records: Sequence[Mapping[str, Any]],
    roles_by_pdb: Dict[str, Dict[str, str]],
) -> pd.DataFrame:
    rows = []
    for record in records:
        structure = parse_pdb(str(record["topology_path"]))
        roles = roles_by_pdb.get(str(record["pdb_id"]), {})
        if len(roles) < len(structure.unique_chains()):
            inferred = infer_roles_from_structure(
                structure, str(record.get("peptide_seq") or "")
            )
            roles = {**inferred, **roles}
            roles_by_pdb[str(record["pdb_id"])] = roles
        interface = interface_mask(structure, roles)
        for residue_index in range(structure.num_residues):
            author_id = structure.author_residue_ids[residue_index]
            digits = "".join(
                character
                for character in author_id
                if character.isdigit() or character == "-"
            )
            insertion = author_id[len(digits) :] if author_id.startswith(digits) else ""
            rows.append(
                {
                    "pdb_id": record["pdb_id"],
                    "chain": structure.chain_ids[residue_index],
                    "residue_number": digits or author_id,
                    "insertion_code": insertion,
                    "amino_acid_3": structure.residue_names[residue_index],
                    "chain_role": roles.get(structure.chain_ids[residue_index], "unknown"),
                    "is_tcr_pmhc_interface_8A": bool(interface[residue_index]),
                }
            )
    return pd.DataFrame(rows)


def common_residue_table(
    consensus: pd.DataFrame,
    flex_path: Path,
    fast_path: Path,
    sample_map: Mapping[str, Tuple[str, str]],
    annotations: pd.DataFrame,
    output_path: Path,
) -> pd.DataFrame:
    key = ["pdb_id", "chain", "residue_number", "insertion_code", "amino_acid_3"]
    flex = pd.read_csv(
        flex_path, dtype={"residue_number": str, "insertion_code": str}
    )
    fast = normalize_fastprotflex_predictions(fast_path, sample_map)
    for frame in (consensus, flex, fast, annotations):
        frame["pdb_id"] = frame["pdb_id"].astype(str).str.lower()
        frame["residue_number"] = frame["residue_number"].astype(str)
        frame["insertion_code"] = frame["insertion_code"].fillna("").astype(str)
        frame["amino_acid_3"] = frame["amino_acid_3"].astype(str).str.upper()
    fast = fast.rename(columns={"fastprotflex_score": "fastprotflex_score"})
    common = (
        consensus.merge(
            flex[key + ["flexppi_rmsf_angstrom", "flexppi_std_angstrom"]],
            on=key,
            how="inner",
        )
        .merge(
            fast[key + ["fastprotflex_score"]],
            on=key,
            how="inner",
        )
        .merge(annotations, on=key, how="left")
    )
    common["region_group"] = common["region"].map(normalize_region)
    common.to_csv(output_path, index=False)
    return common


def subset_mask(frame: pd.DataFrame, subset: str) -> np.ndarray:
    if subset == "overall":
        return np.ones(len(frame), dtype=bool)
    if subset == "interface":
        return frame["is_tcr_pmhc_interface_8A"].fillna(False).to_numpy(bool)
    if subset == "non_interface":
        return ~frame["is_tcr_pmhc_interface_8A"].fillna(False).to_numpy(bool)
    if subset == "tcr":
        return frame["chain_role"].isin(["tcr_alpha", "tcr_beta"]).to_numpy()
    if subset == "pmhc":
        return frame["chain_role"].isin(["mhc", "b2m", "peptide"]).to_numpy()
    if subset in {"mhc", "b2m", "peptide", "tcr_alpha", "tcr_beta"}:
        return (frame["chain_role"] == subset).to_numpy()
    if subset == "cdr_all":
        return frame["region_group"].isin(["cdr1", "cdr2", "cdr3"]).to_numpy()
    if subset in {"cdr1", "cdr2", "cdr3", "tcr_framework"}:
        return (frame["region_group"] == subset).to_numpy()
    raise ValueError(f"Unknown subset: {subset}")


def flatten_metrics(prefix: str, metrics: Mapping[str, Optional[float]]) -> Dict[str, Any]:
    return {f"{prefix}__{name}": value for name, value in metrics.items()}


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg false-discovery-rate adjustment with NaN preservation."""

    result = pd.Series(np.nan, index=values.index, dtype=float)
    finite = values.dropna().astype(float)
    if finite.empty:
        return result
    ordered = finite.sort_values()
    count = len(ordered)
    adjusted = ordered.to_numpy() * count / np.arange(1, count + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[ordered.index] = np.clip(adjusted, 0.0, 1.0)
    return result


def calculate_metrics(
    common: pd.DataFrame,
    reliability: pd.DataFrame,
    contaminated: set[str],
    bootstrap_repetitions: int,
    seed: int,
    output_dir: Path,
) -> Dict[str, Any]:
    subsets = (
        "overall",
        "interface",
        "non_interface",
        "tcr",
        "pmhc",
        "mhc",
        "b2m",
        "peptide",
        "tcr_alpha",
        "tcr_beta",
        "cdr_all",
        "cdr1",
        "cdr2",
        "cdr3",
        "tcr_framework",
    )
    chain_rows: List[Dict[str, Any]] = []
    for label_name, label_column in LABEL_COLUMNS.items():
        for subset in subsets:
            selected = common.loc[subset_mask(common, subset)]
            for (pdb_id, chain), group in selected.groupby(["pdb_id", "chain"], sort=True):
                valid = (
                    np.isfinite(group[label_column])
                    & np.isfinite(group["flexppi_rmsf_angstrom"])
                    & np.isfinite(group["fastprotflex_score"])
                )
                group = group.loc[valid]
                if len(group) < 4:
                    continue
                target = group[label_column].to_numpy(float)
                flex = _system_metrics(
                    group["flexppi_rmsf_angstrom"].to_numpy(float),
                    target,
                    absolute_scale=True,
                )
                fast = _system_metrics(
                    group["fastprotflex_score"].to_numpy(float),
                    target,
                    absolute_scale=False,
                )
                if not flex or not fast:
                    continue
                chain_rows.append(
                    {
                        "pdb_id": pdb_id,
                        "chain": chain,
                        "label_definition": label_name,
                        "subset": subset,
                        "residue_count": len(group),
                        **flatten_metrics("flexppi", flex),
                        **flatten_metrics("fastprotflex", fast),
                    }
                )
    per_chain = pd.DataFrame(chain_rows)
    per_chain.to_csv(output_dir / "per_chain_metrics.csv", index=False)
    metric_columns = [
        column
        for column in per_chain.columns
        if column.startswith("flexppi__") or column.startswith("fastprotflex__")
    ]
    aggregations = {column: "mean" for column in metric_columns}
    aggregations["residue_count"] = "sum"
    per_complex = (
        per_chain.groupby(
            ["pdb_id", "label_definition", "subset"], as_index=False
        )
        .agg(aggregations)
        .sort_values(["label_definition", "subset", "pdb_id"])
    )
    per_complex.to_csv(output_dir / "per_complex_metrics.csv", index=False)

    complex_reliability = reliability.loc[
        reliability["chain"] == "__complex_macro__"
    ].copy()
    stable_chainfit = set(
        complex_reliability.loc[
            (complex_reliability["label_definition"] == "chain_fitted_internal_rmsf")
            & (complex_reliability["mean_pairwise_spearman"] >= 0.5),
            "pdb_id",
        ].astype(str)
    )
    all_ids = set(per_complex["pdb_id"].astype(str))
    clean_ids = all_ids - {value.lower() for value in contaminated}
    cohort_ids = {
        "all_245": all_ids,
        "external_clean": clean_ids,
        "replica_stable_external_clean": clean_ids & stable_chainfit,
    }
    summary_rows: List[Dict[str, Any]] = []
    paired_rows: List[Dict[str, Any]] = []
    for cohort_name, identifiers in cohort_ids.items():
        cohort = per_complex.loc[per_complex["pdb_id"].isin(identifiers)]
        for (label_name, subset), group in cohort.groupby(
            ["label_definition", "subset"], sort=True
        ):
            for method in ("fastprotflex", "flexppi"):
                for metric in METRICS + (
                    ("absolute_mae_angstrom", "absolute_rmse_angstrom")
                    if method == "flexppi"
                    else ()
                ):
                    column = f"{method}__{metric}"
                    if column not in group:
                        continue
                    values = group[column].dropna().to_numpy(float)
                    if len(values) == 0:
                        continue
                    interval = _bootstrap_mean_ci(
                        values.tolist(), bootstrap_repetitions, seed
                    )
                    summary_rows.append(
                        {
                            "cohort": cohort_name,
                            "label_definition": label_name,
                            "subset": subset,
                            "method": method,
                            "metric": metric,
                            "macro_mean": float(np.mean(values)),
                            "macro_median": float(np.median(values)),
                            "ci95_low": interval[0] if interval else np.nan,
                            "ci95_high": interval[1] if interval else np.nan,
                            "complex_count": len(values),
                        }
                    )
            for metric, direction in ADVANTAGE_DIRECTION.items():
                flex_column = f"flexppi__{metric}"
                fast_column = f"fastprotflex__{metric}"
                paired = group[[flex_column, fast_column]].dropna()
                if paired.empty:
                    continue
                differences = direction * (
                    paired[flex_column].to_numpy(float)
                    - paired[fast_column].to_numpy(float)
                )
                interval = _bootstrap_mean_ci(
                    differences.tolist(), bootstrap_repetitions, seed
                )
                nonzero = differences[np.abs(differences) > 1e-12]
                p_value = (
                    float(wilcoxon(nonzero, alternative="two-sided").pvalue)
                    if len(nonzero)
                    else np.nan
                )
                paired_rows.append(
                    {
                        "cohort": cohort_name,
                        "label_definition": label_name,
                        "subset": subset,
                        "metric": metric,
                        "advantage_definition": "positive means FlexPPI is better",
                        "mean_advantage": float(np.mean(differences)),
                        "median_advantage": float(np.median(differences)),
                        "ci95_low": interval[0] if interval else np.nan,
                        "ci95_high": interval[1] if interval else np.nan,
                        "win_rate": float(np.mean(differences > 1e-12)),
                        "wilcoxon_two_sided_p": p_value,
                        "complex_count": len(differences),
                    }
                )
    summary = pd.DataFrame(summary_rows)
    paired = pd.DataFrame(paired_rows)
    if not paired.empty:
        paired["fdr_q_across_subsets_for_metric"] = (
            paired.groupby(
                ["cohort", "label_definition", "metric"],
                group_keys=False,
            )["wilcoxon_two_sided_p"]
            .apply(benjamini_hochberg)
            .sort_index()
        )
        paired["fdr_q_across_all_subset_metric_tests"] = (
            paired.groupby(
                ["cohort", "label_definition"],
                group_keys=False,
            )["wilcoxon_two_sided_p"]
            .apply(benjamini_hochberg)
            .sort_index()
        )
    summary.to_csv(output_dir / "summary_metrics.csv", index=False)
    paired.to_csv(output_dir / "paired_tests.csv", index=False)
    return {
        "all_complexes": len(all_ids),
        "external_clean_complexes": len(clean_ids),
        "replica_stable_external_clean_complexes": len(clean_ids & stable_chainfit),
        "per_chain_rows": len(per_chain),
        "per_complex_rows": len(per_complex),
    }


def write_protocol(
    path: Path,
    records: Sequence[Mapping[str, Any]],
    contaminated: set[str],
    api_failures: Sequence[Mapping[str, str]],
    label_failures: Sequence[Mapping[str, str]],
    metrics_summary: Mapping[str, Any],
    input_length_summary: Mapping[str, Any],
    data_coverage_summary: Mapping[str, Any],
) -> None:
    path.write_text(
        f"""# ITrace external benchmark: FlexPPI versus FastProtFlex

## Dataset and statistical unit

ITrace contributes {len(records)} TCR-pMHC complexes. Each complex normally has
three independent 200 ns trajectories at 310 K. Per-residue labels are first
computed per trajectory and then averaged across replicas. All confidence
intervals and paired tests use the PDB complex, not the trajectory or residue,
as the independent unit.

## Static inputs

Both models receive the same run-1 frame-0 protein topology. FlexPPI receives
the complete multichain complex. FastProtFlex is run independently on every
chain because its published model is a single-chain predictor. FlexPPI uses
the static pretraining checkpoint, temperature 310 K, window 200 ns, and the
DynaRepo-like MD condition code selected before examining ITrace results.

The complete ITrace inputs contain {input_length_summary["minimum"]} to
{input_length_summary["maximum"]} residues (median
{input_length_summary["median"]:.0f}). The pretraining data loader used a
{input_length_summary["training_max_residues"]}-residue maximum; therefore
{input_length_summary["above_training_max_count"]} of {len(records)} complexes
are length-extrapolation tests. Inference itself remains full-complex and does
not crop these structures.

## Ground truth

The primary target is C-alpha RMSF after fitting each chain independently,
matching the FlexPPI training-label definition and excluding inter-chain rigid
motion. The official ITrace whole-complex-fitted RMSF profile is evaluated as a
secondary definition-sensitivity analysis. Three replica profiles are averaged
per residue; replica-to-replica Spearman agreement is reported separately.
Although ITrace documentation describes 1,001 saved frames per trajectory, the
local processed XTC files contain {data_coverage_summary["minimum_frames"]} to
{data_coverage_summary["maximum_frames"]} readable frames (median
{data_coverage_summary["median_frames"]:.0f}); every readable frame is used.
The primary target maps all {data_coverage_summary["common_residues"]} common
residues. The official secondary profile maps
{data_coverage_summary["official_target_residues"]} residues
({100 * data_coverage_summary["official_target_fraction"]:.3f}% coverage);
missing official residues are omitted only from that sensitivity analysis.

## Metrics and aggregation

Spearman, Pearson, within-chain z-MAE/z-RMSE, and top-20% AUROC, average
precision and precision-at-K are calculated on the exact common residue set.
Metrics are calculated within each chain and then averaged with equal chain
weight inside each complex. Dataset summaries are equal-complex macro means
with complex bootstrap 95% confidence intervals. FlexPPI-versus-FastProtFlex
differences use two-sided paired Wilcoxon tests. Absolute Angstrom error is
reported only for FlexPPI because FastProtFlex outputs a unitless standardized
score.

## External-overlap control

Exact PDB identifiers found in FlexPPI training manifests: {", ".join(sorted(contaminated)) or "none"}.
The main external-clean result excludes these complexes; the all-system result
is retained as a sensitivity analysis.

## Fine-grained analyses

Results are stratified by TCR/pMHC component, CDR1/2/3 and TCR framework,
and by an 8 Angstrom C-alpha TCR-pMHC interface mask. Chain roles and CDR
annotations come from the ITrace identity/RMSF-profile metadata rather than
assuming fixed PDB chain letters.

## Execution audit

- API annotation failures: {len(api_failures)}
- Local trajectory-label failures: {len(label_failures)}
- External-clean evaluable complexes: {metrics_summary.get("external_clean_complexes")}
- Replica-stable external-clean complexes: {metrics_summary.get("replica_stable_external_clean_complexes")}
""",
        encoding="utf-8",
    )


def main() -> None:
    arguments = parse_arguments()
    root = Path(arguments.itrace_root)
    output = Path(arguments.output_dir)
    raw = output / "raw"
    metrics_dir = output / "metrics"
    logs = output / "logs"
    for directory in (raw, metrics_dir, logs):
        directory.mkdir(parents=True, exist_ok=True)

    records = load_itrace_records(root, arguments.limit)
    checkpoint = Path(arguments.checkpoint)
    input_residue_counts = np.asarray(
        [
            parse_pdb(str(record["topology_path"])).num_residues
            for record in records
        ],
        dtype=int,
    )
    resolved_config_path = checkpoint.parent / "resolved_config.json"
    training_max_residues = 768
    if resolved_config_path.exists():
        resolved_config = json.loads(resolved_config_path.read_text(encoding="utf-8"))
        training_max_residues = int(
            resolved_config.get("data", {}).get(
                "max_residues",
                training_max_residues,
            )
        )
    input_length_summary = {
        "minimum": int(input_residue_counts.min()),
        "median": float(np.median(input_residue_counts)),
        "maximum": int(input_residue_counts.max()),
        "training_max_residues": training_max_residues,
        "above_training_max_count": int(
            np.sum(input_residue_counts > training_max_residues)
        ),
    }
    manifest_csv = raw / "itrace_complex_manifest.csv"
    flattened_manifest = []
    for record in records:
        flattened_manifest.append(
            {
                key: value
                for key, value in record.items()
                if key != "replicas"
            }
            | {
                "trajectory_ids": ";".join(
                    str(replica["traj_id"]) for replica in record["replicas"]
                ),
                "replica_count": len(record["replicas"]),
            }
        )
    pd.DataFrame(flattened_manifest).to_csv(manifest_csv, index=False)

    training_ids = training_pdb_ids(
        [
            Path("artifacts/data/atlas_states.jsonl"),
            Path("artifacts/data/dynarepo_states.jsonl"),
        ]
    )
    itrace_ids = {str(record["pdb_id"]).upper() for record in records}
    contaminated = itrace_ids & training_ids

    official_path = raw / "itrace_official_rmsf_per_replica.csv"
    roles_path = raw / "itrace_chain_roles.json"
    if arguments.reuse and official_path.exists() and roles_path.exists():
        official = pd.read_csv(official_path, dtype={"residue_number": str})
        roles_by_pdb = json.loads(roles_path.read_text(encoding="utf-8"))
        api_failures: List[Dict[str, str]] = []
    else:
        official, roles_by_pdb, api_failures = collect_api_data(
            records,
            official_path,
            roles_path,
            arguments.api_workers,
            arguments.skip_api,
        )

    chainfit_path = raw / "chainfit_rmsf_per_replica.csv"
    if arguments.reuse and chainfit_path.exists():
        chainfit = pd.read_csv(
            chainfit_path,
            dtype={"residue_number": str, "insertion_code": str},
            low_memory=False,
        )
        label_failures: List[Dict[str, str]] = []
    else:
        chainfit, label_failures = compute_chainfit_labels(
            root, records, chainfit_path
        )

    consensus = residue_consensus(
        chainfit,
        official,
        raw / "rmsf_consensus_per_residue.csv",
    )
    reliability_frames = [
        pairwise_replica_reliability(
            chainfit.assign(
                insertion_code=chainfit["insertion_code"].fillna("").astype(str)
            ),
            "chainfit_rmsf_angstrom",
            "chain_fitted_internal_rmsf",
        )
    ]
    if not official.empty:
        official_for_reliability = official.copy()
        official_for_reliability["pdb_id"] = (
            official_for_reliability["traj_id"].str.split("_run").str[0].str.lower()
        )
        official_for_reliability["insertion_code"] = ""
        reliability_frames.append(
            pairwise_replica_reliability(
                official_for_reliability,
                "official_rmsf_angstrom",
                "itrace_official_complex_fitted_rmsf",
            )
        )
    reliability = pd.concat(
        [frame for frame in reliability_frames if not frame.empty],
        ignore_index=True,
    )
    reliability.to_csv(metrics_dir / "replica_reliability.csv", index=False)

    json_dump(roles_path, roles_by_pdb)
    annotations = static_annotation_rows(records, roles_by_pdb)
    # Add the authoritative region name from the consensus table at merge time.

    fast_manifest = raw / "fastprotflex_chain_manifest.tsv"
    sample_map = write_fastprotflex_manifest(records, fast_manifest)
    fast_predictions = raw / "fastprotflex_chain_predictions.csv"
    fast_timings = raw / "fastprotflex_chain_timings.csv"
    if not (
        arguments.reuse and fast_predictions.exists() and fast_timings.exists()
    ):
        run_fastprotflex(
            fast_manifest,
            fast_predictions,
            fast_timings,
            Path(arguments.rscript),
            Path(arguments.fastprotflex_script),
        )

    flex_predictions = raw / "flexppi_full_complex_predictions.csv"
    flex_timings = raw / "flexppi_full_complex_timings.csv"
    if not (
        arguments.reuse and flex_predictions.exists() and flex_timings.exists()
    ):
        run_flexppi(records, checkpoint, flex_predictions, flex_timings)

    common = common_residue_table(
        consensus,
        flex_predictions,
        fast_predictions,
        sample_map,
        annotations,
        raw / "common_residue_predictions.csv",
    )
    trajectory_frame_counts = (
        chainfit.groupby("traj_id", sort=False)["frame_count"].first().astype(int)
    )
    data_coverage_summary = {
        "minimum_frames": int(trajectory_frame_counts.min()),
        "median_frames": float(trajectory_frame_counts.median()),
        "maximum_frames": int(trajectory_frame_counts.max()),
        "trajectories_with_nominal_1001_frames": int(
            (trajectory_frame_counts == 1001).sum()
        ),
        "trajectories_with_other_frame_counts": int(
            (trajectory_frame_counts != 1001).sum()
        ),
        "common_residues": int(len(common)),
        "primary_target_residues": int(
            common["target_chainfit_rmsf_angstrom"].notna().sum()
        ),
        "official_target_residues": int(
            common["target_official_rmsf_angstrom"].notna().sum()
        ),
        "official_target_fraction": float(
            common["target_official_rmsf_angstrom"].notna().mean()
        ),
    }
    metrics_summary = calculate_metrics(
        common,
        reliability,
        contaminated,
        arguments.bootstrap_repetitions,
        arguments.seed,
        metrics_dir,
    )

    fast_timing_frame = pd.read_csv(fast_timings)
    flex_timing_frame = pd.read_csv(flex_timings)
    audit = {
        "dataset": {
            "name": "ITrace",
            "complexes_requested": len(records),
            "trajectories_requested": sum(len(record["replicas"]) for record in records),
            "temperature_K": 310.0,
            "trajectory_duration_ns": 200.0,
            "primary_statistical_unit": "PDB complex",
            "input_residue_count": input_length_summary,
            "local_data_coverage": data_coverage_summary,
        },
        "models": {
            "flexppi": {
                "checkpoint": str(checkpoint.resolve()),
                "checkpoint_sha256": sha256(checkpoint),
                "input": "full multichain run-1 frame-0 topology",
                "condition": {
                    "temperature_K": 310.0,
                    "window_ns": 200.0,
                    "dataset_code": 1.0,
                },
            },
            "fastprotflex": {
                "input": "each chain independently from the same topology",
                "implementation": str(Path(arguments.fastprotflex_script).resolve()),
            },
        },
        "software": {
            "python": sys.version.split()[0],
            "torch": torch.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "MDAnalysis": mda.__version__,
            "rscript": str(Path(arguments.rscript).resolve()),
        },
        "training_overlap_exact_pdb": sorted(contaminated),
        "failures": {
            "api": api_failures,
            "trajectory_labels": label_failures,
            "fastprotflex": fast_timing_frame.loc[
                fast_timing_frame["status"] != "ok", ["sample_id", "error"]
            ].to_dict("records"),
            "flexppi": flex_timing_frame.loc[
                flex_timing_frame["status"] != "ok", ["pdb_id", "error"]
            ].to_dict("records"),
        },
        "coverage": metrics_summary,
        "runtime": {
            "fastprotflex_chain_total_seconds": float(
                fast_timing_frame["end_to_end_seconds"].dropna().sum()
            ),
            "fastprotflex_chain_median_seconds": float(
                fast_timing_frame["end_to_end_seconds"].dropna().median()
            ),
            "flexppi_complex_total_seconds": float(
                flex_timing_frame["end_to_end_seconds"].dropna().sum()
            ),
            "flexppi_complex_median_seconds": float(
                flex_timing_frame["end_to_end_seconds"].dropna().median()
            ),
        },
        "artifacts": {
            "complex_manifest": str(manifest_csv.resolve()),
            "chainfit_replica_labels": str(chainfit_path.resolve()),
            "official_replica_labels": str(official_path.resolve()),
            "consensus_labels": str(
                (raw / "rmsf_consensus_per_residue.csv").resolve()
            ),
            "common_residue_predictions": str(
                (raw / "common_residue_predictions.csv").resolve()
            ),
            "summary_metrics": str((metrics_dir / "summary_metrics.csv").resolve()),
            "paired_tests": str((metrics_dir / "paired_tests.csv").resolve()),
        },
    }
    json_dump(output / "itrace_benchmark_report.json", audit)
    write_protocol(
        output / "ITRACE_EVALUATION_PROTOCOL.md",
        records,
        contaminated,
        api_failures,
        label_failures,
        metrics_summary,
        input_length_summary,
        data_coverage_summary,
    )
    print(json.dumps(audit["coverage"], indent=2), flush=True)


if __name__ == "__main__":
    main()
