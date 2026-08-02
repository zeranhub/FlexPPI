"""Residue-level multi-metric context analysis for the ITrace TCR-pMHC set.

The same frozen FlexPPI checkpoint is evaluated on four kinds of fixed-coordinate
input:

1. the complete TCR-pMHC complex;
2. the pMHC subcomplex (MHC + beta-2-microglobulin + peptide);
3. the alpha/beta TCR heterodimer; and
4. each of the five chains in isolation.

The pMHC-to-complete and TCR-to-complete comparisons isolate the model response
to adding the binding partner while preserving the experimentally observed
bound coordinates.  They are therefore *binding-context responses*, not true
apo-to-holo measurements.  The complete-to-isolated comparison is retained as
the original context ablation, but isolated short peptides are flagged as
out-of-distribution and are not used for biological binding claims.

Outputs include RMSF, side-chain rotamer entropy (chi_entropy, H_chi), SASA,
heavy-atom contacts, peptide position, TCR CDR annotation, clustered bootstrap
intervals and exact agreement checks against the published RMSF inference run.
Glycine and alanine are explicitly excluded from H_chi summaries because they
do not have a defined side-chain chi torsion in the training labels.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import re
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr
import torch

from flexppi.data.dataset import collate_states
from flexppi.data.structures import ProteinStructure, parse_pdb
from flexppi.inference import _structure_example, load_model
from flexppi.trainer import move_to_device


MODEL_METRICS = ("rmsf", "chi_entropy", "sasa")
PMHC_ROLES = {"mhc", "b2m", "peptide"}
TCR_ROLES = {"tcr_alpha", "tcr_beta"}
CHI_UNDEFINED_AA = {"A", "G"}
CONTACT_CUTOFF_ANGSTROM = 5.0
EXACT_TRAINING_OVERLAP = {"3kxf"}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-dir", default="output/itrace_benchmark")
    parser.add_argument("--context-dir", default="output/itrace_zero_shot_context")
    parser.add_argument("--checkpoint", default="artifacts/pretrain/best.pt")
    parser.add_argument("--output-dir", default="output/itrace_residue_biology")
    parser.add_argument("--device", default="auto")
    parser.add_argument("--temperature", type=float, default=310.0)
    parser.add_argument("--window-ns", type=float, default=200.0)
    parser.add_argument("--dataset-code", type=float, default=1.0)
    parser.add_argument("--bootstrap-repetitions", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--reuse", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_dump(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def _renumber_atom(line: str, serial: int) -> str:
    value = line.rstrip("\r\n").ljust(80)
    return value[:6] + f"{serial:5d}" + value[11:] + "\n"


def _header_before_coordinates(lines: Sequence[str]) -> List[str]:
    header: List[str] = []
    for line in lines:
        record = line[:6].strip()
        if record in {"MODEL", "ATOM", "HETATM"}:
            break
        if record not in {"END", "MASTER", "CONECT"}:
            header.append(line if line.endswith("\n") else line + "\n")
    return header


def extract_chains(source: Path, destination: Path, wanted: Sequence[str]) -> ProteinStructure:
    """Retain selected protein chains without changing any atom coordinate."""

    wanted_set = set(wanted)
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines(True)
    atoms = [
        line
        for line in lines
        if line.startswith("ATOM  ") and line[21:22].strip() in wanted_set
    ]
    if not atoms:
        raise ValueError(f"No requested chains {wanted!r} in {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        handle.writelines(_header_before_coordinates(lines))
        handle.write(
            "REMARK 900 FLEXPPI FIXED-COORDINATE CONTEXT INPUT: RETAINED CHAINS "
            + ",".join(wanted)
            + "\n"
        )
        previous_chain = None
        serial = 1
        for line in atoms:
            chain = line[21:22].strip()
            if previous_chain is not None and chain != previous_chain:
                handle.write(f"TER   {serial:5d}\n")
                serial += 1
            handle.write(_renumber_atom(line, serial))
            serial += 1
            previous_chain = chain
        handle.write("END\n")
    structure = parse_pdb(str(destination))
    observed = set(structure.unique_chains())
    if observed != wanted_set:
        raise ValueError(f"Expected chains {wanted_set!r}; parsed {observed!r}")
    return structure


def author_residue_parts(author_id: str) -> Tuple[str, str]:
    match = re.match(r"^(-?\d+)(.*)$", str(author_id).strip())
    if match is None:
        return str(author_id).strip(), ""
    return match.group(1), match.group(2)


@torch.no_grad()
def predict_condition(
    model: torch.nn.Module,
    model_config: Any,
    device: torch.device,
    pdb_id: str,
    condition: str,
    pdb_path: Path,
    chain_roles: Mapping[str, str],
    temperature: float,
    window_ns: float,
    dataset_code: float,
) -> pd.DataFrame:
    structure = parse_pdb(str(pdb_path))
    example = _structure_example(
        str(pdb_path),
        structure.unique_chains(),
        temperature,
        window_ns,
        dataset_code,
    )
    batch = move_to_device(collate_states([example]), device)
    output = model.forward_state(batch)
    mean = output["mean"][0].detach().cpu().numpy().astype(np.float64)
    std = (
        torch.exp(0.5 * output["log_variance"][0])
        .detach()
        .cpu()
        .numpy()
        .astype(np.float64)
    )
    valid = (batch["target_mask"][0] & output["mask"][0]).detach().cpu().numpy()
    metric_indices = {metric: model_config.metrics.index(metric) for metric in MODEL_METRICS}
    rows: List[Dict[str, Any]] = []
    for residue_index in np.flatnonzero(valid):
        chain = structure.chain_ids[residue_index]
        number, insertion = author_residue_parts(
            structure.author_residue_ids[residue_index]
        )
        row: Dict[str, Any] = {
            "pdb_id": pdb_id,
            "condition": condition,
            "residue_index": int(residue_index),
            "chain": chain,
            "chain_role": chain_roles[chain],
            "author_residue_id": structure.author_residue_ids[residue_index],
            "residue_number": number,
            "insertion_code": insertion,
            "amino_acid_1": structure.sequence[residue_index],
            "amino_acid_3": structure.residue_names[residue_index],
        }
        for metric, metric_index in metric_indices.items():
            row[metric] = float(mean[residue_index, metric_index])
            row[metric + "_std"] = float(std[residue_index, metric_index])
        rows.append(row)
    del example, batch, output
    return pd.DataFrame(rows)


def condition_specs(
    pdb_id: str,
    topology: Path,
    chain_roles: Mapping[str, str],
    prepared_dir: Path,
) -> List[Tuple[str, Path]]:
    ordered_chains = list(chain_roles)
    pmhc_chains = [chain for chain in ordered_chains if chain_roles[chain] in PMHC_ROLES]
    tcr_chains = [chain for chain in ordered_chains if chain_roles[chain] in TCR_ROLES]
    specs: List[Tuple[str, Path]] = [("full", topology)]
    subsets = [("pmhc", pmhc_chains), ("tcr", tcr_chains)]
    subsets.extend((f"isolated_{chain}", [chain]) for chain in ordered_chains)
    for condition, chains in subsets:
        destination = prepared_dir / pdb_id / f"{pdb_id}_{condition}.pdb"
        if not destination.exists():
            extract_chains(topology, destination, chains)
        specs.append((condition, destination))
    return specs


def run_predictions(
    manifest: pd.DataFrame,
    roles: Mapping[str, Mapping[str, str]],
    checkpoint: Path,
    output_dir: Path,
    device_name: str,
    temperature: float,
    window_ns: float,
    dataset_code: float,
    reuse: bool,
) -> Tuple[pd.DataFrame, pd.DataFrame, Sequence[str], str]:
    prepared_dir = output_dir / "prepared"
    prediction_dir = output_dir / "predictions"
    prepared_dir.mkdir(parents=True, exist_ok=True)
    prediction_dir.mkdir(parents=True, exist_ok=True)
    model, model_config, device = load_model(str(checkpoint), device_name)
    model.eval()
    missing_metrics = set(MODEL_METRICS) - set(model_config.metrics)
    if missing_metrics:
        raise ValueError(f"Checkpoint lacks required metrics: {sorted(missing_metrics)}")

    prediction_paths: List[Path] = []
    timing_rows: List[Dict[str, Any]] = []
    task_index = 0
    task_count = len(manifest) * 8
    for record in manifest.to_dict(orient="records"):
        pdb_id = str(record["pdb_id"]).lower()
        topology = Path(str(record["topology_path"]))
        chain_roles = roles[pdb_id]
        specs = condition_specs(pdb_id, topology, chain_roles, prepared_dir)
        for condition, pdb_path in specs:
            task_index += 1
            destination = prediction_dir / pdb_id / f"{condition}.csv"
            prediction_paths.append(destination)
            started = time.perf_counter()
            status, error = "ok", ""
            if reuse and destination.exists():
                frame = pd.read_csv(destination)
                required = set(MODEL_METRICS) | {
                    metric + "_std" for metric in MODEL_METRICS
                }
                if required.issubset(frame.columns):
                    status = "reused"
                else:
                    frame = predict_condition(
                        model, model_config, device, pdb_id, condition, pdb_path,
                        chain_roles, temperature, window_ns, dataset_code,
                    )
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    frame.to_csv(destination, index=False)
            else:
                try:
                    frame = predict_condition(
                        model, model_config, device, pdb_id, condition, pdb_path,
                        chain_roles, temperature, window_ns, dataset_code,
                    )
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    frame.to_csv(destination, index=False)
                except Exception as exception:
                    status, error = "error", repr(exception)
                    frame = pd.DataFrame()
            timing_rows.append(
                {
                    "pdb_id": pdb_id,
                    "condition": condition,
                    "residue_count": len(frame),
                    "end_to_end_seconds": time.perf_counter() - started,
                    "status": status,
                    "error": error,
                }
            )
            if task_index % 20 == 0 or task_index == task_count:
                failures = sum(row["status"] == "error" for row in timing_rows)
                print(
                    f"ITrace multi-metric inference {task_index}/{task_count} "
                    f"(failures={failures})",
                    flush=True,
                )
        gc.collect()
    failures = [row for row in timing_rows if row["status"] == "error"]
    if failures:
        raise RuntimeError(f"Inference failures: {failures[:5]!r}")
    predictions = pd.concat(
        [pd.read_csv(path, keep_default_na=False) for path in prediction_paths],
        ignore_index=True,
    )
    # Per-condition CSVs that contain only integer-looking author residue IDs
    # are otherwise inferred as integers, whereas a full-complex file with
    # insertion codes is inferred as strings.  Normalize before matching.
    for column in ("pdb_id", "condition", "chain", "author_residue_id"):
        predictions[column] = predictions[column].astype(str)
    timings = pd.DataFrame(timing_rows)
    return predictions, timings, list(model_config.metrics), str(device)


def add_prefix(frame: pd.DataFrame, prefix: str) -> pd.DataFrame:
    keys = ["pdb_id", "chain", "author_residue_id"]
    keep = keys + [metric for metric in MODEL_METRICS] + [
        metric + "_std" for metric in MODEL_METRICS
    ]
    output = frame[keep].copy()
    return output.rename(
        columns={column: prefix + column for column in keep if column not in keys}
    )


def build_matched_responses(predictions: pd.DataFrame) -> pd.DataFrame:
    keys = ["pdb_id", "chain", "author_residue_id"]
    full = predictions[predictions["condition"] == "full"].copy()
    full = full.rename(
        columns={
            metric: "full_" + metric
            for metric in MODEL_METRICS
        }
        | {
            metric + "_std": "full_" + metric + "_std"
            for metric in MODEL_METRICS
        }
    )

    reference_frames: List[pd.DataFrame] = []
    for condition in ("pmhc", "tcr"):
        frame = add_prefix(
            predictions[predictions["condition"] == condition], "reference_"
        )
        frame["reference_condition"] = condition
        reference_frames.append(frame)
    reference = pd.concat(reference_frames, ignore_index=True)

    isolated_frames: List[pd.DataFrame] = []
    for chain in sorted(predictions["chain"].unique()):
        condition = f"isolated_{chain}"
        frame = add_prefix(
            predictions[predictions["condition"] == condition], "isolated_"
        )
        isolated_frames.append(frame)
    isolated = pd.concat(isolated_frames, ignore_index=True)

    matched = full.merge(reference, on=keys, how="left", validate="one_to_one")
    matched = matched.merge(isolated, on=keys, how="left", validate="one_to_one")
    expected_reference = np.where(
        matched["chain_role"].isin(PMHC_ROLES), "pmhc", "tcr"
    )
    if not np.all(matched["reference_condition"].to_numpy() == expected_reference):
        raise ValueError("Biological reference condition does not match chain role")
    for metric in MODEL_METRICS:
        matched["binding_context_delta_" + metric] = (
            matched["full_" + metric] - matched["reference_" + metric]
        )
        matched["binding_context_combined_std_" + metric] = np.sqrt(
            matched["full_" + metric + "_std"] ** 2
            + matched["reference_" + metric + "_std"] ** 2
        )
        matched["binding_context_z_" + metric] = (
            matched["binding_context_delta_" + metric]
            / matched["binding_context_combined_std_" + metric]
        )
        matched["isolated_context_delta_" + metric] = (
            matched["full_" + metric] - matched["isolated_" + metric]
        )
    matched["chi_defined"] = ~matched["amino_acid_1"].isin(CHI_UNDEFINED_AA)
    matched["external_clean"] = ~matched["pdb_id"].isin(EXACT_TRAINING_OVERLAP)
    return matched


def common_annotations(context_dir: Path) -> pd.DataFrame:
    path = context_dir / "raw" / "common_context_predictions.csv"
    columns = [
        "pdb_id", "chain", "residue_number", "insertion_code", "region",
        "region_group", "is_tcr_pmhc_interface_8A",
        "target_chainfit_rmsf_angstrom",
    ]
    frame = pd.read_csv(path, usecols=columns, dtype=str, keep_default_na=False)
    frame["pdb_id"] = frame["pdb_id"].str.lower()
    frame["author_residue_id"] = (
        frame["residue_number"].str.strip() + frame["insertion_code"].str.strip()
    )
    frame["is_tcr_pmhc_interface_8A"] = (
        frame["is_tcr_pmhc_interface_8A"].str.lower() == "true"
    )
    frame["target_chainfit_rmsf_angstrom"] = pd.to_numeric(
        frame["target_chainfit_rmsf_angstrom"], errors="coerce"
    )
    return frame.drop(columns=["residue_number", "insertion_code"])


def atom_element(line: str) -> str:
    element = line[76:78].strip().upper()
    if element:
        return element
    atom_name = line[12:16].strip().upper()
    return re.sub(r"^[0-9]+", "", atom_name)[:1]


def heavy_atom_residues(
    pdb_path: Path,
) -> Dict[Tuple[str, str], Dict[str, Any]]:
    residues: MutableMapping[Tuple[str, str], Dict[str, Any]] = {}
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM  "):
            continue
        if line[16:17] not in {" ", "A"}:
            continue
        if atom_element(line) == "H":
            continue
        chain = line[21:22].strip()
        author_id = line[22:26].strip() + line[26:27].strip()
        key = (chain, author_id)
        residue = residues.setdefault(
            key,
            {
                "chain": chain,
                "author_residue_id": author_id,
                "amino_acid_3": line[17:20].strip(),
                "coords": [],
            },
        )
        residue["coords"].append(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
    for residue in residues.values():
        residue["coords"] = np.asarray(residue["coords"], dtype=np.float64)
    return dict(residues)


class TargetAtoms:
    def __init__(
        self,
        residues: Mapping[Tuple[str, str], Mapping[str, Any]],
        keys: Iterable[Tuple[str, str]],
    ) -> None:
        coordinates: List[np.ndarray] = []
        atom_keys: List[Tuple[str, str]] = []
        self.residues = residues
        for key in keys:
            if key not in residues:
                continue
            coords = np.asarray(residues[key]["coords"], dtype=np.float64)
            coordinates.append(coords)
            atom_keys.extend([key] * len(coords))
        self.coords = np.concatenate(coordinates, axis=0) if coordinates else np.empty((0, 3))
        self.atom_keys = atom_keys
        self.tree = cKDTree(self.coords) if len(self.coords) else None

    def describe(
        self, source_coords: np.ndarray, cutoff: float
    ) -> Tuple[float, Tuple[str, str] | None, int]:
        if self.tree is None or not len(source_coords):
            return math.nan, None, 0
        distances, indices = self.tree.query(source_coords, k=1)
        nearest_atom = int(indices[int(np.argmin(distances))])
        nearest_key = self.atom_keys[nearest_atom]
        contacted: set[Tuple[str, str]] = set()
        for atom_indices in self.tree.query_ball_point(source_coords, cutoff):
            contacted.update(self.atom_keys[int(index)] for index in atom_indices)
        return float(np.min(distances)), nearest_key, len(contacted)


def contact_annotations_for_complex(
    pdb_id: str,
    topology: Path,
    chain_roles: Mapping[str, str],
) -> pd.DataFrame:
    residues = heavy_atom_residues(topology)
    by_role: Dict[str, List[Tuple[str, str]]] = {}
    for key in residues:
        by_role.setdefault(chain_roles[key[0]], []).append(key)
    targets = {
        "tcr": TargetAtoms(
            residues, by_role.get("tcr_alpha", []) + by_role.get("tcr_beta", [])
        ),
        "tcr_alpha": TargetAtoms(residues, by_role.get("tcr_alpha", [])),
        "tcr_beta": TargetAtoms(residues, by_role.get("tcr_beta", [])),
        "peptide": TargetAtoms(residues, by_role.get("peptide", [])),
        "mhc": TargetAtoms(residues, by_role.get("mhc", []) + by_role.get("b2m", [])),
    }
    rows: List[Dict[str, Any]] = []
    for key, residue in residues.items():
        role = chain_roles[key[0]]
        row: Dict[str, Any] = {
            "pdb_id": pdb_id,
            "chain": key[0],
            "author_residue_id": key[1],
        }
        requested: Sequence[str]
        if role == "peptide":
            requested = ("tcr", "tcr_alpha", "tcr_beta")
        elif role in TCR_ROLES:
            requested = ("peptide", "mhc")
        else:
            requested = ("tcr",)
        for target_name in requested:
            distance, nearest_key, contact_count = targets[target_name].describe(
                residue["coords"], CONTACT_CUTOFF_ANGSTROM
            )
            row[f"min_distance_to_{target_name}_angstrom"] = distance
            row[f"{target_name}_contact_residue_count_5A"] = contact_count
            if nearest_key is not None:
                nearest = residues[nearest_key]
                row[f"nearest_{target_name}_chain"] = nearest_key[0]
                row[f"nearest_{target_name}_author_residue_id"] = nearest_key[1]
                row[f"nearest_{target_name}_amino_acid_3"] = nearest["amino_acid_3"]
        rows.append(row)
    return pd.DataFrame(rows)


def add_peptide_positions(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    output["peptide_length"] = np.nan
    output["peptide_position_n"] = np.nan
    output["peptide_position_c"] = np.nan
    output["peptide_position_label"] = ""
    peptide = output[output["chain_role"] == "peptide"]
    for (_, _), indices in peptide.groupby(["pdb_id", "chain"], sort=False).groups.items():
        ordered = output.loc[list(indices)].sort_values("residue_index").index.tolist()
        length = len(ordered)
        for position, index in enumerate(ordered, start=1):
            from_c = length - position + 1
            label = f"P{position}" if position <= 3 else f"Pomega-{from_c - 1}"
            if from_c == 1:
                label = "Pomega"
            elif from_c == 2:
                label = "Pomega-1"
            output.at[index, "peptide_length"] = length
            output.at[index, "peptide_position_n"] = position
            output.at[index, "peptide_position_c"] = from_c
            output.at[index, "peptide_position_label"] = label
    for column in ("peptide_length", "peptide_position_n", "peptide_position_c"):
        output[column] = output[column].astype("Int64")
    return output


def cluster_bootstrap_mean(
    frame: pd.DataFrame,
    column: str,
    repetitions: int,
    seed: int,
) -> Tuple[float, float]:
    selected = frame[["pdb_id", column]].replace([np.inf, -np.inf], np.nan).dropna()
    if selected.empty:
        return math.nan, math.nan
    grouped = selected.groupby("pdb_id")[column].agg(["sum", "count"])
    sums = grouped["sum"].to_numpy(dtype=float)
    counts = grouped["count"].to_numpy(dtype=float)
    if len(grouped) == 1:
        value = float(sums[0] / counts[0])
        return value, value
    rng = np.random.default_rng(seed)
    values = np.empty(repetitions, dtype=float)
    for index in range(repetitions):
        sampled = rng.integers(0, len(grouped), size=len(grouped))
        values[index] = sums[sampled].sum() / counts[sampled].sum()
    lower, upper = np.percentile(values, [2.5, 97.5])
    return float(lower), float(upper)


def bootstrap_mean_array(
    values: Sequence[float],
    repetitions: int,
    seed: int,
) -> Tuple[float, float]:
    """Bootstrap a mean when each value already represents one complex."""

    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not len(finite):
        return math.nan, math.nan
    if len(finite) == 1:
        return float(finite[0]), float(finite[0])
    rng = np.random.default_rng(seed)
    sampled = rng.choice(finite, size=(repetitions, len(finite)), replace=True)
    lower, upper = np.percentile(sampled.mean(axis=1), [2.5, 97.5])
    return float(lower), float(upper)


def summarize_md_peptide_positions(
    peptide: pd.DataFrame,
    repetitions: int,
    seed: int,
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for (length, position), group in peptide.groupby(
        ["peptide_length", "peptide_position_n"], sort=True
    ):
        raw = group["target_chainfit_rmsf_angstrom"].dropna().to_numpy(dtype=float)
        z_values = group["md_rmsf_z_within_peptide"].dropna().to_numpy(dtype=float)
        raw_low, raw_high = cluster_bootstrap_mean(
            group, "target_chainfit_rmsf_angstrom", repetitions, seed
        )
        z_low, z_high = cluster_bootstrap_mean(
            group, "md_rmsf_z_within_peptide", repetitions, seed + 1
        )
        rows.append(
            {
                "peptide_length": int(length),
                "peptide_position_n": int(position),
                "complex_count": int(group["pdb_id"].nunique()),
                "mean_md_rmsf_angstrom": float(np.mean(raw)),
                "median_md_rmsf_angstrom": float(np.median(raw)),
                "cluster_bootstrap_ci_low_md_rmsf_angstrom": raw_low,
                "cluster_bootstrap_ci_high_md_rmsf_angstrom": raw_high,
                "mean_md_rmsf_z_within_peptide": float(np.mean(z_values)),
                "median_md_rmsf_z_within_peptide": float(np.median(z_values)),
                "cluster_bootstrap_ci_low_md_rmsf_z": z_low,
                "cluster_bootstrap_ci_high_md_rmsf_z": z_high,
            }
        )
    return pd.DataFrame(rows)


def summarize_group(
    frame: pd.DataFrame,
    group_fields: Sequence[str],
    delta_prefix: str,
    repetitions: int,
    seed: int,
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    grouped = frame.groupby(list(group_fields), dropna=False, sort=True)
    for group_key, group in grouped:
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        row: Dict[str, Any] = dict(zip(group_fields, group_key))
        row["residue_count"] = int(len(group))
        row["complex_count"] = int(group["pdb_id"].nunique())
        for metric_index, metric in enumerate(MODEL_METRICS):
            selected = group
            if metric == "chi_entropy":
                selected = selected[selected["chi_defined"]]
            column = delta_prefix + metric
            values = selected[column].to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            lower, upper = cluster_bootstrap_mean(
                selected, column, repetitions, seed + metric_index
            )
            row[f"{metric}_valid_residue_count"] = int(len(values))
            row[f"mean_delta_{metric}"] = float(np.mean(values)) if len(values) else math.nan
            row[f"median_delta_{metric}"] = float(np.median(values)) if len(values) else math.nan
            row[f"fraction_negative_delta_{metric}"] = (
                float(np.mean(values < 0)) if len(values) else math.nan
            )
            row[f"cluster_bootstrap_ci_low_{metric}"] = lower
            row[f"cluster_bootstrap_ci_high_{metric}"] = upper
        rows.append(row)
    return pd.DataFrame(rows)


def spearman_summary(frame: pd.DataFrame, distance_column: str) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for metric in MODEL_METRICS:
        selected = frame
        if metric == "chi_entropy":
            selected = selected[selected["chi_defined"]]
        pair = selected[[distance_column, "binding_context_delta_" + metric]].dropna()
        statistic, pvalue = spearmanr(pair.iloc[:, 0], pair.iloc[:, 1])
        result[metric] = {
            "residue_count": int(len(pair)),
            "spearman_rho": float(statistic),
            "pvalue_naive_residue_level": float(pvalue),
            "note": "P value is descriptive only because residues are nested within complexes.",
        }
    return result


def within_complex_spearman_summary(
    frame: pd.DataFrame,
    distance_column: str,
    repetitions: int,
    seed: int,
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    """Summarize distance-response association with complexes as replicates."""

    result: Dict[str, Any] = {}
    rows: List[Dict[str, Any]] = []
    for metric_index, metric in enumerate(MODEL_METRICS):
        selected = frame
        if metric == "chi_entropy":
            selected = selected[selected["chi_defined"]]
        column = "binding_context_delta_" + metric
        correlations: List[float] = []
        for pdb_id, group in selected.groupby("pdb_id", sort=True):
            pair = group[[distance_column, column]].dropna()
            if len(pair) < 4 or pair.iloc[:, 0].nunique() < 2 or pair.iloc[:, 1].nunique() < 2:
                continue
            statistic = float(spearmanr(pair.iloc[:, 0], pair.iloc[:, 1]).statistic)
            if np.isfinite(statistic):
                correlations.append(statistic)
                rows.append(
                    {
                        "pdb_id": pdb_id,
                        "metric": metric,
                        "residue_count": int(len(pair)),
                        "spearman_rho": statistic,
                    }
                )
        values = np.asarray(correlations, dtype=float)
        rng = np.random.default_rng(seed + metric_index)
        if len(values):
            bootstrap = np.asarray(
                [
                    np.mean(rng.choice(values, size=len(values), replace=True))
                    for _ in range(repetitions)
                ],
                dtype=float,
            )
            lower, upper = np.percentile(bootstrap, [2.5, 97.5])
        else:
            lower, upper = math.nan, math.nan
        result[metric] = {
            "complex_count": int(len(values)),
            "mean_within_complex_spearman_rho": (
                float(np.mean(values)) if len(values) else math.nan
            ),
            "median_within_complex_spearman_rho": (
                float(np.median(values)) if len(values) else math.nan
            ),
            "complex_bootstrap_ci_low": float(lower),
            "complex_bootstrap_ci_high": float(upper),
        }
    return result, pd.DataFrame(rows)


def paired_contact_effect_summary(
    peptide: pd.DataFrame,
    repetitions: int,
    seed: int,
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    """Within each complex compare contacted and non-contacted peptide residues."""

    rows: List[Dict[str, Any]] = []
    result: Dict[str, Any] = {}
    for metric_index, metric in enumerate(MODEL_METRICS):
        selected = peptide
        if metric == "chi_entropy":
            selected = selected[selected["chi_defined"]]
        column = "binding_context_delta_" + metric
        differences: List[float] = []
        for pdb_id, group in selected.groupby("pdb_id", sort=True):
            contacted = group.loc[group["direct_tcr_contact_5A"], column].dropna()
            noncontacted = group.loc[~group["direct_tcr_contact_5A"], column].dropna()
            if contacted.empty or noncontacted.empty:
                continue
            difference = float(contacted.mean() - noncontacted.mean())
            differences.append(difference)
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "metric": metric,
                    "contacted_residue_count": int(len(contacted)),
                    "noncontacted_residue_count": int(len(noncontacted)),
                    "mean_contacted_minus_noncontacted_delta": difference,
                }
            )
        values = np.asarray(differences, dtype=float)
        rng = np.random.default_rng(seed + metric_index)
        if len(values):
            bootstrap = np.asarray(
                [
                    np.mean(rng.choice(values, size=len(values), replace=True))
                    for _ in range(repetitions)
                ],
                dtype=float,
            )
            lower, upper = np.percentile(bootstrap, [2.5, 97.5])
        else:
            lower, upper = math.nan, math.nan
        result[metric] = {
            "paired_complex_count": int(len(values)),
            "mean_contacted_minus_noncontacted_delta": (
                float(np.mean(values)) if len(values) else math.nan
            ),
            "median_contacted_minus_noncontacted_delta": (
                float(np.median(values)) if len(values) else math.nan
            ),
            "fraction_more_negative_at_contact": (
                float(np.mean(values < 0)) if len(values) else math.nan
            ),
            "complex_bootstrap_ci_low": float(lower),
            "complex_bootstrap_ci_high": float(upper),
        }
    return result, pd.DataFrame(rows)


def rmsf_prediction_fidelity_summary(
    peptide: pd.DataFrame,
    repetitions: int,
    seed: int,
) -> Tuple[Dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Test whether complete-complex predictions recover peptide MD patterns.

    The context deltas answer what the model changes when TCR is removed.  This
    function asks the separate and more important validation question: does the
    complete-complex RMSF prediction agree with the held-out ITrace MD RMSF?
    Complexes, rather than residues, are treated as independent replicates.
    """

    selected = peptide.copy()
    selected["model_rmsf_z_within_peptide"] = (
        selected.groupby(["pdb_id", "chain"])["full_rmsf"].transform(
            lambda values: (
                (values - values.mean()) / values.std(ddof=1)
                if values.std(ddof=1) > 0
                else np.nan
            )
        )
    )
    pair = selected[["full_rmsf", "target_chainfit_rmsf_angstrom"]].dropna()
    pooled_rho, pooled_pvalue = spearmanr(pair.iloc[:, 0], pair.iloc[:, 1])

    complex_rows: List[Dict[str, Any]] = []
    contact_rows: List[Dict[str, Any]] = []
    for pdb_id, group in selected.groupby("pdb_id", sort=True):
        valid = group[[
            "peptide_position_n", "full_rmsf", "target_chainfit_rmsf_angstrom"
        ]].dropna()
        if len(valid) < 3:
            continue
        rho = math.nan
        if valid["full_rmsf"].nunique() > 1 and valid["target_chainfit_rmsf_angstrom"].nunique() > 1:
            rho = float(spearmanr(
                valid["full_rmsf"], valid["target_chainfit_rmsf_angstrom"]
            ).statistic)
        md_top = int(valid.loc[valid["target_chainfit_rmsf_angstrom"].idxmax(), "peptide_position_n"])
        model_top = int(valid.loc[valid["full_rmsf"].idxmax(), "peptide_position_n"])
        md_bottom = int(valid.loc[valid["target_chainfit_rmsf_angstrom"].idxmin(), "peptide_position_n"])
        model_bottom = int(valid.loc[valid["full_rmsf"].idxmin(), "peptide_position_n"])
        complex_rows.append(
            {
                "pdb_id": pdb_id,
                "peptide_length": int(group["peptide_length"].iloc[0]),
                "residue_count": int(len(valid)),
                "spearman_rho": rho,
                "md_most_mobile_position": md_top,
                "model_most_mobile_position": model_top,
                "most_mobile_exact_match": bool(md_top == model_top),
                "most_mobile_within_one_position": bool(abs(md_top - model_top) <= 1),
                "md_least_mobile_position": md_bottom,
                "model_least_mobile_position": model_bottom,
                "least_mobile_exact_match": bool(md_bottom == model_bottom),
                "least_mobile_within_one_position": bool(abs(md_bottom - model_bottom) <= 1),
            }
        )

        contacted = group["direct_tcr_contact_5A"]
        if contacted.any() and (~contacted).any():
            md_difference = float(
                group.loc[contacted, "target_chainfit_rmsf_angstrom"].mean()
                - group.loc[~contacted, "target_chainfit_rmsf_angstrom"].mean()
            )
            model_difference = float(
                group.loc[contacted, "full_rmsf"].mean()
                - group.loc[~contacted, "full_rmsf"].mean()
            )
            contact_rows.append(
                {
                    "pdb_id": pdb_id,
                    "contacted_residue_count": int(contacted.sum()),
                    "noncontacted_residue_count": int((~contacted).sum()),
                    "md_contacted_minus_noncontacted_rmsf_angstrom": md_difference,
                    "model_contacted_minus_noncontacted_rmsf_angstrom": model_difference,
                    "same_direction": bool(np.sign(md_difference) == np.sign(model_difference)),
                }
            )

    complex_frame = pd.DataFrame(complex_rows)
    contact_frame = pd.DataFrame(contact_rows)
    rho_values = complex_frame["spearman_rho"].dropna().to_numpy(dtype=float)
    rho_low, rho_high = bootstrap_mean_array(rho_values, repetitions, seed)
    positive_low, positive_high = bootstrap_mean_array(
        (rho_values > 0).astype(float), repetitions, seed + 1
    )

    nine_mer = selected[selected["peptide_length"] == 9].copy()
    position_rows: List[Dict[str, Any]] = []
    for position, group in nine_mer.groupby("peptide_position_n", sort=True):
        md_low, md_high = cluster_bootstrap_mean(
            group, "md_rmsf_z_within_peptide", repetitions, seed + 2
        )
        model_low, model_high = cluster_bootstrap_mean(
            group, "model_rmsf_z_within_peptide", repetitions, seed + 3
        )
        position_rows.append(
            {
                "peptide_position_n": int(position),
                "complex_count": int(group["pdb_id"].nunique()),
                "mean_md_rmsf_z_within_peptide": float(group["md_rmsf_z_within_peptide"].mean()),
                "cluster_bootstrap_ci_low_md_rmsf_z": md_low,
                "cluster_bootstrap_ci_high_md_rmsf_z": md_high,
                "mean_model_rmsf_z_within_peptide": float(group["model_rmsf_z_within_peptide"].mean()),
                "cluster_bootstrap_ci_low_model_rmsf_z": model_low,
                "cluster_bootstrap_ci_high_model_rmsf_z": model_high,
                "mean_md_rmsf_angstrom": float(group["target_chainfit_rmsf_angstrom"].mean()),
                "mean_model_rmsf_angstrom": float(group["full_rmsf"].mean()),
            }
        )
    position_frame = pd.DataFrame(position_rows)
    position_spearman = spearmanr(
        position_frame["mean_model_rmsf_z_within_peptide"],
        position_frame["mean_md_rmsf_z_within_peptide"],
    )
    position_pearson = pearsonr(
        position_frame["mean_model_rmsf_z_within_peptide"],
        position_frame["mean_md_rmsf_z_within_peptide"],
    )

    md_contact = contact_frame["md_contacted_minus_noncontacted_rmsf_angstrom"].to_numpy(dtype=float)
    model_contact = contact_frame["model_contacted_minus_noncontacted_rmsf_angstrom"].to_numpy(dtype=float)
    md_contact_low, md_contact_high = bootstrap_mean_array(
        md_contact, repetitions, seed + 4
    )
    model_contact_low, model_contact_high = bootstrap_mean_array(
        model_contact, repetitions, seed + 5
    )
    contact_rho, contact_pvalue = spearmanr(md_contact, model_contact)

    result = {
        "scope": (
            "Complete-complex FlexPPI RMSF predictions versus ITrace MD RMSF on "
            "external-clean peptide residues; no context delta is used."
        ),
        "residue_count": int(len(pair)),
        "complex_count": int(complex_frame["pdb_id"].nunique()),
        "pooled_residue_spearman_rho": float(pooled_rho),
        "pooled_residue_pvalue_descriptive_only": float(pooled_pvalue),
        "peptide_rmsf_mae_angstrom": float(
            np.mean(np.abs(pair["full_rmsf"] - pair["target_chainfit_rmsf_angstrom"]))
        ),
        "peptide_rmsf_mean_error_angstrom": float(
            np.mean(pair["full_rmsf"] - pair["target_chainfit_rmsf_angstrom"])
        ),
        "within_complex_ranking": {
            "mean_spearman_rho": float(np.mean(rho_values)),
            "median_spearman_rho": float(np.median(rho_values)),
            "complex_bootstrap_ci_low": rho_low,
            "complex_bootstrap_ci_high": rho_high,
            "fraction_positive": float(np.mean(rho_values > 0)),
            "fraction_positive_bootstrap_ci_low": positive_low,
            "fraction_positive_bootstrap_ci_high": positive_high,
            "fraction_spearman_at_least_0_5": float(np.mean(rho_values >= 0.5)),
            "most_mobile_exact_match_fraction": float(complex_frame["most_mobile_exact_match"].mean()),
            "most_mobile_within_one_position_fraction": float(complex_frame["most_mobile_within_one_position"].mean()),
            "least_mobile_exact_match_fraction": float(complex_frame["least_mobile_exact_match"].mean()),
            "least_mobile_within_one_position_fraction": float(complex_frame["least_mobile_within_one_position"].mean()),
        },
        "nine_mer_mean_position_pattern": {
            "position_count": int(len(position_frame)),
            "complex_count": int(nine_mer["pdb_id"].nunique()),
            "pearson_r": float(position_pearson.statistic),
            "pearson_pvalue": float(position_pearson.pvalue),
            "spearman_rho": float(position_spearman.statistic),
            "spearman_pvalue": float(position_spearman.pvalue),
            "mean_absolute_z_error": float(np.mean(np.abs(
                position_frame["mean_model_rmsf_z_within_peptide"]
                - position_frame["mean_md_rmsf_z_within_peptide"]
            ))),
        },
        "direct_tcr_contact_pattern": {
            "paired_complex_count": int(len(contact_frame)),
            "mean_md_contacted_minus_noncontacted_rmsf_angstrom": float(np.mean(md_contact)),
            "md_complex_bootstrap_ci_low": md_contact_low,
            "md_complex_bootstrap_ci_high": md_contact_high,
            "fraction_md_contacted_lower": float(np.mean(md_contact < 0)),
            "mean_model_contacted_minus_noncontacted_rmsf_angstrom": float(np.mean(model_contact)),
            "model_complex_bootstrap_ci_low": model_contact_low,
            "model_complex_bootstrap_ci_high": model_contact_high,
            "fraction_model_contacted_lower": float(np.mean(model_contact < 0)),
            "across_complex_spearman_rho": float(contact_rho),
            "across_complex_spearman_pvalue": float(contact_pvalue),
            "same_direction_fraction": float(contact_frame["same_direction"].mean()),
        },
    }
    return result, position_frame, complex_frame, contact_frame


def validate_full_rmsf(
    matched: pd.DataFrame, benchmark_dir: Path
) -> Dict[str, Any]:
    path = benchmark_dir / "raw" / "flexppi_full_complex_predictions.csv"
    old = pd.read_csv(path, dtype={"residue_number": str}, keep_default_na=False)
    old["pdb_id"] = old["pdb_id"].str.lower()
    old["author_residue_id"] = (
        old["residue_number"].astype(str).str.strip()
        + old["insertion_code"].astype(str).str.strip()
    )
    old = old[["pdb_id", "chain", "author_residue_id", "flexppi_rmsf_angstrom"]]
    check = matched.merge(
        old, on=["pdb_id", "chain", "author_residue_id"], how="left", validate="one_to_one"
    )
    difference = np.abs(
        check["full_rmsf"].to_numpy(dtype=float)
        - check["flexppi_rmsf_angstrom"].to_numpy(dtype=float)
    )
    finite = difference[np.isfinite(difference)]
    return {
        "matched_residue_count": int(len(finite)),
        "missing_old_prediction_count": int(np.sum(~np.isfinite(difference))),
        "maximum_absolute_difference_angstrom": float(np.max(finite)),
        "mean_absolute_difference_angstrom": float(np.mean(finite)),
        "passes_1e-6_angstrom": bool(np.max(finite) < 1e-6),
    }


def write_analysis(
    matched: pd.DataFrame,
    analysis_dir: Path,
    benchmark_dir: Path,
    repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    analysis_dir.mkdir(parents=True, exist_ok=True)
    clean = matched[matched["external_clean"]].copy()
    peptide = clean[clean["chain_role"] == "peptide"].copy()
    tcr = clean[clean["chain_role"].isin(TCR_ROLES)].copy()
    pmhc = clean[clean["chain_role"].isin({"mhc", "b2m"})].copy()
    peptide.to_csv(analysis_dir / "peptide_residue_context.csv", index=False)
    tcr.to_csv(analysis_dir / "tcr_residue_context.csv", index=False)
    pmhc.to_csv(analysis_dir / "pmhc_residue_context.csv", index=False)

    position = summarize_group(
        peptide,
        ["peptide_length", "peptide_position_n", "amino_acid_1"],
        "binding_context_delta_",
        repetitions,
        seed,
    )
    position.to_csv(analysis_dir / "peptide_length_position_amino_acid_summary.csv", index=False)
    position_only = summarize_group(
        peptide,
        ["peptide_length", "peptide_position_n"],
        "binding_context_delta_",
        repetitions,
        seed + 10,
    )
    position_only.to_csv(analysis_dir / "peptide_length_position_summary.csv", index=False)
    nine_mer = position_only[position_only["peptide_length"] == 9].copy()
    nine_mer.to_csv(analysis_dir / "peptide_9mer_position_summary.csv", index=False)
    md_position = summarize_md_peptide_positions(
        peptide, repetitions, seed + 15
    )
    md_position.to_csv(
        analysis_dir / "peptide_md_rmsf_position_summary.csv", index=False
    )
    md_position[md_position["peptide_length"] == 9].to_csv(
        analysis_dir / "peptide_9mer_md_rmsf_position_summary.csv", index=False
    )
    amino_acid = summarize_group(
        peptide,
        ["amino_acid_1", "amino_acid_3"],
        "binding_context_delta_",
        repetitions,
        seed + 20,
    )
    amino_acid.to_csv(analysis_dir / "peptide_amino_acid_summary.csv", index=False)
    # Remove the mean response of each length/position stratum before asking
    # whether an amino-acid type has an additional response.  This prevents,
    # for example, enrichment of Leu at P2 from masquerading as a Leu effect.
    for metric in MODEL_METRICS:
        source = "binding_context_delta_" + metric
        destination = "position_adjusted_binding_context_delta_" + metric
        if metric == "chi_entropy":
            peptide[destination] = np.nan
            valid = peptide["chi_defined"]
            peptide.loc[valid, destination] = (
                peptide.loc[valid, source]
                - peptide.loc[valid]
                .groupby(["peptide_length", "peptide_position_n"])[source]
                .transform("mean")
            )
        else:
            peptide[destination] = (
                peptide[source]
                - peptide.groupby(["peptide_length", "peptide_position_n"])[source]
                .transform("mean")
            )
    amino_acid_adjusted = summarize_group(
        peptide,
        ["amino_acid_1", "amino_acid_3"],
        "position_adjusted_binding_context_delta_",
        repetitions,
        seed + 25,
    )
    amino_acid_adjusted.to_csv(
        analysis_dir / "peptide_amino_acid_position_adjusted_summary.csv", index=False
    )
    peptide["position_adjusted_md_rmsf_angstrom"] = (
        peptide["target_chainfit_rmsf_angstrom"]
        - peptide.groupby(["peptide_length", "peptide_position_n"])[
            "target_chainfit_rmsf_angstrom"
        ].transform("mean")
    )
    md_amino_rows: List[Dict[str, Any]] = []
    for (aa1, aa3), group in peptide.groupby(
        ["amino_acid_1", "amino_acid_3"], sort=True
    ):
        values = group["position_adjusted_md_rmsf_angstrom"].dropna().to_numpy(dtype=float)
        low, high = cluster_bootstrap_mean(
            group,
            "position_adjusted_md_rmsf_angstrom",
            repetitions,
            seed + 27,
        )
        md_amino_rows.append(
            {
                "amino_acid_1": aa1,
                "amino_acid_3": aa3,
                "residue_count": int(len(values)),
                "complex_count": int(group["pdb_id"].nunique()),
                "mean_position_adjusted_md_rmsf_angstrom": float(np.mean(values)),
                "median_position_adjusted_md_rmsf_angstrom": float(np.median(values)),
                "cluster_bootstrap_ci_low": low,
                "cluster_bootstrap_ci_high": high,
            }
        )
    pd.DataFrame(md_amino_rows).to_csv(
        analysis_dir / "peptide_md_rmsf_amino_acid_position_adjusted_summary.csv",
        index=False,
    )

    peptide["direct_tcr_contact_5A"] = (
        peptide["min_distance_to_tcr_angstrom"] <= CONTACT_CUTOFF_ANGSTROM
    )
    prediction_fidelity, prediction_position, prediction_complex, prediction_contact = (
        rmsf_prediction_fidelity_summary(
            peptide, repetitions, seed + 28
        )
    )
    prediction_position.to_csv(
        analysis_dir / "peptide_9mer_rmsf_prediction_fidelity.csv", index=False
    )
    prediction_complex.to_csv(
        analysis_dir / "peptide_rmsf_prediction_within_complex.csv", index=False
    )
    prediction_contact.to_csv(
        analysis_dir / "peptide_rmsf_contact_prediction_validation.csv", index=False
    )
    contact = summarize_group(
        peptide,
        ["direct_tcr_contact_5A"],
        "binding_context_delta_",
        repetitions,
        seed + 30,
    )
    contact.to_csv(analysis_dir / "peptide_contact_summary.csv", index=False)
    paired_contact, paired_contact_rows = paired_contact_effect_summary(
        peptide, repetitions, seed + 35
    )
    paired_contact_rows.to_csv(
        analysis_dir / "peptide_contact_within_complex_differences.csv", index=False
    )
    within_distance, within_distance_rows = within_complex_spearman_summary(
        peptide,
        "min_distance_to_tcr_angstrom",
        repetitions,
        seed + 37,
    )
    within_distance_rows.to_csv(
        analysis_dir / "peptide_distance_within_complex_spearman.csv", index=False
    )
    isolated = summarize_group(
        clean,
        ["chain_role"],
        "isolated_context_delta_",
        repetitions,
        seed + 40,
    )
    isolated.to_csv(analysis_dir / "full_vs_isolated_multimetric_summary.csv", index=False)
    biological = summarize_group(
        clean,
        ["chain_role", "region_group"],
        "binding_context_delta_",
        repetitions,
        seed + 50,
    )
    biological.to_csv(analysis_dir / "binding_context_region_summary.csv", index=False)

    extreme_rows: List[pd.DataFrame] = []
    identity = [
        "pdb_id", "chain", "chain_role", "region_group", "author_residue_id",
        "amino_acid_1", "peptide_length", "peptide_position_n",
        "min_distance_to_tcr_angstrom", "min_distance_to_peptide_angstrom",
    ]
    for metric in MODEL_METRICS:
        selected = clean
        if metric == "chi_entropy":
            selected = selected[selected["chi_defined"]]
        column = "binding_context_delta_" + metric
        for direction, ascending in (("most_negative", True), ("most_positive", False)):
            columns = [field for field in identity if field in selected.columns] + [column]
            top = selected.nsmallest(50, column) if ascending else selected.nlargest(50, column)
            top = top[columns].copy()
            top.insert(0, "metric", metric)
            top.insert(1, "direction", direction)
            extreme_rows.append(top)
    pd.concat(extreme_rows, ignore_index=True).to_csv(
        analysis_dir / "binding_context_residue_extremes.csv", index=False
    )

    validation = validate_full_rmsf(matched, benchmark_dir)
    result = {
        "complex_count_external_clean": int(clean["pdb_id"].nunique()),
        "residue_count_external_clean": int(len(clean)),
        "peptide_residue_count": int(len(peptide)),
        "peptide_length_distribution": {
            str(key): int(value)
            for key, value in peptide.groupby("pdb_id")["peptide_length"].first().value_counts().sort_index().items()
        },
        "chi_entropy_interpretation": {
            "definition": "Natural-log Shannon entropy over predicted side-chain rotamer states.",
            "negative_delta": "Reduced side-chain conformational diversity in the complete binding context.",
            "excluded_amino_acids": sorted(CHI_UNDEFINED_AA),
            "reason": "Gly and Ala have no defined side-chain chi torsion in the labels.",
        },
        "binding_context_definition": (
            "complete minus pMHC for MHC/b2m/peptide, and complete minus TCR heterodimer "
            "for TCR residues, all at fixed bound coordinates"
        ),
        "peptide_rmsf_prediction_fidelity": prediction_fidelity,
        "peptide_distance_response_spearman": spearman_summary(
            peptide, "min_distance_to_tcr_angstrom"
        ),
        "peptide_distance_response_within_complex": within_distance,
        "peptide_contact_paired_within_complex": paired_contact,
        "rmsf_reproduction_check": validation,
    }
    return result


def main() -> None:
    args = parse_arguments()
    benchmark_dir = Path(args.benchmark_dir)
    context_dir = Path(args.context_dir)
    checkpoint = Path(args.checkpoint)
    output_dir = Path(args.output_dir)
    raw_dir = output_dir / "raw"
    analysis_dir = output_dir / "analysis"
    raw_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(benchmark_dir / "raw" / "itrace_complex_manifest.csv")
    manifest["pdb_id"] = manifest["pdb_id"].astype(str).str.lower()
    if args.limit is not None:
        manifest = manifest.head(args.limit).copy()
    roles = json.loads(
        (benchmark_dir / "raw" / "itrace_chain_roles.json").read_text(encoding="utf-8")
    )
    roles = {str(key).lower(): value for key, value in roles.items()}

    predictions, timings, checkpoint_metrics, device = run_predictions(
        manifest,
        roles,
        checkpoint,
        output_dir,
        args.device,
        args.temperature,
        args.window_ns,
        args.dataset_code,
        args.reuse,
    )
    predictions.to_csv(raw_dir / "all_condition_predictions.csv", index=False)
    timings.to_csv(raw_dir / "condition_timings.csv", index=False)
    matched = build_matched_responses(predictions)

    annotations = common_annotations(context_dir)
    matched = matched.merge(
        annotations,
        on=["pdb_id", "chain", "author_residue_id"],
        how="left",
        validate="one_to_one",
    )
    contact_frames: List[pd.DataFrame] = []
    for index, record in enumerate(manifest.to_dict(orient="records"), start=1):
        pdb_id = str(record["pdb_id"]).lower()
        contact_frames.append(
            contact_annotations_for_complex(
                pdb_id, Path(str(record["topology_path"])), roles[pdb_id]
            )
        )
        if index % 25 == 0 or index == len(manifest):
            print(f"Contact annotation {index}/{len(manifest)}", flush=True)
    contacts = pd.concat(contact_frames, ignore_index=True)
    matched = matched.merge(
        contacts,
        on=["pdb_id", "chain", "author_residue_id"],
        how="left",
        validate="one_to_one",
    )
    matched = add_peptide_positions(matched)
    matched["md_rmsf_z_within_peptide"] = np.nan
    peptide_mask = matched["chain_role"] == "peptide"
    matched.loc[peptide_mask, "md_rmsf_z_within_peptide"] = (
        matched.loc[peptide_mask]
        .groupby(["pdb_id", "chain"])["target_chainfit_rmsf_angstrom"]
        .transform(
            lambda values: (
                (values - values.mean()) / values.std(ddof=1)
                if values.std(ddof=1) > 0
                else np.nan
            )
        )
    )
    matched.to_csv(raw_dir / "matched_residue_context_responses.csv", index=False)

    analysis = write_analysis(
        matched,
        analysis_dir,
        benchmark_dir,
        args.bootstrap_repetitions,
        args.seed,
    )
    report = {
        "protocol": {
            "checkpoint": str(checkpoint.resolve()),
            "checkpoint_sha256": sha256(checkpoint),
            "checkpoint_metrics": checkpoint_metrics,
            "device": device,
            "temperature_kelvin": args.temperature,
            "window_ns": args.window_ns,
            "dataset_code": args.dataset_code,
            "conditions": ["full", "pmhc", "tcr", "five isolated chains"],
            "coordinate_policy": "All retained atom coordinates are unchanged.",
            "causal_boundary": (
                "These are fixed-bound-coordinate context responses, not matched apo-to-holo measurements."
            ),
            "excluded_exact_training_overlap": sorted(EXACT_TRAINING_OVERLAP),
        },
        "timing": {
            "condition_runs": int(len(timings)),
            "computed_runs": int((timings["status"] == "ok").sum()),
            "reused_runs": int((timings["status"] == "reused").sum()),
            "total_seconds": float(timings["end_to_end_seconds"].sum()),
            "median_seconds": float(timings["end_to_end_seconds"].median()),
        },
        "analysis": analysis,
        "artifacts": {
            "all_predictions": str((raw_dir / "all_condition_predictions.csv").resolve()),
            "matched_responses": str((raw_dir / "matched_residue_context_responses.csv").resolve()),
            "analysis_directory": str(analysis_dir.resolve()),
        },
    }
    json_dump(output_dir / "itrace_residue_biology_report.json", json_safe(report))
    print(json.dumps(json_safe(report["timing"]), indent=2), flush=True)


if __name__ == "__main__":
    main()
