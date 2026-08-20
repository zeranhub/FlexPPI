"""Evaluate official AlphaFlow ensembles on their frozen ATLAS overlap.

AlphaFlow is a monomer ensemble generator rather than a direct RMSF regressor.
This script converts each official multi-model PDB into a C-alpha RMSF curve,
using a Kabsch fit of every generated conformation to the corresponding frozen
input structure.  The curve is evaluated with the same per-protein metrics as
the direct FlexPPI baselines.  Results are deliberately written to a separate
overlap table because the released AlphaFlow ensembles cover the authors' 82
protein test split, only six of which overlap this project's ATLAS-191 cohort.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, NamedTuple, Sequence, Tuple

import numpy as np
import pandas as pd

from flexppi.comparison import _method_summary, _system_metrics


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "artifacts" / "evaluation" / "additional_baselines"
STATIC_COMMON = BASE / "results" / "static_common_residue_predictions.csv"
STATIC_MANIFEST = BASE / "manifests" / "static.tsv"
RAW_OUTPUT = BASE / "raw" / "alphaflow_atlas_official_overlap.csv"
PER_SYSTEM_OUTPUT = BASE / "results" / "alphaflow_atlas_overlap_per_system_metrics.csv"
SUMMARY_OUTPUT = BASE / "results" / "alphaflow_atlas_overlap_summary_metrics.csv"
CONVERGENCE_OUTPUT = BASE / "results" / "alphaflow_atlas_overlap_convergence.csv"

MODEL_NAME = "alphaflow_12l_md_templates_distilled_202406"
SOURCE_COMMIT = "0408d7c89dac444a43a9089d7427ce470b0a5e67"

AMINO_ACID_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
}

METHODS = {
    "fastprotflex": ("fastprotflex_score", False),
    "prody_anm": ("anm_sq_fluct", False),
    "prody_gnm": ("gnm_sq_fluct", False),
    "backflip_structure": ("backflip_structure_rmsf_angstrom", True),
    "backflip_seq": ("backflip_seq_rmsf_angstrom", True),
    "rmsf_net_pdb": ("rmsf_net_pdb_score", False),
    "flexppi": ("flexppi_pretrain_rmsf_angstrom", True),
    "alphaflow": ("alphaflow_rmsf_angstrom", True),
}


class CAResidue(NamedTuple):
    chain: str
    number: str
    insertion: str
    name: str
    coord: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ensemble-dir", required=True)
    parser.add_argument("--bootstrap-repetitions", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--sample-counts", default="10,25,50,100,250")
    return parser.parse_args()


def ca_from_line(line: str) -> CAResidue:
    return CAResidue(
        chain=line[21:22].strip() or "_",
        number=line[22:26].strip(),
        insertion=line[26:27].strip(),
        name=line[17:20].strip().upper(),
        coord=np.asarray(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            dtype=np.float64,
        ),
    )


def parse_reference(path: str) -> List[CAResidue]:
    residues: List[CAResidue] = []
    seen: set[Tuple[str, str, str]] = set()
    seen_atoms = False
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ENDMDL") and seen_atoms:
                break
            if not line.startswith("ATOM  ") or len(line) < 54:
                continue
            seen_atoms = True
            if line[12:16].strip() != "CA" or line[16:17] not in (" ", "A"):
                continue
            residue = ca_from_line(line)
            key = residue[:3]
            if key not in seen:
                seen.add(key)
                residues.append(residue)
    if len(residues) < 3:
        raise ValueError(f"{path}: fewer than three reference CA residues")
    return residues


def parse_multimodel_pdb(path: Path) -> Tuple[List[CAResidue], np.ndarray]:
    models: List[List[CAResidue]] = []
    current: List[CAResidue] = []
    seen: set[Tuple[str, str, str]] = set()

    def finish_model() -> None:
        nonlocal current, seen
        if current:
            models.append(current)
        current = []
        seen = set()

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                finish_model()
                continue
            if line.startswith("ENDMDL"):
                finish_model()
                continue
            if not line.startswith("ATOM  ") or len(line) < 54:
                continue
            if line[12:16].strip() != "CA" or line[16:17] not in (" ", "A"):
                continue
            residue = ca_from_line(line)
            key = residue[:3]
            if key not in seen:
                seen.add(key)
                current.append(residue)
    finish_model()
    if len(models) < 2:
        raise ValueError(f"{path}: AlphaFlow ensemble contains fewer than two models")
    identities = [(row.chain, row.number, row.insertion, row.name) for row in models[0]]
    for model_index, model in enumerate(models[1:], start=2):
        observed = [(row.chain, row.number, row.insertion, row.name) for row in model]
        if observed != identities:
            raise ValueError(f"{path}: residue identities differ in MODEL {model_index}")
    coordinates = np.stack([[row.coord for row in model] for model in models])
    return models[0], coordinates


def sequence(residues: Sequence[CAResidue]) -> str:
    return "".join(AMINO_ACID_1.get(residue.name, "X") for residue in residues)


def alignment_pairs(left: str, right: str) -> List[Tuple[int, int]]:
    """Return exact-residue pairs from a deterministic global alignment."""

    rows, columns = len(left) + 1, len(right) + 1
    score = np.empty((rows, columns), dtype=np.int32)
    trace = np.empty((rows, columns), dtype=np.uint8)
    score[:, 0] = -2 * np.arange(rows)
    score[0, :] = -2 * np.arange(columns)
    trace[1:, 0] = 1
    trace[0, 1:] = 2
    trace[0, 0] = 0
    for i in range(1, rows):
        for j in range(1, columns):
            candidates = (
                score[i - 1, j - 1] + (2 if left[i - 1] == right[j - 1] else -1),
                score[i - 1, j] - 2,
                score[i, j - 1] - 2,
            )
            move = int(np.argmax(candidates))
            score[i, j] = candidates[move]
            trace[i, j] = move
    pairs: List[Tuple[int, int]] = []
    i, j = len(left), len(right)
    while i or j:
        move = int(trace[i, j])
        if move == 0:
            i -= 1
            j -= 1
            if left[i] == right[j]:
                pairs.append((i, j))
        elif move == 1:
            i -= 1
        else:
            j -= 1
    pairs.reverse()
    return pairs


def kabsch_align(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    centered_mobile = mobile - mobile_center
    centered_target = target - target_center
    left, _, right_t = np.linalg.svd(centered_mobile.T @ centered_target)
    if np.linalg.det(left @ right_t) < 0:
        left[:, -1] *= -1
    rotation = left @ right_t
    return centered_mobile @ rotation + target_center


def generated_rmsf(
    coordinates: np.ndarray,
    generated_indices: np.ndarray,
    reference_coordinates: np.ndarray,
    sample_count: int,
) -> Tuple[np.ndarray, np.ndarray]:
    frames = coordinates[:sample_count]
    aligned = np.stack(
        [
            kabsch_align(frame[generated_indices], reference_coordinates)
            for frame in frames
        ]
    )
    mean_structure = aligned.mean(axis=0)
    around_mean = np.sqrt(np.mean(np.sum(np.square(aligned - mean_structure), axis=-1), axis=0))
    around_reference = np.sqrt(
        np.mean(np.sum(np.square(aligned - reference_coordinates), axis=-1), axis=0)
    )
    return around_mean, around_reference


def discover_ensembles(directory: Path) -> Dict[str, Path]:
    files: Dict[str, Path] = {}
    for path in directory.rglob("*.pdb"):
        key = path.stem.lower()
        if key in files:
            raise ValueError(f"Duplicate AlphaFlow ensemble stem {key}")
        files[key] = path
    if not files:
        raise FileNotFoundError(f"No PDB ensembles found under {directory}")
    return files


def flatten_metrics(rows: Sequence[Mapping[str, object]]) -> pd.DataFrame:
    flattened: List[Dict[str, object]] = []
    for row in rows:
        item: Dict[str, object] = {
            key: value for key, value in row.items() if not isinstance(value, dict)
        }
        for method in METHODS:
            for metric, value in row.get(method, {}).items():
                item[f"{method}__{metric}"] = value
        flattened.append(item)
    return pd.DataFrame(flattened)


def evaluate_groups(
    common: pd.DataFrame,
    repetitions: int,
    seed: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows: List[Dict[str, object]] = []
    for sample_id, group in common.groupby("sample_id", sort=True):
        target = group["target_rmsf_angstrom"].to_numpy(float)
        metrics: Dict[str, object] = {
            "sample_id": sample_id,
            "residue_count": len(group),
            "frame_count": int(group["alphaflow_frame_count"].iloc[0]),
        }
        for method, (column, absolute) in METHODS.items():
            result = _system_metrics(
                group[column].to_numpy(float), target, absolute_scale=absolute
            )
            if not result:
                raise ValueError(f"{sample_id}: metric calculation failed for {method}")
            metrics[method] = result
        rows.append(metrics)
    summary_rows: List[Dict[str, object]] = []
    for method in METHODS:
        summary = _method_summary(rows, method, repetitions, seed)
        for metric, values in summary.items():
            summary_rows.append(
                {
                    "cohort": "atlas_official_alphaflow_overlap",
                    "sample_count": len(rows),
                    "method": method,
                    "metric": metric,
                    **values,
                }
            )
    return flatten_metrics(rows), pd.DataFrame(summary_rows)


def main() -> None:
    args = parse_args()
    ensembles = discover_ensembles(Path(args.ensemble_dir))
    manifest = pd.read_csv(STATIC_MANIFEST, sep="\t", dtype=str).fillna("")
    manifest = manifest.loc[manifest["cohort"] == "atlas_single_chain"].copy()
    manifest["sample_key"] = manifest["sample_id"].str.lower()
    common = pd.read_csv(STATIC_COMMON, dtype={"sample_id": str, "author_residue_id": str})
    common = common.loc[common["cohort"] == "atlas_single_chain"].copy()
    common["sample_key"] = common["sample_id"].str.lower()
    overlap = sorted(set(manifest["sample_key"]) & set(ensembles))
    if not overlap:
        raise ValueError("The AlphaFlow ensemble directory has no frozen ATLAS overlap")

    requested_counts = sorted(
        {int(value) for value in args.sample_counts.split(",") if int(value) >= 2}
    )
    raw_rows: List[Dict[str, object]] = []
    convergence_rows: List[Dict[str, object]] = []
    for sample_key in overlap:
        manifest_row = manifest.loc[manifest["sample_key"] == sample_key].iloc[0]
        sample_id = str(manifest_row["sample_id"])
        reference = parse_reference(str(manifest_row["pdb_path"]))
        generated, coordinates = parse_multimodel_pdb(ensembles[sample_key])
        pairs = alignment_pairs(sequence(generated), sequence(reference))
        if len(pairs) < 3:
            raise ValueError(f"{sample_id}: fewer than three exact aligned CA residues")
        generated_indices = np.asarray([left for left, _ in pairs], dtype=int)
        reference_indices = np.asarray([right for _, right in pairs], dtype=int)
        reference_coordinates = np.stack(
            [reference[index].coord for index in reference_indices]
        )
        counts = sorted(set([count for count in requested_counts if count <= len(coordinates)] + [len(coordinates)]))
        full_mean: np.ndarray | None = None
        full_reference: np.ndarray | None = None
        for count in counts:
            rmsf_mean, rmsf_reference = generated_rmsf(
                coordinates, generated_indices, reference_coordinates, count
            )
            if count == len(coordinates):
                full_mean, full_reference = rmsf_mean, rmsf_reference
            target_rows = common.loc[common["sample_key"] == sample_key]
            lookup = {
                (
                    reference[reference_index].chain,
                    reference[reference_index].number + reference[reference_index].insertion,
                ): float(value)
                for reference_index, value in zip(reference_indices, rmsf_mean)
            }
            selected = target_rows.copy()
            selected["alphaflow_rmsf_angstrom"] = [
                lookup.get((str(row.chain), str(row.author_residue_id)), np.nan)
                for row in selected.itertuples()
            ]
            selected = selected.dropna(subset=["alphaflow_rmsf_angstrom"])
            result = _system_metrics(
                selected["alphaflow_rmsf_angstrom"].to_numpy(float),
                selected["target_rmsf_angstrom"].to_numpy(float),
                absolute_scale=True,
            )
            convergence_rows.append(
                {
                    "sample_id": sample_id,
                    "frame_count": count,
                    "mapped_residue_count": len(selected),
                    **result,
                }
            )
        assert full_mean is not None and full_reference is not None
        for output_index, (generated_index, reference_index) in enumerate(pairs):
            residue = reference[reference_index]
            raw_rows.append(
                {
                    "cohort": "atlas_official_alphaflow_overlap",
                    "sample_id": sample_id,
                    "residue_index": output_index,
                    "chain": residue.chain,
                    "author_residue_id": residue.number + residue.insertion,
                    "residue_number": residue.number,
                    "insertion_code": residue.insertion,
                    "amino_acid_3": residue.name,
                    "alphaflow_rmsf_angstrom": float(full_mean[output_index]),
                    "alphaflow_reference_displacement_angstrom": float(
                        full_reference[output_index]
                    ),
                    "alphaflow_frame_count": len(coordinates),
                    "model_name": MODEL_NAME,
                    "source_commit": SOURCE_COMMIT,
                }
            )
    raw = pd.DataFrame(raw_rows)
    key = ["sample_id", "chain", "author_residue_id"]
    selected_common = common.loc[common["sample_key"].isin(overlap)].copy()
    merged = selected_common.merge(
        raw[key + ["alphaflow_rmsf_angstrom", "alphaflow_frame_count"]],
        on=key,
        how="inner",
        validate="one_to_one",
    )
    if len(merged) != len(selected_common):
        raise ValueError(
            f"AlphaFlow residue mapping is incomplete: frozen={len(selected_common)}, "
            f"mapped={len(merged)}"
        )
    observed_samples = set(merged["sample_key"])
    if observed_samples != set(overlap):
        raise ValueError(
            f"AlphaFlow residue mapping lost samples: expected={overlap}, "
            f"observed={sorted(observed_samples)}"
        )
    per_system, summary = evaluate_groups(
        merged, args.bootstrap_repetitions, args.seed
    )
    RAW_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    PER_SYSTEM_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    raw.to_csv(RAW_OUTPUT, index=False)
    per_system.to_csv(PER_SYSTEM_OUTPUT, index=False)
    summary.to_csv(SUMMARY_OUTPUT, index=False)
    pd.DataFrame(convergence_rows).to_csv(CONVERGENCE_OUTPUT, index=False)
    print(
        f"AlphaFlow overlap evaluation complete: samples={len(overlap)}, "
        f"residues={len(merged)}, ensembles={args.ensemble_dir}"
    )


if __name__ == "__main__":
    main()
