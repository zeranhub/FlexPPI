"""Create a leakage-safe, difficulty-balanced five-fold FlexPPI split.

The v2 biological/homology groups are treated as indivisible units.  Fold
assignment uses only pre-training dataset metadata and ground-truth label
distribution summaries; no model prediction or benchmark score is consumed.

The script writes the canonical local split directory and a remote-ready ZIP
whose manifest paths are normalised for /ranze/repos/FlexPPI.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import random
import shutil
import statistics
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT
    / "data"
    / "splits"
    / "state_rmsf_hchi_v2_performance_first"
    / "all_ready.jsonl"
)
DEFAULT_OUTPUT = (
    ROOT / "data" / "splits" / "state_rmsf_hchi_v3_difficulty_balanced"
)
DEFAULT_BUNDLE = (
    ROOT
    / "artifacts"
    / "remote_upload"
    / "FlexPPI_v3_difficulty_balanced_5fold.zip"
)
VERSION = "state_rmsf_hchi_v3_difficulty_balanced"
FOLDS = 5
QUANTILE_BINS = 5


def json_dump(path: Path, value: object) -> None:
    write_text_lf(
        path, json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    )


def write_text_lf(path: Path, text: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def jsonl_write(path: Path, records: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_digest(*values: object) -> str:
    text = "\x1f".join(str(value) for value in values)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def remote_path(value: object) -> object:
    """Normalise known local paths to the remote FlexPPI layout."""

    if not isinstance(value, str):
        return value
    normal = value.replace("\\", "/")
    lower = normal.lower()
    replacements = [
        ("e:/project/flexppi", "/ranze/repos/FlexPPI"),
        ("e:/data/itrace", "/ranze/data/ITrace"),
        (
            "e:/data/protein_md_datasets/atlas_extracted",
            "/ranze/data/ATLAS_extracted",
        ),
        (
            "d:/data/dynarepo_non_nucleic_full",
            "/ranze/data/dynarepo/dynarepo_non_nucleic_full",
        ),
    ]
    for source, target in replacements:
        if lower == source or lower.startswith(source + "/"):
            return target + normal[len(source) :]
    return value if not lower.startswith(("e:/", "d:/")) else normal


def local_label_path(record: Mapping[str, object]) -> Path:
    raw = str(record["label_path"])
    direct = Path(raw)
    if direct.is_file():
        return direct
    normal = raw.replace("\\", "/")
    prefixes = [
        ("/ranze/repos/FlexPPI/", ROOT),
        ("E:/Project/FlexPPI/", ROOT),
    ]
    for prefix, target in prefixes:
        if normal.lower().startswith(prefix.lower()):
            candidate = target / normal[len(prefix) :]
            if candidate.is_file():
                return candidate
    dataset = str(record["dataset"]).lower()
    folders = {
        "atlas": "atlas_labels",
        "dynarepo": "dynarepo_labels",
        "itrace": "itrace_performance_first_labels",
    }
    candidate = ROOT / "artifacts" / "data" / folders[dataset] / Path(normal).name
    if candidate.is_file():
        return candidate
    raise FileNotFoundError("Cannot resolve label for {}: {}".format(record["sample_id"], raw))


def metric_values(label_path: Path, metric: str) -> np.ndarray:
    with np.load(label_path, allow_pickle=False) as label:
        names = [str(value) for value in label["metric_names"].tolist()]
        if metric not in names:
            return np.empty(0, dtype=float)
        index = names.index(metric)
        values = np.asarray(label["values"][:, index], dtype=float)
        mask = np.asarray(label["mask"][:, index], dtype=bool) & np.isfinite(values)
        return values[mask]


def describe_values(values: np.ndarray) -> Dict[str, float]:
    if len(values) == 0:
        return {
            "valid": 0,
            "mean": float("nan"),
            "std": float("nan"),
            "iqr": float("nan"),
            "p95": float("nan"),
        }
    q25, q75, p95 = np.quantile(values, [0.25, 0.75, 0.95])
    return {
        "valid": int(len(values)),
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "iqr": float(q75 - q25),
        "p95": float(p95),
    }


def percentile_ranks(values: Sequence[float]) -> List[float]:
    """Return deterministic average ranks in [0, 1], preserving ties."""

    array = np.asarray(values, dtype=float)
    if len(array) <= 1:
        return [0.5] * len(array)
    order = np.argsort(array, kind="mergesort")
    ranks = np.empty(len(array), dtype=float)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and array[order[end]] == array[order[start]]:
            end += 1
        average = (start + end - 1) / 2.0
        ranks[order[start:end]] = average / (len(array) - 1)
        start = end
    return ranks.tolist()


def bin_index(rank: float) -> int:
    return min(QUANTILE_BINS - 1, max(0, int(math.floor(rank * QUANTILE_BINS))))


def load_legacy_cohorts() -> Dict[str, str]:
    path = ROOT / "artifacts" / "evaluation" / "additional_baselines" / "manifests" / "static.tsv"
    if not path.is_file():
        return {}
    result: Dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[str(row["sample_id"]).lower()] = str(row["cohort"])
    return result


def prepare_samples(records: List[MutableMapping[str, object]]) -> List[Dict[str, object]]:
    legacy = load_legacy_cohorts()
    samples: List[Dict[str, object]] = []
    for record in records:
        metadata = dict(record.get("metadata", {}))
        label_path = local_label_path(record)
        rmsf = describe_values(metric_values(label_path, "rmsf"))
        hchi = describe_values(metric_values(label_path, "chi_entropy"))
        residue_count = int(metadata.get("residue_count") or rmsf["valid"])
        chain_count = max(1, len(record.get("target_chains", [])))
        sample_id = str(record["sample_id"])
        samples.append(
            {
                "sample_id": sample_id,
                "sample_key": sample_id.lower(),
                "dataset": str(record["dataset"]),
                "source_group": str(record["source_group"]),
                "group": str(metadata["cluster_id"]),
                "old_fold": int(metadata["fold"]),
                "residue_count": residue_count,
                "chain_count": chain_count,
                "rmsf_valid": int(rmsf["valid"]),
                "rmsf_mean": float(rmsf["mean"]),
                "rmsf_std": float(rmsf["std"]),
                "rmsf_iqr": float(rmsf["iqr"]),
                "rmsf_p95": float(rmsf["p95"]),
                "hchi_valid": int(hchi["valid"]),
                "hchi_mean": float(hchi["mean"]),
                "hchi_std": float(hchi["std"]),
                "hchi_coverage": float(hchi["valid"] / max(rmsf["valid"], 1)),
                "legacy_cohort": legacy.get(sample_id.lower(), "none"),
                "record": record,
                "label_local_path": str(label_path),
            }
        )

    rank_fields = [
        "residue_count",
        "rmsf_mean",
        "rmsf_std",
        "rmsf_p95",
        "hchi_coverage",
    ]
    strata: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    for index, sample in enumerate(samples):
        strata[(str(sample["dataset"]), str(sample["source_group"]))].append(index)
    for indices in strata.values():
        for field in rank_fields:
            values = [float(samples[index][field]) for index in indices]
            finite = [value for value in values if np.isfinite(value)]
            replacement = statistics.median(finite) if finite else 0.0
            values = [value if np.isfinite(value) else replacement for value in values]
            for index, rank in zip(indices, percentile_ranks(values)):
                samples[index][field + "_rank"] = rank
                samples[index][field + "_bin"] = bin_index(rank)
        for index in indices:
            sample = samples[index]
            difficulty = float(
                np.mean(
                    [
                        sample["residue_count_rank"],
                        sample["rmsf_mean_rank"],
                        sample["rmsf_std_rank"],
                        sample["rmsf_p95_rank"],
                    ]
                )
            )
            sample["difficulty_score"] = difficulty
        difficulty_values = [float(samples[index]["difficulty_score"]) for index in indices]
        for index, rank in zip(indices, percentile_ranks(difficulty_values)):
            samples[index]["difficulty_bin"] = bin_index(rank)
    return samples


def feature_weight(name: str) -> float:
    if name == "total_systems":
        return 1.0
    if name.startswith("dataset::"):
        return 6.0
    if name.startswith("source::"):
        return 14.0
    if name.startswith("difficulty::"):
        return 8.0
    if name.startswith("rmsf_mean_bin::") or name.startswith("rmsf_std_bin::"):
        return 6.0
    if name.startswith("length_bin::") or name.startswith("rmsf_p95_bin::"):
        return 4.0
    if name.startswith("chain_count::"):
        return 3.0
    if name.startswith("hchi_coverage_bin::"):
        return 2.0
    if name.startswith("legacy::"):
        return 8.0
    if name.startswith("rmsf_valid_sum::") or name.startswith("hchi_valid_sum::"):
        return 3.0
    return 1.0


def sample_features(sample: Mapping[str, object]) -> Dict[str, float]:
    dataset = str(sample["dataset"])
    source = str(sample["source_group"])
    prefix = dataset + "::" + source
    chain_bin = "1" if int(sample["chain_count"]) == 1 else ("2" if int(sample["chain_count"]) == 2 else "3plus")
    features = {
        "total_systems": 1.0,
        "dataset::" + dataset: 1.0,
        "source::" + prefix: 1.0,
        "difficulty::{}::{}".format(prefix, sample["difficulty_bin"]): 1.0,
        "length_bin::{}::{}".format(prefix, sample["residue_count_bin"]): 1.0,
        "rmsf_mean_bin::{}::{}".format(prefix, sample["rmsf_mean_bin"]): 1.0,
        "rmsf_std_bin::{}::{}".format(prefix, sample["rmsf_std_bin"]): 1.0,
        "rmsf_p95_bin::{}::{}".format(prefix, sample["rmsf_p95_bin"]): 1.0,
        "hchi_coverage_bin::{}::{}".format(prefix, sample["hchi_coverage_bin"]): 1.0,
        "chain_count::{}::{}".format(dataset, chain_bin): 1.0,
        "rmsf_valid_sum::" + dataset: float(sample["rmsf_valid"]),
        "hchi_valid_sum::" + dataset: float(sample["hchi_valid"]),
    }
    legacy = str(sample["legacy_cohort"])
    if legacy != "none":
        features["legacy::" + legacy] = 1.0
        features["legacy_difficulty::{}::{}".format(legacy, sample["difficulty_bin"])] = 1.0
    return features


def add_features(target: MutableMapping[str, float], source: Mapping[str, float]) -> None:
    for name, value in source.items():
        target[name] = target.get(name, 0.0) + float(value)


def build_group_features(samples: Sequence[Mapping[str, object]]) -> Dict[str, Dict[str, float]]:
    groups: Dict[str, Dict[str, float]] = defaultdict(dict)
    for sample in samples:
        add_features(groups[str(sample["group"])], sample_features(sample))
    return dict(groups)


def objective(
    totals: Sequence[Mapping[str, float]],
    targets: Mapping[str, float],
    weights: Mapping[str, float],
) -> float:
    result = 0.0
    for fold in range(FOLDS):
        for name, target in targets.items():
            scale = max(abs(target), 1.0)
            result += weights[name] * ((totals[fold].get(name, 0.0) - target) / scale) ** 2
    return result


def move_delta(
    current: Mapping[str, float],
    features: Mapping[str, float],
    targets: Mapping[str, float],
    weights: Mapping[str, float],
    direction: float,
) -> float:
    result = 0.0
    for name, value in features.items():
        target = targets[name]
        scale = max(abs(target), 1.0)
        before = (current.get(name, 0.0) - target) / scale
        after = (current.get(name, 0.0) + direction * value - target) / scale
        result += weights[name] * (after * after - before * before)
    return result


def assign_groups(
    group_features: Mapping[str, Mapping[str, float]],
    trials: int,
    base_seed: int,
) -> Tuple[Dict[str, int], float, Dict[str, object]]:
    grand: Dict[str, float] = {}
    for features in group_features.values():
        add_features(grand, features)
    targets = {name: value / FOLDS for name, value in grand.items()}
    weights = {name: feature_weight(name) for name in grand}

    def priority(item: Tuple[str, Mapping[str, float]]) -> Tuple[float, int]:
        _, features = item
        hardest = max(
            math.sqrt(weights[name]) * value / max(targets[name], 1.0)
            for name, value in features.items()
        )
        return hardest, int(features.get("total_systems", 0))

    best_assignment: Dict[str, int] = {}
    best_totals: List[Dict[str, float]] = []
    best_score = float("inf")
    for trial in range(trials):
        seed = base_seed + trial * 104729
        ordered = sorted(
            group_features.items(),
            key=lambda item: (
                -priority(item)[0],
                -priority(item)[1],
                stable_digest(seed, item[0]),
            ),
        )
        totals: List[Dict[str, float]] = [dict() for _ in range(FOLDS)]
        assignment: Dict[str, int] = {}
        for group, features in ordered:
            candidates = []
            for fold in range(FOLDS):
                delta = move_delta(totals[fold], features, targets, weights, 1.0)
                candidates.append((delta, stable_digest(seed, group, fold), fold))
            chosen = min(candidates)[2]
            assignment[group] = chosen
            add_features(totals[chosen], features)

        random_order = list(group_features)
        rng = random.Random(seed)
        for _ in range(5):
            rng.shuffle(random_order)
            changed = False
            for group in random_order:
                source_fold = assignment[group]
                features = group_features[group]
                candidates = [(0.0, source_fold)]
                remove = move_delta(totals[source_fold], features, targets, weights, -1.0)
                for target_fold in range(FOLDS):
                    if target_fold == source_fold:
                        continue
                    add = move_delta(totals[target_fold], features, targets, weights, 1.0)
                    candidates.append((remove + add, target_fold))
                delta, chosen = min(candidates)
                if chosen != source_fold and delta < -1e-12:
                    for name, value in features.items():
                        totals[source_fold][name] -= value
                        totals[chosen][name] = totals[chosen].get(name, 0.0) + value
                    assignment[group] = chosen
                    changed = True
            if not changed:
                break

        score = objective(totals, targets, weights)
        if score < best_score:
            best_score = score
            best_assignment = dict(assignment)
            best_totals = [dict(values) for values in totals]

    details = {
        "trials": trials,
        "base_seed": base_seed,
        "feature_count": len(grand),
        "group_count": len(group_features),
    }
    return best_assignment, best_score, details


def assignment_score(
    samples: Sequence[Mapping[str, object]],
    assignment: Mapping[str, int],
    group_features: Mapping[str, Mapping[str, float]],
) -> float:
    grand: Dict[str, float] = {}
    for features in group_features.values():
        add_features(grand, features)
    targets = {name: value / FOLDS for name, value in grand.items()}
    weights = {name: feature_weight(name) for name in grand}
    totals: List[Dict[str, float]] = [dict() for _ in range(FOLDS)]
    for group, features in group_features.items():
        add_features(totals[assignment[group]], features)
    return objective(totals, targets, weights)


def finite_mean(values: Sequence[float]) -> float:
    valid = [float(value) for value in values if np.isfinite(value)]
    return float(np.mean(valid)) if valid else float("nan")


def finite_std(values: Sequence[float]) -> float:
    valid = [float(value) for value in values if np.isfinite(value)]
    return float(np.std(valid)) if valid else float("nan")


def write_assignment_table(path: Path, samples: Sequence[Mapping[str, object]]) -> None:
    fields = [
        "sample_id",
        "dataset",
        "source_group",
        "group",
        "old_fold",
        "fold",
        "residue_count",
        "chain_count",
        "rmsf_valid",
        "rmsf_mean",
        "rmsf_std",
        "rmsf_iqr",
        "rmsf_p95",
        "hchi_valid",
        "hchi_mean",
        "hchi_std",
        "hchi_coverage",
        "difficulty_score",
        "difficulty_bin",
        "legacy_cohort",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for sample in sorted(samples, key=lambda value: str(value["sample_key"])):
            row = {field: sample[field] for field in fields}
            writer.writerow(row)


def write_balance_audit(path: Path, samples: Sequence[Mapping[str, object]]) -> None:
    fields = [
        "dataset",
        "source_group",
        "fold",
        "systems",
        "groups",
        "rmsf_valid_residues",
        "hchi_valid_residues",
        "length_mean",
        "length_std",
        "rmsf_mean_of_system_means",
        "rmsf_mean_of_system_stds",
        "rmsf_mean_of_system_p95",
        "hchi_mean_coverage",
        "difficulty_bin_0",
        "difficulty_bin_1",
        "difficulty_bin_2",
        "difficulty_bin_3",
        "difficulty_bin_4",
    ]
    rows = []
    strata = sorted({(str(s["dataset"]), str(s["source_group"])) for s in samples})
    for dataset, source in strata:
        for fold in range(FOLDS):
            selected = [
                sample
                for sample in samples
                if sample["dataset"] == dataset
                and sample["source_group"] == source
                and sample["fold"] == fold
            ]
            bins = Counter(int(sample["difficulty_bin"]) for sample in selected)
            rows.append(
                {
                    "dataset": dataset,
                    "source_group": source,
                    "fold": fold,
                    "systems": len(selected),
                    "groups": len({str(sample["group"]) for sample in selected}),
                    "rmsf_valid_residues": sum(int(sample["rmsf_valid"]) for sample in selected),
                    "hchi_valid_residues": sum(int(sample["hchi_valid"]) for sample in selected),
                    "length_mean": finite_mean([float(sample["residue_count"]) for sample in selected]),
                    "length_std": finite_std([float(sample["residue_count"]) for sample in selected]),
                    "rmsf_mean_of_system_means": finite_mean([float(sample["rmsf_mean"]) for sample in selected]),
                    "rmsf_mean_of_system_stds": finite_mean([float(sample["rmsf_std"]) for sample in selected]),
                    "rmsf_mean_of_system_p95": finite_mean([float(sample["rmsf_p95"]) for sample in selected]),
                    "hchi_mean_coverage": finite_mean([float(sample["hchi_coverage"]) for sample in selected]),
                    **{"difficulty_bin_{}".format(index): bins[index] for index in range(QUANTILE_BINS)},
                }
            )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def leakage_audit(samples: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    checks: Dict[str, Dict[str, set]] = {
        "biological_group": defaultdict(set),
        "biological_system": defaultdict(set),
        "pdb_id": defaultdict(set),
        "complex_signature": defaultdict(set),
        "label_sha256": defaultdict(set),
        "single_chain_sequence_cluster": defaultdict(set),
    }
    digest_cache: Dict[str, str] = {}
    for sample in samples:
        fold = int(sample["fold"])
        record = sample["record"]
        metadata = dict(record.get("metadata", {}))
        checks["biological_group"][str(sample["group"])].add(fold)
        checks["biological_system"][str(metadata.get("biological_system"))].add(fold)
        for pdb_id in metadata.get("pdb_ids", []):
            checks["pdb_id"][str(pdb_id).upper()].add(fold)
        signature = str(metadata.get("complex_signature") or "")
        if signature:
            checks["complex_signature"][signature].add(fold)
        label_path = str(sample["label_local_path"])
        if label_path not in digest_cache:
            digest_cache[label_path] = sha256(Path(label_path))
        checks["label_sha256"][digest_cache[label_path]].add(fold)
        if int(sample["chain_count"]) == 1:
            for cluster in metadata.get("sequence_clusters", []):
                checks["single_chain_sequence_cluster"][str(cluster)].add(fold)
    violations = {
        name: sorted(key for key, folds in values.items() if len(folds) > 1)
        for name, values in checks.items()
    }
    return {
        "passed": not any(violations.values()),
        "violations": violations,
        "checked_systems": len(samples),
        "note": "Complexes may share individual partner-chain families across folds; complete biological groups and complex signatures remain indivisible.",
    }


def create_configs(bundle_root: Path) -> None:
    config_dir = bundle_root / "configs"
    script_dir = bundle_root / "scripts"
    config_dir.mkdir(parents=True, exist_ok=True)
    script_dir.mkdir(parents=True, exist_ok=True)
    for fold in range(FOLDS):
        text = f"""stage: pretrain
seed: 7
device: auto
output_dir: artifacts/performance_first_v3_difficulty_balanced/fold_{fold}

data:
  state_manifest: data/splits/{VERSION}/fold_{fold}_train.jsonl
  validation_state_manifest: null
  test_state_manifest: data/splits/{VERSION}/fold_{fold}_test.jsonl
  pair_manifest: null
  max_residues: 1024
  batch_size: 1
  num_workers: 0
  path_mappings: {{}}

model:
  sequence_backend: learned
  esm_model_name: facebook/esm2_t6_8M_UR50D
  freeze_esm: true
  single_dim: 128
  pair_dim: 48
  sequence_dim: 128
  ipa_blocks: 4
  ipa_heads: 8
  ipa_scalar_dim: 16
  ipa_qk_points: 4
  ipa_v_points: 8
  rbf_bins: 24
  dropout: 0.1
  use_condition_features: false
  enable_response_head: false
  metrics: [rmsf, chi_entropy]

train:
  epochs: 40
  val_fraction: 0.0
  learning_rate: 0.0003
  weight_decay: 0.0001
  grad_clip: 1.0
  log_every: 25
  checkpoint_every: 1
  sampling_strategy: source_system_balanced
  samples_per_epoch: null
  metric_weights:
    rmsf: 1.0
    chi_entropy: 1.0
"""
        write_text_lf(config_dir / f"performance_first_v3_balanced_fold{fold}.yaml", text)

    launch = f"""#!/usr/bin/env bash
set -euo pipefail
REPO=/ranze/repos/FlexPPI
VENV=\"$REPO/.venv-linux/bin\"
SESSION=flexppi_v3_balanced
LOGDIR=\"$REPO/artifacts/performance_first_v3_difficulty_balanced/logs\"
mkdir -p \"$LOGDIR\"
tmux has-session -t \"$SESSION\" 2>/dev/null || tmux new-session -d -s \"$SESSION\" -n status
for fold in 0 1 2 3 4; do
  LOG=\"$LOGDIR/fold_${{fold}}.log\"
  if [ -f \"$LOG\" ]; then echo \"fold $fold already started; skipping\"; continue; fi
  tmux new-window -t \"$SESSION\" -n \"fold$fold\" \\
    \"cd $REPO && PYTHONUNBUFFERED=1 $VENV/flexppi train --config configs/performance_first_v3_balanced_fold${{fold}}.yaml 2>&1 | tee $LOG\"
done
echo \"Attach with: tmux attach -t $SESSION\"
"""
    write_text_lf(script_dir / "launch_performance_first_v3_balanced.sh", launch)
    evaluate = f"""#!/usr/bin/env bash
set -euo pipefail
REPO=/ranze/repos/FlexPPI
VENV=\"$REPO/.venv-linux/bin\"
OUT=\"$REPO/artifacts/performance_first_v3_difficulty_balanced/evaluation\"
mkdir -p \"$OUT\"
for fold in 0 1 2 3 4; do
  CKPT=\"$REPO/artifacts/performance_first_v3_difficulty_balanced/fold_${{fold}}/epoch_040.pt\"
  TARGET=\"$OUT/fold_${{fold}}_test.json\"
  [ -f \"$CKPT\" ] || {{ echo \"missing $CKPT\"; continue; }}
  [ ! -f \"$TARGET\" ] || {{ echo \"$TARGET exists; skipping\"; continue; }}
  cd \"$REPO\"
  \"$VENV/flexppi\" evaluate \\
    --config \"configs/performance_first_v3_balanced_fold${{fold}}.yaml\" \\
    --checkpoint \"$CKPT\" --split test --bootstrap 200 --output \"$TARGET\" \\
    2>&1 | tee \"$OUT/fold_${{fold}}_eval.log\"
done
"""
    write_text_lf(
        script_dir / "evaluate_performance_first_v3_balanced_epoch40.sh", evaluate
    )


def create_bundle(output_dir: Path, bundle_path: Path) -> None:
    staging = bundle_path.parent / (bundle_path.stem + "_contents")
    if staging.exists():
        shutil.rmtree(staging)
    split_target = staging / "data" / "splits" / VERSION
    split_target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(output_dir, split_target)
    create_configs(staging)
    readme = f"""# FlexPPI v3 difficulty-balanced remote training bundle

Unzip this archive into `/ranze/repos/FlexPPI` so that `data/`, `configs/`,
and `scripts/` merge with the existing repository.  The JSONL manifests use
remote `/ranze/...` paths and the supplied YAML files therefore use an empty
`path_mappings` dictionary.

Verify the split files:

```bash
cd /ranze/repos/FlexPPI/data/splits/{VERSION}
sha256sum -c SHA256SUMS
```

Start training:

```bash
cd /ranze/repos/FlexPPI
bash scripts/launch_performance_first_v3_balanced.sh
```

Evaluate the fixed epoch-40 checkpoints after all folds finish:

```bash
bash scripts/evaluate_performance_first_v3_balanced_epoch40.sh
```

This split was balanced with intrinsic metadata and label-distribution
summaries only.  No FlexPPI or baseline prediction score was used.
"""
    write_text_lf(staging / "REMOTE_README.md", readme)
    bundle_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(bundle_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(staging.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(staging).as_posix())
    shutil.rmtree(staging)
    (bundle_path.with_suffix(bundle_path.suffix + ".sha256")).write_text(
        "{}  {}\n".format(sha256(bundle_path), bundle_path.name), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument("--seed", type=int, default=20260906)
    parser.add_argument("--trials", type=int, default=96)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError("Refusing to overwrite existing split directory: {}".format(args.output_dir))
    if args.bundle.exists():
        raise FileExistsError("Refusing to overwrite existing bundle: {}".format(args.bundle))

    records = [json.loads(line) for line in args.input.read_text(encoding="utf-8").splitlines() if line.strip()]
    sample_keys = [str(record["sample_id"]).lower() for record in records]
    if len(sample_keys) != len(set(sample_keys)):
        raise ValueError("Case-normalised sample IDs are not unique")
    samples = prepare_samples(records)
    group_features = build_group_features(samples)
    assignment, new_score, optimiser = assign_groups(group_features, args.trials, args.seed)
    old_assignment: Dict[str, int] = {}
    for sample in samples:
        group = str(sample["group"])
        old = int(sample["old_fold"])
        if group in old_assignment and old_assignment[group] != old:
            raise ValueError("Existing biological group spans v2 folds: {}".format(group))
        old_assignment[group] = old
        sample["fold"] = int(assignment[group])
    old_score = assignment_score(samples, old_assignment, group_features)

    audit = leakage_audit(samples)
    if not audit["passed"]:
        raise RuntimeError("Leakage audit failed: {}".format(audit["violations"]))

    args.output_dir.mkdir(parents=True)
    output_records = []
    for sample in samples:
        record = copy.deepcopy(sample["record"])
        metadata = dict(record.get("metadata", {}))
        metadata["previous_fold_v2"] = int(sample["old_fold"])
        metadata["fold"] = int(sample["fold"])
        metadata["split_version"] = VERSION
        metadata["difficulty_balance"] = {
            "difficulty_bin": int(sample["difficulty_bin"]),
            "difficulty_score": round(float(sample["difficulty_score"]), 8),
            "length_bin": int(sample["residue_count_bin"]),
            "rmsf_mean_bin": int(sample["rmsf_mean_bin"]),
            "rmsf_std_bin": int(sample["rmsf_std_bin"]),
            "rmsf_p95_bin": int(sample["rmsf_p95_bin"]),
            "hchi_coverage_bin": int(sample["hchi_coverage_bin"]),
        }
        record["metadata"] = metadata
        for field in ["label_path", "structure_path"]:
            record[field] = remote_path(record.get(field))
        record["trajectory_paths"] = [remote_path(path) for path in record.get("trajectory_paths", [])]
        output_records.append(record)

    jsonl_write(args.output_dir / "all_ready.jsonl", output_records)
    for fold in range(FOLDS):
        train = [record for record in output_records if int(record["metadata"]["fold"]) != fold]
        test = [record for record in output_records if int(record["metadata"]["fold"]) == fold]
        jsonl_write(args.output_dir / "fold_{}_train.jsonl".format(fold), train)
        jsonl_write(args.output_dir / "fold_{}_test.jsonl".format(fold), test)
    write_assignment_table(args.output_dir / "fold_assignments.tsv", samples)
    write_balance_audit(args.output_dir / "balance_audit.tsv", samples)
    json_dump(args.output_dir / "leakage_audit.json", audit)

    fold_summary = {}
    for fold in range(FOLDS):
        selected = [sample for sample in samples if int(sample["fold"]) == fold]
        fold_summary[str(fold)] = {
            "systems": len(selected),
            "biological_groups": len({str(sample["group"]) for sample in selected}),
            "datasets": dict(sorted(Counter(str(sample["dataset"]) for sample in selected).items())),
            "source_groups": dict(sorted(Counter(str(sample["source_group"]) for sample in selected).items())),
            "rmsf_valid_residues": sum(int(sample["rmsf_valid"]) for sample in selected),
            "hchi_valid_residues": sum(int(sample["hchi_valid"]) for sample in selected),
            "mean_residue_count": finite_mean([float(sample["residue_count"]) for sample in selected]),
            "mean_system_rmsf": finite_mean([float(sample["rmsf_mean"]) for sample in selected]),
        }
    summary = {
        "version": VERSION,
        "created_from": str(args.input.relative_to(ROOT)).replace("\\", "/"),
        "systems": len(samples),
        "biological_groups": len(group_features),
        "folds": FOLDS,
        "seed": args.seed,
        "optimisation_trials": args.trials,
        "uses_model_predictions": False,
        "difficulty_definition": [
            "within-dataset/source-group residue-count quintile",
            "within-dataset/source-group RMSF mean quintile",
            "within-dataset/source-group RMSF standard-deviation quintile",
            "within-dataset/source-group RMSF 95th-percentile quintile",
            "within-dataset/source-group H-chi label-coverage quintile",
            "target-chain-count category",
        ],
        "indivisible_unit": "v2 biological/homology cluster_id",
        "objective": {
            "v2_assignment_score": old_score,
            "v3_assignment_score": new_score,
            "relative_reduction": 1.0 - new_score / old_score,
            **optimiser,
        },
        "fold_summary": fold_summary,
        "leakage_audit_passed": audit["passed"],
        "remote_path_policy": {
            "repository": "/ranze/repos/FlexPPI",
            "itrace": "/ranze/data/ITrace",
            "atlas": "/ranze/data/ATLAS_extracted",
            "dynarepo": "/ranze/data/dynarepo/dynarepo_non_nucleic_full",
        },
    }
    json_dump(args.output_dir / "summary.json", summary)
    readme = f"""# {VERSION}

This directory contains 3,068 ready ATLAS, DynaRepo, and ITrace systems in a
five-fold, difficulty-balanced, group-safe split.  Existing v2 biological and
homology `cluster_id` values are indivisible.  Fold optimisation uses dataset,
source group, target-chain count, label availability, protein length, and
ground-truth RMSF distribution summaries.  It never uses FlexPPI or baseline
prediction scores.

Training consumes `fold_N_train.jsonl`; evaluation consumes
`fold_N_test.jsonl`.  `fold_assignments.tsv`, `balance_audit.tsv`,
`summary.json`, and `leakage_audit.json` are audit artifacts and should be kept
with the split.  Paths inside all JSONL files are normalised to the remote
`/ranze/...` layout.

Because this split was designed after the v2 result was inspected, report it as
a new, frozen v3 protocol rather than silently replacing the original v2
evaluation.
"""
    write_text_lf(args.output_dir / "README.md", readme)

    checksum_files = sorted(path for path in args.output_dir.iterdir() if path.is_file() and path.name != "SHA256SUMS")
    (args.output_dir / "SHA256SUMS").write_text(
        "".join("{}  {}\n".format(sha256(path), path.name) for path in checksum_files),
        encoding="utf-8",
    )
    create_bundle(args.output_dir, args.bundle)
    print(json.dumps({"output_dir": str(args.output_dir), "bundle": str(args.bundle), "summary": summary}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
