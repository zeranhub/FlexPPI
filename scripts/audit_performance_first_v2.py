"""Independently validate v2 labels, folds, and frozen training configs."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from flexppi.config import load_config  # noqa: E402
from flexppi.constants import METRIC_NAMES  # noqa: E402
from flexppi.data.labels import load_label_file  # noqa: E402
from flexppi.data.records import load_state_records, remap_state_record  # noqa: E402
from flexppi.data.structures import parse_pdb  # noqa: E402


def _record_key(record) -> str:
    return f"{record.dataset}:{record.sample_id}"


def _metadata_values(records, name: str) -> set[str]:
    output = set()
    for record in records:
        value = record.metadata.get(name)
        if isinstance(value, list):
            output.update(str(item) for item in value if item)
        elif value:
            output.add(str(value))
    return output


def audit(split_dir: Path, config_paths: list[Path], output_path: Path) -> dict[str, object]:
    configs = [load_config(str(path)) for path in config_paths]
    mappings = configs[0].data.path_mappings
    all_records = [
        remap_state_record(record, mappings)
        for record in load_state_records(str(split_dir / "all_ready.jsonl"))
    ]
    all_keys = {_record_key(record) for record in all_records}
    violations: dict[str, list[str]] = {
        "config": [],
        "fold_membership": [],
        "fold_overlap": [],
        "missing_file": [],
        "label_schema": [],
    }

    train_appearances: Counter[str] = Counter()
    test_appearances: Counter[str] = Counter()
    fold_summary = {}
    for fold, config in enumerate(configs):
        expected_train = split_dir / f"fold_{fold}_train.jsonl"
        expected_test = split_dir / f"fold_{fold}_test.jsonl"
        config_checks = {
            "state_manifest": Path(str(config.data.state_manifest)).as_posix()
            == expected_train.relative_to(ROOT).as_posix(),
            "test_state_manifest": Path(str(config.data.test_state_manifest)).as_posix()
            == expected_test.relative_to(ROOT).as_posix(),
            "no_outer_validation": config.data.validation_state_manifest is None
            and config.train.val_fraction == 0.0,
            "max_residues": config.data.max_residues == 1024,
            "metrics": config.model.metrics == ["rmsf", "chi_entropy"],
            "structure_only": not config.model.use_condition_features,
            "no_response_head": not config.model.enable_response_head,
            "balanced_sampling": config.train.sampling_strategy == "source_system_balanced",
        }
        failures = [name for name, passed in config_checks.items() if not passed]
        violations["config"].extend(f"fold_{fold}:{name}" for name in failures)

        train_records = load_state_records(str(expected_train))
        test_records = load_state_records(str(expected_test))
        train_keys = {_record_key(record) for record in train_records}
        test_keys = {_record_key(record) for record in test_records}
        for key in train_keys:
            train_appearances[key] += 1
        for key in test_keys:
            test_appearances[key] += 1
        if train_keys & test_keys:
            violations["fold_overlap"].append(f"fold_{fold}:record")
        if train_keys | test_keys != all_keys:
            violations["fold_membership"].append(f"fold_{fold}:incomplete_union")
        for field_name in ("cluster_id", "pdb_ids", "complex_signature"):
            overlap = _metadata_values(train_records, field_name) & _metadata_values(
                test_records, field_name
            )
            if overlap:
                violations["fold_overlap"].append(f"fold_{fold}:{field_name}:{len(overlap)}")
        fold_summary[str(fold)] = {
            "train_records": len(train_records),
            "test_records": len(test_records),
            "config_checks": config_checks,
        }

    for key in sorted(all_keys):
        if train_appearances[key] != 4 or test_appearances[key] != 1:
            violations["fold_membership"].append(
                f"{key}:train={train_appearances[key]},test={test_appearances[key]}"
            )

    metric_counts = Counter()
    for record in all_records:
        structure_path = Path(record.structure_path)
        label_path = Path(record.label_path or "")
        if not structure_path.exists() or not label_path.exists():
            violations["missing_file"].append(_record_key(record))
            continue
        try:
            structure = parse_pdb(str(structure_path))
            values, variances, mask = load_label_file(str(label_path))
            if values.shape[0] != structure.num_residues or structure.num_residues > 1024:
                raise ValueError("structure_label_length_or_limit")
            if values.shape != variances.shape or values.shape != mask.shape:
                raise ValueError("tensor_shape_mismatch")
            for name in ("rmsf", "chi_entropy"):
                index = METRIC_NAMES.index(name)
                valid = mask[:, index] & np.isfinite(values[:, index])
                metric_counts[name] += int(valid.sum())
            if not np.any(mask[:, METRIC_NAMES.index("rmsf")]):
                raise ValueError("no_rmsf")
            if np.any(variances[mask] < 0) or not np.all(np.isfinite(variances[mask])):
                raise ValueError("invalid_variance")
        except (OSError, ValueError, KeyError) as exc:
            violations["label_schema"].append(f"{_record_key(record)}:{exc}")

    result: dict[str, object] = {
        "version": "state_rmsf_hchi_v2_performance_first",
        "passed": not any(violations.values()),
        "ready_records": len(all_records),
        "datasets": dict(sorted(Counter(record.dataset for record in all_records).items())),
        "metric_labels": dict(metric_counts),
        "folds": fold_summary,
        "appearance_rule": (
            "every ready record appears in four train manifests and one test manifest"
        ),
        "violations": violations,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    if not result["passed"]:
        raise RuntimeError(f"v2 readiness audit failed: {violations}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--split-dir",
        default=str(ROOT / "data/splits/state_rmsf_hchi_v2_performance_first"),
    )
    parser.add_argument(
        "--config",
        action="append",
        default=[str(ROOT / f"configs/performance_first_v2_fold{fold}.yaml") for fold in range(5)],
    )
    parser.add_argument(
        "--output",
        default=str(ROOT / "data/splits/state_rmsf_hchi_v2_performance_first/readiness_audit.json"),
    )
    args = parser.parse_args()
    result = audit(
        Path(args.split_dir),
        [Path(path) for path in args.config],
        Path(args.output),
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
