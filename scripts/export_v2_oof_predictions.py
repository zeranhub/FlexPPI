"""Export provenance-checked, final-epoch out-of-fold state predictions.

Every ready sample is inferred only with its held-out fold's epoch_040.pt.
No training, model download, checkpoint selection, or five-model averaging occurs.
The NPZ cache is the reusable source of truth for later analyses and figures.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from flexppi.config import ModelConfig  # noqa: E402
from flexppi.constants import METRIC_NAMES  # noqa: E402
from flexppi.data.dataset import collate_states, load_state_example  # noqa: E402
from flexppi.data.labels import load_label_file  # noqa: E402
from flexppi.data.records import StateRecord  # noqa: E402
from flexppi.data.structures import parse_pdb  # noqa: E402
from flexppi.evaluation import _correlation, _spearman  # noqa: E402
from flexppi.model.network import FlexPPI  # noqa: E402
from flexppi.trainer import move_to_device  # noqa: E402

CACHE_VERSION = 1
METRICS = ("rmsf", "chi_entropy")
SPLIT_DIRECTORY = "data/splits/state_rmsf_hchi_v2_performance_first"
PROTOCOL = {
    "name": "v2_final_epoch40_oof",
    "epoch": 40,
    "context": "full",
    "metrics": list(METRICS),
    "mask": "label_mask & target_mask & geometry_mask & output_mask",
    "scope": "combination-level CV; not zero-shot or unseen-family evaluation",
    "uncertainty": "model_variance and label_variance stored separately",
    "selection": "each sample uses only its held-out fold, no best.pt or ensemble",
}


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
                      allow_nan=False)


def digest_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_jsonl(path: Path) -> list:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def unique_records(records: list, label: str) -> dict:
    result = {}
    for record in records:
        sample = record["sample_id"]
        if sample in result:
            raise ValueError(f"Duplicate sample {sample!r} in {label}")
        result[sample] = record
    return result


def audit_inventory(split_dir: Path, expected_count: int = 3068) -> dict:
    """Audit all folds even if execution later requests only a smoke subset."""
    all_path = split_dir / "all_ready.jsonl"
    ready = unique_records(read_jsonl(all_path), "all_ready")
    if len(ready) != expected_count:
        raise ValueError(f"Expected {expected_count} ready samples, found {len(ready)}")
    assignments, folds, source_hashes = {}, {}, {"all_ready": sha256(all_path)}
    for fold in range(5):
        test_path = split_dir / f"fold_{fold}_test.jsonl"
        train_path = split_dir / f"fold_{fold}_train.jsonl"
        test = unique_records(read_jsonl(test_path), f"fold {fold} test")
        train = unique_records(read_jsonl(train_path), f"fold {fold} train")
        if set(train) & set(test):
            raise ValueError(f"Train/test sample overlap in fold {fold}")
        if set(train) | set(test) != set(ready):
            raise ValueError(f"Incomplete or unexpected fold {fold} membership")
        for key in ("biological_system", "cluster_id"):
            groups = lambda records: {
                r.get("metadata", {}).get(key) for r in records.values()
                if r.get("metadata", {}).get(key)
            }
            if groups(train) & groups(test):
                raise ValueError(f"Train/test {key} overlap in fold {fold}")
        for sample, record in {**train, **test}.items():
            if digest_json(record) != digest_json(ready[sample]):
                raise ValueError(f"Manifest record differs from all_ready: {sample}")
        for sample, record in test.items():
            if sample in assignments:
                raise ValueError(f"Sample is in multiple test folds: {sample}")
            if record.get("metadata", {}).get("fold") != fold:
                raise ValueError(f"Metadata fold disagrees with test membership: {sample}")
            assignments[sample] = fold
        folds[fold] = test
        source_hashes[str(fold)] = {"train": sha256(train_path), "test": sha256(test_path)}
    if set(assignments) != set(ready):
        raise ValueError("Each ready sample must occur in exactly one test fold")
    return {"ready": ready, "folds": folds, "assignments": assignments,
            "source_hashes": source_hashes}


def resolve_existing(path: str, mappings: Mapping[str, str]) -> Path:
    """Prefer existing literal paths, then bidirectional slash-normalised mappings."""
    direct = Path(path)
    if direct.is_file():
        return direct.resolve()
    normal = path.replace("\\", "/").rstrip("/")
    pairs = list(mappings.items()) + [(b, a) for a, b in mappings.items()]
    for source, target in sorted(pairs, key=lambda pair: len(pair[0]), reverse=True):
        prefix = str(source).replace("\\", "/").rstrip("/")
        if normal == prefix or normal.startswith(prefix + "/"):
            suffix = normal[len(prefix):].lstrip("/")
            candidate = Path(target).joinpath(*suffix.split("/"))
            if candidate.is_file():
                return candidate.resolve()
    raise FileNotFoundError(f"No existing input found for {path}")


def atomic_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def atomic_json(path: Path, value: Any) -> None:
    atomic_bytes(path, (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)
                        + "\n").encode("utf-8"))


def cache_name(sample: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", sample)[:100]
    return f"{safe}_{hashlib.sha256(sample.encode()).hexdigest()[:12]}.npz"


def validate_cache_arrays(arrays: dict, expected: dict) -> dict:
    if "metadata_json" not in arrays:
        raise ValueError("Cache has no provenance metadata")
    meta = json.loads(str(arrays["metadata_json"].item()))
    if meta.get("identity") != expected or meta.get("status") != "complete":
        raise ValueError("Cache identity/protocol differs; refusing to reuse or overwrite")
    n = meta["residue_count"]
    for key in ("chain", "author_residue_id", "residue_index", "aatype", "amino_acid",
                "target_mask", "geometry_mask", "output_mask"):
        if key not in arrays or arrays[key].shape != (n,):
            raise ValueError(f"Malformed cache field: {key}")
    for key in ("target", "label_mask", "label_variance", "mean", "model_variance", "valid"):
        if key not in arrays or arrays[key].shape != (n, len(METRICS)):
            raise ValueError(f"Malformed cache field: {key}")
    if arrays.get("metric_names", np.array([])).tolist() != list(METRICS):
        raise ValueError("Cache metric order differs")
    valid = arrays["label_mask"] & arrays["target_mask"][:, None]
    valid &= arrays["geometry_mask"][:, None] & arrays["output_mask"][:, None]
    if not np.array_equal(valid, arrays["valid"]):
        raise ValueError("Cache validity masks disagree")
    for key in ("mean", "model_variance", "target", "label_variance"):
        if not np.isfinite(arrays[key][valid]).all():
            raise ValueError(f"Non-finite valid values in cache: {key}")
    if np.any(arrays["model_variance"][valid] <= 0):
        raise ValueError("Non-positive model variance")
    if np.any(arrays["label_variance"][valid] < 0):
        raise ValueError("Negative label variance")
    return meta


def load_cache(path: Path, expected: dict) -> dict:
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    validate_cache_arrays(arrays, expected)
    return arrays


def save_cache(path: Path, arrays: dict, identity: dict) -> None:
    validate_cache_arrays(arrays, identity)
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    atomic_bytes(path, buffer.getvalue())


def metric_rows(arrays: dict, sample: str, dataset: str, fold: int) -> list:
    rows = []
    for index, metric in enumerate(METRICS):
        selected = arrays["valid"][:, index]
        prediction = arrays["mean"][selected, index].astype(np.float64)
        target = arrays["target"][selected, index].astype(np.float64)
        error = prediction - target
        rows.append({"sample_id": sample, "dataset": dataset, "fold": fold,
                     "metric": metric, "valid_residues": int(selected.sum()),
                     "mae": float(np.abs(error).mean()) if len(error) else None,
                     "rmse": float(np.sqrt(np.square(error).mean())) if len(error) else None,
                     "pearson": _correlation(prediction, target) if len(error) else None,
                     "spearman": _spearman(prediction, target) if len(error) else None})
    return rows


def write_csv(path: Path, rows: list) -> None:
    if not rows:
        return
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    atomic_bytes(path, buffer.getvalue().encode("utf-8"))


def implementation_hash() -> str:
    paths = [Path(__file__)] + sorted((ROOT / "src" / "flexppi").rglob("*.py"))
    return digest_json({str(p.relative_to(ROOT)): sha256(p) for p in paths})


def load_fold(root: Path, fold: int, device: torch.device, infer: bool) -> tuple:
    folder = root / "artifacts" / "performance_first_v2" / f"fold_{fold}"
    checkpoint_path = folder / "epoch_040.pt"
    resolved = json.loads((folder / "resolved_config.json").read_text(encoding="utf-8"))
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    if checkpoint.get("epoch") != 40:
        raise ValueError(f"Fold {fold} checkpoint is not epoch 40")
    config = ModelConfig(**checkpoint["config"]["model"])
    if (list(config.metrics) != list(METRICS) or config.use_condition_features
            or config.enable_response_head or config.sequence_backend != "learned"):
        raise ValueError(f"Fold {fold} is not the frozen two-output, structure-only model")
    if checkpoint["config"]["model"] != resolved["model"]:
        raise ValueError(f"Fold {fold} checkpoint and resolved model configs disagree")
    model = None
    if infer:
        model = FlexPPI(config)
        model.load_state_dict(checkpoint["model"], strict=True)
        model.to(device).eval()
    return model, resolved, checkpoint_path, sha256(checkpoint_path)


def prepare_input(raw: dict, mappings: Mapping[str, str]) -> tuple:
    record = StateRecord.from_dict(raw)
    structure_path = resolve_existing(record.structure_path, mappings)
    label_path = resolve_existing(str(record.label_path), mappings)
    record.structure_path, record.label_path = str(structure_path), str(label_path)
    structure = parse_pdb(str(structure_path))
    values, variances, mask = load_label_file(str(label_path))
    n = structure.num_residues
    if n > 1024 or n != values.shape[0]:
        raise ValueError(f"Length/label mismatch or >1024 residues: {record.sample_id}")
    declared = record.metadata.get("residue_count")
    if declared is not None and int(declared) != n:
        raise ValueError(f"Metadata residue count disagrees: {record.sample_id}")
    if not set(record.target_chains).issubset(set(structure.unique_chains())):
        raise ValueError(f"Requested target chain missing: {record.sample_id}")
    keys = list(zip(structure.chain_ids, structure.author_residue_ids))
    if len(set(keys)) != n:
        raise ValueError(f"Duplicate author residue keys: {record.sample_id}")
    if not np.isfinite(values[mask]).all() or not np.isfinite(variances[mask]).all():
        raise ValueError(f"Non-finite masked labels: {record.sample_id}")
    if np.any(variances[mask] < 0):
        raise ValueError(f"Negative masked label variance: {record.sample_id}")
    return record, structure, (values, variances, mask), structure_path, label_path


def predict_arrays(model: FlexPPI, device: torch.device, prepared: tuple,
                   identity: dict, fold: int) -> dict:
    record, structure, labels, structure_path, label_path = prepared
    batch = move_to_device(collate_states([load_state_example(record)]), device)
    with torch.inference_mode():
        output = model.forward_state(batch)
    indices = [METRIC_NAMES.index(metric) for metric in METRICS]
    target, variance, label_mask = labels
    output_mask = output["mask"][0].cpu().numpy().astype(bool)
    target_mask = batch["target_mask"][0].cpu().numpy().astype(bool)
    valid = label_mask[:, indices] & target_mask[:, None]
    valid &= structure.frame_mask[:, None] & output_mask[:, None]
    metadata = {"identity": identity, "status": "complete", "sample_id": record.sample_id,
                "dataset": record.dataset, "fold": fold, "context": "full",
                "residue_count": structure.num_residues, "structure_path": str(structure_path),
                "label_path": str(label_path), "biological_system": record.metadata.get(
                    "biological_system", record.sample_id), "protocol": PROTOCOL}
    return {"metadata_json": np.asarray(canonical_json(metadata)),
            "metric_names": np.asarray(METRICS), "chain": np.asarray(structure.chain_ids),
            "author_residue_id": np.asarray(structure.author_residue_ids),
            "residue_index": np.arange(structure.num_residues, dtype=np.int64),
            "aatype": structure.aatype, "amino_acid": np.asarray(list(structure.sequence)),
            "target_mask": target_mask, "geometry_mask": structure.frame_mask,
            "output_mask": output_mask, "label_mask": label_mask[:, indices],
            "target": target[:, indices], "label_variance": variance[:, indices],
            "mean": output["mean"][0].cpu().numpy(),
            "model_variance": output["log_variance"][0].exp().cpu().numpy(), "valid": valid}


def run(args: argparse.Namespace) -> dict:
    root, out = Path(args.root).resolve(), Path(args.output_dir).resolve()
    torch.set_num_threads(args.threads)
    torch.set_num_interop_threads(1)
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable; no automatic device substitution")
    inventory = audit_inventory(root / SPLIT_DIRECTORY)
    selected_folds = sorted(set(args.fold if args.fold is not None else range(5)))
    partial = len(selected_folds) != 5 or args.limit is not None
    start = time.monotonic()
    summary = {"cache_version": CACHE_VERSION, "protocol": PROTOCOL,
               "status": "running", "partial": partial, "validate_only": args.validate_only,
               "device": str(device), "torch_version": torch.__version__, "threads": args.threads,
               "expected_samples": 3068, "selected_folds": selected_folds,
               "limit_per_fold": args.limit, "source_hashes": inventory["source_hashes"],
               "manifest_audit": "passed", "processed_samples": 0, "cache_reused": 0,
               "valid_residue_counts": {metric: 0 for metric in METRICS}}
    implementation = implementation_hash()
    manifest = {sample: {"sample_id": sample, "dataset": raw["dataset"],
                        "fold": inventory["assignments"][sample], "status": "not_requested",
                        "record": raw} for sample, raw in inventory["ready"].items()}
    fold_rows = {}

    def publish() -> None:
        summary["elapsed_seconds"] = time.monotonic() - start
        atomic_json(out / "run_summary.json", summary)
        atomic_bytes(out / "prediction_manifest.jsonl", ("\n".join(
            canonical_json(manifest[sample]) for sample in sorted(manifest)
        ) + "\n").encode("utf-8"))
        for fold, rows in fold_rows.items():
            write_csv(out / f"fold_{fold}_per_system.csv", rows)

    try:
        publish()
        for fold in selected_folds:
            model, config, checkpoint_path, checkpoint_sha = load_fold(
                root, fold, device, not args.validate_only)
            fold_rows[fold] = []
            records = list(inventory["folds"][fold].items())
            if args.limit is not None:
                records = records[:args.limit]
            for sample, raw in records:
                prepared = prepare_input(raw, config["data"].get("path_mappings", {}))
                record, structure, labels, structure_path, label_path = prepared
                identity = {"cache_version": CACHE_VERSION, "sample_id": sample, "fold": fold,
                            "checkpoint_sha256": checkpoint_sha,
                            "input_sha256": sha256(structure_path), "label_sha256": sha256(label_path),
                            "record_sha256": digest_json(raw), "protocol_sha256": digest_json(PROTOCOL),
                            "implementation_sha256": implementation,
                            "source_hashes": inventory["source_hashes"],
                            "model_config_sha256": digest_json(config["model"])}
                target = np.ones(structure.num_residues, dtype=bool)
                if record.target_chains:
                    target = np.isin(structure.chain_ids, record.target_chains)
                indices = [METRIC_NAMES.index(metric) for metric in METRICS]
                audit_valid = labels[2][:, indices] & structure.frame_mask[:, None] & target[:, None]
                cache_path = out / f"fold_{fold}" / cache_name(sample)
                if args.validate_only:
                    status, valid = "validated_only", audit_valid
                else:
                    if cache_path.exists():
                        arrays = load_cache(cache_path, identity)
                        summary["cache_reused"] += 1
                    else:
                        arrays = predict_arrays(model, device, prepared, identity, fold)
                        save_cache(cache_path, arrays, identity)
                    if not np.array_equal(audit_valid, arrays["valid"]):
                        raise ValueError(f"Model silently lost valid labelled residues: {sample}")
                    valid = arrays["valid"]
                    status = "complete"
                    fold_rows[fold].extend(metric_rows(arrays, sample, raw["dataset"], fold))
                manifest[sample].update({"status": status, "identity": identity,
                                         "checkpoint": str(checkpoint_path),
                                         "cache_path": str(cache_path) if not args.validate_only else None,
                                         "residue_count": structure.num_residues,
                                         "valid_residue_counts": dict(zip(METRICS, map(int, valid.sum(0))))})
                summary["processed_samples"] += 1
                for index, metric in enumerate(METRICS):
                    summary["valid_residue_counts"][metric] += int(valid[:, index].sum())
                if summary["processed_samples"] % 25 == 0:
                    print(f"Processed {summary['processed_samples']} samples; "
                          f"elapsed {time.monotonic() - start:.1f}s", flush=True)
                    publish()
            del model
            publish()
        if not partial and summary["processed_samples"] != 3068:
            raise ValueError("Full run did not cover all 3068 samples")
        if not partial and summary["valid_residue_counts"] != {
                "rmsf": 959347, "chi_entropy": 431023}:
            raise ValueError("Full-run label counts disagree with frozen training/evaluation audit")
        summary["status"] = ("partial_validation" if partial else "validation_complete") if (
            args.validate_only) else ("partial_smoke" if partial else "complete")
        summary["source_sample_counts"] = dict(Counter(r["dataset"] for r in inventory["ready"].values()))
        publish()
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
        return summary
    except Exception as error:
        summary["status"], summary["error"] = "failed", repr(error)
        publish()
        raise


def main(argv: Optional[list] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--output-dir", default=str(ROOT / "artifacts/evaluation_v2/oof"))
    parser.add_argument("--fold", type=int, choices=range(5), action="append")
    parser.add_argument("--limit", type=int, help="Smoke only: maximum samples per requested fold")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--threads", type=int, default=6)
    args = parser.parse_args(argv)
    if args.threads < 1 or (args.limit is not None and args.limit < 1):
        parser.error("threads and limit must be positive")
    run(args)


if __name__ == "__main__":
    main()
