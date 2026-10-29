"""Export audited epoch-40 OOF predictions for ATLAS-Hchi v1, never old caches.

The frozen v2 exporter's input/mask/cache primitives are reused without editing
that script. This adapter has its own run loop, label counts and provenance.
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import export_v2_oof_predictions as shared
from flexppi.evaluation import _regression_report

ROOT = Path(__file__).resolve().parents[1]
SPLIT_DIRECTORY = "data/splits/state_rmsf_hchi_v2_atlas_hchi_v1"
CHECKPOINT_DIRECTORY = "artifacts/performance_first_v2_atlas_hchi_v1"
PAYLOAD_DIRECTORY = "output/atlas_hchi_retraining_v1_linux/payload"
RELEASE_DIRECTORY = "artifacts/data/atlas_rmsf_hchi_v1"
EXPECTED_COUNTS = {"rmsf": 959347, "chi_entropy": 798385}
EXPECTED_SOURCES = {"atlas": 1929, "dynarepo": 895, "itrace": 244}
PROTOCOL = {**shared.PROTOCOL, "name": "atlas_hchi_v1_final_epoch40_oof",
            "split_version": "state_rmsf_hchi_v2_atlas_hchi_v1",
            "checkpoint_directory": CHECKPOINT_DIRECTORY,
            "valid_residue_counts": EXPECTED_COUNTS,
            "checkpoint_loading": "weights_only=True; explicit standard numpy RNG allowlist"}
DEFAULT_MAPPINGS = {
    "/ranze/repos/FlexPPI/artifacts/data": str(ROOT / "artifacts/data"),
    "/ranze/data/ATLAS_extracted": "E:/data/protein_md_datasets/ATLAS_extracted",
    "/ranze/data/dynarepo/dynarepo_non_nucleic_full": "D:/data/dynarepo_non_nucleic_full",
    "/ranze/data/ITrace": "E:/data/ITrace",
}


def safe_load_checkpoint(path):
    # Numpy 2 pickles use numpy._core; local numpy 1.26 exposes the same standard
    # reconstruct callable under numpy.core. Allow both explicit qualified names.
    reconstruct = np.core.multiarray._reconstruct
    allowed = [np.ndarray, np.dtype, type(np.dtype("uint32")),
               (reconstruct, "numpy._core.multiarray._reconstruct"),
               (reconstruct, "numpy.core.multiarray._reconstruct")]
    permitted = {"numpy.ndarray", "numpy.dtype", "numpy.dtypes.UInt32DType",
                 "numpy._core.multiarray._reconstruct", "numpy.core.multiarray._reconstruct"}
    unknown = set(torch.serialization.get_unsafe_globals_in_checkpoint(path)) - permitted
    if unknown:
        raise ValueError("Unexpected checkpoint globals: " + repr(sorted(unknown)))
    with torch.serialization.safe_globals(allowed):
        return torch.load(path, map_location="cpu", weights_only=True)


def validate_configuration(checkpoint, resolved, fold):
    if checkpoint.get("epoch") != 40:
        raise ValueError("Checkpoint must be epoch 40")
    if checkpoint.get("config") != resolved:
        raise ValueError("Checkpoint and resolved configurations differ")
    data = resolved["data"]
    for field, role in (("state_manifest", "train"), ("test_state_manifest", "test")):
        expected = f"{SPLIT_DIRECTORY}/fold_{fold}_{role}.jsonl"
        actual = str(data.get(field, "")).replace("\\", "/")
        if actual != expected and not actual.endswith("/" + expected):
            raise ValueError("Wrong checkpoint fold/split: " + field)
    if (resolved["stage"] != "pretrain" or resolved["init_checkpoint"] is not None
            or data.get("validation_state_manifest") is not None
            or resolved["train"]["val_fraction"] != 0
            or resolved["train"]["epochs"] != 40):
        raise ValueError("Not the frozen from-scratch, no-validation 40-epoch protocol")
    model = shared.ModelConfig(**resolved["model"])
    if (list(model.metrics) != list(shared.METRICS) or model.use_condition_features
            or model.enable_response_head or model.sequence_backend != "learned"):
        raise ValueError("Not the frozen two-output structure-only model")
    return model


def load_fold(root, fold, device, infer):
    directory = root / CHECKPOINT_DIRECTORY / f"fold_{fold}"
    path = directory / "epoch_040.pt"
    resolved = json.loads((directory / "resolved_config.json").read_text(encoding="utf-8"))
    checkpoint = safe_load_checkpoint(path)
    model_config = validate_configuration(checkpoint, resolved, fold)
    audit = json.loads((root / "output/flexppi_atlas_hchi_v1_five_fold_reference/checkpoint_audit.json").read_text())
    reference = next(row for row in audit["folds"] if row["fold_index"] == fold)
    checksum = shared.sha256(path)
    if checksum != reference["checkpoint_sha256"]:
        raise ValueError("Checkpoint differs from independently audited new run")
    if not all(torch.isfinite(t).all() for t in checkpoint["model"].values()):
        raise ValueError("Nonfinite model tensors")
    model = None
    if infer:
        model = shared.FlexPPI(model_config)
        model.load_state_dict(checkpoint["model"], strict=True)
        model.to(device).eval()
    return model, resolved, path, checksum


def audit_payload(root):
    payload = root / PAYLOAD_DIRECTORY
    checksums = {}
    for line in (payload / RELEASE_DIRECTORY / "PAYLOAD_SHA256SUMS").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        checksums[name.lstrip("*")] = expected
    required = [f"{SPLIT_DIRECTORY}/all_ready.jsonl",
                f"{RELEASE_DIRECTORY}/training_inputs.json"]
    required += [f"{SPLIT_DIRECTORY}/fold_{k}_{role}.jsonl"
                 for k in range(5) for role in ("train", "test")]
    for name in required:
        if name not in checksums or shared.sha256(payload / name) != checksums[name]:
            raise ValueError("Frozen payload checksum mismatch: " + name)
    inventory = shared.audit_inventory(payload / SPLIT_DIRECTORY)
    if dict(Counter(r["dataset"] for r in inventory["ready"].values())) != EXPECTED_SOURCES:
        raise ValueError("Unexpected source composition")
    inputs = shared.unique_records(json.loads(
        (payload / RELEASE_DIRECTORY / "training_inputs.json").read_text()), "training inputs")
    if set(inputs) != set(inventory["ready"]):
        raise ValueError("Training input inventory membership differs")
    for sample, record in inventory["ready"].items():
        if record["dataset"] == "atlas" and "/atlas_rmsf_hchi_v1/labels/" not in record["label_path"]:
            raise ValueError("Old ATLAS labels requested: " + sample)
    return inventory, inputs


def implementation_hash():
    paths = [Path(__file__).resolve(), Path(shared.__file__).resolve()]
    paths += sorted((ROOT / "src/flexppi").rglob("*.py"))
    return shared.digest_json({str(p.relative_to(ROOT)): shared.sha256(p) for p in paths})


def select_records(records, limit):
    if limit is None:
        return records
    # Smoke tests cover all three sources before filling any remaining slots.
    first = {}
    for sample, record in records:
        first.setdefault(record["dataset"], (sample, record))
    chosen = list(first.values())[:limit]
    chosen_ids = {x[0] for x in chosen}
    return (chosen + [x for x in records if x[0] not in chosen_ids])[:limit]


def verify_server(fold, buckets, root, output, tolerance=2e-4):
    reference_path = root / CHECKPOINT_DIRECTORY / f"evaluation/fold_{fold}_test.json"
    reference = json.loads(reference_path.read_text())
    reports, comparisons = {}, []
    for metric, items in buckets.items():
        p = np.concatenate([x[1] for x in items]).astype(np.float64)
        y = np.concatenate([x[2] for x in items]).astype(np.float64)
        groups = np.concatenate([np.repeat(x[0], len(x[1])) for x in items])
        sources = np.concatenate([np.repeat(x[3], len(x[1])) for x in items])
        reports[metric] = {}
        for source in ["all", *EXPECTED_SOURCES]:
            selected = np.ones(len(y), bool) if source == "all" else sources == source
            report = _regression_report(p[selected], y[selected], groups[selected])
            reports[metric][source] = report
            ref = reference["metrics"][metric]["overall"] if source == "all" else reference["metrics"][metric]["by_dataset"][source]
            if report["count"] != ref["count"] or report["system_count"] != ref["system_count"]:
                raise ValueError(f"CPU/server count mismatch: {fold} {metric} {source}")
            for scope in ["pooled", "macro_by_system"]:
                local, remote = (report, ref) if scope == "pooled" else (report[scope], ref[scope])
                for statistic in ["mae", "rmse", "pearson", "spearman"]:
                    if local[statistic] is None or remote[statistic] is None:
                        difference = None
                        passed = local[statistic] == remote[statistic]
                    else:
                        difference = float(local[statistic] - remote[statistic])
                        passed = abs(difference) <= tolerance
                    comparisons.append(dict(fold=fold, metric=metric, source=source,
                                            scope=scope, statistic=statistic,
                                            cpu_minus_server=difference, passed=passed))
    result = dict(fold=fold, server_sha256=shared.sha256(reference_path), tolerance=tolerance,
                  explanation="CPU float32 inference versus server CUDA; tolerance fixed before inference.",
                  passed=all(x["passed"] for x in comparisons), comparisons=comparisons)
    shared.atomic_json(output / f"fold_{fold}_recomputed_metrics.json", reports)
    shared.atomic_json(output / f"fold_{fold}_server_reproduction.json", result)
    if not result["passed"]:
        raise ValueError(f"CPU/server numerical mismatch in fold {fold}; see reproduction report")
    return result


def run(args):
    root, output = Path(args.root).resolve(), Path(args.output_dir).resolve()
    allowed = (root / "artifacts/evaluation_atlas_hchi_v1").resolve()
    if allowed not in output.parents:
        raise ValueError("Output must be a child of artifacts/evaluation_atlas_hchi_v1")
    torch.set_num_threads(args.threads)
    torch.set_num_interop_threads(1)
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA unavailable")
    shared.PROTOCOL = copy.deepcopy(PROTOCOL)
    inventory, frozen_inputs = audit_payload(root)
    implementation = implementation_hash()
    folds = sorted(set(args.fold if args.fold is not None else range(5)))
    partial = len(folds) != 5 or args.limit is not None
    start = time.monotonic()
    summary = dict(cache_version=shared.CACHE_VERSION, protocol=PROTOCOL, status="running",
                   partial=partial, validate_only=args.validate_only, expected_samples=3068,
                   started_at_utc=datetime.now(timezone.utc).isoformat(),
                   device=str(device), torch_version=str(torch.__version__), threads=args.threads,
                   selected_folds=folds, limit_per_fold=args.limit,
                   smoke_selection="first sample from each source, then manifest order",
                   source_hashes=inventory["source_hashes"], manifest_audit="passed",
                   processed_samples=0, cache_reused=0,
                   valid_residue_counts={metric: 0 for metric in shared.METRICS},
                   expected_valid_residue_counts=EXPECTED_COUNTS,
                   source_sample_counts=EXPECTED_SOURCES)
    manifest = {sid: dict(sample_id=sid, dataset=r["dataset"], record=r,
                         fold=inventory["assignments"][sid], status="not_requested")
                for sid, r in inventory["ready"].items()}
    rows, reproduction = {}, []

    def publish():
        summary["elapsed_seconds"] = time.monotonic() - start
        shared.atomic_json(output / "run_summary.json", summary)
        shared.atomic_bytes(output / "prediction_manifest.jsonl", ("\n".join(
            shared.canonical_json(manifest[sid]) for sid in sorted(manifest)) + "\n").encode("utf-8"))
        for fold, entries in rows.items():
            shared.write_csv(output / f"fold_{fold}_per_system.csv", entries)

    try:
        publish()
        for fold in folds:
            model, config, checkpoint_path, checkpoint_sha = load_fold(root, fold, device, not args.validate_only)
            rows[fold], buckets = [], {metric: [] for metric in shared.METRICS}
            mappings = {**config["data"].get("path_mappings", {}), **DEFAULT_MAPPINGS}
            for sample, raw in select_records(list(inventory["folds"][fold].items()), args.limit):
                prepared = shared.prepare_input(raw, mappings)
                record, structure, labels, structure_path, label_path = prepared
                input_sha, label_sha = shared.sha256(structure_path), shared.sha256(label_path)
                frozen = frozen_inputs[sample]
                if (input_sha != frozen["structure_sha256"] or label_sha != frozen["label_sha256"]
                        or frozen["fold"] != fold or frozen["dataset"] != raw["dataset"]):
                    raise ValueError("Input differs from frozen training payload: " + sample)
                identity = dict(cache_version=shared.CACHE_VERSION, sample_id=sample, fold=fold,
                                checkpoint_sha256=checkpoint_sha, input_sha256=input_sha,
                                label_sha256=label_sha, record_sha256=shared.digest_json(raw),
                                protocol_sha256=shared.digest_json(PROTOCOL), implementation_sha256=implementation,
                                source_hashes=inventory["source_hashes"],
                                model_config_sha256=shared.digest_json(config["model"]))
                target = np.ones(structure.num_residues, bool)
                if record.target_chains:
                    target = np.isin(structure.chain_ids, record.target_chains)
                indices = [shared.METRIC_NAMES.index(metric) for metric in shared.METRICS]
                valid = labels[2][:, indices] & structure.frame_mask[:, None] & target[:, None]
                cache_path = output / f"fold_{fold}" / shared.cache_name(sample)
                if args.validate_only:
                    status = "validated_only"
                else:
                    if cache_path.exists():
                        arrays = shared.load_cache(cache_path, identity)
                        summary["cache_reused"] += 1
                    else:
                        arrays = shared.predict_arrays(model, device, prepared, identity, fold)
                        shared.save_cache(cache_path, arrays, identity)
                    if not np.array_equal(valid, arrays["valid"]):
                        raise ValueError("Model silently lost valid residues: " + sample)
                    rows[fold].extend(shared.metric_rows(arrays, sample, raw["dataset"], fold))
                    for j, metric in enumerate(shared.METRICS):
                        keep = arrays["valid"][:, j]
                        buckets[metric].append((sample, arrays["mean"][keep, j], arrays["target"][keep, j], raw["dataset"]))
                    status = "complete"
                manifest[sample].update(status=status, identity=identity, checkpoint=str(checkpoint_path),
                                        cache_path=None if args.validate_only else str(cache_path),
                                        residue_count=structure.num_residues,
                                        valid_residue_counts=dict(zip(shared.METRICS, map(int, valid.sum(0)))))
                summary["processed_samples"] += 1
                for j, metric in enumerate(shared.METRICS):
                    summary["valid_residue_counts"][metric] += int(valid[:, j].sum())
                if summary["processed_samples"] % 25 == 0:
                    print(f"Processed {summary['processed_samples']}/3068; elapsed {time.monotonic()-start:.1f}s", flush=True)
                    publish()
            del model
            if not args.validate_only and args.limit is None:
                reproduction.append(verify_server(fold, buckets, root, output))
                shared.atomic_json(output / "server_reproduction.json", dict(
                    passed=all(x["passed"] for x in reproduction), folds=reproduction))
            publish()
        if not partial and (summary["processed_samples"] != 3068 or summary["valid_residue_counts"] != EXPECTED_COUNTS):
            raise ValueError("Full-run coverage differs from frozen new-run counts")
        summary["status"] = (("partial_validation" if partial else "validation_complete") if args.validate_only
                             else ("partial_smoke" if partial else "complete"))
        summary["server_reproduction_passed"] = all(x["passed"] for x in reproduction) if reproduction else None
        publish()
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
        return summary
    except Exception as error:
        summary.update(status="failed", error=repr(error))
        publish()
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--output-dir", default=str(ROOT / "artifacts/evaluation_atlas_hchi_v1/oof"))
    parser.add_argument("--fold", type=int, choices=range(5), action="append")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--threads", type=int, default=6)
    args = parser.parse_args()
    if args.threads < 1 or (args.limit is not None and args.limit < 1):
        parser.error("threads/limit must be positive")
    run(args)


if __name__ == "__main__":
    main()
