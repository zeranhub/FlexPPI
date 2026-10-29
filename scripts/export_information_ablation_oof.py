"""Export provenance-checked OOF predictions for the information-ablation arms.

One held-out fold checkpoint per sample; no ensemble, no checkpoint selection.
The frozen exporter's cache format, identity and validation primitives are
reused so ablation caches cannot be mixed with the reference run.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

import export_v2_oof_predictions as shared  # noqa: E402
from flexppi.constants import METRIC_NAMES  # noqa: E402
from flexppi.evaluation import _regression_report  # noqa: E402
from flexppi.model.information_ablation import build_information_ablation  # noqa: E402

SPLIT_DIRECTORY = "data/splits/state_rmsf_hchi_v2_atlas_hchi_v1"
EXPECTED_COUNTS = {"rmsf": 959347, "chi_entropy": 798385}
EXPECTED_SOURCES = {"atlas": 1929, "dynarepo": 895, "itrace": 244}
VARIANTS = ("full", "sequence", "geometry", "intrachain")
CHECKPOINT_DIRECTORIES = {
    "full": "artifacts/performance_first_v2_atlas_hchi_v1",
    "sequence": "artifacts/information_ablation/sequence",
    "geometry": "artifacts/information_ablation/geometry",
    "intrachain": "artifacts/information_ablation/intrachain",
}
BASE_PROTOCOL = {
    **shared.PROTOCOL,
    "name": "information_ablation_final_epoch40_oof",
    "split_version": "state_rmsf_hchi_v2_atlas_hchi_v1",
    "selection": "each sample uses only its held-out fold, no best.pt or ensemble",
    "valid_residue_counts": EXPECTED_COUNTS,
    "checkpoint_loading": "weights_only=True; explicit standard numpy RNG allowlist",
}
SMOKE_PROTOCOL = {**BASE_PROTOCOL, "name": "information_ablation_smoke_oof",
                  "selection": "single-epoch smoke checkpoint; never for reported results"}


def protocol_for(variant: str, smoke: bool) -> dict:
    protocol = SMOKE_PROTOCOL if smoke else BASE_PROTOCOL
    return {**protocol, "variant": variant}


def safe_load_checkpoint(path: Path) -> dict:
    reconstruct = np.core.multiarray._reconstruct
    allowed = [np.ndarray, np.dtype, type(np.dtype("uint32")),
               (reconstruct, "numpy._core.multiarray._reconstruct"),
               (reconstruct, "numpy.core.multiarray._reconstruct")]
    permitted = {"numpy.ndarray", "numpy.dtype", "numpy.dtypes.UInt32DType",
                 "numpy._core.multiarray._reconstruct", "numpy.core.multiarray._reconstruct"}
    unknown = set(torch.serialization.get_unsafe_globals_in_checkpoint(str(path))) - permitted
    if unknown:
        raise ValueError("Unexpected checkpoint globals: " + repr(sorted(unknown)))
    with torch.serialization.safe_globals(allowed):
        return torch.load(path, map_location="cpu", weights_only=True)


def validate_configuration(resolved: dict, fold: int, variant: str, epoch: int,
                           smoke: bool) -> None:
    if resolved["model"].get("information_ablation", "full") != variant:
        raise ValueError("Checkpoint was not trained as arm " + variant)
    data = resolved["data"]
    for field, role in (("state_manifest", "train"), ("test_state_manifest", "test")):
        expected = f"{SPLIT_DIRECTORY}/fold_{fold}_{role}.jsonl"
        actual = str(data.get(field, "")).replace("\\", "/")
        if actual != expected and not actual.endswith("/" + expected):
            raise ValueError("Wrong checkpoint fold/split: " + field)
    model = shared.ModelConfig(**resolved["model"])
    if (list(model.metrics) != list(shared.METRICS) or model.use_condition_features
            or model.enable_response_head or model.sequence_backend != "learned"):
        raise ValueError("Not the frozen two-output structure-only model")
    if smoke:
        return
    if (resolved["stage"] != "pretrain" or resolved["init_checkpoint"] is not None
            or data.get("validation_state_manifest") is not None
            or resolved["train"]["val_fraction"] != 0
            or resolved["train"]["epochs"] != epoch):
        raise ValueError("Not the frozen from-scratch, no-validation %d-epoch protocol" % epoch)


def load_fold(root: Path, fold: int, variant: str, epoch: int, smoke: bool,
              device: torch.device, infer: bool, checkpoint_dir: str = None) -> tuple:
    base = Path(checkpoint_dir) if checkpoint_dir else Path(CHECKPOINT_DIRECTORIES[variant])
    if not base.is_absolute():
        base = root / base
    directory = base / ("fold_%d" % fold)
    path = directory / ("epoch_%03d.pt" % epoch)
    resolved = json.loads((directory / "resolved_config.json").read_text(encoding="utf-8"))
    checkpoint = safe_load_checkpoint(path)
    if not smoke and checkpoint.get("epoch") != epoch:
        raise ValueError("Checkpoint is not epoch %d" % epoch)
    validate_configuration(resolved, fold, variant, epoch, smoke)
    if checkpoint.get("config") != resolved:
        raise ValueError("Checkpoint and resolved configurations differ")
    if not all(torch.isfinite(tensor).all() for tensor in checkpoint["model"].values()):
        raise ValueError("Nonfinite model tensors")
    model = None
    if infer:
        model = build_information_ablation(shared.ModelConfig(**resolved["model"]), variant)
        model.load_state_dict(checkpoint["model"], strict=True)
        model.to(device).eval()
    return model, resolved, path, shared.sha256(path)


def predict_arrays(model, device: torch.device, prepared: tuple, identity: dict,
                   fold: int, variant: str, smoke: bool) -> dict:
    record, structure, labels, structure_path, label_path = prepared
    batch = shared.move_to_device(shared.collate_states([shared.load_state_example(record)]), device)
    with torch.inference_mode():
        output = model.forward_state(batch)
    indices = [METRIC_NAMES.index(metric) for metric in shared.METRICS]
    target, variance, label_mask = labels
    output_mask = output["mask"][0].cpu().numpy().astype(bool)
    target_mask = batch["target_mask"][0].cpu().numpy().astype(bool)
    valid = label_mask[:, indices] & target_mask[:, None]
    valid &= structure.frame_mask[:, None] & output_mask[:, None]
    metadata = {"identity": identity, "status": "complete", "sample_id": record.sample_id,
                "dataset": record.dataset, "fold": fold, "context": variant,
                "smoke": smoke, "residue_count": structure.num_residues,
                "structure_path": str(structure_path), "label_path": str(label_path),
                "biological_system": record.metadata.get("biological_system", record.sample_id),
                "protocol": protocol_for(variant, smoke)}
    return {"metadata_json": np.asarray(shared.canonical_json(metadata)),
            "metric_names": np.asarray(shared.METRICS), "chain": np.asarray(structure.chain_ids),
            "author_residue_id": np.asarray(structure.author_residue_ids),
            "residue_index": np.arange(structure.num_residues, dtype=np.int64),
            "aatype": structure.aatype, "amino_acid": np.asarray(list(structure.sequence)),
            "target_mask": target_mask, "geometry_mask": structure.frame_mask,
            "output_mask": output_mask, "label_mask": label_mask[:, indices],
            "target": target[:, indices], "label_variance": variance[:, indices],
            "mean": output["mean"][0].cpu().numpy(),
            "model_variance": output["log_variance"][0].exp().cpu().numpy(), "valid": valid}


def verify_server(fold: int, buckets: dict, root: Path, output: Path, variant: str,
                  tolerance: float = 2e-4) -> dict:
    reference_path = root / CHECKPOINT_DIRECTORIES[variant] / f"evaluation/fold_{fold}_test.json"
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    reports, comparisons = {}, []
    for metric, items in buckets.items():
        prediction = np.concatenate([entry[1] for entry in items]).astype(np.float64)
        target = np.concatenate([entry[2] for entry in items]).astype(np.float64)
        groups = np.concatenate([np.repeat(entry[0], len(entry[1])) for entry in items])
        sources = np.concatenate([np.repeat(entry[3], len(entry[1])) for entry in items])
        reports[metric] = {}
        for source in ["all", *EXPECTED_SOURCES]:
            selected = np.ones(len(target), bool) if source == "all" else sources == source
            report = _regression_report(prediction[selected], target[selected], groups[selected])
            reports[metric][source] = report
            ref = (reference["metrics"][metric]["overall"] if source == "all"
                   else reference["metrics"][metric]["by_dataset"][source])
            if report["count"] != ref["count"] or report["system_count"] != ref["system_count"]:
                raise ValueError(f"CPU/server count mismatch: fold {fold} {metric} {source}")
            for scope in ["pooled", "macro_by_system"]:
                local, remote = (report, ref) if scope == "pooled" else (report[scope], ref[scope])
                for statistic in ["mae", "rmse", "pearson", "spearman"]:
                    if local[statistic] is None or remote[statistic] is None:
                        difference, passed = None, local[statistic] == remote[statistic]
                    else:
                        difference = float(local[statistic] - remote[statistic])
                        passed = abs(difference) <= tolerance
                    comparisons.append(dict(fold=fold, variant=variant, metric=metric,
                                            source=source, scope=scope, statistic=statistic,
                                            cpu_minus_server=difference, passed=passed))
    result = dict(fold=fold, variant=variant, server_sha256=shared.sha256(reference_path),
                  tolerance=tolerance,
                  explanation="CPU float32 inference versus server CUDA; tolerance fixed before inference.",
                  passed=all(entry["passed"] for entry in comparisons), comparisons=comparisons)
    shared.atomic_json(output / f"fold_{fold}_recomputed_metrics.json", reports)
    shared.atomic_json(output / f"fold_{fold}_server_reproduction.json", result)
    if not result["passed"]:
        raise ValueError(f"CPU/server numerical mismatch in fold {fold}; see reproduction report")
    return result


def implementation_hash(script: Path) -> str:
    paths = [script] + sorted((ROOT / "src" / "flexppi").rglob("*.py"))
    return shared.digest_json({str(p.relative_to(ROOT)): shared.sha256(p) for p in paths})


def run(args: argparse.Namespace) -> dict:
    root = Path(args.root).resolve()
    if args.output_dir:
        output = Path(args.output_dir).resolve()
    else:
        suffix = "output/information_ablation_oof_smoke" if args.smoke else "output/information_ablation_oof"
        output = root / suffix / args.variant
    torch.set_num_threads(args.threads)
    torch.set_num_interop_threads(1)
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA unavailable")
    protocol = protocol_for(args.variant, args.smoke)
    inventory = shared.audit_inventory(root / SPLIT_DIRECTORY)
    observed_sources = Counter(record["dataset"] for record in inventory["ready"].values())
    if dict(observed_sources) != EXPECTED_SOURCES:
        raise ValueError("Source composition differs from the frozen run: " + repr(observed_sources))
    implementation = implementation_hash(Path(__file__).resolve())
    folds = sorted(set(args.fold if args.fold is not None else range(5)))
    partial = args.smoke or len(folds) != 5 or args.limit is not None
    start = time.monotonic()
    summary = dict(cache_version=shared.CACHE_VERSION, protocol=protocol, variant=args.variant,
                   status="running", partial=partial, validate_only=args.validate_only,
                   device=str(device), torch_version=str(torch.__version__), threads=args.threads,
                   expected_samples=3068, started_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                   selected_folds=folds, limit_per_fold=args.limit,
                   source_hashes=inventory["source_hashes"], manifest_audit="passed",
                   processed_samples=0, cache_reused=0,
                   valid_residue_counts={metric: 0 for metric in shared.METRICS},
                   expected_valid_residue_counts=EXPECTED_COUNTS,
                   source_sample_counts=EXPECTED_SOURCES)
    manifest = {sample: dict(sample_id=sample, dataset=raw["dataset"],
                             fold=inventory["assignments"][sample], status="not_requested",
                             record=raw)
                for sample, raw in inventory["ready"].items()}
    rows, reproduction = {}, []

    def publish() -> None:
        summary["elapsed_seconds"] = time.monotonic() - start
        shared.atomic_json(output / "run_summary.json", summary)
        shared.atomic_bytes(output / "prediction_manifest.jsonl", ("\n".join(
            shared.canonical_json(manifest[sample]) for sample in sorted(manifest)) + "\n"
        ).encode("utf-8"))
        for fold, entries in rows.items():
            shared.write_csv(output / f"fold_{fold}_per_system.csv", entries)

    try:
        publish()
        for fold in folds:
            model, config, checkpoint_path, checkpoint_sha = load_fold(
                root, fold, args.variant, args.epoch, args.smoke, device,
                not args.validate_only, args.checkpoint_dir)
            rows[fold], buckets = [], {metric: [] for metric in shared.METRICS}
            mappings = {**config["data"].get("path_mappings", {})}
            records = list(inventory["folds"][fold].items())
            if args.limit is not None:
                records = records[:args.limit]
            for sample, raw in records:
                prepared = shared.prepare_input(raw, mappings)
                record, structure, labels, structure_path, label_path = prepared
                identity = dict(cache_version=shared.CACHE_VERSION, sample_id=sample, fold=fold,
                                variant=args.variant,
                                checkpoint_sha256=checkpoint_sha,
                                input_sha256=shared.sha256(structure_path),
                                label_sha256=shared.sha256(label_path),
                                record_sha256=shared.digest_json(raw),
                                protocol_sha256=shared.digest_json(protocol),
                                implementation_sha256=implementation,
                                source_hashes=inventory["source_hashes"],
                                model_config_sha256=shared.digest_json(config["model"]))
                target = np.ones(structure.num_residues, dtype=bool)
                if record.target_chains:
                    target = np.isin(structure.chain_ids, record.target_chains)
                indices = [METRIC_NAMES.index(metric) for metric in shared.METRICS]
                audit_valid = labels[2][:, indices] & structure.frame_mask[:, None] & target[:, None]
                cache_path = output / f"fold_{fold}" / shared.cache_name(sample)
                if args.validate_only:
                    valid, status = audit_valid, "validated_only"
                else:
                    if cache_path.exists():
                        arrays = shared.load_cache(cache_path, identity)
                        summary["cache_reused"] += 1
                    else:
                        arrays = predict_arrays(model, device, prepared, identity, fold,
                                                args.variant, args.smoke)
                        shared.save_cache(cache_path, arrays, identity)
                    if not np.array_equal(audit_valid, arrays["valid"]):
                        raise ValueError("Model silently lost valid residues: " + sample)
                    valid, status = arrays["valid"], "complete"
                    rows[fold].extend(shared.metric_rows(arrays, sample, raw["dataset"], fold))
                    for index, metric in enumerate(shared.METRICS):
                        keep = arrays["valid"][:, index]
                        buckets[metric].append((sample, arrays["mean"][keep, index],
                                                arrays["target"][keep, index], raw["dataset"]))
                manifest[sample].update(status=status, identity=identity,
                                        checkpoint=str(checkpoint_path),
                                        cache_path=None if args.validate_only else str(cache_path),
                                        residue_count=structure.num_residues,
                                        valid_residue_counts=dict(zip(
                                            shared.METRICS, map(int, valid.sum(0)))))
                summary["processed_samples"] += 1
                for index, metric in enumerate(shared.METRICS):
                    summary["valid_residue_counts"][metric] += int(valid[:, index].sum())
                if summary["processed_samples"] % 25 == 0:
                    print(f"Processed {summary['processed_samples']} samples; "
                          f"elapsed {time.monotonic() - start:.1f}s", flush=True)
                    publish()
            del model
            if args.verify and not args.validate_only and args.limit is None:
                reproduction.append(verify_server(fold, buckets, root, output, args.variant))
                shared.atomic_json(output / "server_reproduction.json", dict(
                    passed=all(entry["passed"] for entry in reproduction), folds=reproduction))
            publish()
        if not partial:
            if summary["processed_samples"] != 3068:
                raise ValueError("Full run did not cover all 3068 samples")
            if summary["valid_residue_counts"] != EXPECTED_COUNTS:
                raise ValueError("Full-run label counts disagree with frozen evaluation audit")
        summary["status"] = (("partial_validation" if partial else "validation_complete")
                             if args.validate_only else
                             ("partial_smoke" if partial else "complete"))
        summary["server_reproduction_passed"] = (all(entry["passed"] for entry in reproduction)
                                                 if reproduction else None)
        publish()
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
        return summary
    except Exception as error:
        summary.update(status="failed", error=repr(error))
        publish()
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--variant", required=True, choices=VARIANTS)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--fold", type=int, choices=range(5), action="append")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--epoch", type=int, default=40)
    parser.add_argument("--checkpoint-dir", default=None,
                        help="Override the per-variant checkpoint directory (used by smoke tests)")
    parser.add_argument("--smoke", action="store_true",
                        help="Smoke checkpoints: relaxes epoch/protocol checks, marks partial output")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--verify", action="store_true",
                        help="Compare recomputed fold metrics against the server evaluation JSON")
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--threads", type=int, default=6)
    args = parser.parse_args()
    if args.threads < 1 or (args.limit is not None and args.limit < 1):
        parser.error("threads and limit must be positive")
    run(args)


if __name__ == "__main__":
    main()
