"""Export the six already-planned fold-0 members on the same heldout records.

This collector does not train, tune, change manifests, or use the GPU. All six
epoch-40 checkpoints must be available before the complete export can start.
Run --audit-only --smoke-test to check preparation with available checkpoints.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ["CUDA_VISIBLE_DEVICES"] = ""

import numpy as np
import pandas as pd
import torch
import yaml
from scipy.stats import spearmanr, pearsonr

HERE = Path(__file__).resolve().parent
PAYLOAD = HERE / "training_payload"
sys.path.insert(0, str(PAYLOAD / "src"))
from flexppi.config import ModelConfig
from flexppi.constants import METRIC_NAMES
from flexppi.data.dataset import StateDataset, collate_states
from flexppi.data.structures import parse_pdb
from flexppi.model.information_ablation import build_information_ablation

FREEZE_SHA = "5b994bb59684e3a33c97f4a1801ee9a606e5dd856db9cfe23bf18515697707d3"
REMOTE_PREFIX = "/root/flexppi_revision_20261002/corrected_training_v1"
METRICS = ["rmsf", "chi_entropy"]
MEMBERS = [(v, s, v if s == 7 else f"{v}_seed{s}")
           for v in ["full", "intrachain"] for s in [7, 11, 23]]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def assert_equal(left, right, message):
    if left != right:
        raise ValueError(message)


def json_write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=True), encoding="utf-8")


def check_inputs():
    assert_equal(list(METRIC_NAMES[:2]), METRICS, "Dataset metric order changed")
    assert_equal(sha(PAYLOAD / "training_freeze.json"), FREEZE_SHA, "Freeze hash changed")
    freeze = json.loads((PAYLOAD / "training_freeze.json").read_text(encoding="utf-8"))
    for name, expected in freeze["execution_hashes"].items():
        path = PAYLOAD.joinpath(*name.replace("\\", "/").split("/"))
        assert_equal(sha(path), expected, "Execution or split hash changed: " + name)
    manifest = PAYLOAD / "splits/fold_0_test.jsonl"
    originals = [json.loads(line) for line in manifest.read_text(encoding="utf-8").splitlines()]
    training = [json.loads(line) for line in (PAYLOAD / "splits/fold_0_train.jsonl").read_text(encoding="utf-8").splitlines()]
    test_ids = [row["sample_id"] for row in originals]
    train_ids = [row["sample_id"] for row in training]
    assert_equal(len(test_ids), 612, "Unexpected fold-0 heldout size")
    assert_equal(len(train_ids), 2451, "Unexpected fold-0 training size")
    assert_equal(len(set(test_ids)), len(test_ids), "Duplicate heldout IDs")
    if set(test_ids) & set(train_ids):
        raise ValueError("Fold-0 training and heldout IDs overlap")
    dataset = StateDataset(str(manifest), max_residues=1024,
        path_mappings={REMOTE_PREFIX: str(PAYLOAD)})
    assert_equal([r.sample_id for r in dataset.records], test_ids, "Silent dataset filtering")
    input_hashes = {k.replace("\\", "/"): v for k, v in freeze["input_hashes"].items()}
    provenance = {}
    valid_counts = np.zeros(2, dtype=np.int64)
    schema = json.loads((HERE / "static_corrected/training_schema.json").read_text(encoding="utf-8"))
    static_records = {r["sample_id"]: r for r in schema["records"]}
    assert_equal(len(static_records), len(schema["records"]), "Duplicate static record IDs")
    assert_equal(set(r["sample_id"] for r in schema["records"] if r["fold"] == 0),
                 set(test_ids), "Static fold-0 heldout members differ")
    assert_equal(set(r["sample_id"] for r in schema["records"] if r["fold"] != 0),
                 set(train_ids), "Static fold-0 training members differ")
    static_root = HERE / "static_corrected"
    static_audit = json.loads((static_root / "run_audit.json").read_text(encoding="utf-8"))
    assert_equal(static_audit["status"], "complete", "Static fitting incomplete")
    assert_equal(static_audit["training_freeze_sha256"], FREEZE_SHA, "Static input freeze differs")
    static_source_paths = [static_root / name for name in ["training_features.npz", "training_schema.json",
        "run_audit.json", "static_sequence_complex_oof.npz", "aa_only_oof.npz"]]
    static_source_paths += [HERE / "static_feature_definitions.py", HERE / "run_corrected_static.py"]
    for job in static_audit["jobs"]:
        if job["fold"] == 0:
            checkpoint = static_root / f"static_fold0_{job['metric']}.joblib"
            assert_equal(sha(checkpoint), job["checkpoint_sha256"], "Static member-0 checkpoint changed")
            static_source_paths.append(checkpoint)
    static_predictions = {}
    for name in ["static_sequence_complex", "aa_only"]:
        with np.load(static_root / (name + "_oof.npz")) as z:
            static_predictions[name] = z["prediction"]
    matched_static = {}
    with np.load(HERE / "static_corrected/training_features.npz") as z:
        static_target, static_valid = z["target"], z["valid"]
    for i, record in enumerate(dataset.records):
        identity = {}
        for field in ["structure_path", "label_path"]:
            path = Path(getattr(record, field))
            relative = path.relative_to(PAYLOAD).as_posix()
            identity[field] = sha(path)
            assert_equal(identity[field], input_hashes[relative], "Input hash changed: " + record.sample_id)
        ex = dataset[i]
        target = ex["label_values"][:, :2].numpy()
        valid = (ex["label_mask"][:, :2] & ex["target_mask"][:, None]).numpy()
        reference = static_records[record.sample_id]
        assert_equal(reference["input_sha256"], identity["structure_path"], "Static structure identity differs")
        assert_equal(reference["label_sha256"], identity["label_path"], "Static label identity differs")
        sl = slice(reference["start"], reference["stop"])
        if not np.array_equal(valid, static_valid[sl]):
            raise ValueError("Static/neural heldout masks differ: " + record.sample_id)
        if not np.array_equal(target[valid], static_target[sl].astype(np.float32)[valid]):
            raise ValueError("Static/neural heldout labels differ: " + record.sample_id)
        if not np.isfinite(target[valid]).all():
            raise ValueError("Non-finite valid targets: " + record.sample_id)
        valid_counts += valid.sum(0)
        provenance[record.sample_id] = identity
        matched_static[record.sample_id] = {name: p[sl].copy() for name, p in static_predictions.items()}
        if any(p.shape != target.shape or not np.isfinite(p[valid]).all()
               for p in matched_static[record.sample_id].values()):
            raise ValueError("Invalid matched static predictions")
    return dataset, provenance, matched_static, dict(status="prepared", new_fits=0, expected_members=6,
        seeds=[7, 11, 23], variants=["full", "intrachain"], fold=0, epoch=40,
        records=612, training_records=2451, valid_counts=valid_counts.tolist(),
        test_manifest_sha256=sha(manifest),
        train_manifest_sha256=sha(PAYLOAD / "splits/fold_0_train.jsonl"),
        input_freeze_sha256=FREEZE_SHA, collector_sha256=sha(__file__),
        static_source_sha256={str(p.relative_to(HERE)): sha(p) for p in static_source_paths},
        training_ids=train_ids, heldout_ids=test_ids, device="cpu", threads=2,
        biological_system_unit="Source plus biological_system ID, with records averaged within each such unit.",
        protocol="Fixed fold-0 heldout; separate members, no seed ensemble or score-based selection. Seed sample SD is descriptive, not a family confidence interval.",
        static_match="Fold-0 static and amino-acid OOF member trained on these same 2451 records.")


def load_member(path, variant, seed, member):
    checkpoint = torch.load(path, map_location="cpu", weights_only=False)
    assert_equal(checkpoint["epoch"], 40, "Wrong final epoch")
    config = checkpoint["config"]
    expected = yaml.safe_load((PAYLOAD / f"configs/{member}_fold0.yaml").read_text())
    for key in ["model", "data", "train", "seed", "stage", "init_checkpoint"]:
        assert_equal(config[key], expected[key], "Checkpoint config differs: " + key)
    assert_equal(config["seed"], seed, "Seed mismatch")
    assert_equal(config["model"]["information_ablation"], variant, "Model variant mismatch")
    assert_equal(config["model"]["metrics"], METRICS, "Metric order mismatch")
    model = build_information_ablation(ModelConfig(**config["model"]), variant).cpu().eval()
    model.load_state_dict(checkpoint["model"], strict=True)
    return model


def predict(model, ex):
    with torch.inference_mode():
        result = model.forward_state(collate_states([ex]))
    mean = result["mean"][0].numpy()
    variance = result["log_variance"][0].exp().numpy()
    input_mask = ex["geometry_mask"]
    if not torch.equal(result["mask"][0], input_mask):
        raise ValueError("Model mask changed the prespecified scoring population")
    valid = (ex["label_mask"][:, :2] & ex["target_mask"][:, None]).numpy()
    if mean.shape != ex["label_values"][:, :2].shape or variance.shape != mean.shape:
        raise ValueError("Prediction dimensions differ")
    if not np.isfinite(mean).all() or not np.isfinite(variance).all() or np.any(variance < 0):
        raise ValueError("Non-finite predictions")
    return mean, variance, valid


def statistics(target, prediction):
    error = prediction.astype(np.float64) - target.astype(np.float64)
    result = dict(count=len(error), mae=float(np.mean(abs(error))),
        rmse=float(np.sqrt(np.mean(error ** 2))), bias=float(np.mean(error)),
        spearman=np.nan, pearson=np.nan)
    if len(error) >= 4 and np.ptp(target) > 0 and np.ptp(prediction) > 0:
        result.update(spearman=float(spearmanr(target, prediction).statistic),
                      pearson=float(pearsonr(target, prediction).statistic))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint-root", type=Path, default=HERE / "trained_checkpoints")
    parser.add_argument("--output", type=Path, default=HERE / "fixed_fold_seed_stability")
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()
    if args.smoke_test and not args.audit_only:
        parser.error("--smoke-test requires --audit-only and does not export evaluation scores")
    torch.set_num_threads(2)
    dataset, input_hashes, matched_static, report = check_inputs()
    paths = [(v, s, m, args.checkpoint_root / f"{m}/fold_0/epoch_040.pt") for v, s, m in MEMBERS]
    report["available_members"] = [m for v, s, m, path in paths if path.exists()]
    report["missing_members"] = [m for v, s, m, path in paths if not path.exists()]
    args.output.mkdir(exist_ok=True)
    json_write(args.output / "preflight.json", report)
    if args.audit_only:
        if args.smoke_test:
            smoke = []
            # Exercise short and long inputs plus every source without scoring them.
            smoke_indices = {0, max(range(len(dataset)), key=lambda i: dataset.records[i].metadata["residue_count"])}
            for source in sorted({record.dataset for record in dataset.records}):
                smoke_indices.add(next(i for i, record in enumerate(dataset.records) if record.dataset == source))
            for variant, seed, member, path in paths:
                if path.exists():
                    model = load_member(path, variant, seed, member)
                    for index in sorted(smoke_indices):
                        mean, variance, valid = predict(model, dataset[index])
                        smoke.append(dict(member=member, checkpoint_sha256=sha(path),
                            sample_id=dataset.records[index].sample_id, dataset=dataset.records[index].dataset,
                            finite=True, shape=list(mean.shape), valid_counts=valid.sum(0).tolist()))
                    del model
            json_write(args.output / "smoke_test.json", dict(status="contract_check_only",
                accuracy_evaluation_complete=False, members=smoke, collector_sha256=sha(__file__)))
        print(json.dumps({k: v for k, v in report.items() if not k.endswith("_ids")}))
        return
    if report["missing_members"]:
        raise SystemExit("Missing planned checkpoints; no scores exported: " + ", ".join(report["missing_members"]))
    existing = args.output / "evaluation_audit.json"
    if existing.exists():
        raise SystemExit("Completed output already exists; choose a new output directory for a new audited export.")
    rows, caches, baseline_rows = [], [], []
    started = time.time()
    for i, record in enumerate(dataset.records):
        ex = dataset[i]
        target = ex["label_values"][:, :2].numpy()
        valid = (ex["label_mask"][:, :2] & ex["target_mask"][:, None]).numpy()
        for name, prediction in matched_static[record.sample_id].items():
            for m, metric in enumerate(METRICS):
                if valid[:, m].any():
                    baseline_rows.append(dict(baseline=name, sample_id=record.sample_id,
                        dataset=record.dataset, metric=metric, fold=0,
                        biological_system=record.metadata.get("biological_system", record.source_group),
                        **statistics(target[valid[:, m], m], prediction[valid[:, m], m])))
    for variant, seed, member, path in paths:
        model = load_member(path, variant, seed, member)
        checkpoint_hash = sha(path)
        member_dir = args.output / member
        member_dir.mkdir(exist_ok=True)
        for i, record in enumerate(dataset.records):
            ex = dataset[i]
            target = ex["label_values"][:, :2].numpy()
            expected_valid = (ex["label_mask"][:, :2] & ex["target_mask"][:, None]).numpy()
            structure = parse_pdb(record.structure_path)
            identity = dict(input_hashes[record.sample_id], checkpoint_sha256=checkpoint_hash,
                collector_sha256=report["collector_sha256"], input_freeze_sha256=FREEZE_SHA,
                member=member, seed=seed, epoch=40, sample_id=record.sample_id)
            cache = member_dir / (record.sample_id + ".npz")
            if cache.exists():
                with np.load(cache, allow_pickle=False) as z:
                    assert_equal(json.loads(str(z["provenance"])), identity, "Stale cache provenance")
                    mean, variance, valid = z["mean"], z["model_variance"], z["valid"]
                    if not np.array_equal(z["target"], target, equal_nan=True) or not np.array_equal(valid, expected_valid):
                        raise ValueError("Stale cache target/mask")
                    for name, expected in [("metric_names", np.asarray(METRICS)),
                        ("label_variance", ex["label_variances"][:, :2].numpy()),
                        ("amino_acid", np.asarray(list(structure.sequence))),
                        ("chain", np.asarray(structure.chain_ids)),
                        ("author_residue_id", np.asarray(structure.author_residue_ids))]:
                        saved = z[name]
                        equal = np.array_equal(saved, expected, equal_nan=True) if name == "label_variance" else np.array_equal(saved, expected)
                        if not equal:
                            raise ValueError("Stale cache residue identity or metadata: " + name)
            else:
                mean, variance, valid = predict(model, ex)
                # An incomplete write cannot masquerade as a resumable prediction.
                temporary = cache.with_suffix(".partial.npz")
                np.savez_compressed(temporary, mean=mean, model_variance=variance,
                    metric_names=np.asarray(METRICS),
                    target=target, valid=valid, label_variance=ex["label_variances"][:, :2].numpy(),
                    amino_acid=np.asarray(list(structure.sequence)), chain=np.asarray(structure.chain_ids),
                    author_residue_id=np.asarray(structure.author_residue_ids),
                    provenance=np.asarray(json.dumps(identity, sort_keys=True)))
                temporary.replace(cache)
            if mean.shape != target.shape or variance.shape != target.shape or np.any(variance < 0) or not np.isfinite(mean).all() or not np.isfinite(variance).all():
                raise ValueError("Invalid saved predictions")
            for m, metric in enumerate(METRICS):
                if valid[:, m].any():
                    rows.append(dict(variant=variant, seed=seed, fold=0,
                        sample_id=record.sample_id, dataset=record.dataset,
                        biological_system=record.metadata.get("biological_system", record.source_group),
                        metric=metric, **statistics(target[valid[:, m], m], mean[valid[:, m], m])))
            caches.append(dict(member=member, sample_id=record.sample_id,
                               path=str(cache), sha256=sha(cache)))
            if (i + 1) % 100 == 0:
                print(member, i + 1, len(dataset), flush=True)
        del model
    per = pd.DataFrame(rows)
    per.to_csv(args.output / "per_record_metrics.csv", index=False)
    baseline = pd.DataFrame(baseline_rows)
    baseline.to_csv(args.output / "matched_static_per_record_metrics.csv", index=False)
    columns = ["mae", "rmse", "bias", "spearman", "pearson"]
    keys = ["sample_id", "dataset", "metric", "biological_system", "fold"]
    if per.duplicated(keys + ["seed", "variant"]).any() or baseline.duplicated(keys + ["baseline"]).any():
        raise ValueError("Duplicate paired comparison identities")
    paired = per.merge(baseline, on=keys, suffixes=("_neural", "_static"), validate="many_to_many")
    assert_equal(len(paired), 2 * len(per), "Incomplete matched baseline join")
    if not np.array_equal(paired["count_neural"], paired["count_static"]):
        raise ValueError("Paired static/neural support differs")
    paired["mae_improvement"] = paired["mae_static"] - paired["mae_neural"]
    paired["spearman_improvement"] = paired["spearman_neural"] - paired["spearman_static"]
    paired.to_csv(args.output / "neural_static_paired_metrics.csv", index=False)
    full = per[per.variant == "full"]
    intra = per[per.variant == "intrachain"]
    arm_pairs = full.merge(intra, on=keys + ["seed"], suffixes=("_full", "_intrachain"), validate="one_to_one")
    if not np.array_equal(arm_pairs["count_full"], arm_pairs["count_intrachain"]):
        raise ValueError("Paired full/intrachain support differs")
    arm_pairs["mae_improvement_full"] = arm_pairs["mae_intrachain"] - arm_pairs["mae_full"]
    arm_pairs["spearman_improvement_full"] = arm_pairs["spearman_full"] - arm_pairs["spearman_intrachain"]
    arm_pairs.to_csv(args.output / "same_seed_full_intrachain_paired_metrics.csv", index=False)
    for comparison, frame, group_keys, count_column, correlation_column in [
        ("neural_static", paired, ["variant", "seed", "dataset", "metric", "baseline"], "count_neural", "spearman_improvement"),
        ("same_seed_full_intrachain", arm_pairs, ["seed", "dataset", "metric"], "count_full", "spearman_improvement_full"),
    ]:
        pair_support = []
        frame = frame.assign(source_dataset=frame["dataset"])
        for population in [frame, frame.assign(dataset="all")]:
            for identity, part in population.groupby(group_keys, sort=True):
                system_correlations = part.groupby(["source_dataset", "biological_system"])[correlation_column].mean()
                item = dict(zip(group_keys, identity))
                item.update(records=len(part), valid_residues=int(part[count_column].sum()),
                    biological_systems=len(system_correlations),
                    spearman_both_computable_records=int(part[correlation_column].notna().sum()),
                    spearman_both_computable_systems=int(system_correlations.notna().sum()))
                pair_support.append(item)
        pd.DataFrame(pair_support).to_csv(args.output / (comparison + "_paired_metric_support.csv"), index=False)
    supports = []
    with_source = per.assign(source_dataset=per["dataset"])
    for grouping in [with_source.groupby(["variant", "seed", "dataset", "metric"], sort=True),
                     with_source.assign(dataset="all").groupby(["variant", "seed", "dataset", "metric"], sort=True)]:
        for (variant, seed, source, metric), part in grouping:
            system_correlations = part.groupby(["source_dataset", "biological_system"])[["spearman", "pearson"]].mean()
            supports.append(dict(variant=variant, seed=seed, dataset=source, metric=metric,
                records=len(part), valid_residues=int(part["count"].sum()),
                biological_systems=len(part[["source_dataset", "biological_system"]].drop_duplicates()),
                spearman_computable=int(part["spearman"].notna().sum()),
                pearson_computable=int(part["pearson"].notna().sum()),
                spearman_computable_systems=int(system_correlations["spearman"].notna().sum()),
                pearson_computable_systems=int(system_correlations["pearson"].notna().sum())))
    pd.DataFrame(supports).to_csv(args.output / "metric_support.csv", index=False)
    for estimator in ["record_macro", "biological_system_macro"]:
        data = per
        if estimator == "biological_system_macro":
            data = per.groupby(["variant", "seed", "dataset", "biological_system", "metric"])[columns].mean().reset_index()
        by_source = data.groupby(["variant", "seed", "dataset", "metric"])[columns].mean().reset_index()
        overall = data.groupby(["variant", "seed", "metric"])[columns].mean().reset_index().assign(dataset="all")
        macro = pd.concat([by_source, overall], ignore_index=True)
        macro.to_csv(args.output / (estimator + ".csv"), index=False)
        summary = macro.groupby(["variant", "dataset", "metric"])[columns].agg(["mean", "std", "min", "max"])
        summary.columns = [f"{metric}_{statistic}" for metric, statistic in summary.columns]
        summary["seed_count"] = macro.groupby(["variant", "dataset", "metric"])["seed"].nunique()
        for metric in columns:
            summary[metric + "_finite_seed_count"] = macro.groupby(["variant", "dataset", "metric"])[metric].count()
        summary.to_csv(args.output / (estimator + "_three_seed_descriptive_summary.csv"))
    json_write(args.output / "prediction_manifest.json", caches)
    report.update(status="complete", elapsed_seconds=time.time() - started,
        predictions=len(caches), collected_trained_members=6, checkpoints={m: sha(path) for v, s, m, path in paths},
        python=sys.version, torch=torch.__version__, numpy=np.__version__, platform=platform.platform(),
        inference_cost_scope="CPU six-member export wall clock including model load and file IO, excluding input preflight; not a matched algorithm cost benchmark.",
        statistic_scope="Sample SD (ddof=1) over three training seeds at fixed fold 0. No independent-family CI or p value.")
    report["correlation_population"] = "record_macro averages finite record correlations; biological_system_macro first averages finite record correlations within source dataset plus biological system, then equally averages computable systems. Constant targets or predictions yield NaN, with no imputation. Record and computable system counts are reported per seed/source. Paired gains retain NaN when either member is undefined; paired support tables count records/systems with both correlations computable."
    json_write(existing, report)
    print("Completed six fixed-fold seed members; descriptive seed variation only.", flush=True)


if __name__ == "__main__":
    main()
