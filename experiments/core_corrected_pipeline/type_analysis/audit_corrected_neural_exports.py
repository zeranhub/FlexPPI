"""Independently verify saved corrected predictions without fitting or inference.

--preflight checks real frozen inputs and membership only. A complete audit
requires both successful queue states, all 16 actual checkpoints, and every
saved OOF/external prediction. It never fills missing predictions or changes
the frozen exporters. Reports go to a new directory, preventing stale success.
"""
from pathlib import Path, PurePosixPath
from datetime import datetime, timezone
from functools import lru_cache
import argparse
import hashlib
import json
import os
import sys
import uuid

os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import numpy as np
import pandas as pd
import yaml
from scipy.stats import pearsonr, spearmanr

BASE = Path(__file__).resolve().parents[1]
PIPE = BASE
PAYLOAD = PIPE / "training_payload"
STRICT = PIPE / "strict_external_payload"
EXTERNAL = PIPE / "external_payload/external"
REMOTE = "/root/flexppi_revision_20261002/corrected_training_v1"
sys.path.insert(0, str(PAYLOAD / "src"))
from flexppi.constants import METRIC_NAMES
from flexppi.data.structures import parse_pdb, target_mask
from flexppi.data.labels import load_label_file

METRICS = ["rmsf", "chi_entropy"]
SCORES = ["mae", "rmse", "bias", "spearman", "pearson"]
HASHES = {
    "primary": "5b994bb59684e3a33c97f4a1801ee9a606e5dd856db9cfe23bf18515697707d3",
    "strict": "16be2a49c6cd2b295a8f150de5e0354d3f6e20028de05f530ae8908c26dbacea",
    "external": "4dc8c3d560ac7e11eb2bef1be5bab431c58be32281f2712e4c8019715fe5781b",
    "oof_code": "59be9322393daae296610d5e684cc6d62ea9d815b2c88a3f59db260cd4698fa7",
    "external_code": "d809d798b2b7339d3aa429208da6f2a6f34748c370d9dfa4e19ee43fd467ded5",
    "strict_code": "87623f0466024e92673a2f0b2d622c6763b6b6079d8196d6888c2985367a1c64",
}
PRIMARY_MEMBERS = [(v, k, 7, v) for v in ["full", "intrachain"] for k in range(5)]
PRIMARY_MEMBERS += [(v, 0, s, f"{v}_seed{s}") for v in ["full", "intrachain"] for s in [11, 23]]
STRICT_MEMBERS = [(v, 0, 7, f"strict_{v}") for v in ["full", "intrachain"]]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_lines(path):
    return [json.loads(line) for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]


def portable(root, relative):
    parts = PurePosixPath(relative.replace("\\", "/")).parts
    require(parts and not any(x in ["..", "/"] for x in parts), "Unsafe relative path")
    result = root.joinpath(*parts)
    require(result.resolve().is_relative_to(root.resolve()), "Path escaped intended root")
    return result


def check_hash_table(root, table, evidence):
    for name, digest in table.items():
        path = portable(root, name)
        actual = sha(path)
        require(actual == digest, "Frozen file changed: " + str(path))
        evidence[str(path)] = actual


def dictionary(rows, key="sample_id"):
    result = {r[key]: r for r in rows}
    require(len(result) == len(rows), "Duplicate " + key)
    return result


def local_input_path(row, field):
    remote = str(row[field]).replace("\\", "/")
    require(remote.startswith(REMOTE + "/"), "Input outside frozen execution root")
    relative = remote[len(REMOTE) + 1:]
    if relative.startswith("external/"):
        return portable(EXTERNAL, relative[len("external/"):])
    return portable(PAYLOAD, relative)


@lru_cache(maxsize=96)
def reference(input_path, label_path, target_chains):
    structure = parse_pdb(input_path)
    values, variances, mask = load_label_file(label_path)
    require(values.shape == (structure.num_residues, 3), "Label/structure dimensions differ")
    require(list(METRIC_NAMES) == ["rmsf", "chi_entropy", "sasa"], "Frozen label names changed")
    require(not mask[:, 2].any(), "Unexpected SASA labels; this audit is for two supervised targets")
    require(set(target_chains).issubset(set(structure.chain_ids)), "Target chain absent")
    valid = mask & structure.frame_mask[:, None] & target_mask(structure, target_chains)[:, None]
    require(np.isfinite(values[valid]).all(), "Non-finite valid target")
    require(np.isfinite(variances[valid]).all() and np.all(variances[valid] >= 0), "Invalid label variance")
    identity = list(zip(structure.chain_ids, structure.author_residue_ids))
    require(len(identity) == len(set(identity)), "Duplicate chain/author residue identity")
    return dict(target=values, label_variance=variances, valid=valid,
                amino_acid=np.asarray(list(structure.sequence)), chain=np.asarray(structure.chain_ids),
                author_residue_id=np.asarray(structure.author_residue_ids), aatype=structure.aatype)


def row_reference(row):
    return reference(str(local_input_path(row, "structure_path")), str(local_input_path(row, "label_path")),
                     tuple(row["target_chains"]))


def input_preflight():
    evidence = {}
    freeze_paths = {"primary": PAYLOAD / "training_freeze.json", "strict": STRICT / "strict_external_freeze.json",
                    "external": EXTERNAL / "external_payload_freeze.json"}
    freezes = {}
    for name, path in freeze_paths.items():
        require(sha(path) == HASHES[name], "Freeze changed: " + name)
        evidence[str(path)] = sha(path)
        freezes[name] = read_json(path)
    check_hash_table(PAYLOAD, freezes["primary"]["execution_hashes"], evidence)
    check_hash_table(PAYLOAD, freezes["primary"]["input_hashes"], evidence)
    check_hash_table(STRICT, freezes["strict"]["execution_hashes"], evidence)
    check_hash_table(EXTERNAL, freezes["external"]["data_sha256"], evidence)
    for key, path in [("oof_code", PAYLOAD / "evaluate_corrected_models.py"),
                      ("external_code", PIPE / "external_payload/evaluate_external_models.py")]:
        require(sha(path) == HASHES[key], "Export implementation changed: " + key)
        evidence[str(path)] = sha(path)
    protocol_path = EXTERNAL / "prespecified_protocol.json"
    require(sha(protocol_path) == freezes["external"]["external_protocol_sha256"], "External protocol changed")
    evidence[str(protocol_path)] = sha(protocol_path)
    training_rows = read_lines(PAYLOAD / "splits/all_ready.jsonl")
    records = dictionary(training_rows)
    require(len(records) == 3063, "Unexpected corrected corpus size")
    folds = {}
    heldout = []
    for fold in range(5):
        test = read_lines(PAYLOAD / f"splits/fold_{fold}_test.jsonl")
        train = read_lines(PAYLOAD / f"splits/fold_{fold}_train.jsonl")
        test_dict, train_dict = dictionary(test), dictionary(train)
        require(not set(test_dict) & set(train_dict), "Train/heldout overlap")
        require(set(test_dict) | set(train_dict) == set(records), "Incomplete fold membership")
        for key, row in test_dict.items():
            require(row == records[key], "Heldout record differs from frozen corpus")
            require(int(row["metadata"]["fold"]) == fold, "Record fold differs")
            require(key not in folds, "Record held out more than once")
            folds[key] = fold
        for key, row in train_dict.items():
            require(row == records[key], "Training record differs from frozen corpus")
        heldout.append(len(test))
    require(heldout == [612, 612, 611, 614, 614] and set(folds) == set(records), "Wrong five-fold support")
    external_rows = read_lines(EXTERNAL / "ready_test.jsonl")
    external_records = dictionary(external_rows)
    require(len(external_records) == 32 and not set(records) & set(external_records), "Source-out membership mismatch")
    protocol, extraction = read_json(protocol_path), read_json(EXTERNAL / "extraction_summary.json")
    require(protocol["accepted_count"] == extraction["prespecified"] == extraction["completed"] == 32,
            "Prespecified cohort/extraction differs")
    require(not extraction["failed"], "Unresolved external extraction failure")
    strict_records = dictionary(read_lines(STRICT / "splits/strict_external_train.jsonl"))
    require(len(strict_records) == 2431 and set(strict_records).issubset(records), "Strict training membership differs")
    for key, row in strict_records.items():
        require(row == records[key], "Strict input record differs")
    overlap_path = PIPE / "sequence_audit/record_overlap_audit.csv"
    overlap = pd.read_csv(overlap_path)
    retained = set(overlap.loc[(overlap.role == "training") & ~overlap.purged_training_member, "sample_id"])
    require(retained == set(strict_records), "Strict membership differs from all-chain removal")
    external_overlap = overlap[overlap.role == "external"]
    require(set(external_overlap.sample_id) == set(external_records) and external_overlap.component_id.nunique() == 1,
            "External family dependence changed")
    evidence[str(overlap_path)] = sha(overlap_path)
    counts = {}
    for scope, cohort, expected in [("oof", records, [961321, 824988]),
                                    ("external", external_records, [18748, 16265])]:
        total = np.zeros(2, dtype=np.int64)
        for row in cohort.values():
            ref = row_reference(row)
            n = len(ref["target"])
            require(n == int(row["metadata"]["residue_count"]) and n <= 1024, "Input length metadata differs")
            total += ref["valid"][:, :2].sum(0)
        require(total.tolist() == expected, "Valid population differs: " + scope)
        counts[scope] = total.tolist()
    report = dict(status="prepared_only; final predictions pending", input_sha256=evidence,
                  actual_records={"oof_per_arm": 3063, "external": 32, "strict_training": 2431},
                  heldout_folds=heldout, actual_valid_counts=counts, external_family_components=1,
                  expected_caches={"oof": 6126, "external_individual": 448, "external_ensemble": 64,
                                   "strict_individual": 64, "strict_single_compatibility": 64},
                  supervised_metrics=METRICS, preserved_label_columns=list(METRIC_NAMES),
                  new_fits=0, new_forward_calls=0, final_prediction_audit_complete=False)
    return records, external_records, folds, report


def check_cache(z, ref, oof=False):
    required = {"mean", "model_variance", "target", "valid", "label_variance", "amino_acid", "chain", "author_residue_id"}
    require(required.issubset(z), "Missing prediction fields")
    n = len(ref["target"])
    require(z["mean"].shape == z["model_variance"].shape == (n, 2), "Two-output prediction dimensions differ")
    for name in ["target", "valid", "label_variance", "amino_acid", "chain", "author_residue_id"]:
        require(np.array_equal(z[name], ref[name], equal_nan=True) if z[name].dtype.kind not in "USb"
                else np.array_equal(z[name], ref[name]), "Prediction/source differs: " + name)
    require(z["valid"].dtype.kind == "b", "Scoring mask is not boolean")
    require(np.isfinite(z["mean"]).all(), "Non-finite predictions")
    require(np.all(z["mean"] >= 0), "Negative absolute-state means cannot come from frozen softplus head")
    require(np.isfinite(z["model_variance"]).all() and np.all(z["model_variance"] >= 0), "Invalid model variance")
    if oof:
        require("aatype" in z and np.array_equal(z["aatype"], ref["aatype"]), "OOF amino-acid indices differ")


def check_ensemble(z, members, strict):
    means = np.asarray([m["mean"] for m in members], dtype=np.float64)
    variances = np.asarray([m["model_variance"] for m in members], dtype=np.float64)
    require(len(members) == (1 if strict else 5), "Wrong ensemble membership")
    expected = dict(member_means=means, mean=means.mean(0), within_model_variance=variances.mean(0),
                    between_fold_variance=means.var(0, ddof=0))
    expected["model_variance"] = expected["within_model_variance"] + expected["between_fold_variance"]
    for name, values in expected.items():
        require(name in z and z[name].shape == values.shape and np.allclose(z[name], values, rtol=1e-12, atol=1e-12),
                "Ensemble membership/mean/variance differs: " + name)
    if strict:
        require(np.count_nonzero(z["between_fold_variance"]) == 0, "Single strict member has between-fold variance")


def statistics(y, prediction):
    require(np.isfinite(y).all() and np.isfinite(prediction).all(), "Invalid scoring support")
    difference = prediction.astype(np.float64) - y.astype(np.float64)
    result = dict(count=len(y), mae=float(np.abs(difference).mean()),
                  rmse=float(np.sqrt(np.square(difference).mean())), bias=float(difference.mean()),
                  spearman=np.nan, pearson=np.nan)
    if len(y) >= 4 and np.ptp(y) > 0 and np.ptp(prediction) > 0:
        result.update(spearman=float(spearmanr(y, prediction).statistic),
                      pearson=float(pearsonr(y.astype(np.float64), prediction.astype(np.float64)).statistic))
    return result


def score_rows(z, identifiers):
    rows = []
    for m, metric in enumerate(METRICS):
        valid = z["valid"][:, m]
        if valid.any():
            rows.append(dict(**identifiers, metric=metric, **statistics(z["target"][valid, m], z["mean"][valid, m])))
    return rows


def compare_table(expected, saved_path, keys, values):
    saved = pd.read_csv(saved_path, keep_default_na=False, na_values=[""])
    require(not expected.duplicated(keys).any() and not saved.duplicated(keys).any(), "Duplicate metric row")
    left, right = expected.set_index(keys).sort_index(), saved.set_index(keys).sort_index()
    require(left.index.equals(right.index), "Metric row support differs: " + str(saved_path))
    for name in values:
        a, b = left[name].to_numpy(dtype=np.float64), right[name].to_numpy(dtype=np.float64)
        require(np.array_equal(np.isnan(a), np.isnan(b)), "Correlation validity differs: " + name)
        tol = 1e-10 if name in ["spearman", "count"] else 3e-6
        require(np.allclose(a, b, rtol=tol, atol=tol, equal_nan=True), "Recomputed metric differs: " + name)


def checkpoint_preflight(weight_root, primary_state, strict_state):
    require(primary_state["status"] == strict_state["status"] == "complete", "Both export queues must be complete")
    require(primary_state["evaluation_exit_code"] == primary_state["external_exit_code"] == 0,
            "Primary export did not exit successfully")
    require(strict_state["evaluation_exit_code"] == 0, "Strict export did not exit successfully")
    expected_jobs = {(member, fold) for _, fold, _, member in PRIMARY_MEMBERS}
    jobs = {(j["variant"], j["fold"]): j for j in primary_state["jobs"]}
    require(len(jobs) == len(primary_state["jobs"]) == 14 and set(jobs) == expected_jobs, "Primary queue membership differs")
    for job in jobs.values():
        require(job["status"] in ["complete", "already_complete"], "Incomplete primary fit")
        if job["status"] == "complete":
            require(job["exit_code"] == 0, "Failed primary fit")
    strict_jobs = {j["variant"]: j for j in strict_state["jobs"]}
    require(len(strict_jobs) == len(strict_state["jobs"]) == 2 and set(strict_jobs) == {"full", "intrachain"}, "Strict jobs differ")
    require(all(j["status"] == "complete" for j in strict_jobs.values()), "Incomplete strict fit")
    collected = {}
    for name, field, expected_count, freeze_field, freeze_key in [
        ("completed_members_audit.json", "completed_primary_members", 14, "training_freeze_sha256", "primary"),
        ("completed_strict_members_audit.json", "completed_strict_members", 2, "strict_freeze_sha256", "strict")]:
        report = read_json(weight_root / name)
        require(report[field] == expected_count and report[freeze_field] == HASHES[freeze_key], "Actual checkpoint collection is incomplete")
        require(len(report["checkpoints"]) == expected_count, "Actual collected checkpoint count differs")
        for row in report["checkpoints"]:
            key = (row["member"], row["fold"])
            require(key not in collected and row["epoch"] == 40, "Duplicate or wrong-epoch collected checkpoint")
            require(row["source"] == REMOTE + f"/runs/{row['member']}/fold_{row['fold']}/epoch_040.pt", "Collected checkpoint source differs")
            collected[key] = row["checkpoint_sha256"]
    require(set(collected) == {(m, f) for v, f, s, m in PRIMARY_MEMBERS + STRICT_MEMBERS}, "Collected member identities differ")
    import torch
    torch.set_num_threads(2)
    evidence = {}
    for variant, fold, seed, member in PRIMARY_MEMBERS + STRICT_MEMBERS:
        path = weight_root / f"{member}/fold_{fold}/epoch_040.pt"
        config_path = STRICT / f"configs/{member}.yaml" if member.startswith("strict_") else PAYLOAD / f"configs/{member}_fold{fold}.yaml"
        expected = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
        require(checkpoint["epoch"] == 40 and checkpoint["config"]["seed"] == seed, "Wrong final epoch/seed")
        for key in ["model", "data", "train", "seed", "stage", "init_checkpoint"]:
            require(checkpoint["config"][key] == expected[key], "Checkpoint configuration differs: " + key)
        require(checkpoint["config"]["model"]["information_ablation"] == variant, "Wrong model arm")
        require(checkpoint["config"]["model"]["metrics"] == METRICS, "Wrong output names")
        require(all(torch.isfinite(t).all().item() for t in checkpoint["model"].values() if t.is_floating_point()), "Non-finite model parameters")
        evidence[(member, fold)] = sha(path)
        require(evidence[(member, fold)] == collected[(member, fold)], "Checkpoint differs from actual remote collection")
        del checkpoint
    return evidence


def load_cache(path, hashes):
    digest = sha(path)
    with np.load(path, allow_pickle=False) as archive:
        result = {k: archive[k] for k in archive.files}
    require(sha(path) == digest, "Prediction changed during audit")
    hashes[str(path)] = digest
    return result


def metadata_snapshot(directory, hashes):
    for path in directory.iterdir():
        if path.is_file() and path.suffix in [".csv", ".gz", ".json", ".jsonl"]:
            hashes[str(path)] = sha(path)


def manifest_index(directory, expected_keys, oof):
    rows = read_lines(directory / "prediction_manifest.jsonl")
    index = {}
    for row in rows:
        key = (row["variant"], int(row["fold"]), 7 if oof else int(row["seed"]), row["sample_id"])
        require(key not in index, "Duplicate prediction manifest entry")
        index[key] = row
    require(set(index) == expected_keys, "Incomplete or unexpected prediction manifest")
    return index


def check_manifest(row, record, checkpoint_hash, relative, oof=False):
    required = {"sample_id", "variant", "fold", "checkpoint_sha256", "input_sha256", "label_sha256", "cache_path"}
    required |= {"record", "dataset"} if oof else {"seed", "forward_seconds", "residue_count"}
    require(required.issubset(row), "Incomplete required prediction manifest fields")
    require(row["checkpoint_sha256"] == checkpoint_hash, "Prediction/checkpoint hash differs")
    require(row["input_sha256"] == sha(local_input_path(record, "structure_path")), "Prediction/input hash differs")
    require(row["label_sha256"] == sha(local_input_path(record, "label_path")), "Prediction/label hash differs")
    require(row["cache_path"].replace("\\", "/") == REMOTE + "/" + relative, "Unexpected prediction cache path")
    if oof:
        require(row["record"] == record and row["dataset"] == record["dataset"], "OOF manifest record differs")
    else:
        require(np.isfinite(row["forward_seconds"]) and row["forward_seconds"] > 0, "Invalid forward timing")
        require(row["residue_count"] == record["metadata"]["residue_count"], "Timing record length differs")


def per_residue_audit(path, cache_paths, hashes, oof):
    coverage = {key: np.zeros(len(row_reference(record)["target"]), dtype=bool)
                for key, (cache, record) in cache_paths.items()}
    total = 0
    dtype = {"sample_id": str, "variant": str, "amino_acid": str, "chain": str, "author_residue_id": str}
    for chunk in pd.read_csv(path, compression="gzip", chunksize=100000, dtype=dtype, keep_default_na=False, na_values=[""]):
        for key, part in chunk.groupby(["variant", "sample_id"], sort=False):
            require(key in cache_paths, "Unexpected residue table record")
            cache_path, record = cache_paths[key]
            require(sha(cache_path) == hashes[str(cache_path)], "Cache changed before residue table check")
            with np.load(cache_path, allow_pickle=False) as z:
                indices = part.residue_index.to_numpy(dtype=np.int64)
                require(np.array_equal(indices, part.residue_index.to_numpy()), "Noninteger residue index")
                require(len(set(indices)) == len(indices) and np.all((indices >= 0) & (indices < len(z["mean"]))), "Invalid residue indices")
                require(not coverage[key][indices].any(), "Duplicate residue table position")
                for name in ["amino_acid", "chain", "author_residue_id"]:
                    require(np.array_equal(part[name].to_numpy(), z[name][indices]), "Residue table identity differs: " + name)
                if oof:
                    require((part.fold == record["metadata"]["fold"]).all() and (part.dataset == record["dataset"]).all(), "Residue table fold/source differs")
                for m, metric in enumerate(METRICS):
                    valid = part[f"valid_{metric}"]
                    require(valid.dtype.kind == "b" and np.array_equal(valid.to_numpy(), z["valid"][indices, m]), "Residue table masks differ")
                    comparisons = {f"target_{metric}": z["target"][indices, m], f"pred_{metric}": z["mean"][indices, m],
                                   f"sd_{metric}": np.sqrt(z["model_variance"][indices, m])}
                    if not oof:
                        comparisons[f"between_fold_sd_{metric}"] = np.sqrt(z["between_fold_variance"][indices, m])
                    for name, values in comparisons.items():
                        require(np.allclose(part[name].to_numpy(dtype=np.float64), values, rtol=1e-6, atol=1e-7, equal_nan=True), "Residue table numbers differ: " + name)
                coverage[key][indices] = True
                total += len(indices)
    require(all(x.all() for x in coverage.values()), "Missing residue table positions")
    return total


def audit_oof(root, records, folds, checkpoint_hashes, hashes):
    directory = root / "evaluation"
    metadata_snapshot(directory, hashes)
    audit = read_json(directory / "evaluation_audit.json")
    require(audit["status"] == "complete" and audit["input_freeze_sha256"] == HASHES["primary"]
            and audit["export_code_sha256"] == HASHES["oof_code"], "OOF provenance audit differs")
    expected_keys = {(v, folds[s], 7, s) for v in ["full", "intrachain"] for s in records}
    manifests = manifest_index(directory, expected_keys, True)
    results, cache_paths = [], {}
    for variant, fold, seed, sample in sorted(expected_keys):
        record = records[sample]
        relative = f"evaluation/{variant}/{sample}.npz"
        check_manifest(manifests[(variant, fold, seed, sample)], record, checkpoint_hashes[(variant, fold)], relative, oof=True)
        cache_path = portable(root, relative)
        z = load_cache(cache_path, hashes)
        check_cache(z, row_reference(record), oof=True)
        results += score_rows(z, dict(sample_id=sample, dataset=record["dataset"], variant=variant, fold=fold))
        cache_paths[(variant, sample)] = (cache_path, record)
    require(len(list(directory.rglob("*.npz"))) == 6126, "Unexpected OOF cache count")
    per = pd.DataFrame(results)
    compare_table(per, directory / "per_record_metrics.csv", ["variant", "fold", "dataset", "sample_id", "metric"], ["count"] + SCORES)
    grouped = per.groupby(["variant", "fold", "dataset", "metric"])[SCORES].mean().reset_index()
    all_sources = per.groupby(["variant", "fold", "metric"])[SCORES].mean().reset_index().assign(dataset="all")
    grouped = pd.concat([grouped, all_sources], ignore_index=True)
    compare_table(grouped, directory / "source_fold_macro_metrics.csv", ["variant", "fold", "dataset", "metric"], SCORES)
    summary = grouped.groupby(["variant", "dataset", "metric"])[SCORES].agg(["mean", "std"])
    summary.columns = [f"{name}_{scope}" for name, scope in summary.columns]
    compare_table(summary.reset_index(), directory / "five_fold_mean_sample_sd.csv", ["variant", "dataset", "metric"], list(summary.columns))
    require(audit["records_per_variant"] == 3063, "OOF declared record count differs")
    require(audit["effective_residues_per_variant"] == {"valid_rmsf": {"full": 961321, "intrachain": 961321},
            "valid_chi_entropy": {"full": 824988, "intrachain": 824988}}, "OOF declared label count differs")
    residue_rows = per_residue_audit(directory / "oof_per_residue.csv.gz", cache_paths, hashes, True)
    return dict(status="verified", manifest_entries=len(manifests), cache_count=len(cache_paths),
                recomputed_record_metric_rows=len(per), checked_residue_rows=residue_rows, fold_sd="sample SD; ddof=1")


def audit_external(root, records, checkpoint_hashes, hashes, strict=False):
    name = "strict_external_evaluation" if strict else "external_evaluation"
    directory = root / name
    metadata_snapshot(directory, hashes)
    audit = read_json(directory / "evaluation_audit.json")
    require(audit["status"] == "complete" and audit["epoch"] == 40 and audit["scored_records"] == audit["prespecified_records"] == 32,
            "External declared completion differs")
    require(audit["training_freeze_sha256"] == HASHES["strict" if strict else "primary"]
            and audit["external_payload_sha256"] == HASHES["external"]
            and audit["export_code_sha256"] == HASHES["strict_code" if strict else "external_code"], "External provenance differs")
    require(audit["prespecified_protocol_sha256"] == sha(EXTERNAL / "prespecified_protocol.json") and not audit["extraction_failures"], "External protocol/failures differ")
    members = STRICT_MEMBERS if strict else PRIMARY_MEMBERS
    require(audit["checkpoints"] == len(members), "Declared external checkpoint count differs")
    expected_keys = {(v, f, s, sample) for v, f, s, member in members for sample in records}
    manifests = manifest_index(directory, expected_keys, False)
    results, cache_paths = [], {}
    for variant in ["full", "intrachain"]:
        for sample, record in records.items():
            ensemble_members = []
            for v, fold, seed, member in members:
                if v != variant:
                    continue
                relative = f"{name}/{member}/fold_{fold}/{sample}.npz"
                check_manifest(manifests[(v, fold, seed, sample)], record, checkpoint_hashes[(member, fold)], relative)
                z = load_cache(portable(root, relative), hashes)
                check_cache(z, row_reference(record))
                results += score_rows(z, dict(sample_id=sample, variant=v, fold=fold, seed=seed, role="individual_model"))
                if seed == 7:
                    ensemble_members.append(z)
            ensemble_path = directory / variant / "ensemble" / (sample + ".npz")
            z = load_cache(ensemble_path, hashes)
            check_cache(z, row_reference(record))
            check_ensemble(z, ensemble_members, strict)
            role = "single_strict_model" if strict else "five_fold_ensemble"
            results += score_rows(z, dict(sample_id=sample, variant=variant, fold=-1, seed=7, role=role))
            cache_paths[(variant, sample)] = (ensemble_path, record)
    require(len(list(directory.rglob("*.npz"))) == (128 if strict else 512), "Unexpected external cache count")
    per = pd.DataFrame(results)
    compare_table(per, directory / "per_record_metrics.csv", ["variant", "role", "fold", "seed", "sample_id", "metric"], ["count"] + SCORES)
    macro = per.groupby(["variant", "role", "fold", "seed", "metric"])[SCORES].mean().reset_index()
    compare_table(macro, directory / "macro_metrics.csv", ["variant", "role", "fold", "seed", "metric"], SCORES)
    residue_rows = per_residue_audit(directory / "per_residue.csv.gz", cache_paths, hashes, False)
    return dict(status="verified", manifest_entries=len(manifests), individual_cache_count=len(manifests),
                aggregate_cache_count=64, recomputed_record_metric_rows=len(per), checked_residue_rows=residue_rows,
                aggregate_members=1 if strict else 5, external_family_components=1, family_ci_and_p_values=False,
                strict_aggregate_is_same_prediction_as_individual=strict)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--export-root", type=Path, default=PIPE / "neural_exports")
    parser.add_argument("--weight-root", type=Path, default=PIPE / "trained_checkpoints")
    parser.add_argument("--report-dir", type=Path)
    args = parser.parse_args()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    output = args.report_dir or BASE / "qa/corrected_neural_export_audits" / stamp
    require(not output.exists(), "Refused to overwrite previous audit directory")
    output.mkdir(parents=True)
    report = dict(checked_utc=datetime.now(timezone.utc).isoformat(), auditor_sha256=sha(__file__),
                  new_fits=0, new_forward_calls=0, status="started; not complete", final_prediction_audit_complete=False)
    try:
        records, external, folds, preflight = input_preflight()
        report.update(preflight)
        if not args.preflight:
            primary_path, strict_path = args.export_root / "queue_state.json", args.export_root / "strict_queue_state.json"
            checkpoint_hashes = checkpoint_preflight(args.weight_root, read_json(primary_path), read_json(strict_path))
            report["queue_state_sha256"] = {str(p): sha(p) for p in [primary_path, strict_path]}
            report["checkpoint_sha256"] = {f"{m}/fold_{f}": s for (m, f), s in checkpoint_hashes.items()}
            hashes = {}
            report["scopes"] = dict(oof=audit_oof(args.export_root, records, folds, checkpoint_hashes, hashes),
                external=audit_external(args.export_root, external, checkpoint_hashes, hashes),
                strict=audit_external(args.export_root, external, checkpoint_hashes, hashes, strict=True))
            for name in ["completed_members_audit.json", "completed_strict_members_audit.json"]:
                path = args.weight_root / name
                report["queue_state_sha256"][str(path)] = sha(path)
            for path, digest in {**report["input_sha256"], **report["queue_state_sha256"], **hashes}.items():
                require(sha(path) == digest, "Consumed source/output changed before audit completion: " + path)
            for (member, fold), digest in checkpoint_hashes.items():
                require(sha(args.weight_root / f"{member}/fold_{fold}/epoch_040.pt") == digest, "Checkpoint changed during final audit")
            report.update(status="complete; all actual saved predictions independently verified",
                          final_prediction_audit_complete=True, output_sha256=hashes,
                          metric_tolerance="Per-record Spearman/count 1e-10; MAE/RMSE/bias/Pearson and all summary columns 3e-6 absolute/relative; ensemble 1e-12; residue CSV 1e-6 relative/1e-7 absolute.")
        print(report["status"], flush=True)
    except Exception as error:
        report.update(status="failed or not ready; no final prediction success", final_prediction_audit_complete=False,
                      error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "audit.json").write_text(json.dumps(report, indent=2, ensure_ascii=True), encoding="utf-8")
        print("Audit report:", output / "audit.json", flush=True)


if __name__ == "__main__":
    main()
