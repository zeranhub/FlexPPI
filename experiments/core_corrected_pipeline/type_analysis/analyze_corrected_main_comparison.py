"""Read-only, common-mask main comparison of the frozen corrected cohort.

Before a CURRENT complete neural export audit, use --preflight only. This
program never fits, loads a model/checkpoint, performs a forward pass, creates
figures, bootstraps, or scores the single-family external cohort. Final CSVs
are written only after all audit bindings and actual OOF caches pass.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
import uuid

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

BASE = Path(__file__).resolve().parents[1]
PIPE = BASE
PAYLOAD = PIPE / "training_payload"
STATIC = PIPE / "static_corrected"
BF = PIPE / "backflip2_corrected_scoring"
BF_INPUT = PIPE / "backflip2_corrected_inputs"
BF_PRED = PIPE / "backflip2_corrected_completed/predictions"
SEQ = PIPE / "sequence_audit"
REMOTE = "/root/flexppi_revision_20261002/corrected_training_v1"
QA = BASE / "qa"
PRIMARY_SHA = "5b994bb59684e3a33c97f4a1801ee9a606e5dd856db9cfe23bf18515697707d3"
STRICT_SHA = "16be2a49c6cd2b295a8f150de5e0354d3f6e20028de05f530ae8908c26dbacea"
METRICS = ("rmsf", "chi_entropy")
SOURCES = {"atlas": 1929, "dynarepo": 890, "itrace": 244}
SCORES = ("mae", "rmse", "bias", "spearman", "pearson", "abs_bias")
PRIMARY_MEMBERS = [(v, f, v) for v in ("full", "intrachain") for f in range(5)]
PRIMARY_MEMBERS += [(v, 0, v + "_seed" + str(s)) for v in ("full", "intrachain") for s in (11, 23)]
MEMBERS = PRIMARY_MEMBERS + [(v, 0, "strict_" + v) for v in ("full", "intrachain")]
COMPLETE_STATUS = "complete; all actual saved predictions independently verified"


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
    return [json.loads(s) for s in Path(path).read_text(encoding="utf-8").splitlines() if s.strip()]


def csv(path):
    return pd.read_csv(path, keep_default_na=False, na_values=[""])


def portable(root, name):
    parts = PurePosixPath(str(name).replace("\\", "/")).parts
    require(parts and not any(p in ("..", "/") for p in parts), "Unsafe relative path")
    path = root.joinpath(*parts).resolve()
    require(path.is_relative_to(root.resolve()), "Path escaped root")
    return path


def local_input(row, field):
    remote = row[field].replace("\\", "/")
    require(remote.startswith(REMOTE + "/dataset/"), "OOF input outside frozen training dataset")
    return portable(PAYLOAD, remote[len(REMOTE) + 1:])


def bind(path, evidence, expected=None):
    path = Path(path).resolve()
    actual = sha(path)
    if expected is not None:
        require(actual == expected, "SHA mismatch: " + str(path))
    key = str(path)
    require(key not in evidence or evidence[key] == actual, "File changed while reading: " + key)
    evidence[key] = actual
    return actual


def check_snapshot(evidence):
    for path, digest in evidence.items():
        require(sha(path) == digest, "Consumed file changed during analysis: " + path)


def merge_evidence(evidence, additional):
    """Never replace a prior binding with a later digest of the same file."""
    for path, digest in additional.items():
        require(path not in evidence or evidence[path] == digest,
                "Consumed file changed between audit scopes: " + path)
    evidence.update(additional)


def exact_array(left, right, description):
    require(left.shape == right.shape, "Shape mismatch: " + description)
    equal = np.array_equal(left, right, equal_nan=True) if left.dtype.kind in "fc" else np.array_equal(left, right)
    require(equal, "Value/identity mismatch: " + description)


def statistics(y, prediction):
    """Never remove an invalid prediction from a prespecified valid mask."""
    y, prediction = np.asarray(y, dtype=np.float64), np.asarray(prediction, dtype=np.float64)
    require(y.ndim == prediction.ndim == 1 and y.shape == prediction.shape, "Scoring vector shape differs")
    require(np.isfinite(y).all() and np.isfinite(prediction).all(), "Non-finite value on fixed valid support")
    n = len(y)
    result = dict(count=n, mae=np.nan, rmse=np.nan, bias=np.nan, abs_bias=np.nan,
                  spearman=np.nan, pearson=np.nan, correlation_status="no_valid_residues")
    if n == 0:
        return result
    difference = prediction - y
    result.update(mae=float(np.abs(difference).mean()), rmse=float(np.sqrt(np.square(difference).mean())),
                  bias=float(difference.mean()), abs_bias=float(abs(difference.mean())))
    if n < 4:
        result["correlation_status"] = "fewer_than_4_valid_residues"
    elif np.ptp(y) == 0:
        result["correlation_status"] = "constant_target"
    elif np.ptp(prediction) == 0:
        result["correlation_status"] = "constant_prediction"
    else:
        result.update(spearman=float(spearmanr(y, prediction).statistic),
                      pearson=float(pearsonr(y, prediction).statistic), correlation_status="defined")
        require(np.isfinite(result["spearman"]) and np.isfinite(result["pearson"]), "Unexpected undefined correlation")
    return result


def compare_saved(expected, saved, keys, description):
    require(not expected.duplicated(keys).any() and not saved.duplicated(keys).any(), "Duplicate rows: " + description)
    a, b = expected.set_index(keys).sort_index(), saved.set_index(keys).sort_index()
    require(a.index.equals(b.index), "Saved/recomputed record support differs: " + description)
    for name in ("count", "mae", "rmse", "bias", "spearman", "pearson"):
        x, y = a[name].to_numpy(float), b[name].to_numpy(float)
        require(np.array_equal(np.isnan(x), np.isnan(y)), "Saved NaN support differs: " + description + "/" + name)
        tolerance = 1e-10 if name in ("count", "spearman") else 3e-6
        require(np.allclose(x, y, atol=tolerance, rtol=tolerance, equal_nan=True), "Saved score differs: " + description + "/" + name)


def score_record(info, variant, target, valid, prediction, metric_indices=(0, 1)):
    rows = []
    for m in metric_indices:
        mask = valid[:, m]
        require(mask.dtype.kind == "b", "Valid mask is not boolean")
        p = prediction[:, m] if prediction.ndim == 2 else prediction
        rows.append(dict(sample_id=info["sample_id"], dataset=info["dataset"], fold=int(info["fold"]),
                         biological_system=info["biological_system"], variant=variant, metric=METRICS[m],
                         **statistics(target[mask, m], p[mask])))
    return rows


def baseline_preflight():
    """Verify labels, membership, fitted controls and released BF2; no new fits."""
    evidence = {}
    freeze_path = PAYLOAD / "training_freeze.json"
    bind(freeze_path, evidence, PRIMARY_SHA)
    frozen = read_json(freeze_path)
    require(frozen["records"] == 3063 and frozen["source_records"] == SOURCES, "Primary cohort changed")
    require(frozen["valid_counts"] == [961321, 824988], "Supervised support changed")
    for table_name in ("input_hashes", "execution_hashes"):
        for name, digest in frozen[table_name].items():
            bind(portable(PAYLOAD, name), evidence, digest)
    all_path = PAYLOAD / "splits/all_ready.jsonl"
    bind(all_path, evidence)
    records = read_lines(all_path)
    lookup = {r["sample_id"]: r for r in records}
    require(len(records) == len(lookup) == 3063, "Training record identity count differs")
    bind(STATIC / "training_schema.json", evidence)
    schema = read_json(STATIC / "training_schema.json")
    infos = schema["records"]
    require(len(infos) == 3063 and {i["sample_id"] for i in infos} == set(lookup), "Static schema membership differs")
    require(pd.Series([i["dataset"] for i in infos]).value_counts().to_dict() == SOURCES, "Static source membership differs")
    bind(PIPE / "static_feature_definitions.py", evidence, schema["feature_source_sha256"])
    feature_path = STATIC / "training_features.npz"
    bind(feature_path, evidence)
    with np.load(feature_path, allow_pickle=False) as z:
        target, valid, indices, fold = (z[k].copy() for k in ("target", "valid", "record_index", "fold"))
    require(target.shape == valid.shape == (961321, 2) and valid.dtype.kind == "b", "Static target/mask schema differs")
    require(indices.shape == fold.shape == (len(target),), "Static position membership shape differs")
    require(valid.sum(axis=0).tolist() == schema["valid_counts"] == frozen["valid_counts"], "Valid support count differs")
    require(np.isfinite(target[valid]).all(), "Non-finite valid static target")
    # Frozen scientific parsers are used only to read source identity/masks.
    # Importing their package defines model classes but never instantiates one.
    sys.path.insert(0, str(PAYLOAD / "src"))
    from flexppi.data.structures import parse_pdb, target_mask
    from flexppi.data.labels import load_label_file
    structures = {}
    cursor = 0
    for index, info in enumerate(infos):
        row = lookup[info["sample_id"]]
        require(info["record_index"] == index and info["start"] == cursor and info["stop"] > cursor, "Static record offsets are not contiguous")
        require(info["dataset"] == row["dataset"] and info["fold"] == row["metadata"]["fold"]
                and info["biological_system"] == row["metadata"]["biological_system"], "Static source/fold/system mismatch")
        ip, lp = local_input(row, "structure_path"), local_input(row, "label_path")
        bind(ip, evidence, info["input_sha256"])
        bind(lp, evidence, info["label_sha256"])
        structure = parse_pdb(ip)
        values, variances, mask = load_label_file(lp)
        require(values.shape == variances.shape == mask.shape == (structure.num_residues, 3), "Label/identity dimensions differ")
        require(not mask[:, 2].any(), "Unexpected supervised SASA")
        require(set(row["target_chains"]).issubset(set(structure.chain_ids)), "Unknown target chain")
        reference_mask = mask & structure.frame_mask[:, None] & target_mask(structure, row["target_chains"])[:, None]
        span = slice(info["start"], info["stop"])
        require(structure.num_residues == info["stop"] - info["start"], "Static length mismatch")
        exact_array(target[span], values[:, :2], "static/source target")
        exact_array(valid[span], reference_mask[:, :2], "static/source valid mask")
        require((indices[span] == index).all() and (fold[span] == info["fold"]).all(), "Residue record/fold mapping differs")
        structures[info["sample_id"]] = structure
        cursor = info["stop"]
    require(cursor == len(target), "Trailing static residues are unmapped")
    for f in range(5):
        p = PAYLOAD / f"splits/fold_{f}_test.jsonl"
        bind(p, evidence)
        test = read_lines(p)
        require({r["sample_id"] for r in test} == {i["sample_id"] for i in infos if i["fold"] == f}, "Held-out fold membership differs")
        require(all(r == lookup[r["sample_id"]] for r in test), "Held-out record contents differ")
    system_folds = pd.DataFrame(infos).groupby(["dataset", "biological_system"]).fold.nunique()
    require((system_folds == 1).all(), "A source-specific biological system spans held-out folds")
    run_path = STATIC / "run_audit.json"
    bind(run_path, evidence)
    run = read_json(run_path)
    require(run["status"] == "complete" and run["no_test_tuning"] is True and run["training_freeze_sha256"] == PRIMARY_SHA, "Static audit not complete/current")
    bind(PIPE / "run_corrected_static.py", evidence, run["code_sha256"])
    require(len(run["jobs"]) == 10 and len({(j["fold"], j["metric"]) for j in run["jobs"]}) == 10, "Static fit audit job count differs")
    for job in run["jobs"]:
        f, m = int(job["fold"]), METRICS.index(job["metric"])
        require(f in range(5), "Unexpected static fold")
        require(job["train_count"] == int(((fold != f) & valid[:, m]).sum()) and job["test_count"] == int(((fold == f) & valid[:, m]).sum()), "Static train/test support differs")
        bind(STATIC / f"static_fold{f}_{METRICS[m]}.joblib", evidence, job["checkpoint_sha256"])
    baseline_arrays = {}
    scored = []
    for name in ("static_sequence_complex", "aa_only"):
        p = STATIC / (name + "_oof.npz")
        bind(p, evidence)
        with np.load(p, allow_pickle=False) as z:
            prediction = z["prediction"].copy()
        require(prediction.shape == target.shape and np.isfinite(prediction[valid]).all(), "Invalid static prediction on fixed mask")
        baseline_arrays[name] = prediction
        for info in infos:
            span = slice(info["start"], info["stop"])
            scored += score_record(info, name, target[span], valid[span], prediction[span])
    bind(STATIC / "oof_per_record_metrics.csv", evidence)
    saved = csv(STATIC / "oof_per_record_metrics.csv")
    require(set(saved.role) == {"oof"}, "Unexpected static scoring role")
    compare_saved(pd.DataFrame(scored).query("count > 0"), saved,
                  ["sample_id", "dataset", "fold", "variant", "metric"], "static OOF")
    membership = family_preflight(infos, evidence)
    bf_prediction = bf_preflight(infos, structures, target, valid, evidence)
    baseline_arrays["backflip2_official"] = bf_prediction
    for info in infos:
        span = slice(info["start"], info["stop"])
        scored += score_record(info, "backflip2_official", target[span], valid[span], bf_prediction[span], (0,))
    saved_bf = csv(BF / "per_record_metrics.csv")
    compare_saved(pd.DataFrame(scored).query("variant == 'backflip2_official' and count > 0"),
                  saved_bf.query("role == 'training' and protocol == 'ordinary_oof' and variant == 'backflip2_official'"),
                  ["sample_id", "dataset", "fold", "variant", "metric"], "released BF2 RMSF")
    preflight = dict(status="prepared; final neural audit and complete OOF predictions pending",
                     final_analysis_complete=False, final_neural_scores_emitted=False, figures_created=0,
                     actual_records=3063, source_records=SOURCES, heldout_records=[sum(i["fold"] == f for i in infos) for f in range(5)],
                     valid_residues=valid.sum(axis=0).tolist(), supervised_metrics=list(METRICS),
                     baseline_verified_record_metric_rows=len(scored),
                     baseline_correlation_support=pd.DataFrame(scored).groupby(["variant", "metric"]).agg(
                         records=("sample_id", "size"), spearman_defined=("spearman", "count"), pearson_defined=("pearson", "count")).reset_index().to_dict("records"),
                     family_membership=membership, input_and_baseline_sha256=evidence, new_fits=0, new_forward_calls=0,
                     uncertainty="No bootstrap, CI or p value; external single family is not scored here.")
    return dict(infos=infos, records=lookup, target=target, valid=valid, baseline_arrays=baseline_arrays,
                baseline_rows=scored, structures=structures, evidence=evidence, preflight=preflight)


def family_preflight(infos, evidence):
    for name in ("record_overlap_audit.csv", "chain_overlap_audit.csv", "sequence_audit_summary.json", "strict_overlap_validation.json"):
        bind(SEQ / name, evidence)
    record, chain = csv(SEQ / "record_overlap_audit.csv"), csv(SEQ / "chain_overlap_audit.csv")
    summary, strict = read_json(SEQ / "sequence_audit_summary.json"), read_json(SEQ / "strict_overlap_validation.json")
    require(len(record) == 3095 and not record.duplicated(["role", "sample_id"]).any(), "Family record membership differs")
    training, external = record.query("role == 'training'"), record.query("role == 'external'")
    require(len(training) == 3063 and len(external) == 32 and external.component_id.nunique() == 1, "Family cohort count differs")
    info_df = pd.DataFrame(infos)
    joined = training.merge(info_df[["sample_id", "dataset", "fold", "biological_system"]], on="sample_id", validate="one_to_one", suffixes=("_a", "_b"))
    require(len(joined) == 3063 and (joined.dataset_a == joined.dataset_b).all() and (joined.fold_a == joined.fold_b).all(), "Family/OOF source or fold differs")
    require(len(chain) == 5068 and not chain.duplicated(["role", "sample_id", "chain"]).any(), "Full-chain identity count differs")
    chained = chain.query("role == 'training'").merge(info_df[["sample_id", "dataset", "fold", "biological_system"]], on="sample_id", validate="many_to_one", suffixes=("_a", "_b"))
    require(len(chained) == int((chain.role == "training").sum()), "Chain training record unmatched")
    for name in ("dataset", "fold", "biological_system"):
        require((chained[name + "_a"] == chained[name + "_b"]).all(), "Chain metadata mismatch: " + name)
    require(record.component_id.nunique() == 1664 and summary["components"] == 1664, "Family component count differs")
    require(record.groupby("component_id").size().max() == summary["largest_component_records"] == 664, "Family maximum component differs")
    require(training.purged_training_member.dtype.kind == "b", "Purge membership is not boolean")
    purged = training[training.purged_training_member]
    retained = training[~training.purged_training_member]
    require(len(purged) == summary["purged_training_records"] == 632 and len(retained) == summary["retained_training_records"] == 2431, "Purge count differs")
    require(set(purged.component_id) == set(external.component_id), "Purge is not the external-connected component")
    retained_counts = retained.dataset.value_counts().to_dict()
    require(retained_counts == summary["retained_by_source"] == {"atlas": 1866, "dynarepo": 565}, "Retained source composition differs")
    strict_path = PIPE / "strict_external_payload/strict_external_freeze.json"
    bind(strict_path, evidence, STRICT_SHA)
    strict_freeze = read_json(strict_path)
    bind(SEQ / "all_chain_similarity.tsv", evidence, strict_freeze["sequence_search_sha256"])
    require(strict_freeze["primary_freeze_sha256"] == PRIMARY_SHA and strict_freeze["training_records"] == 2431, "Strict freeze membership differs")
    actual_retained_path = SEQ / "purged_training_local.jsonl"
    bind(actual_retained_path, evidence)
    require({r["sample_id"] for r in read_lines(actual_retained_path)} == set(retained.sample_id), "Retained actual records differ")
    require(strict["status"] == "passed" and strict["strict_training_records"] == 2431 and len(strict["chains"]) == strict["external_chain_count"] == 95, "Strict overlap check incomplete")
    require(strict["maximum_external_chain_hits_after_purge"] == 0 and all(c["purged_training_hits"] == 0 for c in strict["chains"]), "Strict search hit remains")
    ext_chain = chain.query("role == 'external'")
    require(len(ext_chain) == 95 and set(zip(ext_chain.sample_id, ext_chain.chain)) == {(c["sample_id"], c["chain"]) for c in strict["chains"]}, "95 external chain identities differ")
    return dict(training_records=3063, external_records=32, chain_observations=5068, components=1664,
                source_specific_biological_systems=int(info_df.groupby(["dataset", "biological_system"]).ngroups),
                external_components=1, purged_records=632, retained_records=2431, retained_by_source=retained_counts,
                strict_external_chain_count=95, strict_maximum_search_hits=0,
                limitation=summary["limitations"])


def bf_preflight(infos, structures, target, valid, evidence):
    for p in (BF / "scoring_audit.json", BF / "per_record_metrics.csv", BF / "training_aligned.npz", BF_INPUT / "input_freeze.json", BF_INPUT / "record_manifest.json", BF_PRED / "inference_metadata.json"):
        bind(p, evidence)
    audit, freeze = read_json(BF / "scoring_audit.json"), read_json(BF_INPUT / "input_freeze.json")
    require(audit["status"] == "complete" and audit["records"] == 3095 and audit["chains"] == 5068 and audit["inference_failures"] == 0, "BF2 complete-cohort audit missing")
    require(audit["full_identity_and_input_hash_verification"] is True and audit["no_fitted_calibration"] is True and audit["no_error_based_exclusions"] is True, "BF2 identity/protocol audit differs")
    require(audit["training_freeze_sha256"] == PRIMARY_SHA and audit["strict_freeze_sha256"] == STRICT_SHA, "BF2 freeze binding differs")
    bind(PIPE / "score_corrected_backflip2.py", evidence, audit["source_sha256"])
    bind(SEQ / "record_overlap_audit.csv", evidence, audit["search_audit_sha256"])
    bind(BF_INPUT / "residue_mapping.csv", evidence, freeze["mapping_sha256"])
    require(audit["mapping_sha256"] == freeze["mapping_sha256"], "BF2 mapping binding differs")
    bind(BF_INPUT / "input_manifest.csv", evidence, freeze["input_manifest_sha256"])
    metadata = read_json(BF_PRED / "inference_metadata.json")
    require(metadata == audit["official_inference_metadata"], "BF2 actual metadata differs")
    require(metadata["manifest_sha256"] == freeze["input_manifest_sha256"] and metadata["checkpoint_sha256"] == "1af6e13a0e0161b2b92fa6fb9cddd2dee8f3e71c284f5059c77d7a30b9c527b9", "BF2 released input/weight binding differs")
    require(metadata["source_commit"] == "887c6a73b216e3fe2ff9f551ae6e2c983c21b292" and metadata["covariance_unit"] == "nm2", "BF2 released code/unit differs")
    with np.load(BF / "training_aligned.npz", allow_pickle=False) as z:
        prediction = z["prediction"].copy()
        exact_array(z["target"], target[:, 0], "BF2/common target")
        exact_array(z["valid"], valid[:, 0], "BF2/common RMSF mask")
    require(prediction.shape == (len(target),) and np.isfinite(prediction).all(), "Invalid aligned BF2 shape/coverage")
    manifest = pd.read_csv(BF_INPUT / "input_manifest.csv", keep_default_na=False)
    require(len(manifest) == 5068 and manifest.input_id.is_unique, "BF2 chain manifest count differs")
    outputs = {r["input_id"]: r for r in audit["output_hashes"]}
    require(len(outputs) == 5068 and set(outputs) == set(manifest.input_id), "BF2 actual output hash membership differs")
    statuses = pd.DataFrame(read_lines(BF_PRED / "status.jsonl"))
    bind(BF_PRED / "status.jsonl", evidence)
    require(len(statuses) == 5068 and statuses.input_id.is_unique and set(statuses.input_id) == set(outputs) and (statuses.status == "complete").all(), "BF2 inference coverage differs")
    # Reassemble existing chain outputs and verify their physical conversion.
    # External outputs are hash-checked for completed-pass evidence, not scored.
    training_predictions = {}
    for row in manifest.to_dict("records"):
        key = row["input_id"]
        bind(portable(BF_INPUT, row["pdb_file"]), evidence, row["input_sha256"])
        p = BF_PRED / (key + ".npz")
        bind(p, evidence, outputs[key]["sha256"])
        require(outputs[key]["n_residues"] == row["n_residues"], "BF2 recorded chain length differs")
        if row["role"] == "training":
            with np.load(p, allow_pickle=False) as z:
                rmsf, covariance = z["rmsf_angstrom"].copy(), z["covariance_nm2"]
                require(rmsf.shape == (row["n_residues"],) and covariance.shape == (row["n_residues"], 3, 3), "BF2 actual prediction shape differs")
                require(np.isfinite(rmsf).all() and np.isfinite(covariance).all(), "BF2 raw prediction non-finite")
                trace = np.trace(covariance, axis1=-2, axis2=-1)
                require((trace >= 0).all() and np.allclose(rmsf, 10 * np.sqrt(trace), atol=2e-6, rtol=2e-6), "BF2 nm2-to-Angstrom conversion differs")
                training_predictions[key] = rmsf
    by_sample = dict(tuple(manifest.query("role == 'training'").groupby("sample_id", sort=False)))
    require(set(by_sample) == {i["sample_id"] for i in infos}, "BF2 training record membership differs")
    for info in infos:
        structure = structures[info["sample_id"]]
        entries = by_sample[info["sample_id"]]
        require(set(entries.chain) == set(structure.unique_chains()), "BF2 source chains differ")
        assembled = np.full(structure.num_residues, np.nan)
        visits = np.zeros(structure.num_residues, dtype=np.int32)
        for entry in entries.to_dict("records"):
            ids = structure.indices_for_chain(entry["chain"])
            require(len(ids) == entry["n_residues"], "BF2 parsed chain length differs")
            assembled[ids] = training_predictions[entry["input_id"]]
            visits[ids] += 1
        require((visits == 1).all(), "BF2 assembly contains missing/duplicate positions")
        exact_array(prediction[info["start"]:info["stop"]], assembled, "BF2 aligned/raw predictions")
    return prediction


def neural_gate(audit_path, export_root, weight_root, evidence):
    """No successful boolean alone can authorize final analysis."""
    require(audit_path is not None, "Final mode requires an explicit --neural-audit; use --preflight while pending")
    bind(audit_path, evidence)
    audit = read_json(audit_path)
    require(audit.get("status") == COMPLETE_STATUS and audit.get("final_prediction_audit_complete") is True, "Complete actual neural export audit is pending")
    auditor = Path(__file__).resolve().parent / "audit_corrected_neural_exports.py"
    bind(auditor, evidence, audit["auditor_sha256"])
    require(audit["actual_records"] == {"oof_per_arm": 3063, "external": 32, "strict_training": 2431}, "Audited cohort counts differ")
    require(audit["actual_valid_counts"] == {"oof": [961321, 824988], "external": [18748, 16265]}, "Audited target support differs")
    require(set(audit["scopes"]) == {"oof", "external", "strict"} and all(s["status"] == "verified" for s in audit["scopes"].values()), "Neural audit scopes incomplete")
    require(audit["new_fits"] == audit["new_forward_calls"] == 0, "Expected saved-prediction-only audit")
    tables = {}
    for name in ("input_sha256", "queue_state_sha256", "output_sha256"):
        require(isinstance(audit.get(name), dict) and len(audit[name]) > 0, "Audit SHA table missing: " + name)
        normalized = {str(Path(p).resolve()): d for p, d in audit[name].items()}
        require(len(normalized) == len(audit[name]), "Aliased audit paths")
        tables[name] = normalized
        for p, digest in normalized.items():
            bind(p, evidence, digest)
    require(len(tables["input_sha256"]) >= 6190, "Incomplete frozen-input audit bindings")
    for p in (PAYLOAD / "training_freeze.json", PIPE / "strict_external_payload/strict_external_freeze.json"):
        require(str(p.resolve()) in tables["input_sha256"], "Frozen scheme not bound in final audit")
    # All primary frozen inputs and executable files must be part of that audit.
    frozen = read_json(PAYLOAD / "training_freeze.json")
    for table in ("input_hashes", "execution_hashes"):
        for name, digest in frozen[table].items():
            p = str(portable(PAYLOAD, name))
            require(tables["input_sha256"].get(p) == digest, "Final audit does not bind current primary file: " + p)
    state_paths = [export_root / "queue_state.json", export_root / "strict_queue_state.json",
                   weight_root / "completed_members_audit.json", weight_root / "completed_strict_members_audit.json"]
    require(set(tables["queue_state_sha256"]) == {str(p.resolve()) for p in state_paths}, "Actual queue/collected-member audit paths differ")
    expected_members = {f"{m}/fold_{f}" for v, f, m in MEMBERS}
    require(set(audit["checkpoint_sha256"]) == expected_members and len(expected_members) == 16, "Final checkpoint binding membership differs")
    for member, digest in audit["checkpoint_sha256"].items():
        bind(weight_root / member / "epoch_040.pt", evidence, digest)
    directory = export_root / "evaluation"
    for name in ("evaluation_audit.json", "prediction_manifest.jsonl", "per_record_metrics.csv", "source_fold_macro_metrics.csv", "five_fold_mean_sample_sd.csv"):
        require(str((directory / name).resolve()) in tables["output_sha256"], "Consumed OOF metadata not bound by final audit")
    expected_caches = {str((directory / v / (i["sample_id"] + ".npz")).resolve())
                       for v in ("full", "intrachain") for i in read_json(STATIC / "training_schema.json")["records"]}
    bound_caches = {p for p in tables["output_sha256"] if Path(p).is_relative_to(directory.resolve()) and Path(p).suffix == ".npz"}
    require(bound_caches == expected_caches and len(expected_caches) == 6126, "OOF cache bindings incomplete or unexpected")
    return audit, tables["output_sha256"]


def neural_records(context, export_root, audit, output_hashes):
    directory = export_root / "evaluation"
    manifest = read_lines(directory / "prediction_manifest.jsonl")
    by_key = {(r["variant"], r["sample_id"]): r for r in manifest}
    expected = {(v, i["sample_id"]) for v in ("full", "intrachain") for i in context["infos"]}
    require(len(manifest) == len(by_key) == 6126 and set(by_key) == expected, "OOF manifest membership differs")
    rows = []
    for info in context["infos"]:
        span = slice(info["start"], info["stop"])
        reference = context["records"][info["sample_id"]]
        structure = context["structures"][info["sample_id"]]
        for variant in ("full", "intrachain"):
            entry = by_key[(variant, info["sample_id"])]
            require(entry["fold"] == info["fold"] and entry["dataset"] == info["dataset"] and entry["record"] == reference, "OOF record/fold/source differs")
            require(entry["checkpoint_sha256"] == audit["checkpoint_sha256"][f"{variant}/fold_{info['fold']}"]
                    and entry["input_sha256"] == info["input_sha256"] and entry["label_sha256"] == info["label_sha256"], "OOF current input/checkpoint differs")
            require(entry["cache_path"].replace("\\", "/") == REMOTE + f"/evaluation/{variant}/{info['sample_id']}.npz", "Unexpected OOF remote cache path")
            p = directory / variant / (info["sample_id"] + ".npz")
            bind(p, context["evidence"], output_hashes[str(p.resolve())])
            with np.load(p, allow_pickle=False) as z:
                require(z["target"].shape == z["valid"].shape == (structure.num_residues, 3) and z["valid"].dtype.kind == "b", "OOF label dimensions/mask type differ")
                require(z["mean"].shape == (structure.num_residues, 2), "OOF prediction dimensions differ")
                exact_array(z["target"][:, :2], context["target"][span], "OOF/static targets")
                exact_array(z["valid"][:, :2], context["valid"][span], "OOF/static valid mask")
                require(not z["valid"][:, 2].any(), "Unexpected SASA supervision")
                for key, values in (("chain", structure.chain_ids), ("author_residue_id", structure.author_residue_ids), ("amino_acid", list(structure.sequence)), ("aatype", structure.aatype)):
                    exact_array(z[key], np.asarray(values), "OOF residue identity/" + key)
                require(np.isfinite(z["mean"]).all() and (z["mean"] >= 0).all(), "Invalid OOF absolute-state output")
                rows += score_record(info, variant, context["target"][span], context["valid"][span], z["mean"])
    compare_saved(pd.DataFrame(rows).query("count > 0"), csv(directory / "per_record_metrics.csv"),
                  ["variant", "fold", "dataset", "sample_id", "metric"], "independently audited neural OOF")
    return rows


def finite_mean(values):
    values = np.asarray(values, dtype=np.float64)
    return float(values[np.isfinite(values)].mean()) if np.isfinite(values).any() else np.nan


def scopes(frame):
    """Each row participates once in source/all and fold/all descriptive scopes."""
    for dataset, group in list(frame.groupby("dataset", sort=True)) + [("all", frame)]:
        yield dataset, -1, group
        for f, subset in group.groupby("fold", sort=True):
            yield dataset, int(f), subset


def summarize_records(frame):
    macros, system_rows = [], []
    for (variant, metric), group in frame.groupby(["variant", "metric"], sort=True):
        for dataset, fold, part in scopes(group):
            systems = part.groupby(["dataset", "biological_system"], sort=True)[list(SCORES)].mean().reset_index()
            for estimand, units in (("record_macro", part), ("source_system_macro", systems)):
                row = dict(variant=variant, metric=metric, dataset=dataset, fold=fold, estimand=estimand,
                           records=len(part), biological_systems=len(systems), valid_residues=int(part["count"].sum()),
                           records_with_valid_residues=int((part["count"] > 0).sum()),
                           spearman_records_defined=int(np.isfinite(part.spearman).sum()), pearson_records_defined=int(np.isfinite(part.pearson).sum()))
                for score in SCORES:
                    row[score] = finite_mean(units[score])
                    row[score + "_units_defined"] = int(np.isfinite(units[score]).sum())
                macros.append(row)
        for (dataset, system), part in group.groupby(["dataset", "biological_system"], sort=True):
            row = dict(variant=variant, metric=metric, dataset=dataset, biological_system=system,
                       fold=int(part.fold.iloc[0]), records=len(part), valid_residues=int(part["count"].sum()))
            for score in SCORES:
                row[score] = finite_mean(part[score])
                row[score + "_records_defined"] = int(np.isfinite(part[score]).sum())
            system_rows.append(row)
    macro = pd.DataFrame(macros)
    five = []
    for keys, group in macro[macro.fold >= 0].groupby(["variant", "metric", "dataset", "estimand"], sort=True):
        require(set(group.fold) == set(range(5)) and len(group) == 5, "Five-fold support missing")
        row = dict(zip(["variant", "metric", "dataset", "estimand"], keys))
        row["folds"] = 5
        for score in SCORES:
            values = group[score].to_numpy(float)
            finite = values[np.isfinite(values)]
            row[score + "_defined_folds"] = len(finite)
            row[score + "_fold_mean"] = finite_mean(finite)
            row[score + "_fold_sample_sd"] = float(finite.std(ddof=1)) if len(finite) >= 2 else np.nan
        five.append(row)
    return macro, pd.DataFrame(system_rows), pd.DataFrame(five)


def paired_records(frame):
    keys = ["sample_id", "dataset", "fold", "biological_system", "metric"]
    full = frame[frame.variant == "full"]
    results = []
    for comparator in ("intrachain", "static_sequence_complex", "aa_only", "backflip2_official"):
        baseline = frame[frame.variant == comparator]
        pair = full.merge(baseline, on=keys, how="inner", suffixes=("_full", "_comparator"), validate="one_to_one")
        expected = 3063 if comparator == "backflip2_official" else 6126
        require(len(pair) == expected, "Paired record coverage differs: " + comparator)
        require((pair.count_full == pair.count_comparator).all(), "Pair valid residue support differs")
        pair["comparator"] = comparator
        for score in SCORES:
            defined_full = np.isfinite(pair[score + "_full"])
            defined_comparator = np.isfinite(pair[score + "_comparator"])
            common = defined_full & defined_comparator
            pair[score + "_full_defined"] = defined_full
            pair[score + "_comparator_defined"] = defined_comparator
            pair[score + "_both_defined"] = common
            # Signed bias has no monotone improvement direction; abs_bias does.
            if score != "bias":
                delta = pair[score + "_full"] - pair[score + "_comparator"]
                if score in ("mae", "rmse", "abs_bias"):
                    delta = -delta
                pair[score + "_improvement"] = delta.where(common, np.nan)
        results.append(pair)
    return pd.concat(results, ignore_index=True)


def summarize_pairs(pairs):
    rows = []
    for (comparator, metric), group in pairs.groupby(["comparator", "metric"], sort=True):
        for dataset, fold, part in scopes(group):
            for score in ("mae", "rmse", "abs_bias", "spearman", "pearson"):
                f, b = part[score + "_full"].to_numpy(float), part[score + "_comparator"].to_numpy(float)
                common = np.isfinite(f) & np.isfinite(b)
                support = dict(comparator=comparator, metric=metric, dataset=dataset, fold=fold, score=score,
                               records=len(part), full_records_defined=int(np.isfinite(f).sum()),
                               comparator_records_defined=int(np.isfinite(b).sum()), paired_records_defined=int(common.sum()),
                               full_record_mean_own_support=finite_mean(f), comparator_record_mean_own_support=finite_mean(b))
                common_part = part.loc[common].copy()
                system_units = common_part.groupby(["dataset", "biological_system"])[[score + "_full", score + "_comparator", score + "_improvement"]].mean()
                total_systems = int(part.groupby(["dataset", "biological_system"]).ngroups)
                for estimand, units in (("record_macro", common_part), ("source_system_macro", system_units)):
                    rows.append(dict(**support, estimand=estimand, biological_systems=total_systems,
                                     paired_systems_defined=len(system_units), paired_units_defined=len(units),
                                     full_mean_common_support=finite_mean(units[score + "_full"]),
                                     comparator_mean_common_support=finite_mean(units[score + "_comparator"]),
                                     mean_paired_improvement=finite_mean(units[score + "_improvement"])))
    return pd.DataFrame(rows)


def all_record_points(frame, infos):
    points = pd.DataFrame(infos)[["sample_id", "dataset", "fold", "biological_system", "start", "stop", "input_sha256", "label_sha256"]].copy()
    require(len(points) == 3063 and points.sample_id.is_unique, "Main figure point census differs")
    values = ["count", *SCORES, "correlation_status"]
    for (variant, metric), group in frame.groupby(["variant", "metric"], sort=True):
        require(len(group) == 3063 and set(group.sample_id) == set(points.sample_id), "Main figure record lost")
        renamed = group[["sample_id", *values]].rename(columns={s: f"{variant}__{metric}__{s}" for s in values})
        points = points.merge(renamed, on="sample_id", how="left", validate="one_to_one")
    return points


def json_write(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    # allow_nan=False prevents an unnoticed non-finite support count in QA.
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--neural-audit", type=Path, help="Explicit successfully completed audit.json; never inferred from latest")
    parser.add_argument("--export-root", type=Path, default=PIPE / "neural_exports")
    parser.add_argument("--weight-root", type=Path, default=PIPE / "trained_checkpoints")
    parser.add_argument("--output-dir", type=Path, help="New dedicated directory inside qa/corrected_main_comparison/analyses")
    args = parser.parse_args()
    evidence = {}
    report = dict(checked_utc=datetime.now(timezone.utc).isoformat(), analyzer_sha256=sha(__file__),
                  status="pending; no completed final analysis", final_analysis_complete=False,
                  new_fits=0, new_forward_calls=0, figures_created=0, bootstrap_replicates=0)
    try:
        audit = None
        if not args.preflight:
            # The gate precedes score computation and output-directory creation.
            audit, output_hashes = neural_gate(args.neural_audit, args.export_root.resolve(), args.weight_root.resolve(), evidence)
        context = baseline_preflight()
        merge_evidence(evidence, context["evidence"])
        if args.preflight:
            report.update(context["preflight"])
            report["analyzer_sha256"] = sha(__file__)
            check_snapshot(evidence)
            json_write(QA / "corrected_main_comparison_preflight.json", report)
            print(report["status"])
            return
        context["evidence"] = evidence
        rows = context["baseline_rows"] + neural_records(context, args.export_root.resolve(), audit, output_hashes)
        frame = pd.DataFrame(rows)
        require(len(frame) == 3063 * 9 and not frame.duplicated(["sample_id", "variant", "metric"]).any(), "Final complete per-record comparison census differs")
        require(set(frame.dataset) == set(SOURCES) and set(frame.metric) == set(METRICS), "Final source/target names differ")
        macro, systems, five = summarize_records(frame)
        pairs = paired_records(frame)
        pair_macro = summarize_pairs(pairs)
        points = all_record_points(frame, context["infos"])
        check_snapshot(evidence)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
        output_root = (QA / "corrected_main_comparison/analyses").resolve()
        output = args.output_dir.resolve() if args.output_dir else output_root / stamp
        require(output.is_relative_to(output_root) and output != output_root and not output.exists(), "Refuse to overwrite or write outside dedicated final QA analyses")
        output.mkdir(parents=True)
        outputs = {"per_record_metrics.csv": frame, "all_3063_record_points.csv": points,
                   "source_fold_record_and_system_macros.csv": macro, "per_source_biological_system_metrics.csv": systems,
                   "five_fold_mean_sample_sd.csv": five, "paired_record_improvements.csv": pairs,
                   "paired_common_support_macros.csv": pair_macro}
        for name, data in outputs.items():
            data.to_csv(output / name, index=False, na_rep="")
        report.update(status="complete; all actual OOF records compared on identical targets and masks", final_analysis_complete=True,
                      final_neural_scores_emitted=True, neural_audit_path=str(args.neural_audit.resolve()),
                      neural_audit_sha256=sha(args.neural_audit), actual_records=3063, per_record_metric_rows=len(frame),
                      paired_record_metric_rows=len(pairs), input_and_consumed_output_sha256=evidence,
                      generated_csv_sha256={name: sha(output / name) for name in outputs},
                      estimands=["record_macro", "source_system_macro", "five-fold mean and sample SD"],
                      improvement_direction="Positive favors full: comparator-full for MAE/RMSE/absolute bias; full-comparator for correlations. Signed bias is retained without an improvement claim.",
                      uncertainty="Descriptive comparisons only. No bootstrap/CI/p; no external-family independence claim.")
        # Final race check also binds the analyzer and completed audit unchanged.
        check_snapshot(evidence)
        require(report["analyzer_sha256"] == sha(__file__), "Analyzer changed during execution")
        json_write(output / "analysis_audit.json", report)
        print(report["status"])
        print("Dedicated output:", output)
    except Exception as error:
        report.update(status="pending or failed; no final analysis accepted", final_analysis_complete=False,
                      final_neural_scores_emitted=False, error_type=type(error).__name__, error=str(error))
        target = QA / ("corrected_main_comparison_preflight.json" if args.preflight else "corrected_main_comparison_pending.json")
        json_write(target, report)
        raise


if __name__ == "__main__":
    main()
