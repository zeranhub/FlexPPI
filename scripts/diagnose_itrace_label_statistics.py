"""Read-only source-label and held-out amino-acid-prior diagnosis.

No checkpoint is loaded and no model training or inference is performed.
Only this script's output JSON is written; original inputs stay unchanged.
"""
from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from flexppi.constants import RESTYPES
from flexppi.data.labels import load_label_file
from flexppi.data.records import remap_path
from flexppi.data.structures import parse_pdb, target_mask
from flexppi.evaluation import _correlation, _spearman
from package_atlas_hchi_retraining import LOCAL_MAP

SPLIT = ROOT / "output/atlas_hchi_retraining_v1_linux/payload/data/splits/state_rmsf_hchi_v2_atlas_hchi_v1"
EVALUATION = ROOT / "artifacts/performance_first_v2_atlas_hchi_v1/evaluation"
OUT = ROOT / "output/itrace_performance_diagnosis/label_statistics.json"
METRICS = ["rmsf", "chi_entropy"]


def describe(values):
    a = np.asarray(values, dtype=np.float64)
    a = a[np.isfinite(a)]
    if not len(a):
        return {"count": 0}
    q = np.quantile(a, [0.0, 0.25, 0.5, 0.75, 1.0])
    return {"count": int(len(a)), "mean": float(a.mean()),
            "std": float(a.std()), "min": float(q[0]), "q25": float(q[1]),
            "median": float(q[2]), "q75": float(q[3]),
            "iqr": float(q[3] - q[1]), "max": float(q[4])}


def corr_report(prediction, target):
    error = prediction - target
    return {"pearson": _correlation(prediction, target),
            "spearman": _spearman(prediction, target),
            "mae": float(np.abs(error).mean()),
            "rmse": float(np.sqrt(np.square(error).mean()))}


def macro(reports):
    return {key: float(np.mean([row[key] for row in reports if row[key] is not None]))
            for key in ("pearson", "spearman", "mae", "rmse")}


def main():
    start = time.monotonic()
    rows = [json.loads(line) for line in (SPLIT / "all_ready.jsonl").read_text(encoding="utf-8").splitlines() if line]
    assert len(rows) == 3068
    evaluations = [json.loads((EVALUATION / f"fold_{fold}_test.json").read_text()) for fold in range(5)]
    samples = []
    for index, row in enumerate(rows, 1):
        structure = parse_pdb(remap_path(row["structure_path"], LOCAL_MAP))
        values, variances, mask = load_label_file(remap_path(row["label_path"], LOCAL_MAP))
        assert values.shape[0] == structure.num_residues
        valid = mask & np.isfinite(values) & structure.frame_mask[:, None] & target_mask(structure, row["target_chains"])[:, None]
        samples.append({"dataset": row["dataset"], "sample_id": row["sample_id"],
                        "fold": int(row["metadata"]["fold"]), "n_residues": structure.num_residues,
                        "n_chains": len(structure.unique_chains()), "aatype": structure.aatype,
                        "values": values, "variances": variances, "valid": valid})
        if index % 250 == 0 or index == len(rows):
            print(f"Loaded {index}/{len(rows)} systems in {time.monotonic() - start:.1f}s", flush=True)

    datasets = sorted({row["dataset"] for row in samples})
    output = {"scope": "3068 systems in unchanged five-fold pool; geometry and target-chain masks applied",
              "label_units": {"rmsf": "Angstrom", "chi_entropy": "nats"},
              "baseline_definition": "For each held-out fold, predict only the training-fold residue mean for each amino-acid type. Pearson/Spearman calculated separately per system; then equal-system and equal-fold means reported.",
              "variance_warning": "Stored variances arise from dataset-specific replicas/windows and are not directly comparable uncertainty estimates.",
              "datasets": {}, "baselines": {}, "training_hchi_priors_max_absolute_difference": 0.0}
    for dataset in datasets:
        selected = [sample for sample in samples if sample["dataset"] == dataset]
        entry = {"systems": len(selected), "n_residues_total": sum(x["n_residues"] for x in selected),
                 "system_length": describe([x["n_residues"] for x in selected]),
                 "chain_count": describe([x["n_chains"] for x in selected]), "metrics": {}}
        for m, metric in enumerate(METRICS):
            local_values = [x["values"][x["valid"][:, m], m].astype(np.float64) for x in selected]
            local_variances = [x["variances"][x["valid"][:, m], m].astype(np.float64) for x in selected]
            local_aa = [x["aatype"][x["valid"][:, m]] for x in selected]
            values = np.concatenate(local_values)
            variances = np.concatenate(local_variances)
            aatypes = np.concatenate(local_aa)
            per_system = [{"sample_id": x["sample_id"], "fold": x["fold"],
                           "label_count": len(y), "target_mean": float(y.mean()),
                           "target_std": float(y.std()), "target_iqr": float(np.quantile(y, .75) - np.quantile(y, .25)),
                           "label_variance_median": float(np.median(v))}
                          for x, y, v in zip(selected, local_values, local_variances) if len(y)]
            entry["metrics"][metric] = {"pooled_target": describe(values),
                "pooled_label_variance": describe(variances),
                "zero_label_variance_fraction": float(np.mean(variances == 0)),
                "per_system_target_std": describe([x["target_std"] for x in per_system]),
                "per_system_target_iqr": describe([x["target_iqr"] for x in per_system]),
                "per_system_target_mean": describe([x["target_mean"] for x in per_system]),
                "per_system_label_variance_median": describe([x["label_variance_median"] for x in per_system]),
                "amino_acid_targets": {aa: describe(values[aatypes == i]) for i, aa in enumerate(RESTYPES)},
                "per_system": per_system}
        output["datasets"][dataset] = entry

    for scope in ("all_sources_training_priors", "same_source_training_priors"):
        output["baselines"][scope] = {}
        for m, metric in enumerate(METRICS):
            results = defaultdict(list)
            for fold in range(5):
                for dataset in datasets:
                    training = [x for x in samples if x["fold"] != fold and (scope == "all_sources_training_priors" or x["dataset"] == dataset)]
                    sums = np.zeros(len(RESTYPES) + 1, dtype=np.float64)
                    counts = np.zeros(len(RESTYPES) + 1, dtype=np.int64)
                    for sample in training:
                        valid = sample["valid"][:, m]
                        for aa in np.unique(sample["aatype"][valid]):
                            selected = valid & (sample["aatype"] == aa)
                            sums[aa] += float(sample["values"][selected, m].sum())
                            counts[aa] += int(selected.sum())
                    priors = np.divide(sums, counts, out=np.full_like(sums, sums.sum() / counts.sum()), where=counts > 0)
                    if scope == "all_sources_training_priors" and metric == "chi_entropy":
                        stored = evaluations[fold]["metrics"][metric]["amino_acid_corrected"]["training_priors"]
                        difference = max(abs(float(priors[RESTYPES.index(aa)]) - value) for aa, value in stored.items())
                        output["training_hchi_priors_max_absolute_difference"] = max(output["training_hchi_priors_max_absolute_difference"], difference)
                        assert difference < 1e-7, (fold, difference)
                    test = [x for x in samples if x["fold"] == fold and x["dataset"] == dataset]
                    for sample in test:
                        valid = sample["valid"][:, m]
                        target = sample["values"][valid, m].astype(np.float64)
                        if len(target):
                            report = corr_report(priors[sample["aatype"][valid]], target)
                            report.update({"sample_id": sample["sample_id"], "fold": fold})
                            results[dataset].append(report)
            output["baselines"][scope][metric] = {}
            for dataset, reports in results.items():
                folds = [macro([x for x in reports if x["fold"] == fold]) for fold in range(5)]
                output["baselines"][scope][metric][dataset] = {
                    "systems": len(reports), "macro_equal_system": macro(reports),
                    "fold_macro": folds, "equal_fold_mean": macro(folds),
                    "fold_sample_std": {key: float(np.std([x[key] for x in folds], ddof=1)) for key in folds[0]},
                    "per_system": reports}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(output, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}; elapsed {time.monotonic() - start:.1f}s", flush=True)
    for source in datasets:
        print(source, "n", output["datasets"][source]["systems"])
        for metric in METRICS:
            data = output["datasets"][source]["metrics"][metric]
            print(metric, "targets", data["pooled_target"], "within-system target std median", data["per_system_target_std"]["median"], "variance median", data["pooled_label_variance"]["median"])
            print("AA baseline", output["baselines"]["all_sources_training_priors"][metric][source]["equal_fold_mean"])


if __name__ == "__main__":
    main()
