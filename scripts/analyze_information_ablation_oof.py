"""Build the per-residue table for one information-ablation arm.

Reads only the provenance-checked OOF caches written by
``export_information_ablation_oof.py`` and writes the same per-residue schema
the frozen analysis uses, including training-fold-only amino-acid priors.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from export_information_ablation_oof import EXPECTED_COUNTS, EXPECTED_SOURCES, VARIANTS  # noqa: E402
from export_v2_oof_predictions import METRICS, load_cache, sha256  # noqa: E402

ARM_STATUS = ("complete", "partial_smoke")


def read_arm(folder: Path, allow_partial: bool) -> tuple:
    summary = json.loads((folder / "run_summary.json").read_text(encoding="utf-8"))
    if summary["status"] not in ARM_STATUS:
        raise ValueError("Arm export did not finish: " + repr(summary.get("error")))
    if summary["partial"] and not allow_partial:
        raise ValueError("Partial arm export; re-run the exporter or pass --allow-partial")
    manifest = [json.loads(line) for line in (folder / "prediction_manifest.jsonl").read_text().splitlines()]
    if len({row["sample_id"] for row in manifest}) != len(manifest):
        raise ValueError("Duplicate sample in arm manifest")
    frames = []
    for i, row in enumerate(manifest):
        if row["status"] != "complete":
            raise ValueError("Incomplete cache entry: " + row["sample_id"])
        arrays = load_cache(Path(row["cache_path"]), row["identity"])
        meta = json.loads(str(arrays["metadata_json"]))
        if meta["protocol"]["variant"] != summary["protocol"]["variant"]:
            raise ValueError("Mixed arm protocols in caches")
        record = row["record"]
        frame = pd.DataFrame(dict(
            sample_id=row["sample_id"], dataset=record["dataset"], fold=row["fold"],
            chain=arrays["chain"], author_residue_id=arrays["author_residue_id"],
            residue_index=arrays["residue_index"], amino_acid=arrays["amino_acid"],
            biological_system=record["metadata"].get("biological_system", row["sample_id"]),
            cluster_id=record["metadata"].get("cluster_id", ""),
            source_group=record.get("source_group", ""),
            whole_system_length=len(arrays["chain"]), chain_count=len(set(arrays["chain"]))))
        for j, metric in enumerate(METRICS):
            for prefix, array in [("valid", "valid"), ("target", "target"),
                                  ("pred", "mean"), ("model_variance", "model_variance"),
                                  ("label_variance", "label_variance")]:
                frame[f"{prefix}_{metric}"] = arrays[array][:, j]
        frames.append(frame)
        if (i + 1) % 500 == 0:
            print(f"Read {i + 1} arm caches", flush=True)
    data = pd.concat(frames, ignore_index=True)
    if data.duplicated(["sample_id", "chain", "author_residue_id"]).any():
        raise ValueError("Duplicate residue identifiers")
    if not summary["partial"]:
        observed = data.groupby("dataset").sample_id.nunique().to_dict()
        if observed != EXPECTED_SOURCES:
            raise ValueError("Source inventory differs from the frozen run")
        for metric, count in EXPECTED_COUNTS.items():
            if int(data[f"valid_{metric}"].sum()) != count:
                raise ValueError("Full residue coverage mismatch: " + metric)
    return data, manifest, summary


def attach_training_priors(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in METRICS:
        data[f"aa_prior_{metric}"] = np.nan
        for fold in range(5):
            train = data.loc[(data.fold != fold) & data[f"valid_{metric}"]]
            means = train.assign(y=train[f"target_{metric}"].astype(float)).groupby("amino_acid").y.mean()
            test = (data.fold == fold) & data[f"valid_{metric}"]
            data.loc[test, f"aa_prior_{metric}"] = data.loc[test, "amino_acid"].map(means)
            for aa, value in means.items():
                rows.append(dict(metric=metric, fold=fold, amino_acid=aa, training_mean=value))
        if data.loc[data[f"valid_{metric}"], f"aa_prior_{metric}"].isna().any():
            raise ValueError("Held-out amino-acid type has no training prior")
    return pd.DataFrame(rows)


def run(input_dir: Path, output_dir: Path, allow_partial: bool) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    data, manifest, summary = read_arm(input_dir, allow_partial)
    variant = summary["protocol"]["variant"]
    priors = attach_training_priors(data)
    data.to_csv(output_dir / "oof_per_residue.csv.gz", index=False, compression="gzip")
    priors.to_csv(output_dir / "training_fold_aa_priors.csv", index=False, encoding="utf-8")
    report = dict(status="complete", variant=variant, partial=bool(summary["partial"]),
                  n_systems=len(manifest), rows=len(data),
                  valid_residues={metric: int(data[f"valid_{metric}"].sum()) for metric in METRICS},
                  source_manifest_sha256=sha256(input_dir / "prediction_manifest.jsonl"),
                  oof_per_residue_sha256=sha256(output_dir / "oof_per_residue.csv.gz"),
                  analysis_script_sha256=sha256(Path(__file__)),
                  server_reproduction_passed=summary.get("server_reproduction_passed"),
                  baseline="Residue-pooled AA-type means fitted on four training folds only.")
    (output_dir / "analysis_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()
    run(args.input_dir.resolve(), args.output_dir.resolve(), args.allow_partial)


if __name__ == "__main__":
    main()
