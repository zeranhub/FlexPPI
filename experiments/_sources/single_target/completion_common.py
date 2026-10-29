"""Metric-aware utilities for the frozen chemical_identity_corrected_v1 corpus."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

LABEL_NAMES = ("rmsf", "chi_entropy", "sasa")
ARM_METRICS = {"full": ["rmsf", "chi_entropy"], "rmsf_only": ["rmsf"], "chi_entropy_only": ["chi_entropy"]}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def jsonl(path):
    return [json.loads(line) for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def correlation(a, b, ranks=False):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if len(a) < 4 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    if ranks:
        a = pd.Series(a).rank(method="average").to_numpy()
        b = pd.Series(b).rank(method="average").to_numpy()
    return float(np.corrcoef(a, b)[0, 1])


def statistics(target, prediction):
    good = np.isfinite(target) & np.isfinite(prediction)
    target, prediction = np.asarray(target)[good], np.asarray(prediction)[good]
    if not len(target):
        return {key: float("nan") for key in ["mae", "rmse", "bias", "pearson", "spearman"]} | {"count": 0}
    error = prediction - target
    return dict(count=len(target), mae=float(np.mean(abs(error))), rmse=float(np.sqrt(np.mean(error ** 2))),
                bias=float(error.mean()), pearson=correlation(target, prediction), spearman=correlation(target, prediction, True))


def cache_path(evaluation, row):
    """Use the local evaluation root even when the copied manifest has old absolute paths."""
    root = Path(evaluation)
    candidate = root / row["variant"] / (row["sample_id"] + ".npz")
    if candidate.is_file():
        return candidate
    candidate = Path(row["cache_path"])
    if not candidate.is_file():
        raise FileNotFoundError(str(candidate))
    return candidate


def read_cache(evaluation, row, metric):
    with np.load(cache_path(evaluation, row), allow_pickle=False) as cache:
        output_names = cache["metric_names"].tolist() if "metric_names" in cache else ARM_METRICS[row["variant"]]
        label_names = cache["label_names"].tolist() if "label_names" in cache else list(LABEL_NAMES)
        output_column = output_names.index(metric)
        label_column = label_names.index(metric)
        prediction = cache["mean"][:, output_column].astype(float)
        model_variance = cache["model_variance"][:, output_column].astype(float)
        target = cache["target"][:, label_column].astype(float)
        label_variance = cache["label_variance"][:, label_column].astype(float)
        valid = cache["valid"][:, label_column].astype(bool)
        valid &= np.isfinite(prediction) & np.isfinite(target) & np.isfinite(model_variance) & (model_variance > 0)
        return dict(prediction=prediction, model_variance=model_variance, target=target,
                    label_variance=label_variance, valid=valid)


def fold_summary(frame, values, keys):
    """Record macro first; the same record contribution as the original corrected export."""
    by_source = frame.groupby(keys + ["fold", "dataset"], dropna=False)[values].mean().reset_index()
    overall = frame.groupby(keys + ["fold"], dropna=False)[values].mean().reset_index().assign(dataset="all")
    folds = pd.concat([by_source, overall], ignore_index=True)
    summary = folds.groupby(keys + ["dataset"], dropna=False)[values].agg(["mean", "std"])
    summary.columns = [name + "_" + scope for name, scope in summary.columns]
    return folds, summary.reset_index()
