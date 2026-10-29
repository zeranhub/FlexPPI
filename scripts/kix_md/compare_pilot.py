"""Descriptive MD-vs-FlexPPI comparisons of completed first repeats only.

No model selection, significance tests, confidence intervals or fitting.  Called
after each state. Cross-state contrasts appear only when both states are ready.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

STATES = ("apo", "myb", "mll", "both")
METRICS = {"rmsf": ("rmsf_angstrom", "predicted_rmsf_angstrom"),
           "hchi": ("chi_entropy_nats", "predicted_hchi_nats")}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        return
    temporary = path.with_suffix(".csv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def description(x, y):
    x, y = np.asarray(x), np.asarray(y)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < 2:
        raise ValueError("Insufficient matched residues")
    variable = np.ptp(x) > 1e-10 and np.ptp(y) > 1e-10
    return dict(n_residues=len(x), pearson_r=float(np.corrcoef(x, y)[0, 1]) if variable else None,
                spearman_r=float(spearmanr(x, y).statistic) if variable else None,
                mae=float(np.mean(np.abs(y - x))), rmse=float(np.sqrt(np.mean((y - x) ** 2))),
                mean_md=float(x.mean()), mean_model=float(y.mean()))


def compare(root):
    frozen = root / "frozen_predictions"
    meta = json.loads((frozen / "FROZEN.json").read_text())
    assert sha(frozen / "residue_predictions.csv") == meta["prediction_csv_sha256"]
    for state in STATES:
        assert sha(root / "inputs" / f"{state}.pdb") == meta["input_sha256"][state]
    predictions = {}
    for row in read_csv(frozen / "residue_predictions.csv"):
        if row["chain"] == "B":
            key = (row["state"], int(row["author_residue_id"]), int(row["fold_zero_based"]))
            assert key not in predictions
            predictions[key] = row
    assert len(predictions) == 4 * 87 * 5
    matched, summaries, regions, available, provenance, indexed = [], [], [], [], {}, {}
    for state in STATES:
        folder = root / "runs" / state / "rep01"
        path = folder / "analysis" / "summary.json"
        if not (folder / "DONE.json").exists() or not path.exists():
            continue
        info = json.loads(path.read_text())
        done = json.loads((folder / "DONE.json").read_text())
        assert done["repeat"] == 1 and done["production_ns"] == 100 and done["final_step"] == 55000000
        assert info["n_frames"] == 100000 and info["KIX_CA_count"] == 87
        assert info["state"] == state and info["replicate"] == "rep01"
        assert sha(root / "prepared" / state / "protein.pdb") == info["protein_pdb_sha256"]
        table = folder / "analysis" / "residue_metrics.csv"
        assert sha(table) == info["output_csv_sha256"]
        md_rows = [r for r in read_csv(table) if r["chain"] == "B"]
        assert [int(r["author_residue_number"]) for r in md_rows] == list(range(586, 673))
        available.append(state)
        provenance[state] = dict(analysis_summary_sha256=sha(path), metrics_sha256=sha(table),
                                 split_half_sensitivity=info["split_half_sensitivity"])
        for md in md_rows:
            residue = int(md["author_residue_number"])
            for metric, (md_field, model_field) in METRICS.items():
                if metric == "hchi" and md["amino_acid"] in ("A", "G", "P"):
                    continue
                pred = [predictions[(state, residue, fold)] for fold in range(5)]
                assert all(p["amino_acid"] == md["amino_acid"] for p in pred)
                values = [float(p[model_field]) for p in pred]
                for label, value in [(str(f), values[f]) for f in range(5)] + [("ensemble_mean", float(np.mean(values)))]:
                    row = dict(state=state, metric=metric, author_residue_number=residue,
                               amino_acid=md["amino_acid"], model_fold=label,
                               md_value=float(md[md_field]), model_value=value,
                               md_independent_repeats=1, uncertainty="not estimated from one MD repeat")
                    matched.append(row)
                    indexed[(state, metric, residue, label)] = row
    for state in available:
        for metric in METRICS:
            for fold in ["0", "1", "2", "3", "4", "ensemble_mean"]:
                rows = [r for r in matched if r["state"] == state and r["metric"] == metric and r["model_fold"] == fold]
                summaries.append(dict(state=state, metric=metric, model_fold=fold,
                                      **description([r["md_value"] for r in rows], [r["model_value"] for r in rows])))
                for name, low, high in [("L12_614_621", 614, 621), ("KIX_586_672", 586, 672)]:
                    selection = [r for r in rows if low <= r["author_residue_number"] <= high]
                    x, y = np.array([r["md_value"] for r in selection]), np.array([r["model_value"] for r in selection])
                    regions.append(dict(state=state, metric=metric, model_fold=fold, region=name,
                                        n=len(selection), md_mean=float(x.mean()), model_mean=float(y.mean()),
                                        md_sqrt_mean_squared=float(np.sqrt(np.mean(x*x))) if metric == "rmsf" else "",
                                        model_sqrt_mean_squared=float(np.sqrt(np.mean(y*y))) if metric == "rmsf" else ""))
    contrast_rows, contrast_stats = [], []
    for before, after in [("apo", "myb"), ("apo", "mll"), ("apo", "both"), ("myb", "mll")]:
        if before not in available or after not in available:
            continue
        for row in matched:
            if row["state"] != after:
                continue
            key = (before, row["metric"], row["author_residue_number"], row["model_fold"])
            previous = indexed[key]
            contrast_rows.append(dict(contrast=f"{after}_minus_{before}", metric=row["metric"],
                                     author_residue_number=row["author_residue_number"], amino_acid=row["amino_acid"],
                                     model_fold=row["model_fold"], md_delta=row["md_value"]-previous["md_value"],
                                     model_delta=row["model_value"]-previous["model_value"], md_repeats_per_state=1))
        for metric in METRICS:
            for fold in ["0", "1", "2", "3", "4", "ensemble_mean"]:
                sel = [r for r in contrast_rows if r["contrast"] == f"{after}_minus_{before}" and r["metric"] == metric and r["model_fold"] == fold]
                for name, low, high in [("L12_614_621", 614, 621), ("KIX_586_672", 586, 672)]:
                    rs = [r for r in sel if low <= r["author_residue_number"] <= high]
                    contrast_stats.append(dict(contrast=f"{after}_minus_{before}", metric=metric, model_fold=fold,
                                               region=name, **description([r["md_delta"] for r in rs], [r["model_delta"] for r in rs])))
    output = root / "pilot_comparison"
    output.mkdir(exist_ok=True)
    for filename, rows in [("matched_residues.csv", matched), ("absolute_metrics.csv", summaries),
                           ("regions.csv", regions), ("residue_deltas.csv", contrast_rows),
                           ("delta_metrics.csv", contrast_stats)]:
        write_csv(output / filename, rows)
    report = dict(available_states=available, complete_four_state_pilot=len(available)==4,
                  independent_md_repeats_per_state=1, input_provenance=provenance,
                  model_predictions_frozen_sha256=meta["prediction_csv_sha256"],
                  limitations=["Exploratory single-repeat comparison, not a replication or convergence proof.",
                               "MD is a simulation reference, not experimental truth.",
                               "No confidence intervals, p-values or replicate SD are estimated.",
                               "Adjacent residues are not independent experimental repeats.",
                               "Five model checkpoints are not five MD repeats; all are retained without selection.",
                               "Near-identical KIX training sequence prevents interpreting five-fold average as independent validation."],
                  next_step="Review stability, split-half sensitivity, absolute/relative agreement including failures; additional repeats require user approval.")
    temporary = output / "summary.json.tmp"
    temporary.write_text(json.dumps(report, indent=2) + "\n")
    temporary.replace(output / "summary.json")
    print(json.dumps(report), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    compare(parser.parse_args().root.resolve(strict=True))
