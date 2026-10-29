"""Supplementary figure S1: where the unexplained RMSF variance comes from.

Four diagnostics on the frozen out-of-fold predictions:
a) per-system R2 distribution
b) concentration of squared error (top 1/5/10 percent of residues)
c) tail shrinkage: mean prediction versus mean target per target-RMSF bin
d) pooled / within-system / oracle-affine R2 summary

The oracle-affine value is a diagnostic upper bound (per-system linear
recalibration fitted on the test system itself) and is never a model result.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OOF = ROOT / "artifacts/evaluation_v2/analysis/oof_per_residue.csv.gz"
OUTPUT = ROOT / "output/flexppi_paper_figs_v3"


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(OOF, usecols=["sample_id", "dataset", "valid_rmsf",
                                      "target_rmsf", "pred_rmsf"])
    frame = frame[frame.valid_rmsf].reset_index(drop=True)
    target = frame.target_rmsf.to_numpy(float)
    prediction = frame.pred_rmsf.to_numpy(float)
    error = (target - prediction) ** 2
    pooled_r2 = 1 - error.sum() / ((target - target.mean()) ** 2).sum()

    grouped = frame.assign(error=error).groupby("sample_id")
    sse = grouped.error.sum()
    sst = grouped.target_rmsf.var(ddof=0) * grouped.size()
    per_system = (1 - sse / sst).replace([np.inf, -np.inf], np.nan).dropna()

    from scipy.stats import spearmanr
    per_system_rho = []
    for _, group in frame.groupby("sample_id"):
        if len(group) >= 4 and np.ptp(group.target_rmsf.to_numpy(float)) > 0:
            value = spearmanr(group.pred_rmsf, group.target_rmsf).statistic
            if np.isfinite(value):
                per_system_rho.append(value)
    per_system_rho = np.asarray(per_system_rho, dtype=float)

    centered_target = target - frame.groupby("sample_id").target_rmsf.transform("mean")
    centered_prediction = prediction - frame.groupby("sample_id").pred_rmsf.transform("mean")
    within_r2 = 1 - ((centered_target - centered_prediction) ** 2).sum() / (
        centered_target ** 2).sum()

    affine = []
    for _, group in frame.groupby("sample_id"):
        if len(group) < 3:
            continue
        slope, intercept = np.polyfit(group.pred_rmsf, group.target_rmsf, 1)
        fitted = slope * group.pred_rmsf + intercept
        total = ((group.target_rmsf - group.target_rmsf.mean()) ** 2).sum()
        if total > 0:
            affine.append(1 - ((group.target_rmsf - fitted) ** 2).sum() / total)
    affine = np.asarray(affine, dtype=float)

    order = np.sort(error)[::-1]
    cumulative = np.cumsum(order) / order.sum()
    fraction = np.arange(1, len(order) + 1) / len(order)

    edges = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, np.inf]
    rows = []
    for lower, upper in zip(edges[:-1], edges[1:]):
        mask = (target >= lower) & (target < upper)
        if not mask.any():
            continue
        rows.append({"bin": "{:g}-{:g}".format(lower, upper), "residues": int(mask.sum()),
                     "target_mean": float(target[mask].mean()),
                     "prediction_mean": float(prediction[mask].mean()),
                     "variance_ratio": float(prediction[mask].var()
                                             / max(target[mask].var(), 1e-12))})
    bins = pd.DataFrame(rows)

    figure, axes = plt.subplots(2, 3, figsize=(15.5, 7.4))
    axis_a, axis_b, axis_c, axis_d, axis_e, axis_f = axes.ravel()
    axis_f.axis("off")

    axis_a.hist(per_system_rho, bins=40, range=(0, 1), color="#D55E00")
    axis_a.axvline(np.median(per_system_rho), color="black", lw=1.4,
                   label="median {:.2f}".format(np.median(per_system_rho)))
    axis_a.set_xlabel("per-system Spearman")
    axis_a.set_ylabel("systems")
    axis_a.legend(fontsize=8)
    axis_a.text(0.03, 0.95,
                "negative: {} / {} systems ({:.2f}%)\n"
                "ranking is positive in virtually every system".format(
                    int((per_system_rho < 0).sum()), len(per_system_rho),
                    100 * (per_system_rho < 0).mean()),
                transform=axis_a.transAxes, fontsize=7.5, va="top")
    axis_a.set_title("a  system-level ranking (n={})".format(len(per_system_rho)),
                     fontsize=9, loc="left")

    axis_b.hist(per_system.to_numpy(float), bins=40, range=(-2, 1), color="#8A8F98")
    axis_b.axvline(per_system.median(), color="#D55E00", lw=1.4,
                   label="median {:.2f}".format(per_system.median()))
    axis_b.axvline(0, color="grey", lw=0.8)
    axis_b.set_xlabel("per-system R2")
    axis_b.set_ylabel("systems")
    axis_b.legend(fontsize=8)
    below = int((per_system < -2).sum())
    axis_b.text(0.03, 0.95,
                "{} below 0; {} below -2 (clipped)\nout-of-sample R2 turns negative when\n"
                "amplitudes are compressed; it is not anti-correlation".format(
                    int((per_system < 0).sum()), below),
                transform=axis_b.transAxes, fontsize=7, va="top")
    axis_b.set_title("b  system-level squared-error R2", fontsize=9, loc="left")

    axis_c.plot(fraction, cumulative, color="#0072B2", lw=1.6)
    for percent in (0.01, 0.05, 0.10):
        share = cumulative[int(len(cumulative) * percent) - 1]
        axis_c.scatter([percent], [share], color="#D55E00", s=20, zorder=4)
        axis_c.annotate("top {:.0f}% -> {:.0f}%".format(100 * percent, 100 * share),
                        (percent, share), textcoords="offset points", xytext=(6, -10),
                        fontsize=8)
    axis_c.set_xscale("log")
    axis_c.set_xlim(1e-3, 1)
    axis_c.set_ylim(0, 1.02)
    axis_c.set_xlabel("fraction of residues ranked by squared error")
    axis_c.set_ylabel("cumulative fraction of total SSE")
    axis_c.set_title("c  error concentration", fontsize=9, loc="left")

    axis_d.plot([0, 11], [0, 11], color="grey", ls="--", lw=0.9)
    axis_d.plot(bins.target_mean, bins.prediction_mean, "o-", color="#D55E00", lw=1.4)
    for index, (_, row) in enumerate(bins.iterrows()):
        offset = (5, -12) if index % 2 == 0 else (5, 5)
        axis_d.annotate("{:.1f}%".format(100 * row.residues / len(frame)),
                        (row.target_mean, row.prediction_mean),
                        textcoords="offset points", xytext=offset, fontsize=7)
    axis_d.set_xlabel("mean target RMSF in bin (A)")
    axis_d.set_ylabel("mean prediction in bin (A)")
    axis_d.set_xlim(0, 11)
    axis_d.set_ylim(0, 9)
    axis_d.set_title("d  tail shrinkage (labels are residue share)", fontsize=9, loc="left")

    values = [pooled_r2, within_r2, np.nanmean(affine)]
    labels = ["pooled", "within-system\n(system means removed)", "oracle per-system\naffine (diagnostic)"]
    axis_e.bar(np.arange(3), values, color=["#8A8F98", "#0072B2", "#009E73"], alpha=0.85)
    for position, value in enumerate(values):
        axis_e.text(position, value + 0.02, "{:.2f}".format(value), ha="center", fontsize=9)
    axis_e.set_xticks(np.arange(3))
    axis_e.set_xticklabels(labels, fontsize=8)
    axis_e.set_ylim(0, 0.85)
    axis_e.set_ylabel("R2 (959k residues)")
    axis_e.set_title("e  ranking versus scale", fontsize=9, loc="left")

    figure.suptitle("Where the unexplained RMSF variance sits: heavy tails, per-system scale, "
                    "and label noise", fontsize=10)
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    for suffix in ("png", "pdf"):
        figure.savefig(OUTPUT / "figS1_rmsf_error_diagnostics.{}".format(suffix), dpi=300)
    plt.close(figure)
    summary = {"pooled_r2": float(pooled_r2), "within_system_r2": float(within_r2),
               "oracle_affine_r2_mean": float(np.nanmean(affine)),
               "per_system_r2": {"median": float(per_system.median()),
                                 "mean": float(per_system.mean()),
                                 "negative_fraction": float((per_system < 0).mean()),
                                 "n": int(len(per_system))},
               "per_system_spearman": {"median": float(np.median(per_system_rho)),
                                       "negative_fraction": float((per_system_rho < 0).mean()),
                                       "n": int(len(per_system_rho))},
               "error_concentration": {"top1pct": float(cumulative[int(len(order) * 0.01) - 1]),
                                       "top5pct": float(cumulative[int(len(order) * 0.05) - 1]),
                                       "top10pct": float(cumulative[int(len(order) * 0.10) - 1])},
               "bins": bins.to_dict(orient="records"),
               "definition": "squared-error R2 on pooled out-of-fold residues; oracle affine "
                             "fits one slope/intercept per test system and is a diagnostic only"}
    (OUTPUT / "figS1_rmsf_diagnostics.json").write_text(json.dumps(summary, indent=2))
    print("wrote", OUTPUT / "figS1_rmsf_error_diagnostics.png")


if __name__ == "__main__":
    main()
