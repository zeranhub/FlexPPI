"""Paper figure 2 (merged RMSF + Hchi performance overview), draft version.

Panels
------
a) RMSF equal-system macro Spearman per dataset (points are the five outer folds).
b) RMSF residue-pooled MAE per dataset in Angstrom.
c) Common-cohort baseline comparison on the 520-system historical bridge.
d) Hchi pooled Spearman per dataset versus the amino-acid-corrected residual.
e) Top-20% high-flexibility retrieval (macro AUROC / average precision).

Only frozen files are read; nothing is retrained or re-selected.
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
EVALUATION = ROOT / "artifacts/performance_first_v2_atlas_hchi_v1/evaluation"
BRIDGE = ROOT / "artifacts/evaluation_v2/analysis/bridge_macro_metrics.csv"
OUTPUT = ROOT / "output/flexppi_paper_figs_v3"
DATASETS = ("atlas", "dynarepo", "itrace")
DATASET_LABELS = {"atlas": "ATLAS\n(single chain)", "dynarepo": "DynaRepo\n(mixed)",
                  "itrace": "ITrace\n(TCR-pMHC)"}
COLORS = {"flexppi": "#D55E00", "baseline": "#8A8F98", "raw": "#0072B2",
          "corrected": "#009E73", "retrieval": "#CC79A7"}


def load_fold_metrics() -> pd.DataFrame:
    rows = []
    for fold in range(5):
        data = json.loads((EVALUATION / "fold_{}_test.json".format(fold)).read_text())
        for metric in ("rmsf", "chi_entropy"):
            block = data["metrics"][metric]
            for dataset in DATASETS:
                node = block["by_dataset"][dataset]
                rows.append({"fold": fold, "metric": metric, "dataset": dataset,
                             "pooled_spearman": node["spearman"],
                             "pooled_mae": node["mae"],
                             "macro_spearman": node["macro_by_system"]["spearman"],
                             "system_count": node["system_count"]})
            rows.append({"fold": fold, "metric": metric, "dataset": "overall",
                         "pooled_spearman": block["overall"]["spearman"],
                         "pooled_mae": block["overall"]["mae"],
                         "macro_spearman": block["overall"]["macro_by_system"]["spearman"],
                         "system_count": block["overall"]["system_count"]})
            retrieval = block["high_flexibility_retrieval"]
            rows.append({"fold": fold, "metric": metric, "dataset": "retrieval",
                         "macro_auroc": retrieval["macro_auroc"],
                         "macro_average_precision": retrieval["macro_average_precision"],
                         "random_precision": retrieval["random_expected_precision"]})
    return pd.DataFrame(rows)


def load_corrected() -> dict:
    summary = json.loads((EVALUATION / "summary_five_fold.json").read_text())
    corrected = summary["metrics"]["chi_entropy"]["amino_acid_corrected"]
    return {"residual_spearman": corrected["residual_overall"]["spearman"]["mean"],
            "raw_spearman": summary["metrics"]["chi_entropy"]["overall"]["fold_mean"]
            ["spearman"]["mean"]}


def load_hchi_control() -> pd.DataFrame:
    """Per-fold Hchi type-identity control computed from the frozen OOF table.

    raw_model      model versus target
    raw_type       training-fold type-mean lookup versus target
    residual_model (prediction - type prior) versus (target - type prior)
    residual_type  oracle lookup (true test-fold type mean) minus prior, same residual target
    """
    from scipy.stats import spearmanr

    oof = pd.read_csv(ROOT / "artifacts/evaluation_v2/analysis/oof_per_residue.csv.gz",
                      usecols=["fold", "amino_acid", "valid_chi_entropy",
                               "target_chi_entropy", "pred_chi_entropy", "chi_type_prior"])
    oof = oof[oof.valid_chi_entropy].copy()
    rows = []
    for fold, group in oof.groupby("fold"):
        target = group["target_chi_entropy"].to_numpy(float)
        prior = group["chi_type_prior"].to_numpy(float)
        prediction = group["pred_chi_entropy"].to_numpy(float)
        oracle = (group.groupby("amino_acid")["target_chi_entropy"].transform("mean")
                  .to_numpy(float) - prior)
        rows.append({
            "fold": int(fold),
            "raw_model": spearmanr(prediction, target).statistic,
            "raw_type": spearmanr(prior, target).statistic,
            "residual_model": spearmanr(prediction - prior, target - prior).statistic,
            "residual_type": spearmanr(oracle, target - prior).statistic})
    return pd.DataFrame(rows)


def _ridge_cross_fitted_r2(design: np.ndarray, target: np.ndarray, fold: np.ndarray) -> float:
    total = float(((target - target.mean()) ** 2).sum())
    prediction = np.zeros(len(target), dtype=np.float32)
    for k in range(5):
        train = fold != k
        test = ~train
        x_train = design[train].astype(np.float64)
        gram = x_train.T @ x_train
        gram += np.eye(gram.shape[0]) * (1e-4 * np.trace(gram) / gram.shape[0])
        weights = np.linalg.solve(gram, x_train.T @ target[train].astype(np.float64))
        prediction[test] = (design[test].astype(np.float64) @ weights).astype(np.float32)
    return 1.0 - float(((target - prediction) ** 2).sum()) / total


def source_decomposition(metric: str) -> dict:
    """Cross-fitted R2 of nested knowledge sources for RMSF or Hchi.

    Order of addition: amino-acid type, position/termini, sequence neighbours
    (+-3), then the model's remaining gain (structural context and other).
    """
    columns = {"rmsf": ("valid_rmsf", "target_rmsf", "pred_rmsf"),
               "chi": ("valid_chi_entropy", "target_chi_entropy", "pred_chi_entropy")}[metric]
    oof = pd.read_csv(ROOT / "artifacts/evaluation_v2/analysis/oof_per_residue.csv.gz",
                      usecols=["fold", "sample_id", "chain", "residue_index", "amino_acid",
                               columns[0], columns[1], columns[2]])
    oof = oof[oof[columns[0]]].reset_index(drop=True)
    amino_acids = sorted(oof["amino_acid"].unique())
    code = {name: index for index, name in enumerate(amino_acids)}
    onehot = np.zeros((len(oof), len(amino_acids)), dtype=np.float32)
    onehot[np.arange(len(oof)), [code[name] for name in oof["amino_acid"]]] = 1.0
    neighbour_blocks = []
    for shift in (-3, -2, -1, 1, 2, 3):
        block = np.zeros_like(onehot)
        shifted = oof.groupby(["sample_id", "chain"], sort=False)["amino_acid"].shift(shift)
        present = shifted.notna().to_numpy()
        block[np.arange(len(oof)), shifted.fillna(oof["amino_acid"]).map(code).to_numpy()] = (
            present.astype(np.float32))
        neighbour_blocks.append(block)
        neighbour_blocks.append(np.where(present[:, None], 0.0, 1.0).astype(np.float32))
    length = oof.groupby(["sample_id", "chain"])["residue_index"].transform("max")
    relative = oof["residue_index"] / length.clip(lower=1)
    position = np.column_stack([
        relative.to_numpy(np.float32), (relative ** 2).to_numpy(np.float32),
        (oof["residue_index"] <= 1).to_numpy(np.float32),
        (oof["residue_index"] >= length - 1).to_numpy(np.float32),
        np.log1p(length.to_numpy(np.float32))])
    target = oof[columns[1]].to_numpy(np.float32)
    fold = oof["fold"].to_numpy()
    type_r2 = _ridge_cross_fitted_r2(onehot, target, fold)
    position_r2 = _ridge_cross_fitted_r2(np.column_stack([onehot, position]), target, fold)
    sequence_r2 = _ridge_cross_fitted_r2(
        np.column_stack([onehot, position] + neighbour_blocks), target, fold)
    model_r2 = 1.0 - float(((target - oof[columns[2]].to_numpy(np.float32)) ** 2).sum()) / float(
        ((target - target.mean()) ** 2).sum())
    return {"n": int(len(target)), "model_r2": model_r2, "type": type_r2,
            "position": position_r2 - type_r2, "sequence": sequence_r2 - position_r2,
            "structural": model_r2 - sequence_r2, "unexplained": 1.0 - model_r2}


def dot_bar(axis, positions, values, color, label=None, width=0.6):
    means = [np.mean(group) for group in values]
    axis.bar(positions, means, width=width, color=color, alpha=0.55, label=label)
    for position, group in zip(positions, values):
        axis.scatter(np.full(len(group), position) + np.linspace(-0.12, 0.12, len(group)),
                     group, s=16, color=color, zorder=3)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    metrics = load_fold_metrics()
    bridge = pd.read_csv(BRIDGE)

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.4))
    axis_a, axis_b, axis_c, axis_d, axis_e = axes.ravel()[:5]
    axes.ravel()[5].axis("off")

    order = ["atlas", "dynarepo", "itrace", "overall"]
    labels = [DATASET_LABELS.get(name, "overall") for name in order]
    positions = np.arange(len(order))

    subset = metrics[metrics.metric == "rmsf"]
    groups = [subset[subset.dataset == name]["macro_spearman"].to_numpy() for name in order]
    dot_bar(axis_a, positions, groups, COLORS["flexppi"])
    axis_a.set_xticks(positions)
    axis_a.set_xticklabels(labels, fontsize=8)
    axis_a.set_ylabel("RMSF Spearman\n(equal-system macro)")
    axis_a.set_ylim(0, 1)
    axis_a.axhline(0, color="grey", lw=0.6)
    axis_a.set_title("a  ranking within proteins", fontsize=9, loc="left")
    for position, name in zip(positions, order):
        count = int(subset[subset.dataset == name]["system_count"].iloc[0])
        axis_a.text(position, 0.02, "n={}".format(count), ha="center", fontsize=7)

    groups = [subset[subset.dataset == name]["pooled_mae"].to_numpy() for name in order]
    dot_bar(axis_b, positions, groups, COLORS["raw"])
    axis_b.set_xticks(positions)
    axis_b.set_xticklabels(labels, fontsize=8)
    axis_b.set_ylabel("RMSF MAE (A, pooled)")
    axis_b.set_title("b  absolute amplitude", fontsize=9, loc="left")

    cohort_labels = {"atlas_single_chain": "ATLAS single", "dynarepo_complex": "DynaRepo complex",
                     "itrace_oof_bridge": "ITrace"}
    cohorts = ["atlas_single_chain", "dynarepo_complex", "itrace_oof_bridge"]
    methods = ["FlexPPI-v2", "FastProtFlex", "BackFlip-seq", "BackFlip-structure",
               "ProDy-ANM", "RMSF-net-PDB"]
    width = 0.13
    for index, method in enumerate(methods):
        values = []
        for cohort in cohorts:
            row = bridge[(bridge.cohort == cohort) & (bridge.method == method)]
            values.append(float(row["spearman"].iloc[0]) if len(row) else np.nan)
        offset = (index - (len(methods) - 1) / 2) * width
        axis_c.bar(np.arange(len(cohorts)) + offset, values, width=width,
                   color=COLORS["flexppi"] if method == "FlexPPI-v2" else COLORS["baseline"],
                   alpha=0.9 if method == "FlexPPI-v2" else 0.6, label=method)
    axis_c.set_xticks(np.arange(len(cohorts)))
    axis_c.set_xticklabels([cohort_labels[name] for name in cohorts], fontsize=8)
    axis_c.set_ylabel("RMSF Spearman")
    axis_c.set_ylim(0, 1)
    axis_c.legend(fontsize=6.5, ncol=2, loc="lower right")
    axis_c.set_title("c  common-cohort baselines (520 systems)", fontsize=9, loc="left")

    decomposition = {name: source_decomposition(name) for name in ("rmsf", "chi")}
    control = load_hchi_control()
    segments = [("type", "amino-acid type", COLORS["baseline"]),
                ("position", "position / termini", "#56B4E9"),
                ("sequence", "sequence neighbours +/-3", "#E69F00"),
                ("structural", "structural context (model gain)", COLORS["corrected"]),
                ("unexplained", "unexplained", "#D9D9D9")]
    for index, metric in enumerate(("rmsf", "chi")):
        block = decomposition[metric]
        bottom = 0.0
        for key, label, color in segments:
            value = block[key]
            axis_d.bar([index], [value], bottom=bottom, width=0.55, color=color, alpha=0.9,
                       label=label if index == 0 else None)
            if value > 0.02:
                axis_d.text(index, bottom + value / 2, "{:.1f}".format(100 * value),
                            ha="center", va="center", fontsize=8)
            bottom += value
        axis_d.text(index, 1.02, "model R2={:.2f}".format(block["model_r2"]),
                    ha="center", fontsize=8)
    axis_d.set_xticks([0, 1])
    axis_d.set_xticklabels(["RMSF", "Hchi"], fontsize=9)
    axis_d.set_xlim(-0.6, 1.6)
    axis_d.set_ylim(0, 1.1)
    axis_d.set_ylabel("fraction of label variance (pooled OOF)")
    axis_d.set_title("d  what the correlations are made of", fontsize=9, loc="left")
    axis_d.text(1, 0.30,
                "type-only look-up on the\nHchi residual: rho upper bound {:.2f}\n(model {:.2f})".format(
                    control["residual_type"].mean(), control["residual_model"].mean()),
                ha="center", fontsize=7,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))
    handles, labels_ = axis_d.get_legend_handles_labels()
    axes.ravel()[5].legend(handles, labels_, loc="center", fontsize=8, frameon=False,
                           title="variance sources", title_fontsize=8)

    retrieval = metrics[metrics.dataset == "retrieval"]
    for metric, label, color in (("rmsf", "RMSF", COLORS["flexppi"]),
                                 ("chi_entropy", "Hchi", COLORS["retrieval"])):
        group = retrieval[retrieval.metric == metric]
        axis_e.bar([0 + (0.0 if metric == "rmsf" else 0.35)],
                   [group["macro_auroc"].mean()], width=0.3, color=color, alpha=0.75,
                   label="{} AUROC".format(label))
        axis_e.scatter(np.full(5, 0 + (0.0 if metric == "rmsf" else 0.35)) +
                       np.linspace(-0.08, 0.08, 5), group["macro_auroc"], s=12, color=color)
        axis_e.bar([0.7 + (0.0 if metric == "rmsf" else 0.35)],
                   [group["macro_average_precision"].mean()], width=0.3, color=color,
                   alpha=0.4, label="{} AP".format(label))
        axis_e.scatter(np.full(5, 0.7 + (0.0 if metric == "rmsf" else 0.35)) +
                       np.linspace(-0.08, 0.08, 5), group["macro_average_precision"],
                       s=12, color=color)
    axis_e.axhline(0.5, color="grey", ls="--", lw=0.8)
    axis_e.axhline(0.2, color="grey", ls=":", lw=0.8)
    axis_e.set_xticks([0, 0.35, 0.7, 1.05])
    axis_e.set_xticklabels(["AUROC", "", "AP", ""], fontsize=8)
    axis_e.set_ylabel("top-20% retrieval")
    axis_e.set_ylim(0, 1)
    axis_e.legend(fontsize=6.5, loc="lower right")
    axis_e.set_title("e  high-flexibility retrieval", fontsize=9, loc="left")

    fig.suptitle("FlexPPI absolute-state performance on the frozen five-fold outer tests "
                 "(epoch 40, untuned)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for suffix in ("png", "pdf"):
        fig.savefig(OUTPUT / "fig2_rmsf_hchi_performance.{}".format(suffix), dpi=300)
    plt.close(fig)
    print("wrote", OUTPUT / "fig2_rmsf_hchi_performance.png")


if __name__ == "__main__":
    main()
