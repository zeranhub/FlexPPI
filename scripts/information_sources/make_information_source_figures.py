"""Information-source ladder figures for FlexPPI (4.2 review, 2026-09-19).

The five-fold histogram-boosting probes in ``run_static_baselines.py`` are
compared with the frozen five-fold FlexPPI out-of-fold predictions. Every
comparison is paired on the same residue set, and every interval is resampled
over protein sequences rather than over residues or folds.

Scope limit that must travel with every number produced here: the probes and
FlexPPI differ in input information *and* in model family and capacity, so a
probe/FlexPPI gap is not a clean ablation of one architecture. These figures
size the extra value of the frozen model over fixed feature summaries; they do
not prove that one input channel causes the difference.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Academic figure typography baseline, kept verbatim from the project scripts.
import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
    # Academic figure export baseline.
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

import matplotlib.pyplot as plt  # noqa: E402  (needs the rcParams above)
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
SEQUENTIAL = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED = "#B2182B"
GREY = "#999999"
BLACK = "#222222"
SEQ_CMAP = LinearSegmentedColormap.from_list("seq", SEQUENTIAL)

ROOT = Path(__file__).resolve().parents[2]
METRICS = ["rmsf", "chi_entropy"]
UNITS = {"rmsf": "\u00c5", "chi_entropy": "nats"}
LABELS = {"rmsf": "RMSF", "chi_entropy": r"$H_\chi$"}
# (key, display name, description recorded in the audit file)
METHODS = [
    ("aa", "AA prior", "Training-fold amino-acid label mean; residue identity only."),
    ("sequence", "S", "Residue identity, plus/minus 4 sequence offsets, chain position, composition."),
    ("geometry", "G (no identity)", "Backbone angles, intra-chain and partner contact counts, chain position; no explicit residue identity."),
    ("sequence_intra", "S+g", "S plus local structure summaries: phi/psi, N-CA-C angle, intra-chain contact counts."),
    ("sequence_complex", "S+g+p", "S+g plus cross-chain partner contact counts."),
    ("flexppi", "FlexPPI", "Frozen five-fold dual-target model; full 3D geometry with invariant point attention."),
]
METHOD_KEYS = [key for key, _, _ in METHODS]
PROBE_KEYS = ["sequence", "geometry", "sequence_intra", "sequence_complex"]
# Retrained neural arms (same data, folds, seed and 40-epoch budget as the
# frozen model).  "flexppi" above is the frozen full-complex model, which also
# serves as the NN-full arm; these three remove one information source each.
ABLATION_METHODS = [
    ("nn_sequence", "NN-S", "Retrained sequence-only FlexPPI: identity and chain-local scalar attention; no coordinates, no partner."),
    ("nn_geometry", "NN-G (no identity)", "Retrained FlexPPI with residue identity replaced by one unknown token; backbone geometry and partner contacts retained."),
    ("nn_intrachain", "NN (no partner)", "Retrained FlexPPI evaluated one chain at a time; every cross-chain IPA message removed."),
]
ABLATION_KEYS = [key for key, _, _ in ABLATION_METHODS]
ABLATION_FILES = {"nn_sequence": "sequence", "nn_geometry": "geometry", "nn_intrachain": "intrachain"}
ALL_METHODS = METHODS + ABLATION_METHODS
ALL_METHOD_KEYS = [key for key, _, _ in ALL_METHODS]
ATTRITION_ORDER = [
    ("aa", "AA prior"), ("sequence", "S"), ("geometry", "G"),
    ("sequence_intra", "S+g"), ("sequence_complex", "S+g+p"),
    ("nn_sequence", "NN-S"), ("nn_geometry", "NN-G"),
    ("nn_intrachain", "NN no partner"), ("flexppi", "NN full"),
]
GAIN_SPECS = [
    ("identity: NN full \u2212 NN-G", "flexppi", "nn_geometry"),
    ("partner: NN full \u2212 NN no partner", "flexppi", "nn_intrachain"),
    ("geometry+partner: NN full \u2212 NN-S", "flexppi", "nn_sequence"),
    ("3D model vs probes: NN full \u2212 (S+g)", "flexppi", "sequence_intra"),
]
LADDER_STEPS = [
    ("sequence", "aa", "S \u2212 AA prior"),
    ("sequence_intra", "sequence", "(S+g) \u2212 S"),
    ("sequence_complex", "sequence_intra", "(S+g+p) \u2212 (S+g)"),
    ("flexppi", "sequence_intra", "FlexPPI \u2212 (S+g)"),
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    mapped = series.astype(str).map({"True": True, "False": False})
    if mapped.isna().any():
        raise ValueError("Non-boolean validity column")
    return mapped


def bootstrap_mean(values, seed=20260920, repeats=5000):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan, np.nan
    generator = np.random.default_rng(seed)
    means = np.empty(repeats)
    for start in range(0, repeats, 200):
        block = min(200, repeats - start)
        means[start:start + block] = values[generator.integers(0, len(values), size=(block, len(values)))].mean(axis=1)
    return float(values.mean()), float(np.quantile(means, .025)), float(np.quantile(means, .975))


def load_frame(probe_dir: Path, residues: Path, chain_map_path: Path, ablation_dir: Path = None):
    """Rebuild one row per residue carrying labels, masks, every prediction and the protein key."""
    schema = json.loads((probe_dir / "feature_schema.json").read_text(encoding="utf-8"))
    data = np.load(probe_dir / "features.npz")
    records = pd.DataFrame(schema["records"])
    record_index = data["record_index"].astype(np.int64)
    residue_index = data["residue_index"].astype(np.int64)
    if record_index.min() < 0 or record_index.max() >= len(records):
        raise ValueError("record_index out of range for the frozen schema")

    frame = pd.DataFrame({
        "sample_id": records["sample_id"].to_numpy()[record_index],
        "residue_index": residue_index,
        "fold": data["fold"].astype(np.int64),
        "target_rmsf": data["target"][:, 0],
        "target_chi_entropy": data["target"][:, 1],
        "valid_rmsf": data["valid"][:, 0],
        "valid_chi_entropy": data["valid"][:, 1],
        "flexppi_rmsf": data["full"][:, 0],
        "flexppi_chi_entropy": data["full"][:, 1],
    })
    for key in PROBE_KEYS:
        path = probe_dir / (key + "_oof.npz")
        if not path.exists():
            raise FileNotFoundError("Missing probe output: " + str(path))
        prediction = np.load(path)["prediction"]
        if prediction.shape != data["target"].shape:
            raise ValueError("Probe output shape mismatch: " + key)
        frame[key + "_rmsf"] = prediction[:, 0]
        frame[key + "_chi_entropy"] = prediction[:, 1]

    if ablation_dir is not None:
        for key, arm in ABLATION_FILES.items():
            path = ablation_dir / arm / "oof_per_residue.csv.gz"
            if not path.exists():
                raise FileNotFoundError("Missing ablation table: " + str(path))
            columns = ["sample_id", "residue_index", "target_rmsf", "target_chi_entropy",
                       "valid_rmsf", "valid_chi_entropy", "pred_rmsf", "pred_chi_entropy"]
            arm_frame = pd.read_csv(path, dtype={"sample_id": str}, usecols=columns)
            if arm_frame.duplicated(["sample_id", "residue_index"]).any():
                raise ValueError("Duplicated ablation rows: " + arm)
            numeric = [c for c in columns if c.startswith(("target_", "pred_"))]
            for column in numeric:
                arm_frame[column] = pd.to_numeric(arm_frame[column], errors="coerce")
            for metric in METRICS:
                arm_frame["valid_" + metric] = as_bool(arm_frame["valid_" + metric])
            rename = {c: "arm_%s_%s" % (key, c) for c in columns
                      if c not in ("sample_id", "residue_index")}
            arm_frame = arm_frame.rename(columns=rename)
            frame = frame.merge(arm_frame, on=["sample_id", "residue_index"], how="left",
                                validate="one_to_one")
            for metric in METRICS:
                target = frame["arm_%s_target_%s" % (key, metric)]
                valid = frame["arm_%s_valid_%s" % (key, metric)]
                if not np.allclose(frame["target_" + metric], target, rtol=0, atol=1e-6, equal_nan=True):
                    raise ValueError("Ablation %s and feature matrix disagree on %s labels" % (arm, metric))
                if not np.array_equal(frame["valid_" + metric].to_numpy(bool), valid.to_numpy(bool)):
                    raise ValueError("Ablation %s and feature matrix disagree on the %s mask" % (arm, metric))
                frame[key + "_" + metric] = frame["arm_%s_pred_%s" % (key, metric)]
                if frame.loc[frame["valid_" + metric], key + "_" + metric].isna().any():
                    raise ValueError("Missing %s predictions for %s" % (metric, key))
                for suffix in ("target_", "valid_", "pred_"):
                    del frame["arm_%s_%s%s" % (key, suffix, metric)]

    columns = ["sample_id", "residue_index", "chain",
               "aa_prior_rmsf", "aa_prior_chi_entropy",
               "target_rmsf", "target_chi_entropy",
               "pred_rmsf", "pred_chi_entropy",
               "valid_rmsf", "valid_chi_entropy"]
    oof = pd.read_csv(residues, dtype={"sample_id": str, "chain": str},
                      keep_default_na=False, usecols=columns)
    # keep_default_na=False keeps empty label cells as ""; the frozen table
    # stores no value for masked targets, so coerce only the numeric columns.
    numeric = [c for c in columns if c.startswith(("target_", "valid_", "pred_", "aa_"))]
    for column in numeric:
        if column.startswith("valid_"):
            continue
        oof[column] = pd.to_numeric(oof[column], errors="coerce")
    oof["valid_rmsf"] = as_bool(oof["valid_rmsf"])
    oof["valid_chi_entropy"] = as_bool(oof["valid_chi_entropy"])
    if oof.duplicated(["sample_id", "residue_index"]).any():
        raise ValueError("Duplicated (sample_id, residue_index) in the frozen residue table")
    # The frozen table repeats the labels and masks; rename them so the merge
    # cannot silently drop the mask columns the aggregation still needs.
    oof = oof.rename(columns={c: "oof_" + c for c in oof.columns if c.startswith(("target_", "valid_"))})
    merged = frame.merge(oof, on=["sample_id", "residue_index"], how="left", validate="one_to_one")
    if merged["chain"].isna().any():
        raise ValueError("Residue rows without a chain mapping")
    # The two sources must agree about the residue, not merely about the counts.
    for metric in METRICS:
        # Masked residues carry NaN on both sides, so NaN must compare equal here.
        if not np.allclose(merged["target_" + metric], merged["oof_target_" + metric],
                           rtol=0, atol=1e-6, equal_nan=True):
            raise ValueError("Frozen residue table and feature matrix disagree on the " + metric + " labels")
        if not np.array_equal(merged["valid_" + metric].to_numpy(bool), merged["oof_valid_" + metric].to_numpy(bool)):
            raise ValueError("Frozen residue table and feature matrix disagree on the " + metric + " mask")
    if not np.allclose(merged["flexppi_rmsf"], merged["pred_rmsf"], rtol=0, atol=1e-6):
        raise ValueError("Feature matrix 'full' column and frozen prediction differ")
    merged["flexppi_chi_entropy"] = merged["pred_chi_entropy"]
    merged["aa_rmsf"] = merged["aa_prior_rmsf"]
    merged["aa_chi_entropy"] = merged["aa_prior_chi_entropy"]
    merged["valid_rmsf"] = merged["valid_rmsf"].to_numpy(bool)
    merged["valid_chi_entropy"] = merged["valid_chi_entropy"].to_numpy(bool)
    merged = merged.drop(columns=[c for c in merged.columns if c.startswith("oof_") or c.startswith("pred_")
                                  or c.startswith("aa_prior_")])

    chain_map = pd.read_csv(chain_map_path, dtype={"sample_id": str, "chain": str},
                            usecols=["sample_id", "chain", "protein_id", "full_chain_length", "source"])
    merged = merged.merge(chain_map, on=["sample_id", "chain"], how="left", validate="many_to_one")
    if merged["protein_id"].isna().any():
        raise ValueError("Residue rows without a protein identity")
    return merged


def protein_error_table(frame: pd.DataFrame):
    """Chain mean absolute error, then equal weight over chain copies, records and proteins."""
    keys_used = [key for key in ALL_METHOD_KEYS if key + "_rmsf" in frame.columns]
    blocks, pooled, valid_counts = [], {}, {}
    for metric in METRICS:
        valid = frame["valid_" + metric].to_numpy(bool)
        target = frame["target_" + metric].to_numpy(float)[valid]
        keys = frame.loc[valid, ["sample_id", "chain", "protein_id", "full_chain_length", "source"]].copy()
        for key in keys_used:
            prediction = frame[key + "_" + metric].to_numpy(float)[valid]
            if not np.isfinite(prediction).all():
                raise ValueError("Non-finite %s prediction for %s" % (metric, key))
            keys[key] = np.abs(prediction - target)
            pooled.setdefault(metric, {})[key] = float(keys[key].mean())
        per_chain = keys.groupby(["protein_id", "sample_id", "chain"], sort=False)[keys_used].mean()
        per_record = per_chain.groupby(["protein_id", "sample_id"], sort=False)[keys_used].mean()
        per_protein = per_record.groupby("protein_id", sort=False)[keys_used].mean()
        meta = keys.groupby("protein_id", sort=False).agg(
            full_chain_length=("full_chain_length", "first"),
            source=("source", lambda values: "|".join(sorted(set(values)))),
            n_records=("sample_id", "nunique"))
        block = per_protein.join(meta).reset_index()
        block["metric"] = metric
        block["short_chain"] = block["full_chain_length"] < 20
        blocks.append(block)
        valid_counts[metric] = int(valid.sum())
    return pd.concat(blocks, ignore_index=True), pooled, valid_counts


def step_table(table: pd.DataFrame):
    rows = []
    for metric in METRICS:
        block = table.loc[table.metric == metric].set_index("protein_id")
        for higher, baseline, name in LADDER_STEPS:
            gain = (block[baseline] - block[higher]).to_numpy(float)
            finite = gain[np.isfinite(gain)]
            mean, low, high = bootstrap_mean(finite)
            rows.append(dict(metric=metric, step=name, higher=higher, baseline=baseline,
                             n_proteins=len(finite), mean_gain=mean, ci_low=low, ci_high=high,
                             median_gain=float(np.median(finite)),
                             improved_fraction=float((finite > 0).mean()), unit=UNITS[metric]))
    return pd.DataFrame(rows)


def save_figure(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def figure_paired(table: pd.DataFrame, out: Path) -> dict:
    # ``layout="constrained"`` reserves room for the two colorbars, which the
    # default packer did not: they used to sit on top of the right-hand labels.
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), layout="constrained")
    windows = {}
    for ax, metric in zip(axes, METRICS):
        block = table.loc[table.metric == metric].dropna(subset=["sequence", "flexppi"])
        x = block["sequence"].to_numpy(float)
        y = block["flexppi"].to_numpy(float)
        high = float(np.quantile(np.concatenate([x, y]), 0.995)) * 1.05
        handle = ax.scatter(x, y, c=np.log10(block["full_chain_length"].to_numpy(float)), cmap=SEQ_CMAP,
                            s=7, linewidths=0, alpha=0.85, vmin=1.0, vmax=3.2)
        ax.plot([0, high], [0, high], color=GREY, linewidth=0.8, linestyle="--", zorder=0)
        ax.set_xlim(0, high)
        ax.set_ylim(0, high)
        ax.set_xlabel("Sequence probe MAE (%s)" % UNITS[metric])
        ax.set_ylabel("FlexPPI MAE (%s)" % UNITS[metric])
        reduction = x - y
        mean, low, upper = bootstrap_mean(reduction)
        ax.set_title("%s: n=%d proteins, %.0f%% below diagonal"
                     % (LABELS[metric], len(block), 100 * (reduction > 0).mean()))
        ax.text(0.04, 0.95, "mean reduction %.3f (%.3f, %.3f)" % (mean, low, upper),
                transform=ax.transAxes, fontsize=7, va="top", color=BLACK)
        # Long chains reach a sequence-probe MAE far above the window, so the
        # hidden points have to be declared rather than silently dropped.
        outside = int(((x > high) | (y > high)).sum())
        if outside:
            ax.text(0.96, 0.03, "%d proteins above the axis (max %.1f %s); all stay in the statistics"
                                % (outside, float(max(x.max(), y.max())), UNITS[metric]),
                    transform=ax.transAxes, fontsize=6, ha="right", va="bottom", color=GREY)
        windows[metric] = dict(axis_limit=high, proteins_outside_axis=outside,
                               max_sequence_probe_mae=float(x.max()), max_flexppi_mae=float(y.max()))
        fig.colorbar(handle, ax=ax, fraction=0.046, pad=0.02).set_label("log10 chain length", fontsize=7)
    fig.text(0.01, -0.05,
             "One point per exact protein sequence; chain MAE is averaged over chain copies and over records.\n"
             "Sequence probe = fixed-budget boosting on identity, sequence offsets, chain position and composition.\n"
             "The probe has no 3D input, so this panel sizes the total value of the frozen model, not one channel.",
             fontsize=7, color=BLACK)
    save_figure(fig, out / "fig_information_ladder_paired")
    return windows


def figure_steps(table: pd.DataFrame, summary: pd.DataFrame, out: Path) -> None:
    # sharey removes the duplicated step labels on the right panel, which is
    # what used to push the two long x labels on top of each other.
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), sharey=True, layout="constrained")
    generator = np.random.default_rng(7)
    for ax, metric in zip(axes, METRICS):
        block = table.loc[table.metric == metric].set_index("protein_id")
        names = [name for _, _, name in LADDER_STEPS]
        values = [(block[baseline] - block[higher]).to_numpy(float) for higher, baseline, _ in LADDER_STEPS]
        positions = np.arange(len(values))[::-1]
        parts = ax.violinplot([value[np.isfinite(value)] for value in values], positions=positions,
                              vert=False, widths=0.85, showextrema=False, showmedians=False)
        for body in parts["bodies"]:
            body.set_facecolor(SEQUENTIAL[1])
            body.set_alpha(0.35)
            body.set_edgecolor(GREY)
            body.set_linewidth(0.5)
        for position, name, value in zip(positions, names, values):
            finite = value[np.isfinite(value)]
            ax.scatter(finite, position + generator.normal(0, 0.055, size=len(finite)),
                       s=3, color=GREY, alpha=0.35, linewidths=0)
            q1, median, q3 = np.quantile(finite, [.25, .5, .75])
            ax.plot([q1, q3], [position, position], color=BLACK, linewidth=1.1, solid_capstyle="butt")
            ax.plot([median], [position], marker="o", markersize=3.6, color=ACCENT_RED, zorder=4)
            row = summary[(summary.metric == metric) & (summary.step == name)]
            if len(row) == 1:
                ax.plot([row.ci_low.iloc[0], row.ci_high.iloc[0]], [position, position],
                        color=ACCENT_RED, linewidth=0.9)
        ax.axvline(0, color=BLACK, linewidth=0.7)
        ax.set_yticks(positions)
        ax.set_yticklabels([name.replace(" \u2212 ", " \u2212 ") for name in names])
        ax.set_xlabel("Protein-equal MAE reduction (%s)" % UNITS[metric])
        ax.set_title(LABELS[metric])
        ax.set_ylim(-0.6, len(values) - 0.4)
    fig.text(0.01, -0.05,
             "Positive means the higher information level has the lower protein-equal MAE; every protein keeps its weight in every step.\n"
             "Violin = protein density, grey points = individual proteins, bar = IQR, dot and line = mean with a bootstrap\n"
             "95% interval resampled over protein sequences. FlexPPI - (S+g) is not a single-channel ablation: the probes and\n"
             "the frozen model also differ in model family and capacity.",
             fontsize=7, color=BLACK)
    save_figure(fig, out / "fig_information_ladder_steps")


def figure_dual_gain(table: pd.DataFrame, out: Path) -> dict:
    wide = table.pivot_table(index="protein_id", columns="metric",
                             values=["sequence_intra", "flexppi"], aggfunc="first")
    ok = wide[("sequence_intra", "rmsf")].notna() & wide[("sequence_intra", "chi_entropy")].notna()
    wide = wide.loc[ok]
    gain_rmsf = (wide[("sequence_intra", "rmsf")] - wide[("flexppi", "rmsf")]).to_numpy(float)
    gain_hchi = (wide[("sequence_intra", "chi_entropy")] - wide[("flexppi", "chi_entropy")]).to_numpy(float)
    keep = np.isfinite(gain_rmsf) & np.isfinite(gain_hchi)
    gain_rmsf, gain_hchi = gain_rmsf[keep], gain_hchi[keep]
    # The two targets have unrelated units, so one shared symmetric window
    # either hides proteins or wastes most of the panel. Pad each axis to its
    # own full range instead, and say so in the caption.
    def pad(values):
        low, high = float(np.min(values)), float(np.max(values))
        margin = 0.06 * max(high - low, 1e-6)
        return low - margin, high + margin

    x_low, x_high = pad(gain_rmsf)
    y_low, y_high = pad(gain_hchi)
    fig, ax = plt.subplots(figsize=(3.6, 3.4))
    ax.scatter(gain_rmsf, gain_hchi, s=7, color=CATEGORICAL[0], alpha=0.6, linewidths=0)
    ax.axhline(0, color=BLACK, linewidth=0.7)
    ax.axvline(0, color=BLACK, linewidth=0.7)
    ax.set_xlim(x_low, x_high)
    ax.set_ylim(y_low, y_high)
    ax.set_xlabel("RMSF MAE reduction (\u00c5)")
    ax.set_ylabel(r"$H_\chi$ MAE reduction (nats)")
    quadrants = [int(((gain_rmsf > 0) & (gain_hchi > 0)).sum()), int(((gain_rmsf > 0) & (gain_hchi <= 0)).sum()),
                 int(((gain_rmsf <= 0) & (gain_hchi > 0)).sum()), int(((gain_rmsf <= 0) & (gain_hchi <= 0)).sum())]
    ax.set_title("FlexPPI \u2212 (S+g), per protein", fontsize=8)
    ax.text(0.0, 1.02, "n=%d | both %d / RMSF only %d / $H_\\chi$ only %d / neither %d"
                       % (len(gain_rmsf), quadrants[0], quadrants[1], quadrants[2], quadrants[3]),
            transform=ax.transAxes, fontsize=7, va="bottom", ha="left")
    fig.text(0.01, -0.02, "Axes are padded to the full range of the data and scaled independently, so no protein\n"
                          "is dropped; the quadrant counts cover all %d proteins." % len(gain_rmsf),
             fontsize=7, color=BLACK)
    save_figure(fig, out / "fig_information_ladder_dual_gain")
    return dict(proteins=int(len(gain_rmsf)),
                quadrant_counts=dict(both=quadrants[0], rmsf_only=quadrants[1],
                                     hchi_only=quadrants[2], neither=quadrants[3]),
                x_range=[x_low, x_high], y_range=[y_low, y_high])


def figure_attribution_ladder(table: pd.DataFrame, out: Path) -> None:
    """Protein-equal MAE for probes and retrained neural arms, per target."""
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4), layout="constrained")
    rows = []
    positions = np.arange(len(ATTRITION_ORDER))[::-1]
    for ax, metric in zip(axes, METRICS):
        block = table.loc[table.metric == metric]
        for position, (key, display) in zip(positions, ATTRITION_ORDER):
            finite = block[key].dropna().to_numpy(float)
            finite = finite[np.isfinite(finite)]
            mean, low, high = bootstrap_mean(finite)
            color = CATEGORICAL[1] if key in ABLATION_KEYS or key == "flexppi" else CATEGORICAL[0]
            ax.errorbar([mean], [position], xerr=[[mean - low], [high - mean]], fmt="o",
                        markersize=4.2, color=color, ecolor=color, elinewidth=1.0, capsize=2)
            rows.append(dict(metric=metric, method=key, display=display, n_proteins=len(finite),
                             mean_mae=mean, ci_low=low, ci_high=high, unit=UNITS[metric]))
        ax.set_yticks(positions)
        ax.set_title(LABELS[metric])
        ax.set_xlabel("Protein-equal MAE (%s)" % UNITS[metric])
    axes[1].set_yticks([])
    axes[0].set_yticklabels([display for _, display in ATTRITION_ORDER])
    fig.text(0.01, -0.13,
             "Blue = fixed-budget boosting probes and the frozen reference model; red = retrained neural arms (same folds, labels, masks, seed\n"
             "and 40-epoch budget; NN-S has fewer parameters). Dot = mean over proteins, line = bootstrap 95% interval resampled over exact\n"
             "protein sequences. Identity is hidden by one unknown token, so NN-G can still infer residue type from geometry and partners.",
             fontsize=7, color=BLACK)
    save_figure(fig, out / "fig_information_attribution_ladder")
    pd.DataFrame(rows).to_csv(out / "attribution_ladder.csv", index=False)


def figure_attribution_gains(table: pd.DataFrame, out: Path) -> None:
    """Paired per-protein gains between retrained arms, one panel per target."""
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.9), sharey=True, layout="constrained")
    rows = []
    positions = np.arange(len(GAIN_SPECS))[::-1]
    for ax, metric in zip(axes, METRICS):
        block = table.loc[table.metric == metric]
        for position, (name, higher, baseline) in zip(positions, GAIN_SPECS):
            gain = (block[baseline] - block[higher]).to_numpy(float)
            finite = gain[np.isfinite(gain)]
            mean, low, high = bootstrap_mean(finite)
            ax.errorbar([mean], [position], xerr=[[mean - low], [high - mean]], fmt="o",
                        markersize=4.2, color=ACCENT_RED, ecolor=ACCENT_RED, elinewidth=1.0, capsize=2)
            rows.append(dict(metric=metric, comparison=name, n_proteins=len(finite),
                             mean_gain=mean, ci_low=low, ci_high=high,
                             improved_fraction=float((finite > 0).mean()), unit=UNITS[metric]))
        ax.axvline(0, color=BLACK, linewidth=0.7)
        ax.set_title(LABELS[metric])
        ax.set_xlabel("Protein-equal MAE reduction (%s)" % UNITS[metric])
    axes[0].set_yticks(positions)
    axes[0].set_yticklabels([name for name, _, _ in GAIN_SPECS])
    fig.text(0.01, -0.22,
             "Positive means the full model has the lower protein-equal MAE. Each comparison removes one information source from the retrained\n"
             "family, not one input channel from the frozen network: NN-S removes coordinates and partners, NN-G replaces identity with one unknown\n"
             "token, NN no partner forbids cross-chain messages. The no-identity/no-partner corner was not retrained.",
             fontsize=7, color=BLACK)
    save_figure(fig, out / "fig_information_attribution_gains")
    pd.DataFrame(rows).to_csv(out / "attribution_gains.csv", index=False)


def main() -> None:
    # The step names carry a Unicode minus, which the default Windows console
    # code page cannot encode; fail on the report, never on the artifacts.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe-dir", type=Path,
                        default=ROOT / "output/information_sources_20260919/static_probes")
    parser.add_argument("--residues", type=Path,
                        default=ROOT / "artifacts/evaluation_atlas_hchi_v1/analysis/oof_per_residue.csv.gz")
    parser.add_argument("--chain-map", type=Path,
                        default=ROOT / "output/information_sources_20260919/frozen_analysis/chain_map.csv")
    parser.add_argument("--ablation-dir", type=Path,
                        default=ROOT / "output/information_ablation_analysis",
                        help="Per-residue tables of the retrained information-ablation arms")
    parser.add_argument("--output", type=Path, default=ROOT / "output/information_sources_20260919/figures")
    parser.add_argument("--bootstrap-repeats", type=int, default=5000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    frame = load_frame(args.probe_dir, args.residues, args.chain_map, args.ablation_dir)
    print("Aligned %d residue rows, %d records, %d protein identities"
          % (len(frame), frame.sample_id.nunique(), frame.protein_id.nunique()), flush=True)
    table, pooled, valid_counts = protein_error_table(frame)
    summary = step_table(table)
    table.to_csv(args.output / "source_data_protein_error.csv", index=False)
    summary.to_csv(args.output / "information_ladder_steps.csv", index=False)
    pooled_rows = [dict(metric=metric, method=key, display=display,
                        residue_pooled_mae=pooled[metric][key], valid_residues=valid_counts[metric],
                        unit=UNITS[metric])
                   for metric in METRICS for key, display, _ in ALL_METHODS]
    pd.DataFrame(pooled_rows).to_csv(args.output / "residue_pooled_mae.csv", index=False)

    window_paired = figure_paired(table, args.output)
    figure_steps(table, summary, args.output)
    window_dual = figure_dual_gain(table, args.output)
    figure_attribution_ladder(table, args.output)
    figure_attribution_gains(table, args.output)

    audit = dict(
        status="complete", rows=int(len(frame)), records=int(frame.sample_id.nunique()),
        protein_identities=int(frame.protein_id.nunique()),
        short_peptide_identities=int(frame.loc[frame.full_chain_length < 20, "protein_id"].nunique()),
        valid_residues=valid_counts,
        aggregation="chain MAE -> equal weight over chain copies in a record -> equal weight over records of one protein",
        interval_unit="exact protein sequence; homologous proteins may remain dependent",
        method_definitions={key: text for key, _, text in ALL_METHODS},
        probe_interpretation="Fixed-budget histogram boosting (120 iterations, 31 leaves) per fold and metric, sample-weighted to equal source, system and record.",
        neural_ablation_protocol="Same folds, labels, masks, seed 7 and 40-epoch budget as the frozen atlas_hchi_v1 model; NN-S has about two thirds of the other arms' parameters.",
        scope_limit="Probes and FlexPPI differ in input information and in model family; the probe/FlexPPI gap is not a single-channel ablation.",
        neural_scope_limit="Arms remove information sources, not input channels of the frozen network; identity hiding leaves type information inferable from geometry and partners, and the no-identity/no-partner corner was not retrained.",
        label_protocol="Frozen ATLAS-Hchi-v1 labels; existing ILE and His coverage limitations are retained unchanged.",
        figure_scope="Primary panels compare FlexPPI with the sequence probe because the frozen model is the only source of full 3D geometry.",
        figure_axis_windows={"paired_scatter": window_paired, "dual_gain": window_dual},
        file_sha256={str(path): sha256_file(path) for path in [
            args.probe_dir / "features.npz", args.probe_dir / "feature_schema.json",
            args.residues, args.chain_map,
            *[args.ablation_dir / arm / "oof_per_residue.csv.gz" for arm in ABLATION_FILES.values()],
            Path(__file__)]},
    )
    write_json(args.output / "statistics_and_QA.json", audit)
    print(summary.to_string(index=False), flush=True)
    print(pd.DataFrame(pooled_rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
