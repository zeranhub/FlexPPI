# Asset Confirmation Table
# Peptide profiles: LineTrend production asset is cumulative-publication data,
# incompatible with residue differences; inherit clean axes and external legend.
# Complete paired scatter: basic-plots contains no production implementation.
# Regional replicate plots: basic-plots contains no production implementation.
# All panels retain all predefined positions; only chemically undefined Hchi is masked.

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
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
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "dynamics"
positions = pd.read_csv(DATA / "peptide_mutations_complete_per_position.csv")
metrics = pd.read_csv(DATA / "peptide_mutations_predefined_region_change_metrics.csv")
replicas = pd.read_csv(DATA / "predefined_region_means_per_independent_replica.csv")
CONTRASTS = [("4jfd-4jff", "G4A"), ("4jfe-4jff", "I7A")]
counts = positions.groupby(["contrast", "set"]).size().to_dict()
assert all(counts[(c, "peptide")] == 10 and counts[(c, "CDR")] == 45
           and counts[(c, "pMHC_interface")] == 30 for c, _ in CONTRASTS)

fig, axes = plt.subplots(4, 2, figsize=(183 / 25.4, 222 / 25.4))
fig.subplots_adjust(left=0.11, right=0.985, bottom=0.07, top=0.945,
                    wspace=0.36, hspace=0.84)
validation = []

def panel(ax, label, title):
    ax.text(-0.16, 1.10, label, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom")
    ax.set_title(title, loc="left", pad=7)
    ax.tick_params(length=2.5)

peptide_sets = [positions[(positions.contrast == c) & (positions["set"] == "peptide")]
                .sort_values("common_position") for c, _ in CONTRASTS]
lo = min((d.rmsf_MD_delta - d.rmsf_delta_SE).min() for d in peptide_sets)
hi = max((d.rmsf_MD_delta + d.rmsf_delta_SE).max() for d in peptide_sets)
for col, ((contrast, name), d) in enumerate(zip(CONTRASTS, peptide_sets)):
    ax = axes[0, col]
    panel(ax, "ab"[col], f"{name}: complete peptide")
    x = np.arange(1, 11)
    ax.axhline(0, color=GREY, linestyle="--", linewidth=0.65, zorder=0)
    ax.errorbar(x, d.rmsf_MD_delta, yerr=d.rmsf_delta_SE, fmt="o-",
                color=BLACK, ecolor=GREY, capsize=2, markersize=3,
                linewidth=0.9, elinewidth=0.65, label="MD mean ± SE")
    ax.plot(x, d.rmsf_model_delta, "s--", color=CATEGORICAL[0],
            markersize=3, linewidth=1.0, label="FlexPPI")
    ax.set(xlim=(0.5, 10.5), ylim=(lo - 0.035, hi + 0.035),
           xticks=x, xlabel="Peptide position", ylabel="ΔRMSF (Å)")
    ax.tick_params(axis="x", labelsize=6)
    validation.append({"panel": "ab"[col], "contrast": contrast,
                       "positions": 10, "MD_replicas_per_state": 3,
                       "error_bar": "SE of independent means: sqrt((SD_before²+SD_after²)/3)"})
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.54, 0.998),
           ncol=2, frameon=False, columnspacing=1.8, handlelength=2.0)

def paired_scatter(row, set_name, metric, scope, labels):
    datasets = []
    for contrast, _ in CONTRASTS:
        d = positions[(positions.contrast == contrast) &
                      (positions["set"] == set_name) & positions[f"{metric}_eligible"]]
        datasets.append(d)
    arrays = [d[f"{metric}_{source}_delta"].to_numpy(float)
              for d in datasets for source in ["MD", "model"]]
    limit = max(abs(np.concatenate(arrays))) * 1.16
    for col, ((contrast, name), d) in enumerate(zip(CONTRASTS, datasets)):
        ax = axes[row, col]
        x = d[f"{metric}_MD_delta"].to_numpy(float)
        y = d[f"{metric}_model_delta"].to_numpy(float)
        assert np.isfinite(x).all() and np.isfinite(y).all()
        rho = float(spearmanr(x, y).statistic)
        expected = metrics[(metrics.contrast == contrast) &
                           (metrics.scope == scope) & (metrics.metric == metric)].iloc[0]
        assert len(d) == expected["n"]
        assert abs(rho - expected.delta_Spearman) < 1e-12
        assert int((x * y > 0).sum()) == expected.n_direction_agree
        panel(ax, labels[col], f"{name}: {scope.replace('all45_CDR', 'all 45 CDR sites').replace('all30_pMHC_interface', 'pMHC side chains')}")
        ax.set_title(ax.get_title(loc="left"), loc="left", pad=20)
        ax.plot([-limit, limit], [-limit, limit], linestyle="--", color=GREY, linewidth=0.65)
        ax.axhline(0, color=GREY, linewidth=0.45, alpha=0.55)
        ax.axvline(0, color=GREY, linewidth=0.45, alpha=0.55)
        ax.scatter(x, y, s=15, color=CATEGORICAL[0], alpha=0.8, edgecolor="white", linewidth=0.25)
        unit = "ΔRMSF (Å)" if metric == "rmsf" else "ΔHχ (nats)"
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit),
               xlabel=f"MD {unit}", ylabel=f"FlexPPI {unit}")
        # Statistical text is outside the plotting region, so no point is covered.
        ax.text(0, 1.015, f"ρ = {rho:.3f}; direction {int(expected.n_direction_agree)}/{len(d)}",
                transform=ax.transAxes, fontsize=6.5, va="bottom")
        validation.append({"panel": labels[col], "contrast": contrast,
                           "eligible_positions": len(d), "excluded_undefined_or_changed_type":
                           (30 if set_name == "pMHC_interface" else 45) - len(d),
                           "Spearman": rho, "direction_agree": int(expected.n_direction_agree),
                           "all_predefined_positions_retained_in_CSV": True})

paired_scatter(1, "CDR", "rmsf", "all45_CDR", "cd")
paired_scatter(2, "pMHC_interface", "hchi", "all30_pMHC_interface", "ef")

for col, (metric, ylabel) in enumerate([("rmsf", "CDR3β RMSF (Å)"), ("hchi", "CDR3β Hχ (nats)")]):
    ax = axes[3, col]
    panel(ax, "gh"[col], "CDR3β: individual MD replicates")
    for x, (pdb, label, contrast) in enumerate([("4jff", "ELA", "4jfe-4jff"),
                                              ("4jfd", "G4A", "4jfd-4jff"),
                                              ("4jfe", "I7A", "4jfe-4jff")]):
        rep = replicas[(replicas.pdb == pdb) & (replicas.contrast == contrast) &
                       (replicas.scope == "CDR3β") & (replicas.metric == metric)].sort_values("replica")
        assert len(rep) == 3 and rep.replica.nunique() == 3
        ax.scatter(x + np.array([-0.10, 0, 0.10]), rep.value, s=18,
                   color=BLACK, edgecolor="white", linewidth=0.3, zorder=3)
        ax.plot([x - 0.16, x + 0.16], [rep.value.mean()] * 2, color=BLACK, linewidth=1.1)
        m = metrics[(metrics.contrast == contrast) & (metrics.scope == "CDR3β") &
                    (metrics.metric == metric)].iloc[0]
        model = float(m.model_before if pdb == "4jff" else m.model_after)
        ax.scatter([x], [model], s=24, marker="s", color=CATEGORICAL[0], zorder=4)
        validation.append({"panel": "gh"[col], "pdb": pdb, "metric": metric,
                           "independent_MD_values": rep.value.tolist(), "model": model})
    ax.set(xticks=[0, 1, 2], xticklabels=["ELA", "G4A", "I7A"],
           xlim=(-0.45, 2.45), ylabel=ylabel)
axes[3, 0].legend(handles=[Line2D([], [], color=BLACK, marker="o", linestyle="None", label="MD, n = 3 / state"),
                          Line2D([], [], color=CATEGORICAL[0], marker="s", linestyle="None", label="FlexPPI")],
                  loc="upper center", bbox_to_anchor=(1.18, -0.27), ncol=2, fontsize=7)

save_cns_figure(fig, ROOT / "figS6_peptide_dynamic_validation")
fig.savefig(ROOT / "figS6_peptide_dynamic_validation.svg", bbox_inches="tight", dpi=300)
plt.close(fig)
report = {"figure": "figS6_peptide_dynamic_validation", "width_mm": 183, "height_mm": 222,
          "data_selection": "All 10 peptide, 45 CDR, 30 pMHC interface positions per contrast.",
          "missing_Hchi": "Chemically undefined or changed amino-acid type are excluded only from Hchi numeric comparison.",
          "independence": "3 independent MD trajectories per state, not residues/frames/9 pair differences.",
          "p_values": "None; descriptive within-case verification, no biological population inference.",
          "panels": validation,
          "source_sha256": {name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
                             for name in ["peptide_mutations_complete_per_position.csv",
                                          "peptide_mutations_predefined_region_change_metrics.csv",
                                          "predefined_region_means_per_independent_replica.csv"]}}
(ROOT / "dynamic_figure_statistics.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({"status": "rendered", "panels": 8, "validation_rows": len(validation)}, ensure_ascii=False))
