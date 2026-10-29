"""Standalone CaM–CaMKI case figure (publication-oriented draft).

Panels:
a  structure (3Dmol cartoon, CA trace not schematic)
b  per-residue peptide RMSF: MD, Full, Full minus partner, Intrachain
c  interface vs non-interface RMSF MAE
d  mean peptide RMSF across conditions
e  partner-removal controls across systems

All data are observed MD labels or frozen model predictions. The partner-removed
curve is a same-checkpoint input-occlusion inference, documented in the report.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7.5, "axes.titlesize": 7.5, "axes.labelsize": 7.5,
    "xtick.labelsize": 6.6, "ytick.labelsize": 6.6, "legend.fontsize": 6.2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
    "savefig.bbox": "tight", "savefig.dpi": 450,
})
BLUE = "#2166AC"
WARM = "#B87936"
MD = "#62676B"
MID = "#92999F"
INK = "#262626"
LIGHT = "#D9DDE1"
PALE = "#F2F3F4"

HERE = Path(__file__).resolve().parent
DATA = pd.read_csv(HERE / "cam_camki_per_residue.csv", dtype={"author_residue_id": str})
MET = pd.read_csv(HERE / "case_metrics.csv")
CTRL = pd.read_csv(HERE / "partner_removal_controls.csv")
W, H = 183.0, 136.0


def ax_mm(fig, x, y, w, h):
    return fig.add_axes([x / W, y / H, w / W, h / H])


def panel_structure(fig):
    ax = ax_mm(fig, 1, 76, 50, 50)
    img = plt.imread(HERE / "cam_structure.png")
    ax.imshow(img)
    ax.set_axis_off()
    handles = [
        Line2D([], [], color=LIGHT, lw=4.0, label="CaM"),
        Line2D([], [], color="#D98A45", lw=3.4, label="CaMKI peptide"),
        Line2D([], [], marker="s", ls="", color="#8FB3D1", ms=4.0, label="interface"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.015),
              ncol=3, handlelength=1.5, columnspacing=0.7, fontsize=5.8)


def panel_peptide(fig):
    ax = ax_mm(fig, 56, 76, 64, 50)
    q = DATA[DATA.chain == "B"].sort_values("author_residue_id")
    x = q.author_residue_id.astype(int).to_numpy()
    ax.axvspan(148, 166, color=WARM, alpha=0.07, lw=0, zorder=0)
    ax.plot(x, q.full_rmsf, color=BLUE, lw=1.4, label="Full complex", zorder=2)
    ax.plot(x, q.full_no_partner_rmsf, color=BLUE, lw=1.2, ls=(0, (3.0, 1.8)),
            label="Full, partner removed", zorder=1)
    ax.plot(x, q.intra_rmsf, color=WARM, lw=1.2, ls=(0, (4.0, 2.0)), label="Intrachain", zorder=1)
    ax.plot(x, q.target_rmsf, color=MD, lw=2.0, label="MD", zorder=5)
    ax.set_xlim(144.2, 169.8)
    ax.set_ylim(0, 9.3)
    ax.set_xticks([145, 150, 155, 160, 165, 169])
    ax.set_yticks([0, 2, 4, 6, 8])
    ax.set_xlabel("CaMKI peptide residue", labelpad=2)
    ax.set_ylabel("Cα RMSF (Å)", labelpad=2)
    handles = [
        Line2D([], [], color=MD, lw=2.0, label="MD"),
        Line2D([], [], color=BLUE, lw=1.4, label="Full complex"),
        Line2D([], [], color=BLUE, lw=1.2, ls=(0, (3.0, 1.8)), label="Full, partner removed"),
        Line2D([], [], color=WARM, lw=1.2, ls=(0, (4.0, 2.0)), label="Intrachain"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.0),
              ncol=2, handlelength=1.7, columnspacing=0.9, borderaxespad=0.15)
    ax.text(157, 8.55, "bound helix", ha="center", va="center", fontsize=6.0, color=WARM)


def panel_interface(fig):
    ax = ax_mm(fig, 125, 76, 57, 50)
    q = MET[(MET.model == "main_oof") & (MET.metric == "rmsf") &
            MET.subset.isin(["interface_ca8", "outside_interface"])]
    interface = q[q.subset == "interface_ca8"].iloc[0]
    outside = q[q.subset == "outside_interface"].iloc[0]
    x = np.arange(2)
    width = 0.30
    vals_full = [interface.full_mae, outside.full_mae]
    vals_intra = [interface.intra_mae, outside.intra_mae]
    b1 = ax.bar(x - width / 2, vals_full, width, color=BLUE, edgecolor="white", linewidth=0.4, zorder=2)
    b2 = ax.bar(x + width / 2, vals_intra, width, color=WARM, edgecolor="white", linewidth=0.4, zorder=2)
    for bars in (b1, b2):
        for b in bars:
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.018, f"{b.get_height():.2f}",
                    ha="center", va="bottom", fontsize=6.3, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels([f"interface\n(n={int(interface.n)})", f"non-interface\n(n={int(outside.n)})"], fontsize=6.8)
    ax.set_ylabel("RMSF MAE (Å)", labelpad=2)
    ax.set_ylim(0, 1.02)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.text(0.02, 0.97, "Full", color=BLUE, transform=ax.transAxes, fontsize=6.2, va="top", ha="left")
    ax.text(0.02, 0.90, "Intra", color=WARM, transform=ax.transAxes, fontsize=6.2, va="top", ha="left")


def panel_means(fig):
    ax = ax_mm(fig, 1, 10, 87, 56)
    labels = ["MD", "Full\ncomplex", "Full\npartner\nremoved", "Intrachain"]
    vals = [0.789, 0.856, 4.219, 2.003]
    colors = [MD, BLUE, "#A7C4DD", WARM]
    ax.bar(range(4), vals, color=colors, width=0.58, edgecolor="white", linewidth=0.5, zorder=2)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.10, f"{v:.2f}", ha="center", va="bottom", fontsize=6.6, color=INK)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel("mean peptide RMSF (Å)", labelpad=2)
    ax.set_ylim(0, 5.25)
    ax.set_yticks([0, 1, 2, 3, 4, 5])
    ax.annotate("", xy=(1.98, 4.58), xytext=(1.02, 4.58),
                arrowprops=dict(arrowstyle="<->", color="#9A9A9A", lw=0.6, shrinkA=0, shrinkB=0))
    ax.text(1.5, 4.70, "partner removed  +3.36 Å", ha="center", va="bottom", fontsize=6.3, color="#666666")


def panel_controls(fig):
    ax = ax_mm(fig, 93, 10, 89, 56)
    order = ["CaM–CaMKI (1MXE)", "FSH–FSHR (1XWD)", "IL-5–IL-5Rα (3QT2)", "Cdt1–Geminin (2ZXX)"]
    short = ["3WWQ ubiquitin–FAAP20", "2HWN PKA–peptide", "3P8M DLC–GCN4", "6A7V VapBC"]
    ctrl = CTRL.set_index("system").loc[order + short].reset_index()
    labels = ["1MXE", "1XWD", "3QT2", "2ZXX", "3WWQ", "2HWN", "3P8M", "6A7V"]
    y = np.arange(len(labels))
    ax.axhspan(-0.5, 2.5, color="#EAF2F8", lw=0, zorder=0)
    ax.axhspan(3.5, 7.5, color="#F4F4F4", lw=0, zorder=0)
    colors = [WARM] + ["#7FA8C9"] * 3 + ["#C9C9C9"] * 4
    ax.barh(y, ctrl.delta, color=colors, height=0.60, edgecolor="white", linewidth=0.4, zorder=2)
    for yi, v in zip(y, ctrl.delta):
        ax.text(v + 0.08, yi, f"{v:.2f}", va="center", ha="left", fontsize=5.9, color="#555555")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.4)
    ax.invert_yaxis()
    ax.set_xlim(0, 5.5)
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_xlabel("Δ mean RMSF (Å)", labelpad=2)
    ax.text(5.42, 1.0, "folded\npartners", ha="right", va="center", fontsize=5.8, color="#7A8B99")
    ax.text(5.42, 5.5, "short\nchains", ha="right", va="center", fontsize=5.8, color="#8A8A8A")


def main():
    fig = plt.figure(figsize=(W / 25.4, H / 25.4), facecolor="white")
    fig.text(1 / W, 132 / H, "CaM–CaMKI | Partner-dependent folding of a regulatory peptide",
             fontsize=9, fontweight="bold", color=INK)
    fig.text(1 / W, 128 / H, "1MXE / A001K · per-residue flexibility from the static complex · MD reference and frozen OOF predictions",
             fontsize=6.2, color="#777777")

    for x, tag, title in [(1, "a", "Complex and interface"),
                          (56, "b", "CaMKI peptide backbone"),
                          (125, "c", "Interface accuracy"),
                          (1, "d", "Mean peptide mobility"),
                          (93, "e", "Partner-removal controls")]:
        y = 125 if tag in "abc" else 69
        fig.text(x / W, y / H, tag, fontsize=10, fontweight="bold", color=INK)
        fig.text((x + 6) / W, y / H, title, fontsize=7.8, color=INK)

    panel_structure(fig)
    panel_peptide(fig)
    panel_interface(fig)
    panel_means(fig)
    panel_controls(fig)

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(HERE / f"fig_cam_camki_case_v2.{ext}", dpi=450)
    print("wrote", HERE / "fig_cam_camki_case_v2.png")


if __name__ == "__main__":
    main()
