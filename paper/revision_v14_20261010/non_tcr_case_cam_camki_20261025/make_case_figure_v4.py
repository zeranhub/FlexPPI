"""Standalone CaM–CaMKI case figure (v4, free vs bound).

a  two-state schematic: free peptide (coil, model) -> CaM -> bound complex (PyMOL)
b  peptide RMSF: free-state prediction, bound-state MD / Full / Intrachain
c  interface accuracy of the bound-state model
d  mean peptide RMSF: free prediction vs bound MD / Full / Intrachain
e  partner-removal controls

The free-state curve is the same Full checkpoint with CaM removed (input
occlusion). It is a model prediction, not a free-state MD reference; the original
1MXE paper provides only the qualitative free-state disorder.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

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
INK = "#262626"
LIGHT = "#D9DDE1"
FREE = "#8A8A8A"

HERE = Path(__file__).resolve().parent
DATA = pd.read_csv(HERE / "cam_camki_per_residue.csv", dtype={"author_residue_id": str})
MET = pd.read_csv(HERE / "case_metrics.csv")
CTRL = pd.read_csv(HERE / "partner_removal_controls.csv")
W, H = 183.0, 155.0


def ax_mm(fig, x, y, w, h):
    return fig.add_axes([x / W, y / H, w / W, h / H])


def panel_two_state(fig):
    # free peptide: a schematic disordered ensemble
    ax_free = ax_mm(fig, 7, 126, 32, 10)
    t = np.linspace(0, 6 * np.pi, 400)
    for k, alpha in enumerate([0.25, 0.45, 0.75]):
        ax_free.plot(t, np.sin(t + k * 1.1) * 0.26 + (k - 1) * 0.18,
                     color=WARM, alpha=alpha, lw=1.1, solid_capstyle="round")
    ax_free.set_xlim(0, 6 * np.pi)
    ax_free.set_ylim(-0.85, 0.85)
    ax_free.set_axis_off()
    fig.text(7 / W, 137 / H, "free (model)", fontsize=5.8, color=FREE, ha="left")

    # arrow: CaM binding
    ax_arrow = ax_mm(fig, 18, 118, 10, 8)
    ax_arrow.set_xlim(0, 1)
    ax_arrow.set_ylim(0, 1)
    ax_arrow.annotate("", xy=(0.5, 0.05), xytext=(0.5, 0.95),
                      arrowprops=dict(arrowstyle="-|>", color="#9A9A9A", lw=0.9,
                                      mutation_scale=8, shrinkA=0, shrinkB=0))
    ax_arrow.text(0.62, 0.50, "CaM", fontsize=5.8, color="#777777", va="center")
    ax_arrow.set_axis_off()

    # bound complex
    ax_bound = ax_mm(fig, 11, 94, 24, 24)
    ax_bound.imshow(plt.imread(HERE / "cam_structure.png"))
    ax_bound.set_axis_off()
    fig.text(23 / W, 92 / H, "bound (MD)", fontsize=5.8, color=MD, ha="center")


def panel_peptide(fig):
    ax = ax_mm(fig, 52, 94, 66, 46)
    q = DATA[DATA.chain == "B"].sort_values("author_residue_id")
    x = q.author_residue_id.astype(int).to_numpy()
    ax.axvspan(148, 166, color=WARM, alpha=0.06, lw=0, zorder=0)
    ax.plot(x, q.full_no_partner_rmsf, color=FREE, lw=1.2, ls=(0, (2.5, 1.8)), zorder=1)
    ax.plot(x, q.full_rmsf, color=BLUE, lw=1.3, zorder=2)
    ax.plot(x, q.intra_rmsf, color=WARM, lw=1.1, ls=(0, (4.0, 2.0)), zorder=1)
    ax.plot(x, q.target_rmsf, color=MD, lw=1.9, zorder=5)
    ax.set_xlim(144.2, 169.8)
    ax.set_ylim(0, 9.6)
    ax.set_xticks([145, 150, 155, 160, 165, 169])
    ax.set_yticks([0, 2, 4, 6, 8])
    ax.set_xlabel("CaMKI peptide residue", labelpad=2)
    ax.set_ylabel("Cα RMSF (Å)", labelpad=2)
    ax.text(145.3, 8.85, "free (model)", color=FREE, fontsize=5.8, va="bottom")
    ax.text(145.3, 3.55, "bound Intra", color=WARM, fontsize=5.8, va="bottom")
    ax.text(145.3, 2.35, "bound MD / Full", color=INK, fontsize=5.8, va="bottom")


def panel_interface(fig):
    ax = ax_mm(fig, 124, 94, 58, 46)
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
                    ha="center", va="bottom", fontsize=6.2, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels([f"interface\n(n={int(interface.n)})", f"non-interface\n(n={int(outside.n)})"], fontsize=6.7)
    ax.set_ylabel("MAE (Å)", labelpad=2)
    ax.set_ylim(0, 1.02)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])


def panel_means(fig):
    ax = ax_mm(fig, 1, 10, 87, 68)
    labels = ["free\n(model)", "bound\nMD", "bound\nFull", "bound\nIntra"]
    vals = [4.219, 0.789, 0.856, 2.003]
    colors = ["#B8B8B8", MD, BLUE, WARM]
    ax.bar(range(4), vals, color=colors, width=0.58, edgecolor="white", linewidth=0.5, zorder=2)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.10, f"{v:.2f}", ha="center", va="bottom", fontsize=6.6, color=INK)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel("mean peptide RMSF (Å)", labelpad=2)
    ax.set_ylim(0, 5.3)
    ax.set_yticks([0, 1, 2, 3, 4, 5])
    ax.annotate("", xy=(1.98, 4.72), xytext=(1.02, 4.72),
                arrowprops=dict(arrowstyle="<->", color="#9A9A9A", lw=0.6, shrinkA=0, shrinkB=0))
    ax.text(1.5, 4.84, "CaM binding  −3.36 Å", ha="center", va="bottom", fontsize=6.2, color="#666666")


def panel_controls(fig):
    ax = ax_mm(fig, 93, 10, 89, 68)
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
    ax.text(5.42, 1.0, "folded\npartners", ha="right", va="center", fontsize=5.7, color="#7A8B99")
    ax.text(5.42, 5.5, "short\nchains", ha="right", va="center", fontsize=5.7, color="#8A8A8A")


def main():
    fig = plt.figure(figsize=(W / 25.4, H / 25.4), facecolor="white")
    fig.text(1 / W, 151 / H, "CaM–CaMKI | Partner-dependent folding of a regulatory peptide",
             fontsize=9, fontweight="bold", color=INK)
    fig.text(1 / W, 147 / H, "1MXE / A001K · free-state prediction and bound-state MD reference · frozen OOF model",
             fontsize=6.2, color="#777777")

    for x, tag, title in [(1, "a", "Free → bound"),
                          (52, "b", "CaMKI peptide backbone"),
                          (124, "c", "Bound-state interface accuracy"),
                          (1, "d", "Mean peptide mobility"),
                          (93, "e", "Partner-removal controls")]:
        y = 142 if tag in "abc" else 82
        fig.text(x / W, y / H, tag, fontsize=10, fontweight="bold", color=INK)
        fig.text((x + 6) / W, y / H, title, fontsize=7.6, color=INK)

    panel_two_state(fig)
    panel_peptide(fig)
    panel_interface(fig)
    panel_means(fig)
    panel_controls(fig)

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(HERE / f"fig_cam_camki_case_v4.{ext}", dpi=450)
    print("wrote", HERE / "fig_cam_camki_case_v4.png")


if __name__ == "__main__":
    main()
