"""Draft case figure for the CaM–CaMKI peptide (1MXE / A001K).

(a) CA-trace schematic; (b) peptide RMSF with partner-removal ablation;
(c) interface versus non-interface MAE; (d) mean peptide RMSF across conditions;
(e) partner-removal controls. All values are observed MD labels or frozen model
predictions; the partner-removed curve is a same-checkpoint input-occlusion
inference documented in the case report.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.interpolate import splprep, splev

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.4,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
    "savefig.bbox": "tight", "savefig.dpi": 300,
})
BLUE = "#2166AC"
WARM = "#B87936"
MD = "#62676B"
MID = "#92999F"
INK = "#262626"

CASE = Path(r"E:\Project\FlexPPI\paper\revision_v14_20261010\non_tcr_case_cam_camki_20261025")
DATA = pd.read_csv(CASE / "cam_camki_per_residue.csv", dtype={"author_residue_id": str})
MET = pd.read_csv(CASE / "case_metrics.csv")
PDB = Path(r"E:\Project\FlexPPI\experiments\_sources\baselines\cohort\structures\record_aba5e1f3db118b0aba8895de.pdb")
W, H = 183.0, 128.0


def ax_mm(fig, x, y, w, h):
    return fig.add_axes([x / W, y / H, w / W, h / H])


def ca_trace():
    rows = []
    for line in PDB.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            rows.append((line[21], int(line[22:26]), float(line[30:38]), float(line[38:46]), float(line[46:54])))
    df = pd.DataFrame(rows, columns=["chain", "res", "x", "y", "z"])
    xyz = df[["x", "y", "z"]].to_numpy()
    centered = xyz - xyz.mean(0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    proj = centered @ vt[:2].T
    df["px"], df["py"] = proj[:, 0], proj[:, 1]
    return df


def smooth(ax, x, y, color, lw, zorder):
    tck, _ = splprep([x, y], s=0, k=min(3, len(x) - 1))
    u = np.linspace(0, 1, 500)
    xs, ys = splev(u, tck)
    ax.plot(xs, ys, color="white", lw=lw + 2.4, solid_capstyle="round", zorder=zorder - 0.1)
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="round", zorder=zorder)


def panel_structure(fig):
    ax = ax_mm(fig, 1, 76, 52, 42)
    df = ca_trace()
    cam = df[df.chain == "A"].sort_values("res")
    pep = df[df.chain == "B"].sort_values("res")
    smooth(ax, cam.px, cam.py, "#CDD1D5", 5.4, 1)
    smooth(ax, pep.px, pep.py, WARM, 3.2, 3)
    inter = DATA[DATA.interface_ca8.fillna(False)]
    for ch, color, size in [("A", "#5B8DB8", 6), ("B", "#8C3B1E", 8)]:
        q = inter[inter.chain == ch]
        pts = df[df.chain == ch].set_index("res")
        xs = [pts.loc[int(r), "px"] for r in q.author_residue_id if int(r) in pts.index]
        ys = [pts.loc[int(r), "py"] for r in q.author_residue_id if int(r) in pts.index]
        ax.scatter(xs, ys, s=size, color=color, edgecolor="white", linewidth=0.25, zorder=5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-31, 31)
    ax.set_ylim(-28, 27)
    handles = [Line2D([], [], color="#CDD1D5", lw=4, label="CaM"),
               Line2D([], [], color=WARM, lw=3, label="CaMKI peptide"),
               Line2D([], [], marker="o", ls="", color="#5B8DB8", ms=3.2, label="interface")]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=3,
              handlelength=1.6, columnspacing=0.8, fontsize=6.3)


def panel_peptide(fig):
    ax = ax_mm(fig, 59, 76, 60, 42)
    q = DATA[DATA.chain == "B"].sort_values("author_residue_id")
    x = q.author_residue_id.astype(int).to_numpy()
    ax.axvspan(148, 166, color=WARM, alpha=0.06, lw=0)
    ax.plot(x, q.full_rmsf, color=BLUE, lw=1.5, label="Full complex", zorder=2)
    ax.plot(x, q.full_no_partner_rmsf, color=BLUE, lw=1.3, ls=(0, (3, 1.8)),
            label="Full, partner removed", zorder=1)
    ax.plot(x, q.intra_rmsf, color=WARM, lw=1.3, ls="--", label="Intrachain", zorder=1)
    ax.plot(x, q.target_rmsf, color=MD, lw=2.0, label="MD", zorder=5)
    ax.set_xlim(144.2, 169.8)
    ax.set_ylim(0, 11.2)
    ax.set_xticks([145, 150, 155, 160, 165, 169])
    ax.set_xlabel("CaMKI peptide residue")
    ax.set_ylabel("Cα RMSF (Å)")
    ax.set_yticks([0, 2, 4, 6, 8])
    ax.legend(loc="upper left", ncol=2, handlelength=1.5, columnspacing=0.8,
              borderaxespad=0.2, fontsize=6.3)
    ax.text(157, 8.6, "bound helix", ha="center", fontsize=6.2, color=WARM)


def panel_interface(fig):
    ax = ax_mm(fig, 126, 76, 56, 42)
    q = MET[(MET.model == "main_oof") & (MET.metric == "rmsf") & MET.subset.isin(["interface_ca8", "outside_interface"])]
    interface = q[q.subset == "interface_ca8"].iloc[0]
    outside = q[q.subset == "outside_interface"].iloc[0]
    x = np.arange(2)
    width = 0.32
    full = [interface.full_mae, outside.full_mae]
    intra = [interface.intra_mae, outside.intra_mae]
    b1 = ax.bar(x - width / 2, full, width, color=BLUE, edgecolor="white", linewidth=0.4)
    b2 = ax.bar(x + width / 2, intra, width, color=WARM, edgecolor="white", linewidth=0.4)
    for bars in (b1, b2):
        for b in bars:
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.02, f"{b.get_height():.2f}",
                    ha="center", va="bottom", fontsize=6.6)
    ax.set_xticks(x)
    ax.set_xticklabels([f"interface\n(n={int(interface.n)})", f"non-interface\n(n={int(outside.n)})"], fontsize=7)
    ax.set_ylabel("MAE (Å)", fontsize=7.5)
    ax.set_ylim(0, 1.05)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.text(0.02, 0.97, "Full", color=BLUE, transform=ax.transAxes, fontsize=6.4, va="top")
    ax.text(0.02, 0.90, "Intra", color=WARM, transform=ax.transAxes, fontsize=6.4, va="top")


def panel_means(fig):
    ax = ax_mm(fig, 59, 10, 55, 46)
    labels = ["MD", "Full\ncomplex", "Full\npartner\nremoved", "Intrachain"]
    vals = [0.789, 0.856, 4.219, 2.003]
    colors = [MD, BLUE, "#A7C4DD", WARM]
    ax.bar(range(4), vals, color=colors, width=0.60, edgecolor="white", linewidth=0.5)
    for i, v in enumerate(vals):
        if i == 3:
            ax.text(i, v - 0.25, f"{v:.2f}", ha="center", va="top", fontsize=6.6, color="white")
        else:
            ax.text(i, v + 0.12, f"{v:.2f}", ha="center", fontsize=6.6, color=INK)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=6.4)
    ax.set_ylabel("mean peptide RMSF (Å)")
    ax.set_ylim(0, 5.3)
    ax.set_yticks([0, 1, 2, 3, 4, 5])
    ax.annotate("", xy=(1.98, 4.65), xytext=(1.02, 4.65),
                arrowprops=dict(arrowstyle="<->", color="#999999", lw=0.6))
    ax.text(1.5, 4.78, "partner removed: +3.36 Å", ha="center", fontsize=6.2, color="#666666")


def panel_controls(fig):
    ax = ax_mm(fig, 120, 10, 62, 46)
    ctrl = pd.read_csv(CASE / "partner_removal_controls.csv")
    order = ["CaM–CaMKI (1MXE)", "FSH–FSHR (1XWD)", "IL-5–IL-5Rα (3QT2)", "Cdt1–Geminin (2ZXX)"]
    short = ["3WWQ ubiquitin–FAAP20", "2HWN PKA–peptide", "3P8M DLC–GCN4", "6A7V VapBC"]
    ctrl = ctrl.set_index("system").loc[order + short].reset_index()
    labels = ["1MXE", "1XWD", "3QT2", "2ZXX", "3WWQ", "2HWN", "3P8M", "6A7V"]
    y = np.arange(len(labels))
    colors = [WARM] + ["#7FA8C9"] * 3 + ["#CCCCCC"] * 4
    ax.barh(y, ctrl.delta, color=colors, height=0.60, edgecolor="white", linewidth=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.4)
    ax.invert_yaxis()
    ax.set_xlim(0, 5.4)
    ax.set_xlabel("Δ mean RMSF (Å)", fontsize=6.8)
    ax.tick_params(axis="x", labelsize=6.2)
    ax.axhline(3.5, color="#DDDDDD", lw=0.7)


def main():
    fig = plt.figure(figsize=(W / 25.4, H / 25.4), facecolor="white")
    fig.text(1 / W, 125 / H, "CaM–CaMKI | Partner-dependent folding of a regulatory peptide", fontsize=9, fontweight="bold")
    for x, tag, title, y in [(1, "a", "Complex and interface", 120), (59, "b", "CaMKI peptide backbone", 120),
                             (126, "c", "Interface accuracy", 120), (59, "d", "Mean peptide mobility", 60),
                             (120, "e", "Partner-removal controls", 60)]:
        fig.text(x / W, y / H, tag, fontsize=10, fontweight="bold")
        fig.text((x + 6) / W, y / H, title, fontsize=8)
    panel_structure(fig)
    panel_peptide(fig)
    panel_interface(fig)
    panel_means(fig)
    panel_controls(fig)
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(CASE / f"fig_cam_camki_case.{ext}", bbox_inches="tight", dpi=300)
    print("wrote figure")


if __name__ == "__main__":
    main()
