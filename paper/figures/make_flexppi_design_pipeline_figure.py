# -*- coding: utf-8 -*-
"""Generate the FlexPPI design/screening pipeline overview figure.

Figure type: Solution overview (schematic, no data axes except the small
Pareto sketch in panel d). Layout follows the literature patterns collected in
docs/pipeline_figure_literature_review.md:

  - two converging entrance lanes (de novo / mutation)         [RFdiffusion, BindCraft]
  - structure ensemble + fast triage funnel with 3 gates        [BindCraft funnel]
  - FlexPPI dynamics-aware reranking as the visual focus        [BRD3 hybrid AI/physics slot]
  - Pareto selection + validation + single feedback arc         [EVOLVEpro loop, DynamicMPNN syntax]

Graphic grammar: solid arrows = main data flow; grey dashed = evaluation /
consistency checks; one red dashed arc = experimental feedback. Okabe-Ito
colour-blind-safe palette with text/line double encoding.

Box text is wrapped with measured glyph widths (TextPath) so labels never
overflow their boxes. FancyBboxPatch boxstyle pad/rounding_size are in DATA
units (inches on this canvas), so they are converted from mm.

Run from the repo root:
    python paper\\figures\\make_flexppi_design_pipeline_figure.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.textpath import TextPath

# ---------------------------------------------------------------- style ---
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

# Okabe-Ito palette + semantic greys
TARGET = "#0072B2"    # target blue
BINDER = "#E69F00"    # binder orange
PASS = "#009E73"      # surviving candidates green
FLEX = "#CC79A7"      # FlexPPI purple
FEED = "#D55E00"      # experimental feedback red
INK = "#24313D"       # primary text
GREY = "#6B7280"      # secondary / uncertainty
FILL_GREY = "#E5E7EB"
FILL_BLUE = "#D6E8F2"
FILL_ORANGE = "#FBE3C4"
FILL_GREEN = "#D9EEE3"
FILL_PURPLE = "#F0E4F5"
PANEL_BG = "#F7F8FA"
PANEL_EC = "#D1D5DB"

MM = 25.4
W, H = 178 / MM, 95 / MM          # double-column 178 mm x 95 mm

_FP = None


def _measure(s, fontsize):
    """Text width in mm for the configured sans-serif font."""
    global _FP
    if _FP is None:
        _FP = FontProperties(family=plt.rcParams["font.sans-serif"])
    tp = TextPath((0, 0), s, size=fontsize, prop=_FP)
    return tp.get_extents().width / 72.0 * MM


def wrap_to(s, max_w_mm, fontsize):
    """Word-wrap a string so every line fits max_w_mm."""
    words = s.split()
    lines, cur = [], ""
    for w in words:
        test = (cur + " " + w).strip()
        if _measure(test, fontsize) <= max_w_mm:
            cur = test
        else:
            if cur:
                lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    return "\n".join(lines)


# ---------------------------------------------------------------- helpers --
def mm_ax(ax, x, y, w, h, edgecolor=INK, facecolor="white", lw=0.7):
    """Add a rounded module box; coordinates in mm."""
    pad, rnd = 0.8 / MM, 1.2 / MM   # boxstyle args are in DATA units (inches)
    box = FancyBboxPatch(
        (x / MM, y / MM), w / MM, h / MM,
        boxstyle=f"round,pad={pad:.5f},rounding_size={rnd:.5f}",
        linewidth=lw, edgecolor=edgecolor, facecolor=facecolor, zorder=3,
    )
    ax.add_patch(box)
    return box


def mm_txt(ax, x, y, s, size=5.6, color=INK, weight="normal",
           ha="center", va="center", style="normal", zorder=5, rot=0):
    ax.text(x / MM, y / MM, s, fontsize=size, color=color,
            fontweight=weight, ha=ha, va=va, fontstyle=style,
            zorder=zorder, rotation=rot, linespacing=1.15)


def mm_txt_w(ax, x, y, s, max_w_mm, size=5.0, color=INK, weight="normal",
             zorder=5):
    """Centered wrapped text constrained to max_w_mm."""
    mm_txt(ax, x, y, wrap_to(s, max_w_mm, size), size=size, color=color,
           weight=weight, zorder=zorder)


def mm_line(ax, pts, color=INK, lw=1.1, style="-", zorder=4):
    xs = [p[0] / MM for p in pts]
    ys = [p[1] / MM for p in pts]
    ax.plot(xs, ys, color=color, lw=lw, linestyle=style, zorder=zorder,
            solid_capstyle="round")


def mm_arrow(ax, pts, color=INK, lw=1.1, style="-", zorder=4):
    """Polyline through mm points with an arrowhead on the final segment."""
    if len(pts) < 2:
        return
    mm_line(ax, pts, color=color, lw=lw, style=style, zorder=zorder)
    (x1, y1), (x2, y2) = pts[-2], pts[-1]
    ax.annotate(
        "", xy=(x2 / MM, y2 / MM), xytext=(x1 / MM, y1 / MM),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                        linestyle=style, mutation_scale=7,
                        shrinkA=0, shrinkB=0),
        zorder=zorder + 1,
    )


def ribbon(ax, x, y, w, h, colors, gap=0.25):
    """Per-residue style strip built from small vector rectangles."""
    n = len(colors)
    seg = (w - gap * (n - 1)) / n
    for i, c in enumerate(colors):
        ax.add_patch(Rectangle(
            ((x + i * (seg + gap)) / MM, y / MM),
            seg / MM, h / MM, facecolor=c, edgecolor="none", zorder=3,
        ))


def panel(ax, x0, y0, x1, y1):
    pad, rnd = 0.6 / MM, 1.2 / MM   # boxstyle args are in DATA units (inches)
    ax.add_patch(FancyBboxPatch(
        (x0 / MM, y0 / MM), (x1 - x0) / MM, (y1 - y0) / MM,
        boxstyle=f"round,pad={pad:.5f},rounding_size={rnd:.5f}",
        linewidth=0.6, edgecolor=PANEL_EC, facecolor=PANEL_BG, zorder=1,
    ))


# ---------------------------------------------------------------- figure ---
fig = plt.figure(figsize=(W, H), dpi=100)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")

# ---- panels ---------------------------------------------------------------
panel(ax, 2.0, 6.5, 43.0, 90.0)    # (a) candidate space
panel(ax, 45.5, 6.5, 85.0, 90.0)   # (b) structure & fast triage
panel(ax, 88.0, 6.5, 136.0, 90.0)  # (c) dynamics-aware reranking
panel(ax, 139.0, 6.5, 176.0, 90.0) # (d) selection & validation

mm_txt(ax, 22.5, 93.0, "(a) Candidate space", size=8, weight="bold")
mm_txt(ax, 65.0, 93.0, "(b) Structure & fast triage", size=8, weight="bold")
mm_txt(ax, 112.0, 93.0, "(c) Dynamics-aware reranking", size=8, weight="bold")
mm_txt(ax, 157.5, 93.0, "(d) Selection & validation loop", size=8, weight="bold")

# ---- (a) two entrance lanes -------------------------------------------------
LANE1_Y0, LANE1_H = 66.0, 13.0
LANE2_Y0, LANE2_H = 44.0, 13.0
# three 11 mm boxes with 1.5 mm gaps inside panel a (x 2..43)
BX = [(3.5, 14.5), (16.0, 27.0), (28.5, 39.5)]

mm_txt(ax, 21.5, 81.6, "De novo binder design", size=4.8, color=GREY)
mm_txt(ax, 21.5, 59.8, "Mutation screening", size=4.8, color=GREY)

# lane 1: de novo (left -> right)
mm_txt_w(ax, (BX[0][0] + BX[0][1]) / 2, LANE1_Y0 + LANE1_H / 2,
         "target + epitope", 10.5, size=4.3, color=TARGET, weight="bold")
mm_ax(ax, *BX[0], LANE1_Y0, LANE1_H)
mm_txt_w(ax, (BX[1][0] + BX[1][1]) / 2, LANE1_Y0 + LANE1_H / 2,
         "backbone generation (RFdiffusion / BindCraft)", 10.5, size=4.3)
mm_ax(ax, *BX[1], LANE1_Y0, LANE1_H)
mm_txt_w(ax, (BX[2][0] + BX[2][1]) / 2, LANE1_Y0 + LANE1_H / 2,
         "sequence design (ProteinMPNN)", 10.5, size=4.3)
mm_ax(ax, *BX[2], LANE1_Y0, LANE1_H)
mm_arrow(ax, [(14.5, 72.5), (16.0, 72.5)])
mm_arrow(ax, [(27.0, 72.5), (28.5, 72.5)])

# lane 2: mutation (left -> right)
mm_txt_w(ax, (BX[0][0] + BX[0][1]) / 2, LANE2_Y0 + LANE2_H / 2,
         "WT complex", 10.5, size=4.3, color=BINDER, weight="bold")
mm_ax(ax, *BX[0], LANE2_Y0, LANE2_H)
mm_txt_w(ax, (BX[1][0] + BX[1][1]) / 2, LANE2_Y0 + LANE2_H / 2,
         "site selection (interface / hotspot)", 10.5, size=4.3)
mm_ax(ax, *BX[1], LANE2_Y0, LANE2_H)
mm_txt_w(ax, (BX[2][0] + BX[2][1]) / 2, LANE2_Y0 + LANE2_H / 2,
         "single & combinatorial variants", 10.5, size=4.3)
mm_ax(ax, *BX[2], LANE2_Y0, LANE2_H)
mm_arrow(ax, [(14.5, 50.5), (16.0, 50.5)])
mm_arrow(ax, [(27.0, 50.5), (28.5, 50.5)])

# merge box (bottom of panel a)
mm_ax(ax, 8.0, 13.5, 28.0, 11.0, edgecolor=PASS, lw=1.0)
mm_txt_w(ax, 22.0, 19.0, "candidate sequences + complex hypotheses",
         26.0, size=4.7, color=PASS, weight="bold")

# L-shaped orthogonal merge connectors along the right corridor
mm_arrow(ax, [(39.5, 72.5), (42.0, 72.5), (42.0, 21.0), (36.0, 21.0)])   # lane1
mm_arrow(ax, [(39.5, 50.5), (41.0, 50.5), (41.0, 24.0), (36.0, 24.0)])   # lane2
# exit from merge box to panel b
mm_arrow(ax, [(22.0, 13.5), (22.0, 8.2), (45.5, 8.2)])

# ---- (b) structure ensemble + fast triage funnel ----------------------------
mm_ax(ax, 57.0, 66.0, 16.0, 12.0)
mm_txt_w(ax, 65.0, 72.0, "apo + holo ensembles (multi-seed)", 15.0, size=4.5)
mm_arrow(ax, [(65.0, 66.0), (65.0, 62.0)])

mm_ax(ax, 56.0, 50.0, 18.0, 12.0)
mm_txt_w(ax, 65.0, 56.0, "Developability solubility · aggregation · expression",
         16.5, size=4.3)
mm_arrow(ax, [(65.0, 50.0), (65.0, 46.0)])

mm_ax(ax, 58.0, 34.0, 14.0, 12.0)
mm_txt_w(ax, 65.0, 40.0, "Fold & complex confidence pLDDT · PAE · contacts",
         12.5, size=4.3)
mm_arrow(ax, [(65.0, 34.0), (65.0, 30.0)])

mm_ax(ax, 60.0, 18.0, 10.0, 12.0)
mm_txt_w(ax, 65.0, 24.0, "Static interface / ΔΔG (Rosetta · FoldX · Pythia-PPI)",
         8.5, size=4.2)
# funnel exit: route down along the right corridor into the holo input box
mm_arrow(ax, [(65.0, 18.0), (65.0, 10.0), (92.25, 10.0), (92.25, 44.0)])

# ---- (c) FlexPPI dynamics-aware reranking (visual focus) --------------------
apo = mm_ax(ax, 86.0, 63.0, 12.5, 11.0)
mm_txt(ax, 92.25, 68.5, "apo\nstructure", size=4.8, color=TARGET, weight="bold")
holo = mm_ax(ax, 86.0, 44.0, 12.5, 11.0)
mm_txt(ax, 92.25, 49.5, "holo\ncomplex", size=4.8, color=TARGET, weight="bold")

flex = FancyBboxPatch(
    (103.0 / MM, 50.0 / MM), 21.0 / MM, 19.0 / MM,
    boxstyle=f"round,pad={0.8 / MM:.5f},rounding_size={1.2 / MM:.5f}",
    linewidth=1.4, edgecolor=FLEX, facecolor=FILL_PURPLE, zorder=3,
)
ax.add_patch(flex)
mm_txt(ax, 113.5, 65.8, "FlexPPI", size=7.5, weight="bold", color="#8E4FA0")
mm_txt_w(ax, 113.5, 61.2, "absolute states: RMSF · Sχ · SASA", 19.0, size=4.4)
mm_txt_w(ax, 113.5, 55.6, "signed partner-induced Δ + per-residue uncertainty",
         19.0, size=4.4)
mm_txt_w(ax, 113.5, 51.6, "direct vs. subtraction check", 19.0, size=4.2,
         color=GREY)

mm_arrow(ax, [(98.5, 68.5), (103.0, 68.5)])
mm_arrow(ax, [(98.5, 49.5), (103.0, 49.5)])

# per-residue output strips on the right of the FlexPPI box
mm_txt(ax, 130.0, 80.0, "per-residue", size=4.2, color=GREY)
strip_x, strip_w, strip_h = 127.0, 6.0, 2.6
import matplotlib.cm as cm
viri = [cm.viridis(i / 13.0) for i in range(14)]
ribbon(ax, strip_x, 66.0, strip_w, strip_h, viri)
mm_txt(ax, 130.0, 63.2, "RMSF", size=4.2)
ribbon(ax, strip_x, 58.0, strip_w, strip_h, viri)
mm_txt(ax, 130.0, 55.2, "Sχ", size=4.2)
ribbon(ax, strip_x, 50.0, strip_w, strip_h, viri)
mm_txt(ax, 130.0, 47.2, "SASA", size=4.2)

# zero-centred divergent bar for partner-induced change
div = ["#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF",
       "#E0E0E0", "#E0E0E0", "#F4C2B3", "#E08A72", "#D55E00", "#C23B22",
       "#A8321E", "#8B2518"]
ribbon(ax, strip_x, 40.0, strip_w, 3.0, div)
mm_txt(ax, 130.0, 36.0, "partner-\ninduced Δ", size=4.2, color=FEED)
# uncertainty strip
unc = ["#E5E7EB", "#D1D5DB", "#9CA3AF", "#6B7280", "#9CA3AF", "#D1D5DB",
       "#E5E7EB", "#B0B7BF", "#7B838E", "#4B5563", "#9CA3AF", "#D1D5DB",
       "#E5E7EB", "#C4C9D0"]
ribbon(ax, strip_x, 31.0, strip_w, 2.4, unc)
mm_txt(ax, 130.0, 27.8, "uncertainty", size=4.2, color=GREY)

# arrows from FlexPPI box (right edge) to output strips
mm_arrow(ax, [(124.0, 67.3), (127.0, 67.3)])
mm_arrow(ax, [(124.0, 59.3), (127.0, 59.3)])
mm_arrow(ax, [(124.0, 51.3), (127.0, 51.3)])
mm_arrow(ax, [(124.0, 50.0), (127.0, 41.5)])
mm_arrow(ax, [(124.0, 50.0), (127.0, 32.2)])

mm_txt(ax, 111.0, 24.0, "orthogonal dynamics signal · not ΔG",
       size=4.6, color=GREY, style="italic")

# ---- (d) Pareto selection + validation + feedback ---------------------------
pax = fig.add_axes([142.5 / MM, 58.5 / MM, 29.0 / MM, 24.0 / MM])
pax.set_xlim(0, 1)
pax.set_ylim(0, 1)
pax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for s in pax.spines.values():
    s.set_linewidth(0.5)
    s.set_color(GREY)
pax.set_xlabel("structure confidence", fontsize=4.6, color=INK, labelpad=1)
pax.set_ylabel("dynamics match", fontsize=4.6, color=INK, labelpad=1)
pts = [(0.30, 0.25), (0.42, 0.30), (0.55, 0.20), (0.62, 0.45), (0.38, 0.55),
       (0.48, 0.62), (0.70, 0.55), (0.58, 0.78), (0.78, 0.70), (0.88, 0.85)]
for i, (x, y) in enumerate(pts):
    c = PASS if i in (6, 7, 8, 9) else GREY
    m = "*" if i in (6, 7, 8, 9) else "o"
    pax.scatter(x, y, s=13, c=c, marker=m, edgecolors="none", zorder=3)
pax.plot([0.62, 0.70, 0.78, 0.88], [0.45, 0.55, 0.70, 0.85],
         color=PASS, lw=0.8, ls="--", zorder=2)
mm_txt(ax, 157.0, 83.8, "Pareto front", size=4.6, color=PASS, weight="bold")

mm_arrow(ax, [(157.0, 58.5), (157.0, 51.0)])

mm_ax(ax, 144.5, 39.0, 26.0, 12.0)
mm_txt_w(ax, 157.5, 45.0, "MD / expression / binding / specificity / function",
         24.5, size=4.6)

# single feedback arc: validation -> back to candidate generation
mm_arrow(ax, [(146.0, 39.0), (146.0, 4.0), (22.0, 4.0), (22.0, 13.5)],
         color=FEED, lw=1.0, style="--")
mm_txt(ax, 88.0, 2.4, "experimental feedback → ranking / next round "
       "(not an online update of FlexPPI)",
       size=4.4, color=FEED, style="italic")

fig.savefig(os.path.join(os.path.dirname(__file__), "fig_flexppi_design_pipeline.pdf"))
fig.savefig(os.path.join(os.path.dirname(__file__), "fig_flexppi_design_pipeline.svg"))
fig.savefig(os.path.join(os.path.dirname(__file__), "fig_flexppi_design_pipeline.png"), dpi=300)
print("wrote fig_flexppi_design_pipeline.{pdf,svg,png}")
