"""Build Figure 1 of the FlexPPI manuscript (170 mm x 152 mm).

Layout: a wide pipeline row (a) above three detail panels (b, c, d).
Run with D:\\Anaconda\\envs\\pvb_cpu\\python.exe build_figure1.py
"""
from pathlib import Path
import csv
import json

import numpy as np
from PIL import Image
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Circle

import figlib as F
from figlib import (canvas, text, rbox, rect, polyline, arrow, img, fit_img,
                    colorbar, read_rmsf_cmap,
                    INK, BODY, MUTED, FAINT, HAIR, RULE,
                    SCALAR, SCALAR_L, PAIR, PAIR_L, POINT, POINT_L,
                    ENC, ENC_L, READ, READ_L)

plt.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.default": "it",
})

HERE = Path(__file__).resolve().parent
A = HERE / "assets"
R = HERE / "render"
W, H = 170.0, 152.0

receipt = json.loads((R / "render_receipt.json").read_text(encoding="utf-8"))
cmap, (VLO, VHI) = read_rmsf_cmap(A / "all_840_rmsf_colors_original.csv")

curves = list(csv.DictReader(open(A / "curve_local_means.csv", encoding="utf-8")))
num = lambda key: np.array([float(r[key]) if r[key] else np.nan for r in curves])
rmsf_raw, rmsf_mean = num("rmsf_raw"), num("rmsf_local_mean")
hchi_raw, hchi_mean = num("hchi_raw"), num("hchi_local_mean")
x = np.arange(len(curves))

fig, ax = canvas(W, H)
fig.patch.set_facecolor("white")


def panel_letter(x0, y0, letter, title):
    text(ax, x0, y0, letter, 11.5, INK, "bold", va="center")
    text(ax, x0 + 5.0, y0, title, 9, INK, "bold", va="center")


def mmap(box, size):
    bx, by, bw, bh = box
    pw, ph = size
    return lambda px, py: (bx + px / pw * bw, by + (ph - py) / ph * bh)


def strip(x0, x1, y0, y1, raw, mean, vmax, col, raw_col, zbase=2):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="#FDFEFE",
                           edgecolor=RULE, lw=0.5, zorder=zbase))
    xa, xb = x0 + 0.6, x1 - 5.6
    ya, yb = y0 + 0.6, y1 - 0.6
    m = ~np.isnan(raw)
    ax.plot(x[m] / (len(x) - 1) * (xb - xa) + xa,
            ya + raw[m] / vmax * (yb - ya), color=raw_col, lw=0.3, zorder=zbase + 1)
    ax.plot(x / (len(x) - 1) * (xb - xa) + xa,
            ya + mean / vmax * (yb - ya), color=col, lw=0.8, zorder=zbase + 2)
    text(ax, x1 - 0.8, y1 - 1.6, f"{vmax:g}", 8, MUTED, ha="right", va="center")
    text(ax, x1 - 0.8, y0 + 1.6, "0", 8, MUTED, ha="right", va="center")


# =========================================================== panel a: flow ==
panel_letter(3, 146.0, "a", "FlexPPI overview")
STAGES = [("Input structure", 3, 34), ("Geometric encoding", 38, 68),
          ("Shared encoder", 72, 99), ("Joint readout", 103, 130),
          ("Predicted dynamics", 134, 167)]
for label, x0, x1 in STAGES:
    text(ax, (x0 + x1) / 2, 139.4, label, 8, BODY, "bold", ha="center")
for y, xa, xb in ((123.4, 34.2, 37.8), (116.2, 68.2, 71.8),
                  (120.4, 99.2, 102.8), (120.4, 130.2, 133.8)):
    arrow(ax, [(xa, y), (xb, y)], color="#8FA0AC", lw=0.9, head=1.15, z=6, wide=0.52)

# stage 1 ---------------------------------------------------------------------
fit_img(ax, R / "assembly_input.png", 3, 118.6, 31, 14.2)
text(ax, 3, 114.8, "4g9f complex", 8, BODY)
text(ax, 3, 110.4, "5 chains \u00b7 840 residues", 8, MUTED)
CHAIN_COLORS = receipt["chain_colors"]
for i, (chain, col) in enumerate(CHAIN_COLORS.items()):
    cx = 4.1 + i * 6.0
    ax.add_patch(Circle((cx, 105.2), 0.85, facecolor=col, edgecolor="none", zorder=5))
    text(ax, cx + 1.5, 105.2, chain, 8, BODY, va="center")
text(ax, 3, 99.6, "N, C\u03b1, C coordinates", 8, MUTED)
text(ax, 3, 92.0, "sequence", 8, MUTED)
for i, aa in enumerate("GSHSMR"):
    cx = 6.6 + i * 4.4
    ax.add_patch(Circle((cx, 96.4), 1.45, facecolor="none", edgecolor=FAINT, lw=0.55, zorder=5))
    text(ax, cx, 96.4, aa, 8, MUTED, ha="center", va="center", z=6)

# stage 2 ---------------------------------------------------------------------
def bars_glyph(x0, y0, w, h, n=9, color=SCALAR, z=4):
    bw = w / (n * 1.55)
    for i in range(n):
        g = 0.42 + 0.58 * (i % 3) / 2
        rect(ax, x0 + i * bw * 1.55, y0, bw, h * g, color, alpha=0.88, z=z)


def matrix_glyph(x0, y0, size, cells=7, color=PAIR, seed=3, z=4):
    rng = np.random.default_rng(seed)
    vals = rng.random((cells, cells))
    cs = size / cells
    c = matplotlib.colors.to_rgb(color)
    for i in range(cells):
        for j in range(cells):
            ax.add_patch(Rectangle((x0 + j * cs, y0 + (cells - 1 - i) * cs),
                                   cs * 0.92, cs * 0.92, facecolor=c,
                                   alpha=0.18 + 0.8 * vals[i, j], edgecolor="none", zorder=z))


def triad_glyph(x0, y0, s, z=5, lw=0.9):
    o = np.array([x0, y0])
    for ang, col in ((90, POINT), (212, SCALAR), (328, PAIR)):
        v = np.array([np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))])
        arrow(ax, [o, o + v * s], color=col, lw=lw, head=1.0, z=z, wide=0.5)


text(ax, 38, 127.6, "Residue $h_i$", 8, BODY)
bars_glyph(52.6, 123.8, 14.0, 6.4)
text(ax, 38, 112.4, "Pair $b_{ij}$", 8, BODY)
matrix_glyph(50.5, 102.8, 12.4)
text(ax, 38, 97.4, "Frames", 8, BODY)
triad_glyph(50.4, 95.4, 6.0)

# stage 3 ---------------------------------------------------------------------
for k, (dx, dy) in enumerate([(0.0, 0.0), (0.7, 0.75), (1.4, 1.5)][::-1]):
    rbox(ax, 72 + dx, 94.0 + dy, 27.0, 35.0, fill="#F4F7F9", edge=HAIR,
         lw=0.6, r=1.3, z=2 + k)
rbox(ax, 72, 94.0, 27.0, 35.0, fill="#FBFCFD", edge="#9FB0BC", lw=0.8, r=1.3, z=5)
text(ax, 72, 132.6, "\u00d7 4 blocks", 8, ENC, "bold")

text(ax, 85.5, 126.6, "LayerNorm", 8, BODY, ha="center")
rbox(ax, 73.8, 109.6, 23.4, 15.0, fill=POINT_L, edge="#C3B5D6", lw=0.7, r=1.0, z=6)
text(ax, 85.5, 121.8, "IPA \u00b7 8 heads", 8, POINT, "bold", ha="center")
nodes = [(77.4, 117.6), (81.8, 113.8), (89.4, 116.4), (93.6, 112.4)]
for i, j in [(0, 1), (1, 2), (2, 3), (0, 2), (1, 3)]:
    polyline(ax, [nodes[i], nodes[j]], color="#A995C0", lw=0.55, z=7)
for nx, ny in nodes:
    ax.add_patch(Circle((nx, ny), 0.8, facecolor="#FDFCFE", edgecolor=POINT, lw=0.6, zorder=8))
text(ax, 85.5, 111.2, "$q_k$   $b_{ij}$   $\\Vert p_i-p_j\\Vert$", 8, POINT, ha="center")

for yc in (106.6, 96.4):
    ax.add_patch(Circle((85.5, yc), 1.2, facecolor="white", edgecolor=ENC, lw=0.7, zorder=7))
    polyline(ax, [(84.3, yc), (86.7, yc)], color=ENC, lw=0.6, z=8)
    polyline(ax, [(85.5, yc - 1.2), (85.5, yc + 1.2)], color=ENC, lw=0.6, z=8)
# residual routes around the two sub-blocks of one encoder layer
polyline(ax, [(73.9, 94.6), (73.9, 106.6), (84.2, 106.6)], color="#A9B6C0", lw=0.55, z=6)
polyline(ax, [(73.9, 96.4), (84.2, 96.4)], color="#A9B6C0", lw=0.55, z=6)

text(ax, 85.5, 101.4, "Norm + FFN", 8, BODY, ha="center")

# stage 4 ---------------------------------------------------------------------
text(ax, 103, 129.4, "shared MLP", 8, READ, "bold")
for k, name in enumerate(["LayerNorm", "Linear", "SiLU", "Linear"]):
    y = 119.0 - k * 5.8
    rbox(ax, 103, y, 27.0, 4.2, fill=READ_L if k % 2 == 0 else "#F7FAF9",
         edge="#B9D0C7", lw=0.6, r=0.7, z=5)
    text(ax, 116.5, y + 2.1, name, 8, BODY, ha="center")
    if k < 3:
        polyline(ax, [(116.5, y), (116.5, y - 1.1)], color="#9FBDB2", lw=0.6, z=5)
        ax.add_patch(Polygon([(116.5, y - 1.6), (116.0, y - 0.8), (117.0, y - 0.8)],
                             closed=True, facecolor="#9FBDB2", edgecolor="none", zorder=6))
polyline(ax, [(116.5, 101.0), (116.5, 99.4)], color="#9FBDB2", lw=0.6, z=5)
polyline(ax, [(108.2, 99.4), (124.8, 99.4)], color="#9FBDB2", lw=0.6, z=5)
for xc, (name, col) in zip((108.2, 124.8), [("RMSF", SCALAR), ("H$_\\chi$", PAIR)]):
    polyline(ax, [(xc, 99.4), (xc, 98.2)], color=col, lw=0.6, z=5)
    ax.add_patch(Circle((xc, 97.0), 0.75, facecolor=col, edgecolor="none", zorder=6))
    text(ax, xc + 1.4, 97.0, name, 8, BODY, va="center")
text(ax, 116.5, 92.8, "$\\mu$, log $\\sigma^2$", 8, MUTED, ha="center")

# stage 5 ---------------------------------------------------------------------
fit_img(ax, R / "assembly_prediction.png", 134, 118.2, 33, 15.1)
text(ax, 134, 115.2, "RMSF (\u00c5)", 8, BODY)
colorbar(ax, 134, 111.4, 26.0, 2.2, cmap)
text(ax, 134, 109.6, "0", 8, MUTED, ha="left", va="center")
text(ax, 160.4, 109.6, "4", 8, MUTED, ha="left", va="center")
strip(134, 167, 100.6, 107.8, rmsf_raw, rmsf_mean, VHI, SCALAR, "#C8DEEB")
strip(134, 167, 89.0, 96.2, hchi_raw, hchi_mean, 4.1, PAIR, "#F0DEC1")
text(ax, 134, 98.9, "RMSF", 8, SCALAR, "bold")
text(ax, 134, 87.3, "H$_\\chi$ (nats)", 8, PAIR, "bold")
text(ax, 167, 87.3, "residue 1\u2013840", 8, MUTED, ha="right")
ax.plot([3, 167], [84.2, 84.2], color=RULE, lw=0.7, zorder=1)

# ==================================================== panel b: local frames ==
panel_letter(3, 79.8, "b", "Residue frames")
CH = receipt["chain_d"]
ch_size = Image.open(R / "chain_d.png").size
CH_BOX = (3, 46.4, 58.0, 58.0 * ch_size[1] / ch_size[0])
fit_img(ax, R / "chain_d.png", *CH_BOX[:2], CH_BOX[2], CH_BOX[3])
chmap = mmap(CH_BOX, ch_size)
for tag, col in (("i", POINT), ("j", SCALAR)):
    o = chmap(*CH["anchors"][f"{tag}_o"])
    for key in ("x", "y", "z"):
        tip = chmap(*CH["anchors"][f"{tag}_{key}"])
        arrow(ax, [o, tip], color={"x": SCALAR, "y": PAIR, "z": POINT}[key],
              lw=0.7, head=1.0, z=7, wide=0.5)
    text(ax, o[0] + 1.4, o[1] - 1.6, tag, 8.5, col, "bold", z=9)
text(ax, 3, 43.4, "chain D \u00b7 204 residues", 8, MUTED)

Z = receipt["zoom"]
z_size = Image.open(R / "frame_zoom.png").size
ZW = 30.4
Z_BOX = (2.6, 14.2, ZW, ZW * z_size[1] / z_size[0])
fit_img(ax, R / "frame_zoom.png", *Z_BOX[:2], Z_BOX[2], Z_BOX[3])
zmap = mmap(Z_BOX, z_size)
za = Z["anchors"]
o = zmap(*za["ca166"])
for key, col in (("x", SCALAR), ("y", PAIR), ("z", POINT)):
    tip = zmap(*za[key])
    v = (np.array(tip) - np.array(o)) * 0.72
    end = np.array(o) + v
    arrow(ax, [o, end], color=col, lw=0.85, head=1.15, z=8, wide=0.46)
    n = np.hypot(*v)
    text(ax, end[0] + v[0] / n * 1.9, end[1] + v[1] / n * 1.9, key, 8.5, col,
         "bold", ha="center", va="center", z=9, halo=True)
ax.add_patch(Circle(o, 0.95, facecolor="none", edgecolor=FAINT, lw=0.55,
                    ls=(0, (1.4, 1.2)), zorder=7))
# compact backbone legend under the magnification
for i, (lab, col) in enumerate((("N", "#3F6E9E"), ("C\u03b1", "#7D8A92"), ("C", "#B4BEC5"))):
    cx = 4.4 + i * 9.4
    ax.add_patch(Circle((cx, 11.2), 0.8, facecolor=col, edgecolor="none", zorder=6))
    text(ax, cx + 1.4, 11.2, lab, 8, MUTED, va="center")

text(ax, 35.4, 30.6, "$p_i = R_i\\,q_i + t_i$", 9, BODY)
text(ax, 35.4, 26.2, "local \u2192 global", 8, MUTED)

# ================================================ panel c: geometric attention
panel_letter(64, 79.8, "c", "Geometric attention")
SRC = [("Scalar", SCALAR, SCALAR_L), ("Pair", PAIR, PAIR_L), ("Points", POINT, POINT_L)]
for k, (name, col, fill) in enumerate(SRC):
    x0 = 64.0 + k * 15.2
    rbox(ax, x0, 69.0, 14.2, 6.4, fill=fill, edge=col, lw=0.6, r=0.9, z=5)
    text(ax, x0 + 7.1, 72.2, name, 8, col, "bold", ha="center")
    polyline(ax, [(x0 + 7.1, 69.0), (x0 + 7.1, 65.6)], color=col, lw=0.65, z=6)
polyline(ax, [(71.1, 65.6), (101.3, 65.6)], color=HAIR, lw=0.65, z=6)
arrow(ax, [(86.2, 65.6), (86.2, 63.6)], color=FAINT, lw=0.7, head=0.9)
rbox(ax, 64.0, 57.8, 44.4, 5.8, fill="#F3F7F5", edge="#AEC5BA", lw=0.6, r=0.9, z=5)
text(ax, 86.2, 60.7, "$w_{ij} = \\mathrm{softmax}_j\\,(\\cdot)$", 8, BODY, ha="center")

for k, (name, col) in enumerate([("scalar", SCALAR), ("pair", PAIR), ("points", POINT)]):
    y = 54.0 - k * 3.1
    polyline(ax, [(68.4, y), (81.4, y)], color=col, lw=1.0, z=5)
    text(ax, 67.6, y, name, 8, MUTED, ha="right")
    polyline(ax, [(81.4, y), (86.2, 47.0)], color=col, lw=0.5, z=5, ls=":")
ax.add_patch(Circle((86.2, 47.0), 0.8, facecolor="white", edgecolor=ENC, lw=0.6, zorder=6))
arrow(ax, [(86.2, 46.2), (86.2, 43.8)], color=FAINT, lw=0.7, head=0.9)
rbox(ax, 64.0, 37.8, 44.4, 5.6, fill="white", edge="#A994BF", lw=0.6, r=0.9, z=5)
text(ax, 86.2, 40.6, "rotate back to frame $i$", 8, BODY, ha="center")
arrow(ax, [(86.2, 37.8), (86.2, 35.4)], color=FAINT, lw=0.7, head=0.9)
rbox(ax, 64.0, 29.6, 44.4, 5.6, fill="#EDF3F8", edge="#AAC1D2", lw=0.6, r=0.9, z=5)
text(ax, 86.2, 32.4, "merge + linear", 8, BODY, ha="center")
arrow(ax, [(86.2, 29.6), (86.2, 27.2)], color=FAINT, lw=0.7, head=0.9)
rbox(ax, 70.0, 21.4, 32.4, 5.6, fill="white", edge=ENC, lw=0.7, r=0.9, z=5)
text(ax, 86.2, 24.2, "$h_i^{\\prime}$", 9, ENC, "bold", ha="center")
text(ax, 86.2, 17.8, "updated residue representation", 8, MUTED, ha="center")

# ================================================== panel d: MD supervision
panel_letter(116, 79.8, "d", "MD supervision")
SNAP = ["d_snapshot_first.png", "d_snapshot_middle.png", "d_snapshot_last.png"]
for k, (name, t) in enumerate(zip(SNAP, ["0 ns", "100 ns", "200 ns"])):
    x0 = 116.0 + k * 17.4
    fit_img(ax, A / name, x0, 66.4, 15.6, 6.4)
    text(ax, x0 + 7.8, 63.8, t, 8, MUTED, ha="center")
    polyline(ax, [(x0 + 7.8, 62.2), (x0 + 7.8, 59.8)], color=HAIR, lw=0.6, z=4)
polyline(ax, [(123.8, 59.8), (157.2, 59.8)], color=HAIR, lw=0.6, z=4)
arrow(ax, [(140.5, 59.8), (140.5, 57.4)], color=FAINT, lw=0.7, head=0.9)

for k, (name, col) in enumerate([("RMSF", SCALAR), ("H$\\chi$", PAIR)]):
    x0 = 116.0 + k * 26.2
    rbox(ax, x0, 47.0, 24.4, 7.6, fill="white", edge=col, lw=0.6, r=0.9, z=5)
    ax.add_patch(Circle((x0 + 2.5, 52.4), 0.85, facecolor=col, edgecolor="none", zorder=6))
    text(ax, x0 + 4.6, 52.4, name, 8, BODY, va="center")
    text(ax, x0 + 2.5, 49.6, "mean, variance", 8, MUTED, va="center")
polyline(ax, [(128.2, 47.0), (128.2, 44.8), (140.5, 44.8)], color=HAIR, lw=0.6, z=4)
polyline(ax, [(154.4, 47.0), (154.4, 44.8), (140.5, 44.8)], color=HAIR, lw=0.6, z=4)
arrow(ax, [(140.5, 44.8), (140.5, 42.4)], color=FAINT, lw=0.7, head=0.9)
text(ax, 140.5, 40.2, "MD mean $\\pm$ SD", 8, MUTED, ha="center")

rbox(ax, 129.4, 17.6, 37.6, 19.6, fill="#F2F7F4", edge="#A7C1B4", lw=0.7, r=1.2, z=3)
text(ax, 148.2, 33.2, "Gaussian NLL", 8, READ, "bold", ha="center")
text(ax, 148.2, 27.4, "$(y-\\mu)^2/2v + \\log(v)/2$", 8, BODY, ha="center")
text(ax, 148.2, 21.0, "$v = \\sigma^2 + v_{\\mathrm{MD}}$", 8, READ, ha="center")
arrow(ax, [(140.5, 37.2), (140.5, 35.0)], color=FAINT, lw=0.7, head=0.9)
text(ax, 116.0, 26.8, "model", 8, MUTED)
text(ax, 116.0, 23.0, "$\\mu$, $\\sigma^2$", 8, MUTED)
arrow(ax, [(123.6, 24.8), (128.0, 24.8)], color=FAINT, lw=0.7, head=0.9)

fig.savefig(HERE / "figure1_preview.png", dpi=300)

# ------------------------------------------------------------ final output ---
for suffix, kwargs in (("pdf", {"dpi": 600}), ("svg", {"dpi": 600}), ("png", {"dpi": 600})):
    fig.savefig(HERE / f"fig1_framework.{suffix}", facecolor="white", **kwargs)

sizes = sorted({t.get_fontsize() for t in ax.texts})
report = {
    "canvas_mm": [W, H],
    "figures": {s: str(HERE / f"fig1_framework.{s}") for s in ("pdf", "svg", "png")},
    "min_font_pt": min(sizes),
    "font_sizes_pt": sizes,
    "n_text": len(ax.texts),
    "render_assets": {k: receipt["sha256"][k] for k in ("input", "prediction", "chain_d", "zoom")},
    "chain_colors": receipt["chain_colors"],
    "colour_map": "rebuilt from the archived per-residue RMSF colours",
}
(HERE / "figure_audit.json").write_text(json.dumps(report, indent=2, ensure_ascii=False),
                                        encoding="utf-8")
print(json.dumps({"min_font": report["min_font_pt"], "sizes": sizes}))










