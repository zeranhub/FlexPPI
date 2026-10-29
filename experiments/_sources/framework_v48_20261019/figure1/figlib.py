"""Figure 1 drawing library: millimetre canvas, flat vector primitives.

Everything in the figure is positioned in millimetres so that the exported
master is exactly the insertion width and font sizes are the printed sizes.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle, Circle
from PIL import Image

MM = 1 / 25.4

# ------------------------------------------------------------------ palette --
INK = "#1F2C35"
BODY = "#35454F"
MUTED = "#6C7C88"
FAINT = "#93A2AD"
HAIR = "#C6D1D8"
RULE = "#DDE4E9"
PAPER = "#FFFFFF"

SCALAR, SCALAR_L = "#3E7CA6", "#E9F2F8"
PAIR, PAIR_L = "#AF8438", "#F8F0DE"
POINT, POINT_L = "#7A6595", "#EFEAF4"
ENC, ENC_L = "#566F7E", "#EEF3F6"
READ, READ_L = "#4C8173", "#EAF2EF"
DYNA, DYNA_L = "#4C8173", "#EAF2EF"

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 0.6,
})


def canvas(width_mm, height_mm):
    fig = plt.figure(figsize=(width_mm * MM, height_mm * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width_mm)
    ax.set_ylim(0, height_mm)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def text(ax, x, y, s, size=8, color=BODY, weight="normal", ha="left", va="center",
         style="normal", family="Arial", z=6, rotation=0, lsp=1.25, halo=False):
    t = ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha,
                va=va, fontstyle=style, family=family, zorder=z,
                rotation=rotation, linespacing=lsp)
    if halo:
        import matplotlib.patheffects as pe
        t.set_path_effects([pe.withStroke(linewidth=1.6, foreground="white")])
    return t


def rbox(ax, x, y, w, h, fill="none", edge=HAIR, lw=0.7, r=1.1, z=3, ls="-", alpha=1.0):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       facecolor=fill, edgecolor=edge, linewidth=lw, zorder=z,
                       linestyle=ls, alpha=alpha, mutation_aspect=1)
    ax.add_patch(p)
    return p


def rect(ax, x, y, w, h, fill, edge="none", lw=0.0, z=3, alpha=1.0):
    p = Rectangle((x, y), w, h, facecolor=fill, edgecolor=edge,
                  linewidth=lw, zorder=z, alpha=alpha)
    ax.add_patch(p)
    return p


def polyline(ax, pts, color=HAIR, lw=0.7, z=4, ls="-", cap="butt", dash=None):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ln, = ax.plot(xs, ys, color=color, lw=lw, zorder=z, solid_capstyle=cap,
                  linestyle=ls, dash_capstyle="butt")
    if dash:
        ln.set_dashes(dash)
    return ln


def arrow(ax, pts, color=FAINT, lw=0.8, head=1.1, z=4, ls="-", wide=0.55):
    """Polyline with a small solid triangular head at the last point."""
    pts = [np.asarray(p, float) for p in pts]
    end, prev = pts[-1], pts[-2]
    v = end - prev
    L = np.linalg.norm(v)
    if L == 0:
        return
    u = v / L
    base = end - u * head
    n = np.array([-u[1], u[0]])
    polyline(ax, pts[:-1] + [base], color=color, lw=lw, z=z, ls=ls)
    ax.add_patch(Polygon([end, base + n * head * wide, base - n * head * wide],
                         closed=True, facecolor=color, edgecolor="none", zorder=z + .1))


def img(ax, path, x, y, w=None, h=None, z=2, alpha=1.0):
    """Place an RGBA raster inside the box, preserving its aspect ratio."""
    im = Image.open(path).convert("RGBA")
    arr = np.asarray(im)
    iw, ih = im.size
    if w is None:
        w = h * iw / ih
    if h is None:
        h = w * ih / iw
    ax.imshow(arr, extent=[x, x + w, y, y + h], zorder=z, alpha=alpha,
              interpolation="lanczos", aspect="auto")
    return w, h


def fit_img(ax, path, x, y, w, h, z=2, anchor="center", alpha=1.0):
    """Place an image inside a box, scaled to fit, keeping its aspect ratio."""
    im = Image.open(path)
    iw, ih = im.size
    s = min(w / iw, h / ih)
    ww, hh = iw * s, ih * s
    if anchor == "center":
        cx, cy = x + (w - ww) / 2, y + (h - hh) / 2
    elif anchor == "left":
        cx, cy = x, y + (h - hh) / 2
    else:
        cx, cy = x + (w - ww), y + (h - hh) / 2
    return img(ax, path, cx, cy, ww, hh, z=z, alpha=alpha)


def colorbar(ax, x, y, w, h, cmap, z=3, n=256):
    grad = np.linspace(0, 1, n)[None, :]
    ax.imshow(grad, extent=[x, x + w, y, y + h], cmap=cmap, aspect="auto",
              zorder=z, origin="lower", interpolation="bilinear")


def dashed_pair(ax, p0, p1, color=FAINT, lw=0.6, z=2, dash=(2.2, 1.8)):
    ln, = ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, zorder=z)
    ln.set_dashes(dash)


# --------------------------------------------------------------- colormaps --
def read_rmsf_cmap(csv_path):
    """Rebuild the RMSF colour ramp from the archived per-residue colours."""
    import csv as _csv
    rows = list(_csv.DictReader(open(csv_path, encoding="utf-8")))
    vals = np.array([float(r["rmsf_prediction_A"]) for r in rows])
    rgb = np.array([[float(r["red"]), float(r["green"]), float(r["blue"])] for r in rows])
    lo, hi = 0.0, 4.0
    xs = np.clip(vals / hi, 0, 1)
    order = np.argsort(xs)
    xs, rgb = xs[order], rgb[order]
    # average duplicate abscissae, then sample anchors on an even grid
    uniq, idx = np.unique(np.round(xs, 4), return_index=True)
    anchors = []
    for frac in np.linspace(0, 1, 15):
        near = np.argmin(np.abs(uniq - frac))
        anchors.append(tuple(np.clip(rgb[idx[near]], 0, 1)))
    anchors[0] = (0.176, 0.396, 0.553)
    anchors[-1] = (0.607, 0.235, 0.188)
    return LinearSegmentedColormap.from_list("rmsf", anchors), (lo, hi)
