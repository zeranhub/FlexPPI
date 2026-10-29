# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) True-coordinate molecular views -> manuscript PyMOL renders -> param inherit
# (b) Absolute residue profiles -> LineTrend/plot_trend.py -> param inherit
# (c) Absolute paired residue values -> PairedBoxScatter/plot_PairedBoxScatter1.R -> param inherit
# (d) Distinct corroboration/limitation -> case-specific complete data -> param inherit
# RULE: no native-run panels; no observation may be selected by prediction quality.

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
import json
import hashlib
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text

# The manuscript uses 9 pt at the exported physical size; baseline is retained above.
mpl.rcParams.update({key: 9 for key in (
    "font.size", "axes.titlesize", "axes.labelsize", "xtick.labelsize",
    "ytick.labelsize", "legend.fontsize", "figure.titlesize"
)})
FIRST = "#CF8C27"
SECOND = "#762A83"
HLA = "#297B82"
GRAY = "#6B6B6B"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_canvas(height_mm=185):
    assert 175 <= height_mm <= 205
    return plt.figure(figsize=(183 / 25.4, height_mm / 25.4), facecolor="white")


def heading(fig, letter, title, y):
    fig.text(.018, y, letter, fontsize=11, weight="bold", va="baseline")
    fig.text(.061, y, title, fontsize=9.3, va="baseline")


def legend(fig, labels, anchor=(.968, .704), ncol=2):
    handles = [
        Line2D([], [], color=FIRST, marker="o", mfc="white", ms=4, lw=1, label=labels[0]),
        Line2D([], [], color=SECOND, marker="s", mfc=SECOND, ms=3.7, lw=1, label=labels[1]),
    ]
    return fig.legend(handles=handles, loc="upper right", bbox_to_anchor=anchor,
                      ncol=ncol, columnspacing=1.3, handlelength=1.4,
                      handletextpad=.45, borderaxespad=0)


def cropped_image(path, pad=12, threshold=238):
    im = Image.open(path).convert("RGB")
    y, x = np.where(np.any(np.asarray(im) < threshold, axis=2))
    assert len(x), "Molecular panel cannot be a blank image."
    box = (max(0, int(x.min()) - pad), max(0, int(y.min()) - pad),
           min(im.width, int(x.max()) + pad), min(im.height, int(y.max()) + pad))
    return im.crop(box), box


def add_molecule(fig, path, rect, title=None, color=GRAY, crop=True):
    image, cropbox = cropped_image(path) if crop else (Image.open(path).convert("RGB"), None)
    ax = fig.add_axes(rect)
    ax.imshow(image)
    ax.axis("off")
    if title:
        ax.set_title(title, color=color, fontsize=9, pad=3)
    return ax, {"source": str(path), "SHA256": sha256(path), "cropbox": cropbox}


def absolute_profile(ax, frame, x_col, a_col, b_col, ylabel, label_col=None,
                     segment_col=None, highlights=(), ylim=None, all_ticks=True):
    """Plot every row; optional gaps denote predefined regions, never missing chains.

    `highlights` contains x identifiers rather than row indices. The caller supplies
    valid data and records invalid positions separately before invoking this function.
    """
    assert len(frame) and not frame[x_col].isna().any()
    assert np.isfinite(frame[[a_col, b_col]].to_numpy(float)).all()
    assert (frame[[a_col, b_col]].to_numpy(float) >= 0).all()
    groups = list(frame.groupby(segment_col, sort=False)) if segment_col else [(None, frame)]
    plot_rows, tick_x, tick_labels = [], [], []
    xoffset = 0
    data_values = frame[[a_col, b_col]].to_numpy(float)
    ymax = float(data_values.max())
    ytop = ylim[1] if ylim else max(.5, np.ceil(ymax * 1.25 * 5) / 5)
    for group, data in groups:
        x = np.arange(len(data), dtype=float) + xoffset
        for values, color, marker, face in [
            (data[a_col], FIRST, "o", "white"),
            (data[b_col], SECOND, "s", SECOND),
        ]:
            ax.plot(x, values, color=color, marker=marker, ms=3.3, mfc=face,
                    mec=color, lw=1.05, zorder=3)
        if group is not None:
            ax.text(x.mean(), 1.055, str(group), transform=ax.get_xaxis_transform(),
                    ha="center", fontsize=9, color=GRAY)
            ax.plot([x[0]-.2, x[-1]+.2], [1.015, 1.015],
                    transform=ax.get_xaxis_transform(), lw=.6, color="#B5B5B5", clip_on=False)
        for j, (_, row) in enumerate(data.iterrows()):
            xpos = float(x[j])
            if row[x_col] in highlights:
                ax.axvspan(xpos-.4, xpos+.4, color=SECOND, alpha=.06, lw=0)
            label = str(row[label_col]) if label_col else str(row[x_col])
            if all_ticks or j in (0, len(data)-1):
                tick_x.append(xpos)
                tick_labels.append(label)
            plot_rows.append({"plot_x": xpos, "id": str(row[x_col]), "label": label,
                              "segment": str(group), "first": float(row[a_col]),
                              "second": float(row[b_col])})
        xoffset = x[-1] + (2.5 if segment_col else 1)
    ax.set(xlim=(-.6, plot_rows[-1]["plot_x"] + .6), ylim=ylim or (0, ytop),
           ylabel=ylabel, xticks=tick_x)
    ax.set_xticklabels(tick_labels)
    ax.tick_params(length=2.5, pad=3)
    assert ax.get_ylim()[0] <= data_values.min() and ax.get_ylim()[1] >= data_values.max()
    return plot_rows


def paired_residue_values(ax, frame, a_col, b_col, labels, ylabel=r"$H_\chi$ (nats)",
                          highlight_ids=(), id_col=None, mutation_col=None,
                          ylim=None, rotate=55):
    """Absolute pairs; grey connectors join two states of the same site only."""
    assert len(frame) == len(labels)
    values = frame[[a_col, b_col]].to_numpy(float)
    assert np.isfinite(values).all() and (values >= 0).all()
    records = []
    for j, (_, row) in enumerate(frame.iterrows()):
        identifier = row[id_col] if id_col else j
        if identifier in highlight_ids:
            ax.axvspan(j-.4, j+.4, color=HLA, alpha=.07, lw=0)
        ax.plot([j-.14, j+.14], [row[a_col], row[b_col]], color="#AAAAAA", lw=.7, zorder=1)
        ax.plot(j-.14, row[a_col], "o", mfc="white", mec=FIRST, ms=3.8, zorder=3)
        ax.plot(j+.14, row[b_col], "s", mfc=SECOND, mec=SECOND, ms=3.3, zorder=3)
        altered = bool(row[mutation_col]) if mutation_col else False
        if altered:
            ax.plot(j, max(row[a_col], row[b_col]) + .08, marker="v", color=GRAY,
                    ms=3.5, clip_on=False)
        records.append({"plot_x": j, "id": str(identifier), "label": labels[j],
                        "first": float(row[a_col]), "second": float(row[b_col]),
                        "altered_chemistry": altered})
    ymax = max(.5, np.ceil(values.max() * 1.2 * 5) / 5)
    ax.set(xlim=(-.6, len(frame)-.4), ylim=ylim or (0, ymax), ylabel=ylabel,
           xticks=np.arange(len(frame)))
    ax.set_xticklabels(labels, rotation=rotate, ha="right", rotation_mode="anchor")
    ax.tick_params(length=2.5, pad=3)
    assert ax.get_ylim()[0] <= values.min() and ax.get_ylim()[1] >= values.max()
    return records


def delta_comparison(ax, observed, predicted, labels=None, unit="Å", title=None):
    """Independent corroboration panel; a bad relationship is retained unchanged."""
    observed, predicted = np.asarray(observed, float), np.asarray(predicted, float)
    assert len(observed) == len(predicted) and len(observed)
    assert np.isfinite(observed).all() and np.isfinite(predicted).all()
    extent = max(.05, np.max(np.abs(np.r_[observed, predicted])) * 1.18)
    ax.axhline(0, color="#D0D0D0", lw=.5, zorder=0)
    ax.axvline(0, color="#D0D0D0", lw=.5, zorder=0)
    ax.plot([-extent, extent], [-extent, extent], color="#BDBDBD", lw=.6, ls="--", zorder=0)
    ax.scatter(observed, predicted, c=HLA, s=16, alpha=.85, linewidths=0, zorder=3)
    ax.set(xlim=(-extent, extent), ylim=(-extent, extent),
           xlabel=f"Existing MD change ({unit})", ylabel=f"Model change ({unit})")
    # Locator may create off-view ticks. Make the rendered tick set explicit so
    # the audit reports only labels the axes actually show.
    for ticks, setter in [(ax.get_xticks(), ax.set_xticks), (ax.get_yticks(), ax.set_yticks)]:
        setter([tick for tick in ticks if -extent <= tick <= extent])
    if title:
        ax.set_title(title, fontsize=9, pad=5)
    if labels:
        for x, y, text in zip(observed, predicted, labels):
            ax.annotate(text, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=9)
    return {"n_positions": int(len(observed)), "MAE": float(np.mean(np.abs(predicted-observed))),
            "direction_agreement": float(np.mean(np.sign(observed) == np.sign(predicted))),
            "n_definition": "Aligned residue positions; not independent experimental replicas",
            "tests": None, "regression_fit": None, "unit": unit}


def export_audited(fig, stem, extra=None, dpi=450):
    """Export the exact canvas; return an audit bound to the rendered PNG SHA256."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = []
    hidden_axis_text_ids = set()
    for ax in fig.axes:
        if not ax.axison:
            hidden_axis_text_ids.update(id(t) for t in [
                *ax.get_xticklabels(), *ax.get_yticklabels(),
                ax.xaxis.label, ax.yaxis.label, ax.xaxis.offsetText, ax.yaxis.offsetText,
            ])
    for artist in fig.findobj(match=Text):
        if id(artist) in hidden_axis_text_ids or not artist.get_visible() or not artist.get_text():
            continue
        box = artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        texts.append({"text": artist.get_text(), "font_pt": float(artist.get_fontsize()),
                      "bounds": [float(v) for v in box.extents],
                      "inside_canvas": bool(box.x0 >= -.001 and box.y0 >= -.001 and
                                             box.x1 <= 1.001 and box.y1 <= 1.001)})
    assert texts and min(t["font_pt"] for t in texts) >= 9
    assert all(t["inside_canvas"] for t in texts), [t for t in texts if not t["inside_canvas"]]
    assert not any("\u4e00" <= c <= "\u9fff" for t in texts for c in t["text"])
    for ext in ("pdf", "svg", "png"):
        with mpl.rc_context({"savefig.bbox": None}):
            fig.savefig(str(stem) + "." + ext, bbox_inches=None, dpi=dpi)
    audit = {"status": "code_checks_passed_visual_review_pending",
             "dimensions_mm": [float(v * 25.4) for v in fig.get_size_inches()],
             "dpi": dpi, "minimum_font_pt": min(t["font_pt"] for t in texts),
             "all_visible_text_English": True, "all_visible_text_inside_canvas": True,
             "visible_text": texts, "PNG_SHA256": sha256(str(stem) + ".png")}
    if extra:
        audit.update(extra)
    stem.with_suffix(".audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf8")
    return audit
