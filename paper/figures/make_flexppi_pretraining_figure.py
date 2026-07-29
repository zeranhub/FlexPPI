# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) label-construction schematic -> cross-type inherit (PaperGallery grouping) -> param inherit
# (b) single-state encoder schematic -> cross-type inherit (PaperGallery central backbone) -> param inherit
# (c) probabilistic output schematic -> cross-type inherit (PaperGallery split supervision) -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

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
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]  # QA token: DIVERGING = [
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

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


mm_to_inch = 1 / 25.4

BLUE = CATEGORICAL[0]
RED = CATEGORICAL[1]
GREEN = CATEGORICAL[2]
ORANGE = CATEGORICAL[3]
PURPLE = CATEGORICAL[4]
DARK_GREY = "#555555"
MID_GREY = "#C8CDD4"
PALE_BLUE = "#EEF4FA"
PALE_GREEN = "#EEF6F0"
PALE_ORANGE = "#FFF6E8"
PALE_PURPLE = "#F5EFF7"
PALE_RED = "#FAEEEE"
PALE_GREY = "#F5F6F7"
WHITE = "#FFFFFF"


def find_cjk_font(size, bold=False):
    candidates = [
        Path(r"C:\Windows\Fonts\msyh.ttc"),
        Path(r"C:\Windows\Fonts\msyhbd.ttc"),
        Path(r"C:\Windows\Fonts\simhei.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            prop = FontProperties(fname=str(candidate), size=size)
            prop.set_weight("bold" if bold else "normal")
            return prop
    return FontProperties(family="sans-serif", size=size, weight="bold" if bold else "normal")


CJK_TITLE = find_cjk_font(8.8, bold=True)
CJK_PANEL = find_cjk_font(7.4, bold=True)
CJK_HEAD = find_cjk_font(7.0, bold=True)
CJK_TEXT = find_cjk_font(6.1)
CJK_SMALL = find_cjk_font(5.7)
CJK_TINY = find_cjk_font(5.3)


def rounded_box(
    ax,
    x,
    y,
    w,
    h,
    text="",
    edge=MID_GREY,
    fill=WHITE,
    font=CJK_TEXT,
    lw=0.75,
    radius=0.012,
    text_color=BLACK,
    zorder=2,
):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.005,rounding_size={radius}",
        linewidth=lw,
        edgecolor=edge,
        facecolor=fill,
        zorder=zorder,
    )
    ax.add_patch(patch)
    if text:
        ax.text(
            x + w / 2,
            y + h / 2,
            text,
            ha="center",
            va="center",
            color=text_color,
            fontproperties=font,
            linespacing=1.2,
            zorder=zorder + 1,
        )
    return patch


def arrow(
    ax,
    start,
    end,
    color=DARK_GREY,
    lw=0.85,
    linestyle="-",
    connectionstyle="arc3,rad=0",
    zorder=4,
):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=8,
        linewidth=lw,
        color=color,
        linestyle=linestyle,
        connectionstyle=connectionstyle,
        shrinkA=1.2,
        shrinkB=1.2,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def panel_background(ax, x, y, w, h, label, title, fill):
    rounded_box(ax, x, y, w, h, edge="#D2D6DB", fill=fill, lw=0.7, radius=0.018, zorder=0)
    ax.text(
        x + 0.014,
        y + h - 0.028,
        label,
        ha="left",
        va="top",
        color=BLACK,
        fontsize=9,
        fontweight="bold",
        zorder=6,
    )
    ax.text(
        x + 0.048,
        y + h - 0.029,
        title,
        ha="left",
        va="top",
        color=BLACK,
        fontproperties=CJK_PANEL,
        zorder=6,
    )


def chip(ax, x, y, w, text, color, fill, font=CJK_TINY):
    rounded_box(ax, x, y, w, 0.034, text, edge=color, fill=fill, font=font, lw=0.65, radius=0.010)


def draw_traces(ax, x, y, w, h):
    xs = np.linspace(0, 1, 70)
    colors = [BLUE, "#4393C3", "#74ADD1"]
    for idx, (phase, color) in enumerate(zip([0.0, 0.7, 1.35], colors)):
        center = y + h * (0.78 - idx * 0.29)
        ys = center + h * 0.065 * (
            np.sin(xs * 4.2 * np.pi + phase) + 0.32 * np.sin(xs * 11.0 + phase)
        )
        ax.plot(x + xs * w, ys, color=color, lw=0.9, solid_capstyle="round", zorder=4)
        ax.text(
            x - 0.006,
            center,
            f"R{idx + 1}",
            ha="right",
            va="center",
            color=color,
            fontsize=5.0,
            zorder=5,
        )


def draw_windows(ax, x, y, w, h):
    xs = np.linspace(0, 1, 100)
    ys = y + h * (0.50 + 0.18 * np.sin(xs * 5.0 * np.pi) + 0.08 * np.sin(xs * 17.0))
    ax.plot(x + xs * w, ys, color=DARK_GREY, lw=0.75, zorder=3)
    window_colors = [BLUE, ORANGE, PURPLE]
    for idx in range(3):
        left = x + idx * w / 3
        right = x + (idx + 1) * w / 3
        ax.add_patch(
            Rectangle(
                (left, y + 0.002),
                w / 3,
                h - 0.004,
                facecolor=window_colors[idx],
                edgecolor="none",
                alpha=0.08,
                zorder=1,
            )
        )
        if idx:
            ax.plot([left, left], [y, y + h], color=MID_GREY, lw=0.55, ls=(0, (2, 2)), zorder=2)
        ax.text(
            (left + right) / 2,
            y + h + 0.007,
            f"$W_{idx + 1}$",
            ha="center",
            va="bottom",
            fontsize=5.0,
            color=window_colors[idx],
            zorder=5,
        )


def draw_label_tensor(ax, x, y, w, h):
    ax.text(
        x + w / 2,
        y + h - 0.024,
        "逐残基监督张量",
        ha="center",
        va="center",
        fontproperties=CJK_HEAD,
        color=BLACK,
        zorder=6,
    )
    grid_x = x + 0.024
    grid_y = y + 0.040
    cell_w = 0.018
    cell_h = 0.018
    col_colors = [BLUE, PURPLE, GREEN]
    for row in range(5):
        for col in range(3):
            face = col_colors[col] if row not in {1, 4} or col != 2 else WHITE
            alpha = 0.72 if face != WHITE else 1.0
            ax.add_patch(
                Rectangle(
                    (grid_x + col * (cell_w + 0.004), grid_y + row * (cell_h + 0.004)),
                    cell_w,
                    cell_h,
                    facecolor=face,
                    edgecolor=col_colors[col],
                    linewidth=0.55,
                    alpha=alpha,
                    zorder=4,
                )
            )
    labels = ["RMSF", "$H_\\chi$", "SASA"]
    for col, label in enumerate(labels):
        ax.text(
            grid_x + col * (cell_w + 0.004) + cell_w / 2,
            grid_y - 0.010,
            label,
            ha="center",
            va="top",
            fontsize=5.0,
            color=col_colors[col],
            zorder=5,
        )
    ax.text(
        x + 0.112,
        y + 0.102,
        "$\\mathbf{Y}$  标签均值",
        ha="left",
        va="center",
        fontproperties=CJK_SMALL,
        color=BLACK,
        zorder=5,
    )
    ax.text(
        x + 0.112,
        y + 0.070,
        "$\\mathbf{S}^{2}$  观测方差",
        ha="left",
        va="center",
        fontproperties=CJK_SMALL,
        color=BLACK,
        zorder=5,
    )
    ax.text(
        x + 0.112,
        y + 0.038,
        "$\\mathbf{M}$  有效掩码",
        ha="left",
        va="center",
        fontproperties=CJK_SMALL,
        color=BLACK,
        zorder=5,
    )


def draw_token_strip(ax, x, y, w, colors, labels=None, height=0.025):
    count = len(colors)
    gap = 0.004
    cell_w = (w - gap * (count - 1)) / count
    for idx, color in enumerate(colors):
        rounded_box(
            ax,
            x + idx * (cell_w + gap),
            y,
            cell_w,
            height,
            edge=color,
            fill=color,
            lw=0.55,
            radius=0.005,
            zorder=4,
        )
        if labels:
            ax.text(
                x + idx * (cell_w + gap) + cell_w / 2,
                y + height / 2,
                labels[idx],
                ha="center",
                va="center",
                fontsize=5.3,
                color=WHITE,
                fontweight="bold",
                zorder=6,
            )


def draw_backbone_icon(ax, x, y, w, h):
    points = np.asarray(
        [
            [0.05, 0.42],
            [0.22, 0.69],
            [0.40, 0.48],
            [0.58, 0.78],
            [0.76, 0.52],
            [0.93, 0.68],
        ]
    )
    px = x + points[:, 0] * w
    py = y + points[:, 1] * h
    ax.plot(px, py, color=GREEN, lw=1.15, zorder=4)
    for idx, (cx, cy) in enumerate(zip(px, py)):
        ax.add_patch(
            Circle(
                (cx, cy),
                radius=min(w, h) * 0.045,
                facecolor=WHITE,
                edgecolor=GREEN,
                linewidth=0.85,
                zorder=5,
            )
        )
        if idx in {1, 4}:
            ax.plot([cx, cx + w * 0.09], [cy, cy + h * 0.16], color=ORANGE, lw=0.7, zorder=4)
            ax.plot([cx, cx - w * 0.07], [cy, cy + h * 0.12], color=BLUE, lw=0.7, zorder=4)
    arrow(
        ax,
        (px[2], py[2]),
        (px[4], py[4]),
        color=GREEN,
        lw=0.65,
        linestyle=(0, (2, 2)),
        connectionstyle="arc3,rad=-0.34",
        zorder=3,
    )


def draw_ipa_stack(ax, x, y, w, h):
    rounded_box(ax, x, y, w, h, edge=PURPLE, fill=PALE_PURPLE, lw=1.05, radius=0.020, zorder=2)
    ax.text(
        x + w / 2 - 0.008,
        y + h - 0.034,
        "固定框架 IPA",
        ha="center",
        va="center",
        fontproperties=CJK_HEAD,
        color=PURPLE,
        zorder=6,
    )
    ax.text(
        x + w - 0.024,
        y + h - 0.034,
        "× 4",
        ha="center",
        va="center",
        fontsize=6.0,
        color=PURPLE,
        fontweight="bold",
        zorder=6,
    )
    ax.text(
        x + w / 2,
        y + h - 0.070,
        "8 头几何自注意力 · 固定局部框架 $(\\mathbf{R}_i,\\mathbf{t}_i)$",
        ha="center",
        va="center",
        fontproperties=CJK_TINY,
        color=DARK_GREY,
        zorder=6,
    )
    row_y = [y + 0.121, y + 0.083, y + 0.045]
    row_text = [
        "标量通道  $\\mathbf{q}_i^\\top\\mathbf{k}_j$",
        "残基对通道  $b(\\mathbf{p}_{ij})$",
        "三维点通道  $\\|\\mathbf{q}^{pt}-\\mathbf{k}^{pt}\\|^2$",
    ]
    row_colors = [BLUE, GREEN, ORANGE]
    for yy, text, color in zip(row_y, row_text, row_colors):
        rounded_box(
            ax,
            x + 0.028,
            yy,
            w - 0.056,
            0.028,
            text,
            edge=color,
            fill=WHITE,
            font=CJK_TINY,
            lw=0.55,
            radius=0.007,
            zorder=4,
        )
    ax.text(
        x + w / 2,
        y + 0.017,
        "残差 + FFN + LayerNorm",
        ha="center",
        va="center",
        fontproperties=CJK_TINY,
        color=DARK_GREY,
        zorder=6,
    )


def draw_distribution(ax, x, y, w, h, title, color):
    rounded_box(ax, x, y, w, h, edge=color, fill=WHITE, lw=0.75, radius=0.012, zorder=3)
    ax.text(
        x + 0.012,
        y + h - 0.022,
        title,
        ha="left",
        va="center",
        fontproperties=CJK_HEAD,
        color=color,
        zorder=6,
    )
    xs = np.linspace(-2.5, 2.5, 60)
    ys = np.exp(-0.5 * xs**2)
    curve_x = x + 0.050 + (xs + 2.5) / 5.0 * (w - 0.068)
    curve_y = y + 0.018 + ys / ys.max() * (h - 0.047)
    ax.plot(curve_x, curve_y, color=color, lw=0.85, zorder=5)
    ax.fill_between(curve_x, y + 0.018, curve_y, color=color, alpha=0.12, zorder=4)
    mu_x = x + 0.050 + (w - 0.068) / 2
    ax.plot([mu_x, mu_x], [y + 0.018, y + h - 0.030], color=color, lw=0.6, ls=(0, (2, 2)), zorder=5)
    ax.text(
        mu_x,
        y + 0.016,
        "$\\hat{\\mu}$",
        ha="center",
        va="bottom",
        fontsize=5.3,
        color=DARK_GREY,
        zorder=6,
    )


def draw_panel_a(ax):
    x, y, w, h = 0.020, 0.055, 0.258, 0.880
    panel_background(ax, x, y, w, h, "a", "动态监督的构建", PALE_GREY)

    atlas_y = 0.685
    rounded_box(ax, 0.040, atlas_y, 0.218, 0.162, edge=BLUE, fill=WHITE, lw=0.75, radius=0.014)
    ax.text(0.052, atlas_y + 0.137, "ATLAS", ha="left", va="center", fontsize=6.2, color=BLUE, fontweight="bold")
    ax.text(0.103, atlas_y + 0.137, "3 条重复 RMSF", ha="left", va="center", fontproperties=CJK_SMALL, color=DARK_GREY)
    draw_traces(ax, 0.066, atlas_y + 0.034, 0.098, 0.082)
    chip(ax, 0.181, atlas_y + 0.082, 0.055, "$\\overline{y}$", BLUE, PALE_BLUE)
    chip(ax, 0.181, atlas_y + 0.038, 0.055, "$s^2$", BLUE, PALE_BLUE)
    arrow(ax, (0.164, atlas_y + 0.080), (0.181, atlas_y + 0.080), color=BLUE, lw=0.7)

    dyn_y = 0.410
    rounded_box(ax, 0.040, dyn_y, 0.218, 0.230, edge=ORANGE, fill=WHITE, lw=0.75, radius=0.014)
    ax.text(0.052, dyn_y + 0.204, "DynaRepo", ha="left", va="center", fontsize=6.2, color=ORANGE, fontweight="bold")
    ax.text(0.119, dyn_y + 0.204, "非重叠时间窗口", ha="left", va="center", fontproperties=CJK_SMALL, color=DARK_GREY)
    draw_windows(ax, 0.057, dyn_y + 0.125, 0.181, 0.043)
    chip(ax, 0.055, dyn_y + 0.062, 0.050, "RMSF", BLUE, PALE_BLUE)
    chip(ax, 0.112, dyn_y + 0.062, 0.050, "$H_\\chi$", PURPLE, PALE_PURPLE)
    chip(ax, 0.169, dyn_y + 0.062, 0.055, "SASA*", GREEN, PALE_GREEN)
    ax.text(
        0.149,
        dyn_y + 0.025,
        "逐链对齐 · 联合旋转状态 · 表面采样",
        ha="center",
        va="center",
        fontproperties=CJK_TINY,
        color=DARK_GREY,
    )

    label_x, label_y, label_w, label_h = 0.047, 0.136, 0.204, 0.214
    rounded_box(ax, label_x, label_y, label_w, label_h, edge=DARK_GREY, fill=WHITE, lw=0.80, radius=0.014)
    draw_label_tensor(ax, label_x, label_y, label_w, label_h)
    arrow(ax, (0.149, atlas_y), (0.149, 0.640), color=BLUE, lw=0.75)
    arrow(ax, (0.149, dyn_y), (0.149, label_y + label_h), color=ORANGE, lw=0.75)
    ax.text(
        0.149,
        0.096,
        "缺失指标由独立掩码排除",
        ha="center",
        va="center",
        fontproperties=CJK_SMALL,
        color=DARK_GREY,
    )


def draw_panel_b(ax):
    x, y, w, h = 0.296, 0.055, 0.446, 0.880
    panel_background(ax, x, y, w, h, "b", "单结构状态编码器", "#F4F7FA")

    # Sequence and amino-acid type channel
    rounded_box(ax, 0.318, 0.708, 0.124, 0.150, edge=BLUE, fill=WHITE, lw=0.75, radius=0.014)
    ax.text(0.380, 0.829, "残基类别", ha="center", va="center", fontproperties=CJK_HEAD, color=BLUE)
    draw_token_strip(
        ax,
        0.334,
        0.785,
        0.092,
        [BLUE, "#4393C3", BLUE, "#74ADD1", BLUE],
        labels=["A", "R", "G", "L", "X"],
        height=0.025,
    )
    ax.text(0.380, 0.754, "主要氨基酸嵌入 128", ha="center", va="center", fontproperties=CJK_TINY, color=DARK_GREY)
    ax.text(0.380, 0.730, "辅助类型嵌入 32", ha="center", va="center", fontproperties=CJK_TINY, color=DARK_GREY)

    # Condition channel
    rounded_box(ax, 0.457, 0.708, 0.110, 0.150, edge=ORANGE, fill=WHITE, lw=0.75, radius=0.014)
    ax.text(0.512, 0.829, "实验条件", ha="center", va="center", fontproperties=CJK_HEAD, color=ORANGE)
    chip(ax, 0.472, 0.783, 0.025, "$T$", ORANGE, PALE_ORANGE)
    chip(ax, 0.500, 0.783, 0.025, "$t$", ORANGE, PALE_ORANGE)
    chip(ax, 0.528, 0.783, 0.025, "src", ORANGE, PALE_ORANGE)
    ax.text(0.512, 0.742, "条件 MLP", ha="center", va="center", fontproperties=CJK_SMALL, color=DARK_GREY)
    ax.text(0.512, 0.721, "$3\\rightarrow128$", ha="center", va="center", fontsize=5.1, color=DARK_GREY)

    # Backbone geometry channel
    rounded_box(ax, 0.582, 0.688, 0.136, 0.170, edge=GREEN, fill=WHITE, lw=0.75, radius=0.014)
    ax.text(0.650, 0.829, "主链几何", ha="center", va="center", fontproperties=CJK_HEAD, color=GREEN)
    draw_backbone_icon(ax, 0.598, 0.744, 0.104, 0.065)
    ax.text(0.650, 0.713, "N–C$\\alpha$–C 局部坐标系", ha="center", va="center", fontproperties=CJK_TINY, color=DARK_GREY)

    # Fusion ribbon
    rounded_box(
        ax,
        0.365,
        0.590,
        0.244,
        0.073,
        "$[\\mathbf{s}_i\\Vert\\mathbf{a}_i]\\rightarrow128$  +  $\\mathbf{e}^{cond}$"
        "\n初始逐残基表示 $\\mathbf{h}^{(0)}$",
        edge=BLUE,
        fill=PALE_BLUE,
        font=CJK_SMALL,
        lw=0.85,
        radius=0.014,
    )
    arrow(ax, (0.380, 0.708), (0.416, 0.663), color=BLUE, lw=0.75)
    arrow(ax, (0.512, 0.708), (0.512, 0.663), color=ORANGE, lw=0.75)

    # Pair geometry sidecar
    rounded_box(
        ax,
        0.604,
        0.548,
        0.114,
        0.080,
        "残基对几何  $29\\rightarrow48$\nRBF 距离 · 同链\n序列间隔 · 局部方向",
        edge=GREEN,
        fill=PALE_GREEN,
        font=CJK_TINY,
        lw=0.75,
        radius=0.014,
    )
    arrow(ax, (0.650, 0.688), (0.661, 0.628), color=GREEN, lw=0.75)

    # Hero IPA block
    ipa_x, ipa_y, ipa_w, ipa_h = 0.350, 0.285, 0.312, 0.260
    draw_ipa_stack(ax, ipa_x, ipa_y, ipa_w, ipa_h)
    arrow(ax, (0.487, 0.590), (0.487, ipa_y + ipa_h), color=BLUE, lw=0.95)
    pair_row_y = ipa_y + 0.083 + 0.014
    ax.plot([0.681, 0.681], [0.548, pair_row_y], color=GREEN, lw=0.70, zorder=3)
    arrow(ax, (0.681, pair_row_y), (ipa_x + ipa_w - 0.020, pair_row_y), color=GREEN, lw=0.70)

    # Encoded residue ribbon
    rounded_box(ax, 0.375, 0.169, 0.270, 0.072, edge=PURPLE, fill=WHITE, lw=0.75, radius=0.014)
    draw_token_strip(
        ax,
        0.396,
        0.192,
        0.152,
        [PURPLE, "#9970AB", PURPLE, "#B58BC2", PURPLE, "#9970AB"],
        height=0.025,
    )
    ax.text(0.587, 0.205, "$\\mathbf{h}_i\\in\\mathbb{R}^{128}$", ha="center", va="center", fontsize=5.6, color=PURPLE)
    arrow(ax, (0.512, ipa_y), (0.512, 0.241), color=PURPLE, lw=0.95)

    ax.text(
        0.519,
        0.132,
        "全局刚体变换不变",
        ha="center",
        va="center",
        fontproperties=CJK_SMALL,
        color=PURPLE,
    )
    ax.text(
        0.519,
        0.107,
        "坐标固定 · 仅更新特征",
        ha="center",
        va="center",
        fontproperties=CJK_SMALL,
        color=DARK_GREY,
    )


def draw_panel_c(ax):
    x, y, w, h = 0.760, 0.055, 0.220, 0.880
    panel_background(ax, x, y, w, h, "c", "绝对动态分布学习", PALE_RED)

    rounded_box(
        ax,
        0.787,
        0.764,
        0.166,
        0.082,
        "逐残基概率头\nLayerNorm · MLP · $3\\times2$",
        edge=RED,
        fill=WHITE,
        font=CJK_SMALL,
        lw=0.85,
        radius=0.014,
    )
    ax.plot([0.645, 0.748], [0.205, 0.205], color=PURPLE, lw=0.75, zorder=3)
    ax.plot([0.748, 0.748], [0.205, 0.805], color=PURPLE, lw=0.75, zorder=3)
    ax.text(
        0.696,
        0.219,
        "逐残基表示",
        ha="center",
        va="bottom",
        fontproperties=CJK_SMALL,
        color=PURPLE,
        zorder=6,
    )
    arrow(ax, (0.748, 0.805), (0.787, 0.805), color=PURPLE, lw=0.75)

    draw_distribution(ax, 0.787, 0.626, 0.166, 0.104, "RMSF", BLUE)
    draw_distribution(ax, 0.787, 0.494, 0.166, 0.104, "$H_\\chi$", PURPLE)
    draw_distribution(ax, 0.787, 0.362, 0.166, 0.104, "SASA", GREEN)
    arrow(ax, (0.944, 0.764), (0.944, 0.730), color=RED, lw=0.75)
    ax.text(
        0.870,
        0.744,
        "每项输出 $(\\hat{\\mu},\\hat{\\ell})$",
        ha="center",
        va="center",
        fontproperties=CJK_SMALL,
        color=DARK_GREY,
        zorder=6,
    )
    ax.plot([0.870, 0.870], [0.626, 0.598], color=RED, lw=0.65, zorder=3)
    ax.plot([0.870, 0.870], [0.494, 0.466], color=RED, lw=0.65, zorder=3)

    loss_x, loss_y, loss_w, loss_h = 0.787, 0.112, 0.166, 0.188
    rounded_box(ax, loss_x, loss_y, loss_w, loss_h, edge=RED, fill=WHITE, lw=0.95, radius=0.016)
    ax.text(
        loss_x + loss_w / 2,
        loss_y + loss_h - 0.030,
        "噪声感知异方差损失",
        ha="center",
        va="center",
        fontproperties=CJK_HEAD,
        color=RED,
    )
    ax.text(
        loss_x + loss_w / 2,
        loss_y + 0.100,
        "$v_{ik}=\\exp(\\hat{\\ell}_{ik})+s_{ik}^{2}$",
        ha="center",
        va="center",
        fontsize=5.5,
        color=BLACK,
    )
    ax.text(
        loss_x + loss_w / 2,
        loss_y + 0.065,
        "$\\ell_{ik}=\\frac{1}{2}\\left["
        "\\frac{(\\hat{\\mu}_{ik}-y_{ik})^2}{v_{ik}}+\\log v_{ik}\\right]$",
        ha="center",
        va="center",
        fontsize=5.3,
        color=BLACK,
    )
    ax.text(
        loss_x + loss_w / 2,
        loss_y + 0.026,
        "按有效掩码加权平均",
        ha="center",
        va="center",
        fontproperties=CJK_SMALL,
        color=RED,
    )
    arrow(ax, (0.870, 0.362), (0.870, loss_y + loss_h), color=RED, lw=0.80)

    # Supervision route from panel (a) to the loss, deliberately separated from the forward path.
    route_y = 0.073
    route_color = "#777777"
    ax.plot([0.251, 0.251], [0.198, route_y], color=route_color, lw=0.55, ls=(0, (2, 2)), zorder=3)
    ax.plot([0.251, 0.766], [route_y, route_y], color=route_color, lw=0.55, ls=(0, (2, 2)), zorder=3)
    arrow(ax, (0.766, route_y), (loss_x, loss_y + 0.040), color=route_color, lw=0.55, linestyle=(0, (2, 2)))
    ax.text(
        0.690,
        route_y + 0.006,
        "标签监督：$\\mathbf{Y},\\mathbf{S}^{2},\\mathbf{M}$",
        ha="center",
        va="bottom",
        fontproperties=CJK_SMALL,
        color=DARK_GREY,
        zorder=6,
    )


def main():
    fig = plt.figure(figsize=(183 * mm_to_inch, 100 * mm_to_inch), facecolor=WHITE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.020,
        0.974,
        "FlexPPI 预训练：由分子动力学监督学习单结构状态的逐残基动态分布",
        ha="left",
        va="top",
        color=BLACK,
        fontproperties=CJK_TITLE,
    )
    # Direct labeling is used throughout; the schematic intentionally has no legend.
    draw_panel_a(ax)
    draw_panel_b(ax)
    draw_panel_c(ax)

    output_base = Path(__file__).resolve().parent / "fig_flexppi_pretraining"
    save_cns_figure(fig, str(output_base))
    plt.close(fig)
    print(f"Saved {output_base}.pdf and {output_base}.png")


if __name__ == "__main__":
    main()
