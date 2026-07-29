# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) state-encoder schematic -> cross-type inherit -> param inherit
# (b) paired-response schematic -> cross-type inherit -> param inherit
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
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"
mpl.rcParams["savefig.pad_inches"] = 0.0

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
from matplotlib.font_manager import FontProperties
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


MM_TO_INCH = 1 / 25.4
FIGURE_WIDTH_MM = 183
FIGURE_HEIGHT_MM = 104
BLUE = CATEGORICAL[0]
RED = CATEGORICAL[1]
GREEN = CATEGORICAL[2]
ORANGE = CATEGORICAL[3]
PURPLE = CATEGORICAL[4]
DARK_GREY = "#555555"
LIGHT_GREY = "#F3F4F6"
WHITE = "#FFFFFF"


def find_cjk_font():
    candidates = [
        Path(r"C:\Windows\Fonts\msyh.ttc"),
        Path(r"C:\Windows\Fonts\msyhbd.ttc"),
        Path(r"C:\Windows\Fonts\simhei.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return FontProperties(fname=str(candidate), size=7)
    return FontProperties(family="sans-serif", size=7)


CJK = find_cjk_font()
CJK_BOLD = CJK.copy()
CJK_BOLD.set_weight("bold")
CJK_BOLD.set_size(8.5)
CJK_SMALL = CJK.copy()
CJK_SMALL.set_size(6.1)
CJK_TINY = CJK.copy()
CJK_TINY.set_size(5.4)


def rounded_box(ax, x, y, w, h, text, edge, fill, font=CJK_SMALL,
                lw=0.8, radius=0.015, zorder=2, text_color=BLACK):
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        linewidth=lw, edgecolor=edge, facecolor=fill, zorder=zorder,
    )
    ax.add_patch(box)
    ax.text(
        x + w / 2, y + h / 2, text,
        ha="center", va="center", color=text_color,
        fontproperties=font, linespacing=1.22, zorder=zorder + 1,
    )
    return box


def arrow(ax, start, end, color=DARK_GREY, lw=0.8, style="-|>",
          linestyle="-", connectionstyle="arc3,rad=0", zorder=1):
    patch = FancyArrowPatch(
        start, end, arrowstyle=style, mutation_scale=8,
        linewidth=lw, color=color, linestyle=linestyle,
        connectionstyle=connectionstyle, shrinkA=1.5, shrinkB=1.5,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def panel_frame(ax, x, y, w, h, label, title):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=WHITE,
                           edgecolor="#C9CDD3", linewidth=0.7, zorder=0))
    ax.text(x + 0.012, y + h - 0.025, label, fontsize=9, fontweight="bold",
            ha="left", va="top", color=BLACK)
    ax.text(x + 0.043, y + h - 0.026, title, fontproperties=CJK_BOLD,
            ha="left", va="top", color=BLACK)


def draw_state_encoder(ax):
    panel_frame(ax, 0.015, 0.055, 0.49, 0.90, "a", "单个结构状态的逐残基编码")

    rounded_box(ax, 0.038, 0.775, 0.125, 0.105,
                "残基类别编号 $x_i$\n21 类；不是 one-hot", BLUE, "#EAF2FA")
    rounded_box(ax, 0.038, 0.625, 0.125, 0.105,
                "实验/模拟条件 $\\mathbf{c}$\n温度 · 窗口 · 数据源", ORANGE, "#FFF3DF")
    rounded_box(ax, 0.038, 0.420, 0.125, 0.125,
                "主链坐标\nN–$C_\\alpha$–C–O\n$N\\times4\\times3$", GREEN, "#EAF5ED")

    rounded_box(ax, 0.195, 0.790, 0.125, 0.078,
                "序列通道 $\\mathbf{s}_i$\nlearned 或分链 ESM2", BLUE, "#EAF2FA")
    rounded_box(ax, 0.195, 0.675, 0.125, 0.078,
                "氨基酸嵌入 $\\mathbf{e}^{aa}_i$\n独立 32 维查表", BLUE, "#EAF2FA")
    rounded_box(ax, 0.195, 0.565, 0.125, 0.068,
                "条件 MLP\n$3\\rightarrow128$", ORANGE, "#FFF3DF")
    rounded_box(ax, 0.195, 0.405, 0.125, 0.145,
                "固定局部坐标架\n+ 成对几何 $\\mathbf{p}_{ij}$\n距离 · 同链 · 序列间隔 · 方向\n$N\\times N\\times48$", GREEN, "#EAF5ED", font=CJK_TINY)

    rounded_box(ax, 0.350, 0.675, 0.125, 0.150,
                "初始逐残基表示\n拼接 + 投影 + 条件\n默认 learned\n$(128+32)\\rightarrow128$\n$\\mathbf{h}^{(0)}: N\\times128$", BLUE, "#EAF2FA", font=CJK_TINY)
    rounded_box(ax, 0.350, 0.440, 0.125, 0.155,
                "固定坐标架 IPA × 4\n8 个注意力头\n链内 + 链间传播\n只更新特征，不移动坐标", PURPLE, "#F3ECF6", font=CJK_TINY)
    rounded_box(ax, 0.350, 0.205, 0.125, 0.145,
                "绝对状态头\n$128\\rightarrow128\\rightarrow6$\n$\\hat{\\mathbf{y}},\\log\\hat{\\boldsymbol{\\sigma}}^2$\n均为 $N\\times3$", RED, "#F9EAEA", font=CJK_TINY)

    arrow(ax, (0.163, 0.828), (0.195, 0.828), BLUE)
    arrow(ax, (0.163, 0.808), (0.195, 0.714), BLUE, connectionstyle="arc3,rad=0.15")
    arrow(ax, (0.163, 0.678), (0.195, 0.599), ORANGE)
    arrow(ax, (0.163, 0.482), (0.195, 0.482), GREEN)
    arrow(ax, (0.320, 0.829), (0.350, 0.775), BLUE)
    arrow(ax, (0.320, 0.714), (0.350, 0.738), BLUE)
    arrow(ax, (0.320, 0.599), (0.350, 0.698), ORANGE,
          connectionstyle="arc3,rad=-0.18")
    arrow(ax, (0.413, 0.675), (0.413, 0.595), DARK_GREY)
    arrow(ax, (0.320, 0.478), (0.350, 0.500), GREEN)
    arrow(ax, (0.413, 0.440), (0.413, 0.350), DARK_GREY)

    ax.text(0.257, 0.355,
            "全部链参与注意力；目标链掩码仅控制监督与输出读取",
            fontproperties=CJK_TINY, ha="center", va="center", color=DARK_GREY)


def draw_pair_response(ax):
    panel_frame(ax, 0.520, 0.055, 0.465, 0.90, "b", "共享编码器上的 apo/holo 配对响应")

    rounded_box(ax, 0.548, 0.785, 0.105, 0.085,
                "apo 状态\n目标链 ± 伙伴环境", BLUE, "#EAF2FA")
    rounded_box(ax, 0.548, 0.635, 0.105, 0.085,
                "holo 状态\n目标链 + 伙伴环境", GREEN, "#EAF5ED")

    rounded_box(ax, 0.690, 0.770, 0.115, 0.105,
                "状态编码器\n共享同一组参数", PURPLE, "#F3ECF6")
    rounded_box(ax, 0.690, 0.620, 0.115, 0.105,
                "状态编码器\n共享同一组参数", PURPLE, "#F3ECF6")
    ax.plot([0.675, 0.675], [0.628, 0.866], color=PURPLE, lw=0.8, ls="--")
    ax.text(0.666, 0.747, "共享权重", rotation=90, fontproperties=CJK_TINY,
            ha="center", va="center", color=PURPLE)

    rounded_box(ax, 0.838, 0.770, 0.115, 0.105,
                "$\\mathbf{h}^{apo},\\hat{\\mathbf{y}}^{apo}$\n$N_a\\times128,\\ N_a\\times3$", BLUE, "#EAF2FA")
    rounded_box(ax, 0.838, 0.620, 0.115, 0.105,
                "$\\mathbf{h}^{holo},\\hat{\\mathbf{y}}^{holo}$\n$N_h\\times128,\\ N_h\\times3$", GREEN, "#EAF5ED")

    arrow(ax, (0.653, 0.828), (0.690, 0.828), BLUE)
    arrow(ax, (0.653, 0.678), (0.690, 0.678), GREEN)
    arrow(ax, (0.805, 0.828), (0.838, 0.828), PURPLE)
    arrow(ax, (0.805, 0.678), (0.838, 0.678), PURPLE)

    rounded_box(ax, 0.565, 0.440, 0.155, 0.115,
                "可靠残基映射 $\\mathcal{M}$\n按对应行取出 $M$ 个残基\n无映射/缺坐标位置被掩码", DARK_GREY, LIGHT_GREY, font=CJK_TINY)
    rounded_box(ax, 0.748, 0.415, 0.205, 0.165,
                "响应特征拼接 $\\mathbf{z}_i$\n$\\mathbf{h}^{apo},\\mathbf{h}^{holo},\\Delta\\mathbf{h},$\n$\\mathbf{h}^{holo}\\odot\\mathbf{h}^{apo},\\Delta\\hat{\\mathbf{y}}^{sub}$\n$4\\times128+3=515$", RED, "#F9EAEA", font=CJK_TINY)

    arrow(ax, (0.868, 0.770), (0.690, 0.555), BLUE,
          connectionstyle="arc3,rad=0.14")
    arrow(ax, (0.900, 0.620), (0.700, 0.555), GREEN,
          connectionstyle="arc3,rad=-0.10")
    arrow(ax, (0.720, 0.498), (0.748, 0.498), DARK_GREY)

    rounded_box(ax, 0.565, 0.210, 0.155, 0.110,
                "状态相减基线\n$\\Delta\\hat{\\mathbf{y}}^{sub}=\\hat{\\mathbf{y}}^{holo}-\\hat{\\mathbf{y}}^{apo}$\n$M\\times3$", DARK_GREY, LIGHT_GREY, font=CJK_TINY)
    rounded_box(ax, 0.748, 0.225, 0.100, 0.095,
                "响应 MLP\n$515\\rightarrow256\\rightarrow6$", RED, "#F9EAEA")
    rounded_box(ax, 0.875, 0.205, 0.080, 0.135,
                "直接响应\n$\\Delta\\hat{\\mathbf{y}}^{dir}$\n+ 对数方差\n$M\\times3$", RED, "#F9EAEA", font=CJK_TINY)

    arrow(ax, (0.850, 0.415), (0.800, 0.320), RED,
          connectionstyle="arc3,rad=0.08")
    arrow(ax, (0.848, 0.272), (0.875, 0.272), RED)
    arrow(ax, (0.642, 0.440), (0.642, 0.320), DARK_GREY)
    arrow(ax, (0.720, 0.265), (0.875, 0.245), DARK_GREY, lw=0.7,
          linestyle="--", connectionstyle="arc3,rad=-0.12")
    ax.text(0.790, 0.208, "SmoothL1 一致性", fontproperties=CJK_TINY,
            ha="center", va="center", color=DARK_GREY)

    ax.text(0.752, 0.125,
            "训练：apo/holo 状态 NLL + 直接响应 NLL + 一致性约束\n"
            "推理：并列报告直接响应与状态相减，不事后择优",
            fontproperties=CJK_SMALL, ha="center", va="center", color=BLACK,
            linespacing=1.35)


def main():
    fig = plt.figure(figsize=(FIGURE_WIDTH_MM * MM_TO_INCH,
                              FIGURE_HEIGHT_MM * MM_TO_INCH), facecolor=WHITE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    draw_state_encoder(ax)
    draw_pair_response(ax)

    output_base = Path(__file__).resolve().parent / "fig_flexppi_architecture"
    save_cns_figure(fig, str(output_base))
    plt.close(fig)
    print(f"Saved {output_base}.pdf and {output_base}.png")


if __name__ == "__main__":
    main()