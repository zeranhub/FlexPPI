# Asset matching confirmation: a/b inherit MarginalDensity/plot_VariableCorrelation.py;
# c/d inherit LineTrend/plot_trend.py; e inherits heatmap/plot_composition.py;
# f inherits BarComposition/plot_brute_force.py. Source scripts and companion
# heatmap/composition images were inspected. These are cross-type inheritance:
# no regression/KDE, time-series or model-prompt data exist in this label audit.
# All records, empirical CDF steps, eligible amino acids and zero-SD classes
# are used. Palette/fonts/spines/annotated cells/hatch evidence are retained.
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
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch
from matplotlib.ticker import FixedLocator, FixedFormatter
import fitz

HERE = Path(__file__).resolve().parent
SOURCES = ["atlas", "dynarepo", "itrace", "thermopcd"]
DISPLAY = dict(atlas="ATLAS", dynarepo="DynaRepo", itrace="ITrace", thermopcd="ThermoPCD")
COLORS = dict(zip(SOURCES, CATEGORICAL[:4]))
MARKERS = dict(zip(SOURCES, ["o", "^", "s", "D"]))
AA = list("AGCSTVNDHILFPWYQEMRK")
Q = [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def title(ax, letter, text):
    ax.set_title(text, loc="left", pad=9, fontsize=8)
    box = ax.get_position()
    ax.figure.text(.016 if box.x0 < .5 else .51, box.y1 + .010,
                   letter, fontweight="bold", fontsize=10)


def heatmap(fig, ax, values):
    cmap = LinearSegmentedColormap.from_list("label_sequential", SEQUENTIAL)
    cmap.set_bad("#E6E6E6")
    image = ax.imshow(values, aspect="auto", vmin=0, vmax=1, cmap=cmap, interpolation="none")
    ax.set_xticks(range(20), [f"{a}\n{q}" for a, q in zip(AA, Q)])
    ax.set_yticks(range(4), [DISPLAY[s] for s in SOURCES])
    for i in range(4):
        for j in range(20):
            v = values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}".lstrip("0"), ha="center", va="center", fontsize=7.5,
                        color="white" if v > .68 else BLACK)
            else:
                ax.add_patch(Rectangle((j-.5, i-.5), 1, 1, facecolor="none", edgecolor="white", hatch="///", linewidth=0))
    ax.set_xlabel("Amino acid / canonical χ count Q", labelpad=3)
    colorbar = fig.colorbar(image, cax=fig.add_axes([.943, .335, .012, .100]))
    colorbar.set_ticks([0, .5, 1]); colorbar.ax.tick_params(labelsize=7.5)
    return image


def main():
    audit = json.loads((HERE / "label_analysis_audit.json").read_bytes())
    assert audit["actual_labels_analysis"] and not audit["final_neural_results_complete"] and audit["records"] == 3095
    for name, digest in audit["output_sha256"].items():
        assert sha(HERE / name) == digest
    per = pd.read_csv(HERE / "per_record_labels.csv", keep_default_na=False, na_values=[""])
    aa = pd.read_csv(HERE / "by_source_amino_acid_labels.csv", keep_default_na=False, na_values=[""])
    source = pd.read_csv(HERE / "by_source_labels.csv").set_index("dataset")
    assert len(per) == 3095 and per.sample_id.is_unique
    mpl.rcParams.update({"xtick.labelsize": 7.5, "ytick.labelsize": 7.5})
    fig = plt.figure(figsize=(183/25.4, 200/25.4))
    axes = [fig.add_axes(x) for x in [
        [.115, .760, .350, .153], [.605, .760, .350, .153],
        [.115, .525, .350, .145], [.605, .525, .350, .145],
        [.140, .335, .785, .100], [.185, .130, .710, .120]]]
    fig.suptitle("Observed label variation across windows and replicas", y=.99)
    handles = [Line2D([], [], color=COLORS[s], marker=MARKERS[s], linestyle="none", markersize=4,
                     label=f"{DISPLAY[s]} (n={int((per.dataset == s).sum())})") for s in SOURCES]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5,.978), ncol=4, fontsize=7.5,
               handletextpad=.3, columnspacing=.8)
    point_counts, curves = {}, []
    for ax, metric, letter in [(axes[0], "rmsf", "a"), (axes[1], "hchi", "b")]:
        for s in SOURCES:
            d = per[per.dataset == s]; x = d[metric + "_label_mean"]; y = d[metric + "_observed_sd_mean"]
            finite = np.isfinite(x) & np.isfinite(y)
            ax.scatter(x[finite], y[finite], s=5, marker=MARKERS[s], color=COLORS[s], alpha=.52, linewidths=0, rasterized=True)
            point_counts[metric + "_" + s] = dict(plotted=int(finite.sum()), unavailable=int((~finite).sum()))
        if metric == "rmsf":
            ax.set_xscale("symlog", linthresh=.1); ax.set_yscale("symlog", linthresh=.01)
            ax.xaxis.set_major_locator(FixedLocator([.3, 1, 3, 10]))
            ax.xaxis.set_major_formatter(FixedFormatter(["0.3", "1", "3", "10"]))
            ax.yaxis.set_major_locator(FixedLocator([0, .01, .1, 1, 3]))
            ax.yaxis.set_major_formatter(FixedFormatter(["0", "0.01", "0.1", "1", "3"]))
            unit = "Å"; text = "RMSF: label amplitude and variation"
        else:
            ax.set_xlim(left=0); ax.set_ylim(bottom=0)
            unit = "nats"; text = "Hχ: label amplitude and variation"
        ax.set_xlabel(f"Record mean label ({unit})"); ax.set_ylabel(f"Record mean observed SD ({unit})")
        ax.margins(.045)
        title(ax, letter, text)
    for ax, metric, letter in [(axes[2], "rmsf", "c"), (axes[3], "hchi", "d")]:
        for s in SOURCES:
            d = per[per.dataset == s]
            for statistic, style in [("mean", "-"), ("median", "--")]:
                values = d[metric + "_observed_sd_" + statistic].to_numpy(float)
                values = np.sort(values[np.isfinite(values)])
                ecdf = np.arange(1, len(values)+1)/len(values)
                ax.step(values, ecdf, where="post", color=COLORS[s], linestyle=style, linewidth=1)
                curves += [dict(dataset=s, metric=metric, record_sd_statistic=statistic, x=float(x), ecdf=float(y)) for x, y in zip(values,ecdf)]
        ax.set_xscale("symlog", linthresh=.005); ax.set_ylim(0,1.02)
        ax.set_xlim(left=-.00015)
        ax.xaxis.set_major_locator(FixedLocator([0, .005, .05, .5, 5]))
        ax.xaxis.set_major_formatter(FixedFormatter(["0", "0.005", "0.05", "0.5", "5"]))
        ax.set_ylabel("Fraction of records")
        ax.set_xlabel("Record observed SD (" + ("Å" if metric == "rmsf" else "nats") + ")")
        title(ax, letter, ("RMSF" if metric == "rmsf" else "Hχ") + " SD: mean (solid), median (dashed)")
    coverage, normalized = np.full((4,20),np.nan), np.full((4,20),np.nan)
    for i,s in enumerate(SOURCES):
        d=aa[aa.dataset==s].set_index("amino_acid")
        for j,a in enumerate(AA):
            if a in d.index:
                if Q[j]>0:coverage[i,j]=d.loc[a,"hchi_q_eligible_coverage"]
                normalized[i,j]=d.loc[a,"hchi_normalized_label_mean"]
    heatmap(fig, axes[4], normalized)
    title(axes[4], "e", "Mean Hχ / (Q ln 3): all eligible residues, within source and amino acid")
    fractions=[]
    for s in SOURCES:
        d=aa[aa.dataset==s].set_index("amino_acid")
        values=[]
        for a in AA:
            v=d.loc[a,"hchi_valid"] if a in d.index else 0
            zero=d.loc[a,"hchi_valid_zero_label"] if a in d.index else 0
            values.append(zero/v if v else np.nan)
            fractions.append(dict(dataset=s,amino_acid=a,hchi_valid=int(v),hchi_valid_zero_label=int(zero),fraction=values[-1]))
    classes=["positive_sd","zero_n1_proven","zero_n_gt1_saved","zero_n_unknown"]
    fills=["#555555","#D9D9D9","#999999","#FFFFFF"]
    hatches=[None,None,"...","////"]
    bars=[]
    for s in SOURCES:
        d=per[per.dataset==s]
        for metric in ["rmsf","hchi"]:
            valid=int(source.loc[s,metric+"_valid"]); zero=int(source.loc[s,metric+"_valid_zero_variance"])
            single=int(d[metric+"_valid_single_observation_proven"].sum()); unknown=int(d[metric+"_valid_zero_variance_count_unknown"].sum())
            counts=[valid-zero,single,zero-single-unknown,unknown];assert min(counts)>=0 and sum(counts)==valid
            bars.append(dict(dataset=s,metric=metric,valid_residues=valid,**dict(zip(classes,counts))))
    for index,row in enumerate(bars):
        left=0
        for j,c in enumerate(classes):
            width=row[c]/row["valid_residues"]
            axes[5].barh(index,width,left=left,color=fills[j],hatch=hatches[j],edgecolor="#555555",linewidth=.35,height=.72)
            left+=width
        axes[5].text(1.013, index, f'{row["valid_residues"]:,}', va="center", fontsize=7.5, clip_on=False)
    axes[5].set_yticks(range(8),[DISPLAY[x["dataset"]]+" · "+("R" if x["metric"]=="rmsf" else "Hχ") for x in bars])
    axes[5].invert_yaxis();axes[5].set_xlim(0,1);axes[5].set_xlabel("Fraction of valid residue labels")
    title(axes[5],"f","Observed SD and evidence for the per-residue observation count n")
    legend=[Patch(facecolor=fills[i],edgecolor="#555555",hatch=hatches[i],label=x) for i,x in enumerate(["SD > 0","SD = 0; n = 1 proven","SD = 0; n > 1 saved","SD = 0; n not saved"])]
    fig.legend(handles=legend,loc="lower center",bbox_to_anchor=(.50,.025),fontsize=7.5,ncol=2,columnspacing=1.5,handlelength=1.8)
    fig.text(.12,.010,"All records and valid zeros retained. A/G: χ undefined (hatched). SD is not SEM.",fontsize=7.5,ha="left")
    stem=HERE/"figS_corrected_observed_label_variation"
    # Exact final physical size; override tight cropping without shrinking text.
    with mpl.rc_context({"savefig.bbox": None}):
        for extension in ["pdf", "png", "svg"]:
            fig.savefig(str(stem)+"."+extension, bbox_inches=None, dpi=300)
    plt.close(fig)
    pd.DataFrame(curves).to_csv(HERE/"plot_ecdf_all_record_points.csv",index=False)
    pd.DataFrame(fractions).to_csv(HERE/"plot_zero_entropy_amino_acids.csv",index=False,na_rep="")
    pd.DataFrame(bars).to_csv(HERE/"plot_zero_variance_count_evidence.csv",index=False)
    pd.DataFrame(coverage,index=SOURCES,columns=AA).to_csv(HERE/"plot_hchi_coverage_matrix.csv",na_rep="")
    pd.DataFrame(normalized,index=SOURCES,columns=AA).to_csv(HERE/"plot_hchi_normalized_mean_matrix.csv",na_rep="")
    outputs={p.name:sha(p) for p in HERE.glob("figS*.*")}
    doc = fitz.open(str(stem)+".pdf")
    minimum_font = min(span["size"] for block in doc[0].get_text("dict")["blocks"]
                       for line in block.get("lines", []) for span in line["spans"] if span["text"].strip())
    assert minimum_font >= 7.49
    report=dict(status="rendered; visual inspection pending",records=3095,point_counts=point_counts,
                label_analysis_audit_sha256=sha(HERE/"label_analysis_audit.json"),plot_script_sha256=sha(__file__),
                output_sha256=outputs,full_ecdf_rows=len(curves),figsize_mm=[183,200],png_dpi=300,
                vector_minimum_font_points=minimum_font,panels=6,
                final_neural_results_complete=False,subsampled_records=0,variance_exclusions=0,
                heatmap_na="Ala/Gly: χ not defined, not observed zero entropy.")
    (HERE/"figure_render_audit.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
