# Academic Figure Skill Asset Confirmation (previously inspected production assets)
# (A,B) BarAblation/plot_bar_ablation.py -> journal visual adaptation -> param inherit
# (C,D) Forest from BarAblation -> journal visual adaptation -> param inherit
# (supplement) paired ECDF, heatmap, forest -> param inherit
# User requested Bioinformatics / Briefings in Bioinformatics, overriding CNS dimensions.
# All numerical summaries are reused from the already verified deliverables.
"""Compact Bioinformatics / BiB ablation artwork, with captions outside the image."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import fitz
from PIL import Image
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "output/information_ablation_figures/verified"
METHODS = ["flexppi", "nn_geometry", "nn_intrachain", "nn_sequence"]
LABELS = {"flexppi":"Full FlexPPI", "nn_geometry":"Identity masked", "nn_intrachain":"Intrachain only",
          "nn_sequence":"Sequence only", "aa":"AA-type prior", "sequence":"Sequence probe",
          "geometry":"Geometry probe", "sequence_intra":"Sequence + local geometry",
          "sequence_complex":"Sequence + geometry + contacts"}
TICK_LABELS = ["Full\nFlexPPI", "Identity\nmasked", "Intrachain\nonly", "Sequence\nonly"]
# Consistent semantic colors; marker shapes remain useful in greyscale.
COLORS = {"flexppi":"#3C70A4", "nn_geometry":"#D9836C", "nn_intrachain":"#4E948B", "nn_sequence":"#9783AE"}
MARKERS = {"flexppi":"o", "nn_geometry":"s", "nn_intrachain":"D", "nn_sequence":"^"}
METRICS = ["rmsf", "chi_entropy"]
METRIC_TITLES = {"rmsf":"RMSF", "chi_entropy":r"$H_\chi$"}
UNITS = {"rmsf":"Å", "chi_entropy":"nats"}
WIDTH_MM = 178

mpl.rcParams.update({
    "font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","Liberation Sans"],
    "font.size":8, "axes.labelsize":8, "axes.titlesize":8.5,
    "xtick.labelsize":7, "ytick.labelsize":7, "legend.fontsize":7,
    "axes.spines.top":False, "axes.spines.right":False,
    "axes.linewidth":.65, "xtick.direction":"out", "ytick.direction":"out",
    "xtick.major.width":.6, "ytick.major.width":.6,
    "xtick.major.size":2.5, "ytick.major.size":2.5,
    "legend.frameon":False, "pdf.fonttype":42, "svg.fonttype":"none",
    "savefig.facecolor":"white", "figure.facecolor":"white", "axes.facecolor":"white",
    "mathtext.fontset":"dejavusans",
})

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def panel(ax, letter, title):
    ax.set_title(title, loc="center", pad=7, weight="normal")
    ax.text(-.15,1.05,letter,transform=ax.transAxes,ha="left",va="bottom",fontsize=10,weight="bold",color="#222222")

def export(fig, name, out):
    for ext in ("pdf","svg","png"):
        fig.savefig(out/f"{name}.{ext}",dpi=300,bbox_inches=None)
    # Both journals have demanding line-art raster requirements; the vector master
    # stays authoritative and this TIFF is rendered directly from vector artists.
    fig.savefig(out/f"{name}.tif",dpi=1200,bbox_inches=None,
                pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig)

def make_main(performance, contrasts, out):
    fig=plt.figure(figsize=(WIDTH_MM/25.4,112/25.4))
    positions=[(.125,.605,.345,.315),(.615,.605,.345,.315),
               (.125,.13,.345,.305),(.615,.13,.345,.305)]
    for col,metric in enumerate(METRICS):
        ax=fig.add_axes(positions[col])
        values=performance[performance.metric.eq(metric)].set_index("method")
        panel(ax,"AB"[col],METRIC_TITLES[metric])
        for x,method in enumerate(METHODS):
            row=values.loc[method]
            ax.bar(x,row.mean_mae,width=.57,color=COLORS[method],edgecolor="none",zorder=2)
            ax.errorbar(x,row.mean_mae,yerr=[[row.mean_mae-row.ci_low],[row.ci_high-row.mean_mae]],
                        fmt="none",ecolor="#333333",elinewidth=.7,capsize=2,capthick=.7,zorder=3)
            ax.text(x,row.ci_high+.022,f"{row.mean_mae:.3f}",ha="center",va="bottom",fontsize=7)
        ax.set(xticks=range(4),xticklabels=TICK_LABELS,xlim=(-.6,3.6),
               ylim=(0,.78 if col==0 else .74),ylabel=f"MAE ({UNITS[metric]})")
        ax.set_yticks([0,.2,.4,.6])
        ax.tick_params(axis="x",length=0,pad=5)
        ax.set_axisbelow(True)

        ax=fig.add_axes(positions[col+2])
        panel(ax,"CD"[col],f"{METRIC_TITLES[metric]}: paired difference")
        values=contrasts[contrasts.metric.eq(metric)].set_index("method")
        ax.axvline(0,color="#999999",ls="--",lw=.65,zorder=0)
        for y,method in enumerate(METHODS[1:]):
            row=values.loc[method]
            ax.errorbar(row.mean_delta,y,xerr=[[row.mean_delta-row.ci_low],[row.ci_high-row.mean_delta]],
                        fmt=MARKERS[method],ms=4.6,color=COLORS[method],elinewidth=1.1,capsize=2,capthick=.7)
            offset=.012 if col==0 else .025
            ax.text(row.ci_high+offset,y,f"{row.mean_delta:+.4f}",va="center",fontsize=7)
        ax.set(yticks=range(3),yticklabels=[LABELS[m] for m in METHODS[1:]],
               ylim=(2.55,-.55),xlim=(-.025,.225) if col==0 else (-.035,.64),
               xlabel=f"ΔMAE ({UNITS[metric]})")
        ax.set_xticks([0,.05,.10,.15,.20] if col==0 else [0,.2,.4,.6])
        ax.tick_params(axis="y",length=0,pad=5)
        ax.spines["left"].set_visible(False)
    export(fig,"fig_ablation_main",out)

def make_distribution(table,out):
    fig=plt.figure(figsize=(WIDTH_MM/25.4,67/25.4))
    for col,metric in enumerate(METRICS):
        ax=fig.add_axes([.10+col*.49,.255,.36,.625])
        panel(ax,"AB"[col],METRIC_TITLES[metric])
        block=table[table.metric.eq(metric)]
        ax.axvline(0,color="#999999",ls="--",lw=.65)
        for method in METHODS[1:]:
            delta=np.sort((block[method]-block.flexppi).to_numpy())
            ax.step(delta,np.arange(1,len(delta)+1)/len(delta),where="post",lw=1.2,color=COLORS[method])
        all_values=(block[METHODS[1:]].to_numpy()-block.flexppi.to_numpy()[:,None]).ravel()
        ax.set_xscale("symlog",linthresh=.01)
        ax.set(xlim=(all_values.min()*1.2,all_values.max()*1.2),ylim=(0,1),
               xlabel=f"ΔMAE ({UNITS[metric]})",ylabel="Cumulative fraction")
        ax.set_yticks([0,.25,.5,.75,1])
        ax.set_xticks([-1,-.1,0,.1,1] if col==0 else [-.1,0,.1,1])
        ax.set_xticklabels(["−1","−0.1","0","0.1","1"] if col==0 else ["−0.1","0","0.1","1"])
    handles=[Line2D([0],[0],color=COLORS[m],lw=1.3,label=LABELS[m]) for m in METHODS[1:]]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.54,.035),ncol=3,fontsize=7)
    export(fig,"fig_ablation_distribution",out)

def make_controls(performance,out):
    fig=plt.figure(figsize=(WIDTH_MM/25.4,89/25.4))
    order=METHODS+["aa","sequence","geometry","sequence_intra","sequence_complex"]
    ys=[0,1,2,3,5,6,7,8,9]
    for col,metric in enumerate(METRICS):
        ax=fig.add_axes([.30+col*.365,.145,.285,.745])
        panel(ax,"AB"[col],METRIC_TITLES[metric])
        data=performance[performance.metric.eq(metric)].set_index("method")
        for y,method in zip(ys,order):
            row=data.loc[method]
            ax.errorbar(row.mean_mae,y,xerr=[[row.mean_mae-row.ci_low],[row.ci_high-row.mean_mae]],
                        fmt=MARKERS.get(method,"o"),ms=4,color=COLORS.get(method,"#A1A1A1"),elinewidth=1,capsize=2)
            ax.annotate(f"{row.mean_mae:.3f}",xy=(row.ci_high,y),xytext=(5,0),
                        textcoords="offset points",fontsize=6.5,va="center")
        ax.axhline(4,color="#D9D9D9",lw=.6)
        ax.set(yticks=ys,yticklabels=[LABELS[m] for m in order] if col==0 else [""]*9,
               ylim=(9.7,-.7),xlim=(.43,.96) if col==0 else (.13,.82),xlabel=f"MAE ({UNITS[metric]})")
        ax.tick_params(axis="y",length=0,pad=4)
        ax.spines["left"].set_visible(False)
    export(fig,"fig_ablation_controls",out)

def make_heterogeneity(groups,folds,out):
    fig=plt.figure(figsize=(WIDTH_MM/25.4,120/25.4))
    sources=["atlas","dynarepo","itrace","atlas|dynarepo","dynarepo|itrace"]
    source_labels=["ATLAS","DynaRepo","iTrace","ATLAS + DynaRepo","DynaRepo + iTrace"]
    cmap=LinearSegmentedColormap.from_list("paired_delta",["#3C70A4","#F7F7F7","#B95D5D"])
    for col,metric in enumerate(METRICS):
        left=.195+col*.47
        ax=fig.add_axes([left,.665,.30,.235])
        panel(ax,"AB"[col],f"{METRIC_TITLES[metric]}: data source")
        block=groups[groups.metric.eq(metric)]
        values=block.pivot(index="group",columns="model",values="mean_delta").loc[sources,METHODS[1:]]
        limit=float(np.abs(values.to_numpy()).max())
        im=ax.imshow(values,cmap=cmap,norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit),aspect="auto")
        for y in range(5):
            for x in range(3):
                v=float(values.iloc[y,x])
                rgb=np.array(cmap(im.norm(v))[:3])
                linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
                luminance=float(linear @ np.array([.2126,.7152,.0722]))
                # Choose the text with greater measured contrast on each actual cell.
                ink="white" if 1.05/(luminance+.05)>(luminance+.05)/.05 else "black"
                ax.text(x,y,f"{v:+.4f}",ha="center",va="center",fontsize=6.5,color=ink)
        ax.set_xticks(range(3),["Identity\nmasked","Intrachain\nonly","Sequence\nonly"],fontsize=7)
        ax.set_yticks(range(5),source_labels,fontsize=7)
        ax.set_xticks(np.arange(-.5,3),minor=True)
        ax.set_yticks(np.arange(-.5,5),minor=True)
        ax.grid(which="minor",color="white",lw=1)
        ax.tick_params(which="both",length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        cax=fig.add_axes([left,.558,.30,.012])
        cb=fig.colorbar(im,cax=cax,orientation="horizontal",ticks=[-limit,0,limit])
        cb.set_label(f"ΔMAE ({UNITS[metric]})",fontsize=7,labelpad=1.5)
        cb.ax.tick_params(labelsize=6.5,length=2)
        cb.outline.set_visible(False)

        ax=fig.add_axes([left,.16,.30,.265])
        panel(ax,"CD"[col],f"{METRIC_TITLES[metric]}: held-out fold")
        ax.axvline(0,color="#999999",ls="--",lw=.65)
        markers=["o","s","^","D","v"]
        for y,method in enumerate(METHODS[1:]):
            b=folds[folds.metric.eq(metric)&folds.method.eq(method)].set_index("fold")
            for fold in range(5):
                ax.scatter(b.loc[fold,"delta"],y+(fold-2)*.08,marker=markers[fold],s=17,
                           color=COLORS[method],edgecolors="white",linewidths=.4)
        v=folds[folds.metric.eq(metric)].delta
        pad=(v.max()-v.min())*.12
        ax.set(yticks=range(3),yticklabels=[LABELS[m] for m in METHODS[1:]],
               ylim=(2.65,-.65),xlim=(min(v.min()-pad,-.01),v.max()+pad),xlabel=f"ΔMAE ({UNITS[metric]})")
        ax.tick_params(axis="y",length=0,pad=4)
        ax.spines["left"].set_visible(False)
    handles=[Line2D([0],[0],marker=m,ls="none",color="#666666",ms=4,label=f"Fold {i}") for i,m in enumerate(markers)]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.54,.04),ncol=5,fontsize=7,handletextpad=.3,columnspacing=1.1)
    export(fig,"fig_ablation_heterogeneity",out)

def make_dependence(data,out):
    fig=plt.figure(figsize=(WIDTH_MM/25.4,68/25.4))
    for col,metric in enumerate(METRICS):
        ax=fig.add_axes([.165+col*.47,.245,.31,.62])
        panel(ax,"AB"[col],METRIC_TITLES[metric])
        block=data[data.metric.eq(metric)].set_index("method")
        ax.axvline(0,color="#999999",lw=.65,ls="--")
        for y,method in enumerate(METHODS[1:]):
            row=block.loc[method]
            for shift,kind,filled,lw in [(.10,"protein",False,.8),(-.10,"cluster",True,1.3)]:
                lo,hi=row[f"{kind}_ci_low"],row[f"{kind}_ci_high"]
                ax.errorbar(row.mean_delta,y+shift,xerr=[[row.mean_delta-lo],[hi-row.mean_delta]],
                            fmt=MARKERS[method],color=COLORS[method],mfc=COLORS[method] if filled else "white",
                            ms=4.3,elinewidth=lw,capsize=2)
        ax.set(yticks=range(3),yticklabels=[LABELS[m] for m in METHODS[1:]],ylim=(2.6,-.6),
               xlim=(-.04,.21) if col==0 else (-.045,.56),xlabel=f"ΔMAE ({UNITS[metric]})")
        ax.tick_params(axis="y",length=0,pad=4)
        ax.spines["left"].set_visible(False)
    handles=[Line2D([0],[0],marker="o",mfc="white",color="#666666",lw=.8,ms=4,label="Protein bootstrap"),
             Line2D([0],[0],marker="o",color="#666666",lw=1.3,ms=4,label="Connected-group bootstrap")]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.54,.035),ncol=2,fontsize=7)
    export(fig,"fig_ablation_dependence",out)

def captions(out):
    common=("Predictions were evaluated on the same valid residues from five held-out record folds (3,068 records). "
            "Protein-equal summaries include 3,870 exact protein identities for RMSF and 3,867 for Hχ. "
            "MAE denotes mean absolute error; RMSF is root-mean-square fluctuation in ångströms, and Hχ is joint side-chain torsional-state entropy in nats. "
            "ΔMAE = MAE(ablated) − MAE(full), so positive values favour the full model. "
            "All neural models use training seed 7 and a 40-epoch budget. Exact sequences can span folds; these results do not establish strictly unseen-protein generalisation. ")
    text={
        "fig_ablation_main":("Figure 1. Information-source ablations of FlexPPI. (A,B) Protein-equal MAE. "
            "(C,D) Paired MAE differences relative to full FlexPPI. Bars/points are means; whiskers are marginal 95% percentile intervals from 5,000 paired protein resamples. "
            "Identity masked replaces explicit residue identity with an unknown token while retaining structural cues. Intrachain only removes cross-chain messages. "
            "Sequence only removes geometry and partner inputs and uses a smaller scalar-attention architecture (836,676 versus 1,365,812 parameters; 61.26% of full). "+common+
            "Intervals describe test-set variation conditional on the fitted models, not variation between independent training seeds."),
        "fig_ablation_distribution":("Supplementary Figure S1. Distribution of paired protein errors. Every protein pair is retained, including extreme values. "
            "The horizontal axis is symmetric logarithmic, linear within ±0.01 and logarithmic outside that interval. "+common),
        "fig_ablation_heterogeneity":("Supplementary Figure S2. Source heterogeneity and fold consistency. "
            "(A,B) Protein-equal mean paired error differences for all five source-membership categories. "
            "RMSF category sizes are 1,674, 1,536, 640, 12 and 8; Hχ sizes are 1,674, 1,534, 639, 12 and 8, in displayed order. "
            "(C,D) Residue-pooled paired differences in all five complete held-out folds, including sequences spanning folds. "
            "The source and fold panels use different weights. Folds share training data and are not independent repeated experiments. "
            "Source composition and chain context are associated; subgroup patterns are descriptive. "+common),
        "fig_ablation_dependence":("Supplementary Figure S3. Sensitivity to dependence between proteins. "
            "Open symbols show protein-bootstrap intervals; filled symbols show connected-group-bootstrap intervals. "
            "Both estimate the same protein-equal mean using 5,000 paired resamples. "
            "The 2,035 groups merge existing biological split clusters linked by exact sequences; the largest includes 625 proteins. "
            "A sampled group contributes all its proteins, and the estimate divides their summed differences by the sampled protein count. "
            "Full homology clustering was not rerun; grouping does not remove sequence overlap or measure new-seed variability. "+common),
        "fig_ablation_controls":("Supplementary Figure S4. Neural restrictions and feature-probe controls. "
            "The four neural models appear above the separator; the amino-acid prior and four fixed-budget boosting probes appear below it. "
            "Points are protein-equal means and whiskers are marginal 95% protein-bootstrap intervals. "
            "Input representation, architecture and capacity differ across model families, so probe-to-neural differences are descriptive comparisons. "+common),
    }
    (out/"figure_legends.md").write_text("\n\n".join(text.values())+"\n",encoding="utf-8")
    alt={"fig_ablation_main":"Four-panel ablation comparison. Hiding identity sharply increases side-chain entropy error, with little directional change in RMSF. Intrachain and sequence-only models have higher errors for both targets. Paired intervals are shown below the absolute means.",
         "fig_ablation_distribution":"Cumulative distributions of all protein-level paired errors for three restrictions, with a zero reference and an explicitly symmetric logarithmic horizontal scale.",
         "fig_ablation_heterogeneity":"Heatmaps show paired mean differences for all source memberships. Scatter panels show all five residue-pooled fold differences and expose variable RMSF effects of identity masking.",
         "fig_ablation_dependence":"The same paired mean differences are displayed with protein and connected-group resampling intervals. Group intervals are wider; RMSF identity masking continues to span zero.",
         "fig_ablation_controls":"Two forest plots show all nine methods. Neural models and fixed-feature probes are separated visually, with their mean errors and bootstrap intervals."}
    (out/"alt_text.json").write_text(json.dumps(alt,ensure_ascii=False,indent=2),encoding="utf-8")
    return text

def validate_exports(out):
    """Inspect the actual files, including page bounds and TIFF resolution."""
    records=[]
    for path in sorted(out.glob("fig_ablation_*.pdf")):
        with fitz.open(path) as doc:
            page=doc[0]
            width=page.rect.width*25.4/72
            height=page.rect.height*25.4/72
            spans=[s for b in page.get_text("dict")["blocks"] if "lines" in b for line in b["lines"] for s in line["spans"]]
            clipped=[s["text"] for s in spans if s["bbox"][0]<-.1 or s["bbox"][1]<-.1
                     or s["bbox"][2]>page.rect.width+.1 or s["bbox"][3]>page.rect.height+.1]
            embedded=all(len(doc.extract_font(font[0])[3])>0 for font in page.get_fonts())
            min_font=min(s["size"] for s in spans)
            assert abs(width-WIDTH_MM)<.01 and embedded and min_font>=5 and not clipped,(path.name,clipped)
        ET.parse(path.with_suffix(".svg"))
        with Image.open(path.with_suffix(".tif")) as im:
            tiff_dpi=[float(v) for v in im.info["dpi"]]
            size=im.size
            compression=im.info.get("compression")
            assert min(tiff_dpi)>=1199.9
        with Image.open(path.with_suffix(".png")) as im:
            png_dpi=im.info["dpi"]
            assert min(png_dpi)>299.9
        records.append(dict(figure=path.stem,width_mm=width,height_mm=height,
                            minimum_font_pt=min_font,embedded_fonts=embedded,clipped_text=clipped,
                            tiff_dpi=tiff_dpi,tiff_pixels=size,tiff_compression=compression,png_dpi=png_dpi,
                            vector_and_raster_present=True))
    assert len(records)==5
    result={"status":"export_checks_passed","figures":records,"statistics_changed":False,
            "source_rows_retained":7737,"caption_document":"figure_legends.md"}
    (out/"export_QA.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"output/information_ablation_figures/journal_style")
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    paths={k:DATA/v for k,v in {
        "performance":"performance_protein_equal.csv", "contrasts":"paired_ablation_effects.csv",
        "proteins":"source_data_all_proteins.csv", "groups":"source_data_source_groups.csv",
        "folds":"source_data_fold_deltas.csv", "dependence":"source_data_dependence_intervals.csv"}.items()}
    data={k:pd.read_csv(v) for k,v in paths.items()}
    assert len(data["performance"])==18 and len(data["contrasts"])==6 and len(data["proteins"])==7737
    make_main(data["performance"],data["contrasts"],args.output)
    make_distribution(data["proteins"],args.output)
    make_heterogeneity(data["groups"],data["folds"],args.output)
    make_dependence(data["dependence"],args.output)
    make_controls(data["performance"],args.output)
    validate_exports(args.output)
    text=captions(args.output)
    manifest={"profile":"Bioinformatics double-column; compact style also suitable for BiB manuscripts",
              "width_mm":WIDTH_MM,"fonts":"Arial 7–8.5 pt; panel letters 10 pt",
              "vector_formats":["PDF","SVG"],"preview_dpi":300,"tiff_dpi":1200,"tiff_compression":"LZW",
              "statistics_recomputed":False,"source_files":{k:{"path":str(p.resolve()),"sha256":sha256(p)} for k,p in paths.items()},
              "script_sha256":sha256(__file__),"caption_file":"figure_legends.md","alt_text":"alt_text.json",
              "official_sources":["https://academic.oup.com/BIOINFORMATICS/pages/author-guidelines","https://academic.oup.com/bib/pages/msprep_submission"],
              "design_choices":"White background, no figure headline or long in-image footnotes, separate captions, fixed muted model colors, clean axes, compact aligned panels.",
              "statistics_documentation":str((DATA/"statistics_and_reproducibility.json").resolve()),
              "status":"rendered_pending_visual_QA"}
    (args.output/"journal_style_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    lines=["# Bioinformatics / BiB 风格图", "", "以 Bioinformatics 官方双栏宽度 178 mm 为尺寸基准。五张图改为紧凑白底论文面板，去掉整图大标题与图内长段说明，完整统计口径与限制移至独立英文图注。颜色、字号与面板布局属于本次设计选择；期刊并未规定唯一配色模板。现有已核验数据、点估计与置信区间全部复用，没有重算或更改统计结果。", "", "官方依据：[Bioinformatics author guidelines](https://academic.oup.com/BIOINFORMATICS/pages/author-guidelines) 给出 86 / 178 mm 单／双栏宽度、白底及线图高分辨率要求；[BiB manuscript preparation](https://academic.oup.com/bib/pages/msprep_submission) 要求图注独立解释图件，并给出 600 dpi 线图／组合图要求。交付保留矢量 PDF / SVG，另提供 300 dpi PNG 预览与直接渲染的 1200 dpi LZW TIFF。", "", "| 图件 | PDF | SVG | PNG | 1200 dpi TIFF |", "|---|---|---|---|---|"]
    for name in text:
        lines.append("| "+name+" | "+" | ".join(f"[{e.upper()}]({(args.output/(name+'.'+e)).resolve().as_posix()})" for e in ("pdf","svg","png","tif"))+" |")
    lines += ["", "图注："+f"[figure_legends.md]({(args.output/'figure_legends.md').resolve().as_posix()})"+"；统计与数据核验仍见 "+f"[原核验报告]({(DATA/'README_核验与图表说明.md').resolve().as_posix()})"+"。正文建议使用 fig_ablation_main，其余四张作为补充图。", "", "复现：`python -X utf8 scripts/information_sources/make_journal_ablation_figures.py`。"]
    (args.output/"README_期刊风格说明.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print("Saved five journal-style figures to",args.output)

if __name__=="__main__":
    main()
