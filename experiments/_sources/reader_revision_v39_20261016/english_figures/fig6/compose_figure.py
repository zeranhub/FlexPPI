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

"""Compose current-model case evidence; no synthetic points or conformations.

Layout inherits the existing KIX statewise publication figure's paired ranking
maps and shared structure cameras, with the verified fold-1 predictions. The
line-chart styling follows the skill LineTrend asset; scientific data/axes are
protein residue numbers and observed joint-angle occupancy, not example data.
The final canvas is explicit: 183 x 225 mm, avoiding tight-bbox size changes.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/v13_case_refresh_20261009"
BUNDLE = OUT / "figure_bundle"
if (Path(__file__).resolve().parent / "kix").is_dir():
    BUNDLE = Path(__file__).resolve().parent
PROJECT = Path(r"E:\Project\FlexPPI")
BLUE, ORANGE, PURPLE = CATEGORICAL[0], CATEGORICAL[3], CATEGORICAL[4]
MD = "#666666"
W, H = 183., 225.

def ax_mm(fig, x, y, width, height):
    return fig.add_axes([x/W, y/H, width/W, height/H])

def txt(fig, x, y, text, **kwargs):
    kwargs.setdefault("color", BLACK)
    return fig.text(x/W, y/H, text, **kwargs)

def picture(fig, path, x, y, width, height):
    assert path.is_file(), path
    ax = ax_mm(fig, x, y, width, height)
    ax.imshow(plt.imread(path))
    ax.set_axis_off()
    return ax

def source_data():
    data_dir = BUNDLE / "data"
    data_dir.mkdir(exist_ok=True, parents=True)
    portable_files = ["kix_fold1_all_residues.csv", "kix_fold1_case_metrics.csv",
        "ubiquitin_fold0_all_residues.csv", "ubiquitin_Q41_Q62_all27_joint_frequencies.csv",
        "ubiquitin_Q41_Q62_validation_summary.csv"]
    if (BUNDLE/"case_model_protocol.json").is_file():
        protocol=json.loads((BUNDLE/"case_model_protocol.json").read_text(encoding="utf-8"))
        assert protocol["training_protocol"]=="chemical_identity_corrected_v1"
        missing=[name for name in portable_files if not (data_dir/name).is_file()]
        if missing:
            raise FileNotFoundError("Missing authoritative v13 source tables: "+", ".join(missing))
    if all((data_dir/name).is_file() for name in portable_files):
        tables = [pd.read_csv(data_dir/name) for name in portable_files]
        assert len(tables[0].loc[tables[0].metric.eq("rmsf")]) == 4*87
        tables[2] = tables[2].loc[tables[2].chain.eq("A")].sort_values("residue_id").copy()
        assert int(tables[2].valid_rmsf.sum()) == 72
        assert int(tables[2].valid_hchi.sum()) == 66
        assert len(tables[3]) == 54
        assert list(tables[4].states_covering_90_percent) == [7,12]
        return tables
    original = PROJECT / "artifacts/kix_md_20260914/pilot_results_20260916/comparison/matched_residues.csv"
    matched = pd.read_csv(original, dtype={"model_fold":str})
    kix = matched[matched.model_fold.eq("1")].copy()
    states = ["apo", "myb", "mll", "both"]
    rank_rows, metric_rows = [], []
    for state in states:
        for metric in ["rmsf", "hchi"]:
            sub = kix[kix.state.eq(state) & kix.metric.eq(metric)].sort_values("author_residue_number").copy()
            assert len(sub) == (87 if metric == "rmsf" else 74)
            for kind in ["md", "model"]:
                sub[kind+"_percentile"] = 100*(sub[kind+"_value"].rank(method="average")-1)/(len(sub)-1)
            rank_rows.append(sub)
            n_top = int(np.ceil(len(sub)*.2))
            md_top = set(sub.nlargest(n_top, "md_value").author_residue_number)
            pred_top = set(sub.nlargest(n_top, "model_value").author_residue_number)
            loop = sub[sub.author_residue_number.between(614,621)]
            metric_rows.append(dict(state=state, metric=metric, fold=1, n=len(sub),
                spearman=float(spearmanr(sub.md_value, sub.model_value).statistic),
                mae=float(np.abs(sub.md_value-sub.model_value).mean()), top_n=n_top,
                top_overlap=len(md_top & pred_top), l12_n=len(loop),
                l12_md_mean_percentile=float(loop.md_percentile.mean()),
                l12_model_mean_percentile=float(loop.model_percentile.mean())))
    ranks = pd.concat(rank_rows, ignore_index=True)
    ranks.to_csv(data_dir/"kix_fold1_all_residues.csv",index=False)
    kix_metrics = pd.DataFrame(metric_rows)
    kix_metrics.to_csv(data_dir/"kix_fold1_case_metrics.csv",index=False)
    all_res = pd.read_csv(OUT/"named_case_residue_data.csv")
    ubiquitin = all_res[all_res.sample_id.eq("A003R") & all_res.chain.eq("A")].sort_values("residue_id").copy()
    assert int(ubiquitin.valid_rmsf.sum()) == 71
    assert int(ubiquitin.valid_hchi.sum()) == 58
    ubiquitin.to_csv(data_dir/"ubiquitin_fold0_all_residues.csv",index=False)
    frequencies = pd.read_csv(OUT/"rotamers/A003R_Q41_Q62_joint_frequencies_5000_label_frames.csv")
    assert len(frequencies) == 54 and np.allclose(frequencies.groupby("residue").probability.sum(),1.)
    frequencies.to_csv(data_dir/"ubiquitin_Q41_Q62_all27_joint_frequencies.csv",index=False)
    summary = pd.read_csv(OUT/"rotamers/A003R_Q41_Q62_validation_summary.csv")
    summary.to_csv(data_dir/"ubiquitin_Q41_Q62_validation_summary.csv",index=False)
    chain_metrics = pd.read_csv(OUT/"all_chain_case_metrics.csv")
    chain_metrics[chain_metrics.sample_id.eq("A003R")].to_csv(data_dir/"ubiquitin_UBA_fold0_chain_metrics.csv",index=False)
    for name in ["A003R_Q41_Q62_joint_frequencies_per_window.csv", "A003R_Q41_Q62_label_windows.csv", "A003R_Q41_Q62_provenance_and_QA.json"]:
        shutil.copy2(OUT/"rotamers"/name, data_dir/name)
    assert list(kix_metrics[kix_metrics.metric.eq("rmsf")].top_overlap)==[13,15,15,14]
    assert list(summary.states_covering_90_percent)==[7,12]
    return ranks,kix_metrics,ubiquitin,frequencies,summary

def compose():
    ranks, metrics, ub, frequencies, summary = source_data()
    fig = plt.figure(figsize=(W/25.4,H/25.4), facecolor="white")
    txt(fig,1,221,"KIX | Locating a flexible connecting loop",fontsize=9,fontweight="bold")
    txt(fig,1,214,"a",fontsize=10,fontweight="bold")
    txt(fig,7,214,"Ternary complex",fontsize=8)
    picture(fig,BUNDLE/"kix/kix_context.png",1,169,54,43)
    txt(fig,62,214,"b",fontsize=10,fontweight="bold")
    txt(fig,68,214,"L12 close-up (614–621), Both state",fontsize=8)
    picture(fig,BUNDLE/"kix/kix_l12_md.png",61,169,56,42)
    picture(fig,BUNDLE/"kix/kix_l12_fold1.png",122,169,56,42)
    txt(fig,88,210,"MD",ha="center",fontsize=8)
    txt(fig,150,210,"FlexPPI",ha="center",fontsize=8)
    txt(fig,1,164,"c",fontsize=10,fontweight="bold")
    txt(fig,7,164,"Full KIX backbone mobility",fontsize=8)
    states = ["apo","myb","mll","both"]
    names = ["Apo","cMyb","MLL","Both"]
    matrix, labels = [], []
    for state,name in zip(states,names):
        d = ranks[ranks.state.eq(state)&ranks.metric.eq("rmsf")].sort_values("author_residue_number")
        for key,label in [("md_percentile","MD"),("model_percentile","model")]:
            matrix.append(d[key].values)
            labels.append(name+" "+label)
    cmap = LinearSegmentedColormap.from_list("mobility_rank",SEQUENTIAL)
    ax = ax_mm(fig,23,129,83,30)
    im = ax.imshow(matrix, origin="upper", aspect="auto",cmap=cmap,vmin=0,vmax=100,
                   extent=(585.5,672.5,7.5,-.5),interpolation="nearest")
    ax.set_yticks(range(8));ax.set_yticklabels(labels,fontsize=7)
    ax.set_xticks([586,614,621,646,672]);ax.set_xlabel("KIX residue number",labelpad=2)
    ax.tick_params(length=2,pad=2)
    for y in [1.5,3.5,5.5]:ax.axhline(y,color="white",lw=1.1)
    ax.add_patch(Rectangle((613.5,-.5),8,8,fill=False,lw=.9,ec=ORANGE))
    for spine in ax.spines.values():spine.set_visible(False)
    cbax = ax_mm(fig,113,134,2.2,23)
    cb = fig.colorbar(im,cax=cbax,ticks=[0,50,100])
    cb.ax.tick_params(labelsize=7,length=2,pad=2)
    cb.set_label("RMSF percentile",fontsize=7,labelpad=2)
    txt(fig,130,164,"d",fontsize=10,fontweight="bold")
    txt(fig,136,164,"L12 localization",fontsize=8)
    dm = metrics[metrics.metric.eq("rmsf")].set_index("state").loc[states]
    ax = ax_mm(fig,140,136,38,23)
    ax.axhspan(80,100,color=ORANGE,alpha=.12,lw=0)
    ax.plot(range(4),dm.l12_md_mean_percentile,"o-",color=MD,ms=3,lw=1,label="MD")
    ax.plot(range(4),dm.l12_model_mean_percentile,"s--",color=BLUE,ms=3,lw=1,label="model")
    lowest=min(dm.l12_model_mean_percentile.min(),dm.l12_md_mean_percentile.min())
    ax.set_ylim(min(75,5*np.floor((lowest-2)/5)),100);ax.set_yticks([80,90,100]);ax.set_xlim(-.3,3.3)
    ax.set_xticks(range(4));ax.set_xticklabels(names,fontsize=7)
    ax.set_ylabel("L12 mean percentile",fontsize=7,labelpad=1)
    ax.tick_params(length=2,pad=2)
    ax.legend(loc="upper left",bbox_to_anchor=(0,1.01),borderaxespad=.1,
              ncol=2,fontsize=7,handlelength=1,columnspacing=.5)
    txt(fig,133,127,"Top-18 overlap (MD ∩ model)",fontsize=7)
    overlaps="   ".join(f"{int(row.top_overlap)}/{int(row.top_n)}" for _,row in dm.iterrows())
    txt(fig,139,122,overlaps,fontsize=7)
    fig.add_artist(mpl.lines.Line2D([1/W,182/W],[114/H,114/H],transform=fig.transFigure,color="#D9D9D9",lw=.6))
    txt(fig,1,109,"Ubiquitin–UBA | Distinguishing side-chain mobility",fontsize=9,fontweight="bold")
    txt(fig,1,102,"e",fontsize=10,fontweight="bold")
    txt(fig,7,102,"2QHO complex",fontsize=8)
    picture(fig,BUNDLE/"ubiquitin/ubiquitin_context.png",1,60,54,40)
    txt(fig,62,102,"f",fontsize=10,fontweight="bold")
    txt(fig,68,102,"Observed MD side chains and backbone context",fontsize=8)
    picture(fig,BUNDLE/"ubiquitin/q41_zoom.png",61,60,56,40)
    picture(fig,BUNDLE/"ubiquitin/q62_zoom.png",122,60,56,40)
    txt(fig,88,99,"Q41",ha="center",fontsize=8,color=ORANGE)
    txt(fig,150,99,"Q62",ha="center",fontsize=8,color=PURPLE)
    txt(fig,1,56,"g",fontsize=10,fontweight="bold")
    txt(fig,7,56,"Full ubiquitin residue profiles",fontsize=8)
    for key,y,height,label,yticks in [("rmsf",33,17,"RMSF (Å)",[0,.5,1]),("hchi",10,16,r"$H_\chi$ (nats)",[0,1.5,3])]:
        ax = ax_mm(fig,15,y,82,height)
        valid_rows=ub.loc[ub["valid_"+key].eq(True)].sort_values("residue_id")
        assert len(valid_rows)==(72 if key=="rmsf" else 66)
        assert np.isfinite(valid_rows[["md_"+key,"pred_"+key]].values).all()
        # Lines join only observed valid positions in residue-number order.
        # No value is generated at an unavailable position; faint axis rugs
        # identify those positions separately from the numerical observations.
        for source,color,name in [("md",MD,"MD"),("pred",BLUE,"FlexPPI")]:
            ax.plot(valid_rows.residue_id,valid_rows[source+"_"+key],
                color=color,lw=.9,ls="-",label=name,
                marker="o" if source=="md" else None,
                markersize=1.7,markerfacecolor="white",markeredgewidth=.45,
                zorder=3 if source=="md" else 2)
        unavailable=sorted(set(range(1,73))-set(valid_rows.residue_id))
        ax.plot(unavailable,[-.035]*len(unavailable),linestyle="none",marker="|",
            markersize=2.7,markeredgewidth=.6,color="#B8B8B8",
            transform=ax.get_xaxis_transform(),clip_on=False)
        for res,color in [(41,ORANGE),(62,PURPLE)]:ax.axvline(res,color=color,lw=.6,ls=":")
        ax.set_xlim(1,72);ax.set_yticks(yticks);ax.set_ylim(bottom=0)
        ax.set_ylabel(label,fontsize=7,labelpad=1);ax.tick_params(length=2,pad=2)
        ax.set_xticks([1,20,41,62,72])
        if key=="rmsf":
            ax.set_xticklabels([])
            ax.legend(loc="lower left",bbox_to_anchor=(0,1.01),ncol=2,fontsize=7,handlelength=1.4,columnspacing=.8)
            rho=float(spearmanr(valid_rows.md_rmsf,valid_rows.pred_rmsf).statistic)
            ax.text(.98,.87,rf"$\rho={rho:.3f}$",transform=ax.transAxes,ha="right",fontsize=7)
        else:
            ax.set_xlabel("Ubiquitin residue number",fontsize=8,labelpad=2)
            rho=float(spearmanr(valid_rows.md_hchi,valid_rows.pred_hchi).statistic)
            ax.text(.98,.05,rf"$\rho={rho:.3f}$",transform=ax.transAxes,ha="right",fontsize=7)
    txt(fig,110,56,"h",fontsize=10,fontweight="bold")
    txt(fig,116,56,"Observed side-chain state coverage",fontsize=8)
    ax=ax_mm(fig,120,18,57,32)
    ax.axhline(90,color=GREY,ls=":",lw=.8)
    ax.text(27,88,"90%",ha="right",va="top",fontsize=7,color=MD)
    for res,color,marker in [("Q41",ORANGE,"o"),("Q62",PURPLE,"s")]:
        freq=frequencies[frequencies.residue.eq(res)].sort_values(["probability","joint_code"],ascending=[False,True])
        cumulative=freq.probability.cumsum().values*100
        n90=int(np.searchsorted(cumulative,90)+1)
        ax.plot(range(1,28),cumulative,color=color,lw=1.3,marker=marker,ms=2.2,label=res)
        ax.plot([n90,n90],[0,cumulative[n90-1]],color=color,ls="--",lw=.6)
        ax.text(n90+.4,66 if res=="Q41" else 37,f"{res}: {n90} states",fontsize=7,color=color)
    ax.set_xlim(1,27);ax.set_ylim(0,104);ax.set_xticks([1,7,12,20,27]);ax.set_yticks([0,50,100])
    ax.set_xlabel("Joint-angle classes, ranked by occupancy",fontsize=7,labelpad=2)
    ax.set_ylabel("MD frames covered (%)",fontsize=7,labelpad=2)
    ax.tick_params(length=2,pad=2)
    q41=summary.loc[summary.residue.eq("Q41")].iloc[0]
    q62=summary.loc[summary.residue.eq("Q62")].iloc[0]
    txt(fig,113,7,rf"$H_\chi$ MD / model: Q41 {q41.archived_md_hchi:.2f} / {q41.model_hchi:.2f}",fontsize=7)
    txt(fig,139,3,f"Q62 {q62.archived_md_hchi:.2f} / {q62.model_hchi:.2f} nats",fontsize=7)
    destination=BUNDLE/"fig_local_dynamic_cases"
    # Fixed page dimensions preserve the figure's intended print font sizes.
    with mpl.rc_context({"savefig.bbox": None}):
        for extension in ["pdf","svg","png"]:
            fig.savefig(str(destination)+"."+extension,bbox_inches=None,dpi=300)
    plt.close(fig)
    sources=[]
    for path in sorted((BUNDLE/"data").glob("*")):
        sources.append(dict(file=str(path.relative_to(BUNDLE)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    qa=dict(canvas_mm=[W,H], dpi=300, kix_fold=1,ubiquitin_fold=0,
        kix_states=states, kix_all_rmsf_positions_per_state=87,
        ubiquitin_rmsf_positions=int(ub.valid_rmsf.sum()),ubiquitin_hchi_positions=int(ub.valid_hchi.sum()),
        model_protocol="v13 chemical_identity_corrected_v1 full seed7 epoch40; KIX fixed fold1 and ubiquitin OOF fold0",
        rotamer_classes_per_residue=27,structure_rendering="PyMOL ray tracing",
        local_structure_rule="Complete protein cartoon rendered in full frame, then exact pixel crop; full and focus scenes preserve all chains",
        local_crop_manifests=["kix/crop_manifest.json", "ubiquitin/crop_manifest.json"],
        complete_cartoon_CA_counts={"KIX_ternary":131,"ubiquitin_UBA":119},
        side_chain_display="Complete real heavy-atom bonds; each MD side chain has its own whole-ubiquitin CA trace",
        structure_style_reference="Local OPUS-Rota4 Figure 4, IsAb2.0 Figure 3 and Sequeiros-Borja et al. Figure 3A",
        conformation_source="existing MD frames; model predicts scalar mobility only",
        panel_g_line_rule="Solid lines connect available observations at original residue numbers; no imputation",
        panel_g_unavailable_positions={key:sorted(set(range(1,73))-set(ub.loc[ub["valid_"+key].eq(True),"residue_id"])) for key in ["rmsf","hchi"]},
        panel_g_missing_label_marker="faint short ticks below each numerical axis",
        sources=sources)
    (BUNDLE/"figure_provenance_and_QA.json").write_text(json.dumps(qa,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(destination)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Recompose figure from bundled source tables and PyMOL renders.")
    parser.add_argument("--bundle-dir", type=Path)
    args = parser.parse_args()
    if args.bundle_dir:
        BUNDLE = args.bundle_dir.resolve()
    compose()
