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

# This is a reviewable Matplotlib draft, not completed EditaPlot/Origin work.
# User typography and exact physical size supersede the generic baseline above.
from pathlib import Path
import json
import hashlib
import sys
import numpy as np
import pandas as pd
mpl.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data"
OUT = ROOT / "design_draft_preview"
FULL = "#2166AC"
INTRA = "#9A633D"
MD = "#222222"
MIDGREY = "#888888"
SOURCE_ORDER = ["atlas", "dynarepo", "itrace"]
SOURCE_MARKERS = {"atlas": "o", "dynarepo": "s", "itrace": "^"}

mpl.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial"],
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "figure.titlesize": 10, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "text.usetex": False, "axes.unicode_minus": True,
    "savefig.bbox": None,
})

font_path = font_manager.findfont(font_manager.FontProperties(family="Arial"), fallback_to_default=False)
cam = pd.read_csv(DATA / "cam_all25.csv").sort_values("CaMKI_native_position")
ub_all = pd.read_csv(DATA / "ubiquitin_72_position_profile_with_original_mask.csv")
ub = ub_all.loc[ub_all["valid_hchi"].astype(str).str.lower().eq("true")].copy()
gln = pd.read_csv(DATA / "ubiquitin_all6_Gln_Hchi_SASA_plot.csv").sort_values("MD_Hchi_nats")
angles = pd.read_csv(DATA / "defined_angles_long.csv")
chi_stats = pd.read_csv(DATA / "per_chi_correlations.csv")
correlations = pd.read_csv(DATA / "retained_absolute_correlations_all498.csv")
assert len(cam) == 25 and list(cam["CaMKI_native_position"]) == list(range(294, 319))
assert len(ub_all) == 72 and len(ub) == 66 and len(gln) == 6
assert [int((angles.chi == f"chi{i}").sum()) for i in range(1,5)] == [66,52,23,10]
assert len(angles) == 151 and not angles[["predicted_Hchi_nats_AA_centered", "fine_entropy_four_origins_mean_nats_AA_centered"]].isna().any().any()

# Spearman rank correlation without a scipy dependency: average tied ranks.
def spearman(a, b):
    return float(pd.Series(a).rank().corr(pd.Series(b).rank()))

# Pair records before resampling. Preserve original source sample sizes (30,23,30).
pairs = {}
paired_statistics = []
for metric in ["fine_entropy", "lowF_width"]:
    df = correlations.loc[correlations.metric_key.eq(metric)]
    wide = df.pivot(index=["dataset", "sample_id"], columns="predictor", values="rho_type_centered").reset_index()
    wide = wide.sort_values(["dataset", "sample_id"]).reset_index(drop=True)
    assert len(wide) == 83 and wide.groupby("dataset").size().to_dict() == {"atlas":30,"dynarepo":23,"itrace":30}
    wide["delta_rho"] = wide["full_hchi"] - wide["full_rmsf"]
    assert np.isfinite(wide[["full_hchi", "full_rmsf", "delta_rho"]].to_numpy()).all()
    pairs[metric] = wide
    rng = np.random.default_rng(20261026)
    boot = np.zeros(10000)
    for source in SOURCE_ORDER:
        arr = wide.loc[wide.dataset.eq(source), "delta_rho"].to_numpy()
        boot += arr[rng.integers(0, len(arr), size=(10000,len(arr)))].sum(axis=1)
    boot /= len(wide)
    lo, hi = np.quantile(boot,[0.025,0.975])
    paired_statistics.append({"metric_key":metric,"n_records":83,"mean":float(wide.delta_rho.mean()),
                              "bootstrap_95_low":float(lo),"bootstrap_95_high":float(hi),
                              "Full_Hchi_mean":float(wide.full_hchi.mean()),"Full_RMSF_mean":float(wide.full_rmsf.mean()),
                              "positive_records":int(wide.delta_rho.gt(0).sum()),
                              "bootstrap_replicates":10000,"bootstrap_seed":20261026})
stat_f = {r["metric_key"]:r for r in paired_statistics}

fig = plt.figure(figsize=(183/25.4, 210/25.4), facecolor="white")
fig.text(0.035, 0.986, "DRAFT DESIGN PREVIEW", fontsize=10, weight="bold", va="top", color="#555555")
fig.text(0.965, 0.986, "Figure 3 · quantitative panels", fontsize=9, ha="right", va="top", color="#555555")

def heading(letter, title, x, y):
    fig.text(x, y, letter, fontsize=11, weight="bold", va="bottom")
    fig.text(x+0.035, y, title, fontsize=10, va="bottom")

def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(pad=3)
    ax.grid(False)

# b: contiguous sequence, absolute RMSF. All 25 positions for all three series.
ax_b = fig.add_axes([0.100, 0.724, 0.385, 0.207])
heading("b", "CaMKI backbone constraint", 0.045, 0.945)
ax_b.axvspan(302.5, 316.5, color="#EAEAEA", zorder=0)
x = cam.CaMKI_native_position.to_numpy()
ax_b.plot(x, cam.target_rmsf, color=MD, lw=0.9, ls=":", marker="o", markersize=3.1,
          markerfacecolor="white", markeredgewidth=0.7, label="MD", zorder=3)
ax_b.plot(x, cam.pred_rmsf, color=FULL, lw=1.25, ls="-", label="Full", zorder=4)
ax_b.plot(x, cam.intra_rmsf, color=INTRA, lw=1.05, ls="--", label="Intra", zorder=2)
ax_b.set(xlim=(293.5,318.5), ylim=(0,8.15), xticks=[294,300,306,312,318], yticks=[0,2,4,6,8],
         xlabel="CaMKI residue", ylabel="RMSF (Å)")
ax_b.legend(loc="upper left", handlelength=1.7, handletextpad=0.5, borderaxespad=0.25, labelspacing=0.3)
mae_full = float(np.mean(np.abs(cam.pred_rmsf-cam.target_rmsf)))
mae_intra = float(np.mean(np.abs(cam.intra_rmsf-cam.target_rmsf)))
ax_b.text(0.37,0.96,f"MAE (Å)\nFull {mae_full:.3f}\nIntra {mae_intra:.3f}", transform=ax_b.transAxes, va="top", fontsize=9)
ax_b.text(0.63,0.62,"W303–M316", transform=ax_b.transAxes, ha="center", fontsize=9, color="#555555")
clean(ax_b)

# c: equality scatter answers whole-chain agreement, not sequence continuity.
ax_c = fig.add_axes([0.655, 0.719, 0.263, 0.218])
heading("c", "Ubiquitin agreement", 0.555, 0.945)
lim=(-0.08,4.45)
ax_c.plot(lim,lim,color=MIDGREY,lw=0.7,ls="--",zorder=0)
ax_c.scatter(ub.MD_Hchi_nats,ub.Full_Hchi_nats,s=15,c=FULL,alpha=0.72,linewidths=0,zorder=2)
ax_c.set(xlim=lim,ylim=lim,xticks=[0,1,2,3,4],yticks=[0,1,2,3,4],
         xlabel="MD Hχ (nats)",ylabel="Full Hχ (nats)")
ax_c.set_aspect("equal",adjustable="box")
mae_ub = float(np.mean(np.abs(ub.Full_Hchi_nats-ub.MD_Hchi_nats)))
rho_ub = spearman(ub.Full_Hchi_nats,ub.MD_Hchi_nats)
ax_c.text(0.04,0.98,f"n = 66\nMAE = {mae_ub:.3f}\nρ = {rho_ub:.3f}",transform=ax_c.transAxes,va="top",fontsize=9)
for resid,offset in [(40,(14,2)),(41,(-29,-13)),(62,(12,-12))]:
    r=ub.loc[ub.residue_id.eq(resid)].iloc[0]
    ax_c.scatter([r.MD_Hchi_nats],[r.Full_Hchi_nats],s=23,facecolors="white",edgecolors=FULL,lw=0.8,zorder=3)
    ax_c.annotate(f"Q{resid}",(r.MD_Hchi_nats,r.Full_Hchi_nats),xytext=offset,textcoords="offset points",fontsize=9,
                  arrowprops=dict(arrowstyle="-",lw=0.5,color=MIDGREY,shrinkA=1,shrinkB=3))
clean(ax_c)

# d: all six Gln, ordered by observed Hchi; SASA is an explanatory column, no fit.
ax_d = fig.add_axes([0.100, 0.433, 0.305, 0.203])
heading("d", "Same-type side-chain environments", 0.045, 0.655)
y=np.arange(len(gln))
for j,(_,r) in enumerate(gln.iterrows()):
    ax_d.plot([r.MD_Hchi_nats,r.Full_Hchi_nats],[j,j],color=MIDGREY,lw=0.8,zorder=1)
ax_d.scatter(gln.MD_Hchi_nats,y,s=24,facecolors="white",edgecolors=MD,lw=0.8,label="MD",zorder=3)
ax_d.scatter(gln.Full_Hchi_nats,y,s=22,c=FULL,linewidths=0,label="Full",zorder=4)
ax_d.set(xlim=(1.75,3.16),ylim=(5.55,-0.55),yticks=y,yticklabels=gln.residue.tolist(),
         xticks=[1.8,2.2,2.6,3.0],xlabel="Joint Hχ (nats)")
for tick in ax_d.get_yticklabels():
    if tick.get_text() in {"Q41","Q62"}:tick.set_weight("bold")
for j,(_,r) in enumerate(gln.iterrows()):
    ax_d.text(1.13,j,f"{r.sidechain_SASA_A2:.1f}",transform=ax_d.get_yaxis_transform(),ha="center",va="center",fontsize=9)
ax_d.text(1.13,1.04,"SASA\n(Å²)",transform=ax_d.transAxes,ha="center",va="bottom",fontsize=9)
ax_d.legend(loc="lower left",bbox_to_anchor=(-0.04,0.97),ncol=2,handletextpad=0.35,columnspacing=0.8,borderpad=0,handlelength=1.2)
clean(ax_d)

# f: same-record paired deltas make the comparison and analysis unit explicit.
ax_f = fig.add_axes([0.625, 0.433, 0.320, 0.203])
heading("f", "Side-chain occupancy specificity", 0.555, 0.655)
for k,metric in enumerate(["fine_entropy","lowF_width"]):
    w=pairs[metric]
    for source in SOURCE_ORDER:
        a=w.loc[w.dataset.eq(source)].copy()
        # Deterministic spread is display-only; every record is retained once.
        jitter=np.linspace(-0.20,0.20,len(a))
        ax_f.scatter(k+jitter,a.delta_rho,s=13,marker=SOURCE_MARKERS[source],c=FULL,alpha=0.28,linewidths=0,zorder=2)
    s=stat_f[metric]
    ax_f.errorbar(k,s["mean"],yerr=[[s["mean"]-s["bootstrap_95_low"]],[s["bootstrap_95_high"]-s["mean"]]],
                  fmt="o",color=FULL,markerfacecolor=FULL,markersize=5.0,ecolor="#666666",elinewidth=1.05,capsize=3,zorder=5)
    ax_f.text(k,0.556,f"{s['Full_Hchi_mean']:.3f} vs {s['Full_RMSF_mean']:.3f}",ha="center",va="top",fontsize=9,color="#444444")
ax_f.axhline(0,color=MIDGREY,lw=0.7,ls="--",zorder=1)
ax_f.set(xlim=(-0.5,1.5),ylim=(-0.25,0.58),xticks=[0,1],xticklabels=["Fine-bin\nentropy","Low-F* angular\nwidth"],
         yticks=[-0.2,0,0.2,0.4],ylabel="Δρ: Full Hχ − Full RMSF")

clean(ax_f)
# Source glyphs are keyed without assigning the Intra brown to another method.
handles=[Line2D([],[],ls="",marker=SOURCE_MARKERS[s],color=FULL,alpha=0.55,markersize=3.4,label=l)
         for s,l in zip(SOURCE_ORDER,["ATLAS 30","DynaRepo 23","ITrace 30"])]
fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(0.782,0.365),ncol=3,fontsize=9,
           handletextpad=0.2,columnspacing=0.6,borderpad=0)

# e: one shared joint prediction, four separately type-centered single-angle readouts.
heading("e", "Joint prediction and single-angle occupancy", 0.045, 0.318)
axes_e=[]
for i in range(4):
    ax=fig.add_axes([0.100+i*0.220,0.104,0.175,0.177])
    axes_e.append(ax)
    a=angles.loc[angles.chi.eq(f"chi{i+1}")].copy()
    # Ordinary points are drawn first; the two highlighted rows replace their marks.
    highlighted=a.residue.isin(["Q41","Q2"])
    ax.scatter(a.loc[~highlighted,"predicted_Hchi_nats_AA_centered"],
               a.loc[~highlighted,"fine_entropy_four_origins_mean_nats_AA_centered"],
               s=14,c=FULL,alpha=0.55,linewidths=0,zorder=2)
    for residue,marker in [("Q41","D"),("Q2","s")]:
        r=a.loc[a.residue.eq(residue)]
        if len(r):
            ax.scatter(r.predicted_Hchi_nats_AA_centered,r.fine_entropy_four_origins_mean_nats_AA_centered,
                       s=27,marker=marker,facecolors="white",edgecolors=MD,lw=0.8,zorder=4)
    ax.axhline(0,color="#BBBBBB",lw=0.5,zorder=0)
    ax.axvline(0,color="#BBBBBB",lw=0.5,zorder=0)
    ax.set(xlim=(-0.75,0.55),ylim=(-1.05,0.80),xticks=[-0.5,0,0.5],yticks=[-1,-0.5,0,0.5])
    ax.set_title(f"χ{i+1}",fontsize=10,pad=4)
    s=chi_stats.loc[chi_stats.chi.eq(f"chi{i+1}") & chi_stats.predictor.eq("predicted_joint_Hchi") & chi_stats.metric.eq("fine_entropy")].iloc[0]
    ax.text(0.04,0.96,f"ρ = {s.rho_AA_centered:.2f}\nn = {len(a)}",transform=ax.transAxes,fontsize=9,va="top")
    if i>0:ax.tick_params(labelleft=False)
    clean(ax)
fig.text(0.028,0.192,"MD single-χ entropy\nrelative to AA mean (nats/χ)",rotation=90,ha="center",va="center",fontsize=9)
fig.text(0.51,0.057,"Joint Hχ prediction relative to AA mean (nats)",ha="center",fontsize=9)
# Explicit selection legend, all three relevant angles carry the highlighted points.
fig.legend(handles=[Line2D([],[],ls="",marker="D",markerfacecolor="white",markeredgecolor=MD,markersize=4,label="Q41: lowest predicted Gln"),
                    Line2D([],[],ls="",marker="s",markerfacecolor="white",markeredgecolor=MD,markersize=4,label="Q2: highest predicted Gln")],
           loc="lower center",bbox_to_anchor=(0.54,0.018),ncol=2,handletextpad=0.45,columnspacing=1.5,borderpad=0,fontsize=9)
fig.text(0.035,0.003,"Matplotlib draft preview · Origin project and formal readback pending",fontsize=9,color="#555555",va="bottom")

# Guarantee text floors and user-specified font on every Matplotlib Text object.
fig.canvas.draw()
for t in fig.findobj(match=mpl.text.Text):
    if t.get_text():
        assert t.get_fontsize() >= 9, (t.get_text(), t.get_fontsize())
        t.set_fontfamily("Arial")
for ext in ["pdf","png","svg"]:
    meta={"Title":"Figure 3 quantitative panels — draft design preview"}
    if ext=="pdf":meta["Subject"]="Matplotlib design preview only; not verified Origin completion"
    fig.savefig(OUT.with_suffix("."+ext),dpi=300,bbox_inches=None,metadata=meta)

# Precise geometry/text/readout audit supplied beside the review assets.
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
text_audit=[]
for t in fig.findobj(match=mpl.text.Text):
    if t.get_text() and t.get_visible():
        bbox=t.get_window_extent(renderer)
        text_audit.append({"text":t.get_text(),"fontsize_pt":t.get_fontsize(),"fontfamily":t.get_fontfamily(),
                           "bbox_px":[float(x) for x in bbox.extents],
                           "outside_canvas":bool(bbox.x0 < -1 or bbox.y0 < -1 or bbox.x1 > fig.bbox.x1+1 or bbox.y1 > fig.bbox.y1+1)})
file_hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in DATA.glob("*.csv")}
audit={"artifact_status":"draft_preview_not_Origin_completion","font_path":font_path,
       "width_mm":183,"height_mm":210,"raster_dpi":300,"raster_expected_px":[round(183/25.4*300),round(210/25.4*300)],
       "panels":{"b_positions":25,"b_rendered_points":75,"c_valid_positions":66,"c_original_inventory":72,
                 "c_no_chi_positions_unplotted":6,"d_Gln_positions":6,"d_rendered_pairs":6,
                 "e_points_per_chi":[66,52,23,10],"e_total_points":151,"f_records_per_metric":83,"f_delta_points":166},
       "MAE_CaMKI_Full_A":mae_full,"MAE_CaMKI_Intra_A":mae_intra,"MAE_ubiquitin_Full_nats":mae_ub,"rho_ubiquitin":rho_ub,
       "paired_delta_statistics":paired_statistics,"source_sha256":file_hashes,"text":text_audit,
       "caveats":["Not a verified Origin route or editable OPJU.","No structure images or original panel designs were opened.",
                  "Six no-chi ubiquitin inventory rows retained in source; not assigned entropy zero.",
                  "No regression, smoothing, outlier removal, MD resimulation, or model refitting.",
                  "All supplied plotted analysis units are retained; highlighting replaces rather than duplicates angle points."]}
(ROOT/"draft_preview_audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding="utf-8")
pd.DataFrame(paired_statistics).to_csv(ROOT/"draft_preview_f_paired_statistics.csv",index=False,encoding="utf-8-sig")
plt.close(fig)

# Supplemental review asset: same 25 positions at two non-truncating scales.
fig_b = plt.figure(figsize=(89/25.4, 140/25.4), facecolor="white")
fig_b.text(0.07,0.985,"DRAFT · panel b detail",fontsize=10,weight="bold",va="top",color="#555555")
for j,(bottom,upper,title) in enumerate([(0.580,8.15,"All model amplitudes"),(0.145,2.30,"MD and Full · all 25 positions")]):
    ax=fig_b.add_axes([0.19,bottom,0.755,0.300])
    ax.axvspan(302.5,316.5,color="#EAEAEA",zorder=0)
    ax.plot(x,cam.target_rmsf,color=MD,lw=0.9,ls=":",marker="o",markersize=3.0,
            markerfacecolor="white",markeredgewidth=0.7,label="MD",zorder=3)
    ax.plot(x,cam.pred_rmsf,color=FULL,lw=1.25,label="Full",zorder=4)
    if j==0:ax.plot(x,cam.intra_rmsf,color=INTRA,lw=1.05,ls="--",label="Intra",zorder=2)
    ax.set(xlim=(293.5,318.5),ylim=(0,upper),xticks=[294,300,306,312,318],
           xlabel="CaMKI residue",ylabel="RMSF (Å)")
    ax.set_yticks([0,2,4,6,8] if j==0 else [0,0.5,1,1.5,2])
    ax.set_title(title,fontsize=10,loc="left",pad=7)
    clean(ax)
    if j==0:ax.legend(loc="upper left",fontsize=9,handlelength=1.4,borderaxespad=0.1,labelspacing=0.25)
    else:ax.text(0.44,0.88,"W303–M316",transform=ax.transAxes,ha="center",fontsize=9,color="#555555")
fig_b.text(0.07,0.023,"MD and Full use every original position;\nno terminal or core value is clipped.",fontsize=9,color="#555555",va="bottom")
for ext in ["pdf","png","svg"]:
    fig_b.savefig(ROOT/f"b_detail_draft_preview.{ext}",dpi=300,bbox_inches=None)
plt.close(fig_b)
print(json.dumps({"output_prefix":str(OUT),"font":font_path,"width_mm":183,"height_mm":210,
                  "text_below_9pt":0,"text_outside_canvas":sum(t["outside_canvas"] for t in text_audit),
                  "f_statistics":paired_statistics},ensure_ascii=False,indent=2))
