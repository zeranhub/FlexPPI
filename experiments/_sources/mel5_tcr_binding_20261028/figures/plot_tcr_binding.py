# Asset Confirmation Table
# Peptide profiles: LineTrend production asset is cumulative-publication data,
# incompatible with residue differences; inherit clean axes and external legend.
# Complete paired scatter: basic-plots contains no production implementation.
# Regional replicate plots: basic-plots contains no production implementation.
# All panels retain all predefined positions; only chemically undefined Hchi is masked.

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
import pandas as pd
import numpy as np
import json,hashlib
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
ROOT=Path(__file__).resolve().parent
source=ROOT.parent/'combined_binding_scores_wide.csv'
df=pd.read_csv(source)
methods=['ERGO-II McPAS','ERGO-II VDJdb','NetTCR2.2 known-overlap-excluded','NetTCR2.2 all20','ImRex downsample','ImRex full']
labels=['ERGO-II\nMcPAS','ERGO-II\nVDJdb','NetTCR-2.2\n4 models*','NetTCR-2.2\nall 20','ImRex\ndownsample','ImRex\nfull']
ylabels=['MEL5 / ELA (3HG1)','α24β17 / ELA (4JFF)','α24β17 / G4A (4JFD)','α24β17 / I7A (4JFE)','α24β17 / E1A','α24β17 / I5A','α24β17 / G6A','α24β17 / L8A','α24β17 / T9A','α24β17 / ILAKFLHWL','α24β17 / SLLMWITQC','α24β17 / YLEPGPVTA']
m=df[methods].to_numpy(float)
fig=plt.figure(figsize=(183/25.4,132/25.4))
gs=fig.add_gridspec(1,3,width_ratios=[.52,6,.16],left=.27,right=.96,top=.83,bottom=.09,wspace=.09)
ax_truth=fig.add_subplot(gs[0,0]);ax=fig.add_subplot(gs[0,1]);cax=fig.add_subplot(gs[0,2])
cmap=LinearSegmentedColormap.from_list('scientific_blue',SEQUENTIAL)
im=ax.imshow(m,cmap=cmap,vmin=0,vmax=1,aspect='auto',interpolation='nearest')
ax.set_xticks(np.arange(6),labels,fontsize=6.5)
ax.xaxis.tick_top();ax.tick_params(axis='x',length=0,pad=5)
ax.set_yticks([])
ax.set_xticks(np.arange(-.5,6,1),minor=True);ax.set_yticks(np.arange(-.5,12,1),minor=True)
ax.grid(which='minor',color='white',linewidth=.7);ax.tick_params(which='minor',length=0)
for spine in ax.spines.values():spine.set_visible(False)
for i in range(12):
 for j in range(6):
  value=m[i,j]
  label=f'{value:.3f}' if value>=.001 else f'{value:.2e}'
  ax.text(j,i,label,ha='center',va='center',fontsize=6.4,color='white' if value>.58 else BLACK)
for boundary in [3.5,8.5]:ax.axhline(boundary,color=BLACK,lw=.75)
ax_truth.set_xlim(0,1);ax_truth.set_ylim(11.5,-.5)
ax_truth.set_yticks(np.arange(12),ylabels,fontsize=6.7)
ax_truth.tick_params(axis='y',length=0,pad=5)
ax_truth.set_xticks([.5],['SPR\nbinding'],fontsize=6.5);ax_truth.xaxis.tick_top();ax_truth.tick_params(axis='x',length=0,pad=5)
for i,row in df.iterrows():ax_truth.text(.5,i,'+' if row.experimental_binding==1 else 'NB',ha='center',va='center',fontsize=7,fontweight='bold')
for spine in ax_truth.spines.values():spine.set_visible(False)
cb=fig.colorbar(im,cax=cax);cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=6.5,length=2,width=.5);cb.set_label('Binding score',fontsize=7,labelpad=4)
fig.suptitle('TCR–peptide binding predictions across the complete published panel',x=.53,y=.965,fontsize=8)
fig.text(.27,.012,'* Four official CV models exclude confirmed CDR3-pair overlap in partition 0.\nScores are not $K_D$; + = detected SPR binding, NB = reported no binding.',fontsize=6.5,ha='left',va='bottom')
for ext in ['pdf','svg','png']:fig.savefig(ROOT/f'figS8_TCR_binding_scores.{ext}',dpi=300,bbox_inches=None)
(ROOT/'figure_statistics.json').write_text(json.dumps(dict(figure='figS8_TCR_binding_scores',source=source.name,sha256=hashlib.sha256(source.read_bytes()).hexdigest(),configurations=methods,model_families=3,rows=12,points=72,experimental_positive=8,experimental_negative=4,score_scale='Each method raw sigmoid score; not a common calibrated probability or KD scale.',model_selection='All six configurations retained. NetTCR primary four exclude known exact CDR3-pair overlap; all20 shown as sensitivity.',inference_unit='One sequence pair scored per model; CV models are not biological repeats.',threshold='No cutoff used in this figure; 0.5 reference classification is saved only in CSV.',all_rows_retained=True,statistical_tests='None; descriptive single-family case panel, not population performance benchmark.'),ensure_ascii=False,indent=2),encoding='utf-8')
