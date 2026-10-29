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
import hashlib
import json
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'contacts'/'all_peptide_neighbor_contact_contrasts.csv'
data=pd.read_csv(SOURCE)
assert len(data)==1760
fig,axes=plt.subplots(3,2,figsize=(183/25.4,155/25.4))
fig.subplots_adjust(left=.105,right=.985,bottom=.09,top=.95,wspace=.37,hspace=.65)
report=[]
for j,mutation in enumerate(['G4A','I7A']):
    ax=axes[0,j]; d=data[data.mutation==mutation]
    assert len(d)==880
    ax.text(-.16,1.08,'ab'[j],transform=ax.transAxes,fontsize=9,fontweight='bold')
    ax.set_title(f'{mutation}: all 880 peptide–neighbor pairs',loc='left',pad=8)
    ax.plot([0,100],[0,100],color=GREY,linestyle='--',linewidth=.65)
    ax.scatter(100*d.contact4_reference_mean,100*d.contact4_mutant_mean,s=9,color=GREY,alpha=.50,edgecolors='none')
    ax.set(xlim=(-3,103),ylim=(-3,103),xlabel='ELA contact occupancy (%)',ylabel=f'{mutation} contact occupancy (%)')
    report.append({'panel':'ab'[j],'mutation':mutation,'all_pairs':880,'no_zero_contact_pairs_removed':True})
selected=[('I7A',7,'E',96,'I7A: P7–Thrβ96'),('I7A',8,'E',97,'I7A: P8–Glyβ97'),
          ('I7A',8,'E',99,'I7A: P8–Glyβ99'),('G4A',8,'A',97,'G4A: P8–HLA Arg97')]
for k,(mutation,pos,chain,resid,title) in enumerate(selected):
    ax=axes[1+k//2,k%2]
    row=data[(data.mutation==mutation)&(data.peptide_position==pos)&
             (data.neighbor_chain==chain)&(data.paper_common_author_id==resid)]
    assert len(row)==1
    r=row.iloc[0]
    ax.text(-.16,1.08,'cdef'[k],transform=ax.transAxes,fontsize=9,fontweight='bold')
    ax.set_title(title,loc='left',pad=8)
    for q,(prefix,color,label) in enumerate([('contact4',CATEGORICAL[0],'≤4 Å'),('contact4p5',CATEGORICAL[3],'≤4.5 Å')]):
        for x,state in enumerate(['reference','mutant']):
            values=100*np.array([r[f'{prefix}_{state}_rep{i}'] for i in [1,2,3]])
            offset=(-.16 if q==0 else .16)
            ax.scatter(x+offset+np.array([-.045,0,.045]),values,s=19,color=color,
                       marker='o' if q==0 else '^',edgecolors='white',linewidth=.3,zorder=3)
            ax.plot([x+offset-.085,x+offset+.085],[values.mean()]*2,color=color,linewidth=1.15)
    ax.set(xticks=[0,1],xticklabels=['ELA',mutation],xlim=(-.43,1.43),ylim=(-3,103),ylabel='Contact occupancy (%)')
    report.append({'panel':'cdef'[k],'selected_exploratory_pair':title,
                   'values':{prefix:{state:[float(r[f'{prefix}_{state}_rep{i}']) for i in [1,2,3]]
                                    for state in ['reference','mutant']} for prefix in ['contact4','contact4p5']}})
handles=[plt.Line2D([],[],color=CATEGORICAL[0],marker='o',linestyle='None',label='≤4 Å'),
         plt.Line2D([],[],color=CATEGORICAL[3],marker='^',linestyle='None',label='≤4.5 Å')]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.54,.005),ncol=2,fontsize=7.5)
save_cns_figure(fig,ROOT/'figS7_peptide_contact_redistribution')
fig.savefig(ROOT/'figS7_peptide_contact_redistribution.svg',bbox_inches='tight',dpi=300)
plt.close(fig)
(ROOT/'contact_figure_statistics.json').write_text(json.dumps({'panels':report,
    'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'sample_unit':'3 independent trajectories per state; displayed points, no trajectory pairing.',
    'selection':'a/b all 880 pairs; c–f exploratory follow-up from complete screen.',
    'statistics':'Descriptive means and individual points; no p-values, water or free-energy claims.'},ensure_ascii=False,indent=2),encoding='utf-8')
print('Rendered complete contact overview and 4 exploratory replicated pairs')