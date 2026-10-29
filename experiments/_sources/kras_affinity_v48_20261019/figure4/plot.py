# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Actual molecular views -> source_archive/render_structures.py -> param inherit; molecular PNG bytes unchanged
# (b) All10 WT/G12D absolute RMSF profiles -> LineTrend/plot_trend.py -> param inherit
# RULE: no native-run panels; retain all ten observations; six MD replicas are archived, not plotted.
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parent.parent
sys.path.insert(0,str(EXP/'completed_results_20261012'/'plot_runtime'))
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

import json,hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

SOURCE=ROOT/'source_archive'
source_csv=SOURCE/'source_all10_peptide_model_and_MD_changes.csv'
p=pd.read_csv(source_csv,float_precision='round_trip').sort_values('position_a').copy()
assert len(p)==10 and p.rmsf_valid.all()
assert p.position_a.tolist()==list(range(1,11))
assert ''.join(p.aa_a)=='VVVGAGGVGK' and ''.join(p.aa_b)=='VVVGADGVGK'
assert p.loc[~p.same_aa,'position_a'].tolist()==[6]
assert np.isfinite(p[['rmsf_full_a','rmsf_full_b']].to_numpy()).all()
assert (p[['rmsf_full_a','rmsf_full_b']].to_numpy()>=0).all()
assert np.allclose(p.rmsf_full_b-p.rmsf_full_a,p.rmsf_full_delta,atol=1e-14)

def bounds(im):
    pix=np.asarray(im.convert('RGB')); yy,xx=np.where(np.any(pix<235,axis=2))
    return [int(xx.min()),int(yy.min()),int(xx.max()),int(yy.max())]
ims=[Image.open(ROOT/f'{pdb}_binding.png').convert('RGB') for pdb in ['7ow5','7ow6']]
bb=[bounds(im) for im in ims]
shared=(min(v[0] for v in bb)-20,min(v[1] for v in bb)-20,max(v[2] for v in bb)+20,max(v[3] for v in bb)+20)
overlay=Image.open(ROOT/'peptide_superposition.png').convert('RGB');pb=bounds(overlay)
pepcrop=(pb[0]-35,pb[1]-45,pb[2]+35,pb[3]+45)
fig=plt.figure(figsize=(183/25.4,130/25.4))
BLUE='#2166AC';MD='#555555';WT='#CF8C27';MUT='#762A83'
fig.text(.025,.96,'a',weight='bold',fontsize=10)
fig.text(.065,.96,'Single-residue antigen change in the same TCR–HLA complex',fontsize=9)
for i,(im,title,col) in enumerate(zip(ims,['WT · 7OW5','G12D · 7OW6'],[WT,MUT])):
    ax=fig.add_axes([.068+i*.265,.595,.255,.315]);ax.imshow(im.crop(shared));ax.axis('off')
    ax.set_title(title,color=col,pad=1,fontsize=9)
# The peptide is enlarged after the same MHC fit; no additional peptide fit is used.
axp=fig.add_axes([.64,.654,.33,.187]);axp.imshow(overlay.crop(pepcrop));axp.axis('off')
axp.set_title('Peptide superposition',fontsize=9,pad=12)
axp.annotate('G6 / D6',xy=(.57,.52),xytext=(.69,1.02),xycoords='axes fraction',textcoords='axes fraction',
             ha='center',fontsize=9,color=MUT,arrowprops={'arrowstyle':'-','lw':.65,'color':MUT})
fig.text(.806,.602,r'Peptide Cα RMSD = 0.61 Å',fontsize=9,ha='center',color='#444444')
fig.text(.806,.574,'MHC-aligned coordinates',fontsize=9,ha='center',color='#666666')
# WT and G12D sequence labels use precisely aligned letters; mutation is shown without shortening.
for y,name,seq,col,kd in [(.55,'WT','VVVGAGGVGK',WT,r'$K_D$ = 3 μM'),(.517,'G12D','VVVGADGVGK',MUT,r'$K_D$ = 743 pM')]:
    fig.text(.095,y,name,color=col,fontsize=9,ha='right')
    for j,aa in enumerate(seq):
        fig.text(.12+j*.025,y,aa,color=col if j==5 else '#444444',fontsize=9,weight='bold' if j==5 else 'normal',ha='center')
    fig.text(.4,y,kd,fontsize=9,color=col)
fig.text(.574,.533,'Binding measurements from ref. [28]',fontsize=9,color='#666666')

fig.text(.025,.455,'b',weight='bold',fontsize=10)
fig.text(.065,.455,'Predicted peptide RMSF in the WT and G12D complexes',fontsize=9)
ax=fig.add_axes([.095,.135,.87,.28])
x=p.position_a.to_numpy()
wt=p.rmsf_full_a.to_numpy();mut=p.rmsf_full_b.to_numpy()
ax.axvspan(5.6,6.4,color=MUT,alpha=.065,lw=0,zorder=0)
ax.plot(x,wt,marker='o',ms=4,mfc='white',mec=WT,color=WT,lw=1.25,zorder=3)
ax.plot(x,mut,marker='s',ms=3.6,mfc=MUT,mec=MUT,color=MUT,lw=1.25,zorder=4)
ax.plot([10,10.28],[wt[-1],wt[-1]],color=WT,lw=.8)
ax.plot([10,10.28],[mut[-1],mut[-1]],color=MUT,lw=.8)
ax.text(10.38,wt[-1],'WT',fontsize=9,va='center',color=WT)
ax.text(10.38,mut[-1],'G12D',fontsize=9,va='center',color=MUT)
ax.set(xlim=(.6,10.98),ylim=(0,1.4),ylabel='RMSF (Å)',xticks=x,yticks=[0,.4,.8,1.2],xlabel='Peptide position')
ax.xaxis.label.set_size(9);ax.yaxis.label.set_size(9)
ax.set_xticklabels([f'{aa}{i}' if i!=6 else 'G6→D6' for aa,i in zip(p.aa_a,x)])
for tick in ax.get_xticklabels():
    if tick.get_text()=='G6→D6':tick.set_color(MUT);tick.set_fontweight('bold')
ax.tick_params(length=3,pad=4,labelsize=9)
visible=[t.get_text() for t in fig.texts]+[t.get_text() for a in fig.axes for t in [*a.texts,a.xaxis.label,a.yaxis.label,a.title,*a.get_xticklabels(),*a.get_yticklabels()]]
assert not any('\u4e00'<=c<='\u9fff' for s in visible for c in s)
assert not any(s in '\n'.join(visible) for s in ['MD mean','SEM','Δ RMSF','MAE','All 10:','Unchanged 9:'])
with mpl.rc_context({'savefig.bbox':None}):
    for ext in ['pdf','svg','png']:
        fig.savefig(ROOT/f'fig4_kras_affinity.{ext}',dpi=450,bbox_inches=None)
plotted=p[['position_a','aa_a','aa_b','rmsf_full_a','rmsf_full_b']].rename(columns={'position_a':'peptide_position','aa_a':'WT_aa','aa_b':'G12D_aa','rmsf_full_a':'WT_model_RMSF_A','rmsf_full_b':'G12D_model_RMSF_A'})
plotted.to_csv(ROOT/'source_plotted_model_RMSF.csv',index=False)
original_audit=json.loads((SOURCE/'figure6_data_audit.json').read_text(encoding='utf8'))
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
image_assets={name:{'SHA256':sha(ROOT/name),'archive_SHA256':sha(SOURCE/name),'unchanged':sha(ROOT/name)==sha(SOURCE/name)} for name in ['7ow5_binding.png','7ow6_binding.png','peptide_superposition.png']}
assert all(a['unchanged'] for a in image_assets.values())
audit={'status':'passed','dimensions_mm':[183,130],'PNG_dpi':450,'all_peptide_positions_retained':10,'plotted_points':20,'unchanged_positions':9,'substituted_peptide_position':6,
 'WT_sequence':'VVVGAGGVGK','G12D_sequence':'VVVGADGVGK','model':'family_full_no_ITrace, joint two-output model; RMSF head only',
 'model_checkpoint_SHA256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',
 'no_new_training_inference_MD_or_calibration':True,'no_smoothing_or_residue_selection':True,
 'panel_a':{'molecular_PNG_assets':image_assets,'complex_shared_crop':shared,'peptide_crop':pepcrop,'crop_parameters_unchanged':list(shared)==original_audit['complex_shared_crop'] and list(pepcrop)==original_audit['peptide_crop'],'layout_preserved':True,'quantitative_values_preserved':True,'annotation_font_pt':9,'annotation_changes':'Font increased from 7–8 pt to 9 pt; overlay title padding increased to avoid collision.'},
 'panel_b':{'x':'all 10 peptide positions','y':'absolute, untransformed model RMSF in Å','WT_source_column':'rmsf_full_a','G12D_source_column':'rmsf_full_b','two_lines_only':True,'MD_curve_plotted':False,'MD_similarity_metrics_displayed':False,'error_bars':None,'statistical_tests':None,'color_mapping':{'WT':WT,'G12D':MUT},'font_pt':9,'y_limits':[0,1.4]},
 'source_file':str(source_csv),'source_SHA256':sha(source_csv),'source_archive_file_count':len(list(SOURCE.iterdir())),
 'binding_source':'ref. [28], doi:10.1038/s41467-022-32811-1; literature experiment, not model output','visible_text':visible,
 'model_values_all10':plotted.to_dict('records')}
assert audit['panel_a']['crop_parameters_unchanged']
(ROOT/'data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
manifest={path.name:sha(path) for path in SOURCE.iterdir() if path.is_file()}
(ROOT/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
print(json.dumps({'output':str(ROOT),'positions':len(p),'min_WT':float(wt.min()),'max_WT':float(wt.max()),'min_G12D':float(mut.min()),'max_G12D':float(mut.max())},indent=2))

