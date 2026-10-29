# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Molecular views -> v49/render_structures.py -> existing verified real-coordinate raster
# (b) Complete CDR curves -> LineTrend/plot_trend.py -> param inherit
# (c left) Signed residue response -> BarComparison / signed comparison -> param inherit
# (c right) Paired context predictions -> PairedBoxScatter -> param inherit
# Companion previews inspected; all native assets have incompatible biological data structures.
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
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

from pathlib import Path
from datetime import datetime, timezone
import json, hashlib, shutil
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text

ROOT=Path(__file__).resolve().parent
EXP=ROOT.parent.parent
OLD=EXP/'mel5_joint_case_v49_20261019'
WT='#CF8C27'; ENG=CATEGORICAL[4]; TEAL='#297B82'; GRAY='#6B6B6B'
sources={
 'CDR':OLD/'figure4/source_plotted_CDR_RMSF.csv',
 'interface':OLD/'figure4/source_plotted_interface_Hchi.csv',
 'context':OLD/'case_sources/model_context/all_natural_pMHC_context_comparison.csv',
 'molecular_geometry':OLD/'figure4/molecular_view.json',
 'model_anchors':OLD/'case_sources/residue_story/trajectory_link/key_findings.json',
}
cdr=pd.read_csv(sources['CDR'],float_precision='round_trip')
interface=pd.read_csv(sources['interface'],float_precision='round_trip')
context=pd.read_csv(sources['context'],float_precision='round_trip')
key=json.loads(sources['model_anchors'].read_text(encoding='utf8'))
assert len(cdr)==45 and len(interface)==23
interface['delta_Hchi']=interface.hchi_full_b-interface.hchi_full_a
interface=interface.sort_values(['delta_Hchi','chain_a','position_a'],kind='stable').reset_index(drop=True)
assert interface.iloc[0].chain_a=='A' and int(interface.iloc[0].position_a)==72
assert int((interface.delta_Hchi<0).sum())==17
qcontext=context[(context.chain_a=='A') & (context.position_a==72)].iloc[0]
assert qcontext.same_aa and qcontext.hchi_valid
qdelta=key['model']['Q72_Hchi_delta_nats']
gdelta=key['model']['beta_Gly_RMSF_delta_A']

def image_crop(path,pad=12):
    im=Image.open(path).convert('RGB');a=np.asarray(im)
    yy,xx=np.where(np.any(a<238,axis=2))
    box=[max(0,int(xx.min())-pad),max(0,int(yy.min())-pad),min(im.width,int(xx.max())+pad),min(im.height,int(yy.max())+pad)]
    return im.crop(box),box

def heading(fig,letter,title,y):
    fig.text(.018,y,letter,fontsize=11,fontweight='bold')
    fig.text(.062,y,title,fontsize=8.7)

fig=plt.figure(figsize=(183/25.4,157/25.4),facecolor='white')
heading(fig,'a','Q72–βPhe53/βGly54: a shared interface contact network',.973)
overview,bb1=image_crop(OLD/'figure4/4jff_overview.png')
ax=fig.add_axes([.055,.72,.275,.225]);ax.imshow(overview);ax.axis('off')
ax.text(.02,.96,'TCR α',transform=ax.transAxes,color=GRAY,fontsize=8)
ax.text(.72,.96,'TCR β',transform=ax.transAxes,color=ENG,fontsize=8)
fig.text(.06,.713,'HLA-A*02:01',color=TEAL,fontsize=8)
fig.text(.242,.713,'ELA peptide',color=WT,fontsize=8)
fig.text(.065,.69,'MEL5  →  α24β17',fontsize=8)
fig.text(.065,.669,r'$K_D$: 18 μM → 0.6 nM (SPR [28])',fontsize=8)
local,bb2=image_crop(OLD/'figure4/4jff_Q72_F53_G54.png')
loc=fig.add_axes([.35,.738,.40,.208]);loc.imshow(local);loc.axis('off')
loc.annotate('HLA Q72',xy=(.45,.30),xytext=(.02,.38),xycoords='axes fraction',fontsize=8,color=TEAL,weight='bold',arrowprops={'arrowstyle':'-','lw':.6,'color':TEAL})
loc.text(.16,.97,'βPhe53',color=ENG,transform=loc.transAxes,fontsize=8,weight='bold')
loc.text(.79,.87,'βGly54',color=ENG,transform=loc.transAxes,fontsize=8,weight='bold')
loc.text(.44,.71,'N',color=BLACK,transform=loc.transAxes,fontsize=7)
loc.text(.84,.56,'O',color=BLACK,transform=loc.transAxes,fontsize=7)
loc.text(.48,.43,'OE1',color=TEAL,transform=loc.transAxes,fontsize=7)
loc.text(.69,.39,'NE2',color=TEAL,transform=loc.transAxes,fontsize=7)
loc.text(.38,.57,'2.94 Å',color=GRAY,transform=loc.transAxes,fontsize=7,ha='right')
loc.text(.85,.43,'3.12 Å',color=GRAY,transform=loc.transAxes,fontsize=7)
fig.text(.55,.713,'α24β17 · 4JFF',color=GRAY,fontsize=7.3,ha='center')
fig.text(.778,.894,'Model response',fontsize=8,weight='bold')
fig.text(.778,.859,'HLA Q72',color=TEAL,fontsize=8,weight='bold')
fig.text(.778,.831,f'ΔHχ  {qdelta:+.3f} nats',fontsize=8)
fig.text(.778,.79,'TCR βGly54',color=ENG,fontsize=8,weight='bold')
fig.text(.778,.762,f'ΔRMSF  {gdelta:+.3f} Å',fontsize=8)
fig.text(.778,.724,'Both residues unchanged',color=GRAY,fontsize=6.8)

heading(fig,'b','Backbone response across all six recognition loops',.625)
ax=fig.add_axes([.095,.398,.866,.175])
groups=['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
ticks=[];ticklabels=[]
for region in groups:
    q=cdr[cdr.region==region].sort_values('plot_x')
    x=q.plot_x.to_numpy()
    if region=='CDR2β':ax.axvspan(x[0]-.6,x[-1]+.6,color=ENG,alpha=.055,lw=0)
    for field,color,marker,fill in [('MEL5_model_RMSF_A',WT,'o','white'),('alpha24beta17_model_RMSF_A',ENG,'s',ENG)]:
        ax.plot(x,q[field],color=color,marker=marker,ms=2.7,mfc=fill,mec=color,lw=.95,zorder=3)
    ax.text(x.mean(),1.045,region,transform=ax.get_xaxis_transform(),ha='center',fontsize=8)
    ax.plot([x[0]-.2,x[-1]+.2],[1.015,1.015],transform=ax.get_xaxis_transform(),color='#B5B5B5',lw=.6,clip_on=False)
    choose=[0,len(q)-1] if region!='CDR3β' else [0,6,len(q)-1]
    for j in choose:ticks.append(x[j]);ticklabels.append(str(int(q.iloc[j].common_TCR_position)))
    if region=='CDR2β':
        g=q[q.common_TCR_position==54].iloc[0]
        ax.annotate('Gly54',xy=(g.plot_x,g.alpha24beta17_model_RMSF_A),xytext=(g.plot_x+2.4,.99),fontsize=7.5,color=ENG,arrowprops={'arrowstyle':'-','lw':.6,'color':ENG})
ax.set(xlim=(-1,cdr.plot_x.max()+1.5),ylim=(0,1.42),ylabel='RMSF (Å)',xticks=ticks,yticks=[0,.4,.8,1.2])
ax.set_xticklabels(ticklabels)
ax.tick_params(length=2.5,pad=3)
ax.set_xlabel('Common TCR residue numbering',labelpad=3,color=GRAY)
handles=[Line2D([],[],color=WT,marker='o',mfc='white',ms=3.5,lw=1,label='MEL5'),Line2D([],[],color=ENG,marker='s',mfc=ENG,ms=3.3,lw=1,label='α24β17')]
fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(.97,.654),ncol=2,fontsize=8,columnspacing=1.4,handlelength=1.5,handletextpad=.4,borderaxespad=0)

heading(fig,'c','Side-chain response: interface localization and TCR dependence',.306)
left=fig.add_axes([.095,.101,.545,.147])
x=np.arange(23)
colors=[TEAL if (r.chain_a=='A' and int(r.position_a)==72) else '#AEB7BD' for r in interface.itertuples()]
bars=left.bar(x,interface.delta_Hchi,width=.68,color=colors,edgecolor='none',zorder=3)
left.axhline(0,color=GRAY,lw=.6)
left.set(ylim=(-.50,.10),xlim=(-.8,22.8),yticks=[-.4,-.2,0],ylabel='ΔHχ (nats)')
labels=[('' if r.chain_a=='A' else 'p')+f'{r.aa_a}{int(r.position_a)}' for r in interface.itertuples()]
left.set_xticks(x,labels,rotation=60,ha='right',rotation_mode='anchor',fontsize=6.2)
for tick,r in zip(left.get_xticklabels(),interface.itertuples()):
    if r.chain_a=='A' and int(r.position_a)==72:tick.set_color(TEAL);tick.set_fontweight('bold')
left.text(.0,1.10,'All 23 valid interface sites · ranked by change',transform=left.transAxes,fontsize=7.2)
left.text(.02,.03,'Q72',transform=left.transAxes,color=TEAL,fontsize=7,weight='bold')
left.tick_params(length=2.2,pad=3)

right=fig.add_axes([.825,.101,.14,.147])
for row,fields in enumerate([('hchi_full_a','hchi_full_b'),('hchi_pMHC_only_a','hchi_pMHC_only_b')]):
    first,second=[float(qcontext[f]) for f in fields]
    right.plot([first,second],[1-row,1-row],color='#929292',lw=.85,zorder=1)
    right.plot(first,1-row,'o',mfc='white',mec=WT,ms=4,zorder=3)
    right.plot(second,1-row,'s',mfc=ENG,mec=ENG,ms=3.7,zorder=3)
    right.text((first+second)/2,1-row+.24,f'{first:.3f} → {second:.3f}',ha='center',fontsize=6.4)
right.set(xlim=(1.18,2.47),ylim=(-.52,1.54),xticks=[1.5,2.0,2.5],yticks=[1,0],yticklabels=['Full\ncomplex','TCR\nremoved'])
right.tick_params(axis='y',length=0,pad=5,labelsize=7)
right.tick_params(axis='x',length=2.2,pad=3,labelsize=7)
right.spines['left'].set_visible(False)
right.set_xlabel('Q72 Hχ (nats)',fontsize=7.5,labelpad=3)
right.text(-.72,1.10,'Q72 · fixed pMHC coordinates',transform=right.transAxes,fontsize=7.2)
fig.text(.095,.020,'Δ = α24β17 − MEL5; p, peptide. TCR removal changes model input only; pMHC coordinates are retained.',fontsize=6.5,color=GRAY)

fig.canvas.draw();renderer=fig.canvas.get_renderer();texts=[]
for item in fig.findobj(match=Text):
    if not item.get_visible() or not item.get_text():continue
    bb=item.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    texts.append(dict(text=item.get_text(),font_pt=item.get_fontsize(),bounds=list(bb.extents),inside_canvas=bool(bb.x0>=-.001 and bb.y0>=-.001 and bb.x1<=1.001 and bb.y1<=1.001)))
assert all(t['inside_canvas'] for t in texts),[t for t in texts if not t['inside_canvas']]
assert min(t['font_pt'] for t in texts)>=6.2
assert not any('\u4e00'<=char<='\u9fff' for t in texts for char in t['text'])
for ext in ['png','pdf','svg']:
    with mpl.rc_context({'savefig.bbox':None}):fig.savefig(ROOT/f'fig4_mel5_focus.{ext}',dpi=450,bbox_inches=None)
plt.close(fig)
cdr.to_csv(ROOT/'source_CDR_RMSF_all45.csv',index=False,encoding='utf-8-sig')
interface.to_csv(ROOT/'source_interface_Hchi_all23_ranked.csv',index=False,encoding='utf-8-sig')
qcols=['hchi_full_a','hchi_full_b','hchi_pMHC_only_a','hchi_pMHC_only_b','hchi_full_delta','hchi_pMHC_only_delta']
pd.DataFrame([{k:float(qcontext[k]) for k in qcols}]).to_csv(ROOT/'source_Q72_model_context.csv',index=False,encoding='utf-8-sig')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
audit=dict(status='generated_pending_visual_QA',dimensions_mm=[183,157],raster_dpi=450,minimum_font_pt=min(t['font_pt'] for t in texts),model_checkpoint_sha256='ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',sources={k:dict(path=str(p),sha256=sha(p)) for k,p in sources.items()},model_anchors=dict(Q72_Hchi_delta_nats=qdelta,betaG54_RMSF_delta_A=gdelta,Q72_full_context_delta=float(qcontext.hchi_full_delta),Q72_TCR_removed_delta=float(qcontext.hchi_pMHC_only_delta)),all_CDR_positions=45,all_valid_interface_positions=23,Hchi_decrease_n=17,Q72_is_largest_predicted_Hchi_decrease=True,old_panel_d_removed=True,existing_MD_profiles_plotted=False,no_new_inference_or_MD=True,selection_rule='All45CDR and all23same-AA chemically valid interface positions retained; no prediction-quality filtering.',statistics='One deterministic model output per state/context; no error bars or significance tests; source mean three-replica MD occupancy used in prose only.',interpretation='Local joint dynamic response and input-context dependence; no prediction of KD, causal Hchi-to-RMSF direction or real apo state.',visible_text=texts)
(ROOT/'figure_data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'status':audit['status'],'figure':str(ROOT/'fig4_mel5_focus.png'),'interface_positions':23,'CDR_positions':45,'old_d_removed':True},ensure_ascii=False))
