# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Actual molecular overview + one local view -> archived PyMOL render -> native raster placement
# (b,c) Residue profiles -> LineTrend/plot_trend.py -> param inherit
# (d upper/lower) Complete aligned delta maps -> heatmap/plot_composition.py -> param inherit
# Assets and companion PNGs inspected. Monthly counts, grouped fractions and clinical-stage
# counts have different data/axis semantics; preserve clean axes and grouping, use existing
# scientific source mappings. All user-requested values retained; native vectors for all charts.
# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
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
import hashlib, json, shutil, sys
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.text import Text

FIG=Path(__file__).resolve().parent
ROOT=FIG.parent
EXP=ROOT.parent
sys.path.insert(0,str(FIG))
from publication_style import BLUE, WARM, MD, MID, LIGHT, TEAL, SIGNED, signed_cmap, signed_norm
W,H=183,152
WT=WARM; ENG=BLUE; GRAY=MD
G4A=MD; I7A=MID
NA=LIGHT
STATES=['3hg1','4jff','4jfd','4jfe']
NAMES=['MEL5 · ELA','α24β17 · ELA','α24β17 · G4A','α24β17 · I7A']
COLORS=[WT,ENG,G4A,I7A]
MARKERS=['o','s','D','^']
STYLES=['-','-','--',':']
PAIRS=[('3hg1','4jff'),('4jff','4jfd'),('4jff','4jfe')]
PAIR_NAMES=['TCR maturation','Peptide G4A','Peptide I7A']
GROUPS=['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
paths={
    'CDR_absolute':FIG/'source_CDR_absolute.csv',
    'interface_absolute':FIG/'source_interface_absolute.csv',
    'peptide_four_states':FIG/'source_peptide_four_states.csv',
    'CDR_three_contrasts':FIG/'source_CDR_three_contrasts.csv',
    'interface_three_contrasts':FIG/'source_interface_three_contrasts.csv',
    'overview':FIG/'structure/4jff_overview.png',
    'local_structure':FIG/'structure/4jff_Q72_F53_G54.png',
}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
cdr=pd.read_csv(paths['CDR_absolute'],float_precision='round_trip')
iface=pd.read_csv(paths['interface_absolute'],float_precision='round_trip')
pep=pd.read_csv(paths['peptide_four_states'],float_precision='round_trip')
cdelta=pd.read_csv(paths['CDR_three_contrasts'],float_precision='round_trip')
idelta=pd.read_csv(paths['interface_three_contrasts'],float_precision='round_trip')
assert [len(t) for t in [cdr,iface,pep,cdelta,idelta]]==[45,23,40,135,90]

fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
def axes(x,y,w,h): return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def text(x,y,s,**kw): return fig.text(x/W,1-y/H,s,**kw)
def heading(letter,title,y,x=2,title_x=9):
    text(x,y,letter,fontsize=11,fontweight='bold',va='baseline')
    text(title_x,y,title,fontsize=10,va='baseline',fontweight='medium')
def crop(path,pad=12):
    # Freeze the exact original crop so colour changes cannot shift atom labels.
    frozen=json.loads((FIG/'original_figure_data_audit_v59.json').read_text(encoding='utf-8'))['crop_boxes']
    box=frozen['overview' if Path(path).name=='4jff_overview.png' else 'local']
    im=Image.open(path).convert('RGB')
    return im.crop(box),box

# (a) User-selected actual overview and one local molecular view; labels are vectors.
heading('a','Q72–β53/β54 interface contacts',5)
ov,ovbox=crop(paths['overview'])
oa=axes(6,10,43,31.5);oa.imshow(ov);oa.axis('off')
oa.text(.02,.96,'TCR α',transform=oa.transAxes,color=MID,fontsize=8)
oa.text(.72,.96,'TCR β',transform=oa.transAxes,color=ENG,fontsize=8)
text(7,43,'HLA-A*02:01',fontsize=8,color='#92999F')
text(34,43,'ELA peptide',fontsize=8,color=WT)
text(8,48.5,r'MEL5 → α24β17    $K_D$: 18 μM → 0.6 nM',fontsize=8.2)
local,localbox=crop(paths['local_structure'])
la=axes(51,10,62,34);la.imshow(local);la.axis('off')
la.annotate('HLA Q72',xy=(.45,.30),xytext=(.02,.38),xycoords='axes fraction',fontsize=8.2,color=WARM,weight='bold',arrowprops=dict(arrowstyle='-',lw=.6,color=WARM))
for s,x,y,col,fs,weight in [('βPhe53',.16,.97,ENG,8.3,'bold'),('βGly54',.79,.87,ENG,8.3,'bold'),('N',.44,.71,BLACK,8,'normal'),('O',.84,.56,BLACK,8,'normal'),('OE1',.48,.43,WARM,8,'normal'),('NE2',.69,.39,WARM,8,'normal'),('2.94 Å',.27,.57,GRAY,8,'normal'),('3.12 Å',.85,.43,GRAY,8,'normal')]:
    la.text(x,y,s,color=col,transform=la.transAxes,fontsize=fs,fontweight=weight)


# (b) All ten peptide positions for every existing observed bound structure.
heading('b','Peptide dynamics',5,x=119,title_x=126)
b=axes(130,11,48.5,23)
handles=[]
for state,name,col,marker,style in zip(STATES,NAMES,COLORS,MARKERS,STYLES):
    q=pep[pep.pdb.eq(state)].sort_values('id'); assert len(q)==10
    line,=b.plot(q.id,q.rmsf_full,color=col,ls=style,marker=marker,ms=2.2,lw=.9,
                mfc='white' if state=='3hg1' else col,mec=col,mew=.6,label=name)
    handles.append(line)
for pos in [4,10]:
    b.axvline(pos,color='#AAB2B8',lw=.5,ls='--',zorder=0)
    b.text(pos,1.012,f'P{pos}',ha='center',va='bottom',fontsize=8,color=GRAY)
b.set(xlim=(.7,10.3),ylim=(.57,1.02),xticks=range(1,11),yticks=[.6,.7,.8,.9,1],ylabel='Pred. RMSF (Å)',xlabel='Peptide position')
b.tick_params(labelsize=8,length=2,pad=2)
b.xaxis.label.set_size(9); b.yaxis.label.set_size(9)
b.xaxis.labelpad=2;b.yaxis.labelpad=3
state_legend=fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(119/W,1-44/H),ncol=2,
           fontsize=8,handlelength=1.6,handletextpad=.45,columnspacing=.7,labelspacing=.6,borderaxespad=0)

# (c) Complete recognition-loop backbone response. Preserve all45 values.
heading('c','Backbone response across all six recognition loops',54)
c=axes(17,62,161.5,19)
ticks=[];ticklabels=[]
for region in GROUPS:
    q=cdr[cdr.region==region].sort_values('plot_x');x=q.plot_x.to_numpy()
    if region=='CDR2β':c.axvspan(x[0]-.6,x[-1]+.6,color=ENG,alpha=.055,lw=0)
    for field,col,marker,fill in [('MEL5_model_RMSF_A',WT,'o','white'),('alpha24beta17_model_RMSF_A',ENG,'s',ENG)]:
        c.plot(x,q[field],color=col,marker=marker,ms=2.4,mfc=fill,mec=col,lw=.9,mew=.6,zorder=3)
    c.text(x.mean(),1.075,region,transform=c.get_xaxis_transform(),ha='center',fontsize=8)
    c.plot([x[0]-.2,x[-1]+.2],[1.025,1.025],transform=c.get_xaxis_transform(),color='#B5B5B5',lw=.55,clip_on=False)
    choose=[0,len(q)-1] if region!='CDR3β' else [0,6,len(q)-1]
    for j in choose:ticks.append(x[j]);ticklabels.append(str(int(q.iloc[j].common_TCR_position)))
    if region=='CDR2β':
        g=q[q.common_TCR_position==54].iloc[0]
        c.annotate('Gly54',xy=(g.plot_x,g.alpha24beta17_model_RMSF_A),xytext=(g.plot_x+2.3,1.0),fontsize=8,color=ENG,arrowprops=dict(arrowstyle='-',lw=.6,color=ENG))
c.set(xlim=(-1,cdr.plot_x.max()+1.5),ylim=(0,1.42),ylabel='Pred. RMSF (Å)',xticks=ticks,yticks=[0,.4,.8,1.2])
c.set_xticklabels(ticklabels);c.tick_params(length=2,pad=2,labelsize=8)
c.set_xlabel('Common TCR residue numbering',labelpad=2,fontsize=8,color=GRAY)


# (d) Two stacked metrics in one joint-response panel. Existing three contrasts, preserving all alignment and eligibility flags.
rid=lambda x:str(int(x))
ref=cdelta[cdelta.pdb_before.eq('3hg1') & cdelta.pdb_after.eq('4jff')].reset_index(drop=True)
ref_keys=[(r.chain_before,rid(r.id_before)) for r in ref.itertuples()]
engineered_to_ref={(r.chain_after,rid(r.id_after)):key for key,r in zip(ref_keys,ref.itertuples())}
cmat=[];isrc=[];csrc=[]
for k,(before,after) in enumerate(PAIRS):
    q=cdelta[cdelta.pdb_before.eq(before) & cdelta.pdb_after.eq(after)];vals={}
    for r in q.itertuples():
        key=(r.chain_before,rid(r.id_before))
        if before!='3hg1':key=engineered_to_ref[key]
        vals[key]=r.rmsf_change_after_minus_before
        csrc.append(dict(contrast=PAIR_NAMES[k],reference_chain=key[0],reference_position=key[1],before_pdb=before,after_pdb=after,delta_RMSF_A=r.rmsf_change_after_minus_before))
    assert set(vals)==set(ref_keys)
    cmat.append([vals[key] for key in ref_keys])
ic=idelta[idelta.pdb_before.eq('3hg1') & idelta.pdb_after.eq('4jff')].reset_index(drop=True)
ikeys=[(r.chain_before,rid(r.id_before)) for r in ic.itertuples()]
imat=[]
for k,(before,after) in enumerate(PAIRS):
    q=idelta[idelta.pdb_before.eq(before) & idelta.pdb_after.eq(after)];vals={}
    for r in q.itertuples():
        key=(r.chain_before,rid(r.id_before));eligible=bool(r.hchi_valid and r.same_aa)
        val=float(r.hchi_change_after_minus_before) if eligible else np.nan
        vals[key]=val
        isrc.append(dict(contrast=PAIR_NAMES[k],chain=key[0],position=key[1],before_pdb=before,after_pdb=after,eligible=eligible,delta_Hchi_nats=val))
    assert set(vals)==set(ikeys)
    imat.append([vals[key] for key in ikeys])
cmat=np.asarray(cmat);imat=np.asarray(imat)
assert cmat.shape==(3,45) and imat.shape==(3,30)
assert np.isfinite(imat).sum(axis=1).tolist()==[23,23,22]
cmap=signed_cmap().copy();cmap.set_bad(NA)
norm=signed_norm(.45)
assert np.abs(cmat).max()<=.45 and np.nanmax(np.abs(imat))<=.45
heading('d','Joint backbone and side-chain responses',94)
e=axes(26,102,145,11)
im=e.imshow(cmat,cmap=cmap,norm=norm,aspect='auto',interpolation='nearest')
e.set_yticks(range(3),PAIR_NAMES,fontsize=8)
e.set_xticks(range(45),[key[1] for key in ref_keys],rotation=90,fontsize=8)
e.tick_params(length=0,pad=2)
starts=[0,6,12,21,27,32];ends=[6,12,21,27,32,45]
for s,t,name in zip(starts,ends,GROUPS):
    e.text(((s+t)/2)/45,1.13,name,transform=e.transAxes,ha='center',fontsize=8)
    if s:e.axvline(s-.5,color='white',lw=1)
for sp in e.spines.values():sp.set_visible(False)
cbax=axes(172,102,2,11)
cb=fig.colorbar(im,cax=cbax,ticks=[-.45,0,.45])
cb.ax.tick_params(labelsize=8,length=1.5,pad=1)
cb.ax.set_title('ΔRMSF\n(Å)',fontsize=8,pad=4)
cb.outline.set_linewidth(.4)
f=axes(26,126,145,11)
imf=f.imshow(imat,cmap=cmap,norm=norm,aspect='auto',interpolation='nearest')
f.set_yticks(range(3),PAIR_NAMES,fontsize=8)
f.set_xticks(range(30),[f'{r.aa_before}{int(r.id_before)}' for r in ic.itertuples()],rotation=90,fontsize=8)
f.tick_params(length=0,pad=2)
f.axvline(21.5,color='white',lw=1.5)
f.text(11/30,1.13,'HLA-A*02:01',transform=f.transAxes,ha='center',fontsize=8)
f.text(26/30,1.13,'Peptide',transform=f.transAxes,ha='center',fontsize=8)
for sp in f.spines.values():sp.set_visible(False)
q72=ikeys.index(('A','72'))
f.add_patch(Rectangle((q72-.5,-.5),1,3,fill=False,edgecolor=WARM,lw=.8))
f.get_xticklabels()[q72].set_color(WARM);f.get_xticklabels()[q72].set_fontweight('bold')
cbax2=axes(172,126,2,11)
cb2=fig.colorbar(imf,cax=cbax2,ticks=[-.45,0,.45])
cb2.ax.tick_params(labelsize=8,length=1.5,pad=1)
cb2.ax.set_title('ΔHχ\n(nats)',fontsize=8,pad=4)
cb2.outline.set_linewidth(.4)

fig.canvas.draw();renderer=fig.canvas.get_renderer();texts=[]
legend_bb=state_legend.get_window_extent(renderer)
xlabel_bb=b.xaxis.label.get_window_extent(renderer)
assert xlabel_bb.y0 > legend_bb.y1, ('Peptide xlabel overlaps legend',xlabel_bb,legend_bb)

for item in fig.findobj(match=Text):
    if not item.get_visible() or not item.get_text():continue
    bb=item.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    texts.append(dict(text=item.get_text(),font_pt=item.get_fontsize(),bounds=list(bb.extents),inside_canvas=bool(bb.x0>=0 and bb.y0>=0 and bb.x1<=1 and bb.y1<=1)))
assert all(t['inside_canvas'] for t in texts),[t for t in texts if not t['inside_canvas']]
assert min(t['font_pt'] for t in texts)>=8
assert not any('\u4e00'<=ch<='\u9fff' for t in texts for ch in t['text'])
for ext in ['png','pdf','svg']:
    with mpl.rc_context({'savefig.bbox':None}):fig.savefig(FIG/f'fig4_compact_palette.{ext}',dpi=600,bbox_inches=None)
plt.close(fig)
for name,p in paths.items():
    if p.suffix=='.csv' and p.resolve() != (FIG/f'source_{name}.csv').resolve():shutil.copy2(p,FIG/f'source_{name}.csv')
pd.DataFrame(csrc).to_csv(FIG/'source_CDR_heatmap_all135.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(isrc).to_csv(FIG/'source_interface_heatmap_all90.csv',index=False,encoding='utf-8-sig')
np.savez(FIG/'source_heatmap_arrays.npz',CDR_delta_RMSF=cmat,interface_delta_Hchi=imat)
audit=dict(status='generated_pending_visual_QA',dimensions_mm=[W,H],raster_dpi=600,minimum_font_pt=min(t['font_pt'] for t in texts),visible_text=texts,
    source_files={k:dict(path=str(p),sha256=sha(p)) for k,p in paths.items()},
    panels=4,removed_panel="v54 d: ranked Hchi bars",merged_panels="v54 e and f -> v55 d upper and lower",all_requested_remaining_data_retained=True,counts=dict(peptide_values=40,absolute_CDR_values=90,ranked_Hchi=0,CDR_heatmap=135,pMHC_heatmap=90,Hchi_valid_per_contrast=[23,23,22]),
    contrasts=[dict(name=n,before=a,after=b) for n,(a,b) in zip(PAIR_NAMES,PAIRS)],
    colors=dict(MEL5=WT,mature=ENG,G4A=G4A,I7A=I7A,HLA_Q72=WARM,signed=SIGNED,ineligible=NA),
    crop_boxes=dict(overview=ovbox,local=localbox),heatmap_range=[-.45,.45],
    shared_style='Fig2 publication_style.py; Fig3 blue/orange/neutral system',
    heatmap_scale=dict(type='symmetric zero-centred signed differences',negative=BLUE,zero=SIGNED[1],positive=WARM,
        direction='after minus before; negative decrease, positive increase',
        CDR_RMSF=dict(unit='angstrom',actual_min=float(cmat.min()),actual_max=float(cmat.max()),clipped=0),
        pMHC_Hchi=dict(unit='nats',actual_min=float(np.nanmin(imat)),actual_max=float(np.nanmax(imat)),clipped=0),
        interpretation='The matching numerical limits do not equate angstroms with nats. Each strip retains its own unit.'),
    causal_arrows=False,new_inference=0,new_training=0,new_MD=0,
    figure_sha256=sha(FIG/'fig4_compact_palette.png'))
(FIG/'figure_data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','dimensions_mm','minimum_font_pt','counts']},ensure_ascii=False))
