# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Real molecular views -> previous render_structures.py -> param inherit
# (b) Six CDR absolute profiles -> LineTrend/plot_trend.py -> param inherit
# (c) Paired absolute interface values -> PairedBoxScatter/plot_PairedBoxScatter1.R -> param inherit
# (d) Replicate occupancy -> GroupedBarChart/plot_GroupedBarChartv1.py -> param inherit
# RULE: all panels use param inheritance because production data structures are incompatible.
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parent
EXP = Path(r'E:\Project\FlexPPI\experiments\_sources')
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

import json, hashlib
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text

mpl.rcParams.update({'font.size':9,'axes.labelsize':9,'axes.titlesize':9,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,'mathtext.fontset':'dejavusans'})
WT='#CF8C27'; ENG='#762A83'; TEAL='#297B82'; GRAY='#6B6B6B'
V=EXP/'tcr_expanded_screen_v41_20261017'/'mel5_validation'
T=EXP/'tcr_joint_case_search_20261019'/'residue_story'/'trajectory_link'
sources={'cdr':V/'plot_sources'/'mel5_all6_CDR_per_position_source.csv','interface':V/'plot_sources'/'mel5_natural_pMHC_interface_source.csv','contact':T/'per_replica_Q72_beta_contact_and_RMSF.csv','both_hbond':T/'Q72_two_end_backbone_Hbond_same_frame_per_replica.csv'}
cdr=pd.read_csv(sources['cdr'],float_precision='round_trip')
interface_all=pd.read_csv(sources['interface'],float_precision='round_trip')
interface=interface_all[interface_all.hchi_valid].copy()
contact=pd.read_csv(sources['contact'],float_precision='round_trip')
both=pd.read_csv(sources['both_hbond'],float_precision='round_trip')
assert len(cdr)==45 and cdr.rmsf_valid.all()
assert len(interface_all)==30 and len(interface)==23 and interface.same_aa.all()
assert interface[['chain_a','position_a']].duplicated().sum()==0
assert np.isfinite(cdr[['rmsf_full_a','rmsf_full_b']]).all().all()
assert np.isfinite(interface[['hchi_full_a','hchi_full_b']]).all().all()
assert (cdr[['rmsf_full_a','rmsf_full_b']]>0).all().all()
assert (interface[['hchi_full_a','hchi_full_b']]>=0).all().all()
groups=['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
fig=plt.figure(figsize=(183/25.4,185/25.4),facecolor='white')

def heading(letter,title,y):
    fig.text(.018,y,letter,fontsize=11,fontweight='bold',va='baseline')
    fig.text(.06,y,title,fontsize=9.5,va='baseline')

def image_crop(path,pad=15):
    im=Image.open(path).convert('RGB')
    a=np.asarray(im); yy,xx=np.where(np.any(a<238,axis=2))
    box=[max(0,int(xx.min())-pad),max(0,int(yy.min())-pad),min(im.width,int(xx.max())+pad),min(im.height,int(yy.max())+pad)]
    return im.crop(box),box

heading('a','MEL5 affinity maturation remodels a TCR–HLA contact network',.973)
overview,bb1=image_crop(ROOT/'4jff_overview.png')
ax=fig.add_axes([.045,.752,.31,.198]);ax.imshow(overview);ax.axis('off')
ax.text(.01,.97,'TCR α',transform=ax.transAxes,color=GRAY,fontsize=9)
ax.text(.75,.97,'TCR β',transform=ax.transAxes,color=ENG,fontsize=9)
fig.text(.057,.746,'HLA-A*02:01',color=TEAL,fontsize=9)
fig.text(.245,.746,'ELA peptide',color=WT,fontsize=9)
local,bb2=image_crop(ROOT/'4jff_Q72_F53_G54.png')
axloc=fig.add_axes([.408,.745,.55,.215]);axloc.imshow(local);axloc.axis('off')
axloc.annotate('HLA Q72',xy=(.46,.30),xytext=(.04,.38),xycoords='axes fraction',fontsize=9,color=TEAL,weight='bold',arrowprops={'arrowstyle':'-','lw':.6,'color':TEAL})
axloc.text(.17,.97,'βPhe53',color=ENG,transform=axloc.transAxes,fontsize=9,weight='bold')
axloc.text(.79,.87,'βGly54',color=ENG,transform=axloc.transAxes,fontsize=9,weight='bold')
axloc.text(.44,.71,'N',color='#333333',transform=axloc.transAxes,fontsize=9)
axloc.text(.84,.56,'O',color='#333333',transform=axloc.transAxes,fontsize=9)
axloc.text(.48,.43,'OE1',color=TEAL,transform=axloc.transAxes,fontsize=9)
axloc.text(.69,.39,'NE2',color=TEAL,transform=axloc.transAxes,fontsize=9)
axloc.text(.38,.57,'2.94 Å',color=GRAY,transform=axloc.transAxes,fontsize=9,ha='right')
axloc.text(.85,.43,'3.12 Å',color=GRAY,transform=axloc.transAxes,fontsize=9)
fig.text(.068,.724,'MEL5',color=WT,fontsize=9,weight='bold')
fig.text(.145,.724,'→',color=GRAY,fontsize=9)
fig.text(.178,.724,'α24β17',color=ENG,fontsize=9,weight='bold')
fig.text(.068,.705,r'$K_D$: 18 μM → 0.6 nM (SPR)',fontsize=9)
fig.text(.58,.719,'α24β17 complex · 4JFF',fontsize=9,color=GRAY,ha='center')

heading('b','Model: backbone RMSF across all six CDRs',.669)
ax=fig.add_axes([.095,.47,.866,.16])
xpos=0; cdr_records=[]; ticks=[]; ticklabels=[]
for k,region in enumerate(groups):
    d=cdr[cdr.region==region].sort_values('position_a')
    x=np.arange(len(d))+xpos
    if region=='CDR2β':ax.axvspan(x[0]-.6,x[-1]+.6,color=ENG,alpha=.055,lw=0)
    for suffix,col,marker,fill in [('a',WT,'o','white'),('b',ENG,'s',ENG)]:
        ax.plot(x,d[f'rmsf_full_{suffix}'],color=col,marker=marker,ms=3.1,mfc=fill,mec=col,lw=1.05,zorder=3)
    label=region.replace('α','α ').replace('β','β ')
    ax.text(x.mean(),1.04,label,transform=ax.get_xaxis_transform(),ha='center',fontsize=9)
    ax.plot([x[0]-.25,x[-1]+.25],[1.01,1.01],transform=ax.get_xaxis_transform(),color='#B5B5B5',lw=.6,clip_on=False)
    choose=[0,len(d)-1]
    if region=='CDR3β':choose=[0,6,len(d)-1]
    for j in choose:ticks.append(x[j]);ticklabels.append(str(int(d.iloc[j].position_a)))
    if region=='CDR2β':
        g=d[d.position_a==54].iloc[0]; gx=x[list(d.position_a).index(54)]
        ax.annotate('Gly54',xy=(gx,g.rmsf_full_a),xytext=(gx+4.6,1.19),fontsize=9,color=ENG,ha='center',arrowprops={'arrowstyle':'-','lw':.6,'color':ENG})
    for j,(_,r) in enumerate(d.iterrows()):cdr_records.append({'plot_x':int(x[j]),'common_TCR_position':int(r.position_a),'derivative_input_position':int(r.position_b),'region':region,'MEL5_aa':r.aa_a,'alpha24beta17_aa':r.aa_b,'MEL5_model_RMSF_A':r.rmsf_full_a,'alpha24beta17_model_RMSF_A':r.rmsf_full_b})
    xpos=x[-1]+3
ax.set(xlim=(-1,xpos-1.5),ylim=(0,1.45),ylabel='RMSF (Å)',xticks=ticks,yticks=[0,.4,.8,1.2])
ax.set_xticklabels(ticklabels);ax.tick_params(length=2.5,pad=3)
fig.text(.527,.425,'Common TCR residue numbering · segments are separate CDRs',fontsize=9,color=GRAY,ha='center')
handles=[Line2D([],[],color=WT,marker='o',mfc='white',ms=4,lw=1,label='MEL5'),Line2D([],[],color=ENG,marker='s',mfc=ENG,ms=3.7,lw=1,label='α24β17')]
fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(.963,.697),ncol=2,columnspacing=1.4,handlelength=1.5,handletextpad=.45,borderaxespad=0)

heading('c','Model: side-chain entropy at all 23 valid pMHC interface sites',.405)
ax=fig.add_axes([.095,.258,.866,.105])
interface=interface.sort_values(['chain_a','position_a']).reset_index(drop=True)
ix=np.arange(len(interface),dtype=float)
sep=int((interface.chain_a=='A').sum())
ix[sep:]+=1
for i,r in interface.iterrows():
    if r.chain_a=='A' and r.position_a==72:ax.axvspan(ix[i]-.48,ix[i]+.48,color=TEAL,alpha=.08,lw=0)
    ax.plot([ix[i]-.14,ix[i]+.14],[r.hchi_full_a,r.hchi_full_b],color='#9E9E9E',lw=.7,zorder=1)
    ax.plot(ix[i]-.14,r.hchi_full_a,'o',mfc='white',mec=WT,ms=3.6,zorder=3)
    ax.plot(ix[i]+.14,r.hchi_full_b,'s',mfc=ENG,mec=ENG,ms=3.2,zorder=3)
ax.set(xlim=(-.8,ix[-1]+.8),ylim=(0,3.4),yticks=[0,1,2,3],ylabel=r'$H_\chi$ (nats)',xticks=ix)
labels=[f'{r.aa_a}{int(r.position_a)}' for _,r in interface.iterrows()]
ax.set_xticklabels(labels,rotation=60,ha='right',rotation_mode='anchor')
ax.tick_params(length=2.5,pad=3)
for tick,(_,r) in zip(ax.get_xticklabels(),interface.iterrows()):
    if r.chain_a=='A' and r.position_a==72:tick.set_color(TEAL);tick.set_fontweight('bold')
ax.text(ix[:sep].mean(),1.10,'HLA-A*02:01',transform=ax.get_xaxis_transform(),ha='center',fontsize=9,color=GRAY)
ax.text(ix[sep:].mean(),1.10,'Peptide',transform=ax.get_xaxis_transform(),ha='center',fontsize=9,color=GRAY)
ax.axvline((ix[sep-1]+ix[sep])/2,color='#D6D6D6',lw=.6)

heading('d','Existing MD: frame occupancy in three independent replicas',.173)
ax=fig.add_axes([.335,.055,.626,.10])
rep_records=[]
for j,(name,df,colname) in enumerate([('Q72–Gly54 contact',contact[contact.beta_common_position==54],'min_heavy_contact_le4_A_fraction'),('Q72 two-ended bridge',both,'both_same_frame_Hbond_fraction')]):
    ypos=2*(1-j)
    for pdb,col,marker,fill,dy in [('3hg1',WT,'o','white',.33),('4jff',ENG,'s',ENG,-.33)]:
        d=df[df.pdb==pdb].sort_values('replica')
        values=d[colname].to_numpy()*100
        offsets=np.array([-.24,0,.24])
        ax.scatter(values,ypos+dy+offsets,c=fill,edgecolors=col,marker=marker,s=14,linewidths=.8,zorder=3)
        mean=float(values.mean());ax.plot([mean,mean],[ypos+dy-.12,ypos+dy+.12],color=col,lw=1.5,zorder=2)
        for _,r in d.iterrows():rep_records.append({'metric':name,'pdb':pdb,'replica':int(r.replica),'frames':int(r.frames),'frame_occupancy_percent':100*r[colname]})
        ax.text(103,ypos+dy,f'{mean:.1f}%',fontsize=9,va='center',color=col)
ax.set(xlim=(-2,115),ylim=(-.7,2.7),xticks=[0,25,50,75,100],yticks=[2,0],yticklabels=['Q72–Gly54 contact','Q72 two-ended bridge'],xlabel='Frame occupancy (%)')
ax.tick_params(axis='y',length=0,pad=8);ax.tick_params(axis='x',length=2.5,pad=3)
ax.spines['left'].set_visible(False)

fig.canvas.draw()
renderer=fig.canvas.get_renderer()
textaudit=[]
for item in fig.findobj(match=Text):
    if not item.get_visible() or not item.get_text():continue
    bb=item.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    textaudit.append({'text':item.get_text(),'font_pt':float(item.get_fontsize()),'bounds_figure':list(bb.extents),'inside_canvas':bool(bb.x0>=-.001 and bb.y0>=-.001 and bb.x1<=1.001 and bb.y1<=1.001)})
assert min(r['font_pt'] for r in textaudit)>=9
assert all(r['inside_canvas'] for r in textaudit),[r for r in textaudit if not r['inside_canvas']]
assert not any('\u4e00'<=c<='\u9fff' for r in textaudit for c in r['text'])
assert np.asarray(cdr[['rmsf_full_a','rmsf_full_b']]).max()<1.45
assert np.asarray(interface[['hchi_full_a','hchi_full_b']]).max()<3.4
for ext in ['png','svg','pdf']:
    with mpl.rc_context({'savefig.bbox':None}):fig.savefig(ROOT/f'fig4_mel5_joint.{ext}',dpi=450,bbox_inches=None)
pd.DataFrame(cdr_records).to_csv(ROOT/'source_plotted_CDR_RMSF.csv',index=False)
ip=interface[['role','chain_a','position_a','aa_a','hchi_full_a','hchi_full_b']].copy();ip.insert(0,'plot_x',ix)
ip.to_csv(ROOT/'source_plotted_interface_Hchi.csv',index=False)
pd.DataFrame(rep_records).to_csv(ROOT/'source_plotted_existing_MD_occupancy.csv',index=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
audit={'status':'passed','dimensions_mm':[183,185],'PNG_dpi':450,'minimum_font_pt':min(r['font_pt'] for r in textaudit),'font_family':'Arial','source_files':{k:{'path':str(v),'SHA256':sha(v)} for k,v in sources.items()},'all_CDR_RMSF_positions':45,'all_pMHC_interface_positions':30,'all_valid_Hchi_interface_positions':23,'invalid_Hchi_positions_omitted':7,'invalid_Hchi_reason':'No side-chain chi angle; predefined valid flag, not prediction quality. All 30 positions retained in original source.','complete_regions_retained':groups,'model_checkpoint_SHA256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5','model':'joint Full fold 0, trained on 1866 ATLAS + 564 DynaRepo; zero ITrace training systems','model_error_bars':None,'MD_replicate_definition':'Three pre-existing 200 ns trajectories per structure, all saved frames, no new MD. Unweighted replicate means.','MD_occupancy_criteria':{'contact':'Minimum Q72 side-chain heavy to beta Gly54 heavy distance <= 4 Angstrom','two_ended_bridge':'Both Q72 OE1--beta53 N-H and Q72 NE2-H--beta54 O satisfy D-A <=3.5 Angstrom and D-H-A >=150 degrees in same frame'},'statistical_tests':None,'MD_curves_on_model_profiles':False,'between_nonconsecutive_interface_sites_no_lines':True,'CDR_group_gaps_are_not_chain_breaks':True,'new_inference_training_MD':False,'semantic_colors':{'MEL5':WT,'alpha24beta17':ENG,'HLA_Q72_atoms':TEAL},'structure_crops':{'overview':bb1,'local':bb2},'all_plotted_data_inside_axis_limits':True,'RMSF_axis_limits_A':[0,1.45],'RMSF_data_range_A':[float(np.asarray(cdr[['rmsf_full_a','rmsf_full_b']]).min()),float(np.asarray(cdr[['rmsf_full_a','rmsf_full_b']]).max())],'Hchi_axis_limits_nats':[0,3.4],'Hchi_data_range_nats':[float(np.asarray(interface[['hchi_full_a','hchi_full_b']]).min()),float(np.asarray(interface[['hchi_full_a','hchi_full_b']]).max())],'visible_text':textaudit,'PNG_SHA256':sha(ROOT/'fig4_mel5_joint.png')}
(ROOT/'data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
(ROOT/'font_and_bounds_audit.json').write_text(json.dumps({'status':'passed','minimum_font_pt':9,'all_visible_text_inside_canvas':True,'text':textaudit},indent=2),encoding='utf8')
print(json.dumps({'status':'passed','CDR_positions':45,'interface_Hchi_positions':23,'file':str(ROOT/'fig4_mel5_joint.png'),'PNG_SHA256':audit['PNG_SHA256']}))
