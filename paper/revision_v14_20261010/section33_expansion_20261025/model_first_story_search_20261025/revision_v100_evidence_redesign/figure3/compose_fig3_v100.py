# Academic Figure Skill Asset Confirmation (verified against existing native production scripts)
# a -> original v96 PyMOL overview PNG/PSE -> native run, no image editing
# b/c -> v96 absolute residue profiles / LineTrend -> param inherit
# d/e/f -> v96 complete scatter and association plots -> param inherit
# All model and MD arrays are frozen; layout, typography, and annotations only.
# Style target: verified current manuscript Fig2; evidence-focused redesign with horizontal per-chi row.
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




# Reference-informed composition: MAVE-NN Fig4 / ENGINE Fig4 / RMSF-net Fig4.
# The published figures are design references, never source artwork for this figure.
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
ROOT=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent
P=OUT/'source_data';R=ROOT/'pymol'
W,H=183.,222.
BLUE='#2166AC';WARM='#B87936';TEAL='#3A8D8C';PURPLE='#8B72A5';GRAY='#BAC1C7';INK='#222222'
FS=9.0
mpl.rcParams.update({'savefig.bbox':None,'font.size':FS,'axes.labelsize':FS,'xtick.labelsize':FS,'ytick.labelsize':FS,'legend.fontsize':FS,'axes.unicode_minus':False,'axes.linewidth':.6})
cam=pd.read_csv(P/'cam_all25.csv').sort_values('CaMKI_native_position')
ub=pd.read_csv(P/'ubiquitin_72_position_profile_with_original_mask.csv').sort_values('residue_id')
gln=pd.read_csv(P/'ubiquitin_all6_Gln_Hchi_SASA_plot.csv').sort_values('residue_id')
corr=pd.read_csv(P/'retained_absolute_correlations_all498.csv')
stats=pd.read_csv(P/'retained_panel_statistics.csv')
manifest=json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
assert len(cam)==25 and ub.valid_hchi.sum()==66 and len(ub)==72 and len(gln)==6
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');plotted=[]
def txt(x,y,s,size=FS,**kw):return fig.text(x/W,1-y/H,s,fontsize=size,va=kw.pop('va','top'),color=kw.pop('color',INK),**kw)
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def panel(label,title,x,y):
    txt(x,y,label,11,fontweight='bold')
    if title:txt(x+6,y+.25,title,10,fontweight='normal')
def line(x1,y1,x2,y2,**kw):fig.add_artist(Line2D([x1/W,x2/W],[1-y1/H,1-y2/H],transform=fig.transFigure,**kw))
def values(panel,kind,ids,vals,**kw):
    for i,v in zip(ids,vals):
        if np.isfinite(v):plotted.append(dict(panel=panel,kind=kind,identity=str(i),value=float(v),**kw))
def molecule(key,x,y,w,h,labels):
    z=manifest['assets'][key];a=ax(x,y,w,h);a.imshow(Image.open(R/z['png']).convert('RGBA'));a.axis('off')
    # Display viewport only: keep original PNG/PSE and camera/coordinates intact.
    crop={'cam_overview':(200,110,1430,1315),'ub_overview':(90,205,1535,1195)}[key]
    a.set_xlim(crop[0],crop[2]);a.set_ylim(crop[3],crop[1]);a.set_aspect('equal')
    for name,location in labels.items():
        px=z['landmarks'][name]
        a.scatter([px[0]],[px[1]],s=11,facecolors='none',edgecolors=TEAL,linewidth=.85)
        a.annotate(name,xy=px,xytext=location,textcoords='axes fraction',color=INK,fontsize=FS,ha='center',va='center',arrowprops=dict(arrowstyle='-',color=TEAL,lw=.7),bbox=dict(facecolor='white',edgecolor='none',pad=.1,alpha=.92))
    return a

from publication_style import BLUE,WARM,TEAL,INK,MD,LIGHT
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],'font.size':9,'axes.labelsize':9,'axes.titlesize':10,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,'axes.linewidth':.6,'xtick.major.width':.6,'ytick.major.width':.6,'axes.unicode_minus':False,'pdf.fonttype':42,'svg.fonttype':'none'})
axis_records=[]
def style_axis(a,name):
    a.tick_params(length=2.5,pad=2,width=.6);a.xaxis.labelpad=2.5;a.yaxis.labelpad=2.5
    axis_records.append((name,a))
handles=[Line2D([],[],color=c,lw=lw,ls=ls,marker=mark,ms=2.2,mfc='white',mew=.6,label=name) for name,c,lw,ls,mark in [('MD',MD,.9,':','o'),('Full',BLUE,1.25,'-',None)]]

# a: retain the genuine molecular contexts together.
panel('a','',3,2)
txt(27,3,'CaM–CaMKI\n(1MXE)',ha='center',linespacing=1.12)
txt(70,3,'Ubiquitin–UBA\n(2QHO)',ha='center',linespacing=1.12)
molecule('cam_overview',7,14,39,37,{'W303':(.12,.68),'M316':(.91,.08)})
molecule('ub_overview',49,14,39,37,{'Q41':(.64,.96),'Q62':(.91,.10)})
for xx,label,col in [(13,'Analysed chain',PURPLE),(59,'Partner',GRAY)]:
    line(xx,55,xx+4,55,color=col,lw=2.2);txt(xx+5.5,53.4,label)

# b: resolve actual bound-state flexibility, with a separate static geometry track.
panel('b','CaMKI backbone and shielding',98,2)
fig.legend(handles=handles,loc='center',bbox_to_anchor=(145/W,1-10/H),ncol=2,handlelength=1.7,handletextpad=.4,columnspacing=1.3,borderpad=0,frameon=False)
x=cam.CaMKI_native_position.to_numpy()
a=ax(111,18,69,22);style_axis(a,'b_RMSF')
a.axvspan(302.5,316.5,color=TEAL,alpha=.055,lw=0,zorder=0)
for name,color,marker,ls,field,lw in [('MD',MD,'o',':','target_rmsf',.9),('Full',BLUE,None,'-','pred_rmsf',1.25)]:
    y=cam[field].to_numpy();a.plot(x,y,color=color,lw=lw,marker=marker,ms=2.1,ls=ls,mfc='white',mew=.55,zorder=4 if name=='Full' else 3)
    values('b','RMSF_A',x,y,series=name)
mae_full=float(np.mean(abs(cam.pred_rmsf-cam.target_rmsf)));mae_intra=float(np.mean(abs(cam.intra_rmsf-cam.target_rmsf)))
a.text(.47,.965,f'MAE {mae_full:.3f} Å',transform=a.transAxes,ha='center',va='top',fontsize=FS)
a.set(xlim=(293.5,318.5),ylim=(0,2.4),ylabel='RMSF (Å)');a.set_yticks([0,1,2]);a.set_xticks([294,300,306,312,318]);a.tick_params(labelbottom=False)
for pos,label,loc in [(303,'W303',(302.4,1.12)),(316,'M316',(315.3,1.18))]:
    q=cam[cam.CaMKI_native_position.eq(pos)].iloc[0]
    a.annotate(label,xy=(pos,q.pred_rmsf),xytext=loc,fontsize=FS,ha='center',color=INK,arrowprops=dict(arrowstyle='-',color=TEAL,lw=.6))
a=ax(111,44,69,9);style_axis(a,'b_static_shielding')
shield=cam.sidechain_shielding_fraction.to_numpy()*100
a.bar(x,shield,width=.78,color=TEAL,alpha=.60,edgecolor='none')
a.set(xlim=(293.5,318.5),ylim=(0,105),ylabel='CaM\nshielding (%)',xlabel='CaMKI residue')
a.set_yticks([0,100]);a.set_xticks([294,300,306,312,318]);a.tick_params(axis='y',length=2,pad=1)
values('b','static_sidechain_shielding_percent',x,shield,series='starting_structure')

# c/d: complete sequence accuracy, then same chemistry with different environments.
panel('c','Ubiquitin side chains',3,62)
panel('d','Same side chain, different environments',98,62)
a=ax(14,74,70,30);style_axis(a,'c')
for pos in [41,62]:a.axvspan(pos-.55,pos+.55,color=TEAL,alpha=.10,lw=0,zorder=0)
for name,color,marker,ls,field,lw in [('MD',MD,'o',':','MD_Hchi_nats',.8),('Full',BLUE,None,'-','Full_Hchi_nats',1.25)]:
    y=ub[field].where(ub.valid_hchi).to_numpy();a.plot(ub.residue_id,y,color=color,lw=lw,marker=marker,ms=1.9,ls=ls,mfc='white',mew=.45,zorder=3 if name=='Full' else 2)
    values('c','Hchi_nats',ub.residue_id,y,series=name)
mae_ub=float(np.mean(abs(ub.loc[ub.valid_hchi,'Full_Hchi_nats']-ub.loc[ub.valid_hchi,'MD_Hchi_nats'])))
a.text(.035,.98,f'MAE {mae_ub:.3f} nats',transform=a.transAxes,ha='left',va='top',fontsize=FS)
for pos,loc in [(41,(44.5,3.0)),(62,(65.,3.7))]:
    q=ub[ub.residue_id.eq(pos)].iloc[0]
    a.annotate('Q'+str(pos),xy=(pos,q.Full_Hchi_nats),xytext=loc,fontsize=FS,color=INK,arrowprops=dict(arrowstyle='-',color=TEAL,lw=.6))
a.set(xlim=(.5,72.5),ylim=(0,4.75),ylabel='Hχ (nats)',xlabel='Ubiquitin residue');a.set_xticks([1,20,40,60,72]);a.set_yticks([0,2,4])
ordered=gln.sort_values('MD_Hchi_nats');dx=np.arange(6)
a=ax(111,74,69,20);style_axis(a,'d_Hchi')
for i,row in enumerate(ordered.itertuples()):
    a.plot([i-.12,i+.12],[row.MD_Hchi_nats,row.Full_Hchi_nats],c=LIGHT,lw=.75,zorder=1)
    a.scatter(i-.12,row.MD_Hchi_nats,s=18,facecolor='white',edgecolor=MD,linewidth=.75,zorder=3)
    a.scatter(i+.12,row.Full_Hchi_nats,s=14,facecolor=BLUE,edgecolor=BLUE,linewidth=.4,zorder=4)
    for name,val in [('MD',row.MD_Hchi_nats),('Full',row.Full_Hchi_nats)]:values('d','Hchi_nats',[row.residue],[val],series=name,SASA_A2=float(row.sidechain_SASA_A2))
a.set(xlim=(-.6,5.6),ylim=(0,3.5),ylabel='Hχ (nats)');a.set_yticks([0,1,2,3]);a.set_xticks(dx);a.tick_params(labelbottom=False)
a=ax(111,98,69,9);style_axis(a,'d_SASA')
colors=[TEAL if q in ['Q41','Q62'] else '#BBC2C7' for q in ordered.residue]
a.bar(dx,ordered.sidechain_SASA_A2,width=.55,color=colors,edgecolor='none',alpha=.8)
a.set(xlim=(-.6,5.6),ylim=(0,110),ylabel='SASA\n(Å²)');a.set_yticks([0,100]);a.set_xticks(dx,ordered.residue);a.tick_params(axis='x',length=0,pad=2)
for t,q in zip(a.get_xticklabels(),ordered.residue):t.set_color(TEAL if q in ['Q41','Q62'] else INK)
values('d','sidechain_SASA_A2',ordered.residue,ordered.sidechain_SASA_A2,series='starting_structure')

# e: all four angle comparisons in one horizontal row, with identical limits.
panel('e','Per-χ occupancy association',3,115)
eh=[Line2D([],[],ls='',marker=m,ms=3.6,color=c,label=name) for name,m,c in [('Other','o',BLUE),('Q41','D',TEAL),('Q2','s',TEAL)]]
fig.legend(handles=eh,loc='center',bbox_to_anchor=(145/W,1-118/H),ncol=3,handletextpad=.25,columnspacing=.7,borderpad=0,frameon=False)
defined=pd.read_csv(P/'defined_angles_long.csv',float_precision='round_trip')
chi_stats=pd.read_csv(P/'per_chi_correlations.csv',float_precision='round_trip');chi_stats=chi_stats[chi_stats.predictor.eq('predicted_joint_Hchi')&chi_stats.metric.eq('fine_entropy')]
for chi,left in [(1,14),(2,57.5),(3,101),(4,144.5)]:
    a=ax(left,132,35.5,25);style_axis(a,'e_chi'+str(chi))
    z=defined[defined.chi.eq('chi'+str(chi))];st=chi_stats[chi_stats.chi.eq('chi'+str(chi))].iloc[0]
    assert len(z)=={1:66,2:52,3:23,4:10}[chi]
    txt(left,126.5,'χ'+str(chi)+'  n='+str(len(z))+'  ρ='+format(st.rho_AA_centered,'.2f'))
    for name,marker,col in [('Other','o',BLUE),('Q41','D',TEAL),('Q2','s',TEAL)]:
        sub=z[~z.author_residue_id.isin([2,41])] if name=='Other' else z[z.author_residue_id.eq(int(name[1:]))]
        a.scatter(sub.predicted_Hchi_nats_AA_centered,sub.fine_entropy_four_origins_mean_nats_AA_centered,s=8.5 if name=='Other' else 19,marker=marker,color=col,alpha=.44 if name=='Other' else 1,linewidths=.25,edgecolors='white' if name=='Other' else col,zorder=3 if name=='Other' else 4)
    for row in z.itertuples():
        plotted.extend([dict(panel='e',kind='predicted_joint_Hchi_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='Full',value=float(row.predicted_Hchi_nats_AA_centered)),dict(panel='e',kind='MD_single_chi_fine_entropy_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='MD',value=float(row.fine_entropy_four_origins_mean_nats_AA_centered))])
    a.axhline(0,color=LIGHT,lw=.4,ls=':');a.axvline(0,color=LIGHT,lw=.4,ls=':')
    a.set(xlim=(-.7,.5),ylim=(-1,.75));a.set_xticks([-.5,0,.5]);a.set_yticks([-.5,0,.5])
    if chi!=1:a.tick_params(labelleft=False)
txt(3.1,144.5,'MD entropy deviation (nats)',rotation=90,ha='center',va='center')
txt(97,163,'Predicted joint Hχ deviation (nats)',ha='center')

# f: compare two predictors on exactly the same record, not separated clouds.
panel('f','Hχ versus RMSF: matched records',3,173)
paired=[]
for metric,left,title in [('fine_entropy',28,'Fine-bin entropy'),('lowF_width',119,'Low-F* angular width')]:
    z=corr[corr.metric_key.eq(metric)&corr.predictor.isin(['full_hchi','full_rmsf'])]
    q=z.pivot(index='sample_id',columns='predictor',values='rho_type_centered').sort_index();assert len(q)==83 and not q.isna().any().any()
    wins=int((q.full_hchi>q.full_rmsf).sum())
    txt(left+24.5,181,title,ha='center')
    a=ax(left,186,49,25);style_axis(a,'f_'+metric)
    a.plot([-.06,.75],[-.06,.75],color='#999999',lw=.6,ls='--',zorder=0)
    a.scatter(q.full_rmsf,q.full_hchi,s=7.5,facecolor=BLUE,edgecolor='white',linewidth=.2,alpha=.48,zorder=2)
    a.text(.04,.96,f'Hχ higher: {wins}/83',transform=a.transAxes,ha='left',va='top',fontsize=FS)
    a.set(xlim=(-.06,.75),ylim=(-.06,.75));a.set_xticks([0,.4,.75]);a.set_yticks([0,.4,.75])
    a.set_ylabel('Hχ correlation ρ');a.set_xlabel('RMSF correlation ρ')
    for predictor in ['full_hchi','full_rmsf']:
        series=q[predictor];values('f','within_AA_Spearman_rho',series.index,series.values,series=predictor,metric=metric)
    for sid,row in q.iterrows():paired.append(dict(sample_id=sid,metric=metric,Hchi_rho=row.full_hchi,RMSF_rho=row.full_rmsf,Hchi_minus_RMSF=row.full_hchi-row.full_rmsf))

df=pd.DataFrame(plotted);df.to_csv(OUT/'main_plotted_values_v100.csv',index=False)
pd.DataFrame(paired).to_csv(OUT/'paired_record_comparison_v100.csv',index=False)
assert df.groupby('panel').size().to_dict()=={'b':75,'c':132,'d':18,'e':302,'f':332}
fig.canvas.draw();renderer=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in texts if t.get_window_extent(renderer).x0<-.5 or t.get_window_extent(renderer).y0<-.5 or t.get_window_extent(renderer).x1>fig.bbox.width+.5 or t.get_window_extent(renderer).y1>fig.bbox.height+.5]
assert not outside,outside
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'fig3_dual_case_entropy.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
qa=dict(size_mm=[W,H],min_font=min(t.get_fontsize() for t in texts),font_at_167p6mm=9*167.6/W,style='same as verified current Fig2; Arial9/10/11, fine axes, semantic palette',main_plot_counts=df.groupby('panel').size().to_dict(),new_geometry_values=31,original_inventory_853_unchanged=True,Intra_25_values_moved_to_supplement=True,all_masks_preserved=True,all6_Gln_preserved=True,per_chi_all151_points=True,e_four_horizontal=True,f_all83_matched_per_metric=True,f_higher_counts={metric:int(sum(z['Hchi_minus_RMSF']>0 for z in paired if z['metric']==metric)) for metric in ['fine_entropy','lowF_width']},no_joint_CI_constructed=True,no_binding_free_energy_or_apo_holo_claim=True,static_geometry_not_binding_dynamics=True,structure_assets_unchanged=True,annotation_metrics=dict(CaMKI_Full_RMSF_MAE_A=mae_full,CaMKI_Intra_RMSF_MAE_A=mae_intra,Ubiquitin_Full_Hchi_MAE_nats=mae_ub),outside_text=outside,axes_bounds_mm={name:[round(a.get_position().x0*W,3),round((1-a.get_position().y1)*H,3),round(a.get_position().width*W,3),round(a.get_position().height*H,3)] for name,a in axis_records},visual_review='pending',sha256={ext:hashlib.sha256((OUT/f'fig3_dual_case_entropy.{ext}').read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
plt.close(fig)

# Supplement retains the complete context ablation and old distribution/interval view.
supp=plt.figure(figsize=(183/25.4,75/25.4),facecolor='white')
def sa(x,y,w,h):return supp.add_axes([x/183,1-(y+h)/75,w/183,h/75])
supp.text(.02,.96,'a',fontweight='bold',fontsize=11,va='top');supp.text(.085,.96,'Same-bound-input context comparison',fontsize=10,va='top')
a=sa(14,18,70,42)
for name,c,ls,field in [('MD',MD,':','target_rmsf'),('Full',BLUE,'-','pred_rmsf'),('Intra',WARM,'--','intra_rmsf')]:a.plot(x,cam[field],color=c,ls=ls,lw=1.25 if name=='Full' else 1.0,label=name)
a.set(xlim=(293.5,318.5),ylim=(0,8.2),ylabel='RMSF (Å)',xlabel='CaMKI residue');a.set_xticks([294,300,306,312,318]);a.set_yticks([0,4,8]);a.legend(ncol=3,loc='upper left',fontsize=9,handlelength=1.2,columnspacing=.7)
supp.text(.53,.96,'b',fontweight='bold',fontsize=11,va='top');supp.text(.59,.96,'Record distributions and mean intervals',fontsize=10,va='top')
a=sa(112,18,68,42)
supp_handles=[Line2D([],[],ls='',marker=m,ms=3.7,mfc=BLUE if fill else 'white',mec=BLUE,color=BLUE,label=name) for name,m,fill in [('Hχ','o',True),('RMSF','^',False)]]
a.legend(handles=supp_handles,ncol=2,loc='lower center',bbox_to_anchor=(.5,1.04),fontsize=9,handletextpad=.35,columnspacing=1.1,borderpad=0,frameon=False)
for metric,base in [('fine_entropy',0.),('lowF_width',1.)]:
    for predictor,offset,marker,filled in [('full_hchi',-.18,'o',True),('full_rmsf',.18,'^',False)]:
        z=corr[corr.metric_key.eq(metric)&corr.predictor.eq(predictor)].sort_values('sample_id');assert len(z)==83
        a.scatter(base+offset+np.linspace(-.075,.075,83),z.rho_type_centered,s=6,marker=marker,facecolors=BLUE if filled else 'none',edgecolors=BLUE,linewidths=.3,alpha=.24)
        q=stats[stats.panel_kind.eq('correlation')&stats.metric_key.eq(metric)&stats.predictor.eq(predictor)].iloc[0]
        a.errorbar(base+offset,q['mean'],yerr=[[q['mean']-q.bootstrap_95_low],[q.bootstrap_95_high-q['mean']]],fmt=marker,ms=4,color=BLUE,mfc=BLUE if filled else 'white',mew=.75,elinewidth=1.1,capsize=2)
a.set(xlim=(-.48,1.48),ylim=(-.06,.75),ylabel='Within-type ρ');a.set_xticks([0,1],['Fine-bin entropy','Low-F* width']);a.set_yticks([0,.25,.5,.75])
for ext in ['png','pdf','svg']:supp.savefig(OUT/f'figS_context_and_record_intervals.{ext}',dpi=450,bbox_inches=supp.bbox_inches)
plt.close(supp)
print(json.dumps(qa,ensure_ascii=False))
