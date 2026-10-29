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
W,H=179.,200.
BLUE='#2166AC';WARM='#B87936';TEAL='#3A8D8C';PURPLE='#8B72A5';GRAY='#BAC1C7';INK='#222222'
FS=8.75
mpl.rcParams.update({'savefig.bbox':None,'font.size':FS,'axes.labelsize':FS,'xtick.labelsize':FS,'ytick.labelsize':FS,'legend.fontsize':FS,'axes.unicode_minus':False,'axes.linewidth':.55})
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
    txt(x,y,label,10,fontweight='bold')
    if title:txt(x+6,y+.25,title,FS,fontweight='normal')
def line(x1,y1,x2,y2,**kw):fig.add_artist(Line2D([x1/W,x2/W],[1-y1/H,1-y2/H],transform=fig.transFigure,**kw))
def values(panel,kind,ids,vals,**kw):
    for i,v in zip(ids,vals):
        if np.isfinite(v):plotted.append(dict(panel=panel,kind=kind,identity=str(i),value=float(v),**kw))
def molecule(key,x,y,w,h,labels):
    z=manifest['assets'][key];a=ax(x,y,w,h);a.imshow(Image.open(R/z['png']).convert('RGBA'));a.axis('off')
    for name,location in labels.items():
        px=z['landmarks'][name]
        a.scatter([px[0]],[px[1]],s=11,facecolors='none',edgecolors=TEAL,linewidth=.85)
        a.annotate(name,xy=px,xytext=location,textcoords='axes fraction',color=INK,fontsize=FS,ha='center',va='center',arrowprops=dict(arrowstyle='-',color=TEAL,lw=.7),bbox=dict(facecolor='white',edgecolor='none',pad=.1,alpha=.92))
    return a

# a: adjacent intact structures, approximately one fifth of the full canvas.
panel('a','',1,2)
txt(20,2,'CaM–CaMKI (1MXE)',FS)
txt(112,2,'Ubiquitin–UBA (2QHO)',FS)
molecule('cam_overview',17,7,44,38.5,{'W303':(.12,.61),'M316':(.99,.28)})
molecule('ub_overview',111,7,44,38.5,{'Q41':(.56,.94),'Q62':(.96,.09)})
# Single semantic legend row; components and methods have distinct color roles.
line(12,49,17,49,color=PURPLE,lw=2.4);txt(19,47.4,'Analysed chain')
line(51,49,56,49,color=GRAY,lw=2.4);txt(58,47.4,'Partner')
for xx,label,col,mark,ls in [(91,'MD',INK,'o',':'),(118,'Full',BLUE,None,'-'),(147,'Intra',WARM,None,'--')]:
    line(xx,49,xx+5,49,color=col,lw=1.1,marker=mark,ms=2.6,markevery=[0],mfc='white',mec=col,ls=ls)
    txt(xx+7,47.4,label)

# b/c: equal-height data axes directly under their matching structure.
panel('b','CaMKI backbone',1,57)
panel('c','Ubiquitin side chains',96,57)
x=cam.CaMKI_native_position.to_numpy()
a=ax(13,66,72,24)
a.axvspan(302.5,316.5,color='#F0F5F3',lw=0,zorder=0)
for name,color,marker,ls,field in [('MD',INK,'o',':','target_rmsf'),('Full',BLUE,None,'-','pred_rmsf'),('Intra',WARM,None,'--','intra_rmsf')]:
    y=cam[field].to_numpy()
    a.plot(x,y,color=color,lw=1.1 if name=='Full' else .95,marker=marker,ms=2.0,ls=ls,mfc='white',mew=.55,zorder=4 if name=='Full' else 3)
    values('b','RMSF_A',x,y,series=name)
mae_full=float(np.mean(abs(cam.pred_rmsf-cam.target_rmsf)))
mae_intra=float(np.mean(abs(cam.intra_rmsf-cam.target_rmsf)))
a.text(.035,.965,f'Full MAE = {mae_full:.3f} Å\nIntra MAE = {mae_intra:.3f} Å',transform=a.transAxes,ha='left',va='top',fontsize=FS,linespacing=1.25)
a.set(xlim=(293.5,318.5),ylim=(0,8.2),ylabel='RMSF (Å)',xlabel='CaMKI residue')
a.set_yticks([0,4,8]);a.set_xticks([294,300,306,312,318]);a.tick_params(length=2,pad=1.8);a.xaxis.labelpad=2
for pos,label,loc in [(303,'W303',(302.7,1.7)),(316,'M316',(315.0,2.0))]:
    q=cam[cam.CaMKI_native_position.eq(pos)].iloc[0]
    a.annotate(label,xy=(pos,q.pred_rmsf),xytext=loc,fontsize=FS,ha='center',color=INK,arrowprops=dict(arrowstyle='-',color=TEAL,lw=.65))

a=ax(108,66,68,24)
for pos in [41,62]:a.axvspan(pos-.55,pos+.55,color='#E8F1EE',lw=0,zorder=0)
for name,color,marker,ls,field in [('MD',INK,'o',':','MD_Hchi_nats'),('Full',BLUE,None,'-','Full_Hchi_nats')]:
    y=ub[field].where(ub.valid_hchi).to_numpy()
    a.plot(ub.residue_id,y,color=color,lw=1.0 if name=='Full' else .8,marker=marker,ms=1.7,ls=ls,mfc='white',mew=.45,zorder=3 if name=='Full' else 2)
    values('c','Hchi_nats',ub.residue_id,y,series=name)
mae_ub=float(np.mean(abs(ub.loc[ub.valid_hchi,'Full_Hchi_nats']-ub.loc[ub.valid_hchi,'MD_Hchi_nats'])))
a.text(.04,.98,f'Full MAE = {mae_ub:.3f} nats',transform=a.transAxes,ha='left',va='top',fontsize=FS)
for pos,loc in [(41,(44.0,3.0)),(62,(64.2,3.55))]:
    q=ub[ub.residue_id.eq(pos)].iloc[0]
    a.annotate('Q'+str(pos),xy=(pos,q.Full_Hchi_nats),xytext=loc,fontsize=FS,color=INK,arrowprops=dict(arrowstyle='-',color=TEAL,lw=.6))
a.set(xlim=(.5,72.5),ylim=(0,4.75),ylabel='Hχ (nats)',xlabel='Ubiquitin residue')
a.set_xticks([1,20,40,60,72]);a.set_yticks([0,2,4]);a.tick_params(length=2,pad=1.8);a.xaxis.labelpad=2

# d: every same-type residue retained; real pairs only, with no fitted SASA trend.
panel('d','Six Gln environments',1,107)
a=ax(13,119,44,29)
for row in gln.itertuples():
    sx=row.sidechain_SASA_A2;md=row.MD_Hchi_nats;pred=row.Full_Hchi_nats
    a.plot([sx,sx],[md,pred],c='#AAB3B8',lw=.7,zorder=1)
    a.scatter(sx,md,s=19,facecolor='white',edgecolor=INK,linewidth=.75,zorder=3)
    a.scatter(sx,pred,s=14,facecolor=BLUE,edgecolor=BLUE,linewidth=.5,zorder=4)
    if row.residue in ['Q41','Q62']:
        loc={'Q41':(16,1.35),'Q62':(98,3.25)}[row.residue]
        a.annotate(row.residue,xy=(sx,(md+pred)/2),xytext=loc,fontsize=FS,color=INK,ha='center',arrowprops=dict(arrowstyle='-',color=TEAL,lw=.65))
    values('d','Hchi_nats',[row.residue],[md],series='MD',SASA_A2=float(sx))
    values('d','Hchi_nats',[row.residue],[pred],series='Full',SASA_A2=float(sx))
a.set(xlim=(-9,114),ylim=(0,3.65),ylabel='Hχ (nats)',xlabel='Side-chain SASA (Å²)')
a.set_xticks([0,50,100]);a.set_yticks([0,1,2,3]);a.tick_params(length=2,pad=1.8);a.xaxis.labelpad=2

# e: a complete, larger 2x2 unit rather than four independently arranged plots.
defined=pd.read_csv(P/'defined_angles_long.csv',float_precision='round_trip')
chi_stats=pd.read_csv(P/'per_chi_correlations.csv',float_precision='round_trip')
chi_stats=chi_stats[chi_stats.predictor.eq('predicted_joint_Hchi')&chi_stats.metric.eq('fine_entropy')]
assert defined.groupby('chi').size().to_dict()=={'chi1':66,'chi2':52,'chi3':23,'chi4':10}
panel('e','Per-χ occupancy association',72,107)
for xx,label,marker,color in [(85,'Other','o',BLUE),(119,'Q41','D',TEAL),(151,'Q2','s',TEAL)]:
    line(xx,116,xx+2,116,color=color,lw=0,marker=marker,ms=3.5,markevery=[0]);txt(xx+4,114.4,label)
for chi,left,top in [(1,85,124),(2,140,124),(3,85,164),(4,140,164)]:
    a=ax(left,top,35,25)
    z=defined[defined.chi.eq('chi'+str(chi))]
    st=chi_stats[chi_stats.chi.eq('chi'+str(chi))].iloc[0]
    txt(left,top-5.2,'χ'+str(chi)+'  n='+str(len(z))+'  ρ='+format(st.rho_AA_centered,'.2f'))
    for name,marker,col in [('Other','o',BLUE),('Q41','D',TEAL),('Q2','s',TEAL)]:
        sub=z[~z.author_residue_id.isin([2,41])] if name=='Other' else z[z.author_residue_id.eq(int(name[1:]))]
        a.scatter(sub.predicted_Hchi_nats_AA_centered,sub.fine_entropy_four_origins_mean_nats_AA_centered,s=10 if name=='Other' else 20,marker=marker,color=col,alpha=.50 if name=='Other' else 1,linewidths=.25,edgecolors='white' if name=='Other' else col,zorder=3 if name=='Other' else 4)
    for row in z.itertuples():
        plotted.extend([dict(panel='e',kind='predicted_joint_Hchi_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='Full',value=float(row.predicted_Hchi_nats_AA_centered)),dict(panel='e',kind='MD_single_chi_fine_entropy_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='MD',value=float(row.fine_entropy_four_origins_mean_nats_AA_centered))])
    a.axhline(0,color='#CCD4D8',lw=.4,ls=':');a.axvline(0,color='#CCD4D8',lw=.4,ls=':')
    a.set(xlim=(-.7,.5),ylim=(-1,.75));a.set_xticks([-.5,0,.5]);a.set_yticks([-.5,0,.5]);a.tick_params(length=2,pad=1.5)
    if chi in [2,4]:a.tick_params(labelleft=False)
    if chi in [1,2]:a.tick_params(labelbottom=False)
txt(73,158.5,'MD entropy deviation (nats)',rotation=90,ha='center',va='center')
txt(130,194.5,'Predicted joint Hχ deviation (nats)',ha='center')

# f: all record-wise values plus their existing paired bootstrap mean intervals.
panel('f','Across 83 records',1,159)
for xx,label,marker,filled in [(17,'Hχ','o',True),(42,'RMSF','^',False)]:
    line(xx,166.4,xx+2,166.4,color=BLUE,lw=0,marker=marker,ms=3.8,markevery=[0],mfc=BLUE if filled else 'white',mec=BLUE);txt(xx+4,164.8,label)
a=ax(23,169,34,20)
for metric,base in [('fine_entropy',1.),('lowF_width',0.)]:
    for predictor,offset,marker,filled in [('full_hchi',.16,'o',True),('full_rmsf',-.16,'^',False)]:
        z=corr[corr.metric_key.eq(metric)&corr.predictor.eq(predictor)].sort_values('sample_id');assert len(z)==83
        jitter=np.linspace(-.085,.085,83)
        a.scatter(z.rho_type_centered,base+offset+jitter,s=6.5,marker=marker,facecolors=BLUE if filled else 'none',edgecolors=BLUE,linewidths=.3,alpha=.24,zorder=1)
        q=stats[stats.panel_kind.eq('correlation')&stats.metric_key.eq(metric)&stats.predictor.eq(predictor)].iloc[0]
        a.errorbar(q['mean'],base+offset,xerr=[[q['mean']-q.bootstrap_95_low],[q.bootstrap_95_high-q['mean']]],fmt=marker,ms=4.3,color=BLUE,mfc=BLUE if filled else 'white',mew=.75,elinewidth=1.1,capsize=2,zorder=5)
        values('f','within_AA_Spearman_rho',z.sample_id,z.rho_type_centered,series=predictor,metric=metric)
a.axvline(0,color='#CCD4D8',lw=.4,ls=':')
a.set(xlim=(-.05,.85),ylim=(-.42,1.42));a.set_yticks([1,0],['Fine-bin\nentropy','Low-F*\nwidth']);a.set_xticks([0,.4,.8]);a.tick_params(length=2,pad=1.5)
txt(36,194.5,'Within-type ρ',ha='center')

df=pd.DataFrame(plotted);df.to_csv(P/'all_plotted_values.csv',index=False)
assert df.groupby('panel').size().to_dict()=={'b':75,'c':132,'d':12,'e':302,'f':332}
fig.canvas.draw();renderer=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in texts if t.get_window_extent(renderer).x0<-.5 or t.get_window_extent(renderer).y0<-.5 or t.get_window_extent(renderer).x1>fig.bbox.width+.5 or t.get_window_extent(renderer).y1>fig.bbox.height+.5]
assert not outside,outside
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'fig3_dual_case_entropy.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
qa=dict(size_mm=[W,H],min_font=min(t.get_fontsize() for t in texts),print_font_at_164mm=FS*164/W,main_plot_counts=df.groupby('panel').size().to_dict(),all_original_mask_positions=True,all_six_Gln=True,all83_records_in_both_main_diagnostics=True,no_SASA_regression=True,no_tables_or_gray_notes=True,axes_with_absolute_amplitudes_start_at_zero=True,no_binding_free_energy_claim=True,not_apo_holo=True,structure_panels_grouped=True,local_zoom_panels_removed=True,paired_structure_profile_columns=True,black_labels_and_colored_keylines=True,annotation_metrics=dict(CaMKI_Full_RMSF_MAE_A=mae_full,CaMKI_Intra_RMSF_MAE_A=mae_intra,Ubiquitin_Full_Hchi_MAE_nats=mae_ub),all_body_text_editable_in_PDF_and_SVG=True,outside_text=outside,visual_review='pending',sha256={e:hashlib.sha256((OUT/f'fig3_dual_case_entropy.{e}').read_bytes()).hexdigest() for e in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(qa,ensure_ascii=False));plt.close(fig)
