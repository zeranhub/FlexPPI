# Production inheritance: V90 native profiles and V80 full-point occupancy plots.
# All data are frozen model outputs/actual MD references; no generated conformations.
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
ROOT=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent
P=OUT/'source_data';R=ROOT/'pymol'
W,H=183.,225.
BLUE='#2166AC';WARM='#B87936';TEAL='#3A8D8C';PURPLE='#8B72A5';GRAY='#BAC1C7';INK='#222222'
mpl.rcParams.update({'savefig.bbox':None,'font.size':9.0,'axes.labelsize':9.0,'xtick.labelsize':9.0,'ytick.labelsize':9.0,'legend.fontsize':9.0,'axes.unicode_minus':False})
cam=pd.read_csv(P/'cam_all25.csv').sort_values('CaMKI_native_position')
ub=pd.read_csv(P/'ubiquitin_72_position_profile_with_original_mask.csv').sort_values('residue_id')
gln=pd.read_csv(P/'ubiquitin_all6_Gln_Hchi_SASA_plot.csv').sort_values('residue_id')
corr=pd.read_csv(P/'retained_absolute_correlations_all498.csv')
groups=pd.read_csv(P/'retained_prediction_groups_all498.csv')
stats=pd.read_csv(P/'retained_panel_statistics.csv')
manifest=json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
assert len(cam)==25 and ub.valid_hchi.sum()==66 and len(ub)==72 and len(gln)==6
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');plotted=[]
def txt(x,y,s,size=9,**kw):return fig.text(x/W,1-y/H,s,fontsize=size,va=kw.pop('va','top'),color=kw.pop('color',INK),**kw)
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def panel(label,title,x,y):txt(x,y,label,11.5,fontweight='bold');txt(x+6,y,title,9.6,fontweight='bold')
def line(x1,y1,x2,y2,**kw):fig.add_artist(Line2D([x1/W,x2/W],[1-y1/H,1-y2/H],transform=fig.transFigure,**kw))
def add_values(panel,kind,ids,values,**kw):
    for i,v in zip(ids,values):
        if np.isfinite(v):plotted.append(dict(panel=panel,kind=kind,identity=str(i),value=float(v),**kw))
def molecule(key,x,y,w,h,labels=None):
    data=manifest['assets'][key];a=ax(x,y,w,h);a.imshow(Image.open(R/data['png']).convert('RGBA'));a.axis('off')
    for name,location in (labels or {}).items():
        px=data['landmarks'][name]
        a.scatter([px[0]],[px[1]],s=11,facecolors='none',edgecolors=TEAL,linewidth=.8)
        a.annotate(name,xy=px,xytext=location,textcoords='axes fraction',color=TEAL,fontsize=9,ha='center',va='center',arrowprops=dict(arrowstyle='-',color=TEAL,lw=.65),bbox=dict(facecolor='white',edgecolor='none',pad=.1,alpha=.9))
    return a

# a/b: only backbone dynamics are used for the CaMKI main-text example.
panel('a','CaM–CaMKI',1,2)
molecule('cam_overview',1,11,58,50.75,{'W303':(.24,.61),'M316':(.88,.35)})
line(6,62.8,11,62.8,color=PURPLE,lw=2.6);txt(13,61.2,'CaMKI',9)
line(36,62.8,41,62.8,color=GRAY,lw=2.6);txt(43,61.2,'CaM',9)
panel('b','CaMKI · backbone RMSF',65,2)
for xx,label,col,mark,ls in [(116,'MD',INK,'o',':'),(136,'Full',BLUE,'o','-'),(157,'Intra',WARM,'s','--')]:
    line(xx,10.5,xx+4.5,10.5,color=col,lw=1.15,marker=mark,ms=2.8,ls=ls);txt(xx+6,8.9,label,9)
x=cam.CaMKI_native_position.to_numpy();lo,hi=293.5,318.5
for pos,aa in zip(x,cam.amino_acid):txt(75+(pos-lo)/25*105,15.1,aa,9,ha='center',fontweight='bold' if pos in [303,316] else 'normal')
a=ax(75,20,105,27);a.axvspan(302.5,316.5,color='#EDF3F2',zorder=0)
for label,col,mark,ls,field in [('MD',INK,'o',':','target_rmsf'),('Full',BLUE,'o','-','pred_rmsf'),('Intra',WARM,'s','--','intra_rmsf')]:
    vals=cam[field].to_numpy();a.plot(x,vals,color=col,lw=1.05,marker=mark,ms=2.4,ls=ls,mfc='white' if label=='MD' else col,mew=.6,zorder=4 if label=='Full' else 3)
    add_values('b','RMSF_A',x,vals,series=label)
a.set(xlim=(lo,hi),ylim=(0,8.2),ylabel='RMSF (Å)');a.set_yticks([0,4,8]);a.set_xticks([]);a.tick_params(length=2,pad=2)
for pos,label,loc in [(303,'W303',(303.2,1.8)),(316,'M316',(315.5,2.4))]:
    q=cam[cam.CaMKI_native_position.eq(pos)].iloc[0]
    a.annotate(label,xy=(pos,q.pred_rmsf),xytext=loc,fontsize=9,ha='center',color=TEAL,arrowprops=dict(arrowstyle='-',color=TEAL,lw=.6))
a=ax(75,53,105,7);a.bar(x,cam.n_CaM_contact_residues4_5,width=.8,color='#A6B0B6',linewidth=0)
a.set(xlim=(lo,hi),ylim=(0,10.5),ylabel='Contacts');a.set_yticks([0,10]);a.set_xticks([294,298,303,308,312,316,318]);a.tick_params(length=2,pad=1.5)
a.set_xlabel('CaMKI residue number',labelpad=1.5)
add_values('b','CaM_contact_residues_4p5A',x,cam.n_CaM_contact_residues4_5.to_numpy(),series='Static structure')

# c: all 72 positions / 66 original effective side-chain labels; A/G remain gaps.
panel('c','Ubiquitin · complete side-chain profile',1,74)
for xx,label,col,mark,ls in [(138,'MD',INK,'o',':'),(158,'Full',BLUE,'s','-')]:
    line(xx,77,xx+4.5,77,color=col,lw=1.05,marker=mark,ms=2.4,ls=ls);txt(xx+6,75.4,label,9)
a=ax(16,85,164,25)
for pos in [41,62]:a.axvspan(pos-.55,pos+.55,color='#E4EFED',lw=0,zorder=0)
for label,col,mark,ls,field in [('MD',INK,'o',':','MD_Hchi_nats'),('Full',BLUE,'s','-','Full_Hchi_nats')]:
    vals=ub[field].where(ub.valid_hchi).to_numpy();a.plot(ub.residue_id,vals,color=col,lw=.95,marker=mark,ms=2.1,ls=ls,mfc='white' if label=='MD' else col,mew=.55,zorder=3 if label=='Full' else 2)
    add_values('c','Hchi_nats',ub.residue_id,vals,series=label)
for pos,loc in [(41,(42.7,3.05)),(62,(63.7,3.5))]:
    q=ub[ub.residue_id.eq(pos)].iloc[0];a.annotate('Q'+str(pos),xy=(pos,q.Full_Hchi_nats),xytext=loc,fontsize=9,color=TEAL,arrowprops=dict(arrowstyle='-',color=TEAL,lw=.6))
a.set(xlim=(.5,72.5),ylim=(0,4.3),ylabel='Hχ (nats)',xlabel='Ubiquitin residue number');a.set_xticks([1,10,20,30,40,50,60,72]);a.set_yticks([0,2,4]);a.tick_params(length=2,pad=2);a.xaxis.labelpad=2

# d/e/f: true structure and all same-type positions; no SASA regression is implied.
panel('d','Ubiquitin · local environments',1,121);panel('e','All six Gln positions',118,121)
molecule('ub_overview',0,132,45,39.375,{'Q41':(.56,.91),'Q62':(.88,.08)})
txt(7,168,'UBA',9);txt(31,168,'Ub',9)
for key,x0,label in [('ub_Q41_local',48,'Q41'),('ub_Q62_local',81,'Q62')]:
    molecule(key,x0,132,31,27.125)
    txt(x0+15.5,128,label,9,ha='center',color=TEAL,fontweight='bold')
a=ax(129,135,50,28)
label_locations={'Q2':(91,3.4),'Q31':(91,1.9),'Q40':(48,2.2),'Q41':(14,1.35),'Q49':(43,3.35),'Q62':(104,3.1)}
for row in gln.itertuples():
    sx=row.sidechain_SASA_A2;md=row.MD_Hchi_nats;pred=row.Full_Hchi_nats
    a.plot([sx,sx],[md,pred],c='#9DA6AB',lw=.8,zorder=1)
    a.scatter(sx,md,s=18,facecolor='white',edgecolor=INK,linewidth=.8,zorder=3)
    a.scatter(sx,pred,s=15,marker='s',facecolor=BLUE,edgecolor=BLUE,linewidth=.6,zorder=4)
    a.annotate(row.residue,xy=(sx,(md+pred)/2),xytext=label_locations[row.residue],fontsize=9,color=TEAL if row.residue in ['Q41','Q62'] else INK,ha='center',arrowprops=dict(arrowstyle='-',color='#757F85',lw=.5))
    add_values('e','Hchi_nats',[row.residue],[md],series='MD',SASA_A2=float(sx));add_values('e','Hchi_nats',[row.residue],[pred],series='Full',SASA_A2=float(sx))
a.set(xlim=(-9,114),ylim=(0,3.65),ylabel='Hχ (nats)',xlabel='Side-chain SASA (Å²)');a.set_xticks([0,50,100]);a.set_yticks([0,1,2,3]);a.tick_params(length=2,pad=2);a.xaxis.labelpad=2

# g: absolute type-controlled correlations; filled circles / open triangles distinguish readout.
panel('f','Occupancy association',1,181)
for xx,label,mark,filled in [(13,'Hχ','o',True),(45,'RMSF','^',False)]:
    line(xx,190,xx+3,190,color=BLUE,lw=0,marker=mark,ms=4,mfc=BLUE if filled else 'white',mec=BLUE);txt(xx+5,188.3,label,9)
a=ax(29,194,56,20)
for i,(metric,base) in enumerate([('fine_entropy',1.),('lowF_width',0.)]):
    for predictor,offset,marker,filled in [('full_hchi',.16,'o',True),('full_rmsf',-.16,'^',False)]:
        z=corr[corr.metric_key.eq(metric)&corr.predictor.eq(predictor)].sort_values('sample_id');assert len(z)==83
        jitter=np.linspace(-.085,.085,83)
        a.scatter(z.rho_type_centered,base+offset+jitter,s=7,marker=marker,facecolors=BLUE if filled else 'none',edgecolors=BLUE,linewidths=.35,alpha=.18,zorder=1)
        q=stats[stats.panel_kind.eq('correlation')&stats.metric_key.eq(metric)&stats.predictor.eq(predictor)].iloc[0]
        a.errorbar(q['mean'],base+offset,xerr=[[q['mean']-q.bootstrap_95_low],[q.bootstrap_95_high-q['mean']]],fmt=marker,ms=4.6,color=BLUE,mfc=BLUE if filled else 'white',mew=.8,elinewidth=1.2,capsize=2.2,zorder=5)
        add_values('f','within_AA_Spearman_rho',z.sample_id,z.rho_type_centered,series=predictor,metric=metric)
a.axvline(0,color='#BFC5C9',lw=.6,ls=':');a.set(xlim=(-.05,.85),ylim=(-.42,1.42),xlabel='Within-type Spearman ρ');a.set_yticks([1,0],['Fine-bin\nentropy','Low-F*\nwidth']);a.set_xticks([0,.4,.8]);a.tick_params(length=2,pad=2);a.xaxis.labelpad=2

# h: prediction selection is tested against actual MD; all 83 paired records remain visible.
panel('g','Prediction-selected halves',96,181)
for metric,x0,title,ylim,ticks in [('fine_entropy',108,'Entropy (nats/χ)',(0,3.4),[0,1.5,3]),('lowF_width',152,'Width (°/χ)',(0,60),[0,30,60])]:
    txt(x0+13.5,190,title,9,ha='center')
    a=ax(x0,197,27,20.5);z=groups[groups.metric_key.eq(metric)&groups.predictor.eq('full_hchi')].sort_values('sample_id');assert len(z)==83
    for j,row in enumerate(z.itertuples()):
        jitter=((j*37)%83/82-.5)*.19
        a.plot([0+jitter,1+jitter],[row.low_MD_mean,row.high_MD_mean],color='#B6BDC1',lw=.4,alpha=.45,zorder=1)
    for j,field in enumerate(['low_MD_mean','high_MD_mean']):
        jitter=np.array([((i*37)%83/82-.5)*.19 for i in range(83)])
        a.scatter(j+jitter,z[field],color=INK,s=4,alpha=.20,linewidths=0,zorder=2)
        q=stats[stats.panel_kind.eq('prediction_selection')&stats.metric_key.eq(metric)&stats.predictor.eq('full_hchi')&stats.value_key.eq(field)].iloc[0]
        a.errorbar(j,q['mean'],yerr=[[q['mean']-q.bootstrap_95_low],[q.bootstrap_95_high-q['mean']]],fmt='D',ms=4.4,color=INK,mfc='white',mew=.8,elinewidth=1.1,capsize=2.1,zorder=5)
        add_values('g','MD_'+metric,z.sample_id,z[field],series=field)
    assert z[['low_MD_mean','high_MD_mean']].to_numpy().min()>=ylim[0] and z[['low_MD_mean','high_MD_mean']].to_numpy().max()<=ylim[1]
    a.set(xlim=(-.3,1.3),ylim=ylim);a.set_xticks([0,1],['Low','High']);a.set_yticks(ticks);a.tick_params(length=2,pad=2)

df=pd.DataFrame(plotted);df.to_csv(P/'all_plotted_values.csv',index=False)
assert df.groupby('panel').size().to_dict()=={'b':100,'c':132,'e':12,'f':332,'g':332}
fig.canvas.draw();renderer=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in texts if t.get_window_extent(renderer).x0<-.5 or t.get_window_extent(renderer).y0<-.5 or t.get_window_extent(renderer).x1>fig.bbox.width+.5 or t.get_window_extent(renderer).y1>fig.bbox.height+.5]
assert not outside,outside
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'fig3_dual_case_entropy.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
qa=dict(size_mm=[W,H],min_font=min(t.get_fontsize() for t in texts),print_font_at_164mm=min(t.get_fontsize() for t in texts)*164/183,all_original_mask_positions=True,CaM_main_metric='RMSF only',Ub_main_metric='Hchi only',main_plot_counts=df.groupby('panel').size().to_dict(),no_SASA_regression=True,all_six_Gln=True,all83_records_in_both_main_diagnostics=True,no_tables_or_gray_notes=True,axes_with_absolute_amplitudes_start_at_zero=True,no_binding_free_energy_claim=True,not_apo_holo=True,body_text_outside=[],outside_text=outside,visual_review='pending',sha256={e:hashlib.sha256((OUT/f'fig3_dual_case_entropy.{e}').read_bytes()).hexdigest() for e in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(qa,ensure_ascii=False));plt.close(fig)
