from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image

OUT=Path(__file__).resolve().parent; P=OUT/'source_data'; R=OUT/'structure_assets'
W,H=183.,190.; STATE_A='#B87936'; STATE_B='#2166AC'; KEY='#3A8D8C'
N_GREY='#8A9299'; C_GREY='#BAC1C7'; PARTNER='#8B72A5'
STATES=[('NHE1','a',STATE_A,'o'),('IQCG','b',STATE_B,'s')]
mpl.rcParams.update({'font.size':8,'axes.labelsize':7.7,'axes.titlesize':8.5,'xtick.labelsize':7.2,'ytick.labelsize':7.2,'legend.fontsize':7.8,'savefig.bbox':None,'axes.unicode_minus':False})
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white'); axes=[];points=[];checks=[]
def text(x,y,s,fs=8,**kw):return fig.text(x/W,1-y/H,s,fontsize=fs,va=kw.pop('va','top'),color=kw.pop('color',BLACK),**kw)
def axis(x,y,w,h):
    a=fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);axes.append(a);return a
def panel(letter,title,x,y):text(x,y,letter,11,fontweight='bold');text(x+6,y,title,8.5,fontweight='bold')
def snap(key,x,y,w,h):
    info=manifest['assets'][key];info=info.get('paired_trim',info);a=axis(x,y,w,h);a.imshow(Image.open(R/info['file']).convert('RGBA'));a.axis('off')
    if key.startswith('whole'):
        for k,loc in [('M72',(.12,.75)),('M109',(.88,.30))]:
            a.annotate(k,xy=info['landmarks'][k],xytext=loc,textcoords='axes fraction',fontsize=7.2,color=KEY,ha='center',va='center',bbox=dict(fc='white',ec='none',alpha=.85,pad=.4),arrowprops=dict(arrowstyle='-',color=BLACK,lw=.5,shrinkA=1,shrinkB=1))
    return a
def check(a,values):
    z=np.concatenate([np.asarray(x)[np.isfinite(x)] for x in values]);lo,hi=a.get_ylim();assert z.min()>=lo and z.max()<=hi
    checks.append(dict(label=a.get_ylabel(),limits=[lo,hi],range=[float(z.min()),float(z.max())]))

data=pd.read_csv(P/'cam_md_state_contrast_all_common_residues.csv');chi=data[data.metric=='hchi'].sort_values('mature_cam_number');back=data[data.metric=='rmsf'].sort_values('mature_cam_number');mets=chi[chi.amino_acid=='M'];assert (len(back),len(chi),len(mets))==(145,122,9)
manifest=json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
# Displayed predictions are frozen source values. No fit, smoothing or renormalisation.
for _,g in data.groupby('metric'):np.testing.assert_allclose(g.full0_b-g.full0_a,g.full0_delta,atol=1e-12)

panel('a','Recognition structures',1,2)
for i,(state,_,_,_) in enumerate(STATES):
    x=1+i*50;text(x+24.5,10,state+'–CaM',8.4,ha='center',fontweight='bold');snap('whole_'+state,x,14,49,47)
for x,label,c in [(10,'CaM',N_GREY),(38,'Partner',PARTNER),(70,'Key Met',KEY)]:
    fig.add_artist(plt.Line2D([(x-4)/W,(x-1)/W],[1-68.9/H]*2,transform=fig.transFigure,color=c,lw=3));text(x,67.5,label,7.7)
panel('b','Local environments',105,2)
for pos,top in [(72,13),(109,45)]:
    text(106,top,'M'+str(pos),8,fontweight='bold',color=KEY)
    for j,(state,_,_,_) in enumerate(STATES):
        snap('detail_focus_M'+str(pos)+'_'+state,106+j*37,top+5,35,24);text(123.5+j*37,top,state,7.5,ha='center')

# Shared state legend belongs only to the quantitative panels.
for x,(state,s,c,marker) in zip([132,158],STATES):
    fig.add_artist(plt.Line2D([x/W],[1-78.3/H],transform=fig.transFigure,color=c,marker=marker,ms=3.7,lw=0));text(x+3,76.8,state,7.8)
panel('d','Contact groups',1,83);panel('e','All nine Met',64,83);panel('f','Key Met sites',128,83)

masks=[(chi.contact_a>0)&(chi.contact_b==0),(chi.contact_a==0)&(chi.contact_b>0),(chi.contact_a>0)&(chi.contact_b>0),(chi.contact_a==0)&(chi.contact_b==0)]
labels=['NHE1\nonly (8)','IQCG\nonly (21)','Both\n(3)','Neither\n(90)']
a=axis(13,94,43,26);ys=[]
for j,mask in enumerate(masks):
    q=chi.loc[mask];assert len(q)==[8,21,3,90][j];vs=[]
    for dx,(state,s,c,marker) in zip([-.14,.14],STATES):
        val=float(q['full0_'+s].mean());vs.append(val);a.scatter(j+dx,val,s=23,marker=marker,color=c,zorder=3)
        points.append(dict(panel='d',group=j,kind='full0',state=state,value=val))
    a.plot([j-.14,j+.14],vs,color=BLACK,lw=.65,zorder=1);ys.append(vs)
a.set(xlim=(-.45,3.45),ylim=(1.3,2.2),ylabel='Mean Hχ (nats)');a.set_yticks([1.4,1.8,2.2]);a.set_xticks(range(4),labels,fontsize=7);a.tick_params(length=2,pad=2);check(a,ys)

a=axis(75,94,47,26);ys=[]
for state,s,c,marker in STATES:
    vals=mets['full0_'+s].to_numpy();ys.append(vals);a.plot(range(9),vals,color=c,marker=marker,ms=2.7,lw=.9)
    for j,(_,q) in enumerate(mets.iterrows()):
        points.append(dict(panel='e',position=int(q.mature_cam_number),state=state,kind='full0',metric='hchi',value=float(q['full0_'+s])))
        if int(q.mature_cam_number) in [72,109,144]:a.plot(j,q['full0_'+s],marker=marker,ms=3.3,mfc=c,mec=BLACK,mew=.5,zorder=4)
a.set(xlim=(-.4,8.4),ylim=(1.3,2.95),ylabel='Hχ (nats)');a.set_yticks([1.5,2,2.5]);a.set_xticks(range(9),['M'+str(int(x)) for x in mets.mature_cam_number],fontsize=6.9,rotation=55);a.tick_params(length=2,pad=2);check(a,ys)
for lab,pos in zip(a.get_xticklabels(),mets.mature_cam_number):
    if pos in [72,109,144]:lab.set_fontweight('bold')

for metric,x,lims,ticks in [('rmsf',137,(.5,1.62),[.6,1,1.4]),('hchi',166,(1.45,2.35),[1.5,1.9,2.3])]:
    a=axis(x,94,17 if metric=='rmsf' else 16,26);ys=[]
    for j,pos in enumerate([72,109]):
        q=data[(data.metric==metric)&(data.mature_cam_number==pos)].iloc[0];va,vb=q.full0_a,q.full0_b;a.plot([j-.12,j+.12],[va,vb],color=BLACK,lw=.65,zorder=1);ys.append([va,vb])
        for dx,(state,s,c,marker) in zip([-.12,.12],STATES):
            val=q['full0_'+s];a.scatter(j+dx,val,s=22,marker=marker,color=c,zorder=3)
            dy=4 if val==max(va,vb) else -5
            a.annotate(f'{val:.2f}',(j+dx,val),xytext=(0,dy),textcoords='offset points',ha='center',va='bottom' if dy>0 else 'top',fontsize=6.9,color=BLACK)
            points.append(dict(panel='f',position=pos,metric=metric,kind='full0',state=state,value=float(val)))
    a.set(xlim=(-.65,1.65),ylim=lims,ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)');a.set_yticks(ticks);a.set_xticks([0,1],['M72','M109'],fontsize=7);a.tick_params(length=2,pad=2);a.yaxis.labelpad=2;check(a,ys)

panel('c','Predicted residue profiles',1,136)
lo,hi=3.5,148.5;regions=[(68,76,72),(105,113,109),(140,148,144)]
a=axis(20,143,157,2.5);a.set(xlim=(lo,hi),ylim=(0,1));a.axis('off')
for l,h,c,name in [(4,75,N_GREY,'N-lobe'),(76,81,'#D9DDE0',''),(82,148,C_GREY,'C-lobe')]:
    a.add_patch(Rectangle((l-.5,0),h-l+1,1,fc=c,ec='none',alpha=.55))
    if name:a.text((l+h)/2,.5,name,ha='center',va='center',fontsize=7.4,color=BLACK)
for l,h,pos in regions:
    px=20+(pos-lo)/(hi-lo)*157;text(px,147,'M'+str(pos),7.2,ha='center',fontweight='bold')
for metric,y in [('rmsf',151),('hchi',166)]:
    a=axis(20,y,157,10);g=data[data.metric==metric].set_index('mature_cam_number').reindex(range(4,149));ylim=4 if metric=='rmsf' else 4.1
    for l,h,pos in regions:
        a.add_patch(Rectangle((l-.5,.04),h-l+1,ylim-.08,fc='#F5F5F5',ec=BLACK,lw=.6,zorder=0))
    for state,s,c,marker in STATES:
        values=g['full0_'+s];a.plot(g.index,values,color=c,lw=.95)
        for pos in [72,109,144]:a.plot(pos,g.loc[pos,'full0_'+s],marker=marker,ms=3.2,mfc=c,mec=BLACK,mew=.45,zorder=4)
        points.extend(dict(panel='c',metric=metric,position=int(p),kind='full0',state=state,value=float(val)) for p,val in values.dropna().items())
    a.set(xlim=(lo,hi),ylim=(0,ylim),ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)');a.set_yticks([0,2,4]);a.tick_params(axis='y',length=2,pad=2)
    if metric=='rmsf':a.tick_params(axis='x',bottom=False,labelbottom=False);a.spines['bottom'].set_visible(False)
    else:a.set_xticks([4,20,40,60,80,100,120,140,148]);a.tick_params(axis='x',length=2,pad=2);a.set_xlabel('Mature CaM residue number',labelpad=2)
    check(a,[g.full0_a,g.full0_b])

# Audit the full prediction-only panel census and all highlighted values.
plotted=pd.DataFrame(points);assert plotted.groupby('panel').size().to_dict()=={'c':534,'d':8,'e':18,'f':8};assert (plotted.kind=='full0').all()
for _,q in plotted[plotted.panel.isin(['c','e','f'])].iterrows():
    row=data[(data.metric==q.metric)&(data.mature_cam_number==q.position)].iloc[0];s='a' if q.state=='NHE1' else 'b';assert abs(q.value-row['full0_'+s])<1e-12
fig.canvas.draw();renderer=fig.canvas.get_renderer();visible=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in visible if t.get_window_extent(renderer).x0<-.6 or t.get_window_extent(renderer).y0<-.6 or t.get_window_extent(renderer).x1>fig.bbox.width+.6 or t.get_window_extent(renderer).y1>fig.bbox.height+.6];assert not outside,outside
assert not any('Δ' in a.get_ylabel() for a in axes)
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig3_CaM_model_story.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
plotted.to_csv(P/'all_plotted_model_story_quantities.csv',index=False)
metrics={}
for metric,g in data.groupby('metric'):
    metrics[metric]={'positions':len(g),'mean_NHE1':float(g.full0_a.mean()),'mean_IQCG':float(g.full0_b.mean()),'increased':int((g.full0_delta>0).sum()),'decreased':int((g.full0_delta<0).sum())}
(P/'model_story_summary.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2),encoding='utf-8')
qa=dict(canvas_mm=[W,H],prediction_only=True,no_MD_values_in_main=True,no_accuracy_scatter=True,no_numeric_tables=True,no_gray_note_text=True,panel_layout='a/b top; d/e/f middle; c bottom',key_regions=regions,state_palette={'NHE1':STATE_A,'IQCG':STATE_B},structure_palette={'CaM_N':N_GREY,'CaM_C':C_GREY,'partner':PARTNER,'key_Met':KEY},plotted_rows=plotted.groupby('panel').size().to_dict(),outside_text=outside,axis_checks=checks,min_font_pt=min(t.get_fontsize() for t in visible),visual_review='pending',output_sha256={ext:hashlib.sha256((OUT/('fig3_CaM_model_story.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False));plt.close(fig)