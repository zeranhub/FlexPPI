from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image

OUT=Path(__file__).resolve().parent; P=OUT/'source_data'
R=OUT/'structure_assets'
W,H=183.,235.; BLUE='#2166AC'; WARM='#B87936'; TEAL='#3A8D8C'; MD='#626A70'
STATES=[('NHE1','a',WARM,'o'),('IQCG','b',BLUE,'s')]
mpl.rcParams.update({'font.size':8,'axes.labelsize':7.5,'axes.titlesize':8,'xtick.labelsize':7,
    'ytick.labelsize':7,'legend.fontsize':7,'axes.unicode_minus':False,'savefig.bbox':None})
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');axes=[];points=[];checks=[]
def text(x,y,s,fs=8,**kw):return fig.text(x/W,1-y/H,s,fontsize=fs,va=kw.pop('va','top'),color=kw.pop('color',BLACK),**kw)
def axis(x,y,w,h):
    a=fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);axes.append(a);return a
def panel(l,s,x,y):text(x,y,l,11,fontweight='bold');text(x+6,y,s,8.5,fontweight='bold')
def snap(key,x,y,w,h):
    info=manifest['assets'][key];info=info.get('paired_trim',info)
    a=axis(x,y,w,h);a.imshow(Image.open(R/info['file']).convert('RGBA'));a.axis('off')
    for k,loc in [('M72',(.12,.75)),('M109',(.88,.30))] if key.startswith('whole') else [(key.split('_')[2],(.5,.1))]:
        a.annotate(k,xy=info['landmarks'][k],xytext=loc,textcoords='axes fraction',fontsize=6.7,color=TEAL,
            ha='center',va='center',bbox=dict(fc='white',ec='none',alpha=.85,pad=.4),
            arrowprops=dict(arrowstyle='-',color='#777777',lw=.5,shrinkA=1,shrinkB=1))
    return a
def check(a,values):
    v=np.concatenate([np.asarray(z)[np.isfinite(z)] for z in values]);lo,hi=a.get_ylim()
    assert v.min()>=lo and v.max()<=hi
    checks.append(dict(label=a.get_ylabel(),limits=[lo,hi],range=[float(v.min()),float(v.max())]))

data=pd.read_csv(P/'cam_md_state_contrast_all_common_residues.csv')
blocks=pd.read_csv(P/'cam_md_state_contrast_all_residues_five_blocks.csv')
chi=data[data.metric=='hchi'].sort_values('mature_cam_number');back=data[data.metric=='rmsf'].sort_values('mature_cam_number')
mets=chi[chi.amino_acid=='M'];assert len(chi)==122 and len(back)==145 and len(mets)==9
manifest=json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
for _,g in data.groupby('metric'):
    np.testing.assert_allclose(g.md_b-g.md_a,g.md_delta,atol=1e-12)
    np.testing.assert_allclose(g.full0_b-g.full0_a,g.full0_delta,atol=1e-12)

panel('a','Recognition structures',1,2)
for i,(state,_,_,_) in enumerate(STATES):
    x=1+i*50;text(x+24.5,10,state+'–CaM',8.2,ha='center',fontweight='bold')
    snap('whole_'+state,x,14,49,47)

text(10,68,'CaM N-lobe',7.4,color=BLUE);text(43,68,'C-lobe',7.4,color='#85ABD1');text(67,68,'Partner',7.4,color=WARM)

panel('b','Local environments',105,2)
for pos,top in [(72,13),(109,45)]:
    text(106,top,'M'+str(pos),7.7,fontweight='bold',color=TEAL)
    for j,(state,_,c,_) in enumerate(STATES):
        snap('detail_focus_M'+str(pos)+'_'+state,106+j*37,top+5,35,24)
        text(123.5+j*37,top,state,7.3,ha='center',color=c)

panel('c','Backbone and side-chain flexibility',1,80)
text(127,80,'● NHE1',7.5,color=WARM);text(153,80,'■ IQCG',7.5,color=BLUE)
lo,hi=3.5,148.5;a=axis(20,88,157,3);a.set(xlim=(lo,hi),ylim=(0,1));a.axis('off')
for l,h,c,name in [(4,75,BLUE,'N-lobe'),(76,81,'#A6AEB4',''),(82,148,'#85ABD1','C-lobe')]:
    a.add_patch(Rectangle((l-.5,0),h-l+1,1,fc=c,ec='none',alpha=.5))
    if name:a.text((l+h)/2,.5,name,ha='center',va='center',fontsize=7,color='#384E61')
for metric,kind,y in [('rmsf','md',95),('rmsf','full0',111),('hchi','md',130),('hchi','full0',146)]:
    a=axis(20,y,157,12);g=data[data.metric==metric].set_index('mature_cam_number').reindex(range(4,149))
    text(21,y-3,'MD reference' if kind=='md' else 'Full prediction',7.2,color=BLACK)
    for l,h in [(20,31),(56,67),(93,104),(129,140)]:a.axvspan(l-.5,h+.5,color=TEAL,alpha=.065,lw=0)
    for state,s,c,_ in STATES:
        a.plot(g.index,g[kind+'_'+s],color=c,lw=.9)
        points.extend(dict(panel='c',metric=metric,kind=kind,position=int(p),state=state,value=float(v)) for p,v in g[kind+'_'+s].dropna().items())
    a.set(xlim=(lo,hi),ylim=(0,4.1 if metric=='hchi' else 4),ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)')
    a.set_yticks([0,2,4]);a.tick_params(axis='x',labelbottom=False,bottom=False);a.tick_params(axis='y',pad=2,length=2);a.spines['bottom'].set_visible(False)
    check(a,[g[kind+'_a'],g[kind+'_b']])
a=axis(20,162,157,5);a.set(xlim=(lo,hi),ylim=(-.5,1.5))
for j,(col,c) in enumerate([('contact_b',BLUE),('contact_a',WARM)]):
    for p,v in back.set_index('mature_cam_number')[col].items():
        if v:a.add_patch(Rectangle((p-.5,j-.35),1,.70,fc=c,ec='none'))
a.set_yticks([0,1],['IQCG','NHE1']);a.tick_params(axis='y',length=0,pad=4,labelsize=6.8)
a.set_xticks([4,20,40,60,80,100,120,140,148]);a.tick_params(axis='x',length=2,pad=3);a.set_xlabel('Mature CaM residue number',labelpad=2);a.spines['left'].set_visible(False)

panel('d','Contact groups',1,185);panel('e','All nine Met',78,185);panel('f','Hχ prediction',143,185)
masks=[(chi.contact_a>0)&(chi.contact_b==0),(chi.contact_a==0)&(chi.contact_b>0),(chi.contact_a>0)&(chi.contact_b>0),(chi.contact_a==0)&(chi.contact_b==0)]
labels=['NHE1\nonly\nn=8','IQCG\nonly\nn=21','Both\n\nn=3','Neither\n\nn=90']
for kind,x in [('md',13),('full0',49)]:
    text(x+13.5,193,'MD reference' if kind=='md' else 'Full prediction',7.1,ha='center',color=BLACK)
    a=axis(x,198,27,24);ys=[]
    for j,mask in enumerate(masks):
        g=chi.loc[mask];assert len(g)==[8,21,3,90][j];vals=[]
        for offset,(state,s,c,marker) in zip([-.15,.15],STATES):
            v=float(g[kind+'_'+s].mean());vals.append(v)
            a.scatter(j+offset,v,s=18,marker=marker,color=c,zorder=3)
            if kind=='md':
                b=blocks[(blocks.metric=='hchi')&blocks.mature_cam_number.isin(g.mature_cam_number)].groupby('block')['md_'+s].mean()
                a.vlines(j+offset,b.min(),b.max(),color=c,alpha=.4,lw=.8)
            points.append(dict(panel='d',group=j,state=state,kind=kind,value=v))
        a.plot([j-.15,j+.15],vals,color='#B6BEC5',lw=.7,zorder=1);ys.append(vals)
    a.set(xlim=(-.4,3.4),ylim=(1.3,2.2));a.set_yticks([1.4,1.8,2.2]);a.set_xticks(range(4),labels,fontsize=6.2);a.tick_params(pad=2,length=2)
    if kind=='md':a.set_ylabel('Mean Hχ (nats)',labelpad=2)
    else:a.tick_params(labelleft=False)
    check(a,ys)
for kind,x in [('md',88),('full0',118)]:
    text(x+11,193,'MD reference' if kind=='md' else 'Full prediction',7.1,ha='center',color=BLACK)
    a=axis(x,198,22,24);ys=[]
    for state,s,c,marker in STATES:
        values=mets[kind+'_'+s];ys.append(values);a.plot(range(9),values,color=c,marker=marker,ms=2.2,lw=.8)
        for j,(_,q) in enumerate(mets.iterrows()):
            if kind=='md':
                b=blocks[(blocks.metric=='hchi')&(blocks.mature_cam_number==q.mature_cam_number)]['md_'+s]
                a.vlines(j,b.min(),b.max(),color=c,alpha=.35,lw=.8)
            points.append(dict(panel='e',position=int(q.mature_cam_number),state=state,kind=kind,value=float(q[kind+'_'+s])))
    a.set(xlim=(-.4,8.4),ylim=(1.0,3.4));a.set_yticks([1,2,3]);a.set_xticks(range(9),[str(int(v)) for v in mets.mature_cam_number],fontsize=6.1,rotation=60);a.tick_params(pad=2,length=2)
    if kind=='md':a.set_ylabel('Hχ (nats)',labelpad=2)
    else:a.tick_params(labelleft=False)
    check(a,ys)
text(112,232,'Mature CaM Met position',6.8,color=BLACK,ha='center')


a=axis(153,198,27,24);a.plot([0,4.1],[0,4.1],color='#B3B8BC',ls='--',lw=.6)
for state,s,c,marker in STATES:
    a.scatter(chi['md_'+s],chi['full0_'+s],s=6,marker=marker,color=c,alpha=.58,edgecolors='white',linewidths=.18)
    points.extend(dict(panel='f',position=int(q.mature_cam_number),state=state,MD=q['md_'+s],Full=q['full0_'+s]) for _,q in chi.iterrows())
a.set(xlim=(0,4.1),ylim=(0,4.1),xlabel='MD Hχ (nats)',ylabel='Predicted Hχ (nats)');a.set_xticks([0,2,4]);a.set_yticks([0,2,4]);a.tick_params(pad=2,length=2,labelsize=6.5);a.set_xlabel('MD Hχ (nats)',fontsize=6.6,labelpad=2);a.set_ylabel('Predicted Hχ (nats)',fontsize=6.6,labelpad=2)
fig.canvas.draw();renderer=fig.canvas.get_renderer();visible=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in visible if t.get_window_extent(renderer).x0<-.6 or t.get_window_extent(renderer).y0<-.6 or t.get_window_extent(renderer).x1>fig.bbox.width+.6 or t.get_window_extent(renderer).y1>fig.bbox.height+.6]
assert not outside,outside
assert not any('Δ' in a.get_ylabel() for a in axes)
assert not any(t.get_color()==MD for t in visible), 'Gray explanatory text must be in caption, not figure'
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig3_CaM_shared_view.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
pd.DataFrame(points).to_csv(P/'all_plotted_absolute_quantities.csv',index=False)
qa=dict(canvas_mm=[W,H],all_axes_absolute=True,two_state_profiles=True,local_views_right_of_whole=True,no_numeric_tables=True,no_gray_note_text=True,structural_alignment='N-lobe helices I/IV; one whole-complex rigid transform',all_residue_counts=dict(backbone=145,sidechain=122,Met=9,contact_groups=[8,21,3,90]),
    plotted_rows=pd.DataFrame(points).groupby('panel').size().to_dict(),outside_text=outside,min_font_pt=min(t.get_fontsize() for t in visible),axis_checks=checks,visual_review='pending',
    output_sha256={ext:hashlib.sha256((OUT/('fig3_CaM_shared_view.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False));plt.close(fig)
