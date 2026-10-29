from pathlib import Path
import json, hashlib, shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from PIL import Image

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent
P=OUT/'source_data' if (OUT/'source_data').exists() else ROOT.parent/'partner_screen'
R=ROOT/'pymol'
W,H=183.,245.
BLUE='#2166AC'; WARM='#B87936'; TEAL='#3A8D8C'; MD='#626A70'; PALE='#85ABD1'
mpl.rcParams.update({'font.size':8,'axes.labelsize':8,'axes.titlesize':8,
    'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.5,
    'axes.unicode_minus':False,'savefig.bbox':None})
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
all_axes=[]; sources=[]

def text(x,y,s,fs=8,**kw):
    return fig.text(x/W,1-y/H,s,fontsize=fs,va=kw.pop('va','top'),color=kw.pop('color',BLACK),**kw)
def axis(x,y,w,h):
    a=fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);all_axes.append(a);return a
def panel(l,title,x,y):
    text(x,y,l,11,fontweight='bold');text(x+6,y,title,8.6,fontweight='bold')
def snap(key,x,y,w,h,labels=None):
    info=manifest['assets'][key]
    info=info.get('paired_trim',info)
    file=R/info['file'];im=Image.open(file).convert('RGBA')
    a=axis(x,y,w,h);a.imshow(im);a.axis('off')
    for name,loc in (labels or {}).items():
        if name not in info['landmarks']:continue
        a.annotate(name,xy=info['landmarks'][name],xytext=loc,textcoords='axes fraction',
            fontsize=7.1,color=TEAL,ha='center',va='center',
            bbox=dict(fc='white',ec='none',alpha=.88,pad=.5),
            arrowprops=dict(arrowstyle='-',color='#777777',lw=.55,shrinkA=1,shrinkB=1))
    return a

data=pd.read_csv(P/'cam_md_state_contrast_all_common_residues.csv')
blocks=pd.read_csv(P/'cam_md_state_contrast_all_residues_five_blocks.csv')
mets=pd.read_csv(P/'final_cam_all_nine_methionines_with_stability.csv')
groups=pd.read_csv(P/'final_cam_fixed_structural_group_metrics.csv')
manifest=json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
chi=data[data.metric=='hchi'].sort_values('mature_cam_number').copy()
back=data[data.metric=='rmsf'].sort_values('mature_cam_number').copy()
assert len(chi)==122 and len(back)==145
assert back.mature_cam_number.tolist()==list(range(4,149))
assert np.isfinite(chi[['md_delta','full0_delta']].to_numpy()).all()
for metric,g in data.groupby('metric'):
    np.testing.assert_allclose(g.md_b-g.md_a,g.md_delta,atol=1e-12)
    np.testing.assert_allclose(g.full0_b-g.full0_a,g.full0_delta,atol=1e-12)

# a: two actual complete simulation inputs and every side-chain state change.
panel('a','Recognition states and model-supported side-chain redistribution',1,2)
for i,state in enumerate(['NHE1','IQCG']):
    x=1+i*49
    text(x+24,10,state+'–CaM',8.2,ha='center',fontweight='bold')
    snap('whole_'+state,x,15,47,44,{'M72':(.12,.75),'M109':(.88,.30)})
    text(x+24,60,'2YGG-based' if state=='NHE1' else '4LZX-based',7.2,ha='center',color=MD)
text(100,36,'→',12,color=MD,ha='center',va='center')
a=axis(126,18,51,36)
lims=[float(min(chi.md_delta.min(),chi.full0_delta.min())-.08),float(max(chi.md_delta.max(),chi.full0_delta.max())+.08)]
a.plot(lims,lims,color='#B3B8BC',ls='--',lw=.6,zorder=0)
a.axhline(0,color='#D2D6D9',lw=.5);a.axvline(0,color='#D2D6D9',lw=.5)
a.scatter(chi.md_delta,chi.full0_delta,s=11,c=np.where(chi.contact_b>0,TEAL,BLUE),alpha=.7,edgecolors='white',linewidths=.25)
for pos,offset in [(72,(-12,9)),(109,(8,-10)),(144,(8,4))]:
    q=chi[chi.mature_cam_number==pos].iloc[0]
    a.annotate('M'+str(pos),(q.md_delta,q.full0_delta),xytext=offset,textcoords='offset points',fontsize=7,
        arrowprops=dict(arrowstyle='-',lw=.45,color=MD))
a.set(xlim=lims,ylim=lims,xlabel='MD ΔHχ (nats)',ylabel='Predicted ΔHχ (nats)')
a.set_xticks([-1,0,1]);a.set_yticks([-1,0,1]);a.tick_params(pad=2,length=2.5)
text(151.5,10,'All 122 valid side chains',7.8,ha='center')
a.text(.04,.97,'ρ = 0.630',transform=a.transAxes,ha='left',va='top',fontsize=7.5)
text(151.5,64,'MAE: Full 0.171 / Intra 0.234 nats',7.0,ha='center')
text(2,66,'CaM N-/C-lobe',7.0,color=BLUE)
text(31,66,'Partner',7.0,color=WARM)
text(50,66,'Key sites / IQCG contacts',7.0,color=TEAL)

# b: common sequence axis; the complete backbone and side-chain results are kept.
panel('b','A common backbone response accompanies bidirectional side-chain changes',1,75)
text(181,81,'Δ = IQCG − NHE1',7.2,ha='right',color=MD)
lo,hi=3.5,148.5
tracks=axis(19,86,158,3.2);tracks.set_xlim(lo,hi);tracks.set_ylim(0,1);tracks.axis('off')
for l,h,c,name in [(4,75,BLUE,'N-lobe'),(76,81,'#A6AEB4',''),(82,148,PALE,'C-lobe')]:
    tracks.add_patch(Rectangle((l-.5,0),h-l+1,1,fc=c,ec='none',alpha=.5))
    if name:tracks.text((l+h)/2,.5,name,ha='center',va='center',fontsize=7.1,color='#384E61')
for y,metric,height in [(93,'rmsf',17),(115,'hchi',18)]:
    a=axis(19,y,158,height)
    g=data[data.metric==metric].set_index('mature_cam_number').reindex(range(4,149))
    for l,h in [(20,31),(56,67),(93,104),(129,140)]:a.axvspan(l-.5,h+.5,color=TEAL,alpha=.075,lw=0,zorder=0)
    a.axhline(0,color='#9BA2A8',lw=.55,zorder=0)
    a.plot(g.index,g.md_delta,color=MD,lw=1.05,label='MD reference',zorder=2)
    a.plot(g.index,g.full0_delta,color=BLUE,lw=1.05,label='Full prediction',zorder=3)
    a.set_xlim(lo,hi);a.tick_params(axis='x',labelbottom=False,bottom=False)
    a.set_ylabel('ΔRMSF (Å)' if metric=='rmsf' else 'ΔHχ (nats)',labelpad=4)
    a.spines['bottom'].set_visible(False)
    if metric=='rmsf':
        a.set_ylim(min(g.md_delta.min(),g.full0_delta.min())-.12,max(g.md_delta.max(),g.full0_delta.max())+.12)
        text(24,90,'MD',7.0,color=MD);text(35,90,'Full',7.0,color=BLUE)
        a.text(.99,.95,'145 positions · ρ = 0.827',ha='right',va='top',fontsize=7.2,transform=a.transAxes)
    else:
        for pos in [72,109,144]:
            q=chi[chi.mature_cam_number==pos].iloc[0]
            a.scatter([pos],[q.full0_delta],marker='s',s=14,color=TEAL,edgecolors='white',linewidths=.3,zorder=5)
        a.text(.99,.97,'67 increases / 55 decreases in MD',ha='right',va='top',fontsize=7.0,transform=a.transAxes)
    sources.extend([dict(panel='b',metric=metric,position=int(p),md=float(v.md_delta),full=float(v.full0_delta)) for p,v in g.dropna(subset=['md_delta']).iterrows()])
a=axis(19,136,158,5.5);a.set_xlim(lo,hi);a.set_ylim(-.5,1.5)
for j,(col,c) in enumerate([('contact_b',TEAL),('contact_a',WARM)]):
    for p,v in back.set_index('mature_cam_number')[col].items():
        if v:a.add_patch(Rectangle((p-.5,j-.35),1,.70,fc=c,ec='none'))
a.set_yticks([0,1],['IQCG','NHE1']);a.tick_params(axis='y',length=0,pad=4,labelsize=7.2)
a.set_xticks([4,20,40,60,80,100,120,140,148]);a.tick_params(axis='x',length=2,pad=3)
a.set_xlabel('Mature CaM residue number',labelpad=3);a.spines['left'].set_visible(False)
text(19,149,'Pale vertical bands: four EF-hand loops. Contact tracks: interchain Cα distance ≤ 8 Å.',6.8,color=MD)

# c: every residue belongs to exactly one input-defined contact category.
panel('c','Contact groups show opposite responses',1,158)
group_names=['NHE1_contact_only','IQCG_contact_only','contact_in_both_states','outside_both_CA8_interfaces']
group_labels=['NHE1 only','IQCG only','Both','Neither']
masklist=[(chi.contact_a>0)&(chi.contact_b==0),(chi.contact_a==0)&(chi.contact_b>0),
          (chi.contact_a>0)&(chi.contact_b>0),(chi.contact_a==0)&(chi.contact_b==0)]
g=groups[groups.metric=='hchi'].set_index('region').loc[group_names]
assert g.n.tolist()==[8,21,3,90]
a=axis(15,169,72,24);a.axhline(0,color='#B2B7BB',lw=.55)
for j,(name,mask) in enumerate(zip(group_names,masklist)):
    ids=set(chi.loc[mask,'mature_cam_number'])
    q=blocks[(blocks.metric=='hchi')&blocks.mature_cam_number.isin(ids)].groupby('block').md_delta.mean()
    assert len(q)==5
    np.testing.assert_allclose(q.mean(),g.loc[name,'MD_delta'],atol=1e-6)
    a.scatter(j+np.linspace(-.12,.12,5),q,s=10,color='#A2A9AE',edgecolors='none',zorder=2)
    a.plot([j-.17,j+.17],[g.loc[name,'MD_delta'],g.loc[name,'Full_delta']],color='#A6B5C3',lw=.8,zorder=1)
    a.scatter(j-.17,g.loc[name,'MD_delta'],s=22,color=MD,marker='o',edgecolors='white',linewidths=.4,zorder=4)
    a.scatter(j+.17,g.loc[name,'Full_delta'],s=22,color=BLUE,marker='s',edgecolors='white',linewidths=.4,zorder=4)
    sources.extend([dict(panel='c',group=name,block=int(k),MD=float(v)) for k,v in q.items()])
a.set_xticks(range(4),[s+'\nn = '+str(n) for s,n in zip(group_labels,g.n)],fontsize=7.0)
a.set_ylabel('Mean ΔHχ (nats)',labelpad=3);a.set_xlim(-.5,3.5)
a.set_ylim(-.35,.28);a.set_yticks([-.3,0,.3])
text(15,165,'MD mean',7,color=MD);text(37,165,'Full',7,color=BLUE);text(50,165,'Five MD windows',7,color='#8B9399')

# d: all nine Met sites, including inaccurate predictions; temporal windows are not replicates.
panel('d','All nine methionines',97,158)
m=mets[mets.metric=='hchi'].sort_values('mature_cam_number')
assert len(m)==9
a=axis(109,169,68,24);a.axhline(0,color='#B2B7BB',lw=.55)
for j,(_,q) in enumerate(m.iterrows()):
    a.plot([j,j],[q.min_block_delta,q.max_block_delta],color='#BBC0C4',lw=1.15,zorder=1)
    a.plot([j-.12,j+.12],[q.md_delta,q.full0_delta],color='#A6B5C3',lw=.75,zorder=2)
    a.scatter(j-.12,q.md_delta,s=17,color=MD,edgecolors='white',linewidths=.3,zorder=4)
    a.scatter(j+.12,q.full0_delta,s=17,color=BLUE,marker='s',edgecolors='white',linewidths=.3,zorder=4)
    sources.append(dict(panel='d',position=int(q.mature_cam_number),MD=float(q.md_delta),Full=float(q.full0_delta),window_min=float(q.min_block_delta),window_max=float(q.max_block_delta)))
a.set_xticks(range(9),[str(int(v)) for v in m.mature_cam_number],fontsize=7.0)
a.set_xlabel('Mature CaM Met position',labelpad=3);a.set_ylabel('ΔHχ (nats)',labelpad=3)
a.set_ylim(min(m.min_block_delta.min(),m.full0_delta.min())-.07,max(m.max_block_delta.max(),m.full0_delta.max())+.08)
a.set_xlim(-.55,8.55);a.set_yticks([-.4,0,.4])
text(109,165,'MD',7,color=MD);text(121,165,'Full',7,color=BLUE);text(134,165,'MD window range',7,color='#8B9399')

# e: actual local structures sit beside both dynamical changes at the same site.
panel('e','Local environments connect the two readouts at M72 and M109',1,207)
for x,pos in [(1,72),(94,109)]:
    for j,state in enumerate(['NHE1','IQCG']):
        key='detail_focus_M'+str(pos)+'_'+state
        if key not in manifest['assets']:key='detail_M'+str(pos)+'_'+state
        if key not in manifest['assets']:key='M'+str(pos)+'_'+state
        snap(key,x+j*22.5,216,23,25,{'M'+str(pos):(.5,.10)})
        text(x+11.5+j*22.5,213,state,7.0,ha='center',color=MD)
    q=chi[chi.mature_cam_number==pos].iloc[0];b=back[back.mature_cam_number==pos].iloc[0]
    text(x+48,214,'M'+str(pos)+(' · neither' if pos==72 else ' · IQCG contact'),7.3,fontweight='bold',color=TEAL)
    text(x+48,219,'ΔHχ    MD / Full',7.0)
    text(x+48,223,f'{q.md_delta:+.3f} / {q.full0_delta:+.3f} nats',7.1,color=BLUE)
    text(x+48,228,'ΔRMSF  '+f'{b.md_delta:+.3f} / {b.full0_delta:+.3f} Å',7.0,color=MD)
text(181,242,'Fixed fold 0, training-set illustration · One 500 ns trajectory per state; five consecutive 100 ns windows.',6.7,ha='right',color=MD)

fig.canvas.draw();renderer=fig.canvas.get_renderer()
visible=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[]
for t in visible:
    box=t.get_window_extent(renderer)
    if box.x0<-.6 or box.y0<-.6 or box.x1>fig.bbox.width+.6 or box.y1>fig.bbox.height+.6:outside.append({'text':t.get_text(),'bbox':[box.x0,box.y0,box.x1,box.y1]})
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:
    fig.savefig(OUT/('fig3_CaM_integrated.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
sd=OUT/'source_data';sd.mkdir(exist_ok=True)
for file in ['cam_md_state_contrast_all_common_residues.csv','cam_md_state_contrast_all_residues_five_blocks.csv',
    'final_cam_fixed_structural_group_metrics.csv','final_cam_all_nine_methionines_with_stability.csv',
    'cam_md_per_lobe_fitted_rmsf_control_summary.csv','cam_md_mature_number_mapping.csv']:
    if (P/file).resolve()!=(sd/file).resolve():shutil.copy2(P/file,sd/file)
pd.DataFrame(sources).to_csv(sd/'all_plotted_quantities.csv',index=False)
qa=dict(canvas_mm=[W,H],dpi=450,data_positions=dict(backbone=145,sidechain=122,Met=9,contact_groups=[8,21,3,90]),
    outside_text=outside,min_font_pt=min(t.get_fontsize() for t in visible),MD_window_definition='five consecutive100ns blocks of one500ns trajectory/state, not independent replicates',
    palette=dict(Full=BLUE,Intra=WARM,additional=TEAL),visual_review='pending',
    prediction_scope='same frozen fold0 model, both structures present in its training data; interpretive illustration',
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sd.iterdir() if p.is_file()},
    output_sha256={ext:hashlib.sha256((OUT/('fig3_CaM_integrated.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:qa[k] for k in ['canvas_mm','data_positions','outside_text','min_font_pt']},ensure_ascii=False))
plt.close(fig)
