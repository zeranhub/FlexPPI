# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) real PyMOL inputs + 244 absolute-state points -> existing real-structure renderer and correlation style -> param inherit
# (b) shared-sequence profiles + contact tracks -> LineTrend/plot_trend.py and heatmap/plot_comparison.py -> param inherit
# (c) group means + all5 continuous windows -> GroupedBarChart/plot_GroupedBarChartv1.py -> param inherit
# (d) all9 Met absolute states + all5 windows -> same grouped-point visual logic -> param inherit
# (e) actual complete-chain PyMOL local contexts + paired absolute readouts -> image_plate + quant -> param inherit
# Semantic mismatch: timeline and method-comparison heatmap assets do not accept residue tracks.
# Their visual logic is inherited; no template data or downsampling is used.
# User requested manuscript palette overrides after mandatory baselines below.
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

# --- panel code ---
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image

OUT=Path(__file__).resolve().parent; P=OUT/'source_data'
R=OUT/'structure_assets'
W,H=183.,245.; BLUE='#2166AC'; WARM='#B87936'; TEAL='#3A8D8C'; MD='#626A70'
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

panel('a','Recognition states and complete side-chain prediction',1,2)
for i,(state,_,_,_) in enumerate(STATES):
    x=1+i*49;text(x+24,10,state+'–CaM',8.2,ha='center',fontweight='bold')
    snap('whole_'+state,x,14,47,37)
    text(x+24,52,'2YGG-based' if i==0 else '4LZX-based',7,ha='center',color=MD)
text(100,33,'→',12,color=MD,ha='center',va='center')
a=axis(126,17,51,31);a.plot([0,4.1],[0,4.1],color='#B3B8BC',ls='--',lw=.6)
for state,s,c,marker in STATES:
    a.scatter(chi['md_'+s],chi['full0_'+s],s=10,marker=marker,color=c,alpha=.58,edgecolors='white',linewidths=.23)
    points.extend(dict(panel='a',position=int(q.mature_cam_number),state=state,MD=q['md_'+s],Full=q['full0_'+s]) for _,q in chi.iterrows())
a.set(xlim=(0,4.1),ylim=(0,4.1),xlabel='MD Hχ (nats)',ylabel='Predicted Hχ (nats)');a.set_xticks([0,2,4]);a.set_yticks([0,2,4]);a.tick_params(pad=2,length=2)
text(151.5,10,'All 122 side chains in each state',7.4,ha='center')
text(2,57,'Structure: CaM N-/C-lobe',6.8,color=BLUE);text(47,57,'Partner',6.8,color=WARM);text(64,57,'Key sites',6.8,color=TEAL)
text(151.5,57,'Two states shown separately',6.8,ha='center',color=MD)

panel('b','Absolute backbone and side-chain profiles',1,65)
text(127,65,'● NHE1',7.5,color=WARM);text(153,65,'■ IQCG',7.5,color=BLUE)
lo,hi=3.5,148.5;a=axis(20,74,157,3);a.set(xlim=(lo,hi),ylim=(0,1));a.axis('off')
for l,h,c,name in [(4,75,BLUE,'N-lobe'),(76,81,'#A6AEB4',''),(82,148,'#85ABD1','C-lobe')]:
    a.add_patch(Rectangle((l-.5,0),h-l+1,1,fc=c,ec='none',alpha=.5))
    if name:a.text((l+h)/2,.5,name,ha='center',va='center',fontsize=7,color='#384E61')
for metric,kind,y in [('rmsf','md',81),('rmsf','full0',97),('hchi','md',116),('hchi','full0',132)]:
    a=axis(20,y,157,12);g=data[data.metric==metric].set_index('mature_cam_number').reindex(range(4,149))
    text(21,y-3,'MD reference' if kind=='md' else 'Full prediction',6.9,color=MD)
    for l,h in [(20,31),(56,67),(93,104),(129,140)]:a.axvspan(l-.5,h+.5,color=TEAL,alpha=.065,lw=0)
    for state,s,c,_ in STATES:
        a.plot(g.index,g[kind+'_'+s],color=c,lw=.9)
        points.extend(dict(panel='b',metric=metric,kind=kind,position=int(p),state=state,value=float(v)) for p,v in g[kind+'_'+s].dropna().items())
    a.set(xlim=(lo,hi),ylim=(0,4.1 if metric=='hchi' else 4),ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)')
    a.set_yticks([0,2,4]);a.tick_params(axis='x',labelbottom=False,bottom=False);a.tick_params(axis='y',pad=2,length=2);a.spines['bottom'].set_visible(False)
    check(a,[g[kind+'_a'],g[kind+'_b']])
a=axis(20,148,157,5);a.set(xlim=(lo,hi),ylim=(-.5,1.5))
for j,(col,c) in enumerate([('contact_b',BLUE),('contact_a',WARM)]):
    for p,v in back.set_index('mature_cam_number')[col].items():
        if v:a.add_patch(Rectangle((p-.5,j-.35),1,.70,fc=c,ec='none'))
a.set_yticks([0,1],['IQCG','NHE1']);a.tick_params(axis='y',length=0,pad=4,labelsize=6.8)
a.set_xticks([4,20,40,60,80,100,120,140,148]);a.tick_params(axis='x',length=2,pad=3);a.set_xlabel('Mature CaM residue number',labelpad=2);a.spines['left'].set_visible(False)
text(20,161,'145 backbone / 122 side-chain positions · Pale bands: EF-hand loops · Tracks: interchain Cα ≤ 8 Å',6.5,color=MD)

panel('c','Complete contact groups',1,169);panel('d','All nine methionines',96,169)
masks=[(chi.contact_a>0)&(chi.contact_b==0),(chi.contact_a==0)&(chi.contact_b>0),(chi.contact_a>0)&(chi.contact_b>0),(chi.contact_a==0)&(chi.contact_b==0)]
labels=['NHE1\nonly\nn=8','IQCG\nonly\nn=21','Both\n\nn=3','Neither\n\nn=90']
for kind,x in [('md',15),('full0',56)]:
    text(x+15.5,176,'MD reference' if kind=='md' else 'Full prediction',7.1,ha='center',color=MD)
    a=axis(x,181,31,20);ys=[]
    for j,mask in enumerate(masks):
        g=chi.loc[mask];assert len(g)==[8,21,3,90][j];vals=[]
        for offset,(state,s,c,marker) in zip([-.15,.15],STATES):
            v=float(g[kind+'_'+s].mean());vals.append(v)
            a.scatter(j+offset,v,s=18,marker=marker,color=c,zorder=3)
            if kind=='md':
                b=blocks[(blocks.metric=='hchi')&blocks.mature_cam_number.isin(g.mature_cam_number)].groupby('block')['md_'+s].mean()
                a.vlines(j+offset,b.min(),b.max(),color=c,alpha=.4,lw=.8)
            points.append(dict(panel='c',group=j,state=state,kind=kind,value=v))
        a.plot([j-.15,j+.15],vals,color='#B6BEC5',lw=.7,zorder=1);ys.append(vals)
    a.set(xlim=(-.4,3.4),ylim=(1.3,2.2));a.set_yticks([1.4,1.8,2.2]);a.set_xticks(range(4),labels,fontsize=6.2);a.tick_params(pad=2,length=2)
    if kind=='md':a.set_ylabel('Mean Hχ (nats)',labelpad=2)
    else:a.tick_params(labelleft=False)
    check(a,ys)
for kind,x in [('md',108),('full0',149)]:
    text(x+14.5,176,'MD reference' if kind=='md' else 'Full prediction',7.1,ha='center',color=MD)
    a=axis(x,181,29,20);ys=[]
    for state,s,c,marker in STATES:
        values=mets[kind+'_'+s];ys.append(values);a.plot(range(9),values,color=c,marker=marker,ms=2.2,lw=.8)
        for j,(_,q) in enumerate(mets.iterrows()):
            if kind=='md':
                b=blocks[(blocks.metric=='hchi')&(blocks.mature_cam_number==q.mature_cam_number)]['md_'+s]
                a.vlines(j,b.min(),b.max(),color=c,alpha=.35,lw=.8)
            points.append(dict(panel='d',position=int(q.mature_cam_number),state=state,kind=kind,value=float(q[kind+'_'+s])))
    a.set(xlim=(-.4,8.4),ylim=(1.0,3.4));a.set_yticks([1,2,3]);a.set_xticks(range(9),[str(int(v)) for v in mets.mature_cam_number],fontsize=6.1,rotation=60);a.tick_params(pad=2,length=2)
    if kind=='md':a.set_ylabel('Hχ (nats)',labelpad=2)
    else:a.tick_params(labelleft=False)
    check(a,ys)
text(106,211,'Mature CaM Met position',6.7,color=MD)

panel('e','Local environments and both absolute readouts',1,217)
for x,pos in [(1,72),(94,109)]:
    for j,(state,_,c,_) in enumerate(STATES):
        snap('detail_focus_M'+str(pos)+'_'+state,x+j*22.5,224,23,17);text(x+11.5+j*22.5,222,state,6.6,ha='center',color=c)
    q=chi[chi.mature_cam_number==pos].iloc[0];b=back[back.mature_cam_number==pos].iloc[0]
    text(x+48,222,'M'+str(pos)+(' · neither' if pos==72 else ' · IQCG contact'),7,fontweight='bold',color=TEAL)
    for y,s in [(225.5,'                 NHE1    IQCG'),(229,f'Hχ  MD     {q.md_a:.3f}    {q.md_b:.3f}'),(232,f'      Full    {q.full0_a:.3f}    {q.full0_b:.3f}'),(235,f'RMSF MD  {b.md_a:.3f}    {b.md_b:.3f}'),(238,f'      Full    {b.full0_a:.3f}    {b.full0_b:.3f}')]:text(x+48,y,s,6.4)
text(181,242,'Fixed fold 0, training-set illustration · Hχ: nats; RMSF: Å · MD ranges: five consecutive 100 ns windows.',6.3,ha='right',color=MD)
fig.canvas.draw();renderer=fig.canvas.get_renderer();visible=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in visible if t.get_window_extent(renderer).x0<-.6 or t.get_window_extent(renderer).y0<-.6 or t.get_window_extent(renderer).x1>fig.bbox.width+.6 or t.get_window_extent(renderer).y1>fig.bbox.height+.6]
assert not outside,outside
assert not any('Δ' in a.get_ylabel() for a in axes)
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig3_CaM_absolute_profiles.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
pd.DataFrame(points).to_csv(P/'all_plotted_absolute_quantities.csv',index=False)
qa=dict(canvas_mm=[W,H],all_axes_absolute=True,two_state_profiles=True,all_residue_counts=dict(backbone=145,sidechain=122,Met=9,contact_groups=[8,21,3,90]),
    plotted_rows=pd.DataFrame(points).groupby('panel').size().to_dict(),outside_text=outside,min_font_pt=min(t.get_fontsize() for t in visible),axis_checks=checks,visual_review='pending',
    output_sha256={ext:hashlib.sha256((OUT/('fig3_CaM_absolute_profiles.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False));plt.close(fig)
