# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Original3063 ECDF -> original_fig4 production + LineTrend/plot_trend.py -> param inherit
# (b) All83 paired-rho differences -> PairedBoxScatter/plot_PairedBoxScatter1.R -> param inherit
# (c) All83 selection contrasts -> PairedBoxScatter/plot_PairedBoxScatter1.R -> param inherit
# (d) Original66 profiles -> original_fig7/plot_english.py + LineTrend/plot_trend.py -> param inherit
# (e) Original151 angle observations -> original_fig7/plot_english.py -> param inherit
# (f) Original432 curve bins -> original_fig7/plot_english.py + LineTrend/plot_trend.py -> param inherit
# RULE: all native vector plots; no panel image embedding, point sampling, smoothing, or bin changes.
import sys
sys.path.insert(0,r'E:/Project/FlexPPI/experiments/_sources/completed_results_20261012/plot_runtime')
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


import json, hashlib, re, xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from pathlib import Path

OUT=Path(__file__).resolve().parent
DATA=OUT/'source'
STEM='fig5_joint_readouts_conformations'
BLUE=CATEGORICAL[0]; ORANGE=CATEGORICAL[3]; PURPLE=CATEGORICAL[4]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:pd.read_csv(p,float_precision='round_trip')
font_manager.findfont('Arial',fallback_to_default=False)
# The 183 mm source is inserted at 167.64 mm in Word; 9 pt -> 8.24 pt.
mpl.rcParams.update({'font.sans-serif':['Arial','Liberation Sans'], 'font.size':9,
    'axes.labelsize':9,'axes.titlesize':9,'xtick.labelsize':9,'ytick.labelsize':9,
    'legend.fontsize':9,'axes.unicode_minus':False,'savefig.bbox':None})

def filtered(df):
    df=df.copy()
    if 'threshold' in df: df=df[np.isclose(pd.to_numeric(df.threshold),1)]
    if 'scope' in df: df=df[df.scope=='all']
    if 'metric' in df: df=df[df.metric=='fine_entropy_mean_per_chi']
    return df

ecdf=read(DATA/'original_fig4/fig4_joint_complementarity_all_ecdf_points.csv')
assert set(ecdf.role)=={'reference','prediction'}
assert ecdf.groupby('role').size().to_dict()=={'prediction':3063,'reference':3063}
width=filtered(read(DATA/'width/replica1_per_record_correlations.csv'))
entropy=filtered(read(DATA/'entropy/replica1_per_record_correlations.csv'))
selection=filtered(read(DATA/'width/replica1_high_low_per_record.csv'))
width_ci=filtered(read(DATA/'width/replica1_paired_correlation_difference_summary.csv'))
entropy_ci=filtered(read(DATA/'entropy/replica1_paired_correlation_difference_summary.csv'))
selection_ci=filtered(read(DATA/'width/replica1_high_low_summary.csv'))
all_ids=None; paired=[]; pair_stats=[]
for label,df,ci in [('Fine-bin entropy',entropy,entropy_ci),('Low-F* width',width,width_ci)]:
    wide=df.pivot(index='sample_id',columns='predictor',values='rho_type_centered')
    assert wide.shape==(83,3) and np.isfinite(wide.to_numpy()).all()
    if all_ids is None: all_ids=set(wide.index)
    assert set(wide.index)==all_ids
    assert int(df[df.predictor=='full_hchi'].positions.sum())==33142
    meta=df[df.predictor=='full_hchi'].set_index('sample_id')
    delta=wide.full_hchi-wide.full_rmsf
    q=ci[(ci.dataset=='all')&(ci.contrast=='full_hchi_minus_full_rmsf')].iloc[0]
    np.testing.assert_allclose(delta.mean(),q['mean'],rtol=0,atol=1e-14)
    for sid in wide.index:
        paired.append({'metric':label,'sample_id':sid,'dataset':meta.loc[sid,'dataset'],
            'positions':int(meta.loc[sid,'positions']),'full_Hchi_rho':wide.loc[sid,'full_hchi'],
            'full_RMSF_rho':wide.loc[sid,'full_rmsf'],'intrachain_Hchi_rho':wide.loc[sid,'intrachain_hchi'],
            'paired_rho_difference':delta.loc[sid]})
    pair_stats.append({'metric':label,'n':83,'mean_difference':float(q['mean']),
        'CI_low':float(q.bootstrap_95_low),'CI_high':float(q.bootstrap_95_high),
        'positive_records':int((delta>0).sum()),'negative_records':int((delta<0).sum()),
        'mean_full_Hchi':float(wide.full_hchi.mean()),'mean_full_RMSF':float(wide.full_rmsf.mean()),
        'mean_intra_Hchi':float(wide.intrachain_hchi.mean())})
paired=pd.DataFrame(paired)
assert len(selection)==249 and set(selection.sample_id)==all_ids
assert selection.groupby('predictor').size().to_dict()=={'full_hchi':83,'full_rmsf':83,'intrachain_hchi':83}
selection_ci=selection_ci[selection_ci.dataset=='all'].set_index('predictor')
for predictor,group in selection.groupby('predictor'):
    np.testing.assert_allclose(group.high_minus_low_MD.mean(),selection_ci.loc[predictor,'mean'],rtol=0,atol=1e-13)

UB=DATA/'original_fig7'; rows=read(UB/'source/all66_residue_summary.csv')
z=np.load(UB/'source/ubiquitin_all66_chi_501frames.npz',allow_pickle=False)
hist=pd.read_csv(UB/'source/histograms_501frames.csv')
defined=read(UB/'experiment/defined_angles_long.csv')
corr=read(UB/'experiment/per_chi_correlations.csv')
corr=corr[(corr.predictor=='predicted_joint_Hchi')&(corr.metric=='fine_entropy')].copy()
assert len(rows)==66 and len(defined)==151 and len(hist)==432 and len(corr)==4
np.testing.assert_array_equal(rows.predicted_Hchi_nats,z['predicted_Hchi_nats'])
np.testing.assert_array_equal(rows.frozen_MD_joint_Hchi_nats,z['frozen_MD_joint_Hchi_nats'])
rows['predicted_Hchi_nats']=z['predicted_Hchi_nats']; rows['frozen_MD_joint_Hchi_nats']=z['frozen_MD_joint_Hchi_nats']
expected_n={1:66,2:52,3:23,4:10}

# Revised 2026-10-25. Frozen numerical source inherited from original_production.py.
# Layout follows the approved 6-panel organization. No data rows are omitted.
from PIL import Image
BLUE='#2166AC'; WARM='#B87936'; TEAL='#3A8D8C'; ORANGE=WARM; PURPLE=TEAL
BLACK='#222222'; GREY='#999999'
W,H=183.,230.
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
axis=lambda x,y,w,h,**kw:fig.add_axes([x/W,y/H,w/W,h/H],**kw)
text=lambda x,y,s,**kw:fig.text(x/W,y/H,s,**kw)
panel_boxes={}; axes=[]; points=[]
def panel(letter,x,y,title):
    text(x,y,letter,fontsize=12,weight='bold',va='baseline')
    text(x+7,y,title,fontsize=9.5,weight='bold',va='baseline')
def register(name,ax):
    axes.append(ax);panel_boxes[name]=list(ax.get_position().bounds)
def allpoints(panel_name,p,xcol,ycol,kind,**extra):
    for _,v in p.iterrows():
        points.append(dict(panel=panel_name,sample_id=v.get('sample_id',''),kind=kind,x=float(v[xcol]),y=float(v[ycol]),**extra))
# a-c: all records and original statistics.
panel('a',2,223,'Distinct readouts')
text(15,216,'3,063 records',fontsize=9,color=BLACK)
ax=axis(15,182,37,28);register('a',ax)
for role,color,ls,name in [('reference',GREY,'--','MD'),('prediction',BLUE,'-','Model')]:
    p=ecdf[ecdf.role==role].sort_values('y')
    ax.plot(p.x,p.y,color=color,lw=1.15,ls=ls)
    for _,v in p.iterrows():points.append({'panel':'a','sample_id':v.sample_id,'kind':role,'x':v.x,'y':v.y})
    ax.text(.035,.90 if role=='reference' else .70,name+' '+format(float(p.x.mean()),'.3f'),transform=ax.transAxes,color=color,fontsize=9)
ax.set_xlim(-.25,1);ax.set_ylim(0,1.02);ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
ax.set_xlabel('RMSF–Hχ Spearman ρ',labelpad=3);ax.set_ylabel('Fraction of records',labelpad=3)

panel('b',59,223,'Occupancy association')
text(81,216,'83 paired records',fontsize=9)
ax=axis(82,182,36,28);register('b',ax)
for i,st in enumerate(pair_stats):
    y=1-i;p=paired[paired.metric==st['metric']].sort_values('sample_id').reset_index(drop=True)
    jitter=(np.arange(len(p))*37%83/82-.5)*.26
    ax.scatter(p.paired_rho_difference,y+jitter,s=8.5,color=GREY,alpha=.47,linewidths=0,zorder=2)
    ax.errorbar(st['mean_difference'],y,xerr=np.array([[st['mean_difference']-st['CI_low']],[st['CI_high']-st['mean_difference']]]),
        fmt='D',ms=4,color=BLUE,ecolor=BLUE,elinewidth=1.4,capsize=2,zorder=4)
    ax.text(.98,y+.29,format(st['mean_difference'],'.3f'),transform=ax.get_yaxis_transform(),ha='right',color=BLUE,fontsize=9)
    for j,v in p.iterrows():points.append({'panel':'b','sample_id':v.sample_id,'kind':v.metric,'x':v.paired_rho_difference,'y':float(y+jitter[j])})
ax.axvline(0,color=GREY,lw=.6,ls=':',zorder=0)
ax.set_yticks([1,0],['Bin entropy','Low-F* width']);ax.set_ylim(-.5,1.6);ax.set_xlim(-.25,.55);ax.set_xticks([0,.25,.5])
ax.set_xlabel('Δρ: Hχ − RMSF',labelpad=3);ax.tick_params(axis='y',length=0,pad=3)

panel('c',124,223,'Prediction-selected groups')
text(142,216,'83 records',fontsize=9)
ax=axis(142,182,35,28);register('c',ax)
for i,(pred,label,color,marker) in enumerate([('full_hchi','Hχ',BLUE,'D'),('intrachain_hchi','Intra Hχ',WARM,'o'),('full_rmsf','RMSF',BLACK,'s')]):
    y=2-i;p=selection[selection.predictor==pred].sort_values('sample_id').reset_index(drop=True)
    jitter=(np.arange(len(p))*37%83/82-.5)*.26
    ax.scatter(p.high_minus_low_MD,y+jitter,s=8.5,color=GREY,alpha=.43,linewidths=0,zorder=2)
    st=selection_ci.loc[pred];mean=float(st['mean'])
    ax.errorbar(mean,y,xerr=np.array([[mean-st.bootstrap_95_low],[st.bootstrap_95_high-mean]]),
        fmt=marker,ms=4,color=color,mfc=color,ecolor=color,elinewidth=1.3,capsize=2,zorder=4)
    ax.text(.98,y+.29,format(mean,'.2f'),transform=ax.get_yaxis_transform(),ha='right',fontsize=9,color=color)
    for j,v in p.iterrows():points.append({'panel':'c','sample_id':v.sample_id,'kind':pred,'x':v.high_minus_low_MD,'y':float(y+jitter[j])})
ax.axvline(0,color=GREY,lw=.6,ls=':',zorder=0)
ax.set_yticks([2,1,0],['Hχ','Intra Hχ','RMSF']);ax.set_ylim(-.5,2.6);ax.set_xlim(-1.5,17.5);ax.set_xticks([0,8,16])
ax.set_xlabel('MD width difference (°/χ)',labelpad=3);ax.tick_params(axis='y',length=0,pad=3)

# d: complete input structure locator beside all 66 unchanged profile values.
panel('d',2,161,'Ubiquitin: structural positions and complete side-chain profile')
text(62,154,'66 positions  |  MAE 0.080 nats  |  ρ 0.979',fontsize=9)
structure_path=OUT.parent/'pymol/ubiquitin_Q41_Q2_locator_trim.png'
structure_meta=OUT.parent/'pymol/ubiquitin_Q41_Q2_locator.json'
assert structure_path.is_file(),structure_path
im=Image.open(structure_path)
im_ax=axis(3,111,42,42);im_ax.imshow(im);im_ax.axis('off');register('d_structure',im_ax)
# Labels are vector text and target the exact projected residue landmarks.
meta=json.loads(structure_meta.read_text(encoding='utf-8-sig'))
landmarks=meta['asset']['paired_trim']['landmarks']
def xy(name):
    obj=landmarks[name]
    if isinstance(obj,dict):return float(obj['x']),float(obj['y'])
    return float(obj[0]),float(obj[1])
for site,color in [('Q2',TEAL),('Q41',WARM)]:
    xx,yy=xy(site)
    anchor=(.04,.50) if site=='Q2' else (.77,.97)
    im_ax.annotate(site,xy=(xx,yy),xycoords='data',xytext=anchor,textcoords='axes fraction',
        fontsize=9,color=color,ha='left',va='center',weight='bold',
        bbox=dict(facecolor='white',edgecolor='none',pad=.8,alpha=.9), arrowprops=dict(arrowstyle='-',color=color,lw=.7,shrinkA=2,shrinkB=1),annotation_clip=False)
text(24,106,'2QHO input',ha='center',fontsize=9,color=BLACK)
text(24,101,'Ub + UBA',ha='center',fontsize=9,color=GREY)
ax=axis(61,112,116,36);register('d_profile',ax)
ids=np.arange(1,73);model=np.full(72,np.nan);reference=model.copy()
for _,r in rows.iterrows():
    pos=int(r.author_residue_id)-1;model[pos]=r.predicted_Hchi_nats;reference[pos]=r.frozen_MD_joint_Hchi_nats
    for name,value in [('prediction',model[pos]),('MD',reference[pos])]:points.append({'panel':'d','residue_id':pos+1,'kind':name,'x':pos+1,'y':value})
for site,color in [(41,WARM),(2,TEAL)]:ax.axvspan(site-.55,site+.55,color=color,alpha=.18,lw=0,zorder=0)
ax.plot(ids,reference,color=BLACK,lw=.8,ls=':',marker='o',ms=2.5,mfc='white',mew=.6,zorder=2)
ax.plot(ids,model,color=BLUE,lw=1.,marker='s',ms=2.3,mew=0,zorder=3)
ax.set_xlim(.5,72.5);ax.set_ylim(-.15,4.25);ax.set_yticks([0,2,4]);ax.set_ylabel('Joint Hχ (nats)',labelpad=3)
ax.set_xticks([1,10,20,30,40,50,60,72])
for site,color,marker in [(2,TEAL,'s'),(41,WARM,'D')]:
    ax.scatter([site],[model[site-1]],color=color,marker=marker,s=27,zorder=5)
    ax.text(site,4.1,'Q'+str(site),ha='center',va='top',fontsize=9,color=color)
handles=[Line2D([0],[0],color=BLACK,ls=':',lw=.8,marker='o',ms=3,mfc='white',mew=.7),Line2D([0],[0],color=BLUE,lw=1,marker='s',ms=3,mew=0)]
ax.legend(handles,['MD','Model'],loc='upper right',ncol=2,fontsize=9,handlelength=1.5,handletextpad=.4,columnspacing=1.,borderaxespad=.3)
ax.set_xlabel('Ubiquitin residue number',labelpad=3)

# e: four equal axes, all 151 original residue/angle observations.
panel('e',2,92,'One joint output corresponds to occupancy of individual side-chain angles')
for chi,left in [(1,20),(2,63),(3,106),(4,149)]:
    ax=axis(left,55,28,26);register('e_chi'+str(chi),ax)
    p=defined[defined.chi=='chi'+str(chi)];assert len(p)==expected_n[chi]
    r=corr[corr.chi=='chi'+str(chi)].iloc[0]
    ax.set_title('χ'+str(chi)+'  n='+str(len(p))+'\nρ = '+format(float(r.rho_AA_centered),'.2f'),loc='left',fontsize=9,pad=4)
    for label,g in [('other',p[~p.author_residue_id.isin([2,41])]),('Q2',p[p.author_residue_id==2]),('Q41',p[p.author_residue_id==41])]:
        color=TEAL if label=='Q2' else WARM if label=='Q41' else BLUE
        marker='s' if label=='Q2' else 'D' if label=='Q41' else 'o'
        ax.scatter(g.predicted_Hchi_nats_AA_centered,g.fine_entropy_four_origins_mean_nats_AA_centered,
            s=13 if label=='other' else 25,color=color,marker=marker,edgecolors='white' if label=='other' else color,
            linewidths=.25 if label=='other' else .5,alpha=.75 if label=='other' else 1,zorder=3 if label=='other' else 4)
    for _,v in p.iterrows():points.append({'panel':'e','residue_id':int(v.author_residue_id),'chi':chi,'kind':'fine_entropy_pair','x':float(v.predicted_Hchi_nats_AA_centered),'y':float(v.fine_entropy_four_origins_mean_nats_AA_centered)})
    ax.axhline(0,color=GREY,lw=.45,ls=':',zorder=0);ax.axvline(0,color=GREY,lw=.45,ls=':',zorder=0)
    ax.set_xlim(-.70,.50);ax.set_ylim(-1.,.75);ax.set_xticks([-.5,0,.5]);ax.set_yticks([-.5,0,.5]);ax.tick_params(pad=2,length=2)
    if chi!=1:ax.tick_params(labelleft=False)
text(5,68,'MD entropy deviation (nats)',fontsize=9,rotation=90,ha='center',va='center')
text(99,47,'Predicted joint Hχ deviation (nats; both axes amino-acid centered)',fontsize=9,ha='center')

# f: all 432 original bins and original NaN gaps; these are MD-derived curves.
panel('f',2,37,'MD angular occupancy at prediction-selected Gln positions')
text(20,31,'Q41: lowest prediction (2.000 nats)',fontsize=9,color=WARM)
text(100,31,'Q2: highest prediction (2.990 nats)',fontsize=9,color=TEAL)
maxf=float(hist.F_star.max());widths=[]
for chi,left in enumerate([20,75,130],1):
    ax=axis(left,8,43,21);register('f_chi'+str(chi),ax)
    ax.text(.98,.93,'χ'+str(chi),transform=ax.transAxes,ha='right',va='top',fontsize=9,weight='bold')
    for j,(site,color,style) in enumerate([('Q41',WARM,'-'),('Q2',TEAL,'--')]):
        p=hist[(hist.residue==site)&(hist.chi=='chi'+str(chi))].sort_values('angle_center_degrees');assert len(p)==72
        ax.plot(p.angle_center_degrees,p.F_star,color=color,lw=1.2,ls=style,zorder=3)
        good=np.isfinite(p.F_star)&(p.F_star<=1)
        ax.broken_barh([(float(v),5.) for v in p.loc[good,'angle_lower_degrees']],(-.65+j*.3,.16),facecolors=color,edgecolors='none',zorder=3)
        width=5*int(good.sum());widths.append({'residue':site,'chi':chi,'origin0_lowF_width_degrees':width})
        ax.text(.035 if j==0 else .43,.94,site+' '+str(width)+'°',transform=ax.transAxes,color=color,fontsize=9,va='top')
        for _,v in p.iterrows():points.append({'panel':'f','residue_id':int(site[1:]),'chi':chi,'kind':'MD_Fstar','x':v.angle_center_degrees,'y':v.F_star,'low_Fstar_le1':bool(np.isfinite(v.F_star) and v.F_star<=1)})
    ax.axhline(1,color=GREY,lw=.7,ls=':',zorder=0)
    ax.set_xlim(-180,180);ax.set_ylim(-.8,maxf+1.8);ax.set_yticks([0,2,4]);ax.set_xticks([-180,-60,60,180]);ax.set_xlabel('Side-chain angle (°)',fontsize=9,labelpad=1)
    if chi==1:ax.set_ylabel('MD F*',fontsize=9,labelpad=3)
    else:ax.tick_params(labelleft=False)
    ax.tick_params(pad=2,length=2)

fig.canvas.draw();renderer=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[]
for t in texts:
    b=t.get_window_extent(renderer)
    if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
assert min(t.get_fontsize() for t in texts)>=9
assert not [t.get_text() for t in texts if re.search('[\u3400-\u9fff]',t.get_text())]
assert np.array_equal(np.isnan(model),np.isnan(reference)) and np.isnan(model).sum()==6
assert len(points)==6126+166+249+132+151+432
counts=pd.DataFrame(points).groupby('panel').size().to_dict()
# A fixed physical canvas is essential; do not let bbox='tight' alter the page size.
for ext in ['pdf','svg','png']:fig.savefig(OUT/(STEM+'.'+ext),dpi=450,bbox_inches=fig.bbox_inches)
svg_images=ET.parse(OUT/(STEM+'.svg')).findall('.//{http://www.w3.org/2000/svg}image')
assert len(svg_images)==1, 'Only the real PyMOL locator may be raster.'
pd.DataFrame(points).to_csv(OUT/'all_plotted_points.csv',index=False,float_format='%.17g')
paired.to_csv(OUT/'panel_b_all83_paired_rho_differences.csv',index=False,float_format='%.17g')
selection.to_csv(OUT/'panel_c_all83_selection_differences.csv',index=False,float_format='%.17g')
pd.DataFrame(pair_stats).to_csv(OUT/'panel_b_statistics.csv',index=False,float_format='%.17g')
pd.DataFrame(widths).to_csv(OUT/'panel_f_unchanged_display_widths.csv',index=False)
defined.to_csv(OUT/'panel_e_all151_scatter_source.csv',index=False,float_format='%.17g')
(OUT/'visible_text.json').write_text(json.dumps([t.get_text() for t in texts],ensure_ascii=False,indent=2),encoding='utf8')
caption='图5：联合动态读数与MD侧链构象占据的对应。（a）3,063条记录中RMSF与Hχ的同类型相关分布，灰/蓝为MD/预测，标注为记录均值。（b）83条记录、33,142个同类型有效位置中，Hχ与RMSF预测对MD细箱熵、低自由能角区总宽度的相关之差，每点为一条记录。（c）仅按预测划分同类残基高低半组，展示全部83条记录的真实角区总宽度组差；正值表示高预测组更宽。b、c大点及误差线为均值及95%来源分层重采样区间；Intra为链内Hχ，其余为完整模型。（d）PyMOL绘制的泛素–UBA输入结构定位及全部66个Hχ有效位置的预测与MD参考；泛素为蓝色，UBA为灰色，Q41暖赭、Q2青绿，Gly/Ala位置曲线断线。此结构仅用于定位，不是MD代表帧。（e）联合Hχ预测与χ1–χ4的MD单角细箱熵关联，依次66、52、23、10个位点，两轴按氨基酸类型去均值；Q41菱形、Q2方形。（f）完整模型在六个Gln中预测最低的Q41与最高的Q2的MD经验投影自由能F*=−ln(p̂/p̂max)，采用501帧、原点0°的5°箱；零计数箱断线，灰虚线为F*=1，底部色条及数字为低值角区总宽度。曲线由MD统计得到，不是模型预测的自由能。b、c采用第一条共同窗口轨迹及2.7节指标；d参考为原五窗口联合熵均值。'
(OUT/'caption_zh.md').write_text(caption+'\n',encoding='utf8')
audit={'status':'numerical_passed_visual_pending','canvas_mm':[W,H],'source_min_font_pt':min(t.get_fontsize() for t in texts),'word_insertion_width_mm':167.64,'expected_min_word_font_pt':min(t.get_fontsize() for t in texts)*167.64/W,'visible_language':'English','data_rows':counts,'all_source_points_retained':True,'subsampling':False,'panel_a_records_per_curve':3063,'panel_b_record_pairs_per_metric':83,'panel_c_records_per_predictor':83,'common_positions':33142,'panel_b_statistics':pair_stats,'panel_c_statistics':selection_ci.reset_index().to_dict(orient='records'),'panel_d_valid_positions':66,'panel_d_NaN_gaps':[10,28,35,46,47,53],'panel_e_n':expected_n,'panel_e_xlim':[-.70,.50],'panel_e_ylim':[-1.,.75],'panel_f_bins':432,'panel_f_original_float_parser':'pd.read_csv default matches previous production','panel_f_zero_bins':int(hist.F_star.isna().sum()),'panel_f_widths':widths,'panel_f_xlim':[-180,180],'panel_f_ylim':[-.8,maxf+1.8],'panel_f_threshold':1,'text_outside_canvas':outside,'panel_axes':panel_boxes,'raster_elements_in_svg':len(svg_images),'raster_source':str(structure_path),'raster_is_actual_PyMOL_structure':True,'source_sha256':{str(p.relative_to(OUT)):sha(p) for folder in [DATA,OUT/'production_asset_sources'] for p in folder.rglob('*') if p.is_file()},'output_sha256':{ext:sha(OUT/(STEM+'.'+ext)) for ext in ['png','pdf','svg']},'caption_characters':len(caption),'interpretation':'Paired correlation difference is descriptive; panel e uses one joint output, not per-angle predictions; panel f is MD reference selected by original joint predictions.','visual_QA':'pending'}
(OUT/'data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','canvas_mm','source_min_font_pt','expected_min_word_font_pt','data_rows','panel_f_zero_bins']},ensure_ascii=False))

