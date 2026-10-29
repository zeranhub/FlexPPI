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

W,H=183.,202.
fig=plt.figure(figsize=(W/25.4,H/25.4))
axis=lambda x,y,w,h,**kw:fig.add_axes([x/W,y/H,w/W,h/H],**kw)
text=lambda x,y,s,**kw:fig.text(x/W,y/H,s,**kw)
panel_boxes={}; axes=[]; points=[]
def panel(letter,x,y,title):
    text(x,y,letter,fontsize=12,weight='bold',va='baseline')
    text(x+8,y,title,fontsize=9.5,weight='bold',va='baseline')
def register(name,ax):
    axes.append(ax);panel_boxes[name]=list(ax.get_position().bounds)

# a. All 3063 records, unchanged original ECDF vertices: relation, not accuracy.
panel('a',2,195,'Two readouts')
text(15,188,'3,063 records',fontsize=9)
ax=axis(15,152,37,28);register('a',ax)
for role,color,ls,name in [('reference',GREY,'--','MD'),('prediction',BLUE,'-','Model')]:
    p=ecdf[ecdf.role==role].sort_values('y')
    ax.plot(p.x,p.y,color=color,lw=1.05,ls=ls)
    for _,v in p.iterrows():points.append({'panel':'a','sample_id':v.sample_id,'kind':role,'x':v.x,'y':v.y})
    mean=float(p.x.mean())
    ax.text(.03,.92 if role=='reference' else .73,name+' '+format(mean,'.3f'),transform=ax.transAxes,color=color,fontsize=9)
ax.set_xlim(-.25,1);ax.set_ylim(0,1.02);ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
ax.set_xlabel('RMSF–Hχ Spearman ρ',labelpad=3);ax.set_ylabel('Fraction of records',labelpad=3)

# b. Every paired per-record rho difference; no conditional-independence claim.
panel('b',60,195,'Occupancy association')
text(77,188,'83 paired records',fontsize=9)
ax=axis(83,152,36,28);register('b',ax)
for i,st in enumerate(pair_stats):
    y=1-i;p=paired[paired.metric==st['metric']].sort_values('sample_id').reset_index(drop=True)
    jitter=(np.arange(len(p))*37%83/82-.5)*.26
    ax.scatter(p.paired_rho_difference,y+jitter,s=8,color=GREY,alpha=.45,linewidths=0,zorder=2)
    ax.errorbar(st['mean_difference'],y,xerr=np.array([[st['mean_difference']-st['CI_low']],[st['CI_high']-st['mean_difference']]]),
        fmt='D',ms=4,color=BLUE,ecolor=BLUE,elinewidth=1.3,capsize=2,zorder=4)
    ax.text(.98,y+.28,format(st['mean_difference'],'.3f'),transform=ax.get_yaxis_transform(),ha='right',color=BLUE,fontsize=9)
    for j,v in p.iterrows():points.append({'panel':'b','sample_id':v.sample_id,'kind':v.metric,'x':v.paired_rho_difference,'y':float(y+jitter[j])})
ax.axvline(0,color=GREY,lw=.6,ls=':',zorder=0)
ax.set_yticks([1,0],['Bin entropy','Low-F* width']);ax.set_ylim(-.5,1.6);ax.set_xlim(-.25,.55);ax.set_xticks([0,.25,.5])
ax.set_xlabel('Δρ: Hχ − RMSF',labelpad=3);ax.tick_params(axis='y',length=0,pad=3)

# c. All83 predictions determine high/low groups, actual MD width evaluates them.
panel('c',125,195,'Selected angular range')
text(141,188,'83 records',fontsize=9)
ax=axis(141,152,35,28);register('c',ax)
for i,(pred,label,color,marker,filled) in enumerate([
    ('full_hchi','Hχ',BLUE,'D',True),('intrachain_hchi','Intra Hχ',BLUE,'o',False),('full_rmsf','RMSF',BLACK,'s',True)]):
    y=2-i;p=selection[selection.predictor==pred].sort_values('sample_id').reset_index(drop=True)
    jitter=(np.arange(len(p))*37%83/82-.5)*.26
    ax.scatter(p.high_minus_low_MD,y+jitter,s=8,color=GREY,alpha=.4,linewidths=0,zorder=2)
    st=selection_ci.loc[pred];mean=float(st['mean'])
    ax.errorbar(mean,y,xerr=np.array([[mean-st.bootstrap_95_low],[st.bootstrap_95_high-mean]]),
        fmt=marker,ms=4,color=color,mfc=color if filled else 'white',ecolor=color,elinewidth=1.2,capsize=2,zorder=4)
    ax.text(.98,y+.27,format(mean,'.2f'),transform=ax.get_yaxis_transform(),ha='right',fontsize=9,color=color)
    for j,v in p.iterrows():points.append({'panel':'c','sample_id':v.sample_id,'kind':pred,'x':v.high_minus_low_MD,'y':float(y+jitter[j])})
ax.axvline(0,color=GREY,lw=.6,ls=':',zorder=0)
ax.set_yticks([2,1,0],['Hχ','Intra Hχ','RMSF']);ax.set_ylim(-.5,2.6);ax.set_xlim(-1.5,17.5);ax.set_xticks([0,8,16])
ax.set_xlabel('MD width difference (°/χ)',labelpad=3);ax.tick_params(axis='y',length=0,pad=3)

# d. All66 original reference/prediction values and all NaN gaps.
panel('d',2,135,'Ubiquitin: all 66 valid residues')
text(15,129,'MAE 0.080 nats   ρ 0.979',fontsize=9)
ax=axis(15,85,80,34);register('d',ax)
ids=np.arange(1,73);model=np.full(72,np.nan);reference=model.copy()
for _,r in rows.iterrows():
    pos=int(r.author_residue_id)-1;model[pos]=r.predicted_Hchi_nats;reference[pos]=r.frozen_MD_joint_Hchi_nats
    for name,value in [('prediction',model[pos]),('MD',reference[pos])]:points.append({'panel':'d','residue_id':pos+1,'kind':name,'x':pos+1,'y':value})
for site,color in [(41,ORANGE),(2,PURPLE)]:ax.axvspan(site-.5,site+.5,color=color,alpha=.16,lw=0,zorder=0)
ax.plot(ids,reference,color=BLACK,lw=.7,ls=':',marker='o',ms=2.4,mfc='white',mew=.55,zorder=2)
ax.plot(ids,model,color=BLUE,lw=.9,marker='s',ms=2.2,mew=0,zorder=3)
ax.set_xlim(.5,72.5);ax.set_ylim(-.15,4.25);ax.set_yticks([0,2,4]);ax.set_ylabel('Joint Hχ (nats)',labelpad=3)
ax.set_xticks([1,10,20,30,40,50,60,72])
for site,color in [(2,PURPLE),(41,ORANGE)]:
    ax.plot([site,site],[0,-.045],transform=ax.get_xaxis_transform(),color=color,lw=.8,clip_on=False)
    ax.text(site,-.25,'Q'+str(site),transform=ax.get_xaxis_transform(),ha='center',va='center',fontsize=9,color=color,clip_on=False)
handles=[Line2D([0],[0],color=BLACK,ls=':',lw=.7,marker='o',ms=3,mfc='white',mew=.7),Line2D([0],[0],color=BLUE,lw=.9,marker='s',ms=3,mew=0)]
leg=ax.legend(handles,['MD reference','Model prediction'],loc='lower right',bbox_to_anchor=(1.,1.035),ncol=2,
    fontsize=9,handlelength=1.4,handletextpad=.4,columnspacing=1.,borderaxespad=0)

# e. The original151 defined residue-angle observations, no thinning or selection.
panel('e',100,135,'Per-χ occupancy association')
for chi,left,bottom in [(1,114,108),(2,151,108),(3,114,82),(4,151,82)]:
    ax=axis(left,bottom,26,18);register('e_chi'+str(chi),ax)
    p=defined[defined.chi=='chi'+str(chi)]
    assert len(p)==expected_n[chi]
    r=corr[corr.chi=='chi'+str(chi)].iloc[0]
    ax.set_title('χ'+str(chi)+'  n='+str(len(p))+'  ρ='+format(float(r.rho_AA_centered),'.2f'),loc='left',fontsize=9,pad=4)
    for label,g in [('other',p[~p.author_residue_id.isin([2,41])]),('Q2',p[p.author_residue_id==2]),('Q41',p[p.author_residue_id==41])]:
        color=PURPLE if label=='Q2' else ORANGE if label=='Q41' else BLUE
        ax.scatter(g.predicted_Hchi_nats_AA_centered,g.fine_entropy_four_origins_mean_nats_AA_centered,
            s=11 if label=='other' else 20,color=color,edgecolors='white' if label=='other' else color,
            linewidths=.25 if label=='other' else .5,alpha=.75 if label=='other' else 1,zorder=3 if label=='other' else 4)
    for _,v in p.iterrows():points.append({'panel':'e','residue_id':int(v.author_residue_id),'chi':chi,'kind':'fine_entropy_pair',
        'x':float(v.predicted_Hchi_nats_AA_centered),'y':float(v.fine_entropy_four_origins_mean_nats_AA_centered)})
    ax.axhline(0,color=GREY,lw=.45,ls=':',zorder=0);ax.axvline(0,color=GREY,lw=.45,ls=':',zorder=0)
    ax.set_xlim(-.70,.50);ax.set_ylim(-1.,.75);ax.set_xticks([-.5,0,.5]);ax.set_yticks([-.5,0,.5]);ax.tick_params(pad=2,length=2)
    if chi in [1,2]:ax.tick_params(labelbottom=False)
    if chi in [2,4]:ax.tick_params(labelleft=False)
text(102,105,'MD entropy deviation (nats)',fontsize=9,rotation=90,ha='center',va='center')
text(145,72,'Predicted joint Hχ deviation (nats)',fontsize=9,ha='center')

# f. Exactly the previous432 bins, original zero-bin NaNs, F* threshold and bars.
panel('f',2,61,'MD free energy at prediction-selected sites')
text(17,55,'Lowest of six Gln: Q41, Hχ = 2.000',fontsize=9,color=ORANGE)
text(100,55,'Highest: Q2, Hχ = 2.990',fontsize=9,color=PURPLE)
maxf=float(hist.F_star.max());widths=[]
for chi,left in enumerate([20,75,130],1):
    ax=axis(left,13,43,31);register('f_chi'+str(chi),ax)
    ax.set_title('χ'+str(chi),loc='left',pad=5,fontsize=9)
    for j,(site,color,style) in enumerate([('Q41',ORANGE,'-'),('Q2',PURPLE,'--')]):
        p=hist[(hist.residue==site)&(hist.chi=='chi'+str(chi))].sort_values('angle_center_degrees')
        assert len(p)==72
        ax.plot(p.angle_center_degrees,p.F_star,color=color,lw=1.3,ls=style,zorder=3)
        good=np.isfinite(p.F_star)&(p.F_star<=1)
        ax.broken_barh([(float(v),5.) for v in p.loc[good,'angle_lower_degrees']],(-.65+j*.3,.16),facecolors=color,edgecolors='none',zorder=3)
        width=5*int(good.sum());widths.append({'residue':site,'chi':chi,'origin0_lowF_width_degrees':width})
        ax.text(.03,maxf+1.15-j*.90,site+' width '+str(width)+'°',transform=ax.get_yaxis_transform(),color=color,fontsize=9,va='center')
        for _,v in p.iterrows():points.append({'panel':'f','residue_id':int(site[1:]),'chi':chi,'kind':'MD_Fstar','x':v.angle_center_degrees,'y':v.F_star,'low_Fstar_le1':bool(np.isfinite(v.F_star) and v.F_star<=1)})
    ax.axhline(1,color=GREY,lw=.7,ls=':',zorder=0)
    ax.set_xlim(-180,180);ax.set_ylim(-.8,maxf+1.8);ax.set_xticks([-180,-60,60,180]);ax.set_xlabel('Side-chain angle (°)',fontsize=9,labelpad=3)
    if chi==1:ax.set_ylabel('MD F* (dimensionless)',fontsize=9,labelpad=3)
    else:ax.tick_params(labelleft=False)

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
for ext in ['pdf','svg','png']:fig.savefig(OUT/(STEM+'.'+ext),dpi=450,bbox_inches=None)
assert not ET.parse(OUT/(STEM+'.svg')).findall('.//{http://www.w3.org/2000/svg}image')
pd.DataFrame(points).to_csv(OUT/'all_plotted_points.csv',index=False,float_format='%.17g')
paired.to_csv(OUT/'panel_b_all83_paired_rho_differences.csv',index=False,float_format='%.17g')
selection.to_csv(OUT/'panel_c_all83_selection_differences.csv',index=False,float_format='%.17g')
pd.DataFrame(pair_stats).to_csv(OUT/'panel_b_statistics.csv',index=False,float_format='%.17g')
pd.DataFrame(widths).to_csv(OUT/'panel_f_unchanged_display_widths.csv',index=False)
defined.to_csv(OUT/'panel_e_all151_scatter_source.csv',index=False,float_format='%.17g')
(OUT/'visible_text.json').write_text(json.dumps([t.get_text() for t in texts],ensure_ascii=False,indent=2),encoding='utf8')
caption='图5：联合动态读数与MD构象占据的对应。（a）3,063条记录中RMSF与Hχ的同类型相关分布，灰/蓝为MD/预测，标注为记录均值。（b）83条记录、33,142个同类型有效位置中，Hχ与RMSF预测对MD细箱熵、低自由能角区总宽度的相关之差，每点为一条记录。（c）仅按预测划分同类残基高低半组，显示全部83条记录的真实角区总宽度组差；正值表示高预测组更宽。b、c大点及误差线为均值及95%来源分层重采样区间；Intra为链内Hχ，其余为完整模型。（d）泛素全部66个Hχ有效位置的预测与MD参考，Gly/Ala处断线。（e）联合Hχ预测与χ1–χ4的MD单角细箱熵关联，依次66、52、23、10个位点；两轴按氨基酸类型去均值，Q2/41为紫/橙色。（f）预测选出的六个Gln中最低Q41与最高Q2的MD经验投影自由能F*=−ln(p̂/p̂max)，采用501帧、原点0°的5°箱；零计数箱断线，虚线为F*=1，底部色条及数字为低值角区总宽度。b、c采用第一条共同窗口轨迹及2.8节指标；d参考为原五窗口联合熵均值。'
(OUT/'caption_zh.md').write_text(caption+'\n',encoding='utf8')
audit={'status':'numerical_passed_visual_pending','canvas_mm':[W,H],'source_min_font_pt':min(t.get_fontsize() for t in texts),
    'word_insertion_width_mm':167.64,'expected_min_word_font_pt':min(t.get_fontsize() for t in texts)*167.64/W,
    'visible_language':'English','data_rows':counts,'all_source_points_retained':True,'subsampling':False,
    'panel_a_records_per_curve':3063,'panel_b_record_pairs_per_metric':83,'panel_c_records_per_predictor':83,
    'common_positions':33142,'panel_b_statistics':pair_stats,'panel_c_statistics':selection_ci.reset_index().to_dict(orient='records'),
    'panel_d_valid_positions':66,'panel_d_NaN_gaps':[10,28,35,46,47,53],'panel_e_n':expected_n,
    'panel_f_bins':432,'panel_f_original_float_parser':'pd.read_csv default matches previous production',
    'panel_f_zero_bins':int(hist.F_star.isna().sum()),'panel_f_widths':widths,'panel_f_xlim':[-180,180],
    'panel_f_ylim':[-.8,maxf+1.8],'panel_f_threshold':1,'text_outside_canvas':outside,
    'panel_axes':panel_boxes,'source_sha256':{str(p.relative_to(OUT)):sha(p) for folder in [DATA,OUT/'production_asset_sources'] for p in folder.rglob('*') if p.is_file()},
    'output_sha256':{ext:sha(OUT/(STEM+'.'+ext)) for ext in ['png','pdf','svg']},'caption_characters':len(caption),
    'interpretation':'Paired correlation difference is descriptive, not conditional independent contribution; panel e uses one joint output, not per-angle predictions; panel f is MD reference selected by original joint predictions.',
    'visual_QA':'pending'}
(OUT/'data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','canvas_mm','source_min_font_pt','expected_min_word_font_pt','data_rows','panel_b_statistics','caption_characters']},ensure_ascii=False))
