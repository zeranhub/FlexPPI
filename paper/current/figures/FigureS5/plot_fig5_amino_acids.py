# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) GroupedBarChart -> plot_GroupedBarChartv2.py -> param inherit
# (b) GroupedBarChart -> plot_GroupedBarChartv2.py -> param inherit
# Existing asset has hard-coded 2-setting AUROC data; this revision has 18 amino
# acids, two outcomes and record-level bootstrap intervals. Its grouped-column
# and capped-error-bar geometry is inherited, with the project's print style.
# No native-run panel, model inference, trajectory sampling or result selection.
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[4]
sys.path.append(str(ROOT/'experiments/_sources/completed_results_20261012/plot_runtime'))
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

import json,hashlib,shutil
import numpy as np
import pandas as pd
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT=Path(__file__).resolve().parent
DATA=OUT/'source_data';DATA.mkdir(exist_ok=True)
SRC=ROOT/'experiments/_sources'
ENT=SRC/'model_argument_v24_20261013/experiment_predictive_occupancy'
WID=SRC/'free_energy_model_v26_20261013/experiment'
MODEL=['full_hchi','intrachain_hchi','full_rmsf']
ORDER=list('CSTVNDHILFPW YQEMRK'.replace(' ',''))
NAMES={'C':'Cys','S':'Ser','T':'Thr','V':'Val','N':'Asn','D':'Asp','H':'His','I':'Ile','L':'Leu','F':'Phe','P':'Pro','W':'Trp','Y':'Tyr','Q':'Gln','E':'Glu','M':'Met','R':'Arg','K':'Lys'}
assert len(ORDER)==len(set(ORDER))==18
SOURCES={}
def read(path):
    SOURCES[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path,low_memory=False,float_precision='round_trip')

er=read(ENT/'replica1_per_residue.csv.gz')
wr=read(WID/'replica1_per_residue.csv.gz')
keys=['sample_id','residue_index','amino_acid']
merged=er.merge(wr[keys+['eligible_sameAA','width_Fstar_le_1p0_degrees']],on=keys,validate='one_to_one')
assert len(merged)==33493 and merged.sample_id.nunique()==83
eligible=merged[merged.eligible_sameAA].copy()
assert len(eligible)==33142 and set(eligible.amino_acid)==set(ORDER)
eligible.to_csv(DATA/'eligible_residue_values.csv.gz',index=False,compression='gzip')

def rho(x,y):
    a=pd.Series(x).rank(method='average').to_numpy();b=pd.Series(y).rank(method='average').to_numpy()
    assert np.std(a)>0 and np.std(b)>0
    return float(np.corrcoef(a,b)[0,1])

rows=[]
for (sid,aa),g in eligible.groupby(['sample_id','amino_acid'],sort=True):
    assert len(g)>=4
    for metric,col in [('entropy','fine_entropy_mean_per_chi'),('width','width_Fstar_le_1p0_degrees')]:
        for model in MODEL:
            rows.append(dict(sample_id=sid,dataset=g.dataset.iloc[0],amino_acid=aa,metric=metric,predictor=model,n_positions=len(g),value=rho(g[model],g[col])))
cor=pd.DataFrame(rows)
stored=read(ENT/'replica1_per_AA_correlations.csv')
stored=stored[stored.metric=='fine_entropy_mean_per_chi']
check=cor[cor.metric=='entropy'].merge(stored,on=['sample_id','amino_acid','predictor'],validate='one_to_one')
np.testing.assert_allclose(check.value,check.sameAA_rho,atol=2e-14,rtol=0)

eh=read(ENT/'replica1_high_low_per_AA.csv')
eh=eh[(eh.scope=='all')&(eh.metric=='fine_entropy_mean_per_chi')&(eh.predictor=='full_hchi')].copy()
wh=read(WID/'replica1_high_low_per_AA.csv')
wh=wh[(wh.scope=='all')&(wh.threshold==1)&(wh.predictor=='full_hchi')].copy()
assert len(eh)==len(wh)==1296 and eh['rankable'].all() and wh['rankable'].all()
selection=[]
for metric,f in [('entropy',eh),('width',wh)]:
    for _,r in f.iterrows():
        for group,col in [('low','low_MD_mean'),('high','high_MD_mean')]:
            selection.append(dict(sample_id=r.sample_id,dataset=r.dataset,amino_acid=r.amino_acid,metric=metric,predictor=group,n_positions=int(r.positions),value=r[col]))
selection=pd.DataFrame(selection)
assert set(map(tuple,eh[keys[:1]+['amino_acid']].to_numpy()))==set(map(tuple,wh[keys[:1]+['amino_acid']].to_numpy()))
cor.to_csv(DATA/'per_record_amino_acid_correlations.csv',index=False)
selection.to_csv(DATA/'per_record_amino_acid_prediction_groups.csv',index=False)

rng=np.random.default_rng(20261026);summaries=[]
for aa in ORDER:
    one=cor[(cor.amino_acid==aa)&(cor.metric=='entropy')&(cor.predictor=='full_hchi')].sort_values(['dataset','sample_id'])
    ids=one.sample_id.to_numpy();datasets=one.dataset.to_numpy();n=len(ids)
    assert 49<=n<=82
    draws=np.concatenate([rng.choice(np.flatnonzero(datasets==ds),size=(10000,int(sum(datasets==ds))),replace=True) for ds in sorted(set(datasets))],axis=1)
    for panel,frame,predictors in [('a',cor,MODEL),('b',selection,['low','high'])]:
        for metric in ['entropy','width']:
            for pred in predictors:
                f=frame[(frame.amino_acid==aa)&(frame.metric==metric)&(frame.predictor==pred)].set_index('sample_id').loc[ids]
                vals=f.value.to_numpy();assert len(vals)==n and np.isfinite(vals).all()
                lo,hi=np.quantile(vals[draws].mean(axis=1),[.025,.975])
                summaries.append(dict(panel=panel,metric=metric,amino_acid=aa,name=NAMES[aa],predictor=pred,n_records=n,mean=vals.mean(),CI_low=lo,CI_high=hi))
summary=pd.DataFrame(summaries);summary.to_csv(DATA/'plotted_bar_values.csv',index=False)

mpl.rcParams.update({'font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8.5,'axes.unicode_minus':False,'savefig.bbox':None})
fig,axs=plt.subplots(2,2,figsize=(183/25.4,142/25.4))
fig.subplots_adjust(left=.078,right=.99,bottom=.09,top=.825,hspace=.86,wspace=.23)
fig.text(.015,.968,'a',fontsize=11,weight='bold')
fig.text(.052,.968,'Model–MD association within each amino acid',fontsize=10,weight='bold')
fig.text(.015,.489,'b',fontsize=11,weight='bold')
fig.text(.052,.489,'MD observations in model-selected low / high groups',fontsize=10,weight='bold')
fig.text(.99,.968,'83 records; 49–82 records / amino acid',ha='right',va='top',fontsize=7.5)
colors=[CATEGORICAL[0],CATEGORICAL[3],GREY]
fig.legend(handles=[Patch(facecolor=c,edgecolor=BLACK,linewidth=.3,hatch=h,label=l) for c,h,l in zip(colors,['','//','..'],['Full Hχ','Intra Hχ','Full RMSF'])],loc='upper center',bbox_to_anchor=(.51,.938),ncol=3,frameon=False,columnspacing=2.2)
fig.legend(handles=[Patch(facecolor=c,edgecolor=BLACK,linewidth=.3,hatch=h,label=l) for c,h,l in zip([GREY,CATEGORICAL[0]],['//',''],['Low predicted Hχ','High predicted Hχ'])],loc='upper center',bbox_to_anchor=(.51,.463),ncol=2,frameon=False,columnspacing=2.2)
bars=[]
for row,panel in enumerate(['a','b']):
    predictors=MODEL if row==0 else ['low','high']
    for col,metric in enumerate(['entropy','width']):
        ax=axs[row,col];s=summary[(summary.panel==panel)&(summary.metric==metric)]
        width=.235 if row==0 else .32
        for j,pred in enumerate(predictors):
            f=s[s.predictor==pred].set_index('amino_acid').loc[ORDER]
            xx=np.arange(18)+(j-(len(predictors)-1)/2)*width
            color=colors[j] if row==0 else [GREY,CATEGORICAL[0]][j]
            hatch=['','//','..'][j] if row==0 else ['//',''][j]
            bc=ax.bar(xx,f['mean'],width=width*.94,color=color,edgecolor=BLACK,linewidth=.28,hatch=hatch,zorder=2,
                yerr=np.vstack([f['mean']-f.CI_low,f.CI_high-f['mean']]),capsize=1.25,error_kw={'elinewidth':.55,'capthick':.55,'ecolor':BLACK})
            for x,y in zip(xx,f['mean']):bars.append(dict(panel=panel,metric=metric,predictor=pred,x=float(x),height=float(y)))
        ax.set_xlim(-.7,17.7);ax.set_xticks(np.arange(18),[NAMES[x] for x in ORDER],rotation=90)
        ax.tick_params(axis='x',length=0,pad=3)
        lo=float(s.CI_low.min());hi=float(s.CI_high.max())
        ax.set_ylim(min(0,lo*1.15),hi*1.14)
        ax.axhline(0,color=BLACK,lw=.6,zorder=1)
        if row==0:
            ax.set_ylabel('Spearman ρ');ax.set_title('MD fine-bin entropy' if col==0 else 'MD low-F* angular width',pad=8)
        else:
            ax.set_ylabel('MD entropy (nats / χ)' if col==0 else 'MD angular width (° / χ)')
        assert ax.get_ylim()[0]<=lo and ax.get_ylim()[1]>hi
fig.text(.5,.012,'Bars: record means; error bars: 95% source-stratified bootstrap intervals',ha='center',fontsize=8)
stem=OUT/'supplementary_fig5_amino_acid_bars'
fig.savefig(stem.with_suffix('.pdf'),dpi=300)
fig.savefig(stem.with_suffix('.svg'),dpi=300)
fig.savefig(stem.with_suffix('.png'),dpi=300)
plt.close(fig)
audit=dict(source_files=SOURCES,records=83,raw_positions=33493,eligible_positions=33142,excluded_under4_same_type_positions=351,eligible_record_amino_acid_groups=1296,amino_acids=[NAMES[x] for x in ORDER],bars=len(bars),bootstrap={'seed':20261026,'resamples':10000,'unit':'record','stratification':'source within each amino acid','interval':'percentile 2.5–97.5%','same_draws_for_all_compared_predictors':True},new_MD=False,new_model_inference=False,hypothesis_tests='None; descriptive intervals, no p-values or significance claims',archetype='quantitative_grid',dimensions_mm=[183,142],font='Arial',asset_route='GroupedBarChartv2 param inheritance; hard-coded template data/2-setting geometry not reused',all_source_entropy_correlations_reproduced=True,figure_bars=bars)
(OUT/'figure_QA_and_statistics.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
report='''图5统计与复现说明\n\n使用83条原体系五折折外预测记录（ATLAS 30、DynaRepo 23、ITrace 30），每记录第一条已有MD轨迹的0–100 ns、200 ps间隔501帧。全部33,493个侧链位置先应用原方案的记录内同类型至少4个位点条件，保留33,142位、1,296个记录×氨基酸组。351位因同类型数量不足而不参加此项验证；Ala/Gly无χ角，不属于18种可定义侧链读数类型。没有按结果正负或强弱筛选。\n\na：每条记录、每种氨基酸中分别计算模型读数与MD细箱熵或F*≤1角区总宽度的Spearman相关，再对有该氨基酸的记录等权平均。Full Hχ为完整模型侧链读数，Intra Hχ为链内模型侧链读数，Full RMSF为完整模型主链读数。18种类型使用全部49–82条可用记录，计数见plotted_bar_values.csv。此处各氨基酸的相关均值与正文先去类型均值再跨类型计算的83条记录总体相关口径不同，不能简单平均柱高复原正文总体值。\n\nb：先在每条记录、同种氨基酸中仅按Full Hχ预测升序排序，以原残基索引打破并列，最低及最高各取floor(n/2)位，奇数时中间1位不进入两端组；该原有分组规则未改变。分别计算实际MD熵和角区宽度的两组均值，再对记录等权平均。高/低组是模型选择，纵轴是MD观测，不是模型生成的自由能。\n\n误差线为按来源分层重抽记录10,000次的95%百分位区间，种子20261026；同类型内所有比较使用同一批重抽索引。没有新增显著性检验或多重比较星号，残基、χ角及轨迹帧不作为独立重复。细箱熵与角区宽度均先按可定义χ角取平均，原来的5°箱、四个分箱原点及F*阈值1保持不变。\n\nPDF/SVG保留矢量文字和柱图，PNG为300 dpi预览。样式借用分组柱模板的并列柱和带端帽误差线；18种氨基酸与两个观测指标需要自适应布局，未使用模板内的示例AUROC数据。\n'''
(OUT/'统计与复现说明.txt').write_text(report,encoding='utf-8')
print(json.dumps({k:v for k,v in audit.items() if k not in ['source_files','figure_bars']},ensure_ascii=False,indent=2))
