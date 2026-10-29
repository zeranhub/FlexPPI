# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Paired model points and means → PairedBoxScatter/plot_PairedBoxScatter1.R → param inherit
# (b) Paired prediction-defined contrasts → PairedBoxScatter/plot_PairedBoxScatter1.R → param inherit
# RULE: No native-run panels. Three-model paired data use R asset connection/point logic in existing Python figure runtime.
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'completed_results_20261012'/'plot_runtime'))
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

import hashlib,json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
OUT=Path(__file__).resolve().parent
ROOT=OUT.parent;EXP=ROOT/'experiment_predictive_occupancy'
STEM='fig8_model_angular_occupancy'
font_manager.findfont('Microsoft YaHei',fallback_to_default=False)
mpl.rcParams.update({'font.sans-serif':['Microsoft YaHei','Arial'],'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.unicode_minus':False})
metric='fine_entropy_mean_per_chi'
models=['full_hchi','intrachain_hchi','full_rmsf']
names=['Full Hχ','Intra Hχ','Full RMSF']
colors=[CATEGORICAL[0],CATEGORICAL[5],CATEGORICAL[3]]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
sources={}
def read(name):
    p=EXP/name;sources[str(p)]=sha(p);return pd.read_csv(p)
c=read('replica1_per_record_correlations.csv');h=read('replica1_high_low_per_record.csv')
cs=read('replica1_correlation_summary.csv');hs=read('replica1_high_low_summary.csv')
condition=lambda f:(f.scope=='all')&(f.metric==metric)
c=c[condition(c)];h=h[condition(h)]
cs=cs[condition(cs)&(cs.dataset=='all')];hs=hs[condition(hs)&(hs.dataset=='all')]
cw=c.pivot(index='sample_id',columns='predictor',values='rho_type_centered')[models].sort_index()
hw=h.pivot(index='sample_id',columns='predictor',values='high_minus_low_MD')[models].reindex(cw.index)
assert cw.shape==hw.shape==(83,3) and np.isfinite(cw).all().all() and np.isfinite(hw).all().all()
offset=np.linspace(-.08,.08,83)
fig,axes=plt.subplots(1,2,figsize=(183/25.4,108/25.4))
fig.subplots_adjust(left=.095,right=.985,bottom=.27,top=.86,wspace=.38)
points=[];vertices=[];stats=[]
for panel,ax,frame,summary,value,title,ylabel,limits in [
    ('a',axes[0],cw,cs,'rho_type_centered','侧链预测与细粒度占据的排序','同类型去均值 Spearman ρ',(-.04,.72)),
    ('b',axes[1],hw,hs,'high_minus_low_MD','按预测选出的高低组实际差异','MD 细箱占据熵高减低 (nats/角)',(-.025,.47))]:
    for i,(sample,row) in enumerate(frame.iterrows()):
        x=np.arange(3)+offset[i]
        ax.plot(x,row.to_numpy(),color='#A0A0A0',alpha=.28,lw=.45,zorder=1)
        for j,m in enumerate(models):
            points.append({'panel':panel,'sample_id':sample,'predictor':m,'x':x[j],'y':row[m]})
            vertices.append({'panel':panel,'sample_id':sample,'vertex':j,'x':x[j],'y':row[m]})
    labels=[]
    for j,(m,color,name) in enumerate(zip(models,colors,names)):
        ax.scatter(np.full(83,j)+offset,frame[m],s=10,c=color,alpha=.48,linewidths=0,zorder=2)
        s=summary[summary.predictor==m].iloc[0]
        assert abs(frame[m].mean()-s['mean'])<1e-12
        ax.errorbar(j,s['mean'],yerr=[[s['mean']-s.bootstrap_95_low],[s.bootstrap_95_high-s['mean']]],fmt='D',markersize=5.6,color=color,markeredgecolor=BLACK,markeredgewidth=.55,lw=1.4,capsize=3,zorder=5)
        labels.append(name+'\n'+format(s['mean'],'.3f'))
        stats.append({'panel':panel,'predictor':m,'mean':s['mean'],'CI_low':s.bootstrap_95_low,'CI_high':s.bootstrap_95_high,'n_records':83})
    ax.set_xticks(range(3),labels)
    ax.set_xlim(-.36,2.36);ax.set_ylim(*limits)
    ax.set_ylabel(ylabel,fontsize=8)
    ax.set_title(title,fontsize=8.5,pad=10,loc='left')
    ax.text(-.18,1.035,panel,transform=ax.transAxes,fontsize=11,weight='bold',va='bottom')
    ax.axhline(0,color=GREY,lw=.6,ls='--',zorder=0)
    ax.tick_params(axis='x',length=0,pad=7)
fig.legend(handles=[Line2D([],[],color=GREY,marker='o',lw=.45,markersize=3,label='逐记录配对'),Line2D([],[],color=BLACK,marker='D',lw=1.3,markersize=4,label='记录均值与 95% 重采样区间')],loc='upper center',bbox_to_anchor=(.54,.995),ncol=2,fontsize=8)
fig.text(.10,.12,'全部 83 条记录；31,? 个同类型有效位置；每记录第一条轨迹。'.replace('31,?','33,142'),fontsize=7.5)
fig.text(.10,.065,'同一类型至少 4 个位点；高、低各取半组；指标为每角 5°细箱熵的均值。',fontsize=7.5)
fig.canvas.draw()
# Export full fixed canvas so final manuscript dimensions are deterministic.
for extension in ['pdf','svg','png']:
    fig.savefig(OUT/(STEM+'.'+extension),dpi=300,bbox_inches=None)
pd.DataFrame(points).to_csv(OUT/(STEM+'.points.csv'),index=False,float_format='%.17g')
pd.DataFrame(vertices).to_csv(OUT/(STEM+'.paired_vertices.csv'),index=False,float_format='%.17g')
pd.DataFrame(stats).to_csv(OUT/(STEM+'.statistics.csv'),index=False,float_format='%.17g')
assert len(points)==498 and len(vertices)==498 and len(stats)==6
caption='图8：冻结折外预测对细粒度侧链角度占据的群体验证。83条记录统一取首条轨迹0–100 ns的501帧。每个χ角以5°细箱统计占据熵，平均四个箱格原点及该残基全部可定义角，单位为nats/角。相同记录内至少四个共同有效位置的氨基酸类型参与，共33,142个位点。（a）各预测读数与轨迹细箱熵分别减去记录内同类型均值后的Spearman相关。（b）仅按各预测读数在同类型内排序形成高、低半组，比较两组真实细箱熵；先对类型等权平均，再对记录等权平均。小点和灰线覆盖全部83条记录，菱形为记录均值，误差线为按来源分层重抽记录10,000次的95%百分位区间；横轴数字为均值。Full Hχ、Intra Hχ和Full RMSF分别表示完整模型侧链输出、链内模型侧链输出和完整模型主链输出。'
(OUT/(STEM+'.caption_zh.md')).write_text(caption+'\n',encoding='utf8')
contract={'core_claim':'冻结侧链预测识别同类型位点细粒度角度占据差异，联系强于模型主链读数。','panels':{'a':'全部记录的同类型残差排序，比较Full/Intra及RMSF','b':'预测确定位置组后的实际轨迹差值，位置由预测选择'},'archetype':'quantitative_grid','dimensions_mm':[183,108],'all_points':498,'no_downsampling':True,'asset':'PairedBoxScatter1.R param inheritance: paired connections and points; three-method data; statistics use frozen record bootstrap rather than asset paired t-test','sources':sources,'statistics':'Record macro mean; source-stratified 10000 bootstrap percentile95%; no pvalue or multiplicity claims; source cohort83/203 frozen; primary eachrecordrep1; OOF correspondingfivefold weights; Full/Intra separatetrained; units specified','predicted_layout_risks_resolved':['Long methodlabels wrapped into two lines','Panel a negative Intra value retained within ylimit','Global legend and denominator definitions placed outside data areas'],'file_sha256':{(STEM+'.'+ext):sha(OUT/(STEM+'.'+ext)) for ext in ['pdf','svg','png']},'visual_review':'pending_actual_PNG_view'}
(OUT/(STEM+'.audit.json')).write_text(json.dumps(contract,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'points':498,'files':contract['file_sha256']},ensure_ascii=False))
