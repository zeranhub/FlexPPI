# Academic Figure Skill Asset Confirmation
# a-d: PairedBoxScatter/plot_PairedBoxScatter1.R point and connection param inheritance
# Existing validated Fig8 Python implementation extended to four actual model tests.
# No native-run R panel; no outcome-based selection; all83records/threepredictors retained.
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

import hashlib,json,xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from PIL import Image
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent
EXP=ROOT/'experiment'
OLD=ROOT.parent/'model_argument_v24_20261013/experiment_predictive_occupancy'
STEM='fig7_model_free_energy_and_occupancy_validation'
font_manager.findfont('Microsoft YaHei',fallback_to_default=False)
mpl.rcParams.update({'font.sans-serif':['Microsoft YaHei','Arial'],'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.unicode_minus':False})
models=['full_hchi','intrachain_hchi','full_rmsf']
names=['Full Hχ','Intra Hχ','Full RMSF']
colors=[CATEGORICAL[0],CATEGORICAL[5],CATEGORICAL[3]]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
sources={}
def read(base,name):
    p=base/name;sources[str(p)]=sha(p);return pd.read_csv(p)
wc=read(EXP,'replica1_per_record_correlations.csv');wh=read(EXP,'replica1_high_low_per_record.csv')
wcs=read(EXP,'replica1_correlation_summary.csv');whs=read(EXP,'replica1_high_low_summary.csv')
ec=read(OLD,'replica1_per_record_correlations.csv');eh=read(OLD,'replica1_high_low_per_record.csv')
ecs=read(OLD,'replica1_correlation_summary.csv');ehs=read(OLD,'replica1_high_low_summary.csv')
width=lambda f:(f.scope=='all')&(f.threshold==1)
entropy=lambda f:(f.scope=='all')&(f.metric=='fine_entropy_mean_per_chi')
items=[
 ('a',wc[width(wc)],wcs[width(wcs)&(wcs.dataset=='all')],'rho_type_centered','预测排序与 MD 低 F* 范围','同类型去均值 Spearman ρ'),
 ('b',ec[entropy(ec)],ecs[entropy(ecs)&(ecs.dataset=='all')],'rho_type_centered','预测排序与 MD 细箱占据熵','同类型去均值 Spearman ρ'),
 ('c',wh[width(wh)],whs[width(whs)&(whs.dataset=='all')],'high_minus_low_MD','预测高低组的 MD 低 F* 范围差','MD 低 F* 总宽度高减低 (°/角)'),
 ('d',eh[entropy(eh)],ehs[entropy(ehs)&(ehs.dataset=='all')],'high_minus_low_MD','预测高低组的 MD 占据熵差','MD 细箱占据熵高减低 (nats/角)')]
fig,axes=plt.subplots(2,2,figsize=(183/25.4,167/25.4))
fig.subplots_adjust(left=.105,right=.985,bottom=.125,top=.865,wspace=.34,hspace=.59)
points=[];vertices=[];stats=[]
offset=np.linspace(-.08,.08,83)
frames=[]
for letter,data,summary,value,title,ylabel in items:
    frame=data.pivot(index='sample_id',columns='predictor',values=value)[models].sort_index()
    assert frame.shape==(83,3) and np.isfinite(frame).all().all()
    frames.append(frame)
assert all(f.index.equals(frames[0].index) for f in frames)
cor_limits=(min(0,min(f.to_numpy().min() for f in frames[:2]))-.04,max(f.to_numpy().max() for f in frames[:2])+.045)
for (letter,data,summary,value,title,ylabel),frame,ax in zip(items,frames,axes.flat):
    for i,(sample,row) in enumerate(frame.iterrows()):
        x=np.arange(3)+offset[i]
        ax.plot(x,row.to_numpy(),color='#A0A0A0',alpha=.24,lw=.4,zorder=1)
        for j,m in enumerate(models):
            points.append({'panel':letter,'sample_id':sample,'predictor':m,'x':x[j],'y':row[m]})
            vertices.append({'panel':letter,'sample_id':sample,'vertex':j,'x':x[j],'y':row[m]})
    labels=[]
    for j,(m,color,name) in enumerate(zip(models,colors,names)):
        ax.scatter(np.full(83,j)+offset,frame[m],s=9,c=color,alpha=.50,linewidths=0,zorder=2)
        s=summary[summary.predictor==m].iloc[0]
        assert int(s.n_records)==83 and abs(frame[m].mean()-s['mean'])<1e-12
        ax.errorbar(j,s['mean'],yerr=[[s['mean']-s.bootstrap_95_low],[s.bootstrap_95_high-s['mean']]],fmt='D',markersize=5.5,color=color,markeredgecolor=BLACK,markeredgewidth=.55,lw=1.4,capsize=3,zorder=5)
        labels.append(name+'\n'+format(s['mean'],'.3f'))
        stats.append({'panel':letter,'predictor':m,'mean':s['mean'],'CI_low':s.bootstrap_95_low,'CI_high':s.bootstrap_95_high,'n_records':83})
    ax.set_xticks(range(3),labels);ax.set_xlim(-.36,2.36)
    if letter in 'ab':limits=cor_limits
    else:
        lower=min(0,frame.to_numpy().min());upper=frame.to_numpy().max();margin=(upper-lower)*.10
        limits=(lower-margin,upper+margin)
    ax.set_ylim(*limits)
    assert frame.to_numpy().min()>limits[0] and frame.to_numpy().max()<limits[1]
    ax.set_ylabel(ylabel,fontsize=8)
    ax.set_title(title,fontsize=8.5,pad=10,loc='left')
    ax.text(-.19,1.035,letter,transform=ax.transAxes,fontsize=11,weight='bold',va='bottom')
    ax.axhline(0,color=GREY,lw=.6,ls='--',zorder=0)
    ax.tick_params(axis='x',length=0,pad=7)
fig.legend(handles=[Line2D([],[],color=GREY,marker='o',lw=.45,markersize=3,label='83 条记录逐一配对'),Line2D([],[],color=BLACK,marker='D',lw=1.3,markersize=4,label='记录均值与 95% 重采样区间')],loc='upper center',bbox_to_anchor=(.54,.995),ncol=2,fontsize=8)
fig.text(.105,.925,'低自由能取向总宽度',fontsize=9.5,weight='bold')
fig.text(.606,.925,'细粒度角度占据多样性',fontsize=9.5,weight='bold')
fig.text(.105,.052,'83 条记录；33,142 个同类型有效位置；原折外预测决定排序及高低组。',fontsize=7.5)
fig.text(.105,.022,'低 F*：相对最常见角区 ≤ 1；5°箱总长度，四个箱格原点及全部 χ 角平均。',fontsize=7.5)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[]
for t in texts:
    box=t.get_window_extent(renderer)
    if box.x0<-.5 or box.y0<-.5 or box.x1>fig.bbox.width+.5 or box.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
assert min(t.get_fontsize() for t in texts)>=7.5
for ext in ['pdf','svg','png']:fig.savefig(OUT/(STEM+'.'+ext),dpi=300,bbox_inches=None)
assert min(Image.open(OUT/(STEM+'.png')).info['dpi'])>=299
svg=ET.parse(OUT/(STEM+'.svg'));assert not svg.findall('.//{http://www.w3.org/2000/svg}image')
pd.DataFrame(points).to_csv(OUT/(STEM+'.points.csv'),index=False,float_format='%.17g')
pd.DataFrame(vertices).to_csv(OUT/(STEM+'.paired_vertices.csv'),index=False,float_format='%.17g')
pd.DataFrame(stats).to_csv(OUT/(STEM+'.statistics.csv'),index=False,float_format='%.17g')
assert len(points)==len(vertices)==996 and len(stats)==12
caption='图7：冻结折外预测对侧链低自由能取向范围与角度占据的群体验证。83条记录统一使用首条轨迹0–100 ns的501帧，同记录同类型至少四个有效位置，共33,142个位点。低自由能取向总宽度为各χ角5°箱中F*=−ln(p̂/p̂max)≤1的角区总长度；细箱占据熵计入全部非零占据箱。两指标均对四个箱格原点及全部可定义角平均。（a,b）各预测读数与相应MD指标分别减去记录内同类型均值后的Spearman相关。（c,d）仅按各预测读数在同类型内排序形成高、低半组，检验两组MD总宽度或细箱熵之差；先对类型等权平均，再对记录等权平均。小点及灰线完整展示83条记录的配对结果，菱形为记录均值，误差线为来源分层重抽记录10,000次的95%百分位区间；横轴数字为均值。Full Hχ、Intra Hχ和Full RMSF分别为完整模型侧链输出、链内模型侧链输出及完整模型主链输出。'
(OUT/(STEM+'.caption_zh.md')).write_text(caption+'\n',encoding='utf8')
audit={'core_claim':'Frozen sidechainprediction identifiesactualMDprimaryanglecoverage andoveralloccupancydifferences amongsameAApositions','panels':{'a':'PredictionrankingvsMDlowFstarwidth','b':'PredictionrankingvsMDfineentropy','c':'Prediction-onlyhighlowselectionvalidatedwithMDwidth','d':'Prediction-onlyhighlowselectionvalidatedwithMDentropy'},'archetype':'quantitative_grid','sources':sources,'statistics':'83recordunits,source-stratified10000bootstrappercentile95%,all33142sameAAeligiblepositions,predeclaredthreshold1,nofreetenergycurvepredictions,not249independentreplicates','dimensions_mm':[183,167],'min_font_pt':min(t.get_fontsize() for t in texts),'text_outside_canvas':outside,'all_actual_points':996,'no_downsampling':True,'svg_raster_images':0,'file_sha256':{ext:sha(OUT/(STEM+'.'+ext)) for ext in ['png','pdf','svg']},'visual_review':'pending'}
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'points':996,'outputs':audit['file_sha256']},ensure_ascii=False))
