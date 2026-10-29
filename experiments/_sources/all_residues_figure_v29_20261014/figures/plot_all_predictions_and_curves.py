# Academic Figure Skill Asset Confirmation
# a: MarginalDensity/plot_VariableCorrelation.py scatter parameter inheritance without marginal overlays.
# b: PairedBoxScatter/plot_PairedBoxScatter1.R point/connection parameter inheritance.
# c-e: original frozen MD Fstar curves plus LineTrend parameter inheritance.
# All six candidate predictions retained; model picks extrema before MD; all432histrows drawn.
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
CASE=ROOT.parent/'free_energy_curves_v27_20261014'/'case_curves'
STEM='fig7_all_residues_Gln_and_MD_free_energy'
PRED=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\neural_exports\evaluation\full\A003R.npz')
PROFILES=Path(r'E:\Project\FlexPPI\experiments\_case_validation_sources\data\ubiquitin_fold0_all_residues.csv')
METRICS=PROFILES.parent/'ubiquitin_UBA_fold0_chain_metrics.csv'
if not PRED.exists():PRED=ROOT/'source/A003R.npz'
if not PROFILES.exists():PROFILES=ROOT/'source/ubiquitin_fold0_all_residues.csv'
if not METRICS.exists():METRICS=ROOT/'source/ubiquitin_UBA_fold0_chain_metrics.csv'
if not (CASE/'histograms_501frames.csv').exists():CASE=ROOT/'source'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(PRED)=='b8a764739993e3bce18dedb3cff493c50a1a52613a08baac16f1803c4edabb5a'
sources={str(p):sha(p) for p in [PRED,PROFILES,METRICS,CASE/'histograms_501frames.csv',CASE/'selection_inventory_with_MD_outcomes.csv']}
z=np.load(PRED,allow_pickle=False)
indices=np.flatnonzero((z['chain']=='A')&z['valid'][:,1])
allrows=pd.DataFrame({'residue_index':indices,'residue_id':z['author_residue_id'][indices].astype(int),'amino_acid':z['amino_acid'][indices],'MD_Hchi_nats':z['target'][indices,1].astype(float),'model_Hchi_nats':z['mean'][indices,1].astype(float)})
assert len(allrows)==66 and np.isfinite(allrows[['MD_Hchi_nats','model_Hchi_nats']]).all().all()
allrows['absolute_error_nats']=abs(allrows.model_Hchi_nats-allrows.MD_Hchi_nats)
gln=allrows[allrows.amino_acid=='Q'].sort_values('residue_id').copy();assert len(gln)==6
low=gln.loc[gln.model_Hchi_nats.idxmin()];high=gln.loc[gln.model_Hchi_nats.idxmax()]
assert low.residue_id==41 and high.residue_id==2
hist=pd.read_csv(CASE/'histograms_501frames.csv');assert len(hist)==432
profile=pd.read_csv(PROFILES);profile=profile[(profile.chain=='A')&profile.valid_hchi].sort_values('residue_id')
np.testing.assert_allclose(profile.md_hchi,allrows.MD_Hchi_nats,rtol=0,atol=1e-6)
np.testing.assert_allclose(profile.pred_hchi,allrows.model_Hchi_nats,rtol=0,atol=1e-6)
rho=allrows.MD_Hchi_nats.rank().corr(allrows.model_Hchi_nats.rank())
corrected=profile.md_hchi_minus_type_mean.rank().corr(profile.pred_hchi_minus_type_mean.rank())
mae=allrows.absolute_error_nats.mean();gmae=gln.absolute_error_nats.mean()
source_metrics=pd.read_csv(METRICS).query('chain=="A"').iloc[0]
assert abs(rho-source_metrics.hchi_spearman)<1e-12 and abs(corrected-source_metrics.hchi_corrected_spearman)<1e-12
font_manager.findfont('Microsoft YaHei',fallback_to_default=False)
mpl.rcParams.update({'font.sans-serif':['Microsoft YaHei','Arial'],'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.unicode_minus':False})
color=lambda site:CATEGORICAL[3] if int(site)==41 else CATEGORICAL[4] if int(site)==2 else CATEGORICAL[5]
fig=plt.figure(figsize=(183/25.4,167/25.4))
a=fig.add_axes([.08,.61,.27,.28]);b=fig.add_axes([.44,.61,.53,.28])
axes=[fig.add_axes([left,.145,.25,.275],sharey=None) for left in [.07,.40,.73]]
fig.text(.08,.975,'完整留出预测 → 同类型全部位点 → 预测端点的 MD 曲线',fontsize=9,weight='bold')
points=[]
limit=4.2
a.plot([0,limit],[0,limit],color=GREY,ls=':',lw=.7,zorder=0)
for _,r in allrows.iterrows():
    a.scatter(r.MD_Hchi_nats,r.model_Hchi_nats,s=15 if r.residue_id not in [2,41] else 29,c=color(r.residue_id),edgecolors='white',linewidths=.35,zorder=2 if r.residue_id not in [2,41] else 4)
    points.append({'panel':'a','residue_id':int(r.residue_id),'amino_acid':r.amino_acid,'kind':'model_vs_MD','x':r.MD_Hchi_nats,'y':r.model_Hchi_nats})
a.set_xlim(-.1,limit);a.set_ylim(-.1,limit);a.set_aspect('equal',adjustable='box')
a.set_xticks([0,1,2,3,4]);a.set_yticks([0,1,2,3,4])
a.set_xlabel('MD 参考 Hχ (nats)');a.set_ylabel('模型预测 Hχ (nats)')
a.set_title('全部 66 个 Hχ 有效残基',loc='left',fontsize=8.5,pad=10)
a.text(-.17,1.055,'a',transform=a.transAxes,fontsize=11,weight='bold')
a.text(.97,.045,'MAE '+format(mae,'.3f')+' nats\nρ '+format(rho,'.3f')+'\n类型均值校正 ρ '+format(corrected,'.3f'),transform=a.transAxes,ha='right',va='bottom',fontsize=7.5)
for j,(_,r) in enumerate(gln.iterrows()):
    x=np.array([j-.13,j+.13]);y=np.array([r.MD_Hchi_nats,r.model_Hchi_nats]);c=color(r.residue_id)
    b.plot(x,y,color=c,lw=.9,zorder=1)
    b.scatter(x[0],y[0],s=35,marker='o',facecolor='white',edgecolor=c,linewidth=.95,zorder=3)
    b.scatter(x[1],y[1],s=33,marker='s',facecolor=c,edgecolor='white',linewidth=.45,zorder=4)
    for kind,xv,yv in zip(['MD_reference','model_prediction'],x,y):points.append({'panel':'b','residue_id':int(r.residue_id),'amino_acid':'Q','kind':kind,'x':xv,'y':yv})
b.set_xlim(-.55,5.55);b.set_ylim(1.5,3.35)
b.set_xticks(range(6),['Q'+str(int(r.residue_id))+'\n'+format(r.absolute_error_nats,'.3f') for _,r in gln.iterrows()])
for label,site in zip(b.get_xticklabels(),gln.residue_id):label.set_color(color(site))
b.set_ylabel('联合 Hχ (nats)');b.set_xlabel('位点；下方数值为绝对误差 (nats)')
b.set_title('全部 6 个谷氨酰胺位点',loc='left',fontsize=8.5,pad=10)
b.text(-.08,1.055,'b',transform=b.transAxes,fontsize=11,weight='bold')
b.text(.97,.94,'平均绝对误差 '+format(gmae,'.3f')+' nats',transform=b.transAxes,ha='right',va='top',fontsize=7.5)
fig.legend(handles=[Line2D([],[],marker='o',color=BLACK,ls='',markerfacecolor='white',markersize=4,label='MD 参考'),Line2D([],[],marker='s',color=BLACK,ls='',markersize=4,label='模型预测')],loc='upper left',bbox_to_anchor=(.435,.96),ncol=2,fontsize=7.5,handletextpad=.4,columnspacing=1.7)
fig.legend(handles=[Line2D([],[],color=CATEGORICAL[3],lw=1.3,label='模型最低 Q41：Hχ = '+format(low.model_Hchi_nats,'.3f')),Line2D([],[],color=CATEGORICAL[4],lw=1.3,ls='--',label='模型最高 Q2：Hχ = '+format(high.model_Hchi_nats,'.3f'))],loc='upper left',bbox_to_anchor=(.063,.515),ncol=2,fontsize=8,handlelength=1.6,columnspacing=2)
maxf=float(hist.F_star.max());display_widths=[]
for qi,ax in enumerate(axes,1):
    letter='cde'[qi-1]
    ax.set_title('χ'+str(qi)+'：MD 自由能投影',loc='left',pad=10,fontsize=8.5)
    ax.text(-.12,1.035,letter,transform=ax.transAxes,fontsize=11,weight='bold')
    for j,(site,c,style) in enumerate([('Q41',CATEGORICAL[3],'-'),('Q2',CATEGORICAL[4],'--')]):
        p=hist[(hist.residue==site)&(hist.chi=='chi'+str(qi))].sort_values('angle_center_degrees');assert len(p)==72
        ax.plot(p.angle_center_degrees,p.F_star,color=c,lw=1.3,ls=style,zorder=3)
        good=np.isfinite(p.F_star)&(p.F_star<=1)
        ax.broken_barh([(float(x),5.) for x in p.loc[good,'angle_lower_degrees']],(-.65+j*.3,.16),facecolors=c,edgecolors='none',zorder=3)
        width=5*int(good.sum());display_widths.append({'residue':site,'chi':qi,'origin0_lowFstar_width_degrees':width})
        ax.text(.03,maxf+1.15-j*.60,site+' 低值角区 '+str(width)+'°',transform=ax.get_yaxis_transform(),color=c,fontsize=7.5,va='center')
        for _,r in p.iterrows():points.append({'panel':letter,'residue_id':int(site[1:]),'amino_acid':'Q','kind':'MD_Fstar','chi':qi,'x':r.angle_center_degrees,'y':r.F_star,'low_Fstar_le1':bool(np.isfinite(r.F_star) and r.F_star<=1)})
    ax.axhline(1,color=GREY,lw=.7,ls=':',zorder=0)
    ax.set_xlim(-180,180);ax.set_ylim(-.8,maxf+1.8);ax.set_xticks([-180,-60,60,180]);ax.set_xlabel('侧链角度 (°)')
    if qi==1:ax.set_ylabel('MD F* (无量纲)')
    else:ax.tick_params(labelleft=False)
fig.text(.07,.026,'底部色条：MD 中 F* ≤ 1 的实际 5°角区；橙色 Q41，紫色 Q2。',fontsize=7.5)
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
allrows.to_csv(OUT/(STEM+'.all66_residue_predictions.csv'),index=False,float_format='%.17g')
gln.to_csv(OUT/(STEM+'.all6_Gln_predictions.csv'),index=False,float_format='%.17g')
pd.DataFrame(display_widths).to_csv(OUT/(STEM+'.display_widths.csv'),index=False)
assert len(points)==510
stats={'all66_Hchi_MAE_nats':float(mae),'all66_Spearman_rho':float(rho),'all66_training_type_mean_corrected_rho':float(corrected),'all6_Gln_MAE_nats':float(gmae),'all6_errors':{int(r.residue_id):float(r.absolute_error_nats) for _,r in gln.iterrows()},'all66_count':66,'GLN_count':6,'valid_zero_reference_count':int(allrows.MD_Hchi_nats.eq(0).sum())}
(OUT/(STEM+'.statistics.json')).write_text(json.dumps(stats,indent=2),encoding='utf8')
caption='图7：泛素侧链预测的全位点评价、同类型对照及MD自由能验证。（a）泛素链72个输入残基中具有Hχ参考值的全部66个残基，点线为预测等于参考；标注平均绝对误差（MAE）、Spearman相关ρ及减去训练折氨基酸类型均值后的ρ。（b）全部六个有效谷氨酰胺位点，空心圆为MD参考，实心方点为模型预测；位点下方数字为绝对误差，单位nats。a、b均使用第0折冻结完整模型及原主评价联合熵参考，参考为五个连续1,000帧窗口的均值。Q41与Q2分别为六个位点中预测最低与最高的位置，始终以橙色和紫色标识，其余位点为灰色。（c–e）这两个预测端点在0–100 ns、200 ps间隔501帧中的χ1、χ2和χ3经验投影自由能F*=−ln(p̂/p̂max)，采用原点0°的5°箱，零计数箱断线。虚点线为F*=1，底部色条及图内标注为该箱格满足F*≤1的全部角区总长度；±180°两端对应同一圆周接缝。'
(OUT/(STEM+'.caption_zh.md')).write_text(caption+'\n',encoding='utf8')
audit={'core_claim':'Allavailable sidechainpredictions andallsameAAcandidatevalues establish predictionaccuracy andselectionbasis beforeactualMDcurves ofpredictedextrema','archetype':'asymmetric_mixed','dimensions_mm':[183,167],'panel_counts':{'a':66,'b':12,'c':144,'d':144,'e':144},'all_points':510,'all66_available_Hchi_positions_including_valid_zero_values_retained':True,'invalid_positions_not_plotted':[(str(z['author_residue_id'][i]),str(z['amino_acid'][i])) for i in np.flatnonzero((z['chain']=='A')&~z['valid'][:,1])],'all6_Gln_including_largest_error_Q40_retained':True,'curve_data_unchanged':True,'source_sha256':sources,'stats':stats,'minimum_font_pt':min(t.get_fontsize() for t in texts),'text_outside_canvas':outside,'no_sampling_smoothing_pvalues_or_CI':True,'b_axis_range_nats':[1.5,3.35],'semantic_colors':'residueidentity: Q41orange/Q2purple/allothersgrey; circleMD/squaremodel; curveorange/purplearebothMD','output_sha256':{ext:sha(OUT/(STEM+'.'+ext)) for ext in ['png','pdf','svg']},'visual_review':'pending'}
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'points':510,'stats':stats,'output_sha256':audit['output_sha256']},ensure_ascii=False))
