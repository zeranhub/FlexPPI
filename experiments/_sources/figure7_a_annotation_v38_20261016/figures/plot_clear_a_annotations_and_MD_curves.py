# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) All66 frozen joint-Hchi profile -> LineTrend/plot_trend.py -> param inherit
# (b) All151 defined-chi scatter -> cross-type inherit -> param inherit
# (c) Original432 frozen MD Fstar curve rows -> LineTrend/plot_trend.py -> param inherit
# RULE: no native-run panels; retain all source observations without smoothing or selection.
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
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.patches import Rectangle
from PIL import Image

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent;DATA=ROOT/'source';EXP=ROOT/'experiment'
STEM='fig7_clear_a_annotations_and_MD_curves'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
read=lambda p:pd.read_csv(p,float_precision='round_trip')
z=np.load(DATA/'ubiquitin_all66_chi_501frames.npz',allow_pickle=False)
rows=read(DATA/'all66_residue_summary.csv')
stats=readj(DATA/'fig7_all_residues_Gln_and_MD_free_energy.statistics.json')
# Match the v35 production parser exactly for frozen curve float values.
hist=pd.read_csv(DATA/'histograms_501frames.csv')
defined=read(EXP/'defined_angles_long.csv')
corr=read(EXP/'per_chi_correlations.csv')
corr=corr[(corr.predictor=='predicted_joint_Hchi')&corr.metric.isin(['fine_entropy','low_Fstar_width'])].copy()
corr['chi_number']=corr.chi.str.extract(r'(\d+)')[0].astype(int)
assert len(rows)==66 and len(defined)==151 and len(hist)==432 and len(corr)==8
np.testing.assert_allclose(rows.predicted_Hchi_nats,z['predicted_Hchi_nats'],rtol=0,atol=1e-14)
np.testing.assert_allclose(rows.frozen_MD_joint_Hchi_nats,z['frozen_MD_joint_Hchi_nats'],rtol=0,atol=1e-14)
rows['predicted_Hchi_nats']=z['predicted_Hchi_nats']
rows['frozen_MD_joint_Hchi_nats']=z['frozen_MD_joint_Hchi_nats']
expected_n={1:66,2:52,3:23,4:10}
font_manager.findfont('Microsoft YaHei',fallback_to_default=False)
mpl.rcParams.update({'font.sans-serif':['Microsoft YaHei','Arial'],'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.unicode_minus':False})
BLUE=CATEGORICAL[0];ORANGE=CATEGORICAL[3];PURPLE=CATEGORICAL[4]
W,H=183.,147.
fig=plt.figure(figsize=(W/25.4,H/25.4))
axis=lambda x,y,w,h,**kw:fig.add_axes([x/W,y/H,w/W,h/H],**kw)
text=lambda x,y,s,**kw:fig.text(x/W,y/H,s,**kw)
points=[];axis_labels=[];gln_texts=[]

# (a) Original joint-Hchi reference and prediction, all 66 valid positions.
text(3,139,'a',fontsize=11,weight='bold')
text(15,139,'泛素全部 66 个有效位置',fontsize=8.5,weight='bold')
mae_text=text(15,132,'MAE 0.080 nats   ρ 0.979',fontsize=7.5)
type_text=text(15,128,'扣除氨基酸类型均值后 ρ 0.868',fontsize=7.5)
pa=axis(15,82,82,39)
ids=np.arange(1,73);model=np.full(72,np.nan);reference=model.copy()
for _,r in rows.iterrows():
    pos=int(r.author_residue_id)-1
    model[pos]=r.predicted_Hchi_nats;reference[pos]=r.frozen_MD_joint_Hchi_nats
    for kind,value in [('model_prediction',model[pos]),('MD_reference',reference[pos])]:
        points.append(dict(panel='a_Hchi',residue_id=pos+1,kind=kind,x=pos+1,y=value))
for site,color in [(41,ORANGE),(2,PURPLE)]:pa.axvspan(site-.5,site+.5,color=color,alpha=.16,lw=0,zorder=0)
pa.plot(ids,reference,color=BLACK,lw=.7,ls=':',marker='o',markersize=2.4,markerfacecolor='white',markeredgewidth=.55,zorder=2)
pa.plot(ids,model,color=BLUE,lw=.9,marker='s',markersize=2.2,markeredgewidth=0,zorder=3)
pa.set_xlim(.5,72.5);pa.set_ylim(-.15,4.25);pa.set_yticks([0,2,4]);pa.set_ylabel('联合 Hχ\n(nats)',fontsize=7.5)
pa.set_xticks([1,10,20,30,40,50,60,72])
for site,color in [(2,PURPLE),(41,ORANGE)]:
    pa.plot([site,site],[0,-.045],transform=pa.get_xaxis_transform(),color=color,lw=.8,clip_on=False)
    t=pa.text(site,-.22,'Q'+str(site),transform=pa.get_xaxis_transform(),ha='center',va='center',fontsize=7.5,color=color,clip_on=False)
    gln_texts.append(t);axis_labels.append(dict(residue_id=site,anchor_x_data=site,label_x_data=site,y_axes=-.22))
handles=[Line2D([0],[0],color=BLACK,ls=':',lw=.7,marker='o',markersize=3.2,markerfacecolor='white',markeredgewidth=.7),Line2D([0],[0],color=BLUE,lw=.9,marker='s',markersize=3,markeredgewidth=0)]
joint_legend=pa.legend(handles,['MD 参考','模型预测'],loc='lower right',bbox_to_anchor=(1.,1.025),ncol=2,fontsize=7.5,handlelength=1.5,handletextpad=.45,columnspacing=1.2,borderaxespad=0)

# (b) Every defined residue-angle observation, with one common linear scale.
text(100,139,'b',fontsize=11,weight='bold')
text(111,139,'逐χ MD占据与联合预测',fontsize=8.5,weight='bold')
scatter_axes=[]
for chi,left,bottom in [(1,110,105),(2,148,105),(3,110,74),(4,148,74)]:
    ax=axis(left,bottom,30,23)
    scatter_axes.append(ax)
    p=defined[defined.chi=='chi'+str(chi)]
    assert len(p)==expected_n[chi]
    r=corr[(corr.chi_number==chi)&(corr.metric=='fine_entropy')].iloc[0]
    ax.set_title('χ'+str(chi)+'  n='+str(len(p))+'  ρ='+format(float(r.rho_AA_centered),'.2f'),loc='left',fontsize=7.5,pad=4)
    for label,group in [('other',p[~p.author_residue_id.isin([2,41])]),('Q2',p[p.author_residue_id==2]),('Q41',p[p.author_residue_id==41])]:
        color=PURPLE if label=='Q2' else ORANGE if label=='Q41' else BLUE
        ax.scatter(group.predicted_Hchi_nats_AA_centered,group.fine_entropy_four_origins_mean_nats_AA_centered,s=12 if label=='other' else 20,color=color,edgecolors='white' if label=='other' else color,linewidths=.25 if label=='other' else .5,alpha=.75 if label=='other' else 1,zorder=3 if label=='other' else 4)
    for _,v in p.iterrows():
        points.append(dict(panel='b_scatter',residue_id=int(v.author_residue_id),chi=chi,kind='fine_entropy_same_type_deviation_pair',x=float(v.predicted_Hchi_nats_AA_centered),y=float(v.fine_entropy_four_origins_mean_nats_AA_centered),n_defined=len(p),is_proline=bool(v.is_proline)))
    ax.axhline(0,color=GREY,lw=.45,ls=':',zorder=0)
    ax.axvline(0,color=GREY,lw=.45,ls=':',zorder=0)
    ax.set_xlim(-.70,.50);ax.set_ylim(-1.,.75)
    ax.set_xticks([-.5,0,.5]);ax.set_yticks([-.5,0,.5]);ax.tick_params(pad=2,length=2.2)
    if chi in [1,2]:ax.tick_params(labelbottom=False)
    if chi in [2,4]:ax.tick_params(labelleft=False)
text(101,101.5,'MD 细箱熵：同类型偏差 (nats)',fontsize=7.5,rotation=90,ha='center',va='center')
text(144,64,'联合 Hχ 预测：同类型偏差 (nats)',fontsize=7.5,ha='center')

# (c) One panel containing the three unchanged MD projections.
text(3,54,'c',fontsize=11,weight='bold')
text(17,54,'预测端点的逐角 MD 自由能',fontsize=8.5,weight='bold')
text(17,48.5,'六个谷氨酰胺中预测最低 Q41：Hχ = 2.000',fontsize=7.5,color=ORANGE)
text(100,48.5,'预测最高 Q2：Hχ = 2.990',fontsize=7.5,color=PURPLE)
maxf=float(hist.F_star.max());widths=[];curve_axes=[]
for chi,left in enumerate([20,74,128],1):
    ax=axis(left,10,43,31);curve_axes.append(ax)
    ax.set_title('χ'+str(chi),loc='left',pad=5,fontsize=8)
    for j,(site,color,style) in enumerate([('Q41',ORANGE,'-'),('Q2',PURPLE,'--')]):
        p=hist[(hist.residue==site)&(hist.chi=='chi'+str(chi))].sort_values('angle_center_degrees')
        assert len(p)==72
        ax.plot(p.angle_center_degrees,p.F_star,color=color,lw=1.3,ls=style,zorder=3)
        good=np.isfinite(p.F_star)&(p.F_star<=1)
        ax.broken_barh([(float(v),5.) for v in p.loc[good,'angle_lower_degrees']],(-.65+j*.3,.16),facecolors=color,edgecolors='none',zorder=3)
        width=5*int(good.sum());widths.append(dict(residue=site,chi=chi,origin0_lowFstar_width_degrees=width))
        ax.text(.03,maxf+1.15-j*.60,site+' 低值角区 '+str(width)+'°',transform=ax.get_yaxis_transform(),color=color,fontsize=7.5,va='center')
        for _,v in p.iterrows():points.append(dict(panel='c',residue_id=int(site[1:]),chi=chi,kind='MD_Fstar',x=v.angle_center_degrees,y=v.F_star,low_Fstar_le1=bool(np.isfinite(v.F_star) and v.F_star<=1)))
    ax.axhline(1,color=GREY,lw=.7,ls=':',zorder=0)
    ax.set_xlim(-180,180);ax.set_ylim(-.8,maxf+1.8);ax.set_xticks([-180,-60,60,180]);ax.set_xlabel('侧链角度 (°)',fontsize=7.5)
    if chi==1:ax.set_ylabel('MD F* (无量纲)',fontsize=7.5)
    else:ax.tick_params(labelleft=False)

fig.canvas.draw();renderer=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[]
for t in texts:
    box=t.get_window_extent(renderer)
    if box.x0<-.5 or box.y0<-.5 or box.x1>fig.bbox.width+.5 or box.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
legend_box=joint_legend.get_window_extent(renderer)
assert not legend_box.overlaps(type_text.get_window_extent(renderer))
assert not legend_box.overlaps(mae_text.get_window_extent(renderer))
assert legend_box.y0>pa.bbox.y1
assert pa.get_xlabel()==''
missing_positions=[int(i) for i in ids[~np.isfinite(model)]]
assert len(missing_positions)==6 and np.array_equal(np.isnan(model),np.isnan(reference))
assert min(t.get_fontsize() for t in texts)>=7.5
assert len(points)==715
for ext in ['pdf','svg','png']:fig.savefig(OUT/(STEM+'.'+ext),dpi=300,bbox_inches=None)
assert not ET.parse(OUT/(STEM+'.svg')).findall('.//{http://www.w3.org/2000/svg}image')
pd.DataFrame(points).to_csv(OUT/(STEM+'.points.csv'),index=False,float_format='%.17g')
pd.DataFrame(widths).to_csv(OUT/(STEM+'.display_widths.csv'),index=False)
defined.to_csv(OUT/(STEM+'.all151_scatter_source.csv'),index=False,float_format='%.17g')
caption='图7：联合熵预测、逐角占据关联与MD自由能。（a）泛素全部66个Hχ有效残基的联合熵预测与MD参考，仅标出Q2与Q41；空心圆为MD参考，蓝色方点为模型预测。折线为原始Hχ；扣除类型均值后的ρ由预测与MD参考分别减去训练折中该氨基酸类型的同一个均值后计算。全部位置均属泛素A链，无可定义χ角的Gly/Ala位置不绘点，连线在此断开。（b）χ1–χ4全部有效残基的模型联合Hχ与MD单角细箱占据熵散点，n依次为66、52、23、10。每个χ层内按氨基酸类型分别减去预测值与MD指标的均值，散点使用相同坐标范围；Q2和Q41以紫色和橙色标识。各小图标注细箱熵的Spearmanρ。单角细箱熵采用0–100 ns、200 ps间隔的501帧，并对5°箱的四个原点取平均。（c）泛素六个谷氨酰胺中联合Hχ预测最低Q41与最高Q2的χ1、χ2及χ3逐角MD自由能F*=−ln(p̂/p̂max)，采用相同501帧、原点0°的5°箱，零计数箱断线。虚点线为F*=1，底部色条与标注为该箱格低自由能角区总宽度，±180°对应同一圆周接缝。联合Hχ参考为原五个连续1,000帧窗口的均值；模型为第0折冻结完整模型。'
(OUT/(STEM+'.caption_zh.md')).write_text(caption+'\n',encoding='utf8')
audit={'status':'numerical_checks_passed_visual_pending','configured_canvas_mm':[W,H],'panel_letters':['a','b','c'],'a_labels':axis_labels,'all66_joint_values_retained':True,'all151_defined_angles_retained':True,'scatter_n':expected_n,'scatter_xlim':[-.7,.5],'scatter_ylim':[-1.,.75],'b_summary_table_removed':True,'original8correlations_preserved_in_experiment':True,'c_432_Fstar_values_retained':True,'c_curve_xlim':[-180,180],'c_curve_ylim':[-.8,maxf+1.8],'minimum_font_pt':min(t.get_fontsize() for t in texts),'text_outside_canvas':outside,'point_counts':pd.DataFrame(points).groupby('panel').size().to_dict(),'total_points':len(points),'frozen_overall_metrics':stats,'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in list(DATA.iterdir())+list(EXP.iterdir()) if p.is_file()},'output_sha256':{ext:sha(OUT/(STEM+'.'+ext)) for ext in ['png','pdf','svg']},'visual_review':'pending','parent_figure_version':'v37','a_xaxis_title_removed':True,'a_legend_location':'above_curve_axis_upper_right','a_legend_bbox_pixels':list(legend_box.extents),'a_legend_above_axis':True,'a_legend_no_overlap_with_statistics':True,'a_type_mean_label':'扣除氨基酸类型均值后 ρ 0.868','a_missing_positions':missing_positions,'a_NaN_breaks_preserved':True,'a_missing_residue_names':{'10':'GLY','28':'ALA','35':'GLY','46':'ALA','47':'GLY','53':'GLY'},'a_type_mean_correction':'Subtract the same fold0 training reference mean for each amino-acid type from both prediction and MD reference; then compute Spearman rho. This is distinct from panel b within-chi-layer means.'}
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k not in ['source_sha256','frozen_overall_metrics']},ensure_ascii=False))



