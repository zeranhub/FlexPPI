# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a left/right) complete66 Hchi profile and all4chi x2metric correlation summary -> LineTrend/plot_trend.py and cross-type inherit -> param inherit
# (b-d) original frozen MD Fstar curves -> LineTrend/plot_trend.py -> param inherit
# RULE: no native-run panels; production assets have incompatible application-specific data structures.
# Retain all66 joint-Hchi positions, all4chi layers x2 pre-existing MDmetrics, and original432MDcurvepoints; no smoothing or selection.
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
DATA=ROOT/'source';OLD=DATA;CASE=DATA
STEM='fig7_all_chi_associations_and_MD_curves'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
source=DATA/'ubiquitin_all66_chi_501frames.npz'
z=np.load(source,allow_pickle=False)
rows=pd.read_csv(DATA/'all66_residue_summary.csv')
stats=readj(OLD/'fig7_all_residues_Gln_and_MD_free_energy.statistics.json')
gln=rows[rows.amino_acid=='Q'].sort_values('author_residue_id')
assert len(rows)==66 and len(gln)==6 and list(gln.author_residue_id)==[2,31,40,41,49,62]
assert np.array_equal(rows.residue_index,z['residue_index'])
np.testing.assert_allclose(rows.predicted_Hchi_nats,z['predicted_Hchi_nats'],rtol=0,atol=1e-14)
np.testing.assert_allclose(rows.frozen_MD_joint_Hchi_nats,z['frozen_MD_joint_Hchi_nats'],rtol=0,atol=1e-14)
rows['predicted_Hchi_nats']=z['predicted_Hchi_nats']
rows['frozen_MD_joint_Hchi_nats']=z['frozen_MD_joint_Hchi_nats']
gln=rows[rows.amino_acid=='Q'].sort_values('author_residue_id')
hist=pd.read_csv(CASE/'histograms_501frames.csv');assert len(hist)==432
font_manager.findfont('Microsoft YaHei',fallback_to_default=False)
mpl.rcParams.update({'font.sans-serif':['Microsoft YaHei','Arial'],'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.unicode_minus':False})
BLUE=CATEGORICAL[0];ORANGE=CATEGORICAL[3];PURPLE=CATEGORICAL[4]
identity=lambda site: ORANGE if int(site)==41 else PURPLE if int(site)==2 else BLACK
W,H=183.,135.
fig=plt.figure(figsize=(W/25.4,H/25.4))
def axis(x,y,w,h,**kw):return fig.add_axes([x/W,y/H,w/W,h/H],**kw)
def text(x,y,s,**kw):return fig.text(x/W,y/H,s,**kw)
points=[]
text(3,128,'a',fontsize=11,weight='bold')
text(17,128,'泛素全部 66 个有效位置',fontsize=8.5,weight='bold')
text(17,120,'MAE 0.080 nats   ρ 0.979',fontsize=7.5)
text(17,116,'类型均值校正 ρ 0.868',fontsize=7.5)
pa=axis(17,79,100,31)
ids=np.arange(1,73);model=np.full(72,np.nan);reference=model.copy()
for i,r in rows.iterrows():
    pos=int(r.author_residue_id)-1;model[pos]=r.predicted_Hchi_nats;reference[pos]=r.frozen_MD_joint_Hchi_nats
    for kind,y in [('model_prediction',model[pos]),('MD_reference',reference[pos])]:points.append(dict(panel='a_Hchi',residue_id=pos+1,kind=kind,x=pos+1,y=y))
pa.set_xlim(.5,72.5)
for site,c in [(41,ORANGE),(2,PURPLE)]:pa.axvspan(site-.5,site+.5,color=c,alpha=.16,lw=0,zorder=0)
pa.plot(ids,reference,color=BLACK,lw=.7,ls=':',marker='o',markersize=2.4,markerfacecolor='white',markeredgewidth=.55,label='MD 参考',zorder=2)
pa.plot(ids,model,color=BLUE,lw=.9,marker='s',markersize=2.2,markeredgewidth=0,label='模型预测',zorder=3)
pa.set_ylim(-.15,4.25);pa.set_yticks([0,2,4]);pa.set_ylabel('联合 Hχ\n(nats)')
pa.set_xticks([1,10,20,30,40,50,60,72])
axis_labels=[];gln_texts=[]
for site in [2,41]:
    c=identity(site)
    label_x=float(site)
    pa.plot([site,site],[0,-.055],transform=pa.get_xaxis_transform(),color=c,lw=.8,clip_on=False,zorder=5)
    pa.plot([site,label_x],[-.195,-.24],transform=pa.get_xaxis_transform(),color=c,lw=.55,clip_on=False,zorder=5)
    y=-.30
    t=pa.text(label_x,y,'Q'+str(site),transform=pa.get_xaxis_transform(),ha='center',va='center',fontsize=7.5,color=c,clip_on=False)
    gln_texts.append(t)
    axis_labels.append(dict(residue_id=site,label='Q'+str(site),anchor_x_data=float(site),label_x_data=label_x,y_axes=y))
text(67,65,'泛素残基编号',fontsize=7.5,ha='center')
joint_legend_handles=[Line2D([0],[0],color=BLACK,ls=':',lw=.7,marker='o',markersize=3.2,markerfacecolor='white',markeredgewidth=.7),Line2D([0],[0],color=BLUE,lw=.9,marker='s',markersize=3,markeredgewidth=0)]
fig.legend(joint_legend_handles,['MD 参考','模型预测'],loc='lower left',bbox_to_anchor=(72/W,116/H),ncol=2,fontsize=7.5,handlelength=1.5,handletextpad=.45,columnspacing=.8,borderaxespad=0)

CORR=ROOT/'experiment/per_chi_correlations.csv'
correlations=pd.read_csv(CORR,float_precision='round_trip')
chosen=correlations[(correlations.predictor=='predicted_joint_Hchi') & correlations.metric.isin(['fine_entropy','low_Fstar_width'])].copy()
chosen['chi_number']=chosen.chi.astype(str).str.extract(r'(\d+)')[0].astype(int)
assert len(chosen)==8
assert sorted(chosen.chi_number.tolist())==[1,1,2,2,3,3,4,4]
expected_n={1:66,2:52,3:23,4:10}
assert all(int(r.n_defined)==expected_n[int(r.chi_number)] for _,r in chosen.iterrows())
assert np.isfinite(chosen.rho_AA_centered).all()
text(132,128,'逐χ MD占据关联',fontsize=8.5,weight='bold')
text(132,120,'同类型控制 Spearman ρ',fontsize=7.5)
bp=axis(134,79,44,31)
bp.set_ylim(3.5,-.5)
bp.set_xlim(-1,1);bp.set_xticks([-1,0,1]);bp.set_xlabel('Spearman ρ',fontsize=7.5)
bp.set_yticks(range(4),['χ'+str(i)+'  n='+str(expected_n[i]) for i in range(1,5)])
bp.tick_params(axis='y',length=0)
bp.axvline(0,color=GREY,lw=.6,ls=':',zorder=0)
for chi in range(1,5):
    for metric,marker,offset in [('fine_entropy','o',-.12),('low_Fstar_width','D',.12)]:
        row=chosen[(chosen.chi_number==chi)&(chosen.metric==metric)]
        assert len(row)==1
        r=row.iloc[0];value=float(r.rho_AA_centered);y=chi-1+offset
        if metric=='fine_entropy':bp.scatter(value,y,s=24,marker=marker,facecolors=BLUE,edgecolors=BLUE,linewidths=.7,zorder=3)
        else:bp.scatter(value,y,s=24,marker=marker,facecolors='white',edgecolors='#666666',linewidths=.8,zorder=3)
        points.append(dict(panel='a_chi_association',chi=chi,kind=metric,x=value,y=y,n_defined=int(r.n_defined)))
legend_handles=[Line2D([0],[0],color=BLUE,marker='o',lw=0,markersize=4,markerfacecolor=BLUE),Line2D([0],[0],color='#666666',marker='D',lw=0,markersize=4,markerfacecolor='white')]
fig.legend(legend_handles,['细箱占据熵','低F*宽度'],loc='lower left',bbox_to_anchor=(131/W,65/H),ncol=2,fontsize=7.5,handlelength=1,handletextpad=.45,columnspacing=.8,borderaxespad=0)

text(20,53,'谷氨酰胺中预测最低 Q41：Hχ = 2.000',color=ORANGE,fontsize=7.5)
text(100,53,'谷氨酰胺中预测最高 Q2：Hχ = 2.990',color=PURPLE,fontsize=7.5)
maxf=float(hist.F_star.max())
display_widths=[]
for qi,left in enumerate([20,74,128],1):
    ax=axis(left,12,43,33)
    letter='bcd'[qi-1]
    text(left-5,47,letter,fontsize=11,weight='bold')
    ax.set_title('χ'+str(qi)+'：MD 自由能投影',loc='left',pad=7,fontsize=8)
    for j,(site,c,style) in enumerate([('Q41',ORANGE,'-'),('Q2',PURPLE,'--')]):
        p=hist[(hist.residue==site)&(hist.chi=='chi'+str(qi))].sort_values('angle_center_degrees');assert len(p)==72
        ax.plot(p.angle_center_degrees,p.F_star,color=c,lw=1.3,ls=style,zorder=3)
        good=np.isfinite(p.F_star)&(p.F_star<=1)
        ax.broken_barh([(float(v),5.) for v in p.loc[good,'angle_lower_degrees']],(-.65+j*.3,.16),facecolors=c,edgecolors='none',zorder=3)
        width0=5*int(good.sum());display_widths.append(dict(residue=site,chi=qi,origin0_lowFstar_width_degrees=width0))
        ax.text(.03,maxf+1.15-j*.60,site+' 低值角区 '+str(width0)+'°',transform=ax.get_yaxis_transform(),color=c,fontsize=7.5,va='center')
        for _,r in p.iterrows():points.append(dict(panel=letter,residue_id=int(site[1:]),chi=qi,kind='MD_Fstar',x=r.angle_center_degrees,y=r.F_star,low_Fstar_le1=bool(np.isfinite(r.F_star) and r.F_star<=1)))
    ax.axhline(1,color=GREY,lw=.7,ls=':',zorder=0)
    ax.set_xlim(-180,180);ax.set_ylim(-.8,maxf+1.8);ax.set_xticks([-180,-60,60,180]);ax.set_xlabel('侧链角度 (°)',fontsize=7.5)
    if qi==1:ax.set_ylabel('MD F* (无量纲)',fontsize=7.5)
    else:ax.tick_params(labelleft=False)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[]
for t in texts:
    box=t.get_window_extent(renderer)
    if box.x0<-.5 or box.y0<-.5 or box.x1>fig.bbox.width+.5 or box.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
label_boxes=[t.get_window_extent(renderer) for t in gln_texts]
assert all(not a.overlaps(b) for i,a in enumerate(label_boxes) for b in label_boxes[i+1:])
numeric_boxes=[t.get_window_extent(renderer) for t in pa.get_xticklabels()]
assert all(not a.overlaps(b) for a in label_boxes for b in numeric_boxes)
assert len(set(item['y_axes'] for item in axis_labels))==1
assert min(t.get_fontsize() for t in texts)>=7.5
for ext in ['pdf','svg','png']:fig.savefig(OUT/(STEM+'.'+ext),dpi=300,bbox_inches=None)
assert min(Image.open(OUT/(STEM+'.png')).info['dpi'])>=299
assert not ET.parse(OUT/(STEM+'.svg')).findall('.//{http://www.w3.org/2000/svg}image')
pd.DataFrame(points).to_csv(OUT/(STEM+'.points.csv'),index=False,float_format='%.17g')
pd.DataFrame(display_widths).to_csv(OUT/(STEM+'.display_widths.csv'),index=False)
chosen.sort_values(['chi_number','metric']).to_csv(OUT/(STEM+'.all4chi_correlations.csv'),index=False,float_format='%.17g')
assert len(points)==572
caption='图7：侧链联合熵预测与逐角MD占据特征的对应。（a）左侧为泛素全部66个Hχ有效残基的联合熵预测与MD参考，仅标出Q2与Q41；右侧为模型联合Hχ与χ1–χ4逐角MD细箱占据熵、低自由能角区宽度的同类型控制Spearman相关。每个χ层内先按氨基酸类型分别减去预测值与MD指标的均值，再计算该层全部残差的秩相关；n为有该角的残基数，依次为66、52、23、10。逐角MD指标采用与b–d相同的501帧，并对四个箱原点取平均。（b–d）泛素六个谷氨酰胺中联合Hχ预测最低Q41与最高Q2的χ1、χ2及χ3逐角MD自由能F*=−ln(p̂/p̂max)，采用原点0°的5°箱，零计数箱断线。虚点线为F*=1，底部色条与标注为该箱格低自由能角区总宽度，±180°对应同一圆周接缝。联合Hχ参考为原五个连续1,000帧窗口的均值；曲线采用0–100 ns、200 ps间隔501帧，模型为第0折冻结完整模型。'
(OUT/(STEM+'.caption_zh.md')).write_text(caption+'\n',encoding='utf8')
audit={'status':'numerical_checks_passed_visual_pending','configured_canvas_mm':[W,H],'all66_valid_including3zero_Pro':True,'a_left_Gln_labels':axis_labels,'numeric_axis_ticks':[1,10,20,30,40,50,60,72],'single_Gln_label_row':True,'axis_labels_no_overlap':True,'all4_chi_included':[1,2,3,4],'all2_MD_metrics_included':['fine_entropy','low_Fstar_width'],'n_by_chi':expected_n,'association_axis_range':[-1,1],'association_not_per_angle_model_prediction':True,'panel_counts':pd.DataFrame(points).groupby('panel').size().to_dict(),'all_points':len(points),'no_per_angle_prediction_claim':True,'source_sha256':{str(p):sha(p) for p in [source,DATA/'all66_residue_summary.csv',CASE/'histograms_501frames.csv',OLD/'fig7_all_residues_Gln_and_MD_free_energy.statistics.json',CORR]},'output_sha256':{ext:sha(OUT/(STEM+'.'+ext)) for ext in ['png','pdf','svg']},'minimum_font_pt':min(t.get_fontsize() for t in texts),'text_outside_canvas':outside,'frozen_overall_metrics':stats,'visual_review':'pending','parent_figure_version':'v34'}
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k not in ['source_sha256','frozen_overall_metrics']},ensure_ascii=False))
