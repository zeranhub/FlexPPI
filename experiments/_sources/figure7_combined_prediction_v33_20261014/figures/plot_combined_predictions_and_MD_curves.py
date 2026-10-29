# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a left/right) complete66 Hchi profile and all6GLN paired predictions -> LineTrend/plot_trend.py and PairedBoxScatter -> param inherit
# (b-d) original frozen MD Fstar curves -> LineTrend/plot_trend.py -> param inherit
# RULE: no native-run panels; production assets have incompatible application-specific data structures.
# Retain all66 valid Hchi positions, all6GLN predictions, and original432MDcurvepoints; no smoothing or sampling.
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
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;DATA=ROOT.parent/'rich_angle_figure_v30_20261014/experiment'
if not DATA.exists():DATA=ROOT/'source'
OLD=ROOT.parent/'all_residues_figure_v29_20261014/figures'
CASE=ROOT.parent/'free_energy_curves_v27_20261014/case_curves'
if not (OLD/'fig7_all_residues_Gln_and_MD_free_energy.statistics.json').exists():OLD=ROOT/'source'
if not (CASE/'histograms_501frames.csv').exists():CASE=ROOT/'source'
STEM='fig7_combined_predictions_and_MD_curves'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
source=DATA/'ubiquitin_all66_chi_501frames.npz'
z=np.load(source,allow_pickle=False)
rows=pd.read_csv(DATA/'all66_residue_summary.csv')
stats=readj(OLD/'fig7_all_residues_Gln_and_MD_free_energy.statistics.json')
gln=rows[rows.amino_acid=='Q'].sort_values('predicted_Hchi_nats')
assert len(rows)==66 and len(gln)==6 and list(gln.author_residue_id)==[41,31,40,62,49,2]
assert np.array_equal(rows.residue_index,z['residue_index'])
np.testing.assert_allclose(rows.predicted_Hchi_nats,z['predicted_Hchi_nats'],rtol=0,atol=1e-14)
np.testing.assert_allclose(rows.frozen_MD_joint_Hchi_nats,z['frozen_MD_joint_Hchi_nats'],rtol=0,atol=1e-14)
rows['predicted_Hchi_nats']=z['predicted_Hchi_nats']
rows['frozen_MD_joint_Hchi_nats']=z['frozen_MD_joint_Hchi_nats']
gln=rows[rows.amino_acid=='Q'].sort_values('predicted_Hchi_nats')
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
for site in [2,31,40,41,49,62]:
    c=identity(site)
    label_x={40:37.,41:44.}.get(site,float(site))
    pa.plot([site,site],[0,-.055],transform=pa.get_xaxis_transform(),color=c,lw=.8,clip_on=False,zorder=5)
    pa.plot([site,label_x],[-.195,-.24],transform=pa.get_xaxis_transform(),color=c,lw=.55,clip_on=False,zorder=5)
    y=-.30
    t=pa.text(label_x,y,'Q'+str(site),transform=pa.get_xaxis_transform(),ha='center',va='center',fontsize=7.5,color=c,clip_on=False)
    gln_texts.append(t)
    axis_labels.append(dict(residue_id=site,label='Q'+str(site),anchor_x_data=float(site),label_x_data=label_x,y_axes=y))
text(67,65,'泛素残基编号',fontsize=7.5,ha='center')

text(132,128,'全部六个谷氨酰胺',fontsize=8.5,weight='bold')
text(132,120,'MAE 0.105 nats',fontsize=7.5)
bp=axis(134,79,25,31)
be=axis(161,79,17,31)
bp.set_title('联合 Hχ',loc='left',fontsize=8,pad=7)
be.set_title('预测−参考',fontsize=7.5,pad=7)
bp.set_ylim(5.5,-.5);be.set_ylim(5.5,-.5)
bp.set_xlim(1.65,3.15);bp.set_xticks([2,2.5,3]);bp.set_xlabel('Hχ (nats)',fontsize=7.5)
bp.set_yticks(range(6),gln.residue);bp.tick_params(axis='y',length=0)
be.axis('off');be.set_xlim(0,1)
for j,(_,r) in enumerate(gln.iterrows()):
    site=int(r.author_residue_id);c=identity(site)
    for ax in [bp,be]:
        if site in [41,2]:ax.axhspan(j-.48,j+.48,color=c,alpha=.08,lw=0,zorder=0)
    bp.plot([r.frozen_MD_joint_Hchi_nats,r.predicted_Hchi_nats],[j,j],color=c,lw=.8,zorder=1)
    bp.scatter(r.frozen_MD_joint_Hchi_nats,j,s=22,marker='o',facecolors='white',edgecolors=c,linewidths=.8,zorder=3)
    bp.scatter(r.predicted_Hchi_nats,j,s=21,marker='s',facecolors=BLUE,edgecolors=c,linewidths=.4,zorder=4)
    delta=r.predicted_Hchi_nats-r.frozen_MD_joint_Hchi_nats
    be.text(.5,j,format(delta,'+.3f'),ha='center',va='center',fontsize=7.5,color=CATEGORICAL[1] if delta>0 else BLUE)
    for kind,xv in [('model_prediction',r.predicted_Hchi_nats),('MD_reference',r.frozen_MD_joint_Hchi_nats)]:points.append(dict(panel='a_Gln',residue_id=site,kind=kind,x=xv,y=j,error=delta))
for lab,site in zip(bp.get_yticklabels(),gln.author_residue_id):lab.set_color(identity(site))
text(131,65,'○ MD 参考   ■ 模型预测',fontsize=7.5)

text(20,53,'六个位点中预测最低 Q41：Hχ = 2.000',color=ORANGE,fontsize=7.5)
text(100,53,'预测最高 Q2：Hχ = 2.990',color=PURPLE,fontsize=7.5)
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
gln.to_csv(OUT/(STEM+'.all6_prediction_order.csv'),index=False,float_format='%.17g')
assert len(points)==576
caption='图7：侧链联合熵预测与MD自由能取向范围的对应。（a）左侧为泛素全部66个Hχ有效残基的联合熵预测与MD参考，横轴标出全部六个谷氨酰胺；右侧将六个位点按冻结预测从低到高排列，对照联合熵参考、预测及带符号误差。空心圆为MD参考，实心方点为模型预测，Q41与Q2以橙色和紫色标识。（b–d）六个位点中预测最低Q41与最高Q2的χ1、χ2及χ3逐角MD自由能F*=−ln(p̂/p̂max)，采用原点0°的5°箱，零计数箱断线。虚点线为F*=1，底部色条与标注为该箱格低自由能角区总宽度，±180°对应同一圆周接缝。联合Hχ参考为原五个连续1,000帧窗口的均值；曲线采用0–100 ns、200 ps间隔501帧，模型为第0折冻结完整模型。'
(OUT/(STEM+'.caption_zh.md')).write_text(caption+'\n',encoding='utf8')
audit={'status':'numerical_checks_passed_visual_pending','dimensions_mm':[W,H],'all66_valid_including3zero_Pro':True,'a_deleted_error_and_width_tracks':True,'all6_Gln_axis_labels':axis_labels,'numeric_axis_ticks':[1,10,20,30,40,50,60,72],'single_Gln_label_row':True,'axis_labels_no_overlap':True,'occupancy_heatmaps_removed':True,'six_Gln_predictions_moved_to_a_right':True,'all6_order':[int(s) for s in gln.author_residue_id],'all6_Q40_largesterror_included':True,'panel_counts':pd.DataFrame(points).groupby('panel').size().to_dict(),'all_points':len(points),'no_per_angle_prediction_claim':True,'source_sha256':{str(p):sha(p) for p in [source,DATA/'all66_residue_summary.csv',CASE/'histograms_501frames.csv',OLD/'fig7_all_residues_Gln_and_MD_free_energy.statistics.json']},'output_sha256':{ext:sha(OUT/(STEM+'.'+ext)) for ext in ['png','pdf','svg']},'minimum_font_pt':min(t.get_fontsize() for t in texts),'text_outside_canvas':outside,'primary_metrics_unchanged':stats,'visual_review':'pending'}
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k not in ['source_sha256','primary_metrics_unchanged']},ensure_ascii=False))
