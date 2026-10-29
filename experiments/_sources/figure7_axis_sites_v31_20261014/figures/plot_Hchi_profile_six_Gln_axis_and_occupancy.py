# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) complete66 Hchi sequence profile + six GLN axis labels -> LineTrend/plot_trend.py -> param inherit
# (b) six-site Hchi + chi1-chi3 occupancy -> PairedBoxScatter and heatmap/plot_composition.py -> param inherit
# (c-e) original frozen MD Fstar curves -> LineTrend/plot_trend.py -> param inherit
# RULE: no native-run panels; production assets have incompatible application-specific data structures.
# Retain all66 valid Hchi positions and all6GLN18angles; no smoothing, sampling or imputed missing angles.
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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from PIL import Image
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;DATA=ROOT.parent/'rich_angle_figure_v30_20261014/experiment'
if not DATA.exists():DATA=ROOT/'source'
OLD=ROOT.parent/'all_residues_figure_v29_20261014/figures'
CASE=ROOT.parent/'free_energy_curves_v27_20261014/case_curves'
if not (OLD/'fig7_all_residues_Gln_and_MD_free_energy.statistics.json').exists():OLD=ROOT/'source'
if not (CASE/'histograms_501frames.csv').exists():CASE=ROOT/'source'
STEM='fig7_Hchi_profile_six_Gln_axis_and_occupancy'
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
assert z['chi_defined'].sum()==151 and (~z['chi_defined']).sum()==113
hist=pd.read_csv(CASE/'histograms_501frames.csv');assert len(hist)==432
font_manager.findfont('Microsoft YaHei',fallback_to_default=False)
mpl.rcParams.update({'font.sans-serif':['Microsoft YaHei','Arial'],'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.unicode_minus':False})
BLUE=CATEGORICAL[0];ORANGE=CATEGORICAL[3];PURPLE=CATEGORICAL[4]
seq=LinearSegmentedColormap.from_list('occupancy',SEQUENTIAL)
seq.set_bad('#E1E1E1')
identity=lambda site: ORANGE if int(site)==41 else PURPLE if int(site)==2 else BLACK
W,H=183.,184.
fig=plt.figure(figsize=(W/25.4,H/25.4))
def axis(x,y,w,h,**kw):return fig.add_axes([x/W,y/H,w/W,h/H],**kw)
def text(x,y,s,**kw):return fig.text(x/W,y/H,s,**kw)
points=[]
text(3,177,'a',fontsize=11,weight='bold')
text(20,177,'完整泛素：联合熵的模型预测与 MD 参考',fontsize=8.5,weight='bold')
text(20,169,'全部 66 个有效位置   MAE 0.080 nats   ρ 0.979   类型均值校正 ρ 0.868',fontsize=7.5)
pa=axis(20,135,150,27)
ids=np.arange(1,73);model=np.full(72,np.nan);reference=model.copy()
for i,r in rows.iterrows():
    pos=int(r.author_residue_id)-1;model[pos]=r.predicted_Hchi_nats;reference[pos]=r.frozen_MD_joint_Hchi_nats
    for kind,y in [('model_prediction',model[pos]),('MD_reference',reference[pos])]:points.append(dict(panel='a_Hchi',residue_id=pos+1,kind=kind,x=pos+1,y=y))
pa.set_xlim(.5,72.5)
for site,c in [(41,ORANGE),(2,PURPLE)]:pa.axvspan(site-.5,site+.5,color=c,alpha=.16,lw=0,zorder=0)
pa.plot(ids,reference,color=BLACK,lw=.7,ls=':',marker='o',markersize=2.4,markerfacecolor='white',markeredgewidth=.55,label='MD 参考',zorder=2)
pa.plot(ids,model,color=BLUE,lw=.9,marker='s',markersize=2.2,markeredgewidth=0,label='模型预测',zorder=3)
pa.set_ylim(-.15,4.25);pa.set_yticks([0,2,4]);pa.set_ylabel('联合 Hχ\n(nats)')
pa.set_xticks([1,10,20,30,50,60,72])
pa.legend(loc='upper right',bbox_to_anchor=(1.0,1.34),ncol=2,fontsize=7.5,handlelength=1.8,columnspacing=1.4)
axis_labels=[]
for site in [2,31,40,41,49,62]:
    c=identity(site)
    pa.plot([site,site],[0,-.028],transform=pa.get_xaxis_transform(),color=c,lw=1.0,clip_on=False,zorder=5)
    y=-.48 if site==41 else -.33
    t=pa.text(site,y,'Q'+str(site),transform=pa.get_xaxis_transform(),ha='center',va='center',fontsize=7.5,color=c,clip_on=False)
    axis_labels.append(dict(residue_id=site,label='Q'+str(site),x_data=float(site),y_axes=y))
text(20,121,'泛素残基编号',fontsize=7.5)

text(3,116,'b',fontsize=11,weight='bold')
text(20,116,'全部六个谷氨酰胺：按模型预测从低到高排列，逐角展开实际占据',fontsize=8.5,weight='bold')
text(20,109,'联合 Hχ 平均绝对误差 0.105 nats',fontsize=7.5)
bp=axis(20,72,28,30)
be=axis(50,72,15,30)
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
    for kind,xv in [('model_prediction',r.predicted_Hchi_nats),('MD_reference',r.frozen_MD_joint_Hchi_nats)]:points.append(dict(panel='b_Hchi',residue_id=site,kind=kind,x=xv,y=j,error=delta))
for lab,site in zip(bp.get_yticklabels(),gln.author_residue_id):lab.set_color(identity(site))
gi=[int(np.flatnonzero(z['author_residue_id'].astype(int)==int(s))[0]) for s in gln.author_residue_id]
pmax=float(np.ceil(np.max(z['probability_origin0'][gi,:3,:])*20)/20)
for q,xstart in enumerate([70,104,138]):
    ax=axis(xstart,72,30,30)
    mat=z['probability_origin0'][gi,q,:]
    assert mat.shape==(6,72) and np.isfinite(mat).all() and np.allclose(mat.sum(axis=1),1,atol=1e-14)
    hm=ax.pcolormesh(np.arange(-180,181,5),np.arange(-.5,6.5),mat,cmap=seq,vmin=0,vmax=pmax,shading='flat',linewidth=0,rasterized=False)
    ax.set_ylim(5.5,-.5);ax.set_xlim(-180,180);ax.set_yticks([]);ax.set_xticks([-180,0,180]);ax.set_title('MD χ'+str(q+1)+' 占据',loc='left',fontsize=8,pad=7)
    labels=ax.get_xticklabels();labels[0].set_ha('left');labels[-1].set_ha('right')
    for sp in ax.spines.values():sp.set_visible(False)
    ax.tick_params(length=2)
    for j,i in enumerate(gi):
        site=int(z['author_residue_id'][i]);c=identity(site)
        if site in [41,2]:ax.add_patch(Rectangle((-180,j-.5),360,1,fill=False,edgecolor=c,lw=.65,zorder=4))
        low=np.isfinite(z['Fstar_origin0'][i,q])&(z['Fstar_origin0'][i,q]<=1)
        ax.broken_barh([(float(-180+5*k),5.) for k in np.flatnonzero(low)],(j+.34,.095),facecolors=c,edgecolors='none',zorder=5)
        for k in range(72):points.append(dict(panel='b_probability',residue_id=site,chi=q+1,kind='MD_probability',x=-177.5+5*k,y=mat[j,k],count=z['counts_origin0'][i,q,k],low_Fstar_le1=bool(low[k])))
bcax=axis(171,72,2.2,30);bc=fig.colorbar(hm,cax=bcax,ticks=[0,pmax/2,pmax]);bc.outline.set_visible(False);bc.ax.tick_params(labelsize=7.5,length=2,width=.6)
bc.solids.set_rasterized(False)
text(172,108,'占据率',ha='center',fontsize=7.5)
text(119,64,'侧链角度 (°)',ha='center',fontsize=7.5)
text(20,60,'○ MD 参考   ■ 模型预测；占据热图底部细线标出 F* ≤ 1 的实际角区。',fontsize=7.5)

text(20,53,'六个位点中预测最低 Q41：Hχ = 2.000',color=ORANGE,fontsize=7.5)
text(100,53,'预测最高 Q2：Hχ = 2.990',color=PURPLE,fontsize=7.5)
maxf=float(hist.F_star.max())
display_widths=[]
for qi,left in enumerate([20,74,128],1):
    ax=axis(left,12,43,33)
    letter='cde'[qi-1]
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
assert min(t.get_fontsize() for t in texts)>=7.5
for ext in ['pdf','svg','png']:fig.savefig(OUT/(STEM+'.'+ext),dpi=300,bbox_inches=None)
assert min(Image.open(OUT/(STEM+'.png')).info['dpi'])>=299
assert not ET.parse(OUT/(STEM+'.svg')).findall('.//{http://www.w3.org/2000/svg}image')
pd.DataFrame(points).to_csv(OUT/(STEM+'.points.csv'),index=False,float_format='%.17g')
pd.DataFrame(display_widths).to_csv(OUT/(STEM+'.display_widths.csv'),index=False)
gln.to_csv(OUT/(STEM+'.all6_prediction_order.csv'),index=False,float_format='%.17g')
assert len(points)==1872
caption='图7：侧链联合熵预测与完整逐角占据的对应。（a）泛素全部66个Hχ有效残基的联合熵预测与MD参考，按残基编号展开，横轴标出全部六个谷氨酰胺位点。（b）全部六个谷氨酰胺按预测从低到高排列，每行对齐联合熵参考、预测与误差，以及原点0°的χ1–χ3完整5°角度占据。蓝色深浅表示每箱帧占比，三个角共用色标，底部细线标出F*≤1的角区；橙色和紫色行框标识Q41、Q2。（c–e）这两个预测端点的逐角MD自由能F*=−ln(p̂/p̂max)，采用原点0°的5°箱，零计数箱断线。虚点线为F*=1，底部色条与标注为该箱格低自由能角区总宽度，±180°对应同一圆周接缝。联合Hχ参考为原五个连续1,000帧窗口的均值；角度占据、宽度及曲线统一采用0–100 ns、200 ps间隔501帧，模型为第0折冻结完整模型。'
(OUT/(STEM+'.caption_zh.md')).write_text(caption+'\n',encoding='utf8')
audit={'status':'numerical_checks_passed_visual_pending','dimensions_mm':[W,H],'all66_valid_including3zero_Pro':True,'a_deleted_error_and_width_tracks':True,'all6_Gln_axis_labels':axis_labels,'all6_Gln_all3angles_full72bins':True,'all6_order':[int(s) for s in gln.author_residue_id],'all6_Q40_largesterror_included':True,'panel_counts':pd.DataFrame(points).groupby('panel').size().to_dict(),'all_points':len(points),'MD_probability_color_scale':[0,pmax],'no_per_angle_prediction_claim':True,'source_sha256':{str(p):sha(p) for p in [source,DATA/'all66_residue_summary.csv',CASE/'histograms_501frames.csv',OLD/'fig7_all_residues_Gln_and_MD_free_energy.statistics.json']},'output_sha256':{ext:sha(OUT/(STEM+'.'+ext)) for ext in ['png','pdf','svg']},'minimum_font_pt':min(t.get_fontsize() for t in texts),'text_outside_canvas':outside,'primary_metrics_unchanged':stats,'visual_review':'pending'}
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k not in ['source_sha256','primary_metrics_unchanged']},ensure_ascii=False))
