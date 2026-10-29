# Academic Figure Skill Asset Confirmation
# a-c: frozen MD Fstar curve code + LineTrend/plot_trend.py visual-parameter inheritance.
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
STEM='fig7_model_selected_three_MD_free_energy_curves'
font_manager.findfont('Microsoft YaHei',fallback_to_default=False)
mpl.rcParams.update({'font.sans-serif':['Microsoft YaHei','Arial'],'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.unicode_minus':False})
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
sources={}
def read(base,name):
    p=base/name;sources[str(p)]=sha(p);return pd.read_csv(p)
hist=read(CASE,'histograms_501frames.csv')
inventory=read(CASE,'selection_inventory_with_MD_outcomes.csv')
ranked=inventory.sort_values(['predicted_joint_Hchi_nats','residue_index'])
low=ranked.iloc[0];high=ranked.iloc[-1]
assert len(inventory)==6 and low.residue=='Q41' and high.residue=='Q2'
assert low.selection=='lowest_predicted_Hchi' and high.selection=='highest_predicted_Hchi'
assert len(hist)==432 and set(hist.residue)=={low.residue,high.residue} and hist.n_frames.eq(501).all()
for key,p in hist.groupby(['residue','chi']):
    assert p['count'].sum()==501
    occupied=p.probability>0
    np.testing.assert_allclose(p.F_star[occupied],-np.log(p.probability[occupied]/p.probability.max()),rtol=0,atol=1e-12)
    assert p.F_star[~occupied].isna().all()
styles=[(low.residue,CATEGORICAL[3],'-'),(high.residue,CATEGORICAL[4],'--')]
fig,axes=plt.subplots(1,3,figsize=(183/25.4,86/25.4),sharey=True)
fig.subplots_adjust(left=.065,right=.985,bottom=.225,top=.73,wspace=.22)
fig.text(.065,.965,'模型预测排序：六个谷氨酰胺位点 → 选择预测最低与最高位置',fontsize=8.5,weight='bold')
fig.legend(handles=[Line2D([],[],color=CATEGORICAL[3],lw=1.4,ls='-',label='模型低值 Q41：Hχ = '+format(low.predicted_joint_Hchi_nats,'.3f')),Line2D([],[],color=CATEGORICAL[4],lw=1.4,ls='--',label='模型高值 Q2：Hχ = '+format(high.predicted_joint_Hchi_nats,'.3f'))],loc='upper left',bbox_to_anchor=(.056,.90),ncol=2,fontsize=8,handlelength=1.6,columnspacing=2)
points=[];widths=[]
maxf=float(hist.F_star.max())
for qi,ax in enumerate(list(axes.flat)[:3],1):
    letter='abc'[qi-1]
    ax.set_title('χ'+str(qi)+'：MD 自由能投影',loc='left',pad=10,fontsize=8.5)
    ax.text(-.12,1.035,letter,transform=ax.transAxes,fontsize=11,weight='bold',va='bottom')
    for j,(site,color,style) in enumerate(styles):
        p=hist[(hist.residue==site)&(hist.chi=='chi'+str(qi))].sort_values('angle_center_degrees')
        assert len(p)==72
        ax.plot(p.angle_center_degrees,p.F_star,color=color,lw=1.3,ls=style,zorder=3)
        good=np.isfinite(p.F_star)&(p.F_star<=1)
        blocks=[(float(x),5.) for x in p.loc[good,'angle_lower_degrees']]
        ax.broken_barh(blocks,(-.65+j*.3,.16),facecolors=color,edgecolors='none',zorder=3)
        width=5*int(good.sum())
        widths.append({'residue':site,'chi':qi,'display_origin0_lowFstar_width_degrees':width,'four_origin_mean_width_degrees':float(inventory.set_index('residue').loc[site,'MD_width_501_chi'+str(qi)+'_deg'])})
        ax.text(.03,maxf+1.15-j*.60,site+' 低值角区 '+format(width,'.0f')+'°',transform=ax.get_yaxis_transform(),color=color,fontsize=7.5,va='center')
        for _,r in p.iterrows():points.append({'panel':letter,'residue':site,'chi':qi,'x':r.angle_center_degrees,'y':r.F_star,'low_Fstar_le1':bool(np.isfinite(r.F_star) and r.F_star<=1)})
    ax.axhline(1,color=GREY,lw=.7,ls=':',zorder=0)
    ax.set_xlim(-180,180);ax.set_ylim(-.8,maxf+1.8)
    ax.set_xticks([-180,-60,60,180]);ax.set_xlabel('侧链角度 (°)')
    if qi==1:ax.set_ylabel('MD F* (无量纲)')
fig.text(.065,.072,'底部色条：MD 中 F* ≤ 1 的实际 5°角区；黄色为模型低值位点，紫色为高值位点。',fontsize=7.5)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[]
for t in texts:
    b=t.get_window_extent(renderer)
    if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
assert min(t.get_fontsize() for t in texts)>=7.5
for ext in ['pdf','svg','png']:fig.savefig(OUT/(STEM+'.'+ext),dpi=300,bbox_inches=None)
assert min(Image.open(OUT/(STEM+'.png')).info['dpi'])>=299
svg=ET.parse(OUT/(STEM+'.svg'));assert not svg.findall('.//{http://www.w3.org/2000/svg}image')
pd.DataFrame(points).to_csv(OUT/(STEM+'.points.csv'),index=False,float_format='%.17g')
pd.DataFrame(widths).to_csv(OUT/(STEM+'.display_widths.csv'),index=False,float_format='%.17g')
assert len(points)==432
caption='图7：模型预先选择位点的逐角MD自由能验证。冻结第0折完整模型在泛素六个有效谷氨酰胺位点中按Hχ预测排序，选择最低的Q41与最高的Q2；黄色实线和紫色虚线分别对应模型低值与高值。a–c为所选位点在0–100 ns、200 ps间隔501帧中的χ1、χ2和χ3经验投影自由能F*=−ln(p̂/p̂max)，采用原点0°的5°箱，零计数箱断线。虚点线为F*=1；底部色条及图内标注为该箱格中满足F*≤1的全部角区总长度，±180°两端对应同一圆周接缝。'
(OUT/(STEM+'.caption_zh.md')).write_text(caption+'\n',encoding='utf8')
audit={'core_claim':'Frozen model predictions preselect sameAA extrema, MD angular curves verify the selected actual angle ranges','panels':{'a':'chi1 MD','b':'chi2 MD','c':'chi3 MD'},'dimensions_mm':[183,86],'all_points':432,'all_histogram_rows':432,'all_six_prediction_candidates_retained':True,'source_sha256':sources,'minimum_font_pt':min(t.get_fontsize() for t in texts),'text_outside_canvas':outside,'case_selection':'Original frozen model argminQ41/argmaxQ2 unchanged; noMDoutcome selection','removed_panel':'dremoved byuser; originalpopulationdataandstatisticalresultsretained','output_sha256':{ext:sha(OUT/(STEM+'.'+ext)) for ext in ['png','pdf','svg']},'visual_review':'pending'}
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'points':432,'model_selected_sites':[low.residue,high.residue],'outputs':audit['output_sha256']},ensure_ascii=False))
