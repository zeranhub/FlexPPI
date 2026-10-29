from pathlib import Path
import re
ROOT=Path(__file__).resolve().parent
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
SKILL=Path(r'C:\Users\DELL\.codex\skills\academic-figure-skill')
head='''# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Paired scalar prediction/reference → PairedBoxScatter/plot_PairedBoxScatter1.R → param inherit
# (b–d) MD angular Fstar curves → existing plot_case_landscapes.py plus LineTrend/plot_trend.py → param inherit
# RULE: All paired values and all pooled/window histogram points retained. No native-run panel.
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'completed_results_20261012'/'plot_runtime'))
'''
baseline='\n'.join(re.search(r'```python\n(.*?)```',(SKILL/'references'/name).read_text(encoding='utf8'),re.S).group(1) for name in ['typography.md','color-palettes.md','export-specs.md'])
body='''
import hashlib,json,xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager
from PIL import Image
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;TASK=ROOT.parent
CASE=TASK/'energy_landscape_v21_20261012'/'case_analysis'
STEM='fig7_ubiquitin_prediction_and_MD_Fstar'
SOURCES={}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(name,expected):
    p=CASE/name;assert sha(p)==expected;SOURCES[str(p)]=sha(p);return pd.read_csv(p)
v=read('ubiquitin_bin_sensitivity.csv','27ba9b5c933e6789c739e8ec75da58e9d2ebcd3d391a796e3cfe8783b3949df9')
v=v[v.phase==0].set_index('residue').loc[['Q41','Q62']]
f=read('ubiquitin_chi1_histograms_and_Fstar.csv','454d72f79430fb2771912036c8258026c418926ec6b4cb165a54142c92c8f584')
assert len(f)==2592 and set(f.chi)=={1,2,3} and set(f.window)=={0,1,2,3,4,5}
assert f.groupby(['residue','window','chi'])['count'].sum().eq(f.groupby(['residue','window','chi']).frames.first()).all()
for group,part in f.groupby(['residue','window','chi']):
    p=part.probability.to_numpy();expected=np.full(72,np.nan);occupied=p>0
    expected[occupied]=-np.log(p[occupied]/p.max())
    np.testing.assert_allclose(expected,part.F_star,equal_nan=True,rtol=0,atol=1e-12)
font_manager.findfont('Microsoft YaHei',fallback_to_default=False)
mpl.rcParams.update({'font.sans-serif':['Microsoft YaHei','Arial'],'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.unicode_minus':False})
fig=plt.figure(figsize=(183/25.4,146/25.4))
gs=fig.add_gridspec(2,2,left=.085,right=.975,bottom=.115,top=.92,wspace=.31,hspace=.55)
axes=[fig.add_subplot(gs[i//2,i%2]) for i in range(4)]
plot=[];stats=[]
ax=axes[0]
ax.set_title('联合 Hχ 的模型预测与 MD 参考',loc='left',pad=10,fontsize=8.5)
ax.text(-.15,1.04,'a',transform=ax.transAxes,fontsize=11,weight='bold')
method_names=['MD 参考','模型预测'];method_colors=[CATEGORICAL[5],CATEGORICAL[0]]
for i,(site,row) in enumerate(v.iterrows()):
    truth=float(row.frozen_joint_h);prediction=float(row.pred_joint_h)
    assert abs(truth-row.joint_window_mean_h)<1e-12
    x=np.array([i-.16,i+.16]);y=np.array([truth,prediction])
    ax.plot(x,y,color='#A0A0A0',lw=.7,zorder=1)
    for j,method,color,marker in zip(range(2),method_names,method_colors,['o','s']):
        ax.scatter(x[j],y[j],s=37,color=color,marker=marker,edgecolors='white',linewidths=.55,zorder=3)
        ax.annotate(format(y[j],'.3f'),(x[j],y[j]),xytext=(-4 if j==0 else 4,9),textcoords='offset points',ha='right' if j==0 else 'left',va='bottom',fontsize=8,color=color)
        plot.append({'panel':'a','residue':site,'method':method,'x':x[j],'y':y[j]})
    ax.text(i,.39,'绝对误差 '+format(abs(prediction-truth),'.3f'),fontsize=7.5,ha='center',color=BLACK)
    stats.append({'residue':site,'MD_joint_Hchi_nats':truth,'model_joint_Hchi_nats':prediction,'absolute_error_nats':abs(prediction-truth),'MD_label_windows':5,'frames_per_window':1000,'fold':0})
ax.set_xticks([0,1],['Q41','Q62']);ax.set_xlim(-.5,1.5);ax.set_ylim(0,3.7);ax.set_ylabel('联合 Hχ (nats)')
ax.legend(handles=[Line2D([],[],marker='o',color=method_colors[0],ls='',markersize=4,label='MD 参考'),Line2D([],[],marker='s',color=method_colors[1],ls='',markersize=4,label='模型预测')],loc='upper left',bbox_to_anchor=(0,1.02),ncol=2,fontsize=7.5,handletextpad=.4,columnspacing=1)
diff_md=v.frozen_joint_h.iloc[1]-v.frozen_joint_h.iloc[0]
diff_model=v.pred_joint_h.iloc[1]-v.pred_joint_h.iloc[0]
ax.text(.5,.04,'Q62−Q41：MD '+format(diff_md,'.3f')+'；模型 '+format(diff_model,'.3f'),transform=ax.transAxes,ha='center',fontsize=7.5)
for qi,ax in enumerate(axes[1:],1):
    letter='bcd'[qi-1]
    ax.set_title('MD χ'+str(qi)+' 的经验投影自由能',loc='left',pad=10,fontsize=8.5)
    ax.text(-.15,1.04,letter,transform=ax.transAxes,fontsize=11,weight='bold')
    for site,color,style in [('Q41',CATEGORICAL[3],'-'),('Q62',CATEGORICAL[4],'--')]:
        for window in range(1,6):
            p=f[(f.residue==site)&(f.chi==qi)&(f.window==window)].sort_values('angle_center')
            ax.plot(p.angle_center,p.F_star,color=color,lw=.55,alpha=.22,ls=style,zorder=1)
        p=f[(f.residue==site)&(f.chi==qi)&(f.window==0)].sort_values('angle_center')
        ax.plot(p.angle_center,p.F_star,color=color,lw=1.4,ls=style,zorder=3)
    ax.set_xticks([-180,-60,60,180]);ax.set_xlim(-180,180);ax.set_ylim(0,8.5)
    ax.set_xlabel('MD 侧链角度 (°)');ax.set_ylabel('MD F* (无量纲)')
    ax.legend(handles=[Line2D([],[],color=CATEGORICAL[3],ls='-',lw=1.4,label='MD Q41'),Line2D([],[],color=CATEGORICAL[4],ls='--',lw=1.4,label='MD Q62')],loc='upper left',bbox_to_anchor=(0,1.02),ncol=2,fontsize=7.5,handlelength=1.5,columnspacing=1)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[]
for t in texts:
    box=t.get_window_extent(renderer)
    if box.x0<-.5 or box.y0<-.5 or box.x1>fig.bbox.width+.5 or box.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
assert min(t.get_fontsize() for t in texts)>=7.5
for ext in ['pdf','svg','png']:fig.savefig(OUT/(STEM+'.'+ext),dpi=300,bbox_inches=None)
im=Image.open(OUT/(STEM+'.png'));assert min(im.info['dpi'])>=299
svg=ET.parse(OUT/(STEM+'.svg'));ns={'s':'http://www.w3.org/2000/svg'}
assert not svg.findall('.//s:image',ns)
pd.DataFrame(plot).to_csv(OUT/(STEM+'.prediction_points.csv'),index=False,float_format='%.17g')
pd.DataFrame(stats).to_csv(OUT/(STEM+'.statistics.csv'),index=False,float_format='%.17g')
f.to_csv(OUT/(STEM+'.all_MD_curve_points.csv'),index=False,float_format='%.17g')
caption='图7：泛素同类型位点的联合侧链状态熵预测及其MD构象解释。（a）冻结第0折完整模型的联合Hχ均值与MD参考直接对照，灰色圆点为MD，蓝色方点为模型，点旁标注实际读数与绝对误差。MD参考为五个连续1,000帧窗口的联合熵均值；两位点均具有三个χ角。（b–d）同一轨迹前5,000帧（0–499.9 ns，100 ps间隔）的χ1、χ2和χ3经验投影自由能F*=−ln[p̂/p̂max]，以5°等宽箱估计密度，各曲线最低点为零。黄色实线为MD Q41，紫色虚线为MD Q62；粗线合并全部5,000帧，淡线为五个窗口，零计数箱断线。±180°两端对应同一圆周接缝。a检验联合熵的预测幅值与位点差异，b–d展开这些位点在MD中的实际角度占据。'
(OUT/(STEM+'.caption_zh.md')).write_text(caption+'\\n',encoding='utf8')
audit={'core_claim':'模型在该泛素病例中保留同类型位点Hχ幅值及高低差异，MD角度投影解释这些预测对应的占据形态。','evidence_chain':{'a':'模型vs参考直接对照，4actualpoints+2errors+MD/modeldelta','b–d':'onlyMD pooledand5windowFstarcurves explain differences in3angles'},'archetype':'quantitative_grid','statistics':'2 predefinedsameAAcasepositions,1frozenfold0model,MDreference5contiguouswindows(notindependentrepeats),5degreeFhistogramsfrom5000frames,pooledandall5windows,noCIorpvalue,descriptivecasevalidation; generalcohortproofinFigure8','source_sha256':SOURCES,'dimensions_mm':[183,146],'minimum_font_pt':min(t.get_fontsize() for t in texts),'all_curve_rows':2592,'prediction_points':4,'no_rasterized_svg':True,'text_outside_canvas':outside,'semantic_color_contract':{'a':'greyMDreferencecircle/bluemodelsquare; residueidentityonxaxis','b–d':'orangeQ41solid/purpleQ62dashed,allcurvesMD'},'asset_routing':'PairedBoxScatterconnectionlogicandoriginalMDcurvecodeparameterinheritance; no box summary for2sites; no cumulative transform from incompatibleLineTrend asset','all_pooled_and_window_data_preserved':True,'visual_review':'pending','output_sha256':{ext:sha(OUT/(STEM+'.'+ext)) for ext in ['png','pdf','svg']}}
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'prediction_points':4,'MDcurve_points':len(f),'MDdelta':diff_md,'modeldelta':diff_model,'outputs':audit['output_sha256']},ensure_ascii=False))
'''
(FIG/'plot_prediction_and_MD_Fstar.py').write_text(head+baseline+body,encoding='utf8')
print(FIG/'plot_prediction_and_MD_Fstar.py')
