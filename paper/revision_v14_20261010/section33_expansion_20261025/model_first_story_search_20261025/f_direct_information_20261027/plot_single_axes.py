# Academic Figure Skill Asset Confirmation (native Fig2 typography already verified)
# one physical-difference scatter -> frozen matched-record outputs -> param inherit
# All83 records; two physical units, no dual y-axis, no regression or joint ellipse.
from pathlib import Path
from hashlib import sha256
import json,shutil
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],
 'font.size':9,'axes.labelsize':9,'axes.titlesize':10,'xtick.labelsize':9,'ytick.labelsize':9,
 'legend.fontsize':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,
 'xtick.direction':'out','ytick.direction':'out','xtick.major.width':.6,'ytick.major.width':.6,
 'legend.frameon':False,'axes.unicode_minus':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.bbox':None})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

root=Path(__file__).resolve().parent;base=root.parents[2]
records=pd.read_csv(root/'record_means_original_RMSF_tol0.10.csv',float_precision='round_trip')
summary=pd.read_csv(root/'summary_all_prespecified_analyses.csv',float_precision='round_trip')
summary=summary[summary.matching_variable.eq('original_RMSF')&summary.tolerance_A.eq(.1)&summary.value_kind.eq('difference')].set_index('metric')
assert len(records)==83
x=records.original_RMSF_difference.to_numpy();y=records.fine_entropy_mean_per_chi_difference.to_numpy()
assert (y>0).sum()==82
W,H=120.,94.
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
fig.text(.5,1-2/H,'Matched backbone, distinct side-chain occupancy',ha='center',va='top',fontsize=10)
fig.text(.5,1-8/H,'Higher − lower predicted Hχ',ha='center',va='top',fontsize=9)
a=fig.add_axes([18/W,24/H,98/W,54/H])
a.axhline(0,color='#A8B1B7',lw=.6,ls='--',zorder=0)
a.axvline(0,color='#A8B1B7',lw=.6,ls='--',zorder=0)
a.scatter(x,y,s=13,color='#2166AC',alpha=.46,linewidth=.3,edgecolor='white',zorder=2)
sx=summary.loc['original_RMSF'];sy=summary.loc['fine_entropy_mean_per_chi']
a.errorbar(sx['mean'],sy['mean'],xerr=[[sx['mean']-sx.bootstrap95_low],[sx.bootstrap95_high-sx['mean']]],
 yerr=[[sy['mean']-sy.bootstrap95_low],[sy.bootstrap95_high-sy['mean']]],fmt='D',ms=4.6,
 color='#262626',mfc='#262626',mew=.7,capsize=2.,elinewidth=1.0,zorder=6)
a.set(xlim=(-.1,.1),ylim=(-.06,.35),xlabel='MD backbone RMSF difference (Å)',ylabel='MD side-chain entropy\ndifference (nats/χ)')
a.set_xticks([-.1,-.05,0,.05,.1]);a.set_yticks([0,.1,.2,.3])
a.tick_params(length=2.5,pad=2,width=.6);a.xaxis.labelpad=3;a.yaxis.labelpad=3
a.text(.035,.96,'82/83 entropy differences > 0',transform=a.transAxes,ha='left',va='top',fontsize=9)
handles=[Line2D([],[],ls='',marker='o',ms=3.5,color='#2166AC',alpha=.6,label='Records (n=83)'),Line2D([],[],ls='',marker='D',ms=4,color='#262626',label='Mean ± 95% CI')]
fig.legend(handles=handles,loc='center',bbox_to_anchor=(.55,8/H),ncol=2,handletextpad=.4,columnspacing=1.2,frameon=False)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
outside=[]
for t in fig.findobj(mpl.text.Text):
 if t.get_visible() and t.get_text():
  b=t.get_window_extent(renderer)
  if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
stem='figure3_f_v102_20261027_主链侧链差异合一'
for ext in ['png','pdf','svg']:
 p=root/f'f_single_axes_physical_differences.{ext}'
 fig.savefig(p,dpi=450,bbox_inches=fig.bbox_inches)
 shutil.copy2(p,base/(stem+'.'+ext))
plt.close(fig)
plotted=records[['sample_id','dataset','original_RMSF_difference','fine_entropy_mean_per_chi_difference']].copy()
plotted.to_csv(root/'single_axes_all83_points.csv',index=False)
caption='''同一记录内、同一种氨基酸中，将MD参考RMSF差不超过0.10Å的残基作不重复配对，再按模型Hχ预测将每对定为较低/较高两端。每点为一条记录，先等权氨基酸组再汇总；横纵轴均为高预测Hχ端减低预测Hχ端，分别保留Å和nats/χ的物理单位。83条记录均保留，82条的侧链熵差为正；整体MD参考RMSF差0.002239Å，侧链占据熵差0.157966nats/χ。黑菱形为等记录均值，横竖误差线分别为两个均值各自的来源分层bootstrap95%区间，不是二维联合区间。竖直虚线表示参考RMSF差0，水平虚线表示侧链占据熵差0。坐标轴不是两个相关系数，也未做归一化。RMSF近似匹配而非严格相同，不能将图解释为彻底排除主链因素或证明信息独立。原匹配协议、全部阈值敏感性、覆盖审计均沿用同目录冻结分析，无新增筛选、模型推理或MD。'''
(root/'single_axes_caption_zh.md').write_text(caption,encoding='utf-8')
(root/'single_axes_QA.json').write_text(json.dumps({'records':83,'x':'high minus low predicted-Hchi group MD reference RMSF, Angstrom','y':'same direction MD fine-bin side-chain entropy, nats/chi','same_matched_analysis_unchanged':True,'all_record_points':True,'no_dual_axis':True,'no_normalization_or_correlation_of_correlations':True,'zero_reference_lines':True,'marginal_CIs_not_joint':True,'outside_text':outside,'visual_review':'pending','size_mm':[W,H],'sha256':{ext:sha256((root/f'f_single_axes_physical_differences.{ext}').read_bytes()).hexdigest() for ext in ['png','pdf','svg']}},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'image':str(base/(stem+'.png')),'pdf':str(base/(stem+'.pdf')),'svg':str(base/(stem+'.svg')),'records':83,'mean_RMSF_difference_A':float(sx['mean']),'mean_sidechain_entropy_difference_nats_per_chi':float(sy['mean'])},ensure_ascii=False,indent=2))
