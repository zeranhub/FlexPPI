from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(r'E:\Project\FlexPPI\experiments\_sources\antigen_incremental_value_20261028')
allmetrics=[];macros=[];permutations=[]
for variant in ['full','intrachain']:
 for version in ['primary','secondary']:
  P=R/'evaluation'/f'{variant}_{version}'
  assert (P/'evaluation_audit.json').is_file()
  m=pd.read_csv(P/'metrics_by_system_and_fold.csv');m['variant']=variant;m['template_version']=version;allmetrics.append(m)
  a=pd.read_csv(P/'macro_metrics.csv');a['variant']=variant;a['template_version']=version;macros.append(a)
  p=pd.read_csv(P/'shuffled_dynamic_controls_all_seeds.csv');p['variant']=variant;p['template_version']=version;permutations.append(p)
metrics=pd.concat(allmetrics,ignore_index=True);macro=pd.concat(macros,ignore_index=True)
metrics.to_csv(R/'all_settings_metrics.csv',index=False);macro.to_csv(R/'all_settings_macro_metrics.csv',index=False)
models=['B0_binding','B1_binding_static','B2_binding_static_dynamics']
main=metrics[(metrics.variant=='full')&(metrics.template_version=='primary')&(metrics.split=='leave_system_out')&metrics.model.isin(models)].copy()
main.to_csv(R/'primary_three_flow_results.csv',index=False)
primary=macro[(macro.variant=='full')&(macro.template_version=='primary')&(macro.split=='leave_system_out')].set_index('model')
secondary=macro[(macro.variant=='full')&(macro.template_version=='secondary')&(macro.split=='leave_system_out')].set_index('model')
def get(frame,model,col):return float(frame.loc[model,col])
null=pd.concat(permutations,ignore_index=True).groupby(['variant','template_version','seed'])[['top10_hits','average_precision']].mean().reset_index()
nulls=[]
for (v,t),g in null.groupby(['variant','template_version']):
 actual=macro[(macro.variant==v)&(macro.template_version==t)&(macro.split=='leave_system_out')].set_index('model')
 b2=actual.loc['B2_binding_static_dynamics']
 nulls.append({'variant':v,'template_version':t,'n_seeds':len(g),'actual_top10_mean':b2.top10_hits,'shuffled_top10_mean_median':g.top10_hits.median(),
  'fraction_shuffled_top10_ge_actual':float((g.top10_hits>=b2.top10_hits).mean()),'actual_average_precision':b2.average_precision,
  'shuffled_AP_median':g.average_precision.median(),'fraction_shuffled_AP_ge_actual':float((g.average_precision>=b2.average_precision).mean()),
  'interpretation':'Descriptive feature-permutation control, not an IID biological p-value.'})
pd.DataFrame(nulls).to_csv(R/'shuffled_control_summary.csv',index=False)

fig,axes=plt.subplots(1,2,figsize=(10.5,3.7),layout='constrained')
colors=['#666666','#0072B2','#D55E00'];labels=['Binding scores','+ structure/presentation','+ FlexPPI features']
x=np.arange(3);width=.24
for i,(model,color,label) in enumerate(zip(models,colors,labels)):
 values=[int(main[(main.fold==s)&(main.model==model)].top10_hits.iloc[0]) for s in ['nyeso','ebv','mage']]
 bars=axes[0].bar(x+(i-1)*width,values,width,color=color,label=label)
 axes[0].bar_label(bars,fontsize=9,padding=2)
axes[0].set_xticks(x,['NYESO / 1G4','EBV / TK3','MAGE']);axes[0].set_ylim(0,10.8);axes[0].set_ylabel('Known functional positives among top 10')
axes[0].set_title('Primary templates, held-out TCR systems',fontsize=11)
vals=[get(primary,m,'average_precision') for m in models]
bars=axes[1].bar(np.arange(3),vals,color=colors,width=.65);axes[1].bar_label(bars,labels=[f'{v:.3f}' for v in vals],padding=3,fontsize=10)
axes[1].set_xticks(np.arange(3),['Binding','+ structure /\npresentation','+ FlexPPI']);axes[1].set_ylim(0,1);axes[1].set_ylabel('Macro average precision');axes[1].set_title('Same fixed downstream classifier',fontsize=11)
for ax in axes:
 ax.spines[['top','right']].set_visible(False);ax.set_axisbelow(True);ax.yaxis.grid(True,alpha=.18)
fig.legend(*axes[0].get_legend_handles_labels(),frameon=False,fontsize=9,loc='outside lower center',ncol=3)
fig.savefig(R/'primary_screening_comparison.png',dpi=180);plt.close(fig)

base=get(primary,'B1_binding_static','top10_hits')*3;dyn=get(primary,'B2_binding_static_dynamics','top10_hits')*3
report='''# FlexPPI动态表征是否提高抗原候选筛选：实际增量试验

结论：本轮预先固定的区域动态表征没有带来稳定筛选增益，不能据此把FlexPPI写成已经验证的抗原筛选工具。测试的是已发表肽设计面板上的功能候选再排序，不是直接预测亲和力，也不是本轮新做了细胞实验。论文正文未修改。

## 数据与三条流程

使用HERMES研究公开All-designs.xlsx的全部279条设计肽：NYESO/1G4、EBV/TK3、MAGE/工程TCR各93条。作者报告的功能阳性分别为30、13、19条。仅排除每体系的两个WT控制及无肽刺激控制；未按模型分数、响应强度或是否成功删除设计。使用原作者resp/resp_el4标签，未人为设置KD阈值或把未响应称为不结合。三次测量、两个模板和逐残基输出都不增加生物候选数。

B0使用官方NetTCR-2.2统一pan模型与ERGO-II的McPAS、VDJdb版本，共三个结合分数，再以固定L2逻辑回归校准为本面板的功能排序。原始三个结合分数的排序也完整保留。B1在B0上加入作者公开的TCRdock PAE、AF3界面PAE、两种NetMHCpan呈递指示量及与WT的最少替换位数；它是“公开结构预测与呈递特征基线”，不是从本次固定骨架直接计算的纯静态几何控制。B2再加入八个提前固定的FlexPPI特征：TCRα/β第三识别环、α/β第一加第二识别环，各自平均RMSF和侧链Hχ相对参考肽的变化。只比较身份不变的受体残基，侧链熵排除A/G。没有假设柔性越低越好，关系由训练集上的同一种校准器学习。

每次完整留出一个TCR体系的93条候选，用另两个体系的186条训练。标准化、缺失处理和系数只由训练集得到，固定C=1，不根据结果换参数或区域。基础FlexPPI及结合预测器均未重训；训练的是明确列出的下游校准器。主结构使用2BNQ/3MV7/5BRZ，备用结构使用2BNR/4PRP/5BS0；两种FlexPPI分支共生成1128组结构–模型输出，对应564种结构输入。所有候选采用同一参考骨架替换肽身份，主链不移动，替换位旧侧链删除。该测试衡量假设结合姿态下的模型表征用途，没有得到真实突变体平衡构象，也不能判断复合物一定形成。

## 预定义主结果

| 留出的体系 | B0：结合分数 | B1：再加结构与呈递特征 | B2：再加FlexPPI特征 |
|---|---:|---:|---:|
'''
for system,label in [('nyeso','NYESO / 1G4'),('ebv','EBV / TK3'),('mage','MAGE')]:
 hits=[int(main[(main.fold==system)&(main.model==m)].top10_hits.iloc[0]) for m in models]
 report+=f'| {label} | {hits[0]}/10 | {hits[1]}/10 | {hits[2]}/10 |\n'
report+=f'''
每套体系只选前10条，总预算30条。B1选中{base:.0f}条已知功能阳性，B2选中{dyn:.0f}条；NYESO有小幅改善，EBV和MAGE下降。按体系平均的平均精度AP由{get(primary,'B1_binding_static','average_precision'):.3f}降至{get(primary,'B2_binding_static_dynamics','average_precision'):.3f}。AP同时考虑整份候选名单的阳性排序，越高越好；它不是预测亲和力的相关系数。三个独立体系不支持把279个候选当作279次独立跨体系外推检验，本报告不据此给出很小的生物学显著性p值。

![主结果](primary_screening_comparison.png)

## 对照与敏感性

显式加入20种肽氨基酸比例后，B1_AA的30条选择中有{get(primary,'B1_AA','top10_hits')*3:.0f}条阳性，加入动态后为{get(primary,'B2_AA','top10_hits')*3:.0f}条；AP从{get(primary,'B1_AA','average_precision'):.3f}变为{get(primary,'B2_AA','average_precision'):.3f}。加入逐位置氨基酸编码的较强序列对照后，前10命中从平均{get(primary,'B1_position_AA','top10_hits'):.3f}变为{get(primary,'B2_position_AA','top10_hits'):.3f}，AP则从{get(primary,'B1_position_AA','average_precision'):.3f}小幅升至{get(primary,'B2_position_AA','average_precision'):.3f}，指标方向不一致，没有得到稳定的额外收益。序列控制在下游比较指标生成前追加，原始主协议未被替换。

备用模板下，B1与B2的平均前10命中分别为{get(secondary,'B1_binding_static','top10_hits'):.3f}和{get(secondary,'B2_binding_static_dynamics','top10_hits'):.3f}，AP分别为{get(secondary,'B1_binding_static','average_precision'):.3f}和{get(secondary,'B2_binding_static_dynamics','average_precision'):.3f}。完整逐体系及全部模型结果见all_settings_metrics.csv，不用备用结果替换主结果。

链内分支中，改变肽链身份不会改变另两条TCR链的输出，八个受体相对参考特征为零，B1与B2的排序因此相同；这是预期的实现核查，不是新增筛选能力。完整模型的跨链响应确实非零，但“模型产生响应”与“响应帮助选中有效抗原”是不同问题。

同维度打乱对照在每个体系内随机置换整组八维动态特征，共200个固定种子，全部保留。它们仅描述增加相同自由度的无配对特征可能产生什么结果，不能把200个种子当作200次生物实验或正式独立样本p值。主结构实际B2没有表现出优于大多数此类打乱结果的稳定优势，详细数值见shuffled_control_summary.csv。即使优于打乱，也仍需优于B1才能证明实用增量。

次分析按肽序列Hamming距离不超过2的连通分量成组留出。NYESO某训练折只剩一类标签，按预先规则标记不能估计，没有改成随机逐行拆分。可估计的EBV、MAGE在主模板上前10命中分别为4→5、6→5，合计10→10，AP宏平均0.479→0.435；备用模板分别为4→5、6→6，合计10→11，但AP宏平均0.479→0.433下降。这一小幅前10增益也保留，不能因主结论为负而省略；它仍未与整体排序和跨体系检验形成一致改善，不能用这两个体系的次分析冒充三个体系都通过验证。

## 参考、训练覆盖与数值复核

完整训练池3063记录/4973组成链的补充局部与全局比对，没有发现输入长链达到既定阈值的fold0训练命中。但这是透明的Biopython替代核查，不是原MMseqs搜索的严格重跑；短肽近缘和所有基础模型近同源接触没有被完全排除。NetTCR发布训练表没有检出这三组输入的精确配对CDR3，因而按冻结规则保留20模型均值；ERGO发布训练表没有检出这些设计的精确αβ–肽或β–肽配对。它们不等于完全排除同源预训练接触，故本报告不称严格零样本验证。

ERGO所需V/J基因在本轮没有可靠逐项确认，统一保留缺省输入，未猜测基因名；这可能影响它的绝对表现。NetTCR使用经ANARCII核验的六个识别环。基础模型原本预测结合特异性，而本轮终点是细胞功能，因此不能把B0较低的功能排序成绩写成“主流模型无法预测结合”，也不能据此宣称它们普遍不如FlexPPI。

所有候选序列、PDB主链、模型权重SHA、逐残基NPZ、区域映射、每次拟合的训练/测试成员和完整排名均保留。独立代理用另一套逐候选统计程序复算主结果，top10、AP和AUROC一致，未发现体系交叉或标签进入特征；根代理另用独立加权逻辑损失复算校准，并核查所有NPZ与PDB。数值差、优化停止差异和任何top10边界差异见audits/final_numerical_verification.json，不隐去。

另做24次重复前向：同线程下与原输出逐值一致；跨2/4线程的最大单值差约0.0000258，八个区域相对特征的最大差约0.000000954。首次跨线程检查的1e-5阈值曾被超过，初始日志保留，未宣称不同线程逐比特相等；后续将同线程复算与跨线程浮点变化分开记录，见audits/repeated_forward_check.json。

原文补充材料说明实验细胞TCR使用鼠源恒定区，而这里使用公开晶体构建体；因此这是识别区和受体背景近似下的表征试验，不是完整细胞受体动态模拟。原作者部分设计模板来自TCRdock整理坐标，本轮统一使用直接下载的RCSB坐标，也不宣称逐项复现其模型。B1的PAE来自作者的另一套候选结构预测，而FlexPPI使用固定参考骨架；这进一步限制了真实动态机制归因。

## 这轮结果意味着什么

用户的顾虑有实测支持：不能仅靠“柔性参与识别”，就把预测柔性转换成有效筛选能力。总体未稳定改善，当前就没有依据声称筛选能力已被验证。本轮并未证明所有局部动态模式与识别都无关，也没有否定FlexPPI预测动态性质本身；它只表明当前固定骨架、八个预设区域特征和下游校准方式未显示稳定增量。NYESO单项的8到9不能盖过另外两套的下降，也不能重新挑区域或参数把同一实验改写为成功。

对当前论文，适合保留的主张是快速预测残基动态并为识别差异提出可检验的动态假说；抗原筛选应留作尚待验证的潜在应用。已有10月21日旧模型SPR试验也曾显示界面动态均值组合没有稳定改善，但它与本轮模型、数据和终点不同，只作已知背景，不冒称新的独立重复。本轮没有继续搜索新特征直到获得正结果，也未修改论文正文。

数据来源：[HERMES/PNAS 2025](https://doi.org/10.1073/pnas.2504783122)，[作者公开完整候选表](https://github.com/StatPhysBio/tcr_antigen_design/blob/main/peptide_designs/All-designs.xlsx)。代码、协议、源表SHA、完整结果和限制均在本目录。
'''
(R/'实际试跑结果与筛选能力结论.md').write_text(report,encoding='utf-8')
print(main[['fold','model','top10_hits','average_precision']].to_string(index=False))
print('Primary B1/B2 totals',base,dyn)
print('Secondary B1/B2 macro',secondary.loc[['B1_binding_static','B2_binding_static_dynamics'],['top10_hits','average_precision']].to_string())
