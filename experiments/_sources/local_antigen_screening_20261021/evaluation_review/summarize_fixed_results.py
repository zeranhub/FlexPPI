"""Report fixed results without altering any analysis rule or model selection."""
from pathlib import Path
import csv,json,hashlib
HERE=Path(__file__).resolve().parent;RESULTS=HERE/'results_fixed_protocol'
def rows(path):
    with Path(path).open(encoding='utf-8-sig') as h:return list(csv.DictReader(h))
checks=json.loads((HERE/'independent_evaluation_review.json').read_text(encoding='utf8'))
assert checks['status']=='passed_independent_numerical_and_leakage_review'
data=rows(RESULTS/'candidate_level_metrics.csv')
selected={(r['system'],r['model']):r for r in data if r['analysis']=='within_system_leave_position_out' and r['label_variant']=='reported_hybrid' and r['template_pdb']=={'nyeso':'2BNQ','tax':'1AO7'}[r['system']]}
names={
 'zero_change':'全部预测不变','sequence_only':'序列','static_only':'静态结构',
 'sequence_plus_static':'序列＋静态','sequence_plus_static_plus_R':'序列＋静态＋RMSF',
 'sequence_plus_static_plus_H':'序列＋静态＋Hχ','sequence_plus_static_plus_RH':'序列＋静态＋RMSF＋Hχ'}
table=[]
for model,name in names.items():
    item=dict(model=model,method_zh=name)
    for system in ['nyeso','tax']:
        r=selected[system,model]
        item.update({system+'_n':int(r['n_candidates']),system+'_MAE_delta_pKD':float(r['MAE_delta_pKD']),
            system+'_Spearman_rho':float(r['Spearman_rho']) if r['Spearman_rho'] else None,
            system+'_top3_overlap':int(r['top_k_overlap'])})
    table.append(item)
with (HERE/'primary_results_7models.csv').open('w',encoding='utf-8-sig',newline='') as h:
    w=csv.DictWriter(h,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
md='''# 本地抗原筛选pilot：固定评价结果

主检验没有支持“RMSF与Hχ联合读数能够稳定改善抗原候选筛选”。每个体系均为13种单点替换，按六个替换位置留出，使用提前固定的reported-hybrid SPR标签、1G4的2BNQ与A6的1AO7模板。与序列加静态结构相比，联合两项后1G4误差从0.4324增至0.4488 pKD，A6从0.5985增至0.6474 pKD，分别增加约3.8%和8.2%。这里的数值是候选的亲和力变化误差，数值越小越好，不能把它称为残基动态预测误差。完整七种方法如下；ρ为候选排序相关，Top3为实验前三名与预测前三名重合个数。

| 方法 | 1G4 MAE | 1G4 ρ | 1G4 Top3 | A6 MAE | A6 ρ | A6 Top3 |
|---|---:|---:|---:|---:|---:|---:|
'''
for item in table:
    rho=lambda x:'—' if x is None else f'{x:.3f}'
    ntop='无排序' if item['model']=='zero_change' else f"{item['nyeso_top3_overlap']}/3"
    atop='无排序' if item['model']=='zero_change' else f"{item['tax_top3_overlap']}/3"
    md+=f"| {item['method_zh']} | {item['nyeso_MAE_delta_pKD']:.4f} | {rho(item['nyeso_Spearman_rho'])} | {ntop} | {item['tax_MAE_delta_pKD']:.4f} | {rho(item['tax_Spearman_rho'])} | {atop} |\n"
md+='\n全部预测不变产生全体并列，不能提供候选排序。原始CSV按事先规定的编号顺序保留Top3记账值，但它没有筛选含义，因此此表标为“无排序”。\n'
md+='''
1G4有一个有限的正向信号：在同一序列加静态基线上仅加入Hχ，MAE从0.4324降至0.4092，排序相关从0.566升至0.588，前三名重合由0/3升至1/3。但RMSF单项使误差增至0.4751，同时加入两项也不及Hχ单项。这个Hχ信号在三种SPR拟合版本的2BNQ主模板下均为同方向，却未在2BNR敏感性模板中重现；逐study评分也有差异，Pettmann子集改善，而Aleksic子集MAE未改善。因此，不能把它改写成稳健的侧链熵筛选能力，更不能据此说两项读数具有互补收益。A6中序列单独的MAE为0.5299，优于序列加静态的0.5985；加Hχ、RMSF或两项均未改善主模板结果，表明增加特征本身没有形成可重复的候选信息增量。

完整敏感性没有挽救主结论。五模板乘三种SPR拟合版本的15个体系内比较中，联合RMSF/Hχ的MAE有13个比序列加静态更差。两个正值都在A6的1QSE敏感性模板，reported-hybrid仅改善0.0098 pKD、all-fitted仅改善0.0028 pKD，而且联合仍不及其中单项。三种标签版本下，提前定义的“两体系均改善”和“联合又优于两个单项”判断全部为否。不能拿这两个敏感性结果替换主模板，也不能只汇报1G4的Hχ单项。

两方向跨体系检验同样没有一致联合增量。A6训练后预测1G4时，序列加静态MAE为0.3450，联合为0.3475，且排序相关从0.784降至0.659；1G4训练后预测A6时，MAE从0.5433微降至0.5406，但排序相关从0.401降至0.302。A6两种双点替换的主模板附加检验中，联合MAE为1.0340，略低于基线1.0861，但只有两个候选，不能作为一般筛选能力或位置外推的证据；其“Top3重合2/2”只是候选总数为二的必然结果，不具有筛选判别意义。

这些结果限定的是这一次预先固定的实际用法：历史2431联合模型，在固定实验骨架的候选改造输入上，用8 Å环境的两个有符号均值作为亲和力回归特征。基础模型未针对SPR训练，这次也没有采样新骨架或测量游离到结合的熵差。因此，当前结果不支持把这一简单组合直接作为已验证的抗原优化工具；它也不能反过来证明模型在全部残基动态预测任务上无效。可能的后续问题包括局部相反变化被均值抵消、结合前状态缺失或亲和力机制未被两个动态均值覆盖，但这些是新的待检验假设，不能在看完结果后更换特征、模板或候选并沿用本次预先固定验证的名义。

计算与防泄漏核验通过，并不意味着科学收益检验通过。正式评价一次产出2037条分析记录，来源仍为30种生物候选、其中26种主单点变体及两个TCR体系；多模板、拟合版本和模型比较没有增加独立实验数。独立程序没有导入正式评价函数，改用增广最小二乘重建全部630次非零ridge拟合，逐个重算2037个预测及210组MAE、排序相关和Top3；最大预测差为4.7×10⁻¹⁵。各折训练候选、同位置隔离、源体系隔离、仅训练候选计算的RMS尺度、固定λ=1与零截距全部核对一致。协议、U及评价代码SHA没有改变，没有按成绩删候选、改标签、换模板或调整规则。本轮使用历史2431权重，并保留此前已说明的一条间接同家族组件训练记录；没有更改论文。
'''
(HERE/'fixed_results_interpretation_zh.md').write_text(md,encoding='utf8')
print(json.dumps(dict(status='complete_report_without_analysis_changes',main_rows=len(table),independent_review='passed',joint_benefit_primary=False),ensure_ascii=False))
