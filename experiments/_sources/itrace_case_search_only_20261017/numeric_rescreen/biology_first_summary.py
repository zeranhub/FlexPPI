from pathlib import Path
import sys
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/'label_sensitivity/runtime'))
import pandas as pd

bio=pd.read_csv(P.parent/'biology_first/biology_first_ITrace_candidate_series.csv')
v=pd.read_csv(P/'all_models_complete_region_comparison.csv')
rows=[]
for case in bio.itertuples():
    for pair in case.contrasts.split(';'):
        a,b=pair.split('>')
        g=v[((v.a.eq(a)&v.b.eq(b))|(v.a.eq(b)&v.b.eq(a)))&v.model.isin(['full','old_strict'])].drop_duplicates(['a','b','scope','model'])
        assert len(g)>0,(case.case_id,a,b)
        for rec in g.to_dict('records'):
            rec.update(case_id=case.case_id,priority=case.priority,system=case.system,biological_reference=a,biological_variant=b,metric_direction_invariant=True,source_url=case.source_url)
            rows.append(rec)
out=pd.DataFrame(rows)
out.to_csv(P/'biology_first_12_series_17_contrasts_all5regions_2models.csv',index=False,encoding='utf-8-sig')

# This table is a reading aid, not a second selection or evaluation rule. All
# five fixed regions and both outputs for all 17 contrasts remain in the CSV.
focus={'B01':('whole_peptide','hchi'),'B02':('all_CDR','rmsf'),'B03':('whole_peptide','rmsf'),'B04':('all_CDR','rmsf'),'B05':('natural_pMHC_interface','hchi'),'B06':('all_CDR','rmsf'),'B07':('whole_peptide','rmsf'),'B08':('CDR3','rmsf'),'B09':('all_CDR','rmsf'),'B10':('natural_pMHC_interface','hchi'),'B11':('CDR3','rmsf'),'B12':('all_CDR','rmsf')}
table=[]
for (case,a,b),group in out.groupby(['case_id','biological_reference','biological_variant'],sort=False):
    scope,metric=focus[case]
    g=group[group.scope.eq(scope)]
    primary=g[g.model.eq('full')].iloc[0]
    old=g[g.model.eq('old_strict')].iloc[0]
    table.append(dict(case_id=case,a=a,b=b,display_region=scope,display_metric=metric,n=int(primary['n_'+metric]),Full_rho=primary['rho_delta_'+metric],Full_gain=primary['relative_MAE_gain_'+metric],old_strict_rho=old['rho_delta_'+metric],old_strict_gain=old['relative_MAE_gain_'+metric],MD_abs_delta=primary['MD_abs_delta_'+metric],descriptor_above_2SEM_n=int(primary['descriptor_above_2SEM_n_'+metric])))
pd.DataFrame(table).to_csv(P/'biology_first_17_contrasts_reading_table.csv',index=False,encoding='utf-8-sig')
lines=['# ITrace 生物学优先 12 组系列、17 个真实对照：数值复筛判断','',
'本文件只报告搜索与复核，不改动论文。主模型固定为现有 Full family fold0，旧 strict Full 仅平行复测，Intrachain 为既定信息控制；两套 Full 已由根任务分别重新推断。所有模型使用相同 MD 标签与相同配对位置，变化一律是生物对照中后态减前态。以下简表为阅读导引，不能代替完整五区域结果，也没有据此剔除不利子图、位点或模型。相关系数比较各位点变化的排序，MAE 改善比较“一律预测不变”；改善为负说明连这一基线也未超过。表内大于 2SEM 的位点数只是三重复 MD 的描述，不是独立显著性检验。','',
'有价值的单项证据主要有三类。KRAS G12D 的 7OW5→7OW6 对整条 10 肽主链的结果在两个 Full 模型中一致，主模型/旧 strict 的相关为 0.806/0.879，误差改善 38.9%/60.3%；排除发生替换的 G6/D6 后，全部 9 个未变位点仍为 0.733/0.833，误差改善 40.2%/62.8%。9 个位点全部超过 2SEM，方向分别 8/9 与 9/9 正确。这是完整肽段主链环境响应的支持，不能把只有 5 个位点且全部没有超过 2SEM 的侧链高相关一起包装成双输出机制成功。','',
'MEL5 的 IMP2→BST2（7Q99→7Q9A）完整 CDR3 22 个位点主链在两套 Full 中相关 0.713/0.553、误差改善 42.9%/41.8%，4 个超过 2SEM 的位点方向全部正确；Intrachain 相关 0.226、改善 18.2%、方向仅 1/4。全部六环合并也保留改善 14.4%/13.3%，说明并非孤立的一条最佳子环，但侧链结果弱。第三抗原 Melan-A 必须同列：4QOK→7Q99 明显失败，4QOK→7Q9A 主模型虽有 0.824 的完整 CDR3 相关，误差却比零变化基线高 11.0%。因此只能解释 IMP2 与 BST2 这一对的特定主链变化，不能宣称三个抗原整体准确。','',
'PRAME S2 的 A3→A11（8RYM→8RYN）在自然 pMHC 界面 15 个同 AA 侧链上，Full/旧 strict 相关 0.486/0.686、误差改善 14.3%/24.1%，其中 5 个位点超过 2SEM，方向分别 4/5 与 3/5 正确。它可以展示固定 TCR/肽而更换呈递环境后的侧链响应，但 Intrachain 相关 0.671、改善 14.7%，所以此对并未显示 Full 明显胜过 Intrachain。工程对照 8RYM→8RYO 完整 CDR3 主链在主模型为相关 0.559、改善 37.2%，旧 strict 的改善为 -4.1%；其侧链虽有约 0.55–0.59 相关，但 16 个位置无一超过 2SEM。工程对照应保留，并与 HLA 更换分开解释。','',
'DMF4 的 3QDM→3QEQ 自然 pMHC 界面侧链 19 位也是一致单项候选：两 Full 相关 0.779/0.553、改善 25.2%/27.1%，7 位超过 2SEM，主模型方向 7/7、旧 strict 5/7。并行 DMF5 对照 3QDG→3QDJ 在相同区域没有成功，因而不能把它写成所有受体均有的统一响应。1G4 CDR2 工程的 2F54→2F53 全六环主链 49 位在主模型相当强（相关 0.680，改善 59.7%）；43 个未替换位点更强（0.710，改善 63.0%），但旧 strict 误差反而增加 35.9%/45.8%。六环对应完整，不过构建体在环外缺残基，已有按整链拟合得到的 RMSF 可能受恒定区构建差异影响；目前是待构建体/拟合复核的候选，不能当跨模型重复成功。','',
'MAGE-A10、AHIII12.2 四状态、NY-ESO-1 C9V、NRAS K/R、KRAS G12V、MEL5 十肽/九肽均没有比上述候选更强的完整区域结果。MAGE 整肽侧链主模型相关 0.643，但旧 strict 为 0，且只有 1 个位置超过 2SEM；AHIII 的 K66A 全六环主链主模型相关 0.406、误差仅改善 10.7%，旧 strict 未复现；W167A 的 CDR3 侧链相关 0.557，也只有 1 位超过 2SEM。这些不应作为强机制例。','',
'完整文件：`biology_first_12_series_17_contrasts_all5regions_2models.csv`；所有 13 系列与额外对照的 8,160 项原始统计：`all_prespecified_region_model_metrics.csv`；配对逐位点与六环明确边界另存同目录。','',
'|系列|生物对照|阅读导引区域/输出|n|Full ρ / MAE改善|旧 strict ρ / MAE改善|MD平均绝对变化 / >2SEM位点|',
'|---|---|---|---:|---|---|---|']
for row in table:
    units='Å' if row['display_metric']=='rmsf' else 'nats'
    lines.append(f"|{row['case_id']}|{row['a']}→{row['b']}|{row['display_region']} / {row['display_metric']}|{row['n']}|{row['Full_rho']:.3f} / {100*row['Full_gain']:.1f}%|{row['old_strict_rho']:.3f} / {100*row['old_strict_gain']:.1f}%|{row['MD_abs_delta']:.3f} {units} / {row['descriptor_above_2SEM_n']}|")
(P/'生物优先12组_完整区域数值判断.md').write_text('\n'.join(lines)+'\n',encoding='utf8')

for case,group in out.groupby('case_id',sort=False):
    print(case,group.system.iloc[0])
    for (a,b,scope),g in group.groupby(['biological_reference','biological_variant','scope'],sort=False):
        if scope=='natural_MHC_groove':continue
        parts=[]
        for r in g.itertuples():
            parts.append(f'{r.model}: rho {r.rho_delta_rmsf:.3f}/{r.rho_delta_hchi:.3f}; gain {100*r.relative_MAE_gain_rmsf:.1f}/{100*r.relative_MAE_gain_hchi:.1f}%; MD {r.MD_abs_delta_rmsf:.3f}/{r.MD_abs_delta_hchi:.3f}; >2SEM {int(r.descriptor_above_2SEM_n_rmsf)}/{int(r.descriptor_above_2SEM_n_hchi)}')
        print(a+'>'+b,scope,'n='+str(int(g.n_rmsf.iloc[0]))+'/'+str(int(g.n_hchi.iloc[0])), ' | '.join(parts))
