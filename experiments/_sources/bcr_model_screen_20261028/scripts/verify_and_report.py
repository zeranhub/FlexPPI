"""Cross-check stored outputs and write the evidence-limited screening report."""
from pathlib import Path
import csv
import hashlib
import json
import math
import shutil
import statistics
import numpy as np

R = Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028')
def read_csv(name):
    with (R / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def truth(value):
    return value == 'True'
def close(a, b):
    assert math.isclose(float(a), float(b), abs_tol=2e-6, rel_tol=2e-6), (a, b)

states = json.loads((R / 'states.json').read_text(encoding='utf-8'))
audit = json.loads((R / 'inference_audit.json').read_text(encoding='utf-8'))
family = json.loads((R / 'runtime/family_input_audit.json').read_text(encoding='utf-8'))
assert audit['status'] == 'completed' and audit['forward_passes'] == 60
assert sha(R / 'states.json') == audit['frozen_states_sha256']
assert sha(R / 'contrasts.json') == audit['frozen_contrasts_sha256']
assert family['source_sha256'][str(R / 'states.json')] == sha(R / 'states.json')
pred = read_csv('all_states_per_residue_predictions.csv')
lookup = {(p['state_id'], p['variant'], p['chain'], p['author_residue_id']): p for p in pred}
assert len(lookup) == len(pred) == audit['rows']
prediction_groups = {}
for p in pred:
    prediction_groups.setdefault((p['state_id'], p['variant']), []).append(p)
for s in states:
    assert sha(s['path']) == s['sha256']
    for variant in ['full', 'intrachain']:
        rows = prediction_groups[s['id'], variant]
        data = np.load(R / 'predictions' / variant / (s['id'] + '.npz'))
        assert len(rows) == len(data['aa'])
        assert str(data['input_sha256']) == s['sha256']
        entry = next(x for x in audit['predictions'] if x['state_id'] == s['id'] and x['variant'] == variant)
        assert entry['chains'] == list(dict.fromkeys(p['chain'] for p in rows))
        for i, p in enumerate(rows):
            assert (p['chain'], p['author_residue_id'], p['aa']) == (data['chain'][i], data['author'][i], data['aa'][i])
            close(p['rmsf'], data['mean'][i, 0]); close(p['hchi'], data['mean'][i, 1])
            assert truth(p['hchi_chemical_valid']) == (truth(p['model_mask']) and truth(p['geometry_valid']) and p['aa'] not in ['A', 'G'])

contrasts = json.loads((R / 'contrasts.json').read_text(encoding='utf-8'))
by_id = {s['id']: s for s in states}
def backbone(path):
    atoms = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith('ENDMDL'): break
        if line.startswith('ATOM  ') and line[12:16].strip() in ['N', 'CA', 'C', 'O']:
            atoms[line[21], line[22:27].strip(), line[12:16].strip()] = line[30:54]
    return atoms
mutant_backbones_checked = []
for s in states:
    if s['kind'] != 'fixed-backbone sequence counterfactual': continue
    con = next(c for c in contrasts if c['after'] == s['id'])
    assert backbone(s['path']) == backbone(by_id[con['before']]['path'])
    mutant_backbones_checked.append(s['id'])

deltas = read_csv('all_contrasts_complete_per_position.csv')
contrast_map = {c['id']: c for c in contrasts}
for d in deltas:
    c = contrast_map[d['contrast']]
    b = lookup[c['before'], d['variant'], d['chain'], d['author_residue_id']]
    a = lookup[c['after'], d['variant'], d['chain'], d['author_residue_id']]
    for metric in ['rmsf', 'hchi']:
        close(d['delta_' + metric], float(a[metric]) - float(b[metric]))
    if truth(d['valid_hchi']):
        assert b['aa'] == a['aa'] and b['aa'] not in ['A', 'G']

regions = read_csv('all_contrasts_predefined_regions.csv')
checked_regions = 0
for r in regions:
    if r['scope'] not in ['complete_Fv', 'H_CDR1', 'H_CDR3', 'antigen_loop_73_78']: continue
    selected = []
    for d in deltas:
        if (d['contrast'], d['variant']) != (r['contrast'], r['variant']): continue
        if not truth(d['valid_' + r['metric']]): continue
        if r['scope'] == 'complete_Fv':
            keep = bool(d['imgt_number_before']) and bool(d['imgt_number_after'])
        elif r['scope'].startswith('H_'):
            keep = d['antibody_role_before'] == d['antibody_role_after'] == 'H' and d['imgt_region_before'] == d['imgt_region_after'] == r['scope'][2:]
        else:
            keep = d['chain'] == 'A' and 73 <= int(d['author_number_before']) <= 78
        if keep: selected.append(d)
    assert len(selected) == int(r['n'])
    close(r['mean_delta'], statistics.mean(float(d['delta_' + r['metric']]) for d in selected))
    checked_regions += 1

md = read_csv('exact_input_predictions_and_MD_labels.csv')
metrics = read_csv('MD_absolute_prediction_metrics.csv')
kw_h3 = [m for m in md if m['sample_id'] == 'A00KW' and m['variant'] == 'full' and m['antibody_role'] == 'H' and m['imgt_region'] == 'CDR3']
assert len(kw_h3) == 9 and all(float(m['rmsf_pred']) < float(m['rmsf_MD']) for m in kw_h3)
chi_h3 = [m for m in kw_h3 if truth(m['chi_entropy_valid'])]
aa_counts = {aa: sum(m['aa'] == aa for m in chi_h3) for aa in sorted({m['aa'] for m in chi_h3})}
assert aa_counts == {'D': 1, 'N': 1, 'R': 1, 'Y': 2}
with (R / 'reference/D44_1_HCDR3_numbering_and_predictions.csv').open('w', encoding='utf-8-sig', newline='') as f:
    fields = ['chain', 'author_residue_id', 'imgt_number', 'aa', 'rmsf_pred', 'rmsf_MD', 'hchi', 'chi_entropy', 'chi_entropy_valid']
    w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore'); w.writeheader(); w.writerows(kw_h3)

portable = R / 'reference/md_input_structures'
portable.mkdir(exist_ok=True)
for s in states:
    if s['id'].startswith('mdref_'):
        target = portable / (s['id'] + '.pdb')
        shutil.copyfile(s['path'], target)
        assert sha(target) == s['sha256']

verification = {'status': 'passed', 'method': 'Additional CSV/NPZ/source checks by root; not an independent scientific reviewer.',
    'states': len(states), 'forward_passes': 60, 'prediction_rows_checked': len(pred), 'contrast_rows_checked': len(deltas),
    'region_rows_recomputed_with_csv_module': checked_regions, 'fixed_backbone_mutants_verified': len(mutant_backbones_checked),
    'states_sha256': sha(R / 'states.json'), 'contrasts_sha256': sha(R / 'contrasts.json'),
    'MD_HCDR3_all_nine_underpredicted': True, 'MD_HCDR3_valid_Hchi_aa_counts': aa_counts,
    'potential_train_bridge_states': family['states_with_potential_train_bridge'],
    'auditing_note': 'A stale variable originally made inference_audit chains list repeat the final sample chains. Corrected the logging code and reran all 60 forward passes; final CSV and NPZ chain identities checked against source.',
    'files_sha256': {n: sha(R / n) for n in ['all_states_per_residue_predictions.csv', 'all_contrasts_complete_per_position.csv', 'all_contrasts_predefined_regions.csv', 'MD_absolute_prediction_metrics.csv', 'reference/antibody_residue_numbering.csv']}}
(R / 'final_output_verification.json').write_text(json.dumps(verification, ensure_ascii=False, indent=2), encoding='utf-8')

def region(con, scope, metric='rmsf', variant='full'):
    return next(r for r in regions if (r['contrast'],r['scope'],r['metric'],r['variant']) == (con,scope,metric,variant))
def metric(sample, scope='all', target='rmsf', variant='full'):
    return next(m for m in metrics if (m['sample_id'],m['scope'],m['metric'],m['variant']) == (sample,scope,target,variant))
def fmt(row, key, digits=3): return f"{float(row[key]):.{digits}f}"
ag = region('flex_joint_state_change_full', 'antigen_loop_73_78')
h1 = region('flex_joint_state_change_full', 'H_CDR1')
copy_h1 = region('flex_crystal_copy_sensitivity', 'H_CDR1')
copy_ag = region('flex_crystal_copy_sensitivity', 'antigen_loop_73_78')
r71 = region('bound_HR71A', 'H_CDR3')
r45 = region('bound_LR45A', 'H_CDR3')
kw = metric('A00KW'); kw_chi = metric('A00KW',target='hchi'); h3 = metric('A00KW','H_CDR3')
main_rows = []
for con in ['bound_HR71A','bound_LR45A','apo_HR71A','apo_LR45A','flex_joint_state_change_full','flex_joint_state_change_fv','flex_4X_apo_to_bound_full','flex_crystal_copy_sensitivity']:
    main_rows.extend(r for r in regions if r['contrast'] == con and r['variant'] == 'full' and r['scope'] in ['complete_Fv','H_CDR1','H_CDR3','antigen_loop_73_78'])
with (R / 'screening_key_results.csv').open('w', encoding='utf-8-sig', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(regions[0])); w.writeheader(); w.writerows(main_rows)

text = f'''# BCR候选实际试跑与案例筛选结论（2026-10-28）

已经对30种结构输入分别运行完整FlexPPI及仅保留链内信息的对照，共60次实际前向计算；保留17,648条逐残基输出和21组比较的全部共享位置。这里研究的是抗体Fab/Fv与抗原，属于BCR的抗原识别部分，未模拟完整膜上BCR。当前结论是：**尚未找到足以替换现有TCR案例、并能可靠解释抗体单点突变效应的新案例。HyHEL-10暂不推荐作为成功案例；Flex-HEL的抗原柔性环值得保留，但还没有形成抗体侧链机制的证据链。** 本轮未修改论文，没有训练新模型或运行新MD。

## 实际输入与模型

HyHEL-10使用1C08结合态、5AYU游离态以及从1C08去掉抗原的输入，分别测试WT、H:R71A、L:R45A以及近似无亲和力变化的H:Q3A/L:H34A对照，共15种。14种总计的人工突变输入全部只改残基身份、删除被替换位的侧链原子，主链坐标逐项核实未变；模型输入使用主链几何与残基身份，这些输入不是实验测得的突变体结构，也不是重新平衡的构象。Flex-HEL使用6P4A、6P4B、6P4C和6P4D，另保留完整Fab、截取Fv、两个晶体拷贝、去掉抗体的抗原及固定骨架四突变/回复输入。Fv轻链截止作者编号107，重链截止113；识别环通过ANARCII独立映射至IMGT编号，比较保留全部共享位置。另有A00EJ/A00HE/A00KW三个原始MD输入供绝对值核查。

采用当前冻结代码和全组成链家族划分的fold 0、epoch 40检查点。完整模型SHA256为`{audit['predictions'][0]['checkpoint_sha256']}`，链内对照SHA见inference_audit.json。两种模型均严格加载全部参数；先各自预测绝对状态，再相减。本轮的两种模型不是两个外部亲和力预测器，也不输出亲和力。

## HyHEL-10：单点突变输入尚未产生足够有说服力的结果

在结合态固定骨架上，H:R71A的重链第三识别环（H-CDR3）平均主链波动幅度RMSF由{fmt(r71,'mean_before')}变为{fmt(r71,'mean_after')} Å，仅增加{fmt(r71,'mean_delta')} Å；L:R45A由{fmt(r45,'mean_before')}变为{fmt(r45,'mean_after')} Å，降低{abs(float(r45['mean_delta'])):.3f} Å。后者呈现一些跨链响应，但这两种突变没有共同形成清楚、稳健且已有同类动态真值支持的模式。两个近似无亲和力变化的对照也有输出变化；它们不能视作已证明“动态完全不变”的阴性。

游离态NMR正文点名且残基身份可直接对上的8个R71A响应位中，本模型RMSF下降的有0个；R45A可对上的10个位中有2个。这个计数仅记录方向，**不是预测命中率或失败率**：NMR检测毫秒级构象交换，模型的RMSF表示MD尺度下位置摆动的大小，二者不能逐值等同；NMR观察到交换信号消失也不能自动换算成RMSF下降。正文还有H:T56、L:H33两处残基身份/编号不一致，已保留为未解决映射，未改号挑结果。现有材料不能支持“FlexPPI已经解释该突变导致亲和力提高”，也不能据此认定模型对这类突变普遍无效。

## Flex-HEL：抗原环信号清楚，抗体解释仍不充分

以完整Fab输入比较6P4A的WT抗体–rigid抗原与6P4B的4X抗体–flex抗原，预先指定的抗原73–78环平均RMSF从{fmt(ag,'mean_before')}升至{fmt(ag,'mean_after')} Å，增加{fmt(ag,'mean_delta')} Å，6个位置全部升高；Fv输入得到约0.615 Å的同向增幅。同一6P4B的两个晶体拷贝，该环均值差只有{fmt(copy_ag,'mean_delta')} Å，但逐位平均绝对差为{fmt(copy_ag,'mean_abs_delta')} Å，不能用均值差掩盖局部差异。另将预先保留的两种“去掉抗体”输入作描述性诊断，完整模型该环均值约1.080和2.224 Å，方向一致；这些坐标仍来自结合态，不能称为实测游离结构。

抗体的结果弱得多：重链第一识别环（H-CDR1）在上述两状态间平均RMSF只增加{fmt(h1,'mean_delta')} Å，而同一复合物两个晶体拷贝之间增加{fmt(copy_h1,'mean_delta')} Å，后者更大。预先列出的H:W34/I29/Y53、抗原D101/W62侧链网络全部输出已保存；晶体显示的侧链构象重排并不提供侧链熵升降真值，不能挑几个Hχ变化较大的位置就写成机制验证。

原论文的重复MD支持的是**游离抗原**73–78环在flex中更灵活。本轮主要对比的两个复合物同时改变抗体与抗原，因而目前最多得到与该局部现象方向一致的模型信号；不能分离是哪个突变造成，也不能证明“抗体通过侧链柔性改善亲和力”。6P4D的101位为Arg，6P4A/B为Asp，亦不能当作严格相同序列的rigid游离对照。Flex-HEL是继续筛选的候选，现阶段不宜作为新的抗体机制主案例。

## 用本地原始MD标签核查，保留不利结果

这些检查严格使用原始MD输入，未把Fab预测与只有抗原的轨迹混合。下表RMSF误差为逐残基绝对误差的均值，单位Å；Hχ为侧链转动角组合的离散状态熵，单位nats。排序相关越接近1，表示模型给出的高低顺序越接近轨迹；它不代表数值误差很小。

| 样本与实际组成 | RMSF平均误差 | RMSF排序相关 | Hχ平均误差 | Hχ排序相关 |
|---|---:|---:|---:|---:|
'''
for sample, desc in [('A00EJ','A00EJ / 6P4B，仅抗原128位'),('A00HE','A00HE / 1VFB，仅抗原128位'),('A00KW','A00KW / 1MLC，D44.1 Fab–抗原553位')]:
    a=metric(sample); b=metric(sample,target='hchi')
    text += f"| {desc} | {fmt(a,'mae')} | {fmt(a,'spearman')} | {fmt(b,'mae')} | {fmt(b,'spearman')} |\n"
text += f'''
D44.1复合物整体的RMSF误差是{fmt(kw,'mae')} Å，侧链熵排序相关是{fmt(kw_chi,'spearman')}；但是，先扣掉同一种氨基酸各自的平均值、再比较同型残基内部差异，Hχ排序相关降为{fmt(kw_chi,'same_aa_centered_spearman')}，说明不能把较高的总体相关全部解释为捕获局部环境。更关键的是，该抗体H-CDR3的9个位置RMSF全部低估，预测均值{fmt(h3,'mean_prediction')} Å、MD均值{fmt(h3,'mean_MD')} Å，平均误差{fmt(h3,'mae')} Å；链内对照在该环的误差反而较小，约0.265 Å。该环仅5个有效侧链熵位置，只有两个Tyr能在同一种氨基酸内部比较，不能用此处相关等于1来宣称成功。完整模型在全复合物的平均误差小于链内对照，但排序相关和部分链/区域未同时改善，全部结果均保留。

三个样本各只有一条本地轨迹。A00EJ/A00HE的5个时间窗口不是5次独立模拟，A00KW只有1个1000帧窗口且记录时间约0.999 ns，不适合据此声称充分采样的抗体机制。原始标签已从轨迹重算；三个RMSF最大差小于10⁻⁶ Å，A00EJ/A00KW的Hχ完全一致；A00HE有一位Hχ相差0.000134 nats，原因尚未确定，原始及重算值均保留。这些是既有测试记录的绝对状态核查，不是新获得的突变动态真值。

## 训练重叠核查及可用边界

检查当前训练池3063个结构的4973条组成链，而未直接依赖不完整的旧序列文件。A00EJ/A00HE/A00KW都属于当前fold 0留出的同一个633记录家族组件，该组件在该检查点训练清单中的记录数为0。对新候选又做了27条不同输入链序列的穷尽局部/全局比对，使用30%一致性、双向覆盖均至少70%的阈值。由于本地没有MMseqs可执行文件，新查询采用的是公开记录参数的Biopython替代核查，不能称为严格复现原家族搜索。

发现一个边界附近的潜在训练桥接：HyHEL/Hy10未改变的107位轻链Fv与训练样本A003O的B链，局部比对一致性30.49%，输入覆盖72.90%，训练链覆盖75.00%。共影响14种输入。它是**潜在同源重叠**，不是证明模型记忆了该案例；但这些输入也不能被宣称为严格排除该关系的独立泛化证据。完整Fab在整条链覆盖阈值下没有同一命中，不能因此声称其中的可变区已经排除此关系。若要把它们用于论文中的严格新案例，需要用原MMseqs协议确认并排除所有相关训练记录，再用独立动态参考复核；本轮未为追求正结果重新训练。

## 对论文的处理建议与交付

保留现有CaM、泛素、TCR案例；本轮没有实测依据证明新增抗体案例更强，暂不替换或扩写为成功机制案例。HyHEL-10目前只适合保留为探索记录；Flex-HEL可继续作为抗原动态候选，只有补足匹配状态、抗体自身的动态参考和同源排除后，才考虑升级为主案例。此判断针对本轮已测试候选与输入，未测试CH65等其他体系，不能据此否定整个BCR方向。

全部逐残基预测：all_states_per_residue_predictions.csv；全部比较：all_contrasts_complete_per_position.csv；预定义区域：all_contrasts_predefined_regions.csv；便于查看的摘要：screening_key_results.csv；真实MD核查：MD_absolute_prediction_metrics.csv及exact_input_predictions_and_MD_labels.csv；输入、模型、家族及最终一致性审计：states.json、inference_audit.json、runtime/family_input_audit.json、final_output_verification.json。所有人工突变的主链不变、CSV与60份原始NPZ一致、比较差值及主要区域统计均已另行核验。脚本包不包含模型权重和原始XTC，原路径和文件SHA已记录，因此它是结果与核查材料包，并非无需外部数据即可运行全部流程的独立环境。

文献背景：[2017 HyHEL-10 fluctuation editing](https://doi.org/10.1038/s41598-017-10246-9)；[2020 Hy10/Flex-HEL](https://doi.org/10.1073/pnas.2005102117)。结构来源为RCSB，结构原文件及SHA均保留。上文的所有模型数值来自本轮实际输出，而非文献推算。
'''
(R / '结果说明与案例筛选结论.md').write_text(text, encoding='utf-8')
print(json.dumps({'verification': 'passed', 'prediction_rows': len(pred), 'delta_rows': len(deltas), 'checked_region_summaries': checked_regions, 'fixed_backbone_mutants': len(mutant_backbones_checked), 'report': str(R / '结果说明与案例筛选结论.md')}, ensure_ascii=False))
