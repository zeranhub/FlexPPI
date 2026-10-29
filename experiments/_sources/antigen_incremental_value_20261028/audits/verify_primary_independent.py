from pathlib import Path
import csv, hashlib, json, math, statistics
from collections import defaultdict

R = Path(r'E:\Project\FlexPPI\experiments\_sources\antigen_incremental_value_20261028')
E = R/'evaluation/full_primary'
def read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f: return list(csv.DictReader(f))
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def label(value): return int(str(value).lower() in ['1','true'])
def ap(rows):
    ordered=sorted(rows,key=lambda r:-float(r['score']))
    npos=sum(label(r['functional_positive']) for r in ordered)
    if npos==0: return 0.0
    total=hits=0; terms=[]; i=0
    while i<len(ordered):
        j=i+1
        while j<len(ordered) and float(ordered[j]['score'])==float(ordered[i]['score']): j+=1
        positive=sum(label(r['functional_positive']) for r in ordered[i:j])
        total+=j-i; hits+=positive
        terms.append(positive/npos * hits/total)
        i=j
    return sum(terms)
def auc(rows):
    pos=[float(r['score']) for r in rows if label(r['functional_positive'])]
    neg=[float(r['score']) for r in rows if not label(r['functional_positive'])]
    if not pos or not neg: return None
    return sum((a>b)+0.5*(a==b) for a in pos for b in neg)/len(pos)/len(neg)
def metrics(rows):
    ordered=sorted(rows,key=lambda r:(-float(r['score']),r['candidate_id']))
    return {'n':len(rows),'n_positive':sum(label(r['functional_positive']) for r in rows),
        'top10_hits':sum(label(r['functional_positive']) for r in ordered[:10]),
        'average_precision':ap(rows),'auroc':auc(rows),'top10_candidates':'|'.join(r['candidate_id'] for r in ordered[:10])}

candidates=read(R/'functional_candidates.csv'); by_id={c['candidate_id']:c for c in candidates}
assert len(by_id)==279
pred=read(E/'all_heldout_predictions.csv')
members=json.loads((E/'fit_membership.json').read_text())
stored=read(E/'metrics_by_system_and_fold.csv')
components={r['candidate_id']:r['cluster'] for r in read(E/'sequence_components.csv')}
for row in pred:
    c=by_id[row['candidate_id']]
    assert c['system']==row['system'] and label(c['functional_positive'])==label(row['functional_positive'])
    assert math.isclose(float(c['functional_response']),float(row['functional_response']),abs_tol=1e-10)
    assert math.isfinite(float(row['score'])) and 0<=float(row['score'])<=1
groups=defaultdict(list)
for row in pred: groups[row['split'],row['fold'],row['model']].append(row)

for fit in members:
    train=set(fit['training_ids']); test=set(fit['test_ids'])
    assert not train&test
    assert len(train)==len(fit['training_ids']) and len(test)==len(fit['test_ids'])
    assert train|test <= by_id.keys()
    assert set(r['candidate_id'] for r in groups[fit['split'],fit['fold'],fit['model']])==test
    assert set(fit['features']).isdisjoint({'functional_positive','functional_response','candidate_id','system','generation_protocol','wt_for_design','resp','resp_el4'})
    if fit['split']=='leave_system_out':
        assert len(train)==186 and len(test)==93
        assert {by_id[c]['system'] for c in test}=={fit['fold']}
        assert fit['fold'] not in {by_id[c]['system'] for c in train}
        assert train|test==by_id.keys()
    else:
        assert len({by_id[c]['system'] for c in train|test})==1
        assert not {components[c] for c in train}&{components[c] for c in test}

recomputed=[]; checked=0
for key, rows in groups.items():
    assert len({r['candidate_id'] for r in rows})==len(rows)
    calc=metrics(rows)
    s=next(r for r in stored if (r['split'],r['fold'],r['model'])==key)
    for name in ['n','n_positive','top10_hits']: assert int(s[name])==calc[name]
    for name in ['average_precision','auroc']:
        if calc[name] is not None: assert math.isclose(float(s[name]),calc[name],abs_tol=1e-12)
    assert s['top10_candidates']==calc['top10_candidates']
    checked+=1
    if key[0]=='leave_system_out': recomputed.append({'system':key[1],'model':key[2],**calc})

within=defaultdict(list)
for row in pred:
    if row['split']=='within_system_cluster': within[row['system'],row['model']].append(row)
within_calc=[]
for (system,model),rows in within.items():
    assert len(rows)==93 and len({r['candidate_id'] for r in rows})==93
    calc=metrics(rows); within_calc.append({'system':system,'model':model,**calc})
    s=next(r for r in stored if (r['split'],r['fold'],r['model'])==('within_system_cluster_aggregated',system,model))
    assert int(s['top10_hits'])==calc['top10_hits']
    assert math.isclose(float(s['average_precision']),calc['average_precision'],abs_tol=1e-12)

macros={}
for model in sorted({r['model'] for r in recomputed}):
    rows=[r for r in recomputed if r['model']==model]; assert len(rows)==3
    macros[model]={'mean_top10_hits':statistics.mean(r['top10_hits'] for r in rows),
        'mean_average_precision':statistics.mean(r['average_precision'] for r in rows)}

shuffled=read(E/'shuffled_dynamic_controls_all_seeds.csv'); seed_groups=defaultdict(list)
assert len(shuffled)==600
for row in shuffled: seed_groups[int(row['seed'])].append(row)
assert set(seed_groups)==set(range(200))
null=[]
for seed,rows in sorted(seed_groups.items()):
    assert {r['fold'] for r in rows}=={'nyeso','ebv','mage'} and len(rows)==3
    null.append({'seed':seed,'mean_top10_hits':statistics.mean(float(r['top10_hits']) for r in rows),
        'mean_average_precision':statistics.mean(float(r['average_precision']) for r in rows)})
actual=macros['B2_binding_static_dynamics']
null_summary={k:{'median':statistics.median(n[k] for n in null),'min':min(n[k] for n in null),'max':max(n[k] for n in null),
    'fraction_greater_or_equal_observed':sum(n[k]>=actual[k]-1e-14 for n in null)/200} for k in ['mean_top10_hits','mean_average_precision']}

result={'status':'passed','scope':'Full, primary template only; independent CSV-module recomputation, not a new model fit.',
 'blindness':'Pre-output protocol review completed before viewing prediction/metric files; unblinded at explicit parent instruction.',
 'independent_TCR_systems':3,'candidate_rows':279,'fit_membership_checks':len(members),'stored_metric_rows_crosschecked':checked,
 'LOSO_system_separation':'passed','within_system_Hamming_component_separation':'passed',
 'labels_match_original_candidates':'passed','no_label_or_identifying_fields_in_logged_feature_lists':'passed',
 'source_code_manual_review':'Pipeline.fit uses only train rows; median imputation/scaling are inside pipeline; explicit BASE/STATIC/DYNAMIC lists contain no labels. Features-and-labels CSV contains labels for evaluation but these are not fed as model features.',
 'primary_metrics':recomputed,'macro_metrics':macros,'within_system_aggregates':within_calc,
 'permutation_control_summary':null_summary,
 'permutation_caveat':'Only stored per-seed metric summaries are available; distribution checked, each shuffled model was not independently refitted. These fractions are descriptive, not confirmatory significance with three systems.',
 'skipped_analyses':json.loads((E/'skipped_analyses.json').read_text()),
 'conclusion':'No stable incremental functional-screening benefit from the eight frozen FlexPPI regional features on the primary structures. B2 macro top10 and AP are worse than B1; AA-control comparison also worsens. No direct binding/affinity or causal flexibility conclusion follows.',
 'source_sha256':{str(p.relative_to(R)):sha(p) for p in [E/'all_heldout_predictions.csv',E/'metrics_by_system_and_fold.csv',E/'fit_membership.json',E/'features_and_labels.csv',E/'shuffled_dynamic_controls_all_seeds.csv',R/'functional_candidates.csv',R/'scripts/evaluate_incremental_utility.py']}}
(R/'audits/independent_primary_result_review.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')

lines=['# Full主模板结果的独立只读复核','',
'本复核已按明确授权解除结果盲态；此前的协议盲审保持原样。使用Python csv模块从all_heldout_predictions.csv重新排序、按并列分数整体处理独立重算平均精度AP，并从逐对正负比较重算AUROC，未调用原评价脚本的指标函数。逐候选标签与functional_candidates.csv一致，全部可对照的指标及前10候选名单均与保存文件一致。',
'', '| 体系 | B0命中/10 | B1命中/10 | B2命中/10 | B0 AP | B1 AP | B2 AP |','|---|---:|---:|---:|---:|---:|---:|']
for system in ['nyeso','ebv','mage']:
    rows=[next(r for r in recomputed if r['system']==system and r['model']==m) for m in ['B0_binding','B1_binding_static','B2_binding_static_dynamics']]
    lines.append('| '+system+' | '+' | '.join(str(r['top10_hits']) for r in rows)+' | '+' | '.join(f"{r['average_precision']:.6f}" for r in rows)+' |')
lines += ['',f"宏平均B1前10命中{macros['B1_binding_static']['mean_top10_hits']:.3f}个、AP {macros['B1_binding_static']['mean_average_precision']:.6f}；加入动态后的B2为{actual['mean_top10_hits']:.3f}个、AP {actual['mean_average_precision']:.6f}。组成控制B1_AA→B2_AA同样下降，前10由{macros['B1_AA']['mean_top10_hits']:.3f}变为{macros['B2_AA']['mean_top10_hits']:.3f}，AP由{macros['B1_AA']['mean_average_precision']:.6f}变为{macros['B2_AA']['mean_average_precision']:.6f}。逐位置氨基酸控制的AP虽有小幅变化，前10命中下降，不构成一致收益。",
 '',f"200个体系内打乱动态向量的对照，宏平均前10命中中位数{null_summary['mean_top10_hits']['median']:.3f}、AP中位数{null_summary['mean_average_precision']['median']:.6f}；达到或超过实际B2值的比例分别为{null_summary['mean_top10_hits']['fraction_greater_or_equal_observed']:.3f}和{null_summary['mean_average_precision']['fraction_greater_or_equal_observed']:.3f}。这只是描述性噪声/自由度对照，不是三体系下的正式显著性检验。",
 '', 'fit_membership内每个整体系留出折均是186训练/93测试，训练体系与测试体系无交叉，候选无重复；同体系补充划分的训练/测试亦无Hamming连通分量交叉。记录的特征列表不含真实标签、响应强度、候选ID、体系ID或原设计方法。代码的缺失处理与标准化在训练管线内拟合，未发现测试标签或测试特征参与拟合的直接泄漏。NYESO的同体系补充划分因训练折只有单一标签而跳过，不能把仅两体系的补充宏平均当作三体系结果。',
 '', '主模板结论是本轮八个预定义动态特征没有显示稳定筛选增量，且关键汇总指标下降。固定骨架、功能标签、作者其他建模流程的PAE、三个独立受体体系及基模型训练覆盖限制，均沿用先前盲审边界；不能以279候选制造独立外推样本数。本复核未把后续备用模板结果纳入结论，也未修改原预测或协议。']
(R/'audits/independent_primary_result_review.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps({'status':'passed','primary_rows':len(recomputed),'memberships':len(members),'macros':macros,'null':null_summary},ensure_ascii=False))
