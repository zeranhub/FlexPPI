import collections
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
EXP = ROOT.parents[1]
FAMILY = EXP / 'family_folds'
OLD = Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def records(p):
    return [json.loads(x) for x in Path(p).read_text(encoding='utf-8-sig').splitlines() if x.strip()]

remote = json.loads((ROOT / 'remote_checkpoint_manifest_audit.json').read_text(encoding='utf-8-sig'))
train = records(FAMILY / 'splits/fold_0_train.jsonl')
test = records(FAMILY / 'splits/fold_0_test.jsonl')
ready = records(FAMILY / 'splits/all_ready.jsonl')
strict = records(OLD / 'strict_external_payload/splits/strict_external_train.jsonl')
freeze = json.loads((FAMILY / 'training_freeze.json').read_text())
groups = json.loads((FAMILY / 'family_components.json').read_text())
ids = lambda rows: {r['sample_id'] for r in rows}
itrace = [r for r in ready if r['dataset'] == 'itrace']
checks = {
    'family_train_records': len(train),
    'family_train_sources': dict(collections.Counter(r['dataset'] for r in train)),
    'family_test_records': len(test),
    'family_test_sources': dict(collections.Counter(r['dataset'] for r in test)),
    'train_test_sample_overlap': sorted(ids(train) & ids(test)),
    'all_ready_ITrace_count': len(itrace),
    'ITrace_train_overlap': sorted(ids(itrace) & ids(train)),
    'all_ITrace_in_fold0_test': ids(itrace) <= ids(test),
    'train_family_component_ids': len({r['metadata']['family_component_id'] for r in train}),
    'test_family_component_ids': sorted({r['metadata']['family_component_id'] for r in test}),
    'train_test_family_overlap': sorted({r['metadata']['family_component_id'] for r in train} & {r['metadata']['family_component_id'] for r in test}),
    'strict2431_minus_family2430_ids': sorted(ids(strict) - ids(train)),
    'family2430_minus_strict2431_ids': sorted(ids(train) - ids(strict)),
    'strict_ITrace_overlap': sorted(ids(strict) & ids(itrace)),
    '3KXF_train_records': [r['sample_id'] for r in train if '3kxf' in json.dumps(r).lower()],
    '3KXF_test_records': [r['sample_id'] for r in test if '3kxf' in json.dumps(r).lower()],
    '3KXF_strict_train_records': [r['sample_id'] for r in strict if '3kxf' in json.dumps(r).lower()],
    'local_frozen_assets': {name: sha(FAMILY / name) == expected for name, expected in freeze['files'].items()},
}
checks['local_frozen_assets_all_match'] = all(checks['local_frozen_assets'].values())
component = [g for g in groups if g['fold'] == 0]
checks['fold0_component_inventory'] = component

rows = []
for cp in remote['checkpoints']:
    m = cp['manifests']
    c = cp['config']
    trainm = m['state_manifest']
    valm = m.get('validation_state_manifest', {})
    testm = m.get('test_state_manifest', {})
    rows.append(dict(checkpoint_path=cp['path'], checkpoint_sha256=cp['sha256'], epoch=cp['epoch'],
                     seed=c['seed'], outputs=';'.join(c['model']['metrics']), init_checkpoint=c['init_checkpoint'],
                     train_count=trainm['records'], train_ITrace=trainm['sources'].get('itrace', 0),
                     validation_count=valm.get('records', 0), val_fraction=c['train']['val_fraction'],
                     validation_ITrace=valm.get('sources', {}).get('itrace', 0),
                     test_count=testm.get('records', 0), test_ITrace=testm.get('sources', {}).get('itrace', 0),
                     train_3KXF=';'.join(r['sample_id'] for r in trainm['pdb_3kxf_records']),
                     all_ITrace_unused_in_fitting=(trainm['sources'].get('itrace', 0) == 0 and valm.get('sources', {}).get('itrace', 0) == 0)))
    if '/family_folds/' in cp['path']:
        for mk, mv in m.items():
            local = FAMILY / Path(mv['path']).name
            local = FAMILY / 'splits' / Path(mv['path']).name
            assert sha(local) == mv['sha256'], (cp['path'], mk)
with (ROOT / 'checkpoint_inventory.csv').open('w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

fetchpath = EXP / 'tcr_case_v40_20261017/runtime_inventory/family_fetch_audit.json'
fetch = json.loads(fetchpath.read_text())
prediction_manifest = records(EXP / 'tcr_case_v40_20261017/family_predictions/prediction_manifest.jsonl')
pred_checks = []
for variant in ['full', 'intrachain']:
    matched = [cp for cp in rows if '/family_folds/' in cp['checkpoint_path'] and f'/runs/{variant}/fold_0/' in cp['checkpoint_path']][0]
    files = [x for x in fetch['files'] if x['variant'] == variant]
    manifestrows = [r for r in prediction_manifest if r['variant'] == variant]
    for r in files:
        assert Path(r['local_path']).exists()
        assert sha(r['local_path']) == r['sha256'], r['local_path']
        assert r['checkpoint_sha256'] == matched['checkpoint_sha256']
    assert {r['sample_id'] for r in files} == ids(itrace)
    assert len(manifestrows) == len(itrace)
    assert all(r['fold'] == 0 and r['checkpoint_sha256'] == matched['checkpoint_sha256'] for r in manifestrows)
    pred_checks.append(dict(variant=variant, count=len(files), checkpoint_sha256=matched['checkpoint_sha256'],
                            local_directory=str(EXP / f'tcr_case_v40_20261017/family_predictions/{variant}'),
                            all_file_sha256_verified=True, exact_ITrace_set_verified=True))
checks['reusable_predictions'] = pred_checks
checks['remote_checkpoint_count_inspected'] = len(rows)
(ROOT / 'verified_summary.json').write_text(json.dumps(checks, indent=2, ensure_ascii=False), encoding='utf-8')

full = [r for r in rows if '/family_folds/runs/full/fold_0/' in r['checkpoint_path']][0]
old_strict = [r for r in rows if '/corrected_training_v1/runs/strict_full/fold_0/' in r['checkpoint_path']][0]
text = f'''# ITrace整集未参与拟合的既有模型审计

本次只读检查了远程26个实际checkpoint及其中保存的配置，并对照了真实数据清单、原有冻结文件和488个已下载预测文件的SHA256；没有训练模型，没有修改论文或旧结果。

## 推荐主模型

`{full['checkpoint_path']}`，SHA256为`{full['checkpoint_sha256']}`。这是已有的联合RMSF与侧链联合构象熵模型，seed=7、固定epoch40，从随机初始化训练；`init_checkpoint=null`，没有从含ITrace的权重继续训练。实际训练清单共2430条（ATLAS 1866、DynaRepo 564），ITrace为0。`validation_state_manifest=null`且`val_fraction=0.0`，所以验证集为空，ITrace没有参与验证或轮数选择。训练实现对这种配置返回空验证集；test清单仅供导出评估，未进入训练数据加载。

真实fold0测试清单633条，含ITrace 244、ATLAS 63和DynaRepo 326；整份all_ready中全部244条ITrace均在这里。训练与测试sample_id和完整家族component_id都无交集。ITrace所在的关联组为`family_9cd89939c749e38c`，含633个主数据记录及32个只参与连接的外部记录。关联规则包括任意链的精确序列、长链序列identity≥30%且双向覆盖≥70%、短于20残基链的精确连接及相同生物系统，经完整复合物传递连接。此审计复核清单和冻结文件一致，未重新运行全部序列搜索。

DynaRepo中的3KXF是`A005I`。它在fold0测试集，训练集没有该PDB；因此不能把它误当作漏进这个模型的ITrace精确重复。当前主清单是去重后的244个ITrace条目，不能写成对原始ITrace所有可能文件均完成了独立推理。

## 用户所说的更早“不含ITrace”的模型

旧的`{old_strict['checkpoint_path']}`也满足整集ITrace未参与拟合，SHA256为`{old_strict['checkpoint_sha256']}`，训练2431条（ATLAS 1866、DynaRepo 565），验证集为空，epoch40、seed7、联合双输出、随机初始化。真实训练sample_id集合逐条比较后，差异恰好是旧strict多了`A00IL`（PDB 4WV1，抗原单链）；没有其他条目差异。旧同源方案删除632条，新方案因补充全链搜索和短链连接删除633条。历史`single_family_review.md`曾将新增条目描述为A009Z，该说法与这两份真实训练manifest差集不一致，本审计以真实清单差集为准：A009Z在两个训练清单中均不出现。旧strict模型可作预先确定的历史方案复测，但不能按哪个候选效果更好临时选择主模型；主用更严格的2430模型。

相反，旧`corrected_training_v1/runs/full/fold_0..4/epoch_040.pt`的训练集分别包含195、196、195、195、195条ITrace；seed11和23的fold0也各包含195条。它们仅为ITrace部分条目的折外预测，不能称为“整个ITrace没训练过”。新家族fold1至4也都在训练中含244条ITrace，不能混合为整集未用ITrace的集成。完整26个模型的路径、SHA、训练/验证/测试计数见`checkpoint_inventory.csv`。

## 可以直接复用的测量

`{pred_checks[0]['local_directory']}`已有全部244个ITrace双输出预测；每个文件的预测manifest均指向上述主模型的SHA，全部文件SHA256本次复核通过。其intrachain对照在`{pred_checks[1]['local_directory']}`，checkpoint SHA为`{pred_checks[1]['checkpoint_sha256']}`。这意味着上一轮SB27/SB47、MEL5、PIK3CA及JM22的family_predictions筛查本来已经采用整个ITrace未参与拟合的模型，不能把重复读取这些缓存描述为一次新的独立预测实验。候选重新搜索、修正序列对应和从轨迹复算可以继续使用这些一致的输出。

模型训练时未使用ITrace，并不等于现在的候选挑选是盲评：多个条目元数据保留了历史诊断探索标记，而且这次是在查看预测后挑选说明性实例。研究记录应保留全部候选结果与预设区域，不能把挑出的最好实例当作整个数据集能力的无偏估计。

## 审计文件

- `remote_checkpoint_manifest_audit.json`：26个实际保存模型的配置、SHA及其真实manifest来源统计。
- `verified_summary.json`：主模型清单、家族范围、3KXF、旧strict差异、冻结文件hash和复用预测复核。
- `checkpoint_inventory.csv`：完整模型选择表。
- `audit_remote_models.py`与`summarize_audit.py`：可重放只读审计脚本。
'''
(ROOT / 'ITrace未参与训练的既有模型_审计结论.md').write_text(text, encoding='utf-8')
print(json.dumps({k: v for k, v in checks.items() if k not in ['fold0_component_inventory', 'local_frozen_assets']}, ensure_ascii=False, indent=2))
