"""Secondary, audit-derived comparisons on all published ATLAS ID categories.

Membership is frozen before this script reads performance. This is not a new
primary test set, an equal-training contest, or proof of unseen homologs.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
REF=ROOT/'references'
OUT=ROOT/'evaluation'
EXPECTED_MEMBERSHIP_SHA='72521f5abc8609ce7bd4dcf10cf3e1e6c0d16626017f06e9a7be98fdf6942e5b'
CASES=[('DynaProt-M','DynaProt-M'),('backflip-1.0','backflip-1.0 (assembly)'),
       ('backflip-1.0-seq','backflip-1.0-seq (assembly)')]
STRATA=['all_atlas','train','val','test','unlisted']
METRICS=['spearman','top10_ap','mae','z_mae']

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    membership_path=REF/'training_exposure_membership_frozen.csv'
    assert sha(membership_path)==EXPECTED_MEMBERSHIP_SHA
    source_audit=json.loads((REF/'training_exposure_protocol_sources_sha256.json').read_text(encoding='utf-8'))
    assert source_audit['performance_tables_read_in_this_audit']==[]
    assert source_audit['output_hashes'][membership_path.name]==EXPECTED_MEMBERSHIP_SHA
    protocol_path=REF/'TRAINING_EXPOSURE_SENSITIVITY_PROTOCOL.json'
    if not protocol_path.exists():
        protocol=dict(version=1,created_utc=datetime.now(timezone.utc).isoformat(),
            status='audit_derived_secondary_analysis_after_initial_primary_snapshot',
            membership_sha256=EXPECTED_MEMBERSHIP_SHA,
            cases=CASES,strata=STRATA,metrics=METRICS,bootstrap=2000,seed=20261024,
            cohort='All 1929 ATLAS records; official train/val/test/unlisted categories all reported. Main 3063-record cohort unchanged.',
            pairing='Full and each released comparator scored on the identical records within each published ID category.',
            eligibility='Only methods with all 3063 main-cohort records complete. Partial methods are disclosed and omitted until complete.',
            singleton_context='ATLAS records are single-chain; assembly and independent-chain BackFlip outputs are identical and collapsed here, without performance-based context selection.',
            unresolved='BF2 has no bound training ID manifest, so no inferred split. Non-ATLAS original chain-ID resolution is insufficient for these exact-match strata.',
            inference='Paired record-weighted family-cluster bootstrap, 2000 resamples; no independent-family interval for a one-family stratum.',
            limits='ID membership does not prove historical SGD inclusion, identical MD trajectories or absence of homolog/pretraining exposure. BF1 candidate official IDs are weaker than original checkpoint CSVs. Different comparator test IDs must not be used for cross-model test ranking.')
        protocol_path.write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    protocol=json.loads(protocol_path.read_text(encoding='utf-8'))
    assert protocol['membership_sha256']==EXPECTED_MEMBERSHIP_SHA
    assert protocol['strata']==STRATA and protocol['metrics']==METRICS
    # Performance is first read only after membership and the secondary rule are frozen.
    membership=pd.read_csv(membership_path,keep_default_na=False)
    audit_path=OUT/'scoring_audit.json'
    score_path=OUT/'per_record_metrics.csv'
    audit=json.loads(audit_path.read_text(encoding='utf-8'))
    per=pd.read_csv(score_path,keep_default_na=False)
    complete=set(audit['complete_methods'])
    own=membership[membership.method=='Full'].set_index('sample_id')
    assert len(own)==3063 and own.index.is_unique
    full=per[per.method=='Full'].set_index('sample_id')
    assert len(full)==3063 and full.index.is_unique
    rows=[];eligibility=[];rng=np.random.default_rng(protocol['seed'])
    for identifier_method,score_method in CASES:
        eligible=score_method in complete
        eligibility.append(dict(identifier_method=identifier_method,score_method=score_method,
            eligible=eligible,reason='full_common_cohort_complete' if eligible else 'main_cohort_incomplete_not_used'))
        if not eligible: continue
        members=membership[(membership.method==identifier_method)&(membership.dataset=='atlas')].copy()
        assert len(members)==1929 and members.sample_id.is_unique
        assert set(members.documented_split)=={'train','val','test','unlisted'}
        assert members.binding_strength.nunique()==1
        comparison=per[per.method==score_method].set_index('sample_id')
        assert len(comparison)==3063 and comparison.index.is_unique
        if identifier_method.startswith('backflip-'):
            chain_method=score_method.replace('(assembly)','(chains)')
            assert chain_method in complete
            chained=per[per.method==chain_method].set_index('sample_id')
            for metric in METRICS:
                assert np.array_equal(comparison.loc[members.sample_id,metric].to_numpy(),chained.loc[members.sample_id,metric].to_numpy())
        for stratum in STRATA:
            chosen=members if stratum=='all_atlas' else members[members.documented_split==stratum]
            ids=chosen.sample_id.tolist()
            left=full.loc[ids];right=comparison.loc[ids]
            assert left.family.tolist()==right.family.tolist()
            related=sum(own.loc[ids,'own_same_family_in_train'].astype(str)=='True')
            for metric in METRICS:
                a=left[metric].to_numpy(float);b=right[metric].to_numpy(float)
                assert len(a)==len(chosen) and np.isfinite(a).all() and np.isfinite(b).all()
                difference=pd.DataFrame(dict(family=left.family.to_numpy(),delta=a-b))
                family=difference.groupby('family').delta.agg(['sum','count'])
                lo=hi=np.nan
                if len(family)>1:
                    index=rng.integers(0,len(family),(protocol['bootstrap'],len(family)))
                    draws=family['sum'].to_numpy()[index].sum(1)/family['count'].to_numpy()[index].sum(1)
                    lo,hi=np.quantile(draws,[.025,.975])
                rows.append(dict(identifier_method=identifier_method,comparator=score_method,
                    dataset='atlas',stratum=stratum,binding_strength=chosen.binding_strength.iloc[0],
                    n_records=len(chosen),n_families=len(family),full_same_family_in_training_records=int(related),
                    metric=metric,full_mean=float(a.mean()),comparator_mean=float(b.mean()),
                    full_minus_comparator=float((a-b).mean()),ci_low=lo,ci_high=hi,
                    direction='positive_favors_Full' if metric in ['spearman','top10_ap'] else 'negative_favors_Full'))
    result=pd.DataFrame(rows)
    result.to_csv(OUT/'training_exposure_sensitivity.csv',index=False)
    notes=['# 公开训练ID分组的补充敏感性分析','',
        '这是根据公开名单审计派生的次级分析，不改变3,063条主队列或主榜。分组文件及分析规则先冻结，再读取分数。本项仅包含1,929条ATLAS单链记录，并列报告train、val、test、unlisted全部类别；每次比较的Full和外部模型使用完全相同的记录。','',
        '公开test标识不等于与所有预训练数据或同源家族独立。各比较模型的test名单不同，不能据此给模型间排test总榜。BF1系列原checkpoint CSV缺失，此处名单是分割名称与两个官方库一致的候选ID，证据强度单列。BF2没有可绑定的训练ID清单，不猜测其分组。','',
        'Full已确认同记录不在该折训练集，但原折划分不等于序列家族留出；下表同时列出Full在该组中有同家族其他样本参与训练的记录数。家族bootstrap只处理抽样依赖，不会改变训练隔离条件。','']
    for case in eligibility:
        name=case['score_method']
        if not case['eligible']:
            notes += [f'{name}仍未达全队列完整覆盖，本次不计入补充比较。',''];continue
        notes += [f'## {name}','']
        if case['identifier_method'].startswith('backflip-'):
            notes += ['ATLAS均为单链，BackFlip两种上下文的输出已验证逐值一致，此处合并为一份比较。','']
        notes += [
            '| 公开ID组 | n | Full ρ | 对照 ρ | 配对Δρ（95%区间） | Full MAE Å | 对照 MAE Å | Full同家族在训练 |',
            '|---|---:|---:|---:|---|---:|---:|---:|']
        for stratum in STRATA:
            part=result[(result.comparator==name)&(result.stratum==stratum)].set_index('metric')
            rank=part.loc['spearman'];error=part.loc['mae']
            interval=f"{rank.full_minus_comparator:+.4f} [{rank.ci_low:+.4f}, {rank.ci_high:+.4f}]" if np.isfinite(rank.ci_low) else f'{rank.full_minus_comparator:+.4f} [descriptive]'
            notes.append(f'| {stratum} | {int(rank.n_records)} | {rank.full_mean:.3f} | {rank.comparator_mean:.3f} | {interval} | {error.full_mean:.3f} | {error.comparator_mean:.3f} | {int(rank.full_same_family_in_training_records)} |')
        notes += ['','AP和标准化MAE及其全部配对区间见evaluation/training_exposure_sensitivity.csv。','']
    (ROOT/'TRAINING_EXPOSURE_SENSITIVITY_zh.md').write_text('\n'.join(notes)+'\n',encoding='utf-8')
    output_audit=dict(generated_utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=sha(protocol_path),
        membership_sha256=EXPECTED_MEMBERSHIP_SHA,source_scoring_audit_sha256=sha(audit_path),
        source_per_record_metrics_sha256=sha(score_path),script_sha256=sha(__file__),
        eligibility=eligibility,output_rows=len(result),all_published_id_categories_reported=True,
        main_cohort_or_ranking_changed=False,no_refit_no_new_inference=True,
        result_sha256=sha(OUT/'training_exposure_sensitivity.csv'))
    (OUT/'training_exposure_sensitivity_audit.json').write_text(json.dumps(output_audit,indent=2),encoding='utf-8')
    print(json.dumps(dict(rows=len(result),eligibility=eligibility)),flush=True)

if __name__=='__main__':main()
