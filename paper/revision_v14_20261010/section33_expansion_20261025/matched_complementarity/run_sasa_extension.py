"""Quality-filtered, separately frozen SASA extension to fixed matching rules."""
from pathlib import Path
import hashlib
import json
import importlib.util
from datetime import datetime, timezone
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('matched_base', HERE/'analyze_matched_complementarity.py')
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
FEATURES = HERE.parent/'feature_baseline/features_priority83.parquet'
OUT = HERE/'sasa_extension'
RULE = {
    'version': 'sasa_extension_quality_v1',
    'parent_protocol_sha256': BASE.sha(HERE/'protocol.json'),
    'quality': 'Require finite sidechain_sasa_complex_A2, sidechain_complete=True, chi_input_complete=True, nonstandard_residue=False. Identity sample_id/residue_index/amino_acid must match OOF. Quality exclusion is defined before reading feature values.',
    'matching': 'Same record/same amino-acid group >=4 after the quality filter; average rank/n percentiles recomputed in the qualified group. SASA percentile difference<=0.25 plus truth RMSF percentile difference<=0.25 or physical RMSF difference<=0.10 A. Enumerate every i<j pair. Target Hchi; all_close primary and target-percentile-separation>=0.50 secondary. Same aggregation and family-cluster bootstrap as frozen base protocol.',
    'truth_modes': list(BASE.MODES),
    'primary': 'uniform_rep1',
    'caveat': 'Explanatory matching on MD RMSF. All-atom input SASA is not identical to the backbone-only model input, and ranks do not imply equality of absolute SASA; physical SASA balance is reported.',
}


def main():
    OUT.mkdir(exist_ok=True)
    plan_path = OUT/'sasa_extension_protocol.json'
    if plan_path.exists():
        assert json.loads(plan_path.read_text(encoding='utf-8'))['fixed_rules'] == RULE
    else:
        BASE.write_json(plan_path, {'created_utc': datetime.now(timezone.utc).isoformat(),
                                   'fixed_rules': RULE, 'code_sha256': BASE.sha(__file__)})
    # The extension's quality/matching rules above are on disk before feature IO.
    raw = BASE.freeze()
    positions = BASE.assemble(raw)
    f = pd.read_parquet(FEATURES)
    assert f.sample_id.nunique() == 83
    assert not f.duplicated(['sample_id', 'residue_index']).any()
    required = ['sample_id', 'residue_index', 'amino_acid', 'sidechain_sasa_complex_A2',
                'sidechain_complete', 'chi_input_complete', 'nonstandard_residue']
    for col in required:
        assert col in f, col
    take = f[required].rename(columns={'amino_acid': 'feature_AA',
                                      'sidechain_sasa_complex_A2': 'sasa_sidechain_complex_A2'})
    merged = positions.merge(take, on=['sample_id', 'residue_index'], how='left', validate='one_to_one')
    assert merged.feature_AA.notna().all()
    assert (merged.amino_acid == merged.feature_AA).all()
    merged['sasa_finite'] = np.isfinite(merged.sasa_sidechain_complex_A2)
    merged['quality_ok'] = (merged.sasa_finite & merged.sidechain_complete.eq(True) &
                            merged.chi_input_complete.eq(True) & merged.nonstandard_residue.eq(False))
    audit = merged.groupby(['sample_id', 'dataset']).agg(
        common_positions=('residue_index', 'size'), finite_sasa=('sasa_finite', 'sum'),
        sidechain_complete=('sidechain_complete', 'sum'), chi_input_complete=('chi_input_complete', 'sum'),
        nonstandard=('nonstandard_residue', 'sum'), quality_ok=('quality_ok', 'sum')).reset_index()
    BASE.csv(audit, OUT/'sasa_quality_audit_by_record.csv')
    BASE.csv(merged[~merged.quality_ok], OUT/'sasa_quality_excluded_positions.csv')
    finite = merged[merged.quality_ok].copy()
    finite['sameAA_n_quality'] = finite.groupby(['sample_id', 'amino_acid']).residue_index.transform('size')
    BASE.csv(finite, OUT/'common_positions_quality_sasa.csv.gz')
    manifest = dict(path=str(FEATURES), bytes=FEATURES.stat().st_size, sha256=BASE.sha(FEATURES),
                    sasa_column='sidechain_sasa_complex_A2', feature_records=f.sample_id.nunique(),
                    feature_positions=len(f), original_common_positions=len(positions),
                    quality_qualified_positions=len(finite),
                    quality_sameAA4_positions=int((finite.sameAA_n_quality >= 4).sum()),
                    quality_records=finite.sample_id.nunique())
    BASE.write_json(OUT/'sasa_input_manifest.json', manifest)
    print('Quality/SASA identity verified; running all eligible matched pairs.', flush=True)
    per_record, summary = BASE.matched_analysis(finite, OUT, sasa=True)
    primary = summary[(summary.label_mode == 'uniform_rep1')&
                      (summary.target_subset == 'all_close')&(summary.dataset_scope == 'all')]
    BASE.csv(primary, OUT/'primary_matched_summary.csv')
    # Audit physical SASA differences in the exact primary pairs; no new selection.
    pairs = pd.read_csv(OUT/'uniform_rep1_matched_pairs.csv.gz')
    lookup = finite.set_index(['sample_id', 'residue_index']).sasa_sidechain_complex_A2
    ki = pd.MultiIndex.from_frame(pairs[['sample_id', 'residue_index_i']].rename(columns={'residue_index_i': 'residue_index'}))
    kj = pd.MultiIndex.from_frame(pairs[['sample_id', 'residue_index_j']].rename(columns={'residue_index_j': 'residue_index'}))
    pairs['sasa_abs_delta_A2'] = np.abs(lookup.reindex(kj).to_numpy()-lookup.reindex(ki).to_numpy())
    balance = []
    for rule in ('sasa_q25_and_rmsf_q25', 'sasa_q25_and_rmsf_physical_0p10'):
        cond = np.abs(pairs.rmsf_truth_percentile_delta_j_minus_i) <= .25+BASE.EPS if rule.endswith('rmsf_q25') else np.abs(pairs.rmsf_truth_delta_j_minus_i) <= .10+BASE.EPS
        sel = pairs[cond]
        by_type = sel.groupby(['sample_id', 'dataset', 'amino_acid']).agg(
            matched_pairs=('sasa_abs_delta_A2', 'size'), sasa_abs_delta_A2=('sasa_abs_delta_A2', 'mean')).reset_index()
        by_record = by_type.groupby(['sample_id', 'dataset']).agg(
            matched_types=('amino_acid', 'size'), matched_pairs=('matched_pairs', 'sum'),
            sasa_abs_delta_A2=('sasa_abs_delta_A2', 'mean')).reset_index()
        by_record['match_rule'] = rule
        balance.append(by_record)
    BASE.csv(pd.concat(balance, ignore_index=True), OUT/'primary_physical_sasa_balance_by_record.csv')
    lines = ['# SASA与RMSF同时相近条件下的侧链差异', '',
             '质量与匹配规则见sasa_extension_protocol.json；规则在读取特征值前冻结。主结果只使用初始结构侧链完整、χ定义所需原子完整、标准残基且SASA有限的位置。每记录同类型合格位置少于4不进入匹配。SASA相近指同类型百分位差不超过0.25；真实RMSF另采用百分位差不超过0.25或绝对差不超过0.10 Å。匹配选择不依赖模型预测。', '',
             '| 规则 | 指标 | 记录宏平均 [整簇95%区间] | 记录数 | 家族簇数 | 去最大簇 |',
             '|---|---|---:|---:|---:|---:|']
    for r in primary[primary.metric.isin(['target_true_abs_delta_mean', 'full_direction_accuracy', 'surrogate_direction_accuracy', 'full_minus_surrogate_accuracy', 'intra_direction_accuracy', 'full_minus_intra_accuracy'])].itertuples():
        lines.append(f'| {r.match_rule} | {r.metric} | {r.record_macro_mean:.4f} [{r.family_cluster_ci95_low:.4f}, {r.family_cluster_ci95_high:.4f}] | {r.defined_records} | {r.defined_family_clusters} | {r.leave_largest_family_out_mean:.4f} |')
    lines += ['', '该分析可回答：在同一蛋白中、同一种氨基酸、暴露程度和主链活动相近时，模型是否仍能区分侧链占据多样性。它不是去除所有结构混杂的因果证明，也不能将使用真实MD RMSF的分层伪装成部署条件。SASA百分位相近不意味着绝对面积完全相同，故另给primary_physical_sasa_balance_by_record.csv；若残留差较大，仍需结合强结构基线评价。']
    (OUT/'reviewer_interpretation_zh.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    BASE.write_json(OUT/'completion_audit.json', dict(status='complete', input_manifest=manifest,
        protocol_sha256=BASE.sha(plan_path), executed_code_sha256=BASE.sha(__file__),
        base_code_sha256=BASE.sha(HERE/'analyze_matched_complementarity.py'),
        per_record_rows=len(per_record), summary_rows=len(summary), source_files_modified=False,
        word_modified=False, new_training=False, new_inference=False, new_MD=False))
    print(json.dumps(BASE.clean(manifest), ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
