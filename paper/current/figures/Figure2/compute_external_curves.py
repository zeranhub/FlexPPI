"""Reproduce external-method PR/ROC on exactly the frozen Figure 2d cohort.

No model fitting, score calibration, label modification, or source writes.
Run with the benchmark's runtime_py312/Scripts/python.exe.
"""
import os
for thread_variable in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[thread_variable] = '1'

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(r'E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023')
COMPLEX = ROOT / 'protein_complex_benchmark'
OUT = Path(__file__).resolve().parent
TRUTH = ROOT.parent / 'baselines/cohort/truth_residues.csv.gz'
PREDICTIONS = ROOT / 'evaluation/aligned_predictions.npz'
MEMBERSHIP = COMPLEX / 'classification/membership_frozen.csv'
INTERFACES = COMPLEX / 'classification/ppi_interfaces_20261024.csv'
PER_RECORD = COMPLEX / 'evaluation/per_record_metrics.csv'
METHODS = ['Full', 'backflip-2.1 (assembly)', 'backflip-1.0 (assembly)',
           'DynaProt-M', 'FastProtFlex', 'ANM', 'GNM']
SEED = 20261024
BOOTSTRAP = 2000
GRID = np.linspace(0., 1., 1001)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path, **kwargs):
    return pd.read_csv(path, encoding='utf-8-sig', float_precision='round_trip', **kwargs)


def threshold_group_curves(labels, scores):
    """Descending score groups; tied predictions enter together at one threshold."""
    order = np.argsort(-scores, kind='stable')
    ranked_scores = scores[order]
    ranked_labels = labels[order].astype(np.int64)
    group_ends = np.r_[np.flatnonzero(np.diff(ranked_scores) != 0), len(scores) - 1]
    true_positives = np.cumsum(ranked_labels)[group_ends]
    selected = group_ends + 1
    false_positives = selected - true_positives
    n_positive = int(labels.sum())
    n_negative = len(labels) - n_positive
    recall = true_positives / n_positive
    precision = true_positives / selected
    fpr = false_positives / n_negative
    ap = float(np.sum(np.diff(np.r_[0., recall]) * precision))
    # The first threshold group reaching the requested recall is retained;
    # no precision envelope, linear PR interpolation, or curve integration.
    first = np.searchsorted(recall, GRID, side='left')
    pr_grid = precision[np.minimum(first, len(precision) - 1)]
    # Keep the original order of vertical segments. Collapsing a duplicate FPR
    # to its upper TPR before interpolation incorrectly raises the preceding
    # interval. np.interp preserves the lower arrival and upper departure.
    raw_fpr = np.r_[0., fpr]
    raw_tpr = np.r_[0., recall]
    roc_grid = np.interp(GRID, raw_fpr, raw_tpr)
    auroc = float(roc_auc_score(labels, scores))
    raw_trapezoid_auc = float(np.trapz(raw_tpr, raw_fpr))
    assert abs(raw_trapezoid_auc - auroc) <= 2e-14
    checked_ap = float(average_precision_score(labels, scores))
    assert abs(checked_ap - ap) <= 2e-14
    return (pr_grid, roc_grid, ap, auroc,
            dict(threshold=ranked_scores[group_ends], precision=precision,
                 recall=recall, fpr=fpr))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # A vertical rise at FPR=1 must not raise the preceding horizontal segment.
    toy_fpr = np.array([0., 0., 1., 1.])
    toy_tpr = np.array([0., .5, .5, 1.])
    assert np.interp(.5, toy_fpr, toy_tpr) == .5
    assert np.trapz(toy_tpr, toy_fpr) == .5
    scoring_audit = json.loads((ROOT / 'evaluation/scoring_audit.json').read_text(encoding='utf-8'))
    assert sha(TRUTH) == scoring_audit['truth_sha256'], 'Frozen truth hash mismatch'
    members = read_csv(MEMBERSHIP)
    interfaces = read_csv(INTERFACES)
    baseline = read_csv(PER_RECORD)
    selected = baseline[baseline.method.isin(METHODS)].merge(
        members[['sample_id', 'actual_input_n_chains', 'benchmark_group']],
        on='sample_id', how='left', validate='many_to_one')
    selected = selected[(selected.dataset == 'dynarepo') &
        (selected.benchmark_group == 'protein_complex') &
        (selected.actual_input_n_chains >= 2)].merge(
            interfaces[['sample_id', 'n_interfaces']], on='sample_id', validate='many_to_one')
    selected = selected[selected.n_interfaces >= 1].copy()
    assert not selected.duplicated(['method', 'sample_id']).any()
    ids = sorted(selected.sample_id.unique())
    assert len(ids) == 615 and len(selected) == 615 * len(METHODS)
    for method in METHODS:
        assert sorted(selected.loc[selected.method == method, 'sample_id']) == ids
    cohort = selected[selected.method == 'Full'][['sample_id', 'dataset', 'family',
        'n', 'actual_input_n_chains', 'n_interfaces']].set_index('sample_id').loc[ids].reset_index()
    ni = cohort.n_interfaces
    cohort['stratum'] = np.select([ni == 1, ni.between(2, 3), ni >= 4], ['i1', 'i23', 'i4'], default='')
    strata = cohort.stratum.value_counts().to_dict()
    assert strata == {'i1': 414, 'i23': 144, 'i4': 57}, strata
    truth = read_csv(TRUTH, keep_default_na=False,
        usecols=['sample_id', 'residue_index', 'chain', 'author_residue_id', 'target_rmsf', 'valid_rmsf'],
        dtype={'chain': str, 'author_residue_id': str})
    assert len(truth) == 961321
    all_rows = {sample: rows for sample, rows in truth[truth.sample_id.isin(ids)].groupby('sample_id', sort=False)}
    with np.load(PREDICTIONS) as archive:
        assert all(method in archive.files for method in METHODS)
        predictions = {method: archive[method] for method in METHODS}
    for values in predictions.values():
        assert values.shape == (len(truth),)
    print('Loaded frozen arrays and confirmed 615 records / 414,144,57 interface strata.', flush=True)
    n_record = len(ids)
    pr = {method: np.empty((n_record, len(GRID))) for method in METHODS}
    roc = {method: np.empty((n_record, len(GRID))) for method in METHODS}
    metrics = []
    raw_parts = {key: [] for key in ['threshold', 'precision', 'recall', 'fpr']}
    raw_offsets = [0]
    raw_samples, raw_methods = [], []
    prevalence = []
    positivity_cutoffs = []
    reference = selected.set_index(['method', 'sample_id'])
    invalid_truth = 0
    missing = []
    for record_index, sample in enumerate(ids):
        rows = all_rows[sample]
        valid = rows.valid_rmsf.to_numpy(bool)
        positions = rows.index.to_numpy()[valid]
        y = rows.target_rmsf.to_numpy(float)[valid]
        invalid_truth += int((~np.isfinite(y)).sum())
        assert np.isfinite(y).all() and np.ptp(y) > 0
        assert len(y) == int(cohort.loc[record_index, 'n'])
        finite = np.isfinite(np.column_stack([predictions[method][positions] for method in METHODS]))
        if not finite.all():
            for method_index, method in enumerate(METHODS):
                count = int((~finite[:, method_index]).sum())
                if count:
                    missing.append(dict(sample_id=sample, method=method, missing_valid_residues=count))
            continue
        cutoff = float(np.quantile(y, .90))
        labels = y >= cutoff
        assert 0 < labels.sum() < len(labels)
        prevalence.append(float(labels.mean()))
        positivity_cutoffs.append(cutoff)
        for method in METHODS:
            pg, rg, ap, auroc, raw = threshold_group_curves(labels, predictions[method][positions])
            pr[method][record_index] = pg
            roc[method][record_index] = rg
            old_ap = float(reference.loc[(method, sample), 'top10_ap'])
            assert abs(old_ap - ap) <= 2e-13, (method, sample, old_ap, ap)
            metrics.append(dict(sample_id=sample, family=cohort.loc[record_index, 'family'], method=method,
                n_residues=len(y), n_positive=int(labels.sum()), prevalence=float(labels.mean()),
                rmsf_positive_cutoff=cutoff, AP=ap, AUROC=auroc,
                panel_d_AP=old_ap, AP_difference=ap-old_ap, n_score_groups=len(raw['threshold'])))
            for key in raw_parts:
                raw_parts[key].append(raw[key])
            raw_offsets.append(raw_offsets[-1] + len(raw['threshold']))
            raw_samples.append(sample)
            raw_methods.append(method)
    if missing:
        (OUT / 'missing_prediction_audit.json').write_text(json.dumps(missing, indent=2), encoding='utf-8')
        raise RuntimeError('Missing valid-residue predictions; no residues silently dropped. See missing_prediction_audit.json')
    assert invalid_truth == 0 and len(metrics) == n_record * len(METHODS)
    per_record = pd.DataFrame(metrics)
    per_record.to_csv(OUT / 'per_record_metrics.csv', index=False)
    cohort['n_positive'] = per_record[per_record.method == 'Full'].n_positive.to_numpy()
    cohort['prevalence'] = prevalence
    cohort['rmsf_positive_cutoff'] = positivity_cutoffs
    cohort.to_csv(OUT / 'cohort_records.csv', index=False)
    # Resample family clusters, preserving every record in each sampled family.
    names = sorted(cohort.family.unique())
    family_positions = {family: index for index, family in enumerate(names)}
    codes = np.array([family_positions[family] for family in cohort.family], dtype=np.int32)
    sizes = np.bincount(codes, minlength=len(names))
    rng = np.random.default_rng(SEED)
    weights = rng.multinomial(len(names), np.full(len(names), 1 / len(names)), size=BOOTSTRAP)
    denominators = weights @ sizes

    def bootstrap_record_mean(values):
        if values.ndim == 1:
            sums = np.bincount(codes, weights=values, minlength=len(names))
            return (weights @ sums) / denominators
        sums = np.zeros((len(names), values.shape[1]))
        np.add.at(sums, codes, values)
        return (weights @ sums) / denominators[:, None]

    summary = []
    for method in METHODS:
        block = per_record[per_record.method == method].set_index('sample_id').loc[ids]
        for metric in ['AP', 'AUROC']:
            values = block[metric].to_numpy(float)
            draws = bootstrap_record_mean(values)
            lo, hi = np.quantile(draws, [.025, .975])
            summary.append(dict(scope='all', method=method, metric=metric, mean=float(values.mean()),
                ci_low=float(lo), ci_high=float(hi), n_records=n_record, n_families=len(names)))
    pd.DataFrame(summary).to_csv(OUT / 'summary_metrics.csv', index=False)
    curves = []
    for curve_type, values_by_method in [('PR', pr), ('ROC', roc)]:
        for method, values in values_by_method.items():
            lo = hi = np.full(len(GRID), np.nan)
            if method == 'Full':
                draws = bootstrap_record_mean(values)
                lo, hi = np.quantile(draws, [.025, .975], axis=0)
            curves.append(pd.DataFrame(dict(scope='all', method=method, curve_type=curve_type,
                x=GRID, mean=values.mean(axis=0), ci_low=lo, ci_high=hi)))
    pd.concat(curves, ignore_index=True).to_csv(OUT / 'mean_curves.csv', index=False)
    np.savez_compressed(OUT / 'raw_threshold_curves.npz',
        offsets=np.array(raw_offsets, dtype=np.int64), sample_id=np.array(raw_samples), method=np.array(raw_methods),
        **{key: np.concatenate(parts) for key, parts in raw_parts.items()})
    endpoint = float(np.mean(prevalence))
    source_paths = [TRUTH, PREDICTIONS, MEMBERSHIP, INTERFACES, PER_RECORD,
        ROOT / 'evaluate_benchmark.py', COMPLEX / 'aggregate_benchmark.py',
        ROOT / 'evaluation/scoring_audit.json', COMPLEX / 'evaluation/scoring_audit.json']
    audit = dict(generated_utc=datetime.now(timezone.utc).isoformat(), status='complete',
        script_sha256=sha(__file__), sources=[dict(path=str(path), sha256=sha(path)) for path in source_paths],
        cohort='exact Figure 2d cohort: dynarepo; protein_complex; >=2 actual chains; >=1 observed protein interface',
        n_records=n_record, n_families=len(names), n_residues=int(cohort.n.sum()), strata=strata, methods=METHODS,
        invalid_reference_residues=0, missing_valid_predictions=[], discarded_valid_residues=0,
        shared_residue_policy='All frozen valid RMSF residues, identical for every method; all seven predictions must be finite, otherwise fail.',
        label_policy='Per-record reference RMSF >= np.quantile(reference RMSF, 0.90); tied reference values retained.',
        threshold_policy='Descending unique prediction score groups; every tied score group enters together.',
        mean_PR='Per-record first threshold group reaching each common recall, then equal-record mean. r=0 uses first threshold group precision, not an artificial point. r=1 retains first reach; all-selected endpoint is separate.',
        all_selected_precision=endpoint, prevalence_min=float(min(prevalence)), prevalence_max=float(max(prevalence)),
        mean_ROC='Per-record np.interp on original nondecreasing raw FPR/TPR points, preserving duplicate FPR and vertical segment order; then equal-record mean. At exact duplicate FPR, np.interp returns the last upper point; interpolation between distinct FPR values retains arrival at the next vertical lower endpoint. (0,0) may be separately plotted as vertical start.',
        ROC_interpolation_validation='Toy raw ROC (0,0),(0,0.5),(1,0.5),(1,1): interpolated TPR at FPR=0.5 equals 0.5 and raw trapezoid area equals 0.5. Each real raw-threshold trapezoid AUC agrees with sklearn ROC AUC within 2e-14. Duplicate FPR must not be collapsed before interpolation.',
        AP_policy='Sum of recall increments times threshold-group precision, computed on original thresholds; independently checked against sklearn.average_precision_score. Never integral of displayed mean curve.',
        AUROC_policy='Per-record sklearn.roc_auc_score, then equal-record mean; not integral of displayed mean curve.',
        max_abs_AP_difference_from_panel_d=float(per_record.AP_difference.abs().max()),
        bootstrap=dict(seed=SEED, draws=BOOTSTRAP, level='family cluster', percentile=[2.5, 97.5],
            shared_draws=True, denominator='total records in sampled families', pointwise_curve_bands='Full only; not simultaneous confidence bands',
            family_order_sha256=hashlib.sha256('\n'.join(names).encode()).hexdigest(),
            resample_weights_sha256=hashlib.sha256(weights.astype('<i4').tobytes()).hexdigest()),
        raw_curve_format='raw_threshold_curves.npz: one segment per sample_id/method; offsets length = segments+1; thresholds descending; recall=TPR; all-selected endpoint retained; ROC (0,0) not stored as an artificial threshold.',
        interpretation_limits=scoring_audit.get('interpretation_limits', [
            'Published model weights can have mixed training exposure and different target definitions; this is a frozen-cohort comparison, not an equal-training or universal-superiority claim.']))
    (OUT / 'calculation_audit.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding='utf-8')
    (OUT / 'calculation_notes_zh.md').write_text(
        '# 图2外部模型 PR 与 ROC 数据说明\n\n'
        f'沿用图2d的615条DynaRepo复合物记录，涉及{len(names)}个family；界面分层为414、144、57条。'
        f'共{int(cohort.n.sum()):,}个有效RMSF残基，7种方法均完整覆盖，没有删除缺失预测的残基。'
        '每条记录把参考RMSF不低于其90%分位点的残基记为阳性，保留并列值，预测分数相同的残基整组进入。'
        f'实际阳性率的逐记录平均值为{endpoint:.8f}，可作为PR随机参考线。\n\n'
        'PR在共同召回率上取该记录首次达到该召回率的阈值组精确率，然后对记录等权平均。'
        '召回率0处也用首个阈值组精确率；召回率1处保留首次达到1的精确率，并在绘图时额外连到全部残基入选的末端阳性率。'
        'AP来自原始阈值组，不是显示曲线的积分。ROC对每条记录保留原始阈值顺序和重复的假阳性率后线性插值，再等权平均；这样保留垂直段到达下端和离开上端的次序，不会提前抬高前一段曲线。绘图可补上原点到x=0上端的竖线。'
        'AUROC先逐记录计算再取均值。置信区间按family整簇重抽2000次，seed=20261024；Full曲线带为逐点95%区间，不是同时置信带。\n\n'
        f'本次重算AP与图2d逐记录AP最大绝对差为{audit["max_abs_AP_difference_from_panel_d"]:.3g}。'
        'raw_threshold_curves.npz保留每条记录、每种方法全部原始阈值组，calculation_audit.json记录来源哈希和精确规则。'
        '源预测、标签和稿件均未修改。\n', encoding='utf-8')
    print(pd.DataFrame(summary).to_string(index=False), flush=True)
    print(json.dumps({key: audit[key] for key in ['n_records', 'n_families', 'n_residues',
        'all_selected_precision', 'max_abs_AP_difference_from_panel_d']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
