"""Fixed-rule, exhaustive within-record/AA complementarity analysis.

Original files, labels, predictions and the Word manuscript are read-only.
Run --freeze-only before analysis; then run without arguments. Optional --sasa
adds a separately frozen extension after the other branch finishes its features.
Truth-conditioned matching is explanatory, never a deployable predictor.
"""
from __future__ import annotations
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
import gzip
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

ROOT = Path(r'E:\Project\FlexPPI')
HERE = Path(__file__).resolve().parent
OCC = ROOT/'experiments/_sources/model_argument_v24_20261013/experiment_predictive_occupancy'
ENERGY = ROOT/'experiments/_sources/energy_landscape_v21_20261012'
PRED = ROOT/'experiments/core_corrected_pipeline/neural_exports/evaluation'
LABEL = ROOT/'experiments/_sources/label_sensitivity/labels'
FAMILY = ROOT/'experiments/_sources/family_folds/family_components.json'
SOURCES = ('atlas', 'dynarepo', 'itrace')
MODES = ('uniform_rep1', 'uniform_available_replicas', 'original_label')
CLASSES = ('low_RMSF_low_Hchi', 'low_RMSF_high_Hchi',
           'high_RMSF_low_Hchi', 'high_RMSF_high_Hchi', 'unclassified')
B = 10000
SEED = 20261025
EPS = 1e-12

RULES = {
    'version': 'section33_matched_complementarity_v1',
    'primary_cohort': 'Original preselected 83 records, exact original 33493 common positions; same-record same-AA groups >=4 give 33142 positions.',
    'truth': 'Primary: original cached 0-100 ns/200 ps/501-frame labels uniform_rep1. Sensitivity: equal available-replica mean of the same window; original training/reference label. Predictions are the unchanged original OOF predictions.',
    'percentiles': 'Within the exact common valid positions of a record/AA group, average ranks divided by n; ties share their average rank. Matching uses no prediction.',
    'matching': 'Enumerate every unordered pair i<j within each eligible record/AA group. Rank-close rule: absolute truth-conditioning percentile difference <=0.25. Physical rules: absolute truth-conditioning difference <=0.10 primary, <=0.05 and <=0.20 sensitivities. Apply both directions (condition on RMSF, target Hchi; condition on Hchi, target RMSF).',
    'target_subsets': 'all_close contains all matched pairs; separated_target_p50 additionally requires absolute target truth percentile difference >=0.50. This is an explanatory truth-separated stratum, not prediction-led selection.',
    'metrics': 'Actual target and conditioning differences in Angstrom or nats; target pair-order concordance for dedicated Full output, cross-output surrogate, and dedicated Intra output; paired Full-minus-surrogate and Full-minus-Intra differences; dedicated-output difference MAE. True ties are excluded only from ordering scores; predicted ties score 0.5.',
    'aggregation': 'First all pairs equally within a record/AA group; then eligible AA groups equally within each record; then records equally. Counts are descriptive, never the inferential sample size. Accuracy denominators exclude true ties. Undefined groups/records remain in audit tables.',
    'resampling': '10000 percentile bootstrap draws, entire available family_components clusters sampled with replacement. Overall cluster sampling preserves cross-source links; source-specific analyses resample clusters represented in that source. Chunked computation. Report family-equal and leave-largest-family-out descriptive sensitivity.',
    'family_caveat': 'Cluster resampling addresses dependence; it does not convert original OOF into family-held-out predictions or remove train/test family overlap.',
    'full_cohort_review': 'Recompute previous within-record/AA inclusive 25th/75th-quartile five-way classes, confusion and per-record recall on all 3063 original OOF records/824988 positions. Require >=4 and q25<q75 per output. Collapsed and middle groups are unclassified. Verify previous labels and per-record recall numerically.',
    'sasa_extension': 'Predeclared optional extension only: complex sidechain SASA average-rank percentile difference <=0.25 in addition to RMSF percentile-close <=0.25 or physical RMSF-close <=0.10 A. Same-AA counts/ranks use finite SASA positions; missingness audited. Analyze Hchi target for all_close and separated_target_p50; no outcome or prediction-based case picking.',
    'interpretation': 'MD-conditioned matching is explanatory. Relative quartiles are within record/type; they are not absolute physical high/low states. Hchi is coarse joint angular occupancy entropy and not full thermodynamic entropy. Consecutive windows are not independent replicas.',
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(4*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean(v) for v in value]
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(clean(value), ensure_ascii=False, indent=2,
                                    allow_nan=False)+'\n', encoding='utf-8')


def csv(frame, path):
    frame.to_csv(path, index=False, float_format='%.12g', encoding='utf-8',
                 compression='gzip' if str(path).endswith('.gz') else None)


def mean(x):
    a = np.asarray(x, float)
    return float(a[np.isfinite(a)].mean()) if np.isfinite(a).any() else np.nan


def stable_seed(key):
    return (SEED + int(hashlib.sha256(str(key).encode()).hexdigest()[:8], 16)) % (2**32)


def families():
    components = json.loads(FAMILY.read_text(encoding='utf-8'))
    mapping = {}
    for c in components:
        for sid in c['members']:
            assert sid not in mapping, sid
            mapping[sid] = c['component_id']
    return mapping


def cohort_source():
    return pd.read_csv(OCC/'replica1_per_residue.csv.gz', low_memory=False,
                       dtype={'chain': str, 'author_residue_id': str})


def freeze():
    HERE.mkdir(parents=True, exist_ok=True)
    path = HERE/'protocol.json'
    if path.exists():
        old = json.loads(path.read_text(encoding='utf-8'))
        assert old['fixed_rules'] == RULES, 'Protocol rule changes require a new version.'
    else:
        write_json(path, {'created_utc': datetime.now(timezone.utc).isoformat(),
                          'fixed_rules': RULES, 'code_sha256_at_freeze': sha(__file__)})
    d = cohort_source()
    ids = sorted(d.sample_id.unique())
    assert len(ids) == 83 and len(d) == 33493
    paths = [OCC/'replica1_per_residue.csv.gz', FAMILY,
             ENERGY/'oof_analysis/paired_common_residues.csv.gz',
             ENERGY/'oof_analysis/per_record_reference_classes.csv',
             ENERGY/'oof_analysis/per_record_correlations.csv']
    for sid in ids:
        paths += [LABEL/mode/(sid+'.npz') for mode in MODES[:2]]
        paths += [PRED/method/(sid+'.npz') for method in ('full', 'intrachain')]
    manifest = [{'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in paths]
    manifest_path = HERE/'source_manifest.json'
    if manifest_path.exists():
        assert json.loads(manifest_path.read_text(encoding='utf-8'))['sources'] == manifest
    else:
        write_json(manifest_path, {'protocol_sha256': sha(path), 'sources': manifest})
    return d


def assemble(d):
    parts = []
    mapping = families()
    for sid, p in d.groupby('sample_id', sort=True):
        p = p.sort_values('residue_index').copy()
        idx = p.residue_index.to_numpy(int)
        f = np.load(PRED/'full'/(sid+'.npz'), allow_pickle=False)
        intra = np.load(PRED/'intrachain'/(sid+'.npz'), allow_pickle=False)
        for field in ('amino_acid', 'chain', 'author_residue_id'):
            np.testing.assert_array_equal(p[field].astype(str).to_numpy(), f[field][idx].astype(str))
            np.testing.assert_array_equal(f[field], intra[field])
        np.testing.assert_allclose(p.full_rmsf, f['mean'][idx, 0], rtol=0, atol=1e-7)
        np.testing.assert_allclose(p.full_hchi, f['mean'][idx, 1], rtol=0, atol=1e-7)
        np.testing.assert_allclose(p.intrachain_hchi, intra['mean'][idx, 1], rtol=0, atol=1e-7)
        p['pred_full_rmsf'] = f['mean'][idx, 0]
        p['pred_full_hchi'] = f['mean'][idx, 1]
        p['pred_intra_rmsf'] = intra['mean'][idx, 0]
        p['pred_intra_hchi'] = intra['mean'][idx, 1]
        p['family_id'] = mapping[sid]
        for mode in MODES[:2]:
            z = np.load(LABEL/mode/(sid+'.npz'), allow_pickle=False)
            assert z['mask'][idx, :2].all()
            assert np.isfinite(z['values'][idx, :2]).all()
            p[mode+'_rmsf'] = z['values'][idx, 0]
            p[mode+'_hchi'] = z['values'][idx, 1]
        p['original_label_rmsf'] = f['target'][idx, 0]
        p['original_label_hchi'] = f['target'][idx, 1]
        p['sameAA_n'] = p.groupby('amino_acid').residue_index.transform('size')
        p['eligible_sameAA'] = p.sameAA_n >= 4
        parts.append(p)
    result = pd.concat(parts, ignore_index=True)
    assert int(result.eligible_sameAA.sum()) == 33142
    return result


METRICS = ('target_true_abs_delta_mean', 'target_true_abs_delta_median',
           'target_percentile_abs_delta_mean', 'condition_true_abs_delta_mean',
           'condition_percentile_abs_delta_mean',
           'full_direction_accuracy', 'surrogate_direction_accuracy',
           'intra_direction_accuracy', 'full_minus_surrogate_accuracy',
           'full_minus_intra_accuracy', 'full_difference_mae',
           'intra_difference_mae', 'full_minus_intra_difference_mae',
           'full_pred_abs_delta_mean')
COUNT_COLS = ('eligible_positions', 'potential_pairs', 'matched_pairs',
              'matched_positions', 'non_tied_truth_pairs', 'full_pred_tied_pairs')


def percentile(a):
    return pd.Series(a).rank(method='average').to_numpy(float)/len(a)


def pair_stats(g, mode, include_sasa=False):
    # No prediction is used in pair eligibility or target-separation strata.
    n = len(g)
    i, j = np.triu_indices(n, 1)
    truth = {k: g[mode+'_'+k].to_numpy(float) for k in ('rmsf', 'hchi')}
    td = {k: a[j]-a[i] for k, a in truth.items()}
    ranks = {k: percentile(a) for k, a in truth.items()}
    rd = {k: a[j]-a[i] for k, a in ranks.items()}
    pred = {k: g['pred_'+k].to_numpy(float)[j]-g['pred_'+k].to_numpy(float)[i]
            for k in ('full_rmsf', 'full_hchi', 'intra_rmsf', 'intra_hchi')}
    meta = dict(sample_id=g.sample_id.iloc[0], dataset=g.dataset.iloc[0],
                fold=int(g.fold.iloc[0]), family_id=g.family_id.iloc[0],
                amino_acid=g.amino_acid.iloc[0], label_mode=mode)
    if include_sasa:
        sr = percentile(g.sasa_sidechain_complex_A2.to_numpy(float))
        sd = sr[j]-sr[i]
    rows = []
    for cond, target in (('rmsf', 'hchi'), ('hchi', 'rmsf')):
        if include_sasa and target != 'hchi':
            continue
        if include_sasa:
            masks = {'sasa_q25_and_rmsf_q25': (np.abs(sd) <= .25+EPS)&(np.abs(rd[cond]) <= .25+EPS),
                     'sasa_q25_and_rmsf_physical_0p10': (np.abs(sd) <= .25+EPS)&(np.abs(td[cond]) <= .10+EPS)}
        else:
            masks = {'percentile_0p25': np.abs(rd[cond]) <= .25+EPS}
            for caliper in (.05, .10, .20):
                masks['physical_'+format(caliper, '.2f').replace('.', 'p')] = np.abs(td[cond]) <= caliper+EPS
        for rule, close in masks.items():
            for subset in ('all_close', 'separated_target_p50'):
                keep = close.copy()
                if subset == 'separated_target_p50':
                    keep &= np.abs(rd[target]) >= .50-EPS
                count = int(keep.sum())
                r = meta | dict(direction=cond+'_to_'+target, match_rule=rule,
                                target_subset=subset, eligible_positions=n,
                                potential_pairs=len(i), matched_pairs=count,
                                matched_positions=len(np.unique(np.r_[i[keep], j[keep]])),
                                non_tied_truth_pairs=0, full_pred_tied_pairs=0)
                for metric in METRICS:
                    r[metric] = np.nan
                if count:
                    actual = td[target][keep]
                    r['target_true_abs_delta_mean'] = float(np.abs(actual).mean())
                    r['target_true_abs_delta_median'] = float(np.median(np.abs(actual)))
                    r['target_percentile_abs_delta_mean'] = float(np.abs(rd[target][keep]).mean())
                    r['condition_true_abs_delta_mean'] = float(np.abs(td[cond][keep]).mean())
                    r['condition_percentile_abs_delta_mean'] = float(np.abs(rd[cond][keep]).mean())
                    for model in ('full', 'intra'):
                        delta = pred[model+'_'+target][keep]
                        r[model+'_difference_mae'] = float(np.abs(delta-actual).mean())
                    r['full_minus_intra_difference_mae'] = r['full_difference_mae']-r['intra_difference_mae']
                    r['full_pred_abs_delta_mean'] = float(np.abs(pred['full_'+target][keep]).mean())
                    non_tied = np.abs(actual) > EPS
                    r['non_tied_truth_pairs'] = int(non_tied.sum())
                    if non_tied.any():
                        s = np.sign(actual[non_tied])
                        for model, key in (('full', 'full_'+target), ('surrogate', 'full_'+cond), ('intra', 'intra_'+target)):
                            predicted = pred[key][keep][non_tied]
                            score = np.where(np.abs(predicted) <= EPS, .5,
                                             (np.sign(predicted) == s).astype(float))
                            r[model+'_direction_accuracy'] = float(score.mean())
                        r['full_pred_tied_pairs'] = int((np.abs(pred['full_'+target][keep][non_tied]) <= EPS).sum())
                        r['full_minus_surrogate_accuracy'] = r['full_direction_accuracy']-r['surrogate_direction_accuracy']
                        r['full_minus_intra_accuracy'] = r['full_direction_accuracy']-r['intra_direction_accuracy']
                rows.append(r)
    pair_table = None
    if mode == 'uniform_rep1':
        # Union includes every pair eligible under any declared base rule.
        union = ((np.abs(rd['rmsf']) <= .25+EPS)|(np.abs(rd['hchi']) <= .25+EPS)|
                 (np.abs(td['rmsf']) <= .20+EPS)|(np.abs(td['hchi']) <= .20+EPS))
        if include_sasa:
            union = (np.abs(sd) <= .25+EPS)&((np.abs(rd['rmsf']) <= .25+EPS)|(np.abs(td['rmsf']) <= .10+EPS))
        pair_table = pd.DataFrame({
            'sample_id': meta['sample_id'], 'dataset': meta['dataset'],
            'amino_acid': meta['amino_acid'],
            'residue_index_i': g.residue_index.to_numpy(int)[i[union]],
            'residue_index_j': g.residue_index.to_numpy(int)[j[union]],
            **{k+'_truth_delta_j_minus_i': v[union] for k, v in td.items()},
            **{k+'_truth_percentile_delta_j_minus_i': v[union] for k, v in rd.items()},
            **{k+'_pred_delta_j_minus_i': v[union] for k, v in pred.items()},
        })
        if include_sasa:
            pair_table['sasa_percentile_delta_j_minus_i'] = sd[union]
    return rows, pair_table


def matched_analysis(positions, out, sasa=False):
    out.mkdir(parents=True, exist_ok=True)
    frames = []
    with gzip.open(out/'uniform_rep1_matched_pairs.csv.gz', 'wt', encoding='utf-8', newline='') as handle:
        first = True
        for mode in MODES:
            for (sid, aa), g in positions.groupby(['sample_id', 'amino_acid'], sort=True):
                if len(g) < 4:
                    continue
                g = g.sort_values('residue_index')
                rows, pairs = pair_stats(g, mode, include_sasa=sasa)
                frames.extend(rows)
                if pairs is not None and len(pairs):
                    pairs.to_csv(handle, index=False, header=first, float_format='%.10g')
                    first = False
    per_type = pd.DataFrame(frames)
    csv(per_type, out/'matched_per_record_AA.csv.gz')
    config = ['label_mode', 'direction', 'match_rule', 'target_subset',
              'sample_id', 'dataset', 'fold', 'family_id']
    records = []
    for key, g in per_type.groupby(config, sort=True):
        r = dict(zip(config, key))
        r['eligible_types'] = len(g)
        r['matched_types'] = int((g.matched_pairs > 0).sum())
        r['ordering_types'] = int(g.full_direction_accuracy.notna().sum())
        for col in COUNT_COLS:
            r[col] = int(g[col].sum())
        for col in METRICS:
            r[col] = mean(g[col])
        r['below_chance_or_undefined'] = not (r['full_direction_accuracy'] > .5)
        r['dedicated_output_not_better_than_surrogate'] = not (r['full_minus_surrogate_accuracy'] > 0)
        records.append(r)
    per_record = pd.DataFrame(records)
    csv(per_record, out/'matched_per_record.csv')
    audit_cols = config+['eligible_types', 'matched_types', 'ordering_types']+list(COUNT_COLS)
    csv(per_record[per_record.matched_types == 0][audit_cols], out/'no_matched_pair_records.csv')
    failure = per_record[(per_record.full_direction_accuracy <= .5)|
                         (per_record.full_minus_surrogate_accuracy <= 0)|
                         per_record.full_direction_accuracy.isna()]
    csv(failure, out/'weak_or_negative_records_all_rules.csv')
    summary = summarize(per_record, ['label_mode', 'direction', 'match_rule', 'target_subset'], METRICS)
    csv(summary, out/'matched_summary_family_cluster_bootstrap.csv')
    primary = summary[(summary.label_mode == 'uniform_rep1')&
                      summary.match_rule.isin(['percentile_0p25', 'physical_0p10'])&
                      (summary.target_subset == 'all_close')&(summary.dataset_scope == 'all')]
    csv(primary, out/'primary_matched_summary.csv')
    return per_record, summary


def bootstrap_matrix(frame, columns, key):
    """Family-cluster bootstrap; sums/counts preserve record-macro estimand."""
    cols = list(columns)
    a = frame[cols].to_numpy(float)
    mask = np.isfinite(a)
    groups, inv = np.unique(frame.family_id.to_numpy(str), return_inverse=True)
    nf = len(groups)
    sums = np.zeros((nf, len(cols)), float)
    counts = np.zeros_like(sums)
    np.add.at(sums, inv, np.where(mask, a, 0))
    np.add.at(counts, inv, mask.astype(float))
    center = np.divide(sums.sum(0), counts.sum(0), out=np.full(len(cols), np.nan), where=counts.sum(0)>0)
    family_means = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=counts>0)
    family_equal = np.array([mean(family_means[:, k]) for k in range(len(cols))])
    family_sizes = frame.groupby('family_id').size().sort_values(ascending=False, kind='mergesort')
    largest = str(family_sizes.index[0])
    excluded = frame[frame.family_id != largest]
    without = np.array([mean(excluded[c]) for c in cols])
    intervals = np.full((2, len(cols)), np.nan)
    if nf >= 2:
        rng = np.random.default_rng(stable_seed(key))
        draws = np.full((B, len(cols)), np.nan)
        prob = np.full(nf, 1/nf)
        for start in range(0, B, 250):
            n = min(250, B-start)
            w = rng.multinomial(nf, prob, size=n).astype(float)
            num, den = w@sums, w@counts
            draws[start:start+n] = np.divide(num, den, out=np.full_like(num, np.nan), where=den>0)
        for k in range(len(cols)):
            finite = draws[np.isfinite(draws[:, k]), k]
            if len(finite):
                intervals[:, k] = np.quantile(finite, [.025, .975])
    return center, intervals, family_equal, without, largest, int(family_sizes.iloc[0]), nf


def summarize(records, configs, metrics):
    rows = []
    for key, g in records.groupby(configs, sort=True):
        if not isinstance(key, tuple):
            key = (key,)
        meta = dict(zip(configs, key))
        for scope in ('all',)+SOURCES:
            p = g if scope == 'all' else g[g.dataset == scope]
            if not len(p):
                continue
            mu, ci, fm, wo, largest, largest_n, nf = bootstrap_matrix(p, metrics, key+(scope,))
            for k, metric in enumerate(metrics):
                finite = p[np.isfinite(p[metric])]
                ref = .5 if metric.endswith('direction_accuracy') else 0.
                rows.append(meta | dict(dataset_scope=scope, metric=metric,
                    record_macro_mean=mu[k], family_cluster_ci95_low=ci[0, k],
                    family_cluster_ci95_high=ci[1, k], family_equal_mean=fm[k],
                    leave_largest_family_out_mean=wo[k],
                    eligible_records=len(p), defined_records=len(finite),
                    available_family_clusters=nf, defined_family_clusters=finite.family_id.nunique(),
                    largest_family_id=largest, largest_family_records=largest_n,
                    records_above_reference=int((finite[metric] > ref).sum()),
                    records_below_reference=int((finite[metric] < ref).sum()),
                    records_equal_reference=int((finite[metric] == ref).sum()), reference_value=ref,
                    matched_pairs=int(p.matched_pairs.sum()) if 'matched_pairs' in p else np.nan,
                    matched_positions=int(p.matched_positions.sum()) if 'matched_positions' in p else np.nan,
                    matched_record_AA_groups=int(p.matched_types.sum()) if 'matched_types' in p else np.nan,
                    bootstrap_resamples=B))
    return pd.DataFrame(rows)


def tails(x):
    state = np.zeros(len(x), int)
    if len(x) < 4:
        return state
    q25, q75 = np.quantile(x, [.25, .75], method='linear')
    if q25 == q75:
        return state
    state[x <= q25] = -1
    state[x >= q75] = 1
    return state


def joint_classes(frame, truth_prefix='target_', pred_prefix='pred_'):
    actual = np.full(len(frame), 'unclassified', dtype='<U24')
    predicted = actual.copy()
    for _, ix in frame.groupby('amino_acid', sort=True).indices.items():
        g = frame.iloc[ix]
        for prefix, result in ((truth_prefix, actual), (pred_prefix, predicted)):
            r, h = tails(g[prefix+'rmsf'].to_numpy(float)), tails(g[prefix+'hchi'].to_numpy(float))
            for cr, ch, label in ((-1, -1, CLASSES[0]), (-1, 1, CLASSES[1]),
                                   (1, -1, CLASSES[2]), (1, 1, CLASSES[3])):
                result[ix[(r == cr)&(h == ch)]] = label
    return actual, predicted


def class_per_record(frame, cohort, verify_existing=False):
    rows, confusion = [], []
    mapping = families()
    mismatches = 0
    for sid, g in frame.groupby('sample_id', sort=True):
        g = g.reset_index(drop=True)
        truth, pred = joint_classes(g)
        if verify_existing:
            mismatches += int((truth != g.reference_class).sum())+int((pred != g.predicted_class).sum())
        meta = dict(cohort=cohort, sample_id=sid, dataset=g.dataset.iloc[0],
                    fold=int(g.fold.iloc[0]), family_id=mapping[sid])
        for label in CLASSES:
            istrue, ispred = truth == label, pred == label
            nt, npred = int(istrue.sum()), int(ispred.sum())
            hit = int((istrue&ispred).sum())
            rows.append(meta | dict(reference_class=label, count=nt, predicted_class_count=npred,
                intersection_count=hit, reference_fraction=nt/len(g), predicted_fraction=npred/len(g),
                recall=hit/nt if nt and label != 'unclassified' else np.nan,
                precision=hit/npred if npred and label != 'unclassified' else np.nan,
                common_positions=len(g)))
            for predlabel in CLASSES:
                confusion.append(meta | dict(reference_class=label, predicted_class=predlabel,
                                             count=int((istrue&(pred == predlabel)).sum())))
    return pd.DataFrame(rows), pd.DataFrame(confusion), mismatches


def review_classes(positions):
    original = pd.read_csv(ENERGY/'oof_analysis/paired_common_residues.csv.gz', low_memory=False,
                           dtype={'chain': str, 'author_residue_id': str})
    assert len(original) == 824988 and original.sample_id.nunique() == 3063
    full, confusion, mismatches = class_per_record(original, 'original_OOF_3063', verify_existing=True)
    assert mismatches == 0, f'Original quartile class mismatch: {mismatches}'
    previous = pd.read_csv(ENERGY/'oof_analysis/per_record_reference_classes.csv', low_memory=False)
    check = full.merge(previous, on=['sample_id', 'reference_class'], validate='one_to_one', suffixes=('_new', '_old'))
    for c in ('count', 'predicted_class_count', 'intersection_count'):
        np.testing.assert_array_equal(check[c+'_new'], check[c+'_old'])
    np.testing.assert_allclose(check.recall, check.rank_class_recall, atol=1e-12, rtol=0, equal_nan=True)
    csv(full, HERE/'four_combinations_full3063_per_record.csv')
    csv(confusion, HERE/'four_combinations_full3063_confusion.csv.gz')
    summaries = [summarize(full, ['cohort', 'reference_class'], ('reference_fraction', 'predicted_fraction', 'recall', 'precision'))]
    for mode in MODES:
        z = positions.rename(columns={mode+'_rmsf': 'target_rmsf', mode+'_hchi': 'target_hchi',
                                      'pred_full_rmsf': 'pred_rmsf', 'pred_full_hchi': 'pred_hchi'})
        classes, conf, _ = class_per_record(z, 'cohort83_'+mode)
        csv(classes, HERE/('four_combinations_83_'+mode+'_per_record.csv'))
        csv(conf, HERE/('four_combinations_83_'+mode+'_confusion.csv'))
        summaries.append(summarize(classes, ['cohort', 'reference_class'], ('reference_fraction', 'predicted_fraction', 'recall', 'precision')))
    summary = pd.concat(summaries, ignore_index=True)
    csv(summary, HERE/'four_combinations_summary_family_cluster_bootstrap.csv')
    return dict(full_cohort_records=3063, full_cohort_positions=824988,
                class_label_mismatches=mismatches, previous_per_record_counts_reproduced=True,
                previous_per_record_recalls_reproduced=True), summary


def sasa_extension(positions, path):
    path = Path(path).resolve()
    features = pd.read_parquet(path)
    choices = ['sidechain_sasa_complex_A2', 'sasa_sidechain_complex_A2']
    col = next((c for c in choices if c in features), None)
    assert col, f'Complex sidechain SASA column missing; columns={features.columns.tolist()}'
    assert not features.duplicated(['sample_id', 'residue_index']).any()
    feature = features[['sample_id', 'residue_index', 'amino_acid', col]].copy()
    feature = feature.rename(columns={col: 'sasa_sidechain_complex_A2', 'amino_acid': 'feature_AA'})
    merged = positions.merge(feature, on=['sample_id', 'residue_index'], how='left', validate='one_to_one')
    assert (merged.loc[merged.feature_AA.notna(), 'amino_acid'] == merged.loc[merged.feature_AA.notna(), 'feature_AA']).all()
    out = HERE/'sasa_extension'
    out.mkdir(exist_ok=True)
    plan = dict(protocol_sha256=sha(HERE/'protocol.json'), input_path=str(path), input_sha256=sha(path),
                sasa_column=col, rule=RULES['sasa_extension'], created_utc=datetime.now(timezone.utc).isoformat())
    write_json(out/'sasa_input_manifest.json', plan)
    audit = merged.assign(sasa_missing=~np.isfinite(merged.sasa_sidechain_complex_A2)).groupby(['sample_id', 'dataset']).agg(
        positions=('residue_index', 'size'), sasa_missing=('sasa_missing', 'sum')).reset_index()
    csv(audit, out/'sasa_missingness.csv')
    finite = merged[np.isfinite(merged.sasa_sidechain_complex_A2)].copy()
    finite['sameAA_n_sasa_finite'] = finite.groupby(['sample_id', 'amino_acid']).residue_index.transform('size')
    csv(finite, out/'common_positions_with_sasa.csv.gz')
    per_record, summary = matched_analysis(finite, out, sasa=True)
    write_json(out/'completion_audit.json', dict(status='complete', input_manifest=plan,
        positions_before=len(positions), finite_sasa_positions=len(finite),
        sameAA4_positions=int((finite.sameAA_n_sasa_finite >= 4).sum()),
        finite_records=finite.sample_id.nunique(), output_per_record_rows=len(per_record),
        output_summary_rows=len(summary), source_files_modified=False, truth_conditioned_explanation_only=True))


def report(positions, matched, combinations):
    primary = matched[(matched.label_mode == 'uniform_rep1')&(matched.dataset_scope == 'all')&
                      (matched.target_subset == 'all_close')&
                      matched.match_rule.isin(['physical_0p10', 'percentile_0p25'])]
    lines = ['# 同类型匹配与主链/侧链互补：审稿式结果说明', '',
             '本分析锁定原先83条共同窗口记录，使用不变的OOF预测，在同记录同氨基酸类型内穷尽全部符合规则的残基对。主要真实指标来自0—100 ns、每200 ps保存一次、共501帧的统一窗口；可用独立重复均值与原标签作为敏感性分析。规则在结果计算前写入protocol.json；没有依据模型预测或成功程度筛选位置。', '',
             '这是一项以真实MD读数分层的解释分析。它不能被称为用可部署输入构造的新预测器，不能从相关或匹配对直接得出结构约束的因果效应，也不能单独证明结合导致的动态变化。', '',
             '83条记录共有33,493个共同位置，其中同类型至少4个位点的组覆盖33,142个位置。可用家族映射仅包含42个簇，最大簇含38条记录；必须以整簇重抽评价不确定性。家族簇置信区间不能消除原始OOF训练/测试之间可能存在的家族相关。', '',
             '| 匹配方向 | 条件 | 指标 | 记录宏平均 [家族簇95%区间] | 去最大簇 |',
             '|---|---|---|---:|---:|']
    keepmetrics = ['target_true_abs_delta_mean', 'full_direction_accuracy',
                   'surrogate_direction_accuracy', 'full_minus_surrogate_accuracy',
                   'intra_direction_accuracy', 'full_minus_intra_accuracy']
    for r in primary[primary.metric.isin(keepmetrics)].itertuples():
        lines.append(f'| {r.direction} | {r.match_rule} | {r.metric} | {r.record_macro_mean:.4f} [{r.family_cluster_ci95_low:.4f}, {r.family_cluster_ci95_high:.4f}] | {r.leave_largest_family_out_mean:.4f} |')
    lines += ['', '所有方向一致率都先在同类型残基对内求均值，再对类型和记录依次等权。真实差为零的对仅从方向分母剔除，预测平局计0.5。代用输出指用Full的预测RMSF排列真实Hχ，或反向用预测Hχ排列真实RMSF；专属输出优于代用输出，才支持“仅靠另一读数不够”的增量判断。不能把匹配对数量当作独立样本数。', '',
              '四组合分类沿用每记录每类型的25%与75%分位，低于或等于25%分位为低，高于或等于75%分位为高；分位重合、类型少于4位点或处于中段均保留为未分类。全3063记录原分类与原召回均逐项重现。此类高低为相对秩分类，不能解读为绝对物理状态；异向组合的召回必须独立报告，不能用同向组合的较高召回概括全部互补能力。', '',
              '主要表格：primary_matched_summary.csv；全部来源、原标签/重复均值敏感性、失败记录及各规则结果见matched_summary_family_cluster_bootstrap.csv和weak_or_negative_records_all_rules.csv。四组合汇总见four_combinations_summary_family_cluster_bootstrap.csv。matched_per_record_AA.csv.gz、matched_per_record.csv及uniform_rep1_matched_pairs.csv.gz提供逐层源表。所有原始数据和Word均未改动。']
    (HERE/'reviewer_interpretation_zh.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--freeze-only', action='store_true')
    ap.add_argument('--sasa', type=Path)
    args = ap.parse_args()
    raw = freeze()
    print('Frozen protocol and source hashes verified.', flush=True)
    if args.freeze_only:
        return
    positions = assemble(raw)
    if args.sasa:
        sasa_extension(positions, args.sasa)
        print('SASA extension completed.', flush=True)
        return
    csv(positions, HERE/'common_positions_with_uniform_truth_and_predictions.csv.gz')
    records = positions[['sample_id', 'dataset', 'fold', 'family_id']].drop_duplicates()
    records['common_positions'] = records.sample_id.map(positions.groupby('sample_id').size())
    records['sameAA4_positions'] = records.sample_id.map(positions[positions.eligible_sameAA].groupby('sample_id').size())
    csv(records, HERE/'cohort83_record_family_inventory.csv')
    print('83-record uniform labels and identities verified; running exhaustive matches.', flush=True)
    per_record, summary = matched_analysis(positions, HERE)
    print('Matched-pair analysis complete; reproducing full-cohort four-combination classes.', flush=True)
    class_audit, class_summary = review_classes(positions)
    report(positions, summary, class_summary)
    audit = dict(status='complete', completed_utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=sha(HERE/'protocol.json'), source_manifest_sha256=sha(HERE/'source_manifest.json'),
        executed_code_sha256=sha(__file__), python=sys.version, numpy=np.__version__, pandas=pd.__version__,
        records=positions.sample_id.nunique(), common_positions=len(positions),
        sameAA4_positions=int(positions.eligible_sameAA.sum()),
        family_clusters=records.family_id.nunique(), largest_family_records=int(records.groupby('family_id').size().max()),
        matching_per_record_rows=len(per_record), matching_summary_rows=len(summary),
        classification=class_audit, original_sources_modified=False, word_modified=False,
        new_training=False, new_inference=False, new_MD=False,
        outputs={p.name:dict(bytes=p.stat().st_size, sha256=sha(p)) for p in HERE.iterdir()
                 if p.is_file() and p.name not in ('completion_audit.json', '.branch_owner')})
    write_json(HERE/'completion_audit.json', audit)
    print(json.dumps(clean(audit), ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
