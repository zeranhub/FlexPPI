"""Independent read-only reconstruction of low-F* width results.

Does not import the producer implementation. Uses natural logs of normalized
histogram probabilities, rather than the producer's exponential count cutoff.
Only this script and independent_verification.json are written.
"""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import sys
sys.dont_write_bytecode = True
import json
import hashlib
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
PLAN = json.loads((ROOT / 'analysis_plan.json').read_text(encoding='utf8'))
PRED = ['full_hchi', 'intrachain_hchi', 'full_rmsf']
checks = {}
numeric_elements = 0
cache_reports = []

def sha(path):
    value = hashlib.sha256()
    with open(path, 'rb') as stream:
        while True:
            block = stream.read(2**20)
            if not block:
                break
            value.update(block)
    return value.hexdigest()

def check(category, condition):
    if not condition:
        raise AssertionError(category)
    checks[category] = checks.get(category, 0) + 1

def close(category, left, right, atol=1e-12, exact=False):
    global numeric_elements
    left, right = np.asarray(left), np.asarray(right)
    numeric_elements += left.size
    if exact:
        np.testing.assert_array_equal(left, right)
    else:
        np.testing.assert_allclose(left, right, rtol=0, atol=atol, equal_nan=True)
    checks[category] = checks.get(category, 0) + 1

def read(path):
    return pd.read_csv(path, dtype={'sample_id': str, 'chain': str,
                                  'author_residue_id': str},
                       float_precision='round_trip')

def ranks(values):
    unique, inverse, counts = np.unique(values, return_inverse=True, return_counts=True)
    cumulative = np.cumsum(counts)
    average = cumulative - (counts - 1) / 2
    return average[inverse]

def spearman(x, y):
    if len(x) < 4 or np.ptp(x) <= 1e-12 or np.ptp(y) <= 1e-12:
        return np.nan
    rx, ry = ranks(x), ranks(y)
    rx -= rx.mean()
    ry -= ry.mean()
    return (rx @ ry) / np.sqrt((rx @ rx) * (ry @ ry))

def main():
    for path, digest in PLAN['source_hash_files'].items():
        check('frozen_source_hash', sha(path) == digest)
    check('producer_code_hash', sha(ROOT / 'run_low_Fstar_width.py') == PLAN['analysis_code_sha256'])
    frozen = read(PLAN['per_residue_sources']['replica1']['path'])
    idcols = ['sample_id', 'dataset', 'fold', 'residue_index', 'amino_acid',
              'chain', 'author_residue_id', 'chi_count'] + PRED
    rebuilt = {'replica1': [], 'available_replicas': []}
    origins = np.asarray(PLAN['bin_origins_degrees'])
    thresholds = np.asarray(PLAN['thresholds'])
    for number, record in enumerate(PLAN['records'], 1):
        base = frozen.loc[frozen.sample_id == record['sample_id']].copy()
        positions = base.residue_index.to_numpy(int)
        q = base.chi_count.to_numpy(int)
        values = []
        for item in record['caches']:
            check('cache_hash', sha(item['path']) == item['sha256'])
            with np.load(item['path'], allow_pickle=False) as raw:
                for field in ('amino_acid', 'chain', 'author_residue_id', 'chi_count'):
                    close('cache_identity', raw[field][positions], base[field].to_numpy(), exact=True)
                close('cache_time_grid', raw['times_ps'], np.arange(501) * 200, atol=.05)
                check('cache_valid', raw['chi_valid'][positions].all() and raw['ca_valid'][positions].all())
                w = np.full((len(base), 4, 4, 3), np.nan)
                maxima = np.zeros((len(base), 4, 4), np.int16)
                for angle_id in range(4):
                    valid = np.flatnonzero(q > angle_id)
                    if len(valid) == 0:
                        continue
                    angles = raw['chi_degrees'][:, positions[valid], angle_id].astype(float).T
                    check('finite_angles', np.isfinite(angles).all())
                    residue_rows = np.broadcast_to(np.arange(len(valid))[:, None], angles.shape)
                    for origin_id, origin in enumerate(origins):
                        # Explicit probabilities and log ratio; disconnected bins
                        # are counted separately and the periodic seam is wrapped.
                        angle_bins = np.floor(np.remainder(angles - (-180 + origin), 360) / 5).astype(int)
                        counts = np.zeros((len(valid), 72), dtype=np.int64)
                        np.add.at(counts, (residue_rows, angle_bins), 1)
                        check('histogram_501_frames', np.all(counts.sum(axis=1) == 501))
                        probabilities = counts / 501
                        maximum_probability = probabilities.max(axis=1)
                        with np.errstate(divide='ignore'):
                            relative_free_energy = -np.log(probabilities / maximum_probability[:, None])
                        maxima[valid, angle_id, origin_id] = counts.max(axis=1)
                        for threshold_id, threshold in enumerate(thresholds):
                            w[valid, angle_id, origin_id, threshold_id] = 5 * (relative_free_energy <= threshold).sum(axis=1)
                metric = np.nansum(w, axis=(1, 2)) / (4 * q[:, None])
                destination = ROOT / 'per_angle_widths' / f"{record['sample_id']}_rep{item['replica']}.npz"
                with np.load(destination, allow_pickle=False) as saved:
                    close('all_angle_widths', w, saved['width_degrees'], exact=True)
                    close('all_angle_max_counts', maxima, saved['max_count'], exact=True)
                    close('all_residue_metrics', metric, saved['metric_mean_per_chi'], exact=True)
                    close('saved_positions', positions, saved['residue_index'], exact=True)
                    close('saved_chi_counts', q, saved['chi_count'], exact=True)
                    close('saved_origins', origins, saved['bin_origins_degrees'], exact=True)
                    close('saved_thresholds', thresholds, saved['thresholds'], exact=True)
                    check('saved_cache_link', saved['source_cache_sha256'].item() == item['sha256'])
                    check('saved_plan_link', saved['plan_sha256'].item() == sha(ROOT / 'analysis_plan.json'))
                cache_reports.append({'sample_id': record['sample_id'], 'replica': item['replica'],
                                      'positions': len(base), 'passed': True})
                values.append((item['replica'], metric))
        first = next(value for replica, value in values if replica == 1)
        all_replicas = np.mean(np.stack([value for _, value in values]), axis=0)
        for mode, metric in [('replica1', first), ('available_replicas', all_replicas)]:
            part = base[idcols].copy()
            for j, threshold in enumerate(thresholds):
                part['width_Fstar_le_' + str(threshold).replace('.', 'p') + '_degrees'] = metric[:, j]
            part['sameAA_group_positions'] = part.groupby('amino_acid').residue_index.transform('size')
            part['eligible_sameAA'] = part.sameAA_group_positions >= 4
            rebuilt[mode].append(part)
        if number % 20 == 0 or number == 83:
            print('Independent width reconstruction', number, '/ 83', flush=True)

    statistics = {}
    bootstrap_results = {}
    for mode, parts in rebuilt.items():
        frame = pd.concat(parts, ignore_index=True)
        out = read(ROOT / f'{mode}_per_residue.csv.gz')
        check('per_residue_counts', len(frame) == len(out) == 33493)
        check('eligible_counts', int(frame.eligible_sameAA.sum()) == 33142)
        for name in frame.columns:
            close('all_per_residue_columns', frame[name].to_numpy(), out[name].to_numpy(), exact=True)
        correlations = read(ROOT / f'{mode}_per_record_correlations.csv')
        correlations = correlations[(correlations.threshold == 1) & (correlations.scope == 'all')].set_index(['sample_id', 'predictor'])
        group_output = read(ROOT / f'{mode}_high_low_per_AA.csv')
        group_output = group_output[(group_output.threshold == 1) & (group_output.scope == 'all')].set_index(['sample_id', 'predictor', 'amino_acid'])
        highlow_output = read(ROOT / f'{mode}_high_low_per_record.csv')
        highlow_output = highlow_output[(highlow_output.threshold == 1) & (highlow_output.scope == 'all')].set_index(['sample_id', 'predictor'])
        record_values = {'correlation': {}, 'high_low': {}}
        for sid, rawpart in frame.groupby('sample_id', sort=True):
            part = rawpart[rawpart.eligible_sameAA].copy()
            y = part.width_Fstar_le_1p0_degrees
            centered_y = y - part.groupby('amino_acid').width_Fstar_le_1p0_degrees.transform('mean')
            for predictor in PRED:
                centered_x = part[predictor] - part.groupby('amino_acid')[predictor].transform('mean')
                rho = spearman(centered_x.to_numpy(), centered_y.to_numpy())
                close('primary_record_correlations', rho, correlations.loc[(sid, predictor), 'rho_type_centered'])
                record_values['correlation'][(sid, predictor)] = rho
                group_diffs, group_lo, group_hi = [], [], []
                for aa, group in part.groupby('amino_acid', sort=True):
                    order = np.lexsort((group.residue_index.to_numpy(), group[predictor].to_numpy()))
                    ordered = group.iloc[order]
                    half = len(group) // 2
                    is_rankable = np.ptp(group[predictor]) > 1e-12
                    saved = group_output.loc[(sid, predictor, aa)]
                    check('group_rankability', bool(saved['rankable']) == is_rankable)
                    if not is_rankable:
                        continue
                    lo, hi = ordered.iloc[:half], ordered.iloc[-half:]
                    check('prediction_led_low_indices', saved.low_indices == ';'.join(map(str, lo.residue_index)))
                    check('prediction_led_high_indices', saved.high_indices == ';'.join(map(str, hi.residue_index)))
                    low, high = lo.width_Fstar_le_1p0_degrees.mean(), hi.width_Fstar_le_1p0_degrees.mean()
                    close('primary_group_actual_MD_means', [low, high, high-low], [saved.low_MD_mean, saved.high_MD_mean, saved.high_minus_low_MD])
                    group_lo.append(low)
                    group_hi.append(high)
                    group_diffs.append(high-low)
                difference, low, high = np.mean(group_diffs), np.mean(group_lo), np.mean(group_hi)
                saved = highlow_output.loc[(sid, predictor)]
                close('primary_record_high_low', [difference, low, high], [saved.high_minus_low_MD, saved.low_MD_mean, saved.high_MD_mean])
                record_values['high_low'][(sid, predictor)] = difference
        statistics[mode] = {}
        for kind in ('correlation', 'high_low'):
            with np.load(ROOT / f'{mode}_{kind}_primary_bootstrap_draws.npz', allow_pickle=False) as saved:
                indices = saved['record_indices']
                ids = saved['record_ids']
                check('bootstrap_shape', indices.shape == (10000, 83))
                check('bootstrap_bounds', np.all(indices >= 0) and np.all(indices < 83))
                check('bootstrap_seed_resamples', saved['seed'].item() == 20261013 and saved['resamples'].item() == 10000)
                close('bootstrap_record_ids', ids, np.array(sorted(frame.sample_id.unique())), exact=True)
                source_by_id = frame.groupby('sample_id').dataset.first().loc[ids].to_numpy()
                for dataset, expected_count in [('atlas', 30), ('dynarepo', 23), ('itrace', 30)]:
                    check('bootstrap_source_counts', np.all((source_by_id[indices] == dataset).sum(axis=1) == expected_count))
                source_values = np.array([[record_values[kind][(sid, predictor)] for predictor in PRED] for sid in ids])
                draws = np.mean(source_values[indices], axis=1)
                summary = read(ROOT / f'{mode}_{kind}_summary.csv')
                summary = summary[(summary.threshold == 1) & (summary.scope == 'all') & (summary.dataset == 'all')].set_index('predictor')
                paired = read(ROOT / f'{mode}_paired_{kind}_difference_summary.csv')
                paired = paired[(paired.threshold == 1) & (paired.scope == 'all') & (paired.dataset == 'all')].set_index('contrast')
                for j, predictor in enumerate(PRED):
                    close('all_10000_primary_bootstrap_draws', draws[:, j], saved[predictor])
                    ci = np.quantile(draws[:, j], [.025, .975])
                    row = summary.loc[predictor]
                    close('primary_summary_and_CI', [source_values[:, j].mean(), *ci], [row['mean'], row.bootstrap_95_low, row.bootstrap_95_high])
                    check('primary_summary_records', row.n_records == 83)
                    statistics[mode][kind + '_' + predictor] = {'mean': float(source_values[:, j].mean()), 'CI': ci.tolist()}
                for j, predictor in enumerate(PRED[1:], 1):
                    contrast = 'full_hchi_minus_' + predictor
                    contrast_draws = draws[:, 0] - draws[:, j]
                    close('all_10000_paired_bootstrap_draws', contrast_draws, saved[contrast])
                    ci = np.quantile(contrast_draws, [.025, .975])
                    row = paired.loc[contrast]
                    values = source_values[:, 0] - source_values[:, j]
                    close('paired_summary_and_CI', [values.mean(), *ci], [row['mean'], row.bootstrap_95_low, row.bootstrap_95_high])
                    close('paired_positive_negative_counts', [(values > 0).sum(), (values < 0).sum(), (values == 0).sum()], [row.positive_records, row.negative_records, row.zero_records], exact=True)
                    statistics[mode][kind + '_' + contrast] = {'mean': float(values.mean()), 'CI': ci.tolist()}
                bootstrap_results[mode + '_' + kind] = {'draws': 10000, 'records_per_draw': 83, 'source_counts': {'atlas': 30, 'dynarepo': 23, 'itrace': 30}}
    return {'status': 'passed', 'verified_utc': datetime.now(timezone.utc).isoformat(),
            'method': 'Independent raw angle histogram reconstruction via np.add.at; explicit -ln(p/pmax) threshold; independent tied-rank covariance Spearman and np.lexsort prediction-led groups; saved paired record bootstrap index replay',
            'plan_sha256': sha(ROOT / 'analysis_plan.json'),
            'code_sha256': sha(__file__), 'read_only_inputs': True, 'records': 83,
            'trajectories': len(cache_reports), 'check_counts': checks,
            'total_checks': sum(checks.values()), 'numeric_elements_compared': numeric_elements,
            'cache_reports': cache_reports, 'statistics': statistics,
            'bootstrap_replays': bootstrap_results}

if __name__ == '__main__':
    try:
        report = main()
    except Exception as exc:
        report = {'status': 'failed', 'error': repr(exc), 'check_counts': checks,
                  'numeric_elements_compared': numeric_elements}
        (ROOT / 'independent_verification.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf8')
        raise
    (ROOT / 'independent_verification.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf8')
    print(json.dumps({k: v for k, v in report.items() if k not in ('cache_reports', 'statistics')}, indent=2))
