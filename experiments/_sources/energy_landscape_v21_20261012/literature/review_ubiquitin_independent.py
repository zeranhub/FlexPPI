from pathlib import Path
from collections import Counter
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'case_analysis'
BASE = Path(r'E:\Project\FlexPPI\experiments\_case_validation_sources')
ANGLE = BASE / 'ubiquitin/inputs/A003R_Q41_Q62_angles_all_5001_frames.csv'
PROV = BASE / 'data/A003R_Q41_Q62_provenance_and_QA.json'
SUMMARY = BASE / 'data/ubiquitin_Q41_Q62_validation_summary.csv'
MATRIX = BASE / 'data/ubiquitin_Q41_Q62_all27_joint_frequencies.csv'
WINDOWS = BASE / 'data/A003R_Q41_Q62_label_windows.csv'
LABEL_CODE = Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\src\flexppi\data\labels.py')
SCRIPT = ROOT / 'analyze_ubiquitin_landscape.py'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
prov = json.loads(PROV.read_text(encoding='utf-8'))
cache_path = Path(prov['authoritative_numeric_cache'])
label_path = Path(prov['authoritative_label_path'])
cache_sha_before = sha(cache_path)
cache = np.load(cache_path, allow_pickle=False)
labels = np.load(label_path, allow_pickle=False)
source = pd.read_csv(ANGLE)
summary = pd.read_csv(SUMMARY).set_index('residue')
matrix = pd.read_csv(MATRIX)
windows = pd.read_csv(WINDOWS)
reported_stats = pd.read_csv(OUT / 'ubiquitin_bin_sensitivity.csv')
reported_hist = pd.read_csv(OUT / 'ubiquitin_chi1_histograms_and_Fstar.csv')
reported_changes = pd.read_csv(OUT / 'ubiquitin_saved_frame_changes.csv')
reported_states = pd.read_csv(OUT / 'ubiquitin_chi1_state_occupancies.csv')
plan = json.loads((OUT / 'ubiquitin_plan.json').read_text(encoding='utf-8'))
prior_audit = json.loads((OUT / 'ubiquitin_analysis_audit.json').read_text(encoding='utf-8'))


def shannon(categories):
    # A different implementation from the analyzed program: count tuple states directly.
    counts = Counter(categories)
    n = sum(counts.values())
    return -sum((count / n) * np.log(count / n) for count in counts.values())


def original_intervals(degrees):
    return np.where(degrees < -60, 0, np.where(degrees < 60, 1, 2))


def translated_intervals(degrees, phase):
    if phase == 0:
        return original_intervals(degrees)
    wrapped = ((degrees - phase + 180) % 360) - 180
    return original_intervals(wrapped)


def conventional_chi1(degrees):
    theta = ((degrees + 180) % 360) - 180
    return np.where((theta >= -120) & (theta < 0), 0,
                    np.where((theta >= 0) & (theta < 120), 1, 2))


def entropy_joint(states):
    return shannon(map(tuple, states.tolist()))


checks = {
    'source_angle_finite_and_in_closed_180_range': bool(np.isfinite(source[['chi1_degrees', 'chi2_degrees', 'chi3_degrees']]).all().all()
                                                    and source[['chi1_degrees', 'chi2_degrees', 'chi3_degrees']].abs().max().max() <= 180),
    'cache_hash_matches_provenance': cache_sha_before == prov['cache_sha256'],
    'code_hash_matches_analysis_audit': sha(SCRIPT) == prior_audit['code_sha256'],
    'plan_hash_matches_analysis_audit': sha(OUT / 'ubiquitin_plan.json') == prior_audit['plan_sha256'],
    'all_planned_inputs_hash_match': all(sha(path) == expected for path, expected in plan['input_sha256'].items()),
    'full119_label_masks_match_cache': bool(np.array_equal(cache['valid'], labels['mask'])),
    'full119_valid_targets_match_labels': bool(np.allclose(cache['target'][cache['valid']], labels['values'][cache['valid']], atol=1e-7, rtol=0)),
    'full119_valid_variances_match_labels': bool(np.allclose(cache['label_variance'][cache['valid']], labels['variances'][cache['valid']], atol=1e-7, rtol=0)),
}
toy = np.array([-540., -180., -179.999, -120., -119.999, -60., 0., 60., 119.999, 120., 179.999, 180., 540.])
checks['phase60_direct_conventional_intervals_and_periodicity'] = bool(np.array_equal(translated_intervals(toy, 60), conventional_chi1(toy)))
checks['trans_merges_both_sides_of_180_seam'] = bool(conventional_chi1(np.array([-180., -179., 179., 180.])).tolist() == [2, 2, 2, 2])
checks['phase0_preserves_original_exact_endpoint_rule'] = original_intervals(np.array([-180., 180.])).tolist() == [0, 2]

result_rows = []
hist_diagnostics = []
endpoint_details = []
boundary_change_details = []
window_rows = []
for residue, index in [('Q41', 40), ('Q62', 61)]:
    raw = source[source.residue == residue].sort_values('frame_zero_based')
    checks[residue + '_frame_map_complete'] = bool(np.array_equal(raw.frame_zero_based.to_numpy(), np.arange(5001)))
    d = raw[raw.frame_zero_based < 5000].copy()
    a = d[['chi1_degrees', 'chi2_degrees', 'chi3_degrees']].to_numpy()
    checks[residue + '_source_used_csv_matches'] = pd.read_csv(OUT / (residue + '_angles_used.csv')).equals(d.reset_index(drop=True))
    base = original_intervals(a)
    joint_code = (base[:, 0] + 3 * base[:, 1] + 9 * base[:, 2])
    checks[residue + '_original_joint_codes_match'] = bool(np.array_equal(joint_code, d.joint_code.to_numpy()))
    checks[residue + '_record_residue_and_mask_match'] = bool(cache['valid'][index, :2].all()
                                                           and str(cache['chain'][index]) == 'A'
                                                           and str(cache['author_residue_id'][index]) == residue[1:])
    raw_counts = Counter(joint_code.tolist())
    archived = matrix[matrix.residue == residue].sort_values('joint_code')
    checks[residue + '_all27_joint_matrix_match'] = bool(np.array_equal(archived.n_frames.to_numpy(), [raw_counts.get(code, 0) for code in range(27)])
                                                       and np.allclose(archived.probability.to_numpy(), [raw_counts.get(code, 0) / 5000 for code in range(27)], atol=1e-14, rtol=0))
    endpoint_details.append({'residue': residue, 'retained_frame_range': [0, 4999],
                             'retained_time_range_ns': [float(d.time_ns.min()), float(d.time_ns.max())],
                             'excluded_last_frame_time_ns': float(raw[raw.frame_zero_based == 5000].time_ns.iloc[0]),
                             'exact_minus180_count_per_angle': (a == -180).sum(axis=0).tolist(),
                             'exact_plus180_count_per_angle': (a == 180).sum(axis=0).tolist(),
                             'raw_unwrapped_seam_crossings_per_angle': (np.abs(np.diff(a, axis=0)) > 180).sum(axis=0).tolist()})
    for phase in [0, 30, 60, 90]:
        states = translated_intervals(a, phase)
        joint_w = [entropy_joint(states[start:start+1000]) for start in range(0, 5000, 1000)]
        marginals_w = np.array([[shannon(states[start:start+1000, angle].tolist()) for angle in range(3)]
                               for start in range(0, 5000, 1000)])
        marginals_pooled = [shannon(states[:, angle].tolist()) for angle in range(3)]
        joint_mean = float(np.mean(joint_w))
        joint_pooled = entropy_joint(states)
        saved_changes = int(np.count_nonzero(states[1:, 0] != states[:-1, 0]))
        result_rows.append({'residue': residue, 'phase_degrees': phase,
                            'joint_window_mean': joint_mean, 'joint_pooled': joint_pooled,
                            'joint_pooled_minus_window_mean': joint_pooled - joint_mean,
                            'chi1_window_mean': float(marginals_w[:, 0].mean()),
                            'chi2_window_mean': float(marginals_w[:, 1].mean()),
                            'chi3_window_mean': float(marginals_w[:, 2].mean()),
                            'chi1_pooled': marginals_pooled[0], 'chi2_pooled': marginals_pooled[1],
                            'chi3_pooled': marginals_pooled[2],
                            'conditional_chi2_chi3_given_chi1_window_mean': joint_mean - float(marginals_w[:, 0].mean()),
                            'saved_chi1_frame_changes': saved_changes})
        match = reported_stats[(reported_stats.residue == residue) & (reported_stats.phase == phase)].iloc[0]
        recomputed = [joint_mean, joint_pooled, float(marginals_w[:, 0].mean()), marginals_pooled[0], float(np.std(joint_w, ddof=1))]
        printed = match[['joint_window_mean_h', 'joint_pooled_h', 'chi1_window_mean_h', 'chi1_pooled_h', 'joint_window_sd']].to_numpy(dtype=float)
        checks[f'{residue}_phase{phase}_entropy_csv_reproduced'] = bool(np.allclose(recomputed, printed, atol=1e-12, rtol=0))
        reported_change = reported_changes[(reported_changes.residue == residue) & (reported_changes.phase == phase)].iloc[0]
        checks[f'{residue}_phase{phase}_change_csv_reproduced'] = saved_changes == int(reported_change.saved_frame_changes)
        reported_state = reported_states[(reported_states.residue == residue) & (reported_states.phase == phase)].sort_values('state')
        state_counts = Counter(states[:, 0].tolist())
        checks[f'{residue}_phase{phase}_occupancies_reproduced'] = bool(np.allclose([state_counts.get(code, 0) / 5000 for code in range(3)], reported_state.probability.to_numpy(), atol=1e-14, rtol=0))
        if phase == 0:
            checks[residue + '_frozen_Hchi_reproduced'] = abs(joint_mean - float(cache['target'][index, 1])) < 1e-6
            checks[residue + '_frozen_Hchi_variance_reproduced'] = abs(float(np.var(joint_w, ddof=1)) - float(cache['label_variance'][index, 1])) < 1e-6
            checks[residue + '_frozen_predictions_match_case_summary'] = bool(np.allclose(cache['mean'][index, :2], summary.loc[residue, ['model_rmsf', 'model_hchi']].to_numpy(dtype=float), atol=1e-7, rtol=0))
            for w, entropy_value in enumerate(joint_w, start=1):
                archived_window = windows[(windows.residue == residue) & (windows.window == w)].iloc[0]
                checks[f'{residue}_window{w}_archived_entropy_match'] = abs(entropy_value - float(archived_window.hchi_nats)) < 1e-12
                window_rows.append({'residue': residue, 'window': w, 'joint_Hchi': entropy_value,
                                     'start_frame': (w - 1) * 1000, 'stop_frame_exclusive': w * 1000})
    base_changes = base[1:, 0] != base[:-1, 0]
    centred = conventional_chi1(a[:, 0])
    centred_changes = centred[1:] != centred[:-1]
    boundary_change_details.append({'residue': residue,
                                    'original_only_changes': int(np.count_nonzero(base_changes & ~centred_changes)),
                                    'centred_only_changes': int(np.count_nonzero(~base_changes & centred_changes)),
                                    'both_definitions_change': int(np.count_nonzero(base_changes & centred_changes)),
                                    'fraction_original_changes_within_same_conventional_state': float(np.count_nonzero(base_changes & ~centred_changes) / np.count_nonzero(base_changes))})
    for window in range(6):
        aq = a if window == 0 else a[(window - 1) * 1000:window * 1000]
        for chi in range(1, 4):
            # Independent direct indexing of 5-degree histogram, including NumPy's final closed endpoint.
            bucket = np.minimum(np.floor((aq[:, chi-1] + 180) / 5).astype(int), 71)
            counts = Counter(bucket.tolist())
            count_vector = np.array([counts.get(bin_id, 0) for bin_id in range(72)])
            probabilities = count_vector.astype(float) / len(aq)
            visited = probabilities > 0
            expected_f = np.array([-np.log(probabilities[bin_id] / probabilities.max()) if visited[bin_id] else np.nan
                                   for bin_id in range(72)])
            h = reported_hist[(reported_hist.residue == residue) & (reported_hist.window == window)
                              & (reported_hist.chi == chi)].sort_values('angle_lower')
            ok = (len(h) == 72 and np.array_equal(h['count'].to_numpy(), count_vector)
                  and np.allclose(h.probability, probabilities, atol=1e-14, rtol=0)
                  and np.allclose(h.density_per_degree, probabilities / 5, atol=1e-14, rtol=0)
                  and np.allclose(h.F_star, expected_f, equal_nan=True, atol=1e-12, rtol=0)
                  and h.F_star.isna().to_numpy().tolist() == (~visited).tolist())
            checks[f'{residue}_window{window}_chi{chi}_histogram_and_Fstar_match'] = bool(ok)
            hist_diagnostics.append({'residue': residue, 'window': window, 'chi': chi,
                                      'count_sum': int(count_vector.sum()), 'density_integral': float((h.density_per_degree * 5).sum()),
                                      'visited_bins': int(visited.sum()), 'empty_bins_as_NaN': int((~visited).sum()),
                                      'minimum_Fstar': float(np.nanmin(expected_f))})

checks['cache_unchanged_during_review'] = sha(cache_path) == cache_sha_before
checks = {key: bool(value) for key, value in checks.items()}
numeric = pd.DataFrame(result_rows)
numeric.to_csv(OUT / 'ubiquitin_independent_numeric_audit.csv', index=False)
qa = {
    'status': 'pass' if all(checks.values()) else 'requires_inspection',
    'review_type': 'Read-only independent recomputation from existing saved angle CSV; no model, source or label mutation',
    'checks': checks, 'failed_checks': [key for key, value in checks.items() if not value],
    'input_sha256': {str(p): sha(p) for p in [ANGLE, PROV, SUMMARY, MATRIX, WINDOWS, cache_path, label_path, LABEL_CODE, SCRIPT]},
    'cache_protocol': prov['model_protocol'], 'cache_sha256': cache_sha_before,
    'endpoint_and_frame_details': endpoint_details,
    'boundary_change_details': boundary_change_details,
    'independent_entropies': result_rows, 'frozen_window_entropies': window_rows,
    'histogram_diagnostics': hist_diagnostics,
    'interpretation': [
        'Positive phase means wrap(theta-phase), not wrap(theta+phase). At phase60 it agrees exactly with non-Pro conventional chi1 intervals and merges the +/-180 trans seam.',
        'Original phase0 deliberately preserves clip/frozen endpoint handling. The six retained angle series contain no exact +/-180 values, so the phase0 endpoint exception causes no numeric difference for this case.',
        'Five-window mean entropy reproduces frozen Hchi; pooled entropy is a different statistic and is larger for these records by entropy concavity.',
        'Phase60 joint entropy shifts all three angles uniformly and is a partition-sensitivity result. It is not a complete physical rotamer correction for Gln terminal chi3.',
        'Phase60 chi1 marginal entropy is higher for Q41; the complete joint entropy remains higher for Q62. Conditional diversity of chi2/chi3 given chi1 accounts algebraically for that difference.',
        'Saved-frame changes are changes at 100ps observation spacing, not transition rates. Most original changes stay within the same conventional chi1 state.',
        'Every zero-count histogram bin retains NaN Fstar, every visited bin matches -ln(p/pmax), and each curve has minimum zero. Plotting must preserve gaps rather than joining through missing bins.',
        'CSV name ubiquitin_chi1_histograms_and_Fstar.csv is inherited, but its chi column contains all three angles. A three-column figure must filter chi1/chi2/chi3 separately.',
        'The analysis retains frames0..4999 (0..499.9ns) of the existing5001-frame trajectory and matches the five1000-frame label windows; frame5000 at500ns is excluded.'
    ]
}
(OUT / 'ubiquitin_independent_content_QA.json').write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'status': qa['status'], 'failed_checks': qa['failed_checks'], 'cache_sha256': cache_sha_before,
                  'n_checks': len(checks), 'endpoints': endpoint_details, 'boundary_changes': boundary_change_details}, ensure_ascii=False, indent=2))
print(numeric.to_string(index=False))
