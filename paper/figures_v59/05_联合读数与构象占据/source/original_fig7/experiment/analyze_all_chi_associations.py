"""Frozen joint-Hchi readout versus actual per-angle MD observations.

All 66 residues, 151 defined angles, and 113 undefined slots are retained.
No model training, new MD, angle prediction, confidence interval, or threshold
selection occurs here. The four-origin MD summaries are verified from angles.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
PACKAGED_SOURCE = OUT.parent / 'source'
ORIGINAL_SOURCE = PACKAGED_SOURCE
FILES = {name: (PACKAGED_SOURCE / name if (PACKAGED_SOURCE / name).is_file() else ORIGINAL_SOURCE / name) for name in [
    'all66_per_angle_metrics.csv', 'all66_residue_summary.csv',
    'ubiquitin_all66_chi_501frames.npz', 'README_数据定义与复算.md']}
PREDICTORS = {'predicted_joint_Hchi': 'predicted_Hchi_nats',
              'MD_reference_joint_Hchi': 'frozen_MD_joint_Hchi_nats'}
METRICS = {'fine_entropy': 'fine_entropy_four_origins_mean_nats',
           'low_Fstar_width': 'width_four_origins_mean_degrees'}

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def rho(x, y):
    x, y = pd.Series(np.asarray(x, float)), pd.Series(np.asarray(y, float))
    assert np.isfinite(x).all() and np.isfinite(y).all()
    if len(x) < 2 or x.nunique() <= 1 or y.nunique() <= 1:
        return np.nan
    return float(x.rank(method='average').corr(y.rank(method='average')))

def manual_rank(values):
    # Independent average-rank implementation used only for arithmetic QA.
    order = sorted(range(len(values)), key=lambda i: values[i])
    rank = np.empty(len(values), float)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and values[order[stop]] == values[order[start]]:
            stop += 1
        rank[order[start:stop]] = (start + 1 + stop) / 2
        start = stop
    return rank

def emit(frame, name):
    target = OUT / name
    temporary = OUT / (name + '.tmp')
    frame.to_csv(temporary, index=False, float_format='%.17g')
    temporary.replace(target)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    checks = []
    def check(name, passed):
        checks.append({'name': name, 'passed': bool(passed)})
        assert passed, name

    hashes = {str(p): sha(p) for p in FILES.values()}
    provenance = {name: {'original_path': str(ORIGINAL_SOURCE / name),
                         'original_sha256': sha(p), 'read_path': str(p)}
                  for name, p in FILES.items()}
    for name, item in provenance.items():
        original = ORIGINAL_SOURCE / name
        if original.is_file():
            check('packaged source equals original ' + name, sha(original) == item['original_sha256'])
    slots = pd.read_csv(FILES['all66_per_angle_metrics.csv'])
    residues = pd.read_csv(FILES['all66_residue_summary.csv'])
    z = dict(np.load(FILES['ubiquitin_all66_chi_501frames.npz'], allow_pickle=False))
    check('66 distinct valid residue identities', len(residues) == 66 and residues.residue_index.nunique() == 66)
    check('264 distinct residue-angle slots', len(slots) == 264 and not slots.duplicated(['residue_index', 'chi_index']).any())
    check('151 defined plus 113 undefined slots', int(slots.chi_defined.sum()) == 151 and int((~slots.chi_defined).sum()) == 113)
    check('all 151 NPZ defined angles retained', np.array_equal(slots.chi_defined.to_numpy().reshape(66, 4), z['chi_defined']))
    check('501-frame 0-100ns common observations', z['chi_degrees'].shape == (501, 66, 4)
          and np.array_equal(z['source_frame_indices'], np.arange(0, 1001, 2)))
    for name in ['predicted_Hchi_nats', 'frozen_MD_joint_Hchi_nats', 'author_residue_id', 'chi_count']:
        check('residue NPZ field ' + name, np.allclose(residues[name].to_numpy(float), z[name].astype(float), atol=1e-13, rtol=0))
    by_index = {int(r): i for i, r in enumerate(z['residue_index'])}
    max_entropy_error = 0.0
    max_width_error = 0.0
    for _, row in slots.iterrows():
        i, q = by_index[int(row.residue_index)], int(row.chi_index)
        defined = bool(row.chi_defined)
        check('slot identity ' + row.residue + ' ' + row.chi,
              str(row.amino_acid) == str(z['amino_acid'][i]) and int(row.author_residue_id) == int(z['author_residue_id'][i]))
        check('slot definition ' + row.residue + ' ' + row.chi, defined == bool(z['chi_defined'][i, q]))
        if not defined:
            check('undefined metric missing ' + row.residue + ' ' + row.chi,
                  all(pd.isna(row[c]) for c in METRICS.values())
                  and np.isnan(z['chi_degrees'][:, i, q]).all())
            continue
        values = z['chi_degrees'][:, i, q].astype(float)
        check('all 501 angle observations finite ' + row.residue + ' ' + row.chi, np.isfinite(values).all())
        entropies, widths = [], []
        for oi, origin in enumerate([0., 1.25, 2.5, 3.75]):
            bins = np.floor(np.remainder(values + 180 - origin, 360) / 5).astype(int)
            counts = np.bincount(bins, minlength=72)
            p = counts / len(values)
            used = p[p > 0]
            entropy = -math.fsum(float(v * math.log(v)) for v in used)
            width = 5 * np.count_nonzero(used >= used.max() / math.e)
            entropies.append(entropy)
            widths.append(width)
            max_entropy_error = max(max_entropy_error, abs(entropy - z['fine_entropy_per_origin_nats'][i, q, oi]))
            max_width_error = max(max_width_error, abs(width - z['width_per_origin_degrees'][i, q, oi]))
        check('four-origin fine entropy ' + row.residue + ' ' + row.chi,
              abs(np.mean(entropies) - row[METRICS['fine_entropy']]) < 3e-12)
        check('four-origin low-Fstar width ' + row.residue + ' ' + row.chi,
              abs(np.mean(widths) - row[METRICS['low_Fstar_width']]) < 3e-12)
    check('independently reconstructed all604 entropy/width summaries', max_entropy_error < 3e-12 and max_width_error < 3e-12)

    # These derived columns only exist on chemically defined angles. Original
    # undefined rows remain present, with NaN residuals and no imputed zeros.
    numeric = list(PREDICTORS.values()) + list(METRICS.values())
    for col in numeric:
        slots[col + '_AA_centered'] = np.nan
    slots['AA_count_within_chi_layer'] = pd.Series(pd.NA, index=slots.index, dtype='Int64')
    slots['singleton_AA_centered_to_zero'] = pd.Series(pd.NA, index=slots.index, dtype='boolean')
    defined = slots[slots.chi_defined].copy()
    centered = defined.groupby(['chi', 'amino_acid'])[numeric].transform('mean')
    for col in numeric:
        slots.loc[defined.index, col + '_AA_centered'] = defined[col] - centered[col]
    counts = defined.groupby(['chi', 'amino_acid']).amino_acid.transform('size')
    slots.loc[defined.index, 'AA_count_within_chi_layer'] = counts.to_numpy()
    slots.loc[defined.index, 'singleton_AA_centered_to_zero'] = (counts == 1).to_numpy()
    defined = slots[slots.chi_defined].copy()
    expected = {'chi1': 66, 'chi2': 52, 'chi3': 23, 'chi4': 10}
    check('layer sample sizes 66/52/23/10', defined.groupby('chi').size().to_dict() == expected)
    check('3 Pro retained twice each', defined[defined.is_proline].residue.tolist() == ['P19', 'P19', 'P37', 'P37', 'P38', 'P38'])
    check('Pro reference joint entropy zeros retained', (residues.loc[residues.is_proline, 'frozen_MD_joint_Hchi_nats'] == 0).all())
    correlation_rows, repeated_rows, group_rows = [], [], []
    amino_acids = sorted(residues.amino_acid.unique())
    for chi, g in defined.groupby('chi', sort=True):
        for aa in amino_acids:
            t = g[g.amino_acid == aa]
            record = {'chi': chi, 'chi_index': int(g.chi_index.iloc[0]) + 1, 'amino_acid': aa,
                      'n_defined': len(t), 'singleton': len(t) == 1, 'Pro_type': aa == 'P'}
            for key, col in {**PREDICTORS, **METRICS}.items():
                record[key + '_constant_within_AA'] = bool(t[col].nunique() <= 1) if len(t) else None
                record[key + '_range_within_AA'] = float(t[col].max() - t[col].min()) if len(t) else np.nan
            group_rows.append(record)
        for predictor, x in PREDICTORS.items():
            for metric, y in METRICS.items():
                r = rho(g[x + '_AA_centered'], g[y + '_AA_centered'])
                # Independently calculate centered values and tied average ranks.
                xx, yy = [], []
                for _, row in g.iterrows():
                    same = g[g.amino_acid == row.amino_acid]
                    xx.append(float(row[x]) - math.fsum(same[x]) / len(same))
                    yy.append(float(row[y]) - math.fsum(same[y]) / len(same))
                # Centering sums can differ by one rounding ULP across routines.
                # Verify ranks using the saved residuals, and arithmetic means separately.
                xc, yc = g[x + '_AA_centered'].to_numpy(), g[y + '_AA_centered'].to_numpy()
                check('AA-centering exact convention ' + chi + ' ' + predictor + ' ' + metric,
                      np.allclose(xx, xc, atol=2e-14, rtol=0) and np.allclose(yy, yc, atol=2e-14, rtol=0))
                manual = float(np.corrcoef(manual_rank(xc), manual_rank(yc))[0, 1]) if len(set(xc)) > 1 and len(set(yc)) > 1 else np.nan
                check('independent Spearman ranks ' + chi + ' ' + predictor + ' ' + metric,
                      (np.isnan(r) and np.isnan(manual)) or abs(r - manual) < 2e-14)
                correlation_rows.append({'chi': chi, 'chi_index': int(g.chi_index.iloc[0]) + 1,
                                         'predictor': predictor, 'metric': metric, 'n_defined': len(g),
                                         'AA_types': int(g.amino_acid.nunique()),
                                         'n_singleton_AA_positions': int((g.AA_count_within_chi_layer == 1).sum()),
                                         'n_proline_positions': int(g.is_proline.sum()),
                                         'rho_raw': rho(g[x], g[y]), 'rho_AA_centered': r,
                                         'status_AA_centered': 'computed' if np.isfinite(r) else 'undefined_constant_layer',
                                         'centered_zero_x_count': int((g[x + '_AA_centered'] == 0).sum()),
                                         'centered_zero_y_count': int((g[y + '_AA_centered'] == 0).sum())})
                repeats = g[g.AA_count_within_chi_layer >= 2]
                repeated_rows.append({'chi': chi, 'predictor': predictor, 'metric': metric,
                                      'n_repeated_AA_positions': len(repeats),
                                      'n_singletons_not_used_in_this_sensitivity_only': len(g) - len(repeats),
                                      'rho_AA_centered_repeated_AA_only': rho(repeats[x + '_AA_centered'], repeats[y + '_AA_centered'])})
    correlations = pd.DataFrame(correlation_rows)
    groups = pd.DataFrame(group_rows)
    k_rows = []
    for k, g in residues.groupby('chi_count'):
        delta = g.predicted_Hchi_nats - g.frozen_MD_joint_Hchi_nats
        k_rows.append({'chi_count': int(k), 'n_residues': len(g), 'AA_types': int(g.amino_acid.nunique()),
                       'n_proline_residues': int(g.is_proline.sum()), 'MAE_joint_Hchi_nats': float(np.abs(delta).mean()),
                       'RMSE_joint_Hchi_nats': float(np.sqrt(np.mean(delta ** 2))),
                       'mean_signed_error_nats': float(delta.mean()),
                       'rho_joint_Hchi': rho(g.predicted_Hchi_nats, g.frozen_MD_joint_Hchi_nats)})
    k = pd.DataFrame(k_rows)
    check('angle-count strata partition all66 once', int(k.n_residues.sum()) == 66)
    check('weighted stratum MAE matches overall MAE', abs(np.average(k.MAE_joint_Hchi_nats, weights=k.n_residues)
          - np.mean(np.abs(residues.predicted_Hchi_nats - residues.frozen_MD_joint_Hchi_nats))) < 1e-13)
    figure_rows = []
    for chi in ['chi1', 'chi2', 'chi3', 'chi4']:
        e = correlations[(correlations.chi == chi) & (correlations.predictor == 'predicted_joint_Hchi') & (correlations.metric == 'fine_entropy')].iloc[0]
        w = correlations[(correlations.chi == chi) & (correlations.predictor == 'predicted_joint_Hchi') & (correlations.metric == 'low_Fstar_width')].iloc[0]
        figure_rows.append({'chi': chi, 'chi_index': int(e.chi_index), 'n': int(e.n_defined),
                            'rho_entropy_typecentered': e.rho_AA_centered,
                            'rho_width_typecentered': w.rho_AA_centered,
                            'rho_entropy_raw': e.rho_raw, 'rho_width_raw': w.rho_raw})
    figure = pd.DataFrame(figure_rows)
    emit(slots, 'all_angle_slots_long.csv')
    emit(defined, 'defined_angles_long.csv')
    emit(residues, 'all66_joint_predictions.csv')
    emit(correlations, 'per_chi_correlations.csv')
    emit(pd.DataFrame(repeated_rows), 'singleton_sensitivity_correlations.csv')
    emit(groups, 'per_chi_AA_counts_and_constant_groups.csv')
    emit(k, 'joint_prediction_accuracy_by_angle_count.csv')
    emit(figure, 'figure_statistics.csv')
    for p, h in hashes.items():
        check('source unchanged ' + Path(p).name, sha(p) == h)
    handling = {
        'positions': 66, 'defined_angles': 151, 'undefined_angle_slots': 113,
        'defined_per_chi': expected,
        'Pro_residues': residues.loc[residues.is_proline, 'residue'].tolist(),
        'Pro_joint_reference': 'Original zero reference values retained; Pro angle observations remain valid and are included.',
        'undefined_angles': 'Not present in any per-chi correlation; retained as 113 missing slots in the 264-row source table. No zero imputation.',
        'AA_centering': 'Within each chi layer, subtract the arithmetic mean for the same amino-acid type separately from frozen joint Hchi and the MD single-angle metric, then compute Spearman over all defined positions in the layer.',
        'singletons': 'Retained in primary analysis. Both centered values are zero for a one-position amino-acid type. A separate repeated-AA-only sensitivity is provided; no selection by observed result.',
        'constant_type_groups': 'Retained as zero residuals/ties. A whole-layer correlation would be NaN if either centered variable were constant. Average ranks handle all ties.',
        'no_inference_intervals': 'No confidence intervals or p-values: one protein and small chi-layer samples, especially n=10 for chi4.',
    }
    definitions = {
        'model_output': 'One frozen per-residue joint Hchi estimate. There are no model-predicted chi1/2/3/4 angle entropies, widths, or free-energy curves.',
        'primary_MD_metric': 'Per-angle 5-degree-bin Shannon entropy, natural logarithm, averaged over four origins 0,1.25,2.5,3.75 degrees.',
        'secondary_MD_metric': 'Per-angle total length of all 5-degree boxes with Fstar=-ln(p/pmax)<=1, averaged over the same four origins; may be disconnected.',
        'time_window': 'Existing 0-100ns trajectory, 501 frames at 200ps intervals.',
        'joint_reference': 'Original joint 120-degree-state Hchi, mean of five 1000-frame windows. A separate observation definition from the per-angle 501-frame fine-bin metrics.',
        'MD_reference_association': 'A descriptive comparison using observed joint Hchi, not a predictor and not a ceiling for single-angle prediction accuracy.',
        'angle_count_accuracy': 'Ordinary accuracy of the original joint Hchi output, stratified by the residue number of defined chi angles k. It is not the accuracy of predicting the kth angle.',
    }
    (OUT / 'definitions_and_handling.json').write_text(json.dumps({'definitions': definitions, 'handling': handling}, ensure_ascii=False, indent=2), encoding='utf8')
    audit = {'status': 'passed', 'checks_passed': len(checks), 'checks': checks,
             'max_reconstructed_entropy_difference': max_entropy_error,
             'max_reconstructed_width_difference_degrees': max_width_error,
             'source_sha256': hashes, 'original_source_provenance': provenance,
             'all_negative_and_small_associations_retained': True,
             'all_four_chi_layers_retained': True, 'no_new_MD_training_or_inference': True,
             'primary_figure_statistics': figure.to_dict(orient='records'),
             'scientific_boundary': 'Associations of a joint readout with individual-angle MD occupancy; not direct per-angle predictions or MAE. Chi4 is weak for entropy and negative for width in this case.',
             'output_sha256': {p.name: sha(p) for p in sorted(OUT.iterdir()) if p.is_file()
                               and p.suffix in ['.csv', '.json', '.py'] and p.name != 'audit.json'}}
    (OUT / 'audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf8')
    print(figure.to_string(index=False))
    print(k.to_string(index=False))
    print(json.dumps({'status': 'passed', 'checks_passed': len(checks), 'output_directory': str(OUT)}, ensure_ascii=False))

if __name__ == '__main__':
    main()
