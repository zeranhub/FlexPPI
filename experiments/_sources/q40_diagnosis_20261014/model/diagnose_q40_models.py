"""Read-only diagnosis of frozen fold-0 model outputs at ubiquitin GLN sites.

Writes only compact derived CSV/JSON files beside this script. No training,
inference, manuscript edits, trajectory changes, or external services.
"""
from pathlib import Path
from collections import Counter
import hashlib
import json
import math
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
PIPE = Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
FULL = PIPE / 'neural_exports/evaluation/full/A003R.npz'
INTRA = PIPE / 'neural_exports/evaluation/intrachain/A003R.npz'
STATIC = PIPE / 'static_corrected/static_sequence_complex_oof.npz'
AA = PIPE / 'static_corrected/aa_only_oof.npz'
FEATURES = PIPE / 'static_corrected/training_features.npz'
SCHEMA = PIPE / 'static_corrected/training_schema.json'
PROFILE = Path(r'E:\Project\FlexPPI\experiments\_case_validation_sources\data\ubiquitin_fold0_all_residues.csv')
WINDOWS = Path(r'E:\Project\FlexPPI\experiments\_sources\free_energy_curves_v27_20261014\case_curves\joint_entropy_original_five_windows.csv')
SOURCE_FILES = [FULL, INTRA, STATIC, AA, FEATURES, SCHEMA, PROFILE, WINDOWS,
                PIPE / 'run_corrected_static.py', PIPE / 'static_feature_definitions.py',
                PIPE / 'training_payload/evaluate_corrected_models.py',
                PIPE / 'training_payload/src/flexppi/losses.py',
                PIPE / 'training_payload/src/flexppi/model/heads.py',
                PIPE / 'training_payload/src/flexppi/data/labels.py']

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def save_csv(frame, name):
    frame.to_csv(OUT / name, index=False, float_format='%.17g')

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    hashes = {str(p): sha(p) for p in SOURCE_FILES}
    f = dict(np.load(FULL, allow_pickle=False))
    intra = dict(np.load(INTRA, allow_pickle=False))
    schema = json.loads(SCHEMA.read_text(encoding='utf8'))
    record = next(r for r in schema['records'] if r['sample_id'] == 'A003R')
    assert record['fold'] == 0
    sl = slice(record['start'], record['stop'])
    static = np.load(STATIC, allow_pickle=False)['prediction'][sl, 1].astype(float)
    aa_pred = np.load(AA, allow_pickle=False)['prediction'][sl, 1].astype(float)
    for key in ['author_residue_id', 'chain', 'amino_acid', 'valid', 'target', 'label_variance']:
        np.testing.assert_array_equal(f[key], intra[key])
    with np.load(FEATURES, allow_pickle=False) as features:
        x = features['X']
        y = features['target'][:, 1].astype(float)
        valid = features['valid'][:, 1].astype(bool)
        fold = features['fold']
        record_ids = features['record_index']
    names = schema['names']
    np.testing.assert_array_equal(y[sl], f['target'][:, 1])
    np.testing.assert_array_equal(valid[sl], f['valid'][:, 1])
    aa_codes = x[:, names.index('aa_offset_+0')].astype(int)
    np.testing.assert_array_equal(aa_codes[sl], f['aatype'])
    gln_code = int(f['aatype'][1])
    train = (fold != 0) & valid
    train_q = train & (aa_codes == gln_code)
    residue_mean = float(y[train_q].mean())
    # Reconstruct the original source/system/record-balanced AA baseline.
    counts = np.bincount(record_ids[train], minlength=len(schema['records']))
    active = [r for r in schema['records'] if counts[r['record_index']] > 0]
    system_records = Counter((r['dataset'], r['biological_system']) for r in active)
    source_systems = Counter(source for source, system in system_records)
    record_weight = np.zeros(len(schema['records']))
    for r in active:
        i = r['record_index']
        record_weight[i] = 1 / (len(source_systems) * source_systems[r['dataset']]
                                * system_records[(r['dataset'], r['biological_system'])]
                                * counts[i])
    q_weights = record_weight[record_ids[train_q]]
    balanced_mean = float(np.average(y[train_q], weights=q_weights))
    q_record_mean = float(pd.DataFrame({'record': record_ids[train_q], 'target': y[train_q]})
                          .groupby('record').target.mean().mean())
    indices = np.flatnonzero((f['chain'] == 'A') & f['valid'][:, 1])
    table = pd.DataFrame({
        'residue_index': indices,
        'residue_id': f['author_residue_id'][indices].astype(int),
        'amino_acid': f['amino_acid'][indices],
        'MD_Hchi_nats': f['target'][indices, 1].astype(float),
        'full_Hchi_nats': f['mean'][indices, 1].astype(float),
        'intra_Hchi_nats': intra['mean'][indices, 1].astype(float),
        'static_Hchi_nats': static[indices],
        'AA_balanced_Hchi_nats': aa_pred[indices],
        'full_residual_variance_nats2': f['model_variance'][indices, 1].astype(float),
        'intra_residual_variance_nats2': intra['model_variance'][indices, 1].astype(float),
        'label_window_sample_variance_nats2': f['label_variance'][indices, 1].astype(float),
    })
    for name in ['full', 'intra', 'static', 'AA_balanced']:
        table[name + '_signed_error_nats'] = table[name + '_Hchi_nats'] - table.MD_Hchi_nats
        table[name + '_absolute_error_nats'] = table[name + '_signed_error_nats'].abs()
    for name, col in [('full_residual_SD_nats', 'full_residual_variance_nats2'),
                      ('intra_residual_SD_nats', 'intra_residual_variance_nats2'),
                      ('label_window_SD_nats', 'label_window_sample_variance_nats2')]:
        table[name] = np.sqrt(table[col])
    table['full_minus_intra_nats'] = table.full_Hchi_nats - table.intra_Hchi_nats
    q = table[table.amino_acid == 'Q'].copy()
    assert len(table) == 66 and len(q) == 6
    assert abs(balanced_mean - q.AA_balanced_Hchi_nats.iloc[0]) < 2e-7
    assert q.AA_balanced_Hchi_nats.nunique() == 1
    q['unweighted_training_Q_mean_nats'] = residue_mean
    q['full_minus_unweighted_training_Q_mean_nats'] = q.full_Hchi_nats - residue_mean
    q['full_minus_balanced_training_Q_mean_nats'] = q.full_Hchi_nats - balanced_mean
    window = pd.read_csv(WINDOWS)
    window_summary = window.groupby('residue').Hchi_nats.agg(['mean', 'var', 'std', 'min', 'max']).reset_index()
    window_summary['residue_id'] = window_summary.residue.str.removeprefix('Q').astype(int)
    for _, r in q.iterrows():
        observed = window_summary[window_summary.residue_id == r.residue_id].iloc[0]
        assert abs(observed['mean'] - r.MD_Hchi_nats) < 2e-7
        assert abs(observed['var'] - r.label_window_sample_variance_nats2) < 2e-8
    q = q.merge(window_summary.rename(columns={c: 'five_window_' + c for c in ['mean', 'var', 'std', 'min', 'max']}),
                on='residue_id', how='left').sort_values('residue_id')
    feature_table = pd.DataFrame(x[record['start'] + q.residue_index.to_numpy()], columns=names)
    feature_table.insert(0, 'residue_id', q.residue_id.to_numpy())
    for name in names:
        if 'contacts_' in name:
            feature_table[name + '_raw_count'] = np.rint(np.expm1(feature_table[name])).astype(int)
    q40 = q[q.residue_id == 40].iloc[0]
    pairs = []
    for site in [31, 41, 49, 62, 2]:
        other = q[q.residue_id == site].iloc[0]
        pairs.append({'comparison': 'Q40_minus_Q' + str(site),
                      'MD_difference_nats': q40.MD_Hchi_nats - other.MD_Hchi_nats,
                      'full_difference_nats': q40.full_Hchi_nats - other.full_Hchi_nats,
                      'intra_difference_nats': q40.intra_Hchi_nats - other.intra_Hchi_nats,
                      'static_difference_nats': q40.static_Hchi_nats - other.static_Hchi_nats})
    result = {
        'status': 'verified_frozen_outputs_only',
        'record': record,
        'training_Q_means': {
            'fold': 0, 'test_fold_excluded': True,
            'valid_training_Q_residue_count': int(train_q.sum()),
            'unweighted_residue_mean_nats': residue_mean,
            'record_equal_Q_mean_nats': q_record_mean,
            'source_system_record_balanced_Q_mean_nats': balanced_mean,
            'exported_AA_baseline_float32_nats': float(q40.AA_balanced_Hchi_nats),
            'AA_baseline_reconstruction_tolerance': 2e-7,
            'type_corrected_plot_uses': 'unweighted residue mean; the AA predictive baseline uses balanced weights',
        },
        'Q40': {c: (int(q40[c]) if c in ['residue_id', 'residue_index'] else float(q40[c]))
                for c in q.columns if c not in ['amino_acid', 'residue']},
        'Q40_ranks_among_66': {
            'full_absolute_error_descending': int(table.full_absolute_error_nats.rank(method='min', ascending=False)[table.residue_id == 40].iloc[0]),
            'label_variance_descending': int(table.label_window_sample_variance_nats2.rank(method='min', ascending=False)[table.residue_id == 40].iloc[0]),
        },
        'definitions': {
            'model_variance': 'exp(network log_variance), a learned extra residual variance in the heteroscedastic regression objective. Single fold-0 network output; not between-model or posterior uncertainty, nor a measured physical fluctuation.',
            'label_variance': 'Sample variance, ddof=1, of five continuous 1000-frame Hchi observations; not variance of raw angles and not the standard error of their mean.',
            'training_loss': '0.5 * [(prediction-target)^2 / (model_variance + label_variance) + ln(model_variance + label_variance)] on valid positions, then within-record and source/system sampling aggregation.',
            'Full_vs_Intra': 'Separately trained networks with matched architecture; Full minus Intra is an encoding/training comparison, not a same-weight partner-deletion intervention.',
        },
        'inferences': {
            'confirmed': [
                'Q40 is overpredicted by all four exported baselines; Intra, Static, and AA are closer than Full at this one site.',
                'Full is above both unweighted and balanced training-Q means; simple regression to a Q-only mean is insufficient.',
                'Full raises Q40 by about 0.170 nats relative to Intra, whereas the reference is lower; this is a local Full-model error, not overall model performance.',
                'Q40 has the largest five-window label variance among the six GLN sites. The exported residual scale is not the largest and is not a reliable causal explanation.',
            ],
            'not_identified': [
                'Whether partner contacts, a particular intrachain geometry, or missing side-chain/environmental features caused this error.',
                'Whether the static starting state preferentially represents one observed MD window; proximity of a prediction to one window does not establish that mechanism.',
            ],
        },
        'sources_sha256': hashes,
    }
    save_csv(table, 'all66_frozen_model_comparison.csv')
    save_csv(q, 'all6_Gln_frozen_model_comparison.csv')
    save_csv(table[table.residue_id.between(35, 45)], 'sequence_neighborhood_35_45.csv')
    save_csv(feature_table, 'all6_Gln_static_input_features.csv')
    save_csv(pd.DataFrame(pairs), 'Q40_same_type_pairwise_contrasts.csv')
    save_csv(window_summary, 'all6_Gln_five_window_summary.csv')
    (OUT / 'diagnosis_summary.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf8')
    print(json.dumps({'status': result['status'], 'Q40': result['Q40'], 'training_Q_means': result['training_Q_means'],
                      'Q40_ranks_among_66': result['Q40_ranks_among_66']}, ensure_ascii=False))

if __name__ == '__main__':
    main()
