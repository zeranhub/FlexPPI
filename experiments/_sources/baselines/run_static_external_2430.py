"""Re-fit frozen Static/AA controls on new complete-chain external purge2430."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
import joblib
from sklearn.ensemble import HistGradientBoostingRegressor
from scipy.stats import pearsonr, spearmanr
from threadpoolctl import threadpool_limits
from static_feature_definitions import weights


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def statistics(target, prediction):
    result = dict(count=len(target), mae=float(np.mean(abs(target-prediction))), rmse=float(np.sqrt(np.mean((target-prediction)**2))), bias=float(np.mean(prediction-target)), spearman=np.nan, pearson=np.nan)
    if len(target) >= 4 and np.ptp(target) > 0 and np.ptp(prediction) > 0:
        result.update(spearman=float(spearmanr(target, prediction).statistic), pearson=float(pearsonr(target, prediction).statistic))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', required=True)
    parser.add_argument('--family-root', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    reference, family, output = Path(args.reference), Path(args.family_root), Path(args.output)
    if (output / 'run_audit.json').exists():
        raise RuntimeError('Completed fit already exists; use a new directory')
    output.mkdir(parents=True, exist_ok=True)
    allowed = {json.loads(line)['sample_id'] for line in (family / 'splits/fold_0_train.jsonl').read_text().splitlines()}
    freeze = json.loads((family / 'training_freeze.json').read_text())
    assert len(allowed) == 2430 and freeze['largest_training_component'] == 633
    assert sha(family / 'splits/fold_0_train.jsonl') == freeze['files']['splits/fold_0_train.jsonl']
    data = dict(np.load(reference / 'training_features.npz'))
    schema = json.loads((reference / 'training_schema.json').read_text())
    external = dict(np.load(reference / 'external_features.npz'))
    ext_schema = json.loads((reference / 'external_schema.json').read_text())
    assert len(schema['records']) == 3063 and len(ext_schema['records']) == 32
    assert set(allowed) <= {r['sample_id'] for r in schema['records']}
    assert schema['feature_source_sha256'] == sha(Path(__file__).with_name('static_feature_definitions.py'))
    selected = [r['record_index'] for r in schema['records'] if r['sample_id'] in allowed]
    member = np.isin(data['record_index'], selected)
    x, target, valid = data['X'], data['target'], data['valid']
    names = schema['names']
    categorical = [name in schema['categorical'] for name in names]
    configuration = dict(loss='squared_error', learning_rate=.07, max_iter=120, max_leaf_nodes=31, max_depth=None,
        min_samples_leaf=50, l2_regularization=1., max_bins=63, early_stopping=False, random_state=7)
    predictions = {name:np.full(external['target'].shape, np.nan, np.float32) for name in ['static_sequence_complex', 'aa_only']}
    started = time.time()
    audit = dict(status='running', training_records=2430, removed_records=633, external_records=32, configuration=configuration,
        family_freeze_sha256=sha(family / 'training_freeze.json'), training_manifest_sha256=sha(family / 'splits/fold_0_train.jsonl'),
        reference_hashes={name:sha(reference / name) for name in ['training_features.npz','training_schema.json','external_features.npz','external_schema.json']},
        code_sha256=sha(__file__), feature_code_sha256=sha(Path(__file__).with_name('static_feature_definitions.py')),
        no_external_tuning=True, independent_external_family_components=1, confidence_intervals='No independent-family CI or p-value', jobs=[])
    with threadpool_limits(limits=args.threads):
        for m, metric in enumerate(['rmsf', 'chi_entropy']):
            train = member & valid[:, m]
            balance = weights(schema['records'], data['record_index'], train)
            model = HistGradientBoostingRegressor(**configuration, categorical_features=categorical)
            tick = time.time()
            model.fit(x[train], target[train, m], sample_weight=balance)
            joblib.dump(model, output / f'static_{metric}.joblib', compress=3)
            predictions['static_sequence_complex'][:, m] = model.predict(external['X'])
            aa = x[:, names.index('aa_offset_+0')].astype(int)
            centers = np.full(22, np.average(target[train, m], weights=balance))
            ids = np.flatnonzero(train)
            for code in np.unique(aa[train]):
                pick = aa[ids] == code
                centers[code] = np.average(target[ids[pick], m], weights=balance[pick])
            predictions['aa_only'][:, m] = centers[external['X'][:, names.index('aa_offset_+0')].astype(int)]
            np.save(output / f'aa_centers_{metric}.npy', centers)
            audit['jobs'].append(dict(metric=metric, train_residues=int(train.sum()), seconds=time.time()-tick))
            print(metric, int(train.sum()), time.time()-tick, flush=True)
    metrics = []
    for name, prediction in predictions.items():
        assert np.isfinite(prediction[external['valid']]).all()
        np.savez_compressed(output / (name+'_external.npz'), prediction=prediction)
        for record in ext_schema['records']:
            indices = np.arange(record['start'], record['stop'])
            for m, metric in enumerate(['rmsf','chi_entropy']):
                selected = indices[external['valid'][indices,m]]
                metrics.append(dict(sample_id=record['sample_id'], variant=name, metric=metric,
                    **statistics(external['target'][selected,m], prediction[selected,m])))
    per = pd.DataFrame(metrics)
    per.to_csv(output / 'per_record_metrics.csv', index=False)
    per.groupby(['variant','metric'])[['mae','rmse','bias','spearman','pearson']].mean().to_csv(output / 'record_macro_metrics.csv')
    audit.update(status='complete', seconds=time.time()-started, valid_external=external['valid'].sum(0).tolist())
    (output / 'run_audit.json').write_text(json.dumps(audit, indent=2))
    print(per.groupby(['variant','metric'])[['mae','spearman']].mean().to_string(), flush=True)


if __name__ == '__main__':
    main()
