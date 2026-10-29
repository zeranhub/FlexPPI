"""Score all supplied baselines on the same frozen valid residues per record.

Relative scores get rank correlation and independently standardized profile
error, never a fitted Angstrom calibration. Missing predictions are disclosed
before scoring and the common intersection is explicit. Primary full-cohort
results require --require-complete.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def statistics(target, prediction, absolute):
    result = dict(count=len(target), spearman=np.nan, pearson=np.nan, z_mae=np.nan, z_rmse=np.nan,
        mae=np.nan, rmse=np.nan, bias=np.nan, top10_precision=np.nan)
    if absolute:
        result.update(mae=float(np.mean(abs(target-prediction))), rmse=float(np.sqrt(np.mean((target-prediction)**2))), bias=float(np.mean(prediction-target)))
    if len(target) >= 4 and np.ptp(target) > 0 and np.ptp(prediction) > 0:
        result.update(spearman=float(spearmanr(target, prediction).statistic), pearson=float(pearsonr(target, prediction).statistic))
        target_z = (target-target.mean()) / target.std(ddof=0)
        prediction_z = (prediction-prediction.mean()) / prediction.std(ddof=0)
        result.update(z_mae=float(np.mean(abs(target_z-prediction_z))), z_rmse=float(np.sqrt(np.mean((target_z-prediction_z)**2))))
        k = max(1, int(np.ceil(.10 * len(target))))
        target_top = set(np.argsort(-target, kind='stable')[:k])
        prediction_top = set(np.argsort(-prediction, kind='stable')[:k])
        result['top10_precision'] = len(target_top & prediction_top) / k
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', required=True)
    parser.add_argument('--prody')
    parser.add_argument('--rmsf-net')
    parser.add_argument('--planet-full')
    parser.add_argument('--planet-seq')
    parser.add_argument('--flexppi', help='Frozen OOF directory containing sample_id.npz')
    parser.add_argument('--intrachain', help='Frozen OOF directory containing sample_id.npz')
    parser.add_argument('--family-folds', help='Optional record_family_folds.csv for new family OOF results')
    parser.add_argument('--output', required=True)
    parser.add_argument('--require-complete', action='store_true')
    args = parser.parse_args()
    cohort, output = Path(args.cohort), Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    freeze = json.loads((cohort / 'cohort_freeze.json').read_text())
    truth_path = cohort / 'truth_residues.csv.gz'
    assert sha(truth_path) == freeze['files']['truth_residues.csv.gz']
    truth = pd.read_csv(truth_path, keep_default_na=False, float_precision='round_trip', dtype={'chain': str, 'author_residue_id': str, 'amino_acid': str})
    assert truth.sample_id.nunique() == 3063 and int(truth.valid_rmsf.sum()) == 961321
    assert not truth.duplicated(['sample_id', 'chain', 'author_residue_id']).any()
    if args.family_folds:
        folds = pd.read_csv(args.family_folds)
        assert folds.sample_id.is_unique and set(folds.sample_id) == set(truth.sample_id)
        truth['fold'] = truth.sample_id.map(dict(zip(folds.sample_id, folds.fold)))
    methods = {}
    if args.prody:
        methods.update(anm=(Path(args.prody), 'anm', False), gnm=(Path(args.prody), 'gnm', False))
    if args.rmsf_net:
        methods['rmsf_net_pdb'] = (Path(args.rmsf_net), 'score', False)
    for name, directory in [('planet_full', args.planet_full), ('planet_seq', args.planet_seq)]:
        if directory:
            methods[name] = (Path(directory), 'rmsf_angstrom', True)
    for name, directory in [('flexppi', args.flexppi), ('intrachain', args.intrachain)]:
        if directory:
            methods[name] = (Path(directory), 'mean', True)
    assert methods, 'Supply at least one completed prediction directory'
    coverage, metrics, aligned = [], [], []
    for sample, rows in truth.groupby('sample_id', sort=False):
        rows = rows.sort_values('residue_index')
        assert np.array_equal(rows.residue_index.to_numpy(), np.arange(len(rows)))
        target, valid = rows.target_rmsf.to_numpy(float), rows.valid_rmsf.to_numpy(bool)
        predictions = {}
        key = 'record_' + hashlib.sha256(sample.encode()).hexdigest()[:24]
        for name, (directory, variable, absolute) in methods.items():
            prediction = np.full(len(rows), np.nan)
            errors = []
            try:
                if name.startswith('planet_'):
                    for input_id, chain_rows in rows.groupby('input_id', sort=False):
                        path = directory / (input_id + '.npz')
                        with np.load(path) as saved:
                            values = saved[variable]
                            assert values.shape == (len(chain_rows),) and np.isfinite(values).all()
                            assert ''.join(chain_rows.sort_values('input_index').amino_acid) == str(saved['sequence'].item())
                            prediction[chain_rows.residue_index.to_numpy(int)] = values[chain_rows.input_index.to_numpy(int)]
                else:
                    path = directory / ((sample if name in ['flexppi', 'intrachain'] else key) + '.npz')
                    with np.load(path) as saved:
                        values = saved[variable][:, 0] if variable == 'mean' else saved[variable]
                        assert values.shape == (len(rows),)
                        if 'chain' in saved:
                            assert saved['chain'].astype(str).tolist() == rows.chain.tolist()
                            assert saved['author_residue_id'].astype(str).tolist() == rows.author_residue_id.tolist()
                        if variable == 'mean':
                            assert np.array_equal(saved['valid'][:, 0], valid)
                            assert np.array_equal(saved['target'][:, 0], target, equal_nan=True)
                        prediction[:] = values
            except Exception as error:
                errors.append(type(error).__name__ + ': ' + str(error))
            predictions[name] = prediction
            coverage.append(dict(sample_id=sample, method=name, expected_valid=int(valid.sum()),
                covered_valid=int((valid & np.isfinite(prediction)).sum()), complete=bool(np.isfinite(prediction[valid]).all()), error=' | '.join(errors)))
        common = valid.copy()
        for prediction in predictions.values():
            common &= np.isfinite(prediction)
        if common.any():
            for name, prediction in predictions.items():
                metrics.append(dict(sample_id=sample, dataset=rows.dataset.iloc[0], biological_system=rows.biological_system.iloc[0],
                    fold=int(rows.fold.iloc[0]), method=name, expected_valid=int(valid.sum()), common_fraction=float(common.sum()/valid.sum()),
                    absolute_angstrom=methods[name][2], **statistics(target[common], prediction[common], methods[name][2])))
        matched = rows.loc[common, ['sample_id', 'dataset', 'fold', 'residue_index', 'chain', 'author_residue_id', 'target_rmsf']].copy()
        for name, prediction in predictions.items():
            matched[name] = prediction[common]
        aligned.append(matched)
    cover = pd.DataFrame(coverage)
    cover.to_csv(output / 'coverage_audit.csv', index=False)
    coverage_summary = cover.groupby('method').agg(records=('sample_id', 'size'), complete_records=('complete', 'sum'), expected_valid=('expected_valid', 'sum'), covered_valid=('covered_valid', 'sum'))
    coverage_summary.to_csv(output / 'coverage_summary.csv')
    audit = dict(records=3063, expected_rmsf_residues=961321, methods={name:dict(directory=str(value[0]), output_variable=value[1], absolute_angstrom=value[2]) for name,value in methods.items()},
        all_methods_complete=bool(cover.complete.all()), complete_required=args.require_complete,
        mask='Frozen label mask AND frame mask AND target-chain mask, intersected with finite predictions of every supplied method. No error-based filtering.',
        relative_score_policy='Spearman, Pearson and per-record independently standardized profile error; no fitted absolute-unit calibration.',
        aggregation='Each record is one vote; chains are pooled within record. Source and fold summaries retained; family/ITrace source concentration must be disclosed.',
        cohort_freeze_sha256=sha(cohort / 'cohort_freeze.json'), truth_sha256=sha(truth_path), code_sha256=sha(__file__),
        family_fold_sha256=sha(args.family_folds) if args.family_folds else None)
    (output / 'scoring_audit.json').write_text(json.dumps(audit, indent=2))
    if args.require_complete and not cover.complete.all():
        raise SystemExit('Incomplete baseline coverage; inspect coverage_audit.csv. No primary metrics published.')
    per = pd.DataFrame(metrics)
    assert len(per), 'No common valid predictions'
    per.to_csv(output / 'per_record_metrics.csv', index=False)
    numeric = ['mae', 'rmse', 'bias', 'spearman', 'pearson', 'z_mae', 'z_rmse', 'top10_precision']
    source = per.groupby(['method', 'dataset'])[numeric].mean().reset_index()
    overall = per.groupby('method')[numeric].mean().reset_index().assign(dataset='all')
    pd.concat([source, overall]).to_csv(output / 'record_macro_metrics.csv', index=False)
    per.groupby(['method', 'dataset', 'fold'])[numeric].mean().to_csv(output / 'source_fold_macro_metrics.csv')
    pd.concat(aligned, ignore_index=True).to_csv(output / 'common_residue_predictions.csv.gz', index=False, compression='gzip')
    print(coverage_summary.to_string())
    print(overall.to_string(index=False))


if __name__ == '__main__':
    main()
