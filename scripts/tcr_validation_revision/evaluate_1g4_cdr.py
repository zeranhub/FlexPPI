"""Evaluate all IMGT CDRs on frozen 1G4 predictions and matched MD references.

This script does not train, refit, select snapshots, or change old artifacts.
Whole-run bootstrap preserves each prediction/reference pair, while apo and
holo ensembles are sampled independently. Residues are never bootstrapped.
"""
from pathlib import Path
import hashlib
import json
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'output/tcr_prediction_reference_revision_20260915'
OUT = ROOT / 'output/tcr_cdr_reassessment_20260917'
ANNOTATION = ROOT / 'artifacts/tcr_cdr_reassessment_20260917/annotation_agent'
KEY = ['chain', 'author_id']
SEED = 20260917
N_BOOT = 10000


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def score(reference, prediction):
    reference = np.asarray(reference, float)
    prediction = np.asarray(prediction, float)
    assert np.isfinite(reference).all() and np.isfinite(prediction).all()
    error = prediction - reference
    # Two points necessarily give +/-1; do not report that as performance.
    rho = float('nan')
    if len(reference) >= 3 and np.ptp(reference) > 0 and np.ptp(prediction) > 0:
        rho = float(spearmanr(reference, prediction).statistic)
    return dict(n=len(reference), spearman=rho,
                mae=float(np.abs(error).mean()),
                rmse=float(np.sqrt(np.mean(error ** 2))),
                bias=float(error.mean()),
                reference_mean=float(reference.mean()),
                prediction_mean=float(prediction.mean()),
                small_region_warning=len(reference) <= 5)


def subsets(mapping):
    region = mapping.region.to_numpy()
    chain = mapping.chain.to_numpy()
    result = [('all_TCR', np.ones(len(mapping), bool)),
              ('all_CDR', np.isin(region, ['CDR1', 'CDR2', 'CDR3'])),
              ('CDR1', region == 'CDR1'), ('CDR2', region == 'CDR2'),
              ('CDR3', region == 'CDR3'),
              ('non_CDR', ~np.isin(region, ['CDR1', 'CDR2', 'CDR3'])),
              ('V_framework', region == 'V_framework'),
              ('non_variable_remainder', region == 'non_variable_remainder')]
    for ch in ['alpha', 'beta']:
        result.append((ch + '_all_CDR', (chain == ch) & np.isin(region, ['CDR1', 'CDR2', 'CDR3'])))
        for r in ['CDR1', 'CDR2', 'CDR3']:
            result.append((ch + '_' + r, (chain == ch) & (region == r)))
    return result


def bootstrap_intervals(md, pred):
    gain = np.abs(md).mean(axis=1) - np.abs(pred - md).mean(axis=1)
    result = {}
    series = {'mae_gain': gain, 'reference_mean': md.mean(axis=1),
              'prediction_mean': pred.mean(axis=1),
              'mae': np.abs(pred - md).mean(axis=1)}
    if md.shape[1] >= 3:
        a = rankdata(md, axis=1)
        b = rankdata(pred, axis=1)
        a -= a.mean(axis=1, keepdims=True)
        b -= b.mean(axis=1, keepdims=True)
        denominator = np.sqrt((a*a).sum(axis=1) * (b*b).sum(axis=1))
        rho = np.divide((a*b).sum(axis=1), denominator,
                        out=np.full(len(md), np.nan), where=denominator > 0)
        series['spearman'] = rho
    for name, values in series.items():
        finite = values[np.isfinite(values)]
        lo, hi = np.quantile(finite, [.025, .975])
        result[name + '_ci95_low'] = float(lo)
        result[name + '_ci95_high'] = float(hi)
    result['mae_gain_interval_positive'] = bool(result['mae_gain_ci95_low'] > 0)
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mapping_audit_path = ANNOTATION / 'independent_mapping_audit.json'
    recovery_audit_path = OUT / 'ile_reference_recovery_audit.json'
    mapping_audit = json.loads(mapping_audit_path.read_text(encoding='utf-8'))
    recovery_audit = json.loads(recovery_audit_path.read_text(encoding='utf-8'))
    assert mapping_audit['status'] == recovery_audit['status'] == 'passed'
    mapping_path = ANNOTATION / 'independent_imgt_author_mapping.csv'
    mapping = pd.read_csv(mapping_path, dtype={'author_id': str})
    mapping.author_id = mapping.author_id.astype(str)
    assert len(mapping) == 444 and not mapping.duplicated(KEY).any()
    assert mapping.region.value_counts().to_dict() == {
        'non_variable_remainder': 220, 'V_framework': 175,
        'CDR3': 25, 'CDR2': 13, 'CDR1': 11}
    mapping.to_csv(OUT / 'cdr_residue_mapping.csv', index=False)
    run_path = SOURCE / '1g4_run_paired_residues.csv'
    assert sha(run_path) == mapping_audit['md_residue_table_sha256']
    protected_hashes = recovery_audit['protected_old_files_sha256_before']
    assert protected_hashes == recovery_audit['protected_old_files_sha256_after']
    for protected_path, digest in protected_hashes.items():
        assert sha(protected_path) == digest, protected_path
    original = pd.read_csv(run_path, dtype={'author_id': str})
    assert len(original) == 17760
    run_key = ['state', 'run', 'chain', 'author_id', 'metric']
    assert not original.duplicated(run_key).any()
    corrected = original.copy()
    recovery_path = OUT / 'ile_reference_recovery.csv'
    assert sha(recovery_path) == recovery_audit['recovered_csv_sha256']
    assert recovery_audit['max_control_absolute_difference'] < 1e-10
    recovered = pd.read_csv(recovery_path, dtype={'author_id': str})
    assert len(recovered) == 280 and not recovered.duplicated(run_key[:-1]).any()
    assert set(recovered.amino_acid) == {'I'}
    corrected = corrected.merge(recovered.rename(columns={'reference': 'recovered_reference'})[
        run_key[:-1] + ['recovered_reference']], on=run_key[:-1], how='left', validate='many_to_one')
    fix = (corrected.metric == 'hchi') & corrected.recovered_reference.notna()
    assert fix.sum() == 280 and not corrected.loc[fix, 'evaluation_mask'].any()
    assert corrected.loc[fix, 'reference'].isna().all()
    corrected.loc[fix, 'reference'] = corrected.loc[fix, 'recovered_reference']
    corrected.loc[fix, 'evaluation_mask'] = True
    corrected = corrected.drop(columns='recovered_reference')
    np.testing.assert_array_equal(corrected.prediction, original.prediction)
    np.testing.assert_allclose(corrected.loc[~fix, 'reference'], original.loc[~fix, 'reference'], equal_nan=True)
    assert corrected[corrected.evaluation_mask].groupby(['state','run','metric']).size().value_counts().to_dict() == {444:20, 362:20}
    annotated = corrected.merge(mapping[KEY + ['region', 'pdb_chain', 'imgt_number', 'imgt_insertion', 'amino_acid']].rename(
        columns={'amino_acid':'annotation_amino_acid'}), on=KEY, validate='many_to_one')
    assert (annotated.amino_acid == annotated.annotation_amino_acid).all()
    annotated.to_csv(OUT / 'cdr_annotated_run_predictions_and_references.csv', index=False)
    priors_path = ROOT / 'artifacts/evaluation_atlas_hchi_v1/analysis/training_fold_aa_priors.csv'
    priors = pd.read_csv(priors_path)
    priors = priors[priors.fold == 0].copy()
    priors.metric = priors.metric.replace({'chi_entropy': 'hchi'})
    rng = np.random.default_rng(SEED)
    weights = {state:rng.multinomial(10, [.1]*10, size=N_BOOT)/10 for state in ['apo','holo']}
    delta_metrics, absolute_metrics, individual_metrics = [], [], []
    per_residue, coverage = [], []
    versions = [('legacy_as_reported', original), ('corrected_ILE_alias', corrected)]
    for version, data in versions:
        for metric in ['rmsf', 'hchi']:
            full = data[data.metric == metric]
            valid = full[full.evaluation_mask].copy()
            arrays = {}
            for state in ['apo','holo']:
                for col in ['reference','prediction']:
                    table = valid[valid.state == state].pivot(index='run', columns=KEY, values=col).sort_index(axis=1)
                    assert len(table) == 10 and not table.isna().any().any()
                    arrays[state,col] = table
            keys = arrays['apo','reference'].columns
            for table in arrays.values():
                assert table.columns.equals(keys)
            mm = mapping.set_index(KEY).loc[keys].reset_index()
            md = arrays['holo','reference'].mean(0).to_numpy() - arrays['apo','reference'].mean(0).to_numpy()
            pred = arrays['holo','prediction'].mean(0).to_numpy() - arrays['apo','prediction'].mean(0).to_numpy()
            boot_md = weights['holo'] @ arrays['holo','reference'].to_numpy() - weights['apo'] @ arrays['apo','reference'].to_numpy()
            boot_pred = weights['holo'] @ arrays['holo','prediction'].to_numpy() - weights['apo'] @ arrays['apo','prediction'].to_numpy()
            for name, mask in subsets(mm):
                assert mask.any(), (version, metric, name)
                q = score(md[mask], pred[mask])
                zero_mae = float(np.abs(md[mask]).mean())
                q.update(reference_version=version, metric=metric, region=name,
                         zero_change_mae=zero_mae, mae_gain=zero_mae-q['mae'],
                         mae_reduction_percent=100*(1-q['mae']/zero_mae),
                         sign_correct_count=int(np.sum(md[mask]*pred[mask] > 0)),
                         sign_agreement=float(np.mean(md[mask]*pred[mask] > 0)),
                         always_rigidify_sign_agreement=float(np.mean(md[mask] < 0)),
                         always_soften_sign_agreement=float(np.mean(md[mask] > 0)))
                q.update(bootstrap_intervals(boot_md[:,mask], boot_pred[:,mask]))
                delta_metrics.append(q)
                if version == 'legacy_as_reported':
                    continue
                for state in ['apo','holo']:
                    md_state = arrays[state,'reference'].to_numpy()[:,mask]
                    pr_state = arrays[state,'prediction'].to_numpy()[:,mask]
                    aa = mm.loc[mask, 'amino_acid']
                    prior = priors[priors.metric==metric].set_index('amino_acid').loc[aa, 'training_mean'].to_numpy()
                    base = dict(metric=metric, region=name, state=state)
                    state_score = score(md_state.mean(0), pr_state.mean(0))
                    state_score['aa_only_mae'] = float(abs(prior-md_state.mean(0)).mean())
                    absolute_metrics.append(dict(base, **state_score))
                    for index, run in enumerate(arrays[state,'reference'].index):
                        row = dict(base, run=int(run), **score(md_state[index], pr_state[index]))
                        row['aa_only_mae'] = float(abs(prior-md_state[index]).mean())
                        individual_metrics.append(row)
            if version == 'corrected_ILE_alias':
                residue = mm.copy()
                residue['metric'] = metric
                residue['md_delta'] = md
                residue['prediction_delta'] = pred
                residue['signed_error'] = pred-md
                residue['absolute_error'] = abs(pred-md)
                for state in ['apo','holo']:
                    for col in ['reference','prediction']:
                        residue[state+'_'+col+'_mean'] = arrays[state,col].mean(0).to_numpy()
                        residue[state+'_'+col+'_sd'] = arrays[state,col].std(0,ddof=1).to_numpy()
                per_residue.append(residue)
            for name, mask in subsets(mapping):
                ids = pd.MultiIndex.from_frame(mapping.loc[mask,KEY])
                n = int(keys.isin(ids).sum())
                coverage.append(dict(reference_version=version, metric=metric, region=name,
                                     structural_residues=int(mask.sum()), evaluated_residues=n,
                                     excluded_residues=int(mask.sum())-n))
    metrics = pd.DataFrame(delta_metrics)
    metrics.to_csv(OUT / 'cdr_delta_metrics_with_bootstrap.csv', index=False)
    pd.DataFrame(absolute_metrics).to_csv(OUT / 'cdr_absolute_state_mean_metrics.csv', index=False)
    individuals = pd.DataFrame(individual_metrics)
    individuals.to_csv(OUT / 'cdr_individual_run_absolute_metrics.csv', index=False)
    individuals.groupby(['metric','region','state']).agg(
        n_runs=('run','size'), mean_mae=('mae','mean'), sd_mae=('mae','std'),
        mean_spearman=('spearman','mean'), sd_spearman=('spearman','std'),
        mean_aa_only_mae=('aa_only_mae','mean')).reset_index().to_csv(OUT / 'cdr_individual_run_summary.csv', index=False)
    pd.concat(per_residue, ignore_index=True).to_csv(OUT / 'cdr_residue_delta_comparison.csv', index=False)
    pd.DataFrame(coverage).to_csv(OUT / 'cdr_coverage.csv', index=False)
    # Recover old global scores exactly before accepting any new analysis.
    old = pd.read_csv(SOURCE / '1g4_validation_metrics.csv')
    for metric in ['rmsf','hchi']:
        a = metrics[(metrics.reference_version=='legacy_as_reported') & (metrics.metric==metric) & (metrics.region=='all_TCR')].iloc[0]
        b = old[(old.task=='apo_to_holo_delta') & (old.metric==metric)].iloc[0]
        np.testing.assert_allclose(a[['mae','rmse','spearman']].to_numpy(float), b[['mae','rmse','spearman']].to_numpy(float), atol=1e-12)
    rmsf = metrics[metrics.metric=='rmsf']
    np.testing.assert_allclose(rmsf[rmsf.reference_version=='legacy_as_reported'].mae,
                               rmsf[rmsf.reference_version=='corrected_ILE_alias'].mae, atol=0)
    for protected_path, digest in protected_hashes.items():
        assert sha(protected_path) == digest, protected_path
    audit = dict(status='passed', system='1G4 only; no CD3', n_independent_protein_systems=1,
        n_apo_runs=10, n_holo_runs=10, seed=SEED, n_bootstrap=N_BOOT,
        source_hashes={str(path.relative_to(ROOT)):sha(path) for path in [mapping_path,run_path,recovery_path,priors_path,mapping_audit_path,recovery_audit_path,Path(__file__)]},
        protected_old_files_sha256=protected_hashes,
        no_new_training=True, no_prediction_changes=True, old_files_unchanged=True,
        recovered_ile_rows=280, all_reference_counts={'rmsf':444,'hchi':362},
        cdr_reference_counts={'rmsf':49,'hchi':37},
        primary_definition='IMGT via STCRDab 2BNR numbered structure and V-domain FASTA',
        source_urls=['https://opig.stats.ox.ac.uk/webapps/stcrdab-stcrpred/StrViewer?pdb=2bnr',
                     'https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html'],
        limitations=['Post hoc single-system regional re-evaluation, not unseen-family validation.',
            'Inputs remain fixed 10 ns MD snapshots, reference 20-100 ns; not MD-free experimental PDB evaluation.',
            'Apo/holo are not paired by run number; bootstrap independent state ensembles and whole runs.',
            'No independent-residue p values, no multiple-comparison significance claim.',
            'Bootstrap intervals are descriptive, not corrected for region selection/multiple comparisons.',
            'Absolute state means use ten predictions; individual-run metrics separately supplied.',
            'Small CDR samples make rank correlations unstable; n<3 not reported.',
            'Outside the annotated variable domain is not called framework.',
            'No arbitrary biological-importance weights and no new test-set calibration.'])
    (OUT / 'cdr_evaluation_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    show = metrics[(metrics.reference_version=='corrected_ILE_alias') & metrics.region.isin(
        ['all_TCR','all_CDR','CDR1','CDR2','CDR3','alpha_CDR3','beta_CDR3','V_framework','non_variable_remainder'])]
    print(show[['metric','region','n','spearman','mae','zero_change_mae','mae_reduction_percent',
                'reference_mean','prediction_mean','mae_gain_ci95_low','mae_gain_ci95_high']].to_string(index=False))


if __name__ == '__main__':
    warnings.filterwarnings('ignore', message='invalid value encountered in cast')
    main()
