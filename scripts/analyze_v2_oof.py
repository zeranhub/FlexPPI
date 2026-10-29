"""Build full-coverage OOF source data and a separately labelled cached-baseline bridge."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from export_v2_oof_predictions import ROOT, METRICS, load_cache, sha256, atomic_json
from flexppi.comparison import _system_metrics
from flexppi.evaluation import _correlation, _spearman


def regression(prediction, target):
    p, y = np.asarray(prediction, float), np.asarray(target, float)
    if not len(y) or not np.isfinite(p).all() or not np.isfinite(y).all():
        raise ValueError("Empty/nonfinite regression input")
    e = p - y
    return dict(mae=float(np.abs(e).mean()), rmse=float(np.sqrt((e*e).mean())),
                bias=float(e.mean()), pearson=_correlation(p, y), spearman=_spearman(p, y))


def read_oof(folder):
    summary = json.loads((folder / 'run_summary.json').read_text())
    if summary['status'] != 'complete' or summary['partial'] or summary['processed_samples'] != 3068:
        raise ValueError('Refusing figures/analysis from incomplete OOF run')
    manifest = [json.loads(x) for x in (folder / 'prediction_manifest.jsonl').read_text().splitlines()]
    frames = []
    for row in manifest:
        if row['status'] != 'complete':
            raise ValueError('Incomplete sample: '+row['sample_id'])
        a = load_cache(Path(row['cache_path']), row['identity'])
        r = row['record']
        frame = pd.DataFrame(dict(sample_id=row['sample_id'], dataset=r['dataset'],
                                  fold=row['fold'], chain=a['chain'], author_residue_id=a['author_residue_id'],
                                  residue_index=a['residue_index'], amino_acid=a['amino_acid'],
                                  biological_system=r['metadata']['biological_system'],
                                  cluster_id=r['metadata']['cluster_id'], source_group=r['source_group']))
        for i, metric in enumerate(METRICS):
            frame['valid_'+metric] = a['valid'][:, i]
            frame['target_'+metric] = a['target'][:, i]
            frame['pred_'+metric] = a['mean'][:, i]
            frame['model_variance_'+metric] = a['model_variance'][:, i]
            frame['label_variance_'+metric] = a['label_variance'][:, i]
        frames.append(frame)
    data = pd.concat(frames, ignore_index=True)
    keys = ['sample_id', 'chain', 'author_residue_id']
    if data.duplicated(keys).any():
        raise ValueError('Duplicate OOF residue identities')
    if data.valid_rmsf.sum() != 959347 or data.valid_chi_entropy.sum() != 431023:
        raise ValueError('OOF valid residue coverage mismatch')
    return data, manifest


def analyze_states(data, out):
    priors = []
    data['chi_type_prior'] = np.nan
    for fold in range(5):
        train = data.loc[(data.fold != fold) & data.valid_chi_entropy]
        means = train.groupby('amino_acid').target_chi_entropy.mean()
        # Float64 summation matches the training-fold evaluation definition.
        means = train.assign(y=train.target_chi_entropy.astype(float)).groupby('amino_acid').y.mean()
        holdout = (data.fold == fold) & data.valid_chi_entropy
        data.loc[holdout, 'chi_type_prior'] = data.loc[holdout, 'amino_acid'].map(means)
        for aa, value in means.items():
            priors.append(dict(fold=fold, amino_acid=aa, training_mean=value,
                               training_residues=int((train.amino_acid == aa).sum())))
        server = json.loads((ROOT / f'artifacts/performance_first_v2/evaluation/fold_{fold}_test.json').read_text())
        recorded = server['metrics']['chi_entropy']['amino_acid_corrected']['training_priors']
        if any(abs(float(means[aa])-v) > 1e-6 for aa, v in recorded.items()):
            raise ValueError('Recomputed training-only amino-acid priors disagree with server')
    if data.loc[data.valid_chi_entropy, 'chi_type_prior'].isna().any():
        raise ValueError('Missing training prior for a test amino-acid type')
    pd.DataFrame(priors).to_csv(out/'training_fold_aa_priors.csv', index=False)
    rows, chain_rows = [], []
    for sid, group in data.groupby('sample_id', sort=False):
        base = {k:group.iloc[0][k] for k in ['sample_id','dataset','fold','biological_system','cluster_id']}
        for metric in METRICS:
            valid = group.loc[group['valid_'+metric]]
            if valid.empty:
                continue
            p, y = valid['pred_'+metric], valid['target_'+metric]
            rows.append(dict(**base, metric=metric, n_residues=len(valid), **regression(p,y)))
            for chain, part in valid.groupby('chain', sort=False):
                chain_rows.append(dict(**base, chain=chain, metric=metric, n_residues=len(part),
                                       **regression(part['pred_'+metric],part['target_'+metric])))
            if metric == 'chi_entropy':
                prior = valid.chi_type_prior
                for name, pp, yy in [('chi_type_prior',prior,y),('chi_residual',p-prior,y-prior)]:
                    rows.append(dict(**base, metric=name,n_residues=len(valid), **regression(pp, yy)))
    systems = pd.DataFrame(rows)
    systems.to_csv(out/'per_system_metrics.csv',index=False)
    pd.DataFrame(chain_rows).to_csv(out/'per_chain_metrics.csv',index=False)
    stat = ['mae','rmse','bias','pearson','spearman']
    folds = systems.groupby(['dataset','metric','fold'])[stat].mean().reset_index()
    folds.to_csv(out/'fold_metrics.csv',index=False)
    byfold = folds.groupby(['dataset','metric'])[stat].agg(['mean','std'])
    byfold.to_csv(out/'five_fold_mean_sd_ddof1.csv')
    aa_rows = []
    for (fold, aa), part in data.loc[data.valid_chi_entropy].groupby(['fold','amino_acid']):
        aa_rows.append(dict(fold=fold,amino_acid=aa,n_residues=len(part),
                            **regression(part.pred_chi_entropy,part.target_chi_entropy)))
    pd.DataFrame(aa_rows).to_csv(out/'within_aa_metrics.csv',index=False)
    # Verify the independently exported CPU predictions against server fold reports.
    checks = []
    for fold in range(5):
        report = json.loads((ROOT/f'artifacts/performance_first_v2/evaluation/fold_{fold}_test.json').read_text())
        for metric in METRICS:
            s = systems.loc[(systems.fold==fold)&(systems.metric==metric)]
            ref = report['metrics'][metric]['overall']['macro_by_system']
            for key in ['mae','rmse','pearson','spearman']:
                diff = float(s[key].mean()-ref[key])
                checks.append(dict(fold=fold,metric=metric,statistic=key,cpu_minus_server=diff))
                if abs(diff)>2e-4:
                    raise ValueError(f'CPU/server mismatch: {fold} {metric} {key}: {diff}')
    pd.DataFrame(checks).to_csv(out/'cpu_server_reproduction.csv',index=False)
    data.to_csv(out/'oof_per_residue.csv.gz', index=False, compression='gzip')
    return systems


BASELINE_COLUMNS = {
    'FlexPPI-v2':'pred_rmsf', 'FlexPPI-old':'flexppi_pretrain_rmsf_angstrom',
    'FastProtFlex':'fastprotflex_score', 'BackFlip-seq':'backflip_seq_rmsf_angstrom',
    'BackFlip-structure':'backflip_structure_rmsf_angstrom',
    'ProDy-ANM':'anm_sq_fluct','ProDy-GNM':'gnm_sq_fluct','RMSF-net-PDB':'rmsf_net_pdb_score'}


def baseline_bridge(data, manifest, out):
    lookup = {r['sample_id'].lower():r for r in manifest}
    if len(lookup)!=3068:
        raise ValueError('Case-normalised sample ID collision')
    source = ROOT/'artifacts/evaluation/additional_baselines'
    old_inputs = {}
    for name in ['static','itrace']:
        m = pd.read_csv(source/f'manifests/{name}.tsv',sep='\t')
        for _, r in m.iterrows():
            sid = str(r.sample_id).lower()
            old_inputs['itrace_'+sid if name=='itrace' else sid] = str(r.pdb_path)
    input_checks = []
    for sid, path in old_inputs.items():
        actual=sha256(Path(path)); expected=lookup[sid]['identity']['input_sha256']
        if actual!=expected:
            raise ValueError('Cached-baseline input differs from v2 input: '+sid)
        input_checks.append(dict(sample_id=sid,input_sha256=actual))
    comparisons=[]; audits=[]; tables=[]
    columns=['sample_id','chain','author_residue_id','fold','biological_system','cluster_id',
             'target_rmsf','pred_rmsf','valid_rmsf']
    new=data[columns].copy();new.sample_id=new.sample_id.str.lower()
    for name in ['static','itrace']:
        path=source/f'results/{name}_common_residue_predictions.csv'
        old=pd.read_csv(path,dtype={'sample_id':str,'pdb_id':str,'author_residue_id':str,
                                  'residue_number':str,'insertion_code':str},keep_default_na=False)
        if name=='itrace':
            old['sample_id']='itrace_'+old.pdb_id.str.lower()
            old['author_residue_id']=old.residue_number+old.insertion_code
            old['cohort']='itrace_oof_bridge'
            target='target_chainfit_rmsf_angstrom'
        else:
            old.sample_id=old.sample_id.str.lower();target='target_rmsf_angstrom'
        old.chain=old.chain.replace('','_')
        keys=['sample_id','chain','author_residue_id']
        merged=old.merge(new,on=keys,how='left',validate='one_to_one',indicator=True)
        if not (merged._merge=='both').all() or not merged.valid_rmsf.all():
            raise ValueError('Baseline bridge misses a valid OOF residue')
        required=list(BASELINE_COLUMNS.values())+['target_rmsf']
        if not np.isfinite(merged[required].to_numpy(float)).all():
            raise ValueError('Baseline bridge has missing/nonfinite predictions')
        delta=np.abs(merged[target].to_numpy(float)-merged.target_rmsf.to_numpy(float))
        audits.append(dict(cohort=name,samples=int(merged.sample_id.nunique()),residues=len(merged),
                           old_target_differences_above_1e_6=int((delta>1e-6).sum()),
                           old_target_max_difference=float(delta.max()),source_sha256=sha256(path)))
        merged.drop(columns='_merge').to_csv(out/f'{name}_common_predictions.csv.gz',index=False,compression='gzip')
        for (cohort,sid), group in merged.groupby(['cohort','sample_id'],sort=False):
            chain_groups = [part for _,part in group.groupby('chain',sort=False)] if (
                cohort in ['dynarepo_complex','itrace_oof_bridge']) else [group]
            for method,col in BASELINE_COLUMNS.items():
                chain_metrics=[]
                for part in chain_groups:
                    values=_system_metrics(part[col].to_numpy(float),part.target_rmsf.to_numpy(float),
                                           absolute_scale=method.startswith(('FlexPPI','BackFlip')))
                    if not values:
                        raise ValueError(f'Undefined benchmark metrics: {sid} {method}')
                    chain_metrics.append(values)
                keys_metrics=chain_metrics[0].keys()
                result={key:float(np.mean([r[key] for r in chain_metrics if r.get(key) is not None]))
                        for key in keys_metrics if any(r.get(key) is not None for r in chain_metrics)}
                comparisons.append(dict(cohort=cohort,sample_id=sid,fold=int(group.fold.iloc[0]),
                                         biological_system=group.biological_system.iloc[0],method=method,
                                         n_residues=len(group),n_chains=len(chain_groups),**result))
    comp=pd.DataFrame(comparisons)
    comp.to_csv(out/'bridge_per_system_metrics.csv',index=False)
    metrics=[c for c in comp if c not in ['cohort','sample_id','fold','biological_system','method','n_residues','n_chains']]
    comp.groupby(['cohort','method'])[metrics].mean().to_csv(out/'bridge_macro_metrics.csv')
    comp.groupby(['cohort','method','fold'])[metrics].mean().to_csv(out/'bridge_fold_metrics.csv')
    atomic_json(out/'bridge_audit.json',dict(status='complete_legacy_cohort_bridge',cohorts=audits,
                input_hash_checks=input_checks, total_systems=520, full_pool_systems=3068,
                limitation='Existing baseline predictions reused, not rerun. Current source structures match v2 inputs; historical per-input hashes were not universally recorded. Full 3068 baseline coverage is pending.'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=ROOT/'artifacts/evaluation_v2/oof')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'artifacts/evaluation_v2/analysis')
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    data,manifest=read_oof(args.input_dir)
    print('Loaded full OOF:',len(data),'residues',flush=True)
    systems=analyze_states(data,args.output_dir)
    print('Full OOF statistics verified against server; updating legacy-cohort bridge',flush=True)
    baseline_bridge(data,manifest,args.output_dir)
    atomic_json(args.output_dir/'analysis_summary.json',dict(status='complete',samples=3068,
                rmsf_residues=int(data.valid_rmsf.sum()),hchi_residues=int(data.valid_chi_entropy.sum()),
                bridge_samples=520,baseline_full_pool_complete=False,
                fold_sd='sample SD, ddof=1; descriptive variation, not independent repetitions or CI',
                estimands='whole-system full OOF; chain-then-system legacy bridge, reported separately',
                input_run_sha256=sha256(args.input_dir/'run_summary.json')))
    print('Analysis complete',flush=True)


if __name__=='__main__':
    main()
