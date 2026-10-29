"""Identity-checked, record-equal benchmark. Never calibrate on test labels."""
import os
for name in ['OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS']:
    os.environ[name] = '2'
import argparse
import hashlib
import json
from pathlib import Path
import time
import warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'baselines'
OUT = ROOT / 'evaluation'
OUT.mkdir(exist_ok=True)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def stats(y, p, absolute):
    result = dict(spearman=np.nan, pearson=np.nan, z_mae=np.nan, top10_ap=np.nan,
                  mae=np.nan, rmse=np.nan, bias=np.nan)
    if absolute:
        result.update(mae=float(np.mean(abs(y-p))), rmse=float(np.sqrt(np.mean((y-p)**2))), bias=float(np.mean(p-y)))
    if len(y) >= 4 and np.ptp(y) > 0 and np.ptp(p) > 0:
        result.update(spearman=float(spearmanr(y, p).statistic), pearson=float(pearsonr(y, p).statistic),
            z_mae=float(np.mean(abs((y-y.mean())/y.std(ddof=0)-(p-p.mean())/p.std(ddof=0)))))
    if np.ptp(y) > 0:
        result['top10_ap'] = float(average_precision_score(y >= np.quantile(y, .90), p))
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bootstrap', type=int, default=2000)
    args = parser.parse_args()
    records = pd.read_csv(ROOT/'cohort/records.csv', keep_default_na=False)
    truth_path = SOURCE/'cohort/truth_residues.csv.gz'
    frozen = json.loads((SOURCE/'cohort/cohort_freeze.json').read_text())
    assert sha(truth_path) == frozen['files']['truth_residues.csv.gz']
    truth = pd.read_csv(truth_path, keep_default_na=False, float_precision='round_trip',
        dtype={'chain':str, 'author_residue_id':str, 'amino_acid':str})
    assert len(records)==3063 and len(truth)==961321
    assert records.main_group_candidate.value_counts().to_dict()=={'single_chain':2193,'other_multichain':626,'tcr_pmhc_itrace':244}
    assert not truth.duplicated(['sample_id','chain','author_residue_id']).any()
    grouped = {name: rows for name,rows in truth.groupby('sample_id', sort=False)}
    catalog = pd.read_csv(ROOT/'cohort/existing_prediction_manifest.csv', keep_default_na=False)
    methods, metadata, source_audit = {}, {}, []
    cache = {}
    def register(name, absolute, provenance, context):
        methods[name]=np.full(len(truth),np.nan)
        metadata[name]=dict(absolute_angstrom=absolute,provenance=provenance,input_context=context)
    for name in ['Full','Intra','Static','AA']:
        register(name,True,'frozen_five_fold_oof_reuse', 'complete_structure' if name in ['Full','Static'] else 'within_chain' if name=='Intra' else 'sequence_residue_type')
        for row in catalog[catalog.method==name].itertuples():
            rows=grouped[row.sample_id]
            # A source hash is checked exactly once; immutable source arrays are then reused.
            if row.path not in cache:
                assert sha(row.path)==row.file_sha256
                with np.load(row.path) as f:
                    cache[row.path]={key:f[key] for key in [row.array_key,'target','valid'] if key in f}
            obj=cache[row.path]
            arr=obj[row.array_key][int(row.row_start):int(row.row_stop),int(row.rmsf_column)]
            assert arr.shape==(len(rows),)
            if 'target' in obj:
                y=obj['target'][int(row.row_start):int(row.row_stop),0]
                assert np.array_equal(y,rows.target_rmsf.to_numpy(),equal_nan=True), row.sample_id
            methods[name][rows.index]=arr
        cache.clear()
        print('Collected',name,flush=True)
    # New local independent-chain DynaProt inference, exact input sequence and output hash.
    register('DynaProt-M',True,'new_local_cpu','independent_chain')
    directory=ROOT/'dynaprot/predictions'
    for ident,rows in truth.groupby('input_id',sort=False):
        path=directory/(ident+'.npz'); status=directory/(ident+'.json')
        if path.exists() and status.exists():
            audit=json.loads(status.read_text())
            assert audit['status']=='complete' and sha(path)==audit['output_sha256']
            with np.load(path) as f:
                assert str(f['sequence'].item())==''.join(rows.sort_values('input_index').amino_acid)
                value=f['rmsf_angstrom']
                assert value.shape==(len(rows),) and np.isfinite(value).all()
                methods['DynaProt-M'][rows.index]=value[rows.input_index.to_numpy(int)]
    print('Collected DynaProt-M',flush=True)
    # Original ANM/GNM full-position arrays: retained explicitly as verified historical caches.
    path=SOURCE/'prody_common_evaluation/common_residue_predictions.csv.gz'
    if path.exists():
        old=pd.read_csv(path,keep_default_na=False,float_precision='round_trip',dtype={'chain':str,'author_residue_id':str})
        assert old[['sample_id','residue_index','chain','author_residue_id']].equals(truth[['sample_id','residue_index','chain','author_residue_id']])
        assert np.array_equal(old.target_rmsf,truth.target_rmsf)
        for name,key in [('ANM','anm'),('GNM','gnm')]:
            register(name,False,'verified_frozen_official_prody_cache','complete_structure')
            methods[name][:]=old[key]
        source_audit.append(dict(path=str(path),sha256=sha(path)))
        del old
    # Fresh BackFlip run, with all original author residues including insertion codes retained.
    back=ROOT/'backflip/results'
    for variant in ['backflip-2.1','backflip-1.0','backflip-1.0-seq']:
        directory=back/variant
        if not directory.exists():
            continue
        for context,suffix in [('whole_complex','assembly'),('independent_chain','chains')]:
            name=f'{variant} ({suffix})'
            register(name,True,'new_local_cpu',context)
            for sample,rows in grouped.items():
                path=directory/(sample+'__'+context+'.npz')
                if not path.exists(): continue
                with np.load(path) as f:
                    assert f['chain'].astype(str).tolist()==rows.chain.tolist(), sample
                    assert f['author_residue_id'].astype(str).tolist()==rows.author_residue_id.tolist(),sample
                    assert np.array_equal(f['target'],rows.target_rmsf.to_numpy(),equal_nan=True),sample
                    value=f['prediction']; assert value.shape==(len(rows),)
                    methods[name][rows.index]=value
    coverage,per=[],[]
    for name,p in methods.items():
        for r in records.itertuples():
            rows=grouped[r.sample_id]
            valid=rows.valid_rmsf.to_numpy(bool)
            values=p[rows.index]
            complete=bool(np.isfinite(values[valid]).all())
            coverage.append(dict(method=name,sample_id=r.sample_id,group=r.main_group_candidate,expected=int(valid.sum()),covered=int(np.isfinite(values[valid]).sum()),complete=complete))
            if not complete: continue
            per.append(dict(method=name,sample_id=r.sample_id,group=r.main_group_candidate,dataset=r.dataset,
                family=r.family_component_id,n=int(valid.sum()),**stats(rows.target_rmsf.to_numpy(float)[valid],values[valid],metadata[name]['absolute_angstrom'])))
        print('Scored',name,flush=True)
    cover=pd.DataFrame(coverage)
    scores=pd.DataFrame(per)
    cover.to_csv(OUT/'coverage_per_record.csv',index=False)
    csum=cover.groupby(['method','group'],sort=False).agg(records=('sample_id','size'),complete_records=('complete','sum'),expected_residues=('expected','sum'),covered_residues=('covered','sum')).reset_index()
    csum.to_csv(OUT/'coverage_summary.csv',index=False)
    scores.to_csv(OUT/'per_record_metrics.csv',index=False)
    complete_methods=[n for n in methods if cover.loc[cover.method==n,'complete'].all()]
    summaries=[]
    rng=np.random.default_rng(20261023)
    for group in ['single_chain','other_multichain','tcr_pmhc_itrace','all_multichain','all']:
        rows=scores if group=='all' else scores[scores.group!='single_chain'] if group=='all_multichain' else scores[scores.group==group]
        for name in methods:
            part=rows[rows.method==name]
            if part.empty: continue
            for metric in ['spearman','pearson','z_mae','top10_ap','mae','rmse','bias']:
                ok=part[np.isfinite(part[metric])]
                if ok.empty: continue
                point=float(ok[metric].mean());lo=hi=np.nan
                fam=ok.groupby('family')[metric].agg(['sum','count'])
                if len(fam)>1 and args.bootstrap:
                    index=rng.integers(0,len(fam),size=(args.bootstrap,len(fam)))
                    draws=fam['sum'].to_numpy()[index].sum(1)/fam['count'].to_numpy()[index].sum(1)
                    lo,hi=np.quantile(draws,[.025,.975])
                summaries.append(dict(method=name,group=group,metric=metric,mean=point,ci_low=lo,ci_high=hi,n_records=len(ok),n_families=len(fam),primary_complete_cohort=name in complete_methods,
                    uncertainty='record_weighted_family_cluster_bootstrap' if len(fam)>1 else 'descriptive_single_family_no_independent_CI'))
    pd.DataFrame(summaries).to_csv(OUT/'group_summary.csv',index=False)
    # Store aligned full-length outputs for reanalysis, never only the model-wise favorable intersection.
    np.savez_compressed(OUT/'aligned_predictions.npz',**methods)
    audit=dict(records=len(records),residues=len(truth),group_confirmed=True,complete_methods=complete_methods,methods=metadata,
        truth_sha256=sha(truth_path),script_sha256=sha(__file__),additional_source_hashes=source_audit,
        metric_policy='All frozen valid RMSF residues; score only complete records; equal record mean; no fitted calibration; AP uses >=90th reference percentile and retains ties.',
        ci_policy='2000 family-cluster bootstrap draws, ratios of record sums/counts. One-family strata get no independent CI. Published weights may have training overlap.',
        pointwise_hchi_policy='External flexibility models lack the same joint side-chain state entropy output; not scored as Hchi.',
        no_primary_incomplete_methods=True)
    (OUT/'scoring_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print('Complete full-cohort methods:',complete_methods,flush=True)
    print(pd.DataFrame(summaries).query("group=='all' and metric=='spearman'").to_string(index=False),flush=True)

if __name__=='__main__':
    main()
