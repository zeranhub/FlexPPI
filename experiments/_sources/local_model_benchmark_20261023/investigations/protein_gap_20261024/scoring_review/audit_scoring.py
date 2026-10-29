"""Read-only raw-output/scoring audit; no model forward or target modification."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from pathlib import Path
import datetime
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):
    return pd.read_csv(path, encoding='utf-8-sig', keep_default_na=False, float_precision='round_trip')
def result_dict(y,p):
    return dict(spearman=float(spearmanr(y,p).statistic), pearson=float(pearsonr(y,p).statistic), mae=float(np.mean(np.abs(y-p))))

truth_path = ROOT.parent / 'baselines/cohort/truth_residues.csv.gz'
truth = pd.read_csv(truth_path, keep_default_na=False, float_precision='round_trip',
    usecols=['sample_id','residue_index','amino_acid','chain','author_residue_id','target_rmsf','valid_rmsf'],
    dtype={'amino_acid':str,'chain':str,'author_residue_id':str})
records = read(ROOT/'cohort/records.csv')
members = read(ROOT/'protein_complex_benchmark/classification/membership_frozen.csv')
stored = read(ROOT/'protein_complex_benchmark/evaluation/per_record_metrics.csv')
catalog = read(ROOT/'cohort/existing_prediction_manifest.csv')
groups = {k:v for k,v in truth.groupby('sample_id',sort=False)}
frozen_proof = json.loads((ROOT/'evaluation/frozen_predictions_audit.json').read_text(encoding='utf-8'))
assert sha(ROOT/'evaluation/frozen_predictions.npz') == frozen_proof['cache_sha256']
assert sha(truth_path) == frozen_proof['truth_sha256']
assert sha(ROOT/'cohort/existing_prediction_manifest.csv') == frozen_proof['original_prediction_manifest_sha256']
with np.load(ROOT/'evaluation/frozen_predictions.npz') as f:
    frozen_full = f['Full'].copy()
with np.load(ROOT/'evaluation/aligned_predictions.npz') as f:
    arrays = {m:f[m].copy() for m in ['Full','backflip-2.1 (assembly)','backflip-2.1 (chains)']}
assert np.array_equal(frozen_full,arrays['Full'],equal_nan=True)

raw_checks = []
for cr in catalog[catalog.method=='Full'].itertuples():
    rows = groups[cr.sample_id]
    assert sha(cr.path)==cr.file_sha256
    with np.load(cr.path) as f:
        sl=slice(int(cr.row_start),int(cr.row_stop))
        assert int(cr.rmsf_column)==0 and int(cr.hchi_column)==1 and cr.array_key=='mean'
        assert f['mean'].shape[1]==2
        for k in ['amino_acid','chain','author_residue_id']:
            assert f[k][sl].astype(str).tolist()==rows[k].tolist(),(cr.sample_id,k)
        assert np.array_equal(f['target'][sl,0],rows.target_rmsf.values,equal_nan=True)
        assert np.array_equal(f['valid'][sl,0],rows.valid_rmsf.to_numpy(bool))
        assert np.array_equal(f['mean'][sl,0],arrays['Full'][rows.index],equal_nan=True)
        raw_checks.append({'sample_id':cr.sample_id,'raw_sha256':cr.file_sha256,'n':len(rows),'same_raw_sequence_chain_author_target_valid_and_mean':True})
assert len(raw_checks)==3063
(OUT/'full_raw_identity_checks.csv').write_text(pd.DataFrame(raw_checks).to_csv(index=False),encoding='utf-8')

back_counts={}
bf_root=ROOT/'backflip/results/backflip-2.1'
statuses={}
for line in (bf_root/'status.jsonl').read_text(encoding='utf-8').splitlines():
    s=json.loads(line);statuses[s['sample_id'],s['context']]=s
bf_checks=[]
for context,method in [('whole_complex','backflip-2.1 (assembly)'),('independent_chain','backflip-2.1 (chains)')]:
    count=0;residues=0
    for sample,rows in groups.items():
        path=bf_root/(sample+'__'+context+'.npz')
        status=statuses[sample,context]
        assert status['status']=='complete' and sha(path)==status['output_sha256']
        with np.load(path) as f:
            for k in ['chain','author_residue_id']:
                assert f[k].astype(str).tolist()==rows[k].tolist(),(sample,context,k)
            assert np.array_equal(f['target'],rows.target_rmsf.values,equal_nan=True)
            assert np.array_equal(f['valid'],rows.valid_rmsf.to_numpy(bool))
            assert np.array_equal(f['prediction'],arrays[method][rows.index],equal_nan=True)
            assert np.isfinite(f['prediction'][f['valid']]).all()
            residues+=int(f['valid'].sum())
        count+=1
        bf_checks.append({'sample_id':sample,'context':context,'sha256':status['output_sha256'],'same_chain_author_target_valid_and_prediction':True})
    back_counts[context]={'records':count,'valid_residues':residues}
pd.DataFrame(bf_checks).to_csv(OUT/'backflip2_raw_identity_checks.csv',index=False)

recomputed=[];max_metric_diff=0.
for sample,rows in groups.items():
    valid=rows.valid_rmsf.to_numpy(bool);y=rows.target_rmsf.to_numpy()[valid]
    for method,p in arrays.items():
        z=result_dict(y,p[rows.index][valid]);ref=stored[(stored.method==method)&(stored.sample_id==sample)].iloc[0]
        for metric in ['spearman','pearson','mae']:
            max_metric_diff=max(max_metric_diff,abs(z[metric]-float(ref[metric])))
        recomputed.append({'sample_id':sample,'method':method,'n':int(valid.sum()),**z})
recomputed=pd.DataFrame(recomputed)
recomputed.to_csv(OUT/'independently_recomputed_record_scores.csv',index=False)

member_group=members.set_index('sample_id').benchmark_group
record_meta=records.set_index('sample_id')
recomputed['group']=recomputed.sample_id.map(member_group)
recomputed['dataset']=recomputed.sample_id.map(record_meta.dataset)
recomputed['n_chains']=recomputed.sample_id.map(record_meta.n_chains)
truth['group']=truth.sample_id.map(member_group)
truth['dataset']=truth.sample_id.map(record_meta.dataset)
aggregation=[]
for group in ['protein','protein_complex']:
    for dataset in ['all','atlas','dynarepo','itrace']:
        mask=(truth.group==group)&(truth.valid_rmsf.astype(bool))
        per=recomputed[recomputed.group==group]
        if dataset!='all': mask &= truth.dataset==dataset; per=per[per.dataset==dataset]
        if not mask.any():continue
        for method,p in arrays.items():
            part=per[per.method==method]
            pooled=result_dict(truth.loc[mask,'target_rmsf'].to_numpy(),p[mask.values])
            aggregation.append({'group':group,'dataset':dataset,'method':method,'records':len(part),'residues':int(mask.sum()),
                'record_mean_spearman':float(part.spearman.mean()),'record_mean_pearson':float(part.pearson.mean()),'record_mean_mae':float(part.mae.mean()),
                'residue_weighted_record_spearman':float(np.average(part.spearman,weights=part.n)),
                'pooled_residue_spearman':pooled['spearman'],'pooled_residue_pearson':pooled['pearson'],'pooled_residue_mae':pooled['mae']})
pd.DataFrame(aggregation).to_csv(OUT/'aggregation_definition_comparison.csv',index=False)
composition=records.merge(members[['sample_id','benchmark_group']],on='sample_id',validate='one_to_one')
composition['chain_count_class']=np.where(composition.n_chains==1,'one_input_chain','multiple_input_chains')
composition=composition.groupby(['benchmark_group','dataset','chain_count_class'],sort=False).agg(records=('sample_id','size'),valid_residues=('n_valid_rmsf','sum')).reset_index()
composition.to_csv(OUT/'group_input_chain_composition.csv',index=False)

sources=[Path(__file__),ROOT/'evaluate_benchmark.py',ROOT/'evaluate_benchmark_initial.py',ROOT/'prepare_frozen_cache.py',
    ROOT/'evaluation/frozen_predictions_audit.json',ROOT/'evaluation/frozen_predictions.npz',ROOT/'evaluation/aligned_predictions.npz',
    ROOT/'cohort/existing_prediction_manifest.csv',ROOT/'cohort/records.csv',truth_path,
    ROOT/'backflip/run_benchmark.py',ROOT/'backflip/validate_outputs.py',ROOT/'backflip/prepare_inputs.py',
    bf_root/'status.jsonl',bf_root/'run_metadata.json',ROOT/'protein_complex_benchmark/evaluation/scoring_audit.json',
    ROOT/'protein_complex_benchmark/evaluation/per_record_metrics.csv',ROOT/'protein_complex_benchmark/classification/membership_frozen.csv']
audit={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'no_model_forward':True,'frozen_files_modified':False,
    'full_raw_records_checked':len(raw_checks),'backflip_raw_contexts':back_counts,'recomputed_scores':len(recomputed),
    'max_abs_difference_from_published_record_metrics':max_metric_diff,
    'truth_valid_value_counts':{str(k):int(v) for k,v in truth.valid_rmsf.value_counts().items()},
    'sources_sha256':{str(p):sha(p) for p in sources if p.exists()},
    'output_sha256':{p.name:sha(p) for p in OUT.glob('*.csv')}}
(OUT/'scoring_identity_audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'status':'read_only_scoring_audit_complete','raw_full':len(raw_checks),'raw_backflip':back_counts,'max_abs_metric_difference':max_metric_diff,'composition':composition.to_dict('records'),'aggregation':aggregation},ensure_ascii=True),flush=True)
