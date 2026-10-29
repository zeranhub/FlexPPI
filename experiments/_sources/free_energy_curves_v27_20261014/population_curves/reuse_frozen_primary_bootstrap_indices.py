"""Reuse the already frozen v26 primary paired record indices for all W(E).

The original v27 seed20261014 plan and outputs are preserved. The same prior
evidence chain supplies indices, rather than choosing a new random seed by CI.
"""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
import numpy as np
import pandas as pd

P=Path(__file__).resolve().parent
S=P.parents[1]/'free_energy_model_v26_20261013'/'experiment'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):h.update(b)
    return h.hexdigest()

def write_json(path,obj):Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')

def main():
    delivery=json.loads((S/'experiment_delivery_audit.json').read_text(encoding='utf8'))
    frozen_npz=S/'replica1_high_low_primary_bootstrap_draws.npz'
    original_record=S/'replica1_high_low_per_record.csv'
    original_summary=S/'replica1_high_low_summary.csv'
    for f in (frozen_npz,original_record,original_summary):assert sha(f)==delivery['files'][f.name]['sha256']
    source_files={str(f):sha(f) for f in (frozen_npz,original_record,original_summary,S/'analysis_plan.json',S/'experiment_delivery_audit.json',P/'analysis_plan.json',P/'per_record_curves.csv')}
    addendum=dict(protocol='Reuse_v26_original_primary_paired_indices_for_v27_all31_curves',created_utc=datetime.now(timezone.utc).isoformat(),
        source_hash_files=source_files,code_sha256=sha(__file__),
        rationale='Keep figure curve uncertainty on the same frozen v26 primary evidence-chain indices; no selection by interval quality; original v27 seed20261014 plan and results retained untouched',
        bootstrap='Exactly copy v26 record_ids and record_indices (seed20261013); apply paired indices to both MD group curves across all31 thresholds; pointwise percentile95%, not simultaneous',
        E1_anchor='Verify new curve low/high/difference per record equal old v26 within1e-12; use frozen v26 high_minus_low_MD at E=1 so its original AA-difference averaging and floating-point ordering are preserved',
        no_new_inputs='No caches re-read, MD, training, inference, regrouping or threshold selection')
    addendum_path=P/'frozen_primary_indices_addendum.json'
    if addendum_path.exists():
        previous=json.loads(addendum_path.read_text(encoding='utf8'))
        assert previous['source_hash_files']==source_files and previous['code_sha256']==sha(__file__)
        addendum=previous
    else:write_json(addendum_path,addendum)
    frozen=np.load(frozen_npz,allow_pickle=False)
    ids=frozen['record_ids'];indices=frozen['record_indices']
    curves=pd.read_csv(P/'per_record_curves.csv',float_precision='round_trip',dtype={'sample_id':str})
    energies=np.arange(31)/10
    low=curves.pivot(index='sample_id',columns='E',values='low_width_degrees').loc[ids,energies].to_numpy()
    high=curves.pivot(index='sample_id',columns='E',values='high_width_degrees').loc[ids,energies].to_numpy()
    difference=high-low
    old=pd.read_csv(original_record,float_precision='round_trip',dtype={'sample_id':str})
    old=old[(old.threshold==1)&(old.scope=='all')&(old.predictor=='full_hchi')].set_index('sample_id').loc[ids]
    np.testing.assert_allclose(low[:,10],old.low_MD_mean,rtol=0,atol=1e-12)
    np.testing.assert_allclose(high[:,10],old.high_MD_mean,rtol=0,atol=1e-12)
    np.testing.assert_allclose(difference[:,10],old.high_minus_low_MD,rtol=0,atol=1e-12)
    anchor_rounding_error=float(np.max(abs(difference[:,10]-old.high_minus_low_MD.to_numpy())))
    difference[:,10]=old.high_minus_low_MD.to_numpy()
    assert indices.shape==(10000,83)
    sources=old.dataset.to_numpy(str)
    for ds,n in [('atlas',30),('dynarepo',23),('itrace',30)]:
        np.testing.assert_array_equal((sources[indices]==ds).sum(axis=1),np.full(10000,n))
    low_draws=low[indices].mean(axis=1);high_draws=high[indices].mean(axis=1);difference_draws=difference[indices].mean(axis=1)
    np.testing.assert_allclose(difference_draws[:,10],frozen['full_hchi'],rtol=0,atol=1e-12)
    replay_rounding_error=float(np.max(abs(difference_draws[:,10]-frozen['full_hchi'])))
    # Preserve the exact already published primary floating-point draw values.
    difference_draws[:,10]=frozen['full_hchi']
    np.savez_compressed(P/'frozen_primary_paired_bootstrap_curves.npz',energy_cutoffs=energies,record_ids=ids,record_sources=sources,record_indices=indices,
        low_mean_draws=low_draws,high_mean_draws=high_draws,difference_mean_draws=difference_draws,
        seed=int(frozen['seed']),resamples=int(frozen['resamples']),source_bootstrap_sha256=np.array(sha(frozen_npz)),
        frozen_primary_anchor_difference_by_record=difference[:,10],addendum_sha256=np.array(sha(addendum_path)))
    rows=[]
    for ei,e in enumerate(energies):
        row=dict(E=e,n_records=83,positive_difference_records=int((difference[:,ei]>0).sum()),negative_difference_records=int((difference[:,ei]<0).sum()),zero_difference_records=int((difference[:,ei]==0).sum()))
        for name,values,draws in [('low',low,low_draws),('high',high,high_draws),('difference',difference,difference_draws)]:
            ci=np.quantile(draws[:,ei],[.025,.975])
            row[name+'_mean']=float(values[:,ei].mean());row[name+'_bootstrap_95_low']=float(ci[0]);row[name+'_bootstrap_95_high']=float(ci[1])
        rows.append(row)
    summary=pd.DataFrame(rows)
    summary.to_csv(P/'population_curve_summary_frozen_primary_indices.csv',index=False,float_format='%.17g')
    old_summary=pd.read_csv(original_summary,float_precision='round_trip')
    old_summary=old_summary[(old_summary.threshold==1)&(old_summary.scope=='all')&(old_summary.predictor=='full_hchi')&(old_summary.dataset=='all')].iloc[0]
    anchor=summary[summary.E==1].iloc[0]
    np.testing.assert_array_equal([anchor.difference_mean,anchor.difference_bootstrap_95_low,anchor.difference_bootstrap_95_high],
        [old_summary['mean'],old_summary.bootstrap_95_low,old_summary.bootstrap_95_high])
    for f,digest in source_files.items():assert sha(f)==digest
    outputs={f.name:dict(sha256=sha(f),bytes=f.stat().st_size) for f in [addendum_path,P/'frozen_primary_paired_bootstrap_curves.npz',P/'population_curve_summary_frozen_primary_indices.csv',Path(__file__)]}
    write_json(P/'frozen_primary_indices_audit.json',dict(status='complete',finished_utc=datetime.now(timezone.utc).isoformat(),
        code_sha256=sha(__file__),addendum_sha256=sha(addendum_path),source_hashes_checked_before_and_after=True,
        all_original_indices_exact=True,all_original_record_ids_exact=True,source_counts_each_draw=[30,23,30],all31cutoffs=31,
        E1_pre_anchor_record_rounding_max_abs=anchor_rounding_error,E1_pre_anchor_draw_rounding_max_abs=replay_rounding_error,
        E1_original_delta_draws_exact=True,E1_original_delta_mean_and_CI_exact=True,
        E1=dict(low_mean=float(anchor.low_mean),high_mean=float(anchor.high_mean),difference_mean=float(anchor.difference_mean),
            difference_bootstrap_95_low=float(anchor.difference_bootstrap_95_low),difference_bootstrap_95_high=float(anchor.difference_bootstrap_95_high)),outputs=outputs))
    print('FROZENPRIMARY E1',anchor.to_json(),flush=True)

if __name__=='__main__':main()
