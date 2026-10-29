"""Replay frozen v26 indices against every curve point, without producer imports."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib
from datetime import datetime,timezone
import numpy as np
import pandas as pd

P=Path(__file__).resolve().parent
S=P.parents[1]/'free_energy_model_v26_20261013'/'experiment'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    a=json.loads((P/'frozen_primary_indices_addendum.json').read_text(encoding='utf8'))
    for f,digest in a['source_hash_files'].items():assert sha(f)==digest
    assert sha(P/'reuse_frozen_primary_bootstrap_indices.py')==a['code_sha256']
    source=np.load(S/'replica1_high_low_primary_bootstrap_draws.npz',allow_pickle=False)
    saved=np.load(P/'frozen_primary_paired_bootstrap_curves.npz',allow_pickle=False)
    for key in ('record_ids','record_indices','seed','resamples'):np.testing.assert_array_equal(saved[key],source[key])
    curves=pd.read_csv(P/'per_record_curves.csv',float_precision='round_trip',dtype={'sample_id':str})
    ids=source['record_ids'];indices=source['record_indices'];e=np.arange(31)/10
    low=curves.pivot(index='sample_id',columns='E',values='low_width_degrees').loc[ids,e].to_numpy()
    high=curves.pivot(index='sample_id',columns='E',values='high_width_degrees').loc[ids,e].to_numpy()
    diff=high-low
    original=pd.read_csv(S/'replica1_high_low_per_record.csv',float_precision='round_trip',dtype={'sample_id':str})
    original=original[(original.threshold==1)&(original.scope=='all')&(original.predictor=='full_hchi')].set_index('sample_id').loc[ids]
    np.testing.assert_allclose(diff[:,10],original.high_minus_low_MD.to_numpy(),rtol=0,atol=1e-12)
    diff[:,10]=original.high_minus_low_MD.to_numpy()
    np.testing.assert_array_equal(saved['frozen_primary_anchor_difference_by_record'],diff[:,10])
    np.testing.assert_array_equal(saved['record_sources'],original.dataset.to_numpy(str))
    for ds,n in [('atlas',30),('dynarepo',23),('itrace',30)]:
        np.testing.assert_array_equal((saved['record_sources'][indices]==ds).sum(axis=1),np.full(10000,n))
    summary=pd.read_csv(P/'population_curve_summary_frozen_primary_indices.csv',float_precision='round_trip')
    np.testing.assert_array_equal(summary.E.to_numpy(),e)
    for name,values in [('low',low),('high',high),('difference',diff)]:
        replay=values[indices].mean(axis=1)
        if name=='difference':
            np.testing.assert_allclose(replay[:,10],source['full_hchi'],rtol=0,atol=1e-12)
            replay[:,10]=source['full_hchi']
        np.testing.assert_array_equal(saved[name+'_mean_draws'],replay)
        np.testing.assert_array_equal(summary[name+'_mean'].to_numpy(),np.array([values[:,ei].mean() for ei in range(31)]))
        ci=np.quantile(replay,[.025,.975],axis=0)
        np.testing.assert_array_equal(summary[name+'_bootstrap_95_low'].to_numpy(),ci[0])
        np.testing.assert_array_equal(summary[name+'_bootstrap_95_high'].to_numpy(),ci[1])
    np.testing.assert_array_equal(saved['difference_mean_draws'][:,10],source['full_hchi'])
    old=pd.read_csv(S/'replica1_high_low_summary.csv',float_precision='round_trip')
    old=old[(old.threshold==1)&(old.scope=='all')&(old.predictor=='full_hchi')&(old.dataset=='all')].iloc[0]
    anchor=summary[summary.E==1].iloc[0]
    np.testing.assert_array_equal([anchor.difference_mean,anchor.difference_bootstrap_95_low,anchor.difference_bootstrap_95_high],
        [old['mean'],old.bootstrap_95_low,old.bootstrap_95_high])
    for name,condition in [('positive',diff>0),('negative',diff<0),('zero',diff==0)]:
        np.testing.assert_array_equal(summary[name+'_difference_records'].to_numpy(),condition.sum(axis=0))
    result=dict(status='passed',finished_utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        source_original_bootstrap_sha256=sha(S/'replica1_high_low_primary_bootstrap_draws.npz'),
        all_indices_and_record_ids_exact=True,source_counts_each_draw=[30,23,30],all31points_replayed=True,
        paired_bootstrap_arrays=3,draws_per_curve=10000,compared_draw_elements=930000,
        E1_original_delta_draws_and_mean_and_CI_exact=True,source_hashes_checked=True,
        verification_scope='Separate replay implementation by same experiment agent; no imports from producer/adapter; original v27 plan and outputs preserved')
    (P/'frozen_primary_indices_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
