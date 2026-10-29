"""Second implementation: explicit negative-log F* and independent histograms.

Does not import the producer. Reconstructs every raw MD width and replayed curve.
"""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib
from datetime import datetime,timezone
import numpy as np
import pandas as pd

P=Path(__file__).resolve().parent
S=P.parents[1]/'free_energy_model_v26_20261013'/'experiment'
checks=0
elements=0

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):h.update(b)
    return h.hexdigest()

def equal(a,b,atol=0):
    global checks,elements
    a,b=np.asarray(a),np.asarray(b)
    assert a.shape==b.shape
    if atol:np.testing.assert_allclose(a,b,rtol=0,atol=atol)
    else:np.testing.assert_array_equal(a,b)
    checks+=1;elements+=a.size

def main():
    global checks
    plan=json.loads((P/'analysis_plan.json').read_text(encoding='utf8'))
    complete=json.loads((P/'completion_audit.json').read_text(encoding='utf8'))
    assert complete['status']=='complete'
    assert sha(P/'run_population_width_curves.py')==plan['analysis_code_sha256'];checks+=1
    for f,digest in plan['source_hash_files'].items():assert sha(f)==digest;checks+=1
    for f,item in complete['outputs'].items():assert sha(P/f)==item['sha256'];checks+=1
    data=np.load(P/'all_per_residue_widths.npz',allow_pickle=False)
    energies=np.array(plan['energy_cutoffs'])
    equal(energies,np.arange(31)/10)
    equal(data['energy_cutoffs'],energies)
    old=pd.read_csv(S/'replica1_per_residue.csv.gz',float_precision='round_trip',dtype={'sample_id':str,'chain':str,'author_residue_id':str})
    assert len(old)==33493 and int(old.eligible_sameAA.sum())==33142;checks+=1
    for name in ('sample_id','dataset','residue_index','amino_acid','chain','author_residue_id','chi_count','eligible_sameAA','full_hchi'):
        equal(data[name],old[name].to_numpy())
    equal(data['width_degrees'][:,10],old.width_Fstar_le_1p0_degrees.to_numpy())
    reconstructed=[]
    for r in plan['records']:
        selected=np.flatnonzero(data['sample_id']==r['sample_id'])
        source=[c for c in r['caches'] if c['replica']==1];assert len(source)==1
        assert sha(source[0]['path'])==source[0]['sha256'];checks+=1
        c=np.load(source[0]['path'],allow_pickle=False)
        idx=data['residue_index'][selected];q=data['chi_count'][selected]
        equal(c['chi_count'][idx],q)
        equal(c['times_ps'],np.arange(501)*200,atol=.05)
        for name in ('amino_acid','chain','author_residue_id'):equal(c[name][idx],data[name][selected])
        per_angle=[]
        for qi in range(4):
            valid=np.flatnonzero(q>qi)
            if not len(valid):continue
            angle=c['chi_degrees'][:,idx[valid],qi].astype(float)
            row_indices=np.broadcast_to(np.arange(len(valid)),angle.shape)
            for origin in plan['bin_origins_degrees']:
                code=np.floor(np.mod(angle+180-origin,360)/5).astype(int)
                histogram=np.zeros((len(valid),72),int)
                np.add.at(histogram,(row_indices.ravel(),code.ravel()),1)
                equal(histogram.sum(axis=1),np.full(len(valid),501))
                probability=histogram/501.
                peak=probability.max(axis=1)
                f=np.full(probability.shape,np.inf)
                occupied=probability>0
                ratio=probability/peak[:,None]
                f[occupied]=-np.log(ratio[occupied])
                vectors=5*np.sum(f[:,:,None]<=energies[None,None,:],axis=1)
                expanded=np.full((len(idx),31),np.nan,float);expanded[valid]=vectors
                per_angle.append(expanded)
        values=np.nanmean(np.stack(per_angle,axis=1),axis=1)
        equal(values,data['width_degrees'][selected])
        reconstructed.append((r['sample_id'],values))
        print('CACHE VERIFIED',r['sample_id'],flush=True)
    # The npz row order must match plan concatenation, not merely each identity.
    equal(np.concatenate([v for _,v in reconstructed]),data['width_degrees'])
    groups=pd.read_csv(P/'prediction_defined_groups.csv',float_precision='round_trip',dtype={'sample_id':str})
    original_groups=pd.read_csv(S/'replica1_high_low_per_AA.csv',float_precision='round_trip',dtype={'sample_id':str})
    original_groups=original_groups[(original_groups.threshold==1)&(original_groups.scope=='all')&(original_groups.predictor=='full_hchi')].reset_index(drop=True)
    # A one-threshold subset CSV infers threshold1 as int rather than float;
    # exact values/indices are required, serialization-only dtype is ignored.
    pd.testing.assert_frame_equal(groups,original_groups,check_exact=True,check_dtype=False);checks+=1
    record=pd.read_csv(P/'per_record_curves.csv',float_precision='round_trip',dtype={'sample_id':str})
    aa=pd.read_csv(P/'per_AA_group_curves.csv',float_precision='round_trip',dtype={'sample_id':str})
    ids=sorted(data['sample_id'].tolist());ids=sorted(set(ids))
    lows,highs=[],[]
    for sid in ids:
        loc=np.flatnonzero(data['sample_id']==sid)
        rowmap={int(ri):i for ri,i in zip(data['residue_index'][loc],loc)}
        lo_parts,hi_parts=[],[]
        for _,g in groups[groups.sample_id==sid].sort_values('amino_acid').iterrows():
            low=[int(v) for v in g.low_indices.split(';')];high=[int(v) for v in g.high_indices.split(';')]
            eligible=old[(old.sample_id==sid)&(old.amino_acid==g.amino_acid)&old.eligible_sameAA].sort_values(['full_hchi','residue_index'],kind='mergesort')
            equal(eligible.residue_index.iloc[:len(eligible)//2].to_numpy(),low)
            equal(eligible.residue_index.iloc[-(len(eligible)//2):].to_numpy(),high)
            lo=np.mean(data['width_degrees'][[rowmap[x] for x in low]],axis=0)
            hi=np.mean(data['width_degrees'][[rowmap[x] for x in high]],axis=0)
            saved=aa[(aa.sample_id==sid)&(aa.amino_acid==g.amino_acid)].sort_values('E')
            equal(saved.E.to_numpy(),energies)
            equal(saved.low_width_degrees.to_numpy(),lo)
            equal(saved.high_width_degrees.to_numpy(),hi)
            equal(saved.difference_degrees.to_numpy(),hi-lo)
            lo_parts.append(lo);hi_parts.append(hi)
        low=np.mean(lo_parts,axis=0);high=np.mean(hi_parts,axis=0)
        saved=record[record.sample_id==sid].sort_values('E')
        equal(saved.low_width_degrees.to_numpy(),low)
        equal(saved.high_width_degrees.to_numpy(),high)
        equal(saved.difference_degrees.to_numpy(),high-low)
        lows.append(low);highs.append(high)
    low=np.array(lows);high=np.array(highs);difference=high-low
    boot=np.load(P/'primary_paired_bootstrap_curves.npz',allow_pickle=False)
    equal(boot['record_ids'],ids)
    metadata=record[['sample_id','dataset']].drop_duplicates().set_index('sample_id').loc[ids]
    equal(boot['record_sources'],metadata.dataset.to_numpy(str))
    index=boot['record_indices'];assert index.shape==(10000,83) and index.min()==0 and index.max()==82;checks+=1
    rng=np.random.default_rng(plan['bootstrap']['seed']);index_parts=[]
    for ds,n in [('atlas',30),('dynarepo',23),('itrace',30)]:
        loc=np.flatnonzero(boot['record_sources']==ds)
        assert len(loc)==n;checks+=1
        equal((boot['record_sources'][index]==ds).sum(axis=1),np.full(10000,n))
        index_parts.append(loc[rng.integers(0,n,size=(10000,n))])
    equal(index,np.concatenate(index_parts,axis=1))
    summary=pd.read_csv(P/'population_curve_summary.csv',float_precision='round_trip')
    equal(summary.E.to_numpy(),energies)
    for label,values in [('low',low),('high',high),('difference',difference)]:
        replay=values[index].mean(axis=1)
        equal(replay,boot[label+'_mean_draws'])
        equal(summary[label+'_mean'].to_numpy(),values.mean(axis=0),atol=1e-12)
        ci=np.quantile(replay,[.025,.975],axis=0)
        equal(summary[label+'_bootstrap_95_low'].to_numpy(),ci[0])
        equal(summary[label+'_bootstrap_95_high'].to_numpy(),ci[1])
    equal(summary.positive_difference_records.to_numpy(),(difference>0).sum(axis=0))
    equal(summary.negative_difference_records.to_numpy(),(difference<0).sum(axis=0))
    equal(summary.zero_difference_records.to_numpy(),(difference==0).sum(axis=0))
    assert np.all(np.diff(low,axis=1)>=0) and np.all(np.diff(high,axis=1)>=0);checks+=1
    result=dict(status='passed',finished_utc=datetime.now(timezone.utc).isoformat(),verification_code_sha256=sha(__file__),
        producer_code_sha256=sha(P/'run_population_width_curves.py'),plan_sha256=sha(P/'analysis_plan.json'),checks=checks,compared_array_elements=elements,
        raw_caches_reconstructed=83,per_residue_positions=33493,cutoffs=31,all_model_defined_half_groups_exact=True,
        all_raw_width_vectors_exact=True,all_AA_and_record_curves_exact=True,all_saved_bootstrap_draws_replayed=True,
        bootstrap_source_counts=[30,23,30],all_pointwise_intervals_match=True,
        verification_scope='Second implementation by the experiment agent; uses explicit -ln(p/pmax) and np.add.at; no producer imports. This is not a separate-agent audit.',
        claim_scope='MD widths conditioned on frozen model ranking, not model-predicted F* curves, joint angular coupling, barriers or kinetics')
    (P/'second_implementation_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,ensure_ascii=False),flush=True)

if __name__=='__main__':main()
