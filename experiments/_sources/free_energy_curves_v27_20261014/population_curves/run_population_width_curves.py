"""Prediction-defined population curves of MD low-F* angular width.

Only this directory is written. Fixed v26 Full Hchi same-AA halves are reused.
The curves are measured from MD, not predicted free-energy curves.
"""
from __future__ import annotations
import sys
sys.dont_write_bytecode=True
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,argparse
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'free_energy_model_v26_20261013'/'experiment'
E=np.arange(31,dtype=float)/10
ORIGINS=(0.,1.25,2.5,3.75)
SEED=20261014
B=10000

def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as h:
        for block in iter(lambda:h.read(4*1024*1024),b''):result.update(block)
    return result.hexdigest()

def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.generic):return clean(x.item())
    if isinstance(x,float) and not np.isfinite(x):return None
    return x

def write_json(path,x):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(clean(x),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')

def utc():return datetime.now(timezone.utc).isoformat()

def freeze():
    delivery=json.loads((SOURCE/'experiment_delivery_audit.json').read_text(encoding='utf8'))
    assert delivery['status']=='complete'
    inputs={str(SOURCE/'experiment_delivery_audit.json'):sha(SOURCE/'experiment_delivery_audit.json')}
    names=['analysis_plan.json','independent_verification.json','replica1_per_residue.csv.gz',
        'replica1_high_low_per_AA.csv','replica1_high_low_per_record.csv','replica1_high_low_summary.csv']
    for name in names:
        assert sha(SOURCE/name)==delivery['files'][name]['sha256']
        inputs[str(SOURCE/name)]=delivery['files'][name]['sha256']
    old=json.loads((SOURCE/'analysis_plan.json').read_text(encoding='utf8'))
    for path,digest in old['source_hash_files'].items():
        assert sha(path)==digest
        inputs[path]=digest
    plan=dict(protocol='Fixed_prediction_halves_MD_low_Fstar_width_population_curves_v1',created_utc=utc(),
        analysis_code_sha256=sha(__file__),source_plan_sha256=sha(SOURCE/'analysis_plan.json'),source_hash_files=inputs,
        records=old['records'],sample_counts=old['sample_counts'],
        energy_cutoffs=E.tolist(),main_anchor_cutoff=1.,bin_width_degrees=5.,bin_origins_degrees=list(ORIGINS),histogram_bins=72,
        metric='W(E): total measure in degrees of occupied per-angle bins satisfying F*=-ln(p/pmax)<=E; average four circular origins and all chemically defined chi angles including Pro',
        metric_semantics='Union of qualifying bins, possibly disconnected; measured MD marginal angular support; not a connected basin width or predicted free-energy function',
        hist_rule='No smoothing, interpolation or pseudocount; p=counts/501; select count>0 and count>=countmax*exp(-E)',
        trajectory_grid=dict(start_ns=0,end_ns=100,step_ps=200,frames=501,replica=1),
        common_positions=33493,eligible_sameAA_positions=33142,
        grouping='Reuse exact v26 threshold1/scopeall/predictorfull_hchi same-AA prediction-defined low/high indices; each eligible type has>=4 positions; top/bottomfloor(n/2), tie canonical index; odd middle excluded only from halves; no regrouping using MD or cutoff',
        expected_half_positions=dict(low=16250,high=16250,omitted_middle=642),
        weights='Within each record, average equally across original eligible amino-acid groups; then equal record weights across all83 records for each E',
        bootstrap=dict(seed=SEED,resamples=B,unit='record',method='Source-stratified paired resampling with original counts30atlas/23dynarepo/30itrace; same record indices for high/low/difference and all31 cutoff points; pointwise percentile95% intervals, not simultaneous band'),
        preserve_all='Keep all31 cutoffs, all83 records and all33493 per-residue width vectors; no threshold selection based on curve result',
        required_anchor_validation='At E=1 reproduce each v26 per-residue width, each original AA group low/high, each record low/high/difference, and population means within1e-12 absolute serialization tolerance',
        models='Full Hchi scores define groups; v26 Intra and RMSF comparisons remain frozen elsewhere, no new model inference',
        freeze_boundary='Analysis plan written before any new W(E) or population result; all sources checked before and after; no MD/training/inference')
    dest=HERE/'analysis_plan.json'
    if dest.exists():
        previous=json.loads(dest.read_text(encoding='utf8'))
        assert previous['analysis_code_sha256']==sha(__file__)
        assert previous['source_hash_files']==inputs and previous['records']==plan['records']
        plan=previous
    else:write_json(dest,plan)
    return plan

def width_vectors(cache_path,base):
    cache=np.load(cache_path,allow_pickle=False)
    idx=base.residue_index.to_numpy(int);q=base.chi_count.to_numpy(int)
    np.testing.assert_allclose(cache['times_ps'],np.arange(501)*200,rtol=0,atol=.05)
    np.testing.assert_array_equal(q,cache['chi_count'][idx])
    for name in ('amino_acid','chain','author_residue_id'):
        np.testing.assert_array_equal(base[name].astype(str).to_numpy(),cache[name][idx].astype(str))
    assert cache['chi_valid'][idx].all() and cache['ca_valid'][idx].all()
    w=np.full((len(idx),4,4,31),np.nan,float)
    for qi in range(4):
        valid=np.flatnonzero(q>qi)
        if not len(valid):continue
        angles=cache['chi_degrees'][:,idx[valid],qi].astype(float)
        assert np.isfinite(angles).all()
        offsets=np.arange(len(valid),dtype=np.int64)*72
        for oi,origin in enumerate(ORIGINS):
            code=np.floor(((angles+180-origin)%360)/5).astype(np.int64)
            assert code.min()>=0 and code.max()<72
            counts=np.bincount((code+offsets[None,:]).ravel(),minlength=len(valid)*72).reshape(len(valid),72)
            assert np.all(counts.sum(axis=1)==501)
            bound=counts.max(axis=1)[:,None,None]*np.exp(-E)[None,None,:]
            selected=(counts[:,:,None]>0)&(counts[:,:,None]>=bound)
            w[valid,qi,oi,:]=5*selected.sum(axis=1)
    values=np.nanmean(w.reshape(len(idx),16,31),axis=1)
    assert np.isfinite(values).all() and np.all(np.diff(values,axis=1)>=0)
    assert values.min()>=5 and values.max()<=360
    return values

def compute(plan):
    residue=pd.read_csv(SOURCE/'replica1_per_residue.csv.gz',dtype={'sample_id':str,'chain':str,'author_residue_id':str},float_precision='round_trip')
    old_group=pd.read_csv(SOURCE/'replica1_high_low_per_AA.csv',dtype={'sample_id':str},float_precision='round_trip')
    group=old_group[(old_group.threshold==1)&(old_group.scope=='all')&(old_group.predictor=='full_hchi')].copy()
    assert group.rankable.all() and group.sample_id.nunique()==83
    assert len(residue)==33493 and int(residue.eligible_sameAA.sum())==33142
    assert group.low_positions.sum()==16250 and group.high_positions.sum()==16250 and group.omitted_middle_positions.sum()==642
    group.to_csv(HERE/'prediction_defined_groups.csv',index=False,float_format='%.17g')
    old_record=pd.read_csv(SOURCE/'replica1_high_low_per_record.csv',dtype={'sample_id':str},float_precision='round_trip')
    old_record=old_record[(old_record.threshold==1)&(old_record.scope=='all')&(old_record.predictor=='full_hchi')].set_index('sample_id')
    arrays,identities,record_rows,aa_rows=[],[],[],[]
    anchor_max_error=0.
    for number,record in enumerate(plan['records'],start=1):
        sid=record['sample_id'];base=residue[residue.sample_id==sid].copy().reset_index(drop=True)
        caches=[c for c in record['caches'] if c['replica']==1];assert len(caches)==1
        values=width_vectors(caches[0]['path'],base)
        expected=base.width_Fstar_le_1p0_degrees.to_numpy()
        np.testing.assert_allclose(values[:,10],expected,rtol=0,atol=1e-12)
        anchor_max_error=max(anchor_max_error,float(np.max(abs(values[:,10]-expected))))
        arrays.append(values);identities.append(base)
        lookup={int(index):i for i,index in enumerate(base.residue_index)}
        aa_low,aa_high=[],[]
        selected_groups=group[group.sample_id==sid].sort_values('amino_acid')
        for _,g in selected_groups.iterrows():
            low=[int(x) for x in g.low_indices.split(';')];high=[int(x) for x in g.high_indices.split(';')]
            assert len(low)==g.low_positions and len(high)==g.high_positions and not set(low)&set(high)
            # Verify frozen identities and model-only ranking independently.
            part=base[(base.amino_acid==g.amino_acid)&base.eligible_sameAA].sort_values(['full_hchi','residue_index'],kind='mergesort')
            half=len(part)//2
            assert low==part.residue_index.iloc[:half].tolist() and high==part.residue_index.iloc[-half:].tolist()
            lo=values[[lookup[i] for i in low]].mean(axis=0);hi=values[[lookup[i] for i in high]].mean(axis=0)
            np.testing.assert_allclose([lo[10],hi[10],hi[10]-lo[10]],[g.low_MD_mean,g.high_MD_mean,g.high_minus_low_MD],rtol=0,atol=1e-12)
            aa_low.append(lo);aa_high.append(hi)
            for ei,e in enumerate(E):
                aa_rows.append(dict(sample_id=sid,dataset=record['dataset'],fold=record['fold'],amino_acid=g.amino_acid,chi_count=int(g.chi_count),E=e,
                    low_width_degrees=lo[ei],high_width_degrees=hi[ei],difference_degrees=hi[ei]-lo[ei],low_positions=len(low),high_positions=len(high)))
        lo=np.mean(aa_low,axis=0);hi=np.mean(aa_high,axis=0)
        orig=old_record.loc[sid]
        np.testing.assert_allclose([lo[10],hi[10],hi[10]-lo[10]],[orig.low_MD_mean,orig.high_MD_mean,orig.high_minus_low_MD],rtol=0,atol=1e-12)
        for ei,e in enumerate(E):
            record_rows.append(dict(sample_id=sid,dataset=record['dataset'],fold=record['fold'],E=e,
                low_width_degrees=lo[ei],high_width_degrees=hi[ei],difference_degrees=hi[ei]-lo[ei],AA_groups=len(aa_low),
                low_positions=int(selected_groups.low_positions.sum()),high_positions=int(selected_groups.high_positions.sum()),
                omitted_middle_positions=int(selected_groups.omitted_middle_positions.sum())))
        print('RECORD',number,83,sid,'E1MATCH',flush=True)
    all_frame=pd.concat(identities,ignore_index=True);all_values=np.concatenate(arrays,axis=0)
    np.savez_compressed(HERE/'all_per_residue_widths.npz',energy_cutoffs=E,width_degrees=all_values,
        sample_id=all_frame.sample_id.to_numpy(str),dataset=all_frame.dataset.to_numpy(str),fold=all_frame.fold.to_numpy(int),
        residue_index=all_frame.residue_index.to_numpy(int),amino_acid=all_frame.amino_acid.to_numpy(str),chain=all_frame.chain.to_numpy(str),
        author_residue_id=all_frame.author_residue_id.to_numpy(str),chi_count=all_frame.chi_count.to_numpy(int),
        eligible_sameAA=all_frame.eligible_sameAA.to_numpy(bool),full_hchi=all_frame.full_hchi.to_numpy(float),
        plan_sha256=np.array(sha(HERE/'analysis_plan.json')))
    aa=pd.DataFrame(aa_rows);rec=pd.DataFrame(record_rows)
    aa.to_csv(HERE/'per_AA_group_curves.csv',index=False,float_format='%.17g')
    rec.to_csv(HERE/'per_record_curves.csv',index=False,float_format='%.17g')
    record_ids=sorted(rec.sample_id.unique());record_meta=rec[['sample_id','dataset']].drop_duplicates().set_index('sample_id').loc[record_ids]
    low=rec.pivot(index='sample_id',columns='E',values='low_width_degrees').loc[record_ids,E].to_numpy()
    high=rec.pivot(index='sample_id',columns='E',values='high_width_degrees').loc[record_ids,E].to_numpy()
    difference=high-low
    rng=np.random.default_rng(SEED);index_parts=[]
    for ds in ('atlas','dynarepo','itrace'):
        loc=np.flatnonzero(record_meta.dataset.to_numpy()==ds)
        assert len(loc)=={'atlas':30,'dynarepo':23,'itrace':30}[ds]
        index_parts.append(loc[rng.integers(0,len(loc),size=(B,len(loc)))])
    indices=np.concatenate(index_parts,axis=1)
    low_draw=low[indices].mean(axis=1);high_draw=high[indices].mean(axis=1);diff_draw=difference[indices].mean(axis=1)
    np.testing.assert_allclose(high_draw-low_draw,diff_draw,rtol=0,atol=1e-12)
    np.savez_compressed(HERE/'primary_paired_bootstrap_curves.npz',energy_cutoffs=E,record_ids=np.array(record_ids),record_sources=record_meta.dataset.to_numpy(str),
        record_indices=indices,low_mean_draws=low_draw,high_mean_draws=high_draw,difference_mean_draws=diff_draw,seed=SEED,resamples=B)
    summary=[]
    for ei,e in enumerate(E):
        row=dict(E=e,n_records=83,positive_difference_records=int((difference[:,ei]>0).sum()),negative_difference_records=int((difference[:,ei]<0).sum()),zero_difference_records=int((difference[:,ei]==0).sum()))
        for label,v,d in [('low',low,low_draw),('high',high,high_draw),('difference',difference,diff_draw)]:
            ci=np.quantile(d[:,ei],[.025,.975])
            row[label+'_mean']=v[:,ei].mean();row[label+'_bootstrap_95_low']=ci[0];row[label+'_bootstrap_95_high']=ci[1]
        summary.append(row)
    summary=pd.DataFrame(summary)
    summary.to_csv(HERE/'population_curve_summary.csv',index=False,float_format='%.17g')
    anchor=summary[summary.E==1].iloc[0].to_dict()
    np.testing.assert_allclose([anchor['low_mean'],anchor['high_mean'],anchor['difference_mean']],
        [old_record.low_MD_mean.mean(),old_record.high_MD_mean.mean(),old_record.high_minus_low_MD.mean()],rtol=0,atol=1e-12)
    for path,digest in plan['source_hash_files'].items():assert sha(path)==digest
    outputs={str(p.relative_to(HERE)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in HERE.rglob('*') if p.is_file() and p.name!='completion_audit.json'}
    write_json(HERE/'completion_audit.json',dict(status='complete',finished_utc=utc(),plan_sha256=sha(HERE/'analysis_plan.json'),code_sha256=sha(__file__),
        records=83,trajectories_read=83,all31cutoffs_saved=True,all33493positions_saved=True,eligible_sameAA_positions=33142,
        original_group_indices_exact=True,source_hashes_checked_before_and_after=True,E1_per_residue_max_abs_error=anchor_max_error,E1_record_groups_match=True,
        E1_population_match=True,E1=anchor,bootstrap_pointwise_only=True,outputs=outputs))
    print('COMPLETE E1',json.dumps(clean(anchor)),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=['freeze','compute'],default='compute');args=parser.parse_args()
    plan=freeze();print('FROZEN',sha(HERE/'analysis_plan.json'),plan['created_utc'],flush=True)
    if args.phase=='compute':compute(plan)

if __name__=='__main__':main()
