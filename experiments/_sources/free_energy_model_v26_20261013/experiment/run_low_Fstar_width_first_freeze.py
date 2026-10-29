"""Frozen, OOF-prediction-led validation against MD low-F* angular width.

Writes only this experiment directory; all previous outputs and sources are read-only.
No MD, training, inference, smoothing or outcome-based subset selection.
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
import json,hashlib,time,argparse
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'model_argument_v24_20261013'/'experiment_predictive_occupancy'
ORIGINS=(0.,1.25,2.5,3.75)
THRESHOLDS=(.5,1.,2.)
PREDICTORS=('full_hchi','intrachain_hchi','full_rmsf')
SCOPES=('all','q1','q2','q3','q4')
SEED=20261013
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
    old_plan_path=SOURCE/'analysis_plan.json'
    old=json.loads(old_plan_path.read_text(encoding='utf8'))
    source_files={str(old_plan_path):sha(old_plan_path)}
    per_residue_sources={}
    for mode in ('replica1','available_replicas'):
        audit_path=SOURCE/(mode+'_analysis_audit.json')
        audit=json.loads(audit_path.read_text(encoding='utf8'))
        assert audit['status']=='complete' and audit['records']==83
        item=audit['outputs']['per_residue']
        assert sha(item['path'])==item['sha256'] and item['rows']==33493
        source_files[str(audit_path)]=sha(audit_path)
        source_files[item['path']]=item['sha256']
        per_residue_sources[mode]=item
    for path,digest in old['source_hash_files'].items():
        assert sha(path)==digest
        source_files[path]=digest
    # Hash verification and identity inventory happen before any new metric.
    for r in old['records']:
        for field in ('full_prediction','intra_prediction','old_label'):
            assert sha(r[field])==r[field+'_sha256']
            source_files[r[field]]=r[field+'_sha256']
        for c in r['caches']:
            assert sha(c['path'])==c['sha256']
            source_files[c['path']]=c['sha256']
    plan=dict(protocol='OOF_sameAA_low_Fstar_width_v1',created_utc=utc(),analysis_code_sha256=sha(__file__),
        parent_plan_sha256=sha(old_plan_path),records=old['records'],per_residue_sources=per_residue_sources,
        source_hash_files=source_files,sample_counts=old['sample_counts'],
        primary_metric='Mean per-angle angular width in degrees of all occupied 5-degree bins satisfying F*=-ln(p/pmax)<=1; bin-union measure, not one connected basin width; each angle averaged over four circular origins then averaged over all chemically defined chi angles',
        primary_threshold=1.,secondary_thresholds=[.5,2.],thresholds=list(THRESHOLDS),
        exact_bin_rule='count>0 and -ln(count/countmax)<=threshold, equivalently count>=countmax*exp(-threshold); no pseudocount, smoothing or joining of disconnected bins',
        bin_width_degrees=5,bin_origins_degrees=list(ORIGINS),histogram_bins=72,
        all_angles_definition='All chemically defined chi angles including Pro; periodic +/-180 seam wrapped; no three-rotamer assumption',
        primary_replicas='replica1 from each of 83 records; 501 frames, 0-100 ns inclusive at 200 ps',
        repeat_sensitivity='Compute metric independently in all 203 cached trajectories and average metric over available replicas per original residue (ATLAS/ITrace 3, DynaRepo 1)',
        common_positions='Exact 33493 original v24 common positions, source identities/predictions/validity preserved; original same-AA>=4 rule retains 33142',
        type_control='Within each record subtract each eligible amino-acid group mean separately from MD width and predictor; Spearman over all retained residual positions; retain all original eligible types/residues',
        predictors=list(PREDICTORS),scopes=list(SCOPES),summary_unit='Equal-weight arithmetic mean across defined per-record quantities; record is resampling unit',
        high_low_rule='Sort each predictor within each eligible same-AA group ascending, tie by canonical residue_index; bottom/top floor(n/2), omit odd middle only here; record average of eligible-AA-equal high mean minus low mean in actual MD width; constant prediction group excluded as unrankable with count',
        high_low_scope='Overall and q1/q2/q3/q4; actual MD high and low means also reported',
        paired_comparisons=['full_hchi_minus_intrachain_hchi','full_hchi_minus_full_rmsf'],
        amino_acid_type_baseline='Constant within same-AA group, so centered rho is undefined and not zero; no same-AA ranking possible',
        bootstrap=dict(seed=SEED,resamples=B,unit='record',method='Source-stratified paired resampling retaining original source counts for each finite record set; same sampled record indices applied to all compared predictors; percentile 2.5/97.5%; main threshold=1/all-scope/all-source draws preserved'),
        frozen_boundary='Plan written before any width/correlation/high-low result; thresholds and metric fixed; no new MD/training/inference; weak and negative results retained',
        output_columns=dict(summary=['replica_mode','threshold','scope','predictor','dataset','n_records','mean','bootstrap_95_low','bootstrap_95_high','positive_records','negative_records','zero_records']))
    destination=HERE/'analysis_plan.json'
    if destination.exists():
        existing=json.loads(destination.read_text(encoding='utf8'))
        assert existing['analysis_code_sha256']==sha(__file__)
        assert existing['records']==plan['records'] and existing['source_hash_files']==plan['source_hash_files']
        plan=existing
    else:write_json(destination,plan)
    return plan

def spearman(x,y):
    x,y=np.asarray(x,float),np.asarray(y,float)
    assert np.isfinite(x).all() and np.isfinite(y).all()
    if len(x)<4 or np.ptp(x)<=1e-12 or np.ptp(y)<=1e-12:return np.nan
    return float(np.corrcoef(pd.Series(x).rank(method='average'),pd.Series(y).rank(method='average'))[0,1])

def width_cache(c,base,record):
    cache=np.load(c['path'],allow_pickle=False)
    np.testing.assert_allclose(cache['times_ps'],np.arange(501)*200,atol=.05,rtol=0)
    idx=base.residue_index.to_numpy(int)
    q=base.chi_count.to_numpy(int)
    for name in ('amino_acid','chain','author_residue_id'):
        np.testing.assert_array_equal(base[name].astype(str).to_numpy(),cache[name][idx].astype(str))
    np.testing.assert_array_equal(q,cache['chi_count'][idx])
    assert cache['chi_valid'][idx].all() and cache['ca_valid'][idx].all()
    widths=np.full((len(idx),4,4,3),np.nan,float)
    maximum=np.zeros((len(idx),4,4),np.int16)
    for qi in range(4):
        valid=np.flatnonzero(q>qi)
        if not len(valid):continue
        angles=cache['chi_degrees'][:,idx[valid],qi].astype(float)
        assert np.isfinite(angles).all()
        offsets=np.arange(len(valid),dtype=np.int64)*72
        for oi,origin in enumerate(ORIGINS):
            bins=np.floor(((angles+180-origin)%360)/5).astype(np.int64)
            assert bins.min()>=0 and bins.max()<72
            counts=np.bincount((bins+offsets[None,:]).ravel(),minlength=len(valid)*72).reshape(len(valid),72)
            assert np.all(counts.sum(axis=1)==501)
            maxcount=counts.max(axis=1)
            maximum[valid,qi,oi]=maxcount
            for ti,threshold in enumerate(THRESHOLDS):
                keep=(counts>0)&(counts>=maxcount[:,None]*np.exp(-threshold))
                widths[valid,qi,oi,ti]=5*keep.sum(axis=1)
    assert np.all(widths[...,0]<=widths[...,1]) or np.isnan(widths).any()
    for qi in range(4):
        valid=q>qi
        assert np.isfinite(widths[valid,qi,:,:]).all()
        assert np.all(np.diff(widths[valid,qi,:,:],axis=-1)>=0)
    means=np.nanmean(widths.reshape(len(idx),16,3),axis=1)
    assert np.isfinite(means).all() and np.all(means>=5) and np.all(means<=360)
    dest=HERE/'per_angle_widths'/f"{record['sample_id']}_rep{c['replica']}.npz"
    dest.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(dest,residue_index=idx,chi_count=q,width_degrees=widths,max_count=maximum,
        bin_origins_degrees=ORIGINS,thresholds=THRESHOLDS,metric_mean_per_chi=means,
        source_cache_sha256=np.array(c['sha256']),plan_sha256=np.array(sha(HERE/'analysis_plan.json')))
    return means,dict(path=str(dest),sha256=sha(dest),positions=len(idx),replica=c['replica'])

def record_stats(frame,threshold,mode):
    idcols=dict(sample_id=frame.sample_id.iloc[0],dataset=frame.dataset.iloc[0],fold=int(frame.fold.iloc[0]),replica_mode=mode,threshold=threshold)
    eligible=frame.groupby('amino_acid').size()
    part=frame[frame.amino_acid.isin(eligible[eligible>=4].index)].copy()
    corr,groups,masks=[],[],[]
    for scope in SCOPES:
        p=part if scope=='all' else part[part.chi_count==int(scope[-1])]
        masks.append(idcols|dict(scope=scope,eligible_positions=len(p),AA_types=p.amino_acid.nunique()))
        if len(p)<4:continue
        cy=p.width_degrees-p.groupby('amino_acid').width_degrees.transform('mean')
        for pred in PREDICTORS:
            cx=p[pred]-p.groupby('amino_acid')[pred].transform('mean')
            corr.append(idcols|dict(scope=scope,predictor=pred,positions=len(p),AA_types=p.amino_acid.nunique(),rho_type_centered=spearman(cx,cy)))
            for aa,a in p.groupby('amino_acid',sort=True):
                ordered=a.sort_values([pred,'residue_index'],kind='mergesort')
                n=len(ordered)//2
                constant=np.ptp(ordered[pred].to_numpy())<=1e-12
                lo,hi=ordered.iloc[:n],ordered.iloc[-n:]
                groups.append(idcols|dict(scope=scope,predictor=pred,amino_acid=aa,chi_count=int(a.chi_count.iloc[0]),positions=len(a),low_positions=n,high_positions=n,
                    omitted_middle_positions=len(a)-2*n,rankable=not constant,
                    low_MD_mean=float(lo.width_degrees.mean()) if not constant else np.nan,
                    high_MD_mean=float(hi.width_degrees.mean()) if not constant else np.nan,
                    high_minus_low_MD=float(hi.width_degrees.mean()-lo.width_degrees.mean()) if not constant else np.nan,
                    low_prediction_mean=float(lo[pred].mean()) if not constant else np.nan,
                    high_prediction_mean=float(hi[pred].mean()) if not constant else np.nan,
                    low_indices=';'.join(map(str,lo.residue_index)) if not constant else '',
                    high_indices=';'.join(map(str,hi.residue_index)) if not constant else ''))
    return corr,groups,masks

def summary_with_pairs(frame,value,mode,kind):
    rng=np.random.default_rng(SEED)
    rows,pairs,paired_rows,draw_data=[],[],[],{}
    for (threshold,scope),part in frame.groupby(['threshold','scope'],sort=True):
        wide=part.pivot(index=['sample_id','dataset'],columns='predictor',values=value)
        for source in ('all','atlas','dynarepo','itrace'):
            subset=wide if source=='all' else wide[wide.index.get_level_values('dataset')==source]
            # All compared predictors share eligibility; assert equal finite record masks.
            masks=subset[list(PREDICTORS)].notna()
            assert (masks.eq(masks.iloc[:,0],axis=0)).all().all()
            subset=subset.dropna(subset=list(PREDICTORS)).sort_index()
            if not len(subset):continue
            idx_parts=[]
            for ds in sorted(subset.index.get_level_values('dataset').unique()):
                loc=np.flatnonzero(subset.index.get_level_values('dataset')==ds)
                idx_parts.append(loc[rng.integers(0,len(loc),size=(B,len(loc)))])
            indices=np.concatenate(idx_parts,axis=1)
            arr=subset[list(PREDICTORS)].to_numpy()
            draws=arr[indices].mean(axis=1)
            for pi,pred in enumerate(PREDICTORS):
                v=arr[:,pi];d=draws[:,pi];low,high=np.quantile(d,[.025,.975])
                rows.append(dict(replica_mode=mode,threshold=threshold,scope=scope,predictor=pred,dataset=source,n_records=len(v),mean=v.mean(),
                    bootstrap_95_low=low,bootstrap_95_high=high,record_sample_sd=pd.Series(v).std(),median=np.median(v),min=v.min(),max=v.max(),
                    positive_records=int((v>0).sum()),negative_records=int((v<0).sum()),zero_records=int((v==0).sum())))
                if threshold==1 and scope=='all' and source=='all':draw_data[pred]=d
            for bi,baseline in enumerate(PREDICTORS[1:],start=1):
                contrast='full_hchi_minus_'+baseline
                v=arr[:,0]-arr[:,bi];d=draws[:,0]-draws[:,bi];low,high=np.quantile(d,[.025,.975])
                paired_rows.append(dict(replica_mode=mode,threshold=threshold,scope=scope,contrast=contrast,dataset=source,n_records=len(v),mean=v.mean(),
                    bootstrap_95_low=low,bootstrap_95_high=high,positive_records=int((v>0).sum()),negative_records=int((v<0).sum()),zero_records=int((v==0).sum())))
                if source=='all':
                    for (sid,ds),vv in zip(subset.index,v):pairs.append(dict(replica_mode=mode,threshold=threshold,scope=scope,contrast=contrast,sample_id=sid,dataset=ds,difference=vv))
                if threshold==1 and scope=='all' and source=='all':draw_data[contrast]=d
            if threshold==1 and scope=='all' and source=='all':
                draw_data['record_indices']=indices
                draw_data['record_ids']=subset.index.get_level_values('sample_id').to_numpy(str)
    pd.DataFrame(rows).to_csv(HERE/f'{mode}_{kind}_summary.csv',index=False,float_format='%.17g')
    pd.DataFrame(paired_rows).to_csv(HERE/f'{mode}_paired_{kind}_difference_summary.csv',index=False,float_format='%.17g')
    pd.DataFrame(pairs).to_csv(HERE/f'{mode}_paired_{kind}_differences.csv',index=False,float_format='%.17g')
    np.savez_compressed(HERE/f'{mode}_{kind}_primary_bootstrap_draws.npz',**draw_data,seed=SEED,resamples=B)
    return pd.DataFrame(rows),pd.DataFrame(paired_rows)

def compute(plan):
    old_frames={mode:pd.read_csv(item['path'],dtype={'sample_id':str,'chain':str,'author_residue_id':str},float_precision='round_trip') for mode,item in plan['per_residue_sources'].items()}
    for frame in old_frames.values():assert len(frame)==33493 and frame.sample_id.nunique()==83
    identities=['sample_id','dataset','fold','residue_index','amino_acid','chain','author_residue_id','chi_count','full_hchi','intrachain_hchi','full_rmsf']
    pd.testing.assert_frame_equal(old_frames['replica1'][identities],old_frames['available_replicas'][identities],check_exact=True)
    tables={mode:[] for mode in old_frames}
    cache_outputs=[]
    for number,record in enumerate(plan['records'],start=1):
        base=old_frames['replica1'][old_frames['replica1'].sample_id==record['sample_id']].copy()
        full=np.load(record['full_prediction'],allow_pickle=False)
        intra=np.load(record['intra_prediction'],allow_pickle=False)
        idx=base.residue_index.to_numpy(int)
        np.testing.assert_array_equal(base.full_hchi.to_numpy(),full['mean'][idx,1].astype(float))
        np.testing.assert_array_equal(base.intrachain_hchi.to_numpy(),intra['mean'][idx,1].astype(float))
        np.testing.assert_array_equal(base.full_rmsf.to_numpy(),full['mean'][idx,0].astype(float))
        assert full['valid'][idx,:2].all() and intra['valid'][idx,:2].all()
        values=[]
        for c in record['caches']:
            metric,report=width_cache(c,base,record)
            values.append((c['replica'],metric));cache_outputs.append(report)
        primary=[v for replica,v in values if replica==1]
        assert len(primary)==1
        averaged=np.mean([v for _,v in values],axis=0)
        for mode,metric in [('replica1',primary[0]),('available_replicas',averaged)]:
            x=base[identities].copy();x['replica_mode']=mode;x['replicas']=1 if mode=='replica1' else len(values)
            for ti,threshold in enumerate(THRESHOLDS):x['width_Fstar_le_'+str(threshold).replace('.','p')+'_degrees']=metric[:,ti]
            count=x.groupby('amino_acid').residue_index.transform('size')
            x['sameAA_group_positions']=count;x['eligible_sameAA']=count>=4
            tables[mode].append(x)
        print('CACHE',number,83,record['sample_id'],len(base),len(values),flush=True)
    for mode,parts in tables.items():
        frame=pd.concat(parts,ignore_index=True)
        assert len(frame)==33493 and int(frame.eligible_sameAA.sum())==33142
        frame.to_csv(HERE/f'{mode}_per_residue.csv.gz',index=False,compression='gzip',float_format='%.17g')
        cr,gr,mr=[],[],[]
        for threshold in THRESHOLDS:
            col='width_Fstar_le_'+str(threshold).replace('.','p')+'_degrees'
            for _,record_frame in frame.groupby('sample_id',sort=True):
                temp=record_frame.copy();temp['width_degrees']=temp[col]
                c,g,m=record_stats(temp,threshold,mode);cr.extend(c);gr.extend(g);mr.extend(m)
        correlations=pd.DataFrame(cr);groups=pd.DataFrame(gr);masks=pd.DataFrame(mr)
        correlations.to_csv(HERE/f'{mode}_per_record_correlations.csv',index=False,float_format='%.17g')
        groups.to_csv(HERE/f'{mode}_high_low_per_AA.csv',index=False,float_format='%.17g')
        masks.to_csv(HERE/f'{mode}_eligibility_masks.csv',index=False)
        selected=groups[groups.rankable]
        high_record=selected.groupby(['sample_id','dataset','fold','replica_mode','threshold','scope','predictor']).agg(
            high_minus_low_MD=('high_minus_low_MD','mean'),low_MD_mean=('low_MD_mean','mean'),high_MD_mean=('high_MD_mean','mean'),
            AA_groups=('amino_acid','size'),positions=('positions','sum'),low_positions=('low_positions','sum'),high_positions=('high_positions','sum'),
            omitted_middle_positions=('omitted_middle_positions','sum')).reset_index()
        high_record.to_csv(HERE/f'{mode}_high_low_per_record.csv',index=False,float_format='%.17g')
        cs,cd=summary_with_pairs(correlations,'rho_type_centered',mode,'correlation')
        hs,hd=summary_with_pairs(high_record,'high_minus_low_MD',mode,'high_low')
        primary=lambda x:x[(x.threshold==1)&(x.scope=='all')&(x.dataset=='all')].to_dict(orient='records')
        write_json(HERE/f'{mode}_analysis_audit.json',dict(status='complete',finished_utc=utc(),plan_sha256=sha(HERE/'analysis_plan.json'),
            records=83,common_positions=len(frame),eligible_sameAA_positions=int(frame.eligible_sameAA.sum()),excluded_under4_sameAA_positions=int((~frame.eligible_sameAA).sum()),
            primary_correlations=primary(cs),primary_paired_correlations=primary(cd),primary_high_low=primary(hs),primary_paired_high_low=primary(hd),
            constant_prediction_groups=int((~groups.rankable).sum()),unreported_failed_records=0,
            exact_predictions_and_identities=True,scope='Scalar OOF scores associated with MD marginal low-F* bin-union width; no model-predicted F* curve or barrier/kinetics'))
        print(mode,'PRIMARY',json.dumps(dict(correlation=primary(cs),high_low=primary(hs),correlation_differences=primary(cd)),ensure_ascii=False),flush=True)
    outputs={}
    for path in HERE.rglob('*'):
        if path.is_file() and path.name not in ('completion_audit.json','run.log','run.stderr.log'):
            outputs[str(path.relative_to(HERE))]=dict(sha256=sha(path),bytes=path.stat().st_size)
    # Verify hashes again at the end so no input mutation is hidden.
    for path,digest in plan['source_hash_files'].items():assert sha(path)==digest
    write_json(HERE/'completion_audit.json',dict(status='complete',finished_utc=utc(),plan_sha256=sha(HERE/'analysis_plan.json'),code_sha256=sha(__file__),
        records=83,trajectories=203,source_hashes_checked_before_and_after=True,caches=cache_outputs,outputs=outputs))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=['freeze','compute'],default='compute');args=parser.parse_args()
    plan=freeze();print('FROZEN',sha(HERE/'analysis_plan.json'),plan['created_utc'],flush=True)
    if args.phase=='compute':compute(plan)

if __name__=='__main__':main()
