"""Frozen OOF-driven, same-amino-acid continuous angular-occupancy experiment.

Only this experiment_* directory is written. All v21/v22 caches, labels,
predictions, models and manuscript files are read-only inputs. No GPU or MD.
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
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse
import hashlib
import json
import time
import traceback
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
TASK=HERE.parents[1]
SENS=TASK/'label_sensitivity'
TRAJ=TASK/'energy_landscape_v21_20261012'/'trajectory_analysis'
OLD=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
BIN_ORIGINS=(0,1.25,2.5,3.75)
SEED=20261013
BOOTSTRAPS=10000
MIN_SAME_AA=4
PREDICTORS=('full_hchi','intrachain_hchi','full_rmsf')
SCOPES=('all','q1','q2','q3','q4')


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''):
            result.update(block)
    return result.hexdigest()


def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def write_json(path,value):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(clean(value),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def utc():return datetime.now(timezone.utc).isoformat()


def freeze():
    source_path=SENS/'sensitivity_plan.json'
    source=json.loads(source_path.read_text(encoding='utf-8'))
    audit_path=TRAJ/'completion_audit.json'
    audit=json.loads(audit_path.read_text(encoding='utf-8'))
    assert len(source['records'])==83 and audit['status']=='complete' and audit['trajectories']==203
    manifest_path=OLD/'neural_exports/evaluation/prediction_manifest.jsonl'
    manifest=[json.loads(x) for x in manifest_path.read_text(encoding='utf-8').splitlines() if x.strip()]
    intra={row['sample_id']:row for row in manifest if row['variant']=='intrachain'}
    frozen=[]
    for row in source['records']:
        assert sha(row['prediction'])==row['prediction_sha256']
        assert sha(row['old_label'])==row['old_label_sha256']
        im=intra[row['sample_id']]
        assert im['input_sha256']==row['structure_sha256'] and im['label_sha256']==row['old_label_sha256']
        assert im['fold']==row['fold'] and im['dataset']==row['dataset']
        path=OLD/'neural_exports/evaluation/intrachain'/(row['sample_id']+'.npz')
        assert path.is_file()
        caches=[]
        for item in row['trajectories']:
            name=row['sample_id']+'_rep'+str(item['replica_index'])+'.npz'
            cache=TRAJ/'cache'/name
            assert sha(cache)==audit['ca_and_chi_cache_sha256'][name]
            caches.append(dict(replica=item['replica_index'],path=str(cache),sha256=audit['ca_and_chi_cache_sha256'][name]))
        frozen.append(dict(sample_id=row['sample_id'],dataset=row['dataset'],fold=row['fold'],
            full_prediction=row['prediction'],full_prediction_sha256=row['prediction_sha256'],
            full_checkpoint_sha256=row['old_prediction_manifest']['checkpoint_sha256'],
            old_label=row['old_label'],old_label_sha256=row['old_label_sha256'],
            intra_prediction=str(path),intra_prediction_sha256=sha(path),intra_checkpoint_sha256=im['checkpoint_sha256'],
            input_structure_sha256=row['structure_sha256'],caches=caches))
    plan=dict(protocol='OOF_sameAA_fine_angular_occupancy_v1',created_utc=utc(),analysis_code_sha256=sha(__file__),
        records=frozen,source_plan_sha256=sha(source_path),trajectory_completion_sha256=sha(audit_path),
        prediction_manifest_sha256=sha(manifest_path),source_hash_files={str(source_path):sha(source_path),str(audit_path):sha(audit_path),str(manifest_path):sha(manifest_path)},
        sample_counts=dict(records=83,trajectories=203,atlas=30,dynarepo=23,itrace=30),
        primary_metric='Mean per-angle5degree-bin Shannon entropy: for each definedχ angle, mean entropy over4circular bin origins0/1.25/2.5/3.75degrees; arithmetic mean overtheqdefined angles; unitsnats/angle',
        auxiliary_metric='Per-angle inverseSimpson participation1/sum(p_bin²), mean overthe same4origins andqangles; auxiliary only, cannot replace thepredeclared primary metric',
        bin_width_degrees=5,bin_origins_degrees=list(BIN_ORIGINS),histogram_bins=72,all_angles_definition='Allchemically definedχangles,includingPro,withoutphysicalthree-rotamer interpretation; periodic180seamwrapped; no smoothing/pseudocount',
        primary_replicas='replica1onlyforeachsource',repeat_sensitivity='Arithmetic metricmean across allavailable independent cached replicas: ATLAS/ITrace3,DynaRepo1; per-anglemetric remains same',
        trajectory_grid=dict(start_ns=0,end_ns=100,step_ps=200,frames=501),
        common_positions='Anglecacheχvalid andCAvalid,Full/Intra original validRMSF/χlabels, bothmodelHchi andFullRMSF predictions finite; source identities/targets/masks asserted equal; exclude incompleteχangles with explicitcounts',
        type_control='Within each record keep amino-acidtypeswith>=4commonresidues; subtract thattypegroupmean separately fromMDmetric andeachprediction; Spearman overALLretained residualpositions in thatrecord; do not select favorable sites',
        group_scope='Overall,andχanglecountq1/2/3/4separately; AA-qmapping assertedunique; samepositions/types reusedforallpredictors',
        primary_predictors=list(PREDICTORS),amino_acid_mean_baseline='SameAAtype-only score isconstant withinAA; residual0 so rho undefined, not reported asrho0; high-low predicted difference0',
        summary_unit='One record/protein isunit; arithmeticmean ofdefinedper-recordrho; reportnumberofrecordsandpositions; no poolingresiduesasreplicates',
        bootstrap=dict(seed=SEED,resamples=BOOTSTRAPS,unit='record',method='Source-stratified resampling withoriginalsourcecounts,pairedindexresamples forcomparedpredictors,percentile2.5/97.5%; descriptiverecordresamplinginterval'),
        paired_difference='Withinpairedfinite recordset,FullHchi rho minusIntraHchi orFullRMSF rho; same pairedbootstrap',
        all_sameAA_pair_secondary='Allunordered pairs insideeligibleAA groups; symmetric+/-differences removearbitrarysequence-indexorientation; descriptiveSpearman ΔpredictionvsΔMD,recordmacro; no pair treatedasindependent',
        high_low_secondary='Within everyeligibleAAgroup, sorteachpredictor ascending(tiebycanonicalresidueindex), bottomfloor(n/2) andtopfloor(n/2); oddmiddle1omittedonlyfromsecondary; groupconstantprediction(<=1e-12span)unrankable andreported; actualMD highmean-lowmean, thenAAequalmean withinrecord; usesallrankablegroups, nooutcome-based selection',
        q_sensitivity='AllprimarypositionsincludingPro; qgroupresultsreporteligiblecounts; optionalsameAA-group-equalrho sensitivitypredeclared',
        freeze_boundary='Primarymetric fixedbeforeanycorrelationresult; originalcohort/caches/predictions/labels unchanged; nonewMD/training/inference; resultsreportedregardlessdirection')
    destination=HERE/'analysis_plan.json'
    if destination.exists():
        existing=json.loads(destination.read_text(encoding='utf-8'))
        assert existing['analysis_code_sha256']==sha(__file__)
        assert existing['records']==frozen and existing['source_plan_sha256']==sha(source_path)
        plan=existing
    else:write_json(destination,plan)
    return plan,sha(destination)


def occupancy(angles):
    assert angles.ndim==2 and angles.shape[0]==501 and 1<=angles.shape[1]<=4 and np.isfinite(angles).all()
    h,p,per_angle_h,per_angle_p=[],[],[],[]
    for qi in range(angles.shape[1]):
        entropies,participation=[],[]
        a=angles[:,qi].astype(float)
        for origin in BIN_ORIGINS:
            code=np.floor(((a+180-origin)%360)/5).astype(int)
            assert code.min()>=0 and code.max()<72
            counts=np.bincount(code,minlength=72).astype(float)
            assert counts.sum()==501
            probability=counts/501
            visited=probability[probability>0]
            entropies.append(float(-np.sum(visited*np.log(visited))))
            participation.append(float(1/np.sum(probability**2)))
        per_angle_h.append(float(np.mean(entropies)))
        per_angle_p.append(float(np.mean(participation)))
        h.extend(entropies)
        p.extend(participation)
    result=dict(fine_entropy_mean_per_chi=float(np.mean(h)),fine_entropy_sum_over_chi=float(np.sum(per_angle_h)),
        participation_mean_per_chi=float(np.mean(p)))
    for qi in range(4):
        result['chi'+str(qi+1)+'_fine_entropy']=per_angle_h[qi] if qi<len(per_angle_h) else np.nan
        result['chi'+str(qi+1)+'_participation']=per_angle_p[qi] if qi<len(per_angle_p) else np.nan
    assert 0<=result['fine_entropy_mean_per_chi']<=np.log(72)+1e-12
    assert 1<=result['participation_mean_per_chi']<=72+1e-12
    return result


def spearman(x,y):
    x,y=np.asarray(x,dtype=float),np.asarray(y,dtype=float)
    assert x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all()
    if len(x)<4 or np.ptp(x)<=1e-12 or np.ptp(y)<=1e-12:return np.nan
    rx=pd.Series(x).rank(method='average').to_numpy()
    ry=pd.Series(y).rank(method='average').to_numpy()
    return float(np.corrcoef(rx,ry)[0,1])


def record_statistics(frame,mode):
    identity=dict(sample_id=frame.sample_id.iloc[0],dataset=frame.dataset.iloc[0],fold=int(frame.fold.iloc[0]),replica_mode=mode)
    counts=frame.groupby('amino_acid').size()
    eligible=counts[counts>=MIN_SAME_AA].index
    all_rows=frame[frame.amino_acid.isin(eligible)].copy()
    group_rows=[]
    for aa,part in all_rows.groupby('amino_acid',sort=True):
        assert part.chi_count.nunique()==1
        for metric in ['fine_entropy_mean_per_chi','participation_mean_per_chi']:
            for predictor in PREDICTORS:
                x,y=part[predictor].to_numpy(),part[metric].to_numpy()
                group_rows.append(identity|dict(amino_acid=aa,chi_count=int(part.chi_count.iloc[0]),metric=metric,predictor=predictor,
                    positions=len(part),sameAA_rho=spearman(x,y)))
    groups=pd.DataFrame(group_rows)
    rows=[]
    highlow=[]
    pairs=[]
    mask=[]
    for scope in SCOPES:
        part=all_rows if scope=='all' else all_rows[all_rows.chi_count==int(scope[-1])]
        mask.append(identity|dict(scope=scope,common_positions=len(frame if scope=='all' else frame[frame.chi_count==int(scope[-1])]),
            eligible_positions=len(part),AA_types=part.amino_acid.nunique(),type_minimum=MIN_SAME_AA,
            excluded_under4_positions=len(frame if scope=='all' else frame[frame.chi_count==int(scope[-1])])-len(part)))
        if len(part)<4:continue
        for metric in ['fine_entropy_mean_per_chi','participation_mean_per_chi']:
            centered_y=part[metric]-part.groupby('amino_acid')[metric].transform('mean')
            for predictor in PREDICTORS:
                centered_x=part[predictor]-part.groupby('amino_acid')[predictor].transform('mean')
                g=groups[(groups.metric==metric)&(groups.predictor==predictor)]
                if scope!='all':g=g[g.chi_count==int(scope[-1])]
                rows.append(identity|dict(scope=scope,metric=metric,predictor=predictor,positions=len(part),AA_types=part.amino_acid.nunique(),
                    rho_type_centered=spearman(centered_x,centered_y),rho_AA_group_equal=float(g.sameAA_rho.mean()),
                    AA_types_defined_rho=int(g.sameAA_rho.notna().sum()),MD_centered_sample_SD=float(centered_y.std()),prediction_centered_sample_SD=float(centered_x.std())))
        for aa,aa_part in part.groupby('amino_acid',sort=True):
            for predictor in PREDICTORS:
                ordered=aa_part.sort_values([predictor,'residue_index'],kind='mergesort')
                half=len(ordered)//2
                constant=np.ptp(ordered[predictor].to_numpy())<=1e-12
                low,high=ordered.iloc[:half],ordered.iloc[-half:]
                for metric in ['fine_entropy_mean_per_chi','participation_mean_per_chi']:
                    highlow.append(identity|dict(scope=scope,amino_acid=aa,chi_count=int(ordered.chi_count.iloc[0]),predictor=predictor,metric=metric,
                        positions=len(ordered),low_positions=half,high_positions=half,omitted_middle_positions=len(ordered)-2*half,
                        rankable=not constant,excluded_reason='sameAApredictionconstant' if constant else '',
                        low_MD_mean=float(low[metric].mean()) if not constant else np.nan,
                        high_MD_mean=float(high[metric].mean()) if not constant else np.nan,
                        high_minus_low_MD=float(high[metric].mean()-low[metric].mean()) if not constant else np.nan,
                        low_prediction_mean=float(low[predictor].mean()) if not constant else np.nan,
                        high_prediction_mean=float(high[predictor].mean()) if not constant else np.nan,
                        low_indices=';'.join(str(i) for i in low.residue_index.tolist()) if not constant else '',
                        high_indices=';'.join(str(i) for i in high.residue_index.tolist()) if not constant else '',AA_mean_score_high_minus_low=0.0))
        # Pair contrasts are supplementary and computed only for the primary
        # entropy metric; both signs ensure sequence numbering cannot define
        # the direction of a pair. Protein remains the summary unit.
        for predictor in PREDICTORS:
            dx,dy=[],[]
            for aa,aa_part in part.groupby('amino_acid',sort=True):
                ids=np.triu_indices(len(aa_part),k=1)
                x=aa_part[predictor].to_numpy()
                y=aa_part.fine_entropy_mean_per_chi.to_numpy()
                dx.append(x[ids[1]]-x[ids[0]])
                dy.append(y[ids[1]]-y[ids[0]])
            x,y=np.concatenate(dx),np.concatenate(dy)
            pairs.append(identity|dict(scope=scope,predictor=predictor,unordered_sameAA_pairs=len(x),
                rho_symmetric_pair_contrasts=spearman(np.concatenate([x,-x]),np.concatenate([y,-y]))))
    return pd.DataFrame(rows),groups,pd.DataFrame(highlow),pd.DataFrame(pairs),pd.DataFrame(mask)


def process_record(arguments):
    record,mode,plan_sha=arguments
    start=time.monotonic()
    identity=dict(sample_id=record['sample_id'],dataset=record['dataset'],fold=record['fold'])
    try:
        for key in ['full_prediction','intra_prediction','old_label']:
            assert sha(record[key])==record[key+'_sha256']
        with np.load(record['full_prediction'],allow_pickle=False) as z:
            full={k:z[k] for k in z.files}
        with np.load(record['intra_prediction'],allow_pickle=False) as z:
            intra={k:z[k] for k in z.files}
        with np.load(record['old_label'],allow_pickle=False) as z:
            np.testing.assert_array_equal(full['target'],z['values'])
            np.testing.assert_array_equal(full['valid'],z['mask'])
        assert full['mean'].shape[1]==intra['mean'].shape[1]==2
        assert full['target'].shape[1]==3
        for name in ['target','valid','amino_acid','chain','author_residue_id']:
            np.testing.assert_array_equal(full[name],intra[name])
        selected=record['caches'][:1] if mode=='replica1' else record['caches']
        caches=[]
        for entry in selected:
            assert sha(entry['path'])==entry['sha256']
            with np.load(entry['path'],allow_pickle=False) as z:
                # Avoid loadingCAcoordinate arrays, which are not needed.
                caches.append({k:z[k] for k in ['chi_degrees','chi_count','chi_valid','ca_valid','times_ps','chain','author_residue_id','amino_acid','uniform_hchi']})
        common=full['valid'][:,:2].all(axis=1)&intra['valid'][:,:2].all(axis=1)
        common &= np.isfinite(full['mean']).all(axis=1)&np.isfinite(intra['mean']).all(axis=1)
        for cache in caches:
            for name in ['chain','author_residue_id','amino_acid']:
                np.testing.assert_array_equal(cache[name],full[name])
            np.testing.assert_allclose(cache['times_ps'],np.arange(501)*200,rtol=0,atol=.05)
            common &= cache['chi_valid']&cache['ca_valid']
            np.testing.assert_array_equal(cache['chi_count'],caches[0]['chi_count'])
        rows=[]
        original_H_maxerror=0.0
        for ri in np.flatnonzero(common):
            q=int(caches[0]['chi_count'][ri])
            results=[]
            for cache in caches:
                angle=cache['chi_degrees'][:,ri,:q]
                assert np.isfinite(angle).all()
                results.append(occupancy(angle))
                codes=np.clip(np.floor((angle+np.float32(180))/np.float32(120)).astype(int),0,2)
                joint=np.sum(codes*(3**np.arange(q)),axis=1)
                _,counts=np.unique(joint,return_counts=True)
                p=counts/501
                original_H=-np.sum(p*np.log(p))
                error=abs(original_H-cache['uniform_hchi'][ri])
                assert error<=1e-12
                original_H_maxerror=max(original_H_maxerror,error)
            values={key:float(np.mean([r[key] for r in results])) for key in results[0]}
            rows.append(identity|dict(replica_mode=mode,residue_index=int(ri),amino_acid=full['amino_acid'][ri],chain=full['chain'][ri],
                author_residue_id=full['author_residue_id'][ri],chi_count=q,replicas=len(caches),
                full_hchi=float(full['mean'][ri,1]),intrachain_hchi=float(intra['mean'][ri,1]),full_rmsf=float(full['mean'][ri,0]),
                original_joint_Hchi=float(full['target'][ri,1]),original_RMSF=float(full['target'][ri,0]),**values))
        frame=pd.DataFrame(rows)
        assert len(frame)>0
        assert frame.groupby('amino_acid').chi_count.nunique().le(1).all()
        correlations,groups,highlow,pairs,masks=record_statistics(frame,mode)
        part=HERE/'table_parts'/mode/record['sample_id']
        part.mkdir(parents=True,exist_ok=True)
        tables={}
        for name,table in [('per_residue',frame),('per_record_correlations',correlations),('per_AA_correlations',groups),
                           ('high_low_per_AA',highlow),('sameAA_pair_contrasts',pairs),('eligibility_masks',masks)]:
            path=part/(name+'.csv')
            table.to_csv(path,index=False,float_format='%.17g')
            tables[name]=dict(path=str(path),sha256=sha(path),rows=len(table))
        result=identity|dict(status='complete',replica_mode=mode,replicas=len(caches),common_positions=len(frame),
            original_joint_H_reproduction_max_abs=original_H_maxerror,targets_masks_identities_exact=True,
            plan_sha256=plan_sha,tables=tables,elapsed_seconds=time.monotonic()-start)
    except Exception as error:
        result=identity|dict(status='failed',replica_mode=mode,plan_sha256=plan_sha,error=repr(error),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-start)
    write_json(HERE/'record_audit'/mode/(record['sample_id']+'.json'),result)
    return result


def bootstrap(frame,value,groups,output_name):
    """Each estimate is arithmeticrecordmean; identical resamples comparemodels."""
    rng=np.random.default_rng(SEED)
    summaries=[]
    draws=[]
    for groupkey,part in frame.groupby(groups,dropna=False,sort=True):
        if not isinstance(groupkey,tuple):groupkey=(groupkey,)
        group=dict(zip(groups,groupkey))
        for source in ['all','atlas','dynarepo','itrace']:
            sample=part if source=='all' else part[part.dataset==source]
            sample=sample[sample[value].notna()].sort_values('sample_id')
            if len(sample)==0:continue
            assert not sample.sample_id.duplicated().any()
            index_parts=[]
            for dataset,subset in sample.groupby('dataset',sort=True):
                v=subset[value].to_numpy()
                indexes=rng.integers(0,len(v),size=(BOOTSTRAPS,len(v)))
                index_parts.append(v[indexes])
            values=np.concatenate(index_parts,axis=1).mean(axis=1)
            low,high=np.quantile(values,[.025,.975])
            summaries.append(group|dict(dataset=source,n_records=len(sample),mean=float(sample[value].mean()),record_sample_sd=float(sample[value].std()),
                median=float(sample[value].median()),min=float(sample[value].min()),max=float(sample[value].max()),bootstrap_95_low=float(low),bootstrap_95_high=float(high),
                positive_records=int((sample[value]>0).sum()),negative_records=int((sample[value]<0).sum()),zero_records=int((sample[value]==0).sum())))
            draws.append(pd.DataFrame(dict(bootstrap_index=np.arange(BOOTSTRAPS),value=values,**group,dataset=source)))
    summary=pd.DataFrame(summaries)
    summary.to_csv(HERE/(output_name+'.csv'),index=False,float_format='%.17g')
    pd.concat(draws,ignore_index=True).to_csv(HERE/(output_name+'.bootstrap_draws.csv.gz'),index=False,compression='gzip',float_format='%.17g')
    return summary


def analyze(plan,mode):
    reports=[json.loads((HERE/'record_audit'/mode/(r['sample_id']+'.json')).read_text(encoding='utf-8')) for r in plan['records']]
    assert all(r['status']=='complete' for r in reports) and len(reports)==83
    outputs={}
    assembled={}
    for name in ['per_residue','per_record_correlations','per_AA_correlations','high_low_per_AA','sameAA_pair_contrasts','eligibility_masks']:
        parts=[]
        for report in reports:
            item=report['tables'][name]
            assert sha(item['path'])==item['sha256']
            parts.append(pd.read_csv(item['path'],dtype={'sample_id':str,'chain':str,'author_residue_id':str}))
        table=pd.concat(parts,ignore_index=True)
        path=HERE/(mode+'_'+name+('.csv.gz' if name=='per_residue' else '.csv'))
        table.to_csv(path,index=False,compression='gzip' if path.suffix=='.gz' else None,float_format='%.17g')
        outputs[name]=dict(path=str(path),sha256=sha(path),rows=len(table))
        assembled[name]=table
    residues=assembled['per_residue']
    assert residues.sample_id.nunique()==83 and len(residues)==33493
    all_scopes=assembled['eligibility_masks']
    correlations=assembled['per_record_correlations']
    summary=bootstrap(correlations,'rho_type_centered',['replica_mode','scope','metric','predictor'],mode+'_correlation_summary')
    contrasts=[]
    for (scope,metric),part in correlations.groupby(['scope','metric'],sort=True):
        w=part.pivot(index=['sample_id','dataset'],columns='predictor',values='rho_type_centered')
        for baseline in ['intrachain_hchi','full_rmsf']:
            pairs=w[['full_hchi',baseline]].dropna()
            for (sample,dataset),row in pairs.iterrows():
                contrasts.append(dict(replica_mode=mode,scope=scope,metric=metric,sample_id=sample,dataset=dataset,
                    contrast='full_hchi_minus_'+baseline,difference=float(row.full_hchi-row[baseline])))
    contrast=pd.DataFrame(contrasts)
    contrast.to_csv(HERE/(mode+'_paired_correlation_differences.csv'),index=False,float_format='%.17g')
    diff_summary=bootstrap(contrast,'difference',['replica_mode','scope','metric','contrast'],mode+'_paired_correlation_difference_summary')
    highlow=assembled['high_low_per_AA']
    selected=highlow[highlow.rankable].copy()
    high_record=selected.groupby(['sample_id','dataset','fold','replica_mode','scope','metric','predictor']).agg(
        high_minus_low_MD=('high_minus_low_MD','mean'),AA_groups=('amino_acid','size'),
        positions=('positions','sum'),low_positions=('low_positions','sum'),high_positions=('high_positions','sum'),
        omitted_middle_positions=('omitted_middle_positions','sum')).reset_index()
    high_record['AA_mean_baseline_difference']=0.0
    high_record.to_csv(HERE/(mode+'_high_low_per_record.csv'),index=False,float_format='%.17g')
    high_summary=bootstrap(high_record,'high_minus_low_MD',['replica_mode','scope','metric','predictor'],mode+'_high_low_summary')
    pair_summary=bootstrap(assembled['sameAA_pair_contrasts'],'rho_symmetric_pair_contrasts',
        ['replica_mode','scope','predictor'],mode+'_pair_contrast_summary')
    group_summary=bootstrap(correlations,'rho_AA_group_equal',['replica_mode','scope','metric','predictor'],mode+'_AA_group_equal_correlation_summary')
    for path in HERE.glob(mode+'_*'):
        if path.is_file():outputs[path.name]=dict(path=str(path),sha256=sha(path),bytes=path.stat().st_size)
    primary=summary[(summary.scope=='all')&(summary.metric=='fine_entropy_mean_per_chi')&(summary.dataset=='all')]
    primary_diff=diff_summary[(diff_summary.scope=='all')&(diff_summary.metric=='fine_entropy_mean_per_chi')&(diff_summary.dataset=='all')]
    primary_high=high_summary[(high_summary.scope=='all')&(high_summary.metric=='fine_entropy_mean_per_chi')&(high_summary.dataset=='all')]
    write_json(HERE/(mode+'_analysis_audit.json'),dict(status='complete',replica_mode=mode,plan_sha256=sha(HERE/'analysis_plan.json'),
        records=83,common_positions=len(residues),eligible_sameAA_positions=int(all_scopes[all_scopes.scope=='all'].eligible_positions.sum()),
        excluded_under4_sameAA_positions=int(all_scopes[all_scopes.scope=='all'].excluded_under4_positions.sum()),
        all_original_targets_masks_identities_exact=True,source_hashes_checked=True,primary_metric_unchanged=True,
        primary_correlations=primary.to_dict(orient='records'),primary_paired_differences=primary_diff.to_dict(orient='records'),
        primary_high_low=primary_high.to_dict(orient='records'),outputs=outputs,finished_utc=utc(),
        inference_scope='Finite-record resamplingCI; outcomesareindependent finerangular descriptors,notmodel-predictedFcurve; sameAAcontrolremovesconstanttypevalue; n-recordsnotpair/residuecount'))
    print(mode,'PRIMARYCORRELATIONS',primary.to_string(index=False),flush=True)
    print(mode,'PRIMARYPAIREDDIFF',primary_diff.to_string(index=False),flush=True)
    print(mode,'PRIMARYHIGHLOW',primary_high.to_string(index=False),flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['freeze','replica1','available_replicas','all'],default='all')
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    HERE.mkdir(parents=True,exist_ok=True)
    plan,plan_sha=freeze()
    print(json.dumps(dict(plan_sha256=plan_sha,code_sha256=sha(__file__),records=83,phase=args.phase)),flush=True)
    modes=['replica1','available_replicas'] if args.phase=='all' else ([] if args.phase=='freeze' else [args.phase])
    for mode in modes:
        jobs=[(row,mode,plan_sha) for row in plan['records']]
        results=[]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(process_record,x) for x in jobs]
            for index,future in enumerate(as_completed(futures),start=1):
                result=future.result()
                results.append(result)
                print(mode,index,83,result['sample_id'],result['status'],round(result['elapsed_seconds'],2),flush=True)
                if result['status']!='complete':print(result['traceback'],flush=True)
        write_json(HERE/(mode+'_compute_audit.json'),dict(status='complete' if all(x['status']=='complete' for x in results) else 'partial',
            plan_sha256=plan_sha,records=len(results),failures=[r for r in results if r['status']!='complete'],finished_utc=utc()))
        assert all(r['status']=='complete' for r in results),'No failures silentlyexcluded'
        analyze(plan,mode)
    if len(modes)==2:
        write_json(HERE/'completion_audit.json',dict(status='complete',plan_sha256=plan_sha,code_sha256=sha(__file__),finished_utc=utc(),
            replica1_analysis_sha256=sha(HERE/'replica1_analysis_audit.json'),available_replicas_analysis_sha256=sha(HERE/'available_replicas_analysis_audit.json')))


if __name__=='__main__':
    main()
