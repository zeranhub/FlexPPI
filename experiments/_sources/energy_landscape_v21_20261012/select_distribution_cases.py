from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'label_sensitivity'
OUT=ROOT/'case_analysis'
PLAN=OLD/'sensitivity_plan.json'
METRIC=OLD/'prediction_effect_per_record.csv'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
plan=json.loads(PLAN.read_text(encoding='utf8'))
records={r['sample_id']:r for r in plan['records']}
score=pd.read_csv(METRIC)
score=score[score.variant=='uniform_rep1'].pivot(index=['sample_id','dataset'],columns='metric',values='old_mae').reset_index()
choices=[]
for source in ['atlas','dynarepo','itrace']:
    d=score[score.dataset==source].copy()
    # Percentile ranks are computed within source because the targets have
    # different physical units and the sources have different error scales.
    d['mean_error_rank']=(d.rmsf.rank(method='average',pct=True)+d.chi_entropy.rank(method='average',pct=True))/2
    d['median_rank_distance']=abs(d.mean_error_rank-.5)
    d=d.sort_values(['median_rank_distance','sample_id'])
    row=d.iloc[0]
    choices.append({'sample_id':row.sample_id,'dataset':source,'mean_error_rank':float(row.mean_error_rank),
                    'rmsf_mae_A':float(row.rmsf),'hchi_mae_nats':float(row.chi_entropy)})
protocol={'status':'selection_frozen_before_selected_angle_distributions',
          'selection':'One record per source, average of within-source original OOF RMSF/Hchi MAE percentile ranks closest to 0.5; tie by sample_id',
          'candidate_records':83,'input_sha256':{str(PLAN):sha(PLAN),str(METRIC):sha(METRIC)},
          'choices':choices,'additional_existing_case':'2QHO ubiquitin Q41/Q62 from v20',
          'residue_rule':'Among nonPro residues with both valid targets and chi1, group by AA with at least4; within each group use original reference RMSF/Hchi percentile-rank difference. Choose same-AA pair maximizing its largest minus smallest difference, tie by AA then input index; no prediction-score selection',
          'angular_window':'0–100ns,200ps,501 actual frames,first trajectory; all three repeats shown when available',
          'F_star':'5degree equal-width bins; -ln(probability/max probability); zero count missing; no smoothing'}
(OUT/'representative_cases_selection.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding='utf8')
profiles=[];selected_rows=[];time_rows=[];hist_rows=[];sources=[]
for choice in choices:
    sample=choice['sample_id'];row=records[sample]
    prediction=Path(row['prediction'])
    assert sha(prediction)==row['prediction_sha256']
    with np.load(prediction,allow_pickle=False) as z:
        y=z['target'][:,:2];p=z['mean'][:,:2];valid=z['valid'][:,:2].astype(bool)
    file=ROOT/'trajectory_analysis/cache'/f'{sample}_rep1.npz'
    with np.load(file,allow_pickle=False) as z:
        aa=z['amino_acid'].astype(str);chain=z['chain'].astype(str);ids=z['author_residue_id'].astype(str)
        q=z['chi_degrees'];chi_ok=z['chi_atom_available'][:,0].astype(bool)
    assert y.shape==(len(aa),2)
    common=valid.all(axis=1)&chi_ok&(aa!='P')&np.isfinite(y).all(axis=1)
    candidates=[]
    for typ in sorted(set(aa[common])):
        index=np.flatnonzero(common&(aa==typ))
        if len(index)<4:continue
        rrank=pd.Series(y[index,0]).rank(method='average',pct=True).to_numpy()
        hrank=pd.Series(y[index,1]).rank(method='average',pct=True).to_numpy()
        gap=hrank-rrank
        low=np.lexsort((index,gap))[0];high=np.lexsort((index,-gap))[0]
        if low==high:continue
        candidates.append((float(gap[high]-gap[low]),typ,int(index[high]),int(index[low]),float(gap[high]),float(gap[low])))
    assert candidates,sample
    gap,typ,hi,lo,gh,gl=sorted(candidates,key=lambda x:(-x[0],x[1],x[2],x[3]))[0]
    selected=[hi,lo]
    for i in range(len(aa)):
        profiles.append({'sample_id':sample,'dataset':row['dataset'],'input_index':i,'chain':chain[i],
                         'author_residue_id':ids[i],'amino_acid':aa[i],
                         'rmsf_valid':int(valid[i,0]),'hchi_valid':int(valid[i,1]),
                         'rmsf_reference':float(y[i,0]) if valid[i,0] else np.nan,
                         'rmsf_prediction':float(p[i,0]) if valid[i,0] else np.nan,
                         'hchi_reference':float(y[i,1]) if valid[i,1] else np.nan,
                         'hchi_prediction':float(p[i,1]) if valid[i,1] else np.nan})
    for i in selected:
        selected_rows.append({'sample_id':sample,'dataset':row['dataset'],'input_index':i,'chain':chain[i],
                              'author_residue_id':ids[i],'amino_acid':aa[i],
                              'reference_H_minus_R_rank':gh if i==hi else gl,
                              'rmsf_reference':float(y[i,0]),'rmsf_prediction':float(p[i,0]),
                              'hchi_reference':float(y[i,1]),'hchi_prediction':float(p[i,1])})
    for rep in range(1,len(row['trajectories'])+1):
        path=ROOT/'trajectory_analysis/cache'/f'{sample}_rep{rep}.npz'
        sources.append({'path':str(path),'sha256':sha(path)})
        with np.load(path,allow_pickle=False) as z:
            angles=z['chi_degrees'];times=z['times_ps']/1000
        for i in selected:
            a=angles[:,i,0]
            assert np.isfinite(a).all()
            for t,val in zip(times,a):time_rows.append({'sample_id':sample,'replica':rep,'input_index':i,
                                                       'residue':aa[i]+ids[i],'time_ns':float(t),'chi1_degrees':float(val)})
            counts,edges=np.histogram(a,bins=np.arange(-180,180.001,5))
            prob=counts/counts.sum();f=np.full(len(prob),np.nan)
            f[prob>0]=-np.log(prob[prob>0]/prob.max())
            for j in range(len(prob)):hist_rows.append({'sample_id':sample,'replica':rep,'input_index':i,
                    'residue':aa[i]+ids[i],'chi1_degrees':float((edges[j]+edges[j+1])/2),'count':int(counts[j]),
                    'probability':float(prob[j]),'F_star':float(f[j])})
for name,rows in [('representative_full_profiles',profiles),('representative_selected_residues',selected_rows),
                  ('representative_chi1_timeseries',time_rows),('representative_chi1_Fstar',hist_rows)]:
    pd.DataFrame(rows).to_csv(OUT/(name+'.csv'),index=False)
audit={'status':'complete','proteins':3,'additional_ubiquitin_protein':1,'selected_residues':6,
       'cache_sources':sources,'selection_sha256':sha(OUT/'representative_cases_selection.json'),
       'code_sha256':sha(__file__),'no_frames_subsampled':True,'source_targets_and_predictions_unchanged':True,
       'choices':choices,'residues':selected_rows}
(OUT/'representative_cases_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'choices':choices,'residues':selected_rows},ensure_ascii=False,indent=2))
