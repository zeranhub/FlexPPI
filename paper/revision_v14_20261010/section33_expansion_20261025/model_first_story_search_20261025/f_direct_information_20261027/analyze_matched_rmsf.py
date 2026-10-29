from pathlib import Path
from hashlib import sha256
import json
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parent;root.mkdir(parents=True,exist_ok=True)
ent=Path(r'E:\Project\FlexPPI\experiments\_sources\model_argument_v24_20261013\experiment_predictive_occupancy\replica1_per_residue.csv.gz')
wid=Path(r'E:\Project\FlexPPI\experiments\_sources\free_energy_model_v26_20261013\experiment\replica1_per_residue.csv.gz')
cohort=root.parent/'revision_v100_evidence_redesign'/'figure3'/'source_data'/'retained_absolute_correlations_all498.csv'
digest=lambda p:sha256(p.read_bytes()).hexdigest()
plan={'status':'fixed_before_matched_outcome_analysis','primary_matching_variable':'original_RMSF (frozen MD reference target, not claimed identical to replica1 ensemble)',
 'primary_tolerance_A':.10,'sensitivity_tolerances_A':[.05,.20],'additional_matching_variable_check':'full_rmsf',
 'pairing':'within same record and amino acid; sort matching RMSF then residue_index; greedy adjacent disjoint matching; if gap exceeds tolerance, skip lower index and try next; no residue reused',
 'orientation':'after matching, lower/higher frozen predicted full_hchi; no Hchi difference threshold; exclude exact Hchi ties only because direction undefined',
 'outcomes':['fine_entropy_mean_per_chi','width_Fstar_le_1p0_degrees'],'control_outcome':'original_RMSF',
 'weights':'equal pairs within amino-acid group, equal amino-acid groups within record, equal records in grand mean',
 'uncertainty':'4000 whole-record resamples stratified by source; seed20261027; identical indices for low/high/difference within each analysis',
 'no_success_filtering':True,'no_model_training_inference_or_new_MD':True,
 'inputs':{str(p):digest(p) for p in [ent,wid,cohort]}}
planpath=root/'analysis_plan_before_results.json'
if planpath.exists():assert json.loads(planpath.read_text(encoding='utf-8'))==plan
else:planpath.write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')

a=pd.read_csv(ent,float_precision='round_trip',low_memory=False)
b=pd.read_csv(wid,float_precision='round_trip',low_memory=False)
ids=set(pd.read_csv(cohort).sample_id)
assert set(a.sample_id)==set(b.sample_id)==ids and len(ids)==83
keys=['sample_id','residue_index']
assert not a.duplicated(keys).any() and not b.duplicated(keys).any()
z=a.merge(b[keys+['amino_acid','full_hchi','full_rmsf','width_Fstar_le_1p0_degrees']],on=keys,validate='one_to_one',suffixes=('','_width'))
for col in ['amino_acid','full_hchi','full_rmsf']:np.testing.assert_array_equal(z[col],z[col+'_width'])
z=z.drop(columns=[c for c in z if c.endswith('_width')])
metrics=['original_RMSF','full_rmsf','fine_entropy_mean_per_chi','width_Fstar_le_1p0_degrees']
summaries=[];audits=[]
for match_col in ['original_RMSF','full_rmsf']:
    for tolerance in [.05,.10,.20]:
        tag=f'{match_col}_tol{tolerance:.2f}'
        pairs=[];ties=0;eligible_pairs=0
        for (sid,aa),part in z.groupby(['sample_id','amino_acid'],sort=True):
            part=part.sort_values([match_col,'residue_index'],kind='stable').reset_index(drop=True)
            i=0
            while i+1<len(part):
                left,right=part.iloc[i],part.iloc[i+1]
                if right[match_col]-left[match_col]>tolerance:
                    i+=1;continue
                eligible_pairs+=1;i+=2
                if left.full_hchi==right.full_hchi:ties+=1;continue
                lo,hi=(left,right) if left.full_hchi<right.full_hchi else (right,left)
                row={'sample_id':sid,'dataset':lo.dataset,'amino_acid':aa,'low_residue_index':int(lo.residue_index),'high_residue_index':int(hi.residue_index),
                     'low_predicted_Hchi':lo.full_hchi,'high_predicted_Hchi':hi.full_hchi,'matching_variable':match_col,'tolerance_A':tolerance,
                     'absolute_matching_RMSF_gap_A':abs(hi[match_col]-lo[match_col])}
                for metric in metrics:row.update({metric+'_low':float(lo[metric]),metric+'_high':float(hi[metric]),metric+'_difference':float(hi[metric]-lo[metric])})
                pairs.append(row)
        pairs=pd.DataFrame(pairs);assert not pairs.empty
        for sid,p in pairs.groupby('sample_id'):
            used=np.r_[p.low_residue_index,p.high_residue_index];assert len(used)==len(set(used))
        assert pairs.absolute_matching_RMSF_gap_A.max()<=tolerance+1e-12
        paircols=[metric+suffix for metric in metrics for suffix in ['_low','_high','_difference']]
        aa=pairs.groupby(['sample_id','dataset','amino_acid'],as_index=False)[paircols].mean()
        rec=aa.groupby(['sample_id','dataset'],as_index=False)[paircols].mean().sort_values('sample_id').reset_index(drop=True)
        paircounts=pairs.groupby('sample_id').size();aacounts=aa.groupby('sample_id').size()
        rec['n_pairs']=rec.sample_id.map(paircounts);rec['n_AA_groups']=rec.sample_id.map(aacounts)
        pairs.to_csv(root/f'pairs_{tag}.csv.gz',index=False,compression='gzip')
        aa.to_csv(root/f'AAmeans_{tag}.csv',index=False)
        rec.to_csv(root/f'record_means_{tag}.csv',index=False)
        rng=np.random.default_rng(20261027);bootparts=[]
        for source in ['atlas','dynarepo','itrace']:
            pos=np.flatnonzero(rec.dataset.to_numpy()==source)
            if len(pos):bootparts.append(pos[rng.integers(0,len(pos),size=(4000,len(pos)))])
        indices=np.concatenate(bootparts,axis=1)
        for metric in metrics:
            for suffix in ['_low','_high','_difference']:
                values=rec[metric+suffix].to_numpy();lo,hi=np.quantile(values[indices].mean(axis=1),[.025,.975])
                summaries.append(dict(matching_variable=match_col,tolerance_A=tolerance,metric=metric,value_kind=suffix[1:],n_records=len(rec),n_pairs=len(pairs),
                    mean=float(values.mean()),bootstrap95_low=float(lo),bootstrap95_high=float(hi),positive_records=int((values>0).sum()) if suffix=='_difference' else None))
        audits.append(dict(matching_variable=match_col,tolerance_A=tolerance,raw_positions=len(z),pairs_without_Hchi_tie=len(pairs),exact_Hchi_ties=ties,
            included_positions=2*len(pairs),not_included_positions=len(z)-2*len(pairs),records=len(rec),missing_records=sorted(ids-set(rec.sample_id)),
            record_source_counts=rec.dataset.value_counts().to_dict(),min_pairs_per_record=int(rec.n_pairs.min()),max_pairs_per_record=int(rec.n_pairs.max()),
            AA_groups=len(aa),max_matching_gap=float(pairs.absolute_matching_RMSF_gap_A.max()),no_reused_residue=True))
summary=pd.DataFrame(summaries);summary.to_csv(root/'summary_all_prespecified_analyses.csv',index=False)
(root/'coverage_audit.json').write_text(json.dumps(audits,ensure_ascii=False,indent=2),encoding='utf-8')
show=summary[summary.value_kind.eq('difference')]
print(show[['matching_variable','tolerance_A','metric','n_records','n_pairs','mean','bootstrap95_low','bootstrap95_high','positive_records']].to_string(index=False))
print('PRIMARY COVERAGE '+json.dumps([a for a in audits if a['matching_variable']=='original_RMSF' and a['tolerance_A']==.1][0],ensure_ascii=False))
