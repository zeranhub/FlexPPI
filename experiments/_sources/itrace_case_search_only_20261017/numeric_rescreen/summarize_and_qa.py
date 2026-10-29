from pathlib import Path
import sys,json
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/'label_sensitivity/runtime'))
import numpy as np,pandas as pd

f=pd.read_csv(P/'all_prespecified_region_model_metrics.csv')
main=f[f.region_interpretable&f.scope.isin(['all_CDR','CDR3','whole_peptide','natural_MHC_groove','natural_pMHC_interface'])]
main=main[((main.metric.eq('rmsf'))&(main.chemistry.eq('including_substitutions')|main.scope.str.startswith('natural_')))|(main.metric.eq('hchi')&main.chemistry.eq('same_AA'))]
vals=['n','rho_delta','MAE_delta','zero_change_MAE','relative_MAE_gain','MD_abs_delta','model_abs_delta','descriptor_above_2SEM_n','descriptor_above_2SEM_direction']
v=main.pivot(index=['series','a','b','scope','model'],columns='metric',values=vals)
v.columns=['_'.join(c) for c in v.columns];v=v.reset_index()
v.to_csv(P/'all_models_complete_region_comparison.csv',index=False,encoding='utf-8-sig')
v[v.model.eq('full')].to_csv(P/'primary_case_region_comparison.csv',index=False,encoding='utf-8-sig')
v[v.model.eq('old_strict')].to_csv(P/'old_strict_case_region_comparison.csv',index=False,encoding='utf-8-sig')
g=v[(v.relative_MAE_gain_rmsf>0)&(v.relative_MAE_gain_hchi>0)]
g.to_csv(P/'dual_positive_gain_regions_all_models.csv',index=False,encoding='utf-8-sig')

# Independent numerical spot-checks from complete position-level rows.
positions=pd.read_csv(P/'all_aligned_position_changes.csv.gz',low_memory=False)
al=pd.read_csv(P/'complete_chain_alignment_audit.csv')
bound=pd.read_csv(P/'all_CDR_boundaries_and_sequences.csv')
short=al[al.role.eq('peptide')]
assert (short.n_aligned==short[['length_a','length_b']].min(axis=1)).all(), 'Peptide truncated by alignment'
checks=[]
for a,b,scope in [('3hg1','4jff','CDR3'),('7q99','7q9a','whole_peptide'),('8rym','8ryo','CDR3'),('2ak4','4jry','natural_pMHC_interface')]:
    p=positions[positions.a.eq(a)&positions.b.eq(b)]
    p=p.drop_duplicates(['role','index_a','index_b'])
    if scope=='CDR3':mask=p.role.isin(['tcr_alpha','tcr_beta'])&p.region_a.str.contains('CDR3',na=False)
    elif scope=='whole_peptide':mask=p.role.eq('peptide')
    else:mask=p.role.isin(['mhc','peptide'])&p.interface_union&p.same_aa
    for metric in ['rmsf','hchi']:
        select=mask&p[metric+'_valid']
        if metric=='hchi':select &= p.same_aa
        x=p.loc[select,metric+'_md_delta'].to_numpy()
        for model in ['full','old_strict','intrachain']:
            y=p.loc[select,f'{metric}_{model}_delta'].to_numpy()
            rr=main[main.a.eq(a)&main.b.eq(b)&main.scope.eq(scope)&main.metric.eq(metric)&main.model.eq(model)].iloc[0]
            assert len(x)==rr.n
            mae=np.mean(abs(x-y));zero=np.mean(abs(x))
            assert abs(mae-rr.MAE_delta)<1e-12
            assert abs(zero-rr.zero_change_MAE)<1e-12
            checks.append(dict(a=a,b=b,scope=scope,metric=metric,model=model,n=len(x),MAE_recalculated=mae,zero_MAE_recalculated=zero))
mel=bound[bound.a.eq('3hg1')&bound.b.eq('4jff')&bound.CDR.str.startswith('CDR3')]
assert mel.aligned_n.sum()==22
qa=dict(status='passed',n_structures=len(set(al.a)|set(al.b)),unique_pairs=len(al[['a','b']].drop_duplicates()),series_pair_memberships=len(al[['series','a','b']].drop_duplicates()),all_peptide_positions_retained=True,whole_MEL5_engineering_CDR3_RMSF_n=22,all_models_metrics_rows=len(f),numerical_spotchecks=checks,CDR_boundaries='All explicit sequences and indices retained in all_CDR_boundaries_and_sequences.csv')
(P/'independent_numeric_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
print('COMPLETED',qa['n_structures'],'structures',qa['unique_pairs'],'unique pairs',qa['series_pair_memberships'],'series pair memberships')
cols=['series','a','b','scope','model','n_rmsf','n_hchi','rho_delta_rmsf','rho_delta_hchi','relative_MAE_gain_rmsf','relative_MAE_gain_hchi','descriptor_above_2SEM_n_rmsf','descriptor_above_2SEM_n_hchi']
print('Primary Full both-positive gain fixed complete regions, no post hoc N/rho/amplitude filters:')
print(g[g.model.eq('full')][cols].to_string(index=False))
print('Added MEL5 series and preserved SB27, all five complete regions, both no-ITrace Full models:')
print(v[v.series.str.startswith(('MEL5','preserved'))&v.model.isin(['full','old_strict'])][cols].to_string(index=False))
