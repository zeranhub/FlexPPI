from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'label_sensitivity/runtime'))
import numpy as np
import pandas as pd
z=pd.read_csv(ROOT/'numeric_rescreen/all_aligned_position_changes.csv.gz',dtype={'id_a':str,'id_b':str},low_memory=False)
pairs=[('7ow5','7ow6'),('7ow6','7pb2'),('7q99','7q9a'),('4qok','7q99'),('4qok','7q9a'),('8rym','8ryn'),('8rym','8ryo'),('3qdg','3qdj'),('3qdm','3qeq')]
OUT=ROOT/'promising_cases_replica_validation';OUT.mkdir(exist_ok=True)
allstats=[];files=[]
def rho(a,b):
    return float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1]) if len(a)>=4 else None
for a,b in pairs:
    f=z[z.a.eq(a)&z.b.eq(b)].drop_duplicates(['role','id_a','id_b']).copy()
    assert len(f), (a,b)
    for state,pdb in [('a',a),('b',b)]:
        with np.load(ROOT/'replica_recheck'/pdb/f'{pdb}_all_residue_replica_labels.npz') as q:
            indices=f['index_'+state].to_numpy(dtype=int)
            np.testing.assert_array_equal(f['chain_'+state].to_numpy(dtype=str),q['chain'][indices])
            np.testing.assert_array_equal(f['id_'+state].to_numpy(dtype=str),q['author_residue_id'][indices])
            np.testing.assert_array_equal(f['aa_'+state].to_numpy(dtype=str),q['amino_acid'][indices])
            for j,m in enumerate(['rmsf','hchi']):
                valid=f[m+'_valid'].to_numpy(dtype=bool)
                np.testing.assert_allclose(f.loc[valid,m+'_md_'+state],q['mean'][indices[valid],j],atol=3e-6,rtol=0)
                f[m+'_md_'+state]=q['mean'][indices,j]
                for rep in range(3):f[f'{m}_rep{rep+1}_{state}']=q['values'][rep,indices,j]
                f[m+'_variance_'+state]=q['sample_variance'][indices,j]
    for m in ['rmsf','hchi']:
        f[m+'_md_delta']=f[m+'_md_b']-f[m+'_md_a']
        f[m+'_delta_sem']=np.sqrt((f[m+'_variance_a']+f[m+'_variance_b'])/3)
        av=f[[f'{m}_rep{i}_a' for i in [1,2,3]]].to_numpy()
        bv=f[[f'{m}_rep{i}_b' for i in [1,2,3]]].to_numpy()
        f[m+'_replica_ranges_nonoverlap']=(np.max(av,axis=1)<np.min(bv,axis=1))|(np.max(bv,axis=1)<np.min(av,axis=1))
    f.to_csv(OUT/f'{a}_{b}_all_positions_three_replicas.csv',index=False)
    files.append(OUT/f'{a}_{b}_all_positions_three_replicas.csv')
    scopes={'whole_peptide_all':f.role.eq('peptide'),'whole_peptide_unchanged_AA':f.role.eq('peptide')&f.same_aa,
      'all_CDR':f.role.str.startswith('tcr')&f.region_a.str.startswith('CDR'),
      'CDR3':f.role.str.startswith('tcr')&f.region_a.str.startswith('CDR3'),
      'natural_pMHC_interface':f.role.isin(['mhc','peptide'])&f.interface_union&f.same_aa}
    for scope,mask in scopes.items():
        for m in ['rmsf','hchi']:
            ff=f[mask&f[m+'_valid']&(f.same_aa if m=='hchi' else True)]
            if not len(ff):continue
            y=ff[m+'_md_delta'].to_numpy()
            for model in ['full','old_strict','intrachain']:
                p=ff[f'{m}_{model}_delta'].to_numpy()
                large=abs(y)>2*ff[m+'_delta_sem'].to_numpy()
                allstats.append(dict(a=a,b=b,scope=scope,metric=m,model=model,n=len(ff),
                    rho=rho(y,p),MAE=float(np.mean(abs(y-p))),zero_MAE=float(np.mean(abs(y))),
                    MAE_improvement=float(1-np.mean(abs(y-p))/np.mean(abs(y))),
                    MD_mean_delta=float(np.mean(y)),model_mean_delta=float(np.mean(p)),
                    direction_correct_n=int((y*p>0).sum()),n_abs_delta_above_2SEM=int(large.sum()),
                    direction_correct_among_above_2SEM_n=int((y[large]*p[large]>0).sum()),
                    nonoverlapping_replica_ranges_n=int(ff[m+'_replica_ranges_nonoverlap'].sum())))
stats=pd.DataFrame(allstats);stats.to_csv(OUT/'recomputed_replica_complete_region_summary.csv',index=False)
audit={'status':'passed','cases':pairs,'recomputed_MD_means_and_variances_match_ITrace_cache':True,
       'regions':'Complete predeclared regions; all Hchi comparisons use unchanged AA. Full backbone retains substitutions, and unchanged-AA peptide is separately reported.',
       'source':'Each state has3 existing independent200ns MD runs, all frames retained; replicas not residues or frames are the repeat unit.',
       'statistics':'Descriptive rho and MAE only. abs(delta)>2SEM and disjoint replica ranges are descriptive checks, not formal significance tests.',
       'source_SHA256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(OUT/'replica_comparison_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
selection=((stats.a.eq('7ow5')&stats.b.eq('7ow6')&stats.scope.eq('whole_peptide_unchanged_AA')&stats.metric.eq('rmsf'))|
           (stats.a.eq('7q99')&stats.b.eq('7q9a')&stats.scope.eq('CDR3')&stats.metric.eq('rmsf'))|
           (stats.a.eq('8rym')&stats.b.eq('8ryn')&stats.scope.eq('natural_pMHC_interface')&stats.metric.eq('hchi'))|
           (stats.a.eq('3qdm')&stats.b.eq('3qeq')&stats.scope.eq('natural_pMHC_interface')&stats.metric.eq('hchi')))
print(stats[selection].to_string(index=False))
