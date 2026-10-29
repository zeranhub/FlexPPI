"""Evaluate all frozen candidates without residue pseudo-replication."""
from pathlib import Path
import json, itertools, hashlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
OUT=Path(__file__).resolve().parent
ANALYSIS=OUT.parent
KEY=['sample_id','residue_index']
METHODS=['full','intra','static47_sasa','aa_sasa_contacts','aa_mean']


def rho(a,b):
    return float(spearmanr(a,b).statistic) if len(a)>=3 and np.ptp(a)>0 and np.ptp(b)>0 else np.nan


def merged_table():
    f=pd.read_parquet(OUT/'candidate_with_baselines56.parquet')
    original=pd.read_parquet(ANALYSIS/'shared/oof_candidates56.parquet')
    for source,method in [('full','full'),('intrachain','intra')]:
        x=original[original.variant.eq(source)].copy()
        keep=KEY+['chain','author_residue_id','amino_acid','fold','target_rmsf','target_chi_entropy','valid_rmsf','valid_chi_entropy','pred_rmsf','pred_chi_entropy']
        x=x[keep].rename(columns={c:'check_'+c for c in keep if c not in KEY})
        f=f.merge(x,on=KEY,how='left',validate='one_to_one')
        assert f.check_chain.eq(f.chain).all()
        assert f.check_author_residue_id.astype(str).eq(f.author_residue_id.astype(str)).all()
        assert f.check_amino_acid.eq(f.amino_acid).all()
        assert f.check_fold.eq(f.fold).all()
        for metric,rawmetric in [('rmsf','rmsf'),('hchi','chi_entropy')]:
            valid=f['valid_'+metric].to_numpy()
            assert np.array_equal(valid,f['check_valid_'+rawmetric])
            assert np.allclose(f.loc[valid,'md_'+metric],f.loc[valid,'check_target_'+rawmetric],rtol=0,atol=1e-5)
            f['pred_'+method+'_'+metric]=f['check_pred_'+rawmetric].to_numpy()
        f=f.drop(columns=[c for c in f if c.startswith('check_')])
    for metric in ['rmsf','hchi']:
        for method in METHODS:
            assert np.isfinite(f.loc[f['valid_'+metric],'pred_'+method+'_'+metric]).all()
    f.to_parquet(OUT/'all56_features_and_family_predictions.parquet',index=False)
    return f


def record_metrics(f):
    records=[]
    for (sid,role),part in f.groupby(['sample_id','role'],sort=True):
        common={k:part[k].iloc[0] for k in ['pdb_id','component_id','fold','chain']}
        common.update(sample_id=sid,role=role)
        sets={'all':np.ones(len(part),bool),'interface_ca8':part.interface_ca8.to_numpy(),'noninterface_ca8':~part.interface_ca8.to_numpy(),'complete_sidechains':part.sidechain_complete.to_numpy(),'buried_sidechain_le5_A2':part.sidechain_sasa_complex_A2.le(5).to_numpy(),'exposed_sidechain_ge40_A2':part.sidechain_sasa_complex_A2.ge(40).to_numpy()}
        for subset,mask in sets.items():
            for metric in ['rmsf','hchi']:
                keep=mask&part['valid_'+metric].to_numpy()
                y=part.loc[keep,'md_'+metric].to_numpy()
                for method in METHODS:
                    pred=part.loc[keep,'pred_'+method+'_'+metric].to_numpy()
                    error=pred-y
                    records.append(dict(common,subset=subset,metric=metric,method=method,n_valid=len(y),mae=float(np.mean(np.abs(error))) if len(y) else np.nan,rmse=float(np.sqrt(np.mean(error**2))) if len(y) else np.nan,bias=float(np.mean(error)) if len(y) else np.nan,rho=rho(y,pred)))
    result=pd.DataFrame(records)
    result.to_csv(OUT/'metrics_per_record_all56.csv',index=False)
    return result


def cluster_interval(pdb,value_col='gain',draws=10000):
    # Resample family clusters; retain all PDBs per selected family. Estimand stays PDB macro.
    family=pdb.groupby('component_id')[value_col].agg(['sum','count','mean'])
    rng=np.random.default_rng(20261026)
    ids=rng.integers(0,len(family),size=(draws,len(family)))
    sums=family['sum'].to_numpy()[ids].sum(axis=1)
    counts=family['count'].to_numpy()[ids].sum(axis=1)
    low,high=np.quantile(sums/counts,[.025,.975])
    flo,fhi=np.quantile(family['mean'].to_numpy()[ids].mean(axis=1),[.025,.975])
    return dict(pdb_macro_gain=float(pdb[value_col].mean()),pdb_gain_ci_lo=float(low),pdb_gain_ci_hi=float(high),family_macro_gain=float(family['mean'].mean()),family_gain_ci_lo=float(flo),family_gain_ci_hi=float(fhi),n_pdb=len(pdb),n_family=len(family),n_pdb_positive=int((pdb[value_col]>0).sum()),n_pdb_negative=int((pdb[value_col]<0).sum()),n_pdb_ties=int((pdb[value_col]==0).sum()),bootstrap_draws=draws)


def summarize(record):
    levels=[]
    group=['role','subset','metric','method']
    for keys,part in record.groupby(group):
        pdb=part[part.n_valid.gt(0)].groupby(['pdb_id','component_id'])[['mae','rmse','bias','rho']].mean().reset_index()
        if len(pdb)==0: continue
        family=pdb.groupby('component_id')[['mae','rmse','bias','rho']].mean()
        levels.append(dict(zip(group,keys),n_records=int(part.n_valid.gt(0).sum()),n_pdb=len(pdb),n_family=len(family),pdb_macro_mae=pdb.mae.mean(),pdb_macro_rho=pdb.rho.mean(),family_macro_mae=family.mae.mean(),family_macro_rho=family.rho.mean()))
    summary=pd.DataFrame(levels)
    summary.to_csv(OUT/'summary_pdb_family_macro.csv',index=False)
    intervals=[]
    for (role,subset,metric),part in record.groupby(['role','subset','metric']):
        w=part.pivot(index=['sample_id','pdb_id','component_id'],columns='method',values='mae')
        for comparator in [m for m in METHODS if m!='full']:
            p=w[['full',comparator]].dropna().copy()
            if len(p)==0: continue
            p['gain']=p[comparator]-p['full']
            pdb=p.reset_index().groupby(['pdb_id','component_id'])[['gain','full',comparator]].mean().reset_index()
            intervals.append(dict(role=role,subset=subset,metric=metric,comparison='full_vs_'+comparator,positive='comparator_MAE_minus_Full_MAE',n_records=len(p),full_pdb_macro_mae=float(pdb['full'].mean()),comparator_pdb_macro_mae=float(pdb[comparator].mean()),**cluster_interval(pdb)))
    ci=pd.DataFrame(intervals)
    ci.to_csv(OUT/'paired_full_vs_baselines_family_cluster_ci.csv',index=False)
    return summary,ci


def matched_pairs(f):
    pair_rows=[]; coverage=[]
    for (sid,role),part in f.groupby(['sample_id','role'],sort=True):
        common={k:part[k].iloc[0] for k in ['pdb_id','component_id']};common.update(sample_id=sid,role=role)
        p=part[part.valid_hchi & part.valid_rmsf]
        for tolerance in [.10,.25,.50]:
            count=0
            for aa,a in p.groupby('amino_acid'):
                for i,j in itertools.combinations(a.index,2):
                    v=f.loc[i];w=f.loc[j]
                    if abs(v.md_rmsf-w.md_rmsf)>tolerance: continue
                    count+=1
                    rawdiff=w.md_hchi-v.md_hchi
                    row=dict(common,amino_acid=aa,rmsf_tolerance_A=tolerance,residue_i=int(v.residue_index),residue_j=int(w.residue_index),target_hchi_difference=float(rawdiff),target_rmsf_difference=float(w.md_rmsf-v.md_rmsf),sasa_difference_A2=float(w.sidechain_sasa_complex_A2-v.sidechain_sasa_complex_A2))
                    for method in METHODS:
                        d=w['pred_'+method+'_hchi']-v['pred_'+method+'_hchi']
                        row[method+'_difference']=float(d);row[method+'_absolute_error']=float(abs(d-rawdiff))
                    pair_rows.append(row)
            coverage.append(dict(common,rmsf_tolerance_A=tolerance,n_pairs=count))
    pairs=pd.DataFrame(pair_rows);cover=pd.DataFrame(coverage)
    pairs.to_parquet(OUT/'sameAA_backbone_matched_all_pairs.parquet',index=False)
    cover.to_csv(OUT/'sameAA_backbone_matched_coverage_all56.csv',index=False)
    metrics=[]
    for keys,p in pairs.groupby(['sample_id','pdb_id','component_id','role','rmsf_tolerance_A']):
        common=dict(zip(['sample_id','pdb_id','component_id','role','rmsf_tolerance_A'],keys))
        for method in METHODS:
            metrics.append(dict(common,method=method,n_pairs=len(p),pair_difference_mae=p[method+'_absolute_error'].mean(),pair_difference_rho=rho(p.target_hchi_difference,p[method+'_difference'])))
    metrics=pd.DataFrame(metrics)
    metrics.to_csv(OUT/'sameAA_backbone_matched_record_metrics.csv',index=False)
    intervals=[]
    for (role,tol),p in metrics.groupby(['role','rmsf_tolerance_A']):
        w=p.pivot(index=['sample_id','pdb_id','component_id'],columns='method',values='pair_difference_mae')
        for comparator in [m for m in METHODS if m!='full']:
            a=w[['full',comparator]].dropna().copy();a['gain']=a[comparator]-a['full']
            pdb=a.reset_index().groupby(['pdb_id','component_id'])[['gain','full',comparator]].mean().reset_index()
            intervals.append(dict(role=role,rmsf_tolerance_A=tol,comparison='full_vs_'+comparator,full_pair_mae=pdb.full.mean(),comparator_pair_mae=pdb[comparator].mean(),n_records=len(a),**cluster_interval(pdb)))
    intervals=pd.DataFrame(intervals)
    intervals.to_csv(OUT/'sameAA_backbone_matched_paired_ci.csv',index=False)
    return pairs,cover,metrics,intervals


def stratified_aa(f):
    rows=[]
    for keys,p in f[f.valid_hchi].groupby(['sample_id','pdb_id','component_id','role','amino_acid']):
        common=dict(zip(['sample_id','pdb_id','component_id','role','amino_acid'],keys))
        for subset,mask in [('all',np.ones(len(p),bool)),('buried_le5',p.sidechain_sasa_complex_A2.le(5)),('exposed_ge40',p.sidechain_sasa_complex_A2.ge(40))]:
            a=p.loc[mask]
            for method in METHODS:
                rows.append(dict(common,subset=subset,method=method,n_valid=len(a),mae=float(abs(a['pred_'+method+'_hchi']-a.md_hchi).mean()) if len(a) else np.nan,target_mean=float(a.md_hchi.mean()) if len(a) else np.nan,pred_mean=float(a['pred_'+method+'_hchi'].mean()) if len(a) else np.nan))
    pd.DataFrame(rows).to_csv(OUT/'sameAA_environment_record_metrics.csv',index=False)


def main():
    f=merged_table();r=record_metrics(f);s,ci=summarize(r)
    pairs,cover,m,pci=matched_pairs(f);stratified_aa(f)
    complete=dict(records=f.sample_id.nunique(),pdb=f.pdb_id.nunique(),family=f.component_id.nunique(),rows=len(f),merged_key_1_to_1=True,old_fold_predictions_used=False,feature_and_family_oof_targets_identical_with_atol_1e_5=True,all_methods_valid_finite=True,all56_included=True,pair_count=len(pairs),note='same-AA pairs are correlated within a record; aggregate record then PDB/family; no per-residue p-values')
    (OUT/'evaluation_audit.json').write_text(json.dumps(complete,indent=2,ensure_ascii=False),encoding='utf-8')
    print('MAIN ALL')
    print(ci[ci.subset.eq('all')].to_string(index=False))
    print('SAME AA .25')
    print(pci[pci.rmsf_tolerance_A.eq(.25)].to_string(index=False))
if __name__=='__main__':main()
