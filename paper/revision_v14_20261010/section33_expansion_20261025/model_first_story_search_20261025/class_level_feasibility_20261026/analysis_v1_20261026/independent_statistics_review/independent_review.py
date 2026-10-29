"""Independent read-only audit. Writes only next to this script; never fits models."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OUT=Path(__file__).resolve().parent
R=OUT.parent
ROOT=Path(r'E:\Project\FlexPPI')
B=ROOT/'paper/revision_v14_20261010/section33_expansion_20261025/feature_baseline'
KEY=['sample_id','residue_index']
METHODS=['full','intra','static47_sasa','aa_sasa_contacts','aa_mean']

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()

def js(p,obj):
    p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')

def weighted_means(rows,target,mask):
    active=rows.loc[mask,['record_index','dataset','biological_system']].groupby('record_index').agg(dataset=('dataset','first'),system=('biological_system','first'),n=('dataset','size')).reset_index()
    nsystem=active[['dataset','system']].drop_duplicates().groupby('dataset').size()
    nrecord=active.groupby(['dataset','system']).size()
    weights={int(x.record_index):1/(len(nsystem)*nsystem[x.dataset]*nrecord[(x.dataset,x.system)]*x.n) for x in active.itertuples()}
    a=rows.loc[mask,['record_index','aatype',target]].copy()
    a['w']=a.record_index.map(weights)
    a['wy']=a.w*a[target]
    return a.groupby('aatype').wy.sum()/a.groupby('aatype').w.sum()

def record_maes(frame,label='original'):
    result=[]
    for (sid,role),a in frame.groupby(['sample_id','role']):
        for metric in ['rmsf','hchi']:
            m=a['valid_'+metric].to_numpy(bool)&np.isfinite(a['md_'+metric])
            for method in METHODS:
                yy=a.loc[m,'md_'+metric].to_numpy()
                pp=a.loc[m,'pred_'+method+'_'+metric].to_numpy()
                assert np.isfinite(pp).all()
                result.append(dict(sample_id=sid,role=role,pdb_id=a.pdb_id.iloc[0],component_id=a.component_id.iloc[0],fold=int(a.fold.iloc[0]),metric=metric,method=method,label=label,n=int(m.sum()),mae=float(np.abs(pp-yy).mean()) if len(yy) else np.nan))
    return pd.DataFrame(result)

def cluster(p):
    # Independently formed PDB means and cluster-ratio interval.
    grp=[a.gain.to_numpy() for _,a in p.groupby('component_id')]
    cnt=np.array([len(x) for x in grp]); sums=np.array([x.sum() for x in grp]); means=sums/cnt
    rng=np.random.default_rng(20261026)
    draw=rng.integers(0,len(grp),size=(10000,len(grp)))
    boot=sums[draw].sum(1)/cnt[draw].sum(1)
    return dict(gain=float(p.gain.mean()),lo=float(np.quantile(boot,.025)),hi=float(np.quantile(boot,.975)),family_gain=float(means.mean()),n_pdb=len(p),n_family=len(grp),n_positive=int(p.gain.gt(0).sum()))

def comparisons(rec,label):
    rows=[]
    for (role,metric),a in rec.groupby(['role','metric']):
        w=a.pivot(index=['sample_id','pdb_id','component_id'],columns='method',values='mae')
        for comparator in METHODS[1:]:
            s=w[['full',comparator]].dropna().copy()
            s['gain']=s[comparator]-s.full
            p=s.reset_index().groupby(['pdb_id','component_id'])[['full',comparator,'gain']].mean().reset_index()
            leave=[]
            for fam in p.component_id.unique():
                z=p[~p.component_id.eq(fam)]
                leave.append((fam,float(z.gain.mean())))
            leave.sort(key=lambda x:x[1])
            rows.append(dict(label=label,role=role,metric=metric,comparison='full_vs_'+comparator,full_mae=float(p.full.mean()),comparator_mae=float(p[comparator].mean()),**cluster(p),leave_one_family_gain_min=leave[0][1],leave_one_family_gain_min_excluded=leave[0][0],leave_one_family_gain_max=leave[-1][1],leave_one_family_gain_max_excluded=leave[-1][0]))
    return pd.DataFrame(rows)

def main():
    OUT.mkdir(exist_ok=True)
    f=pd.read_parquet(R/'controls/candidate_with_baselines56.parquet')
    old=pd.read_parquet(R/'shared/oof_candidates56.parquet')
    fam=pd.read_csv(ROOT/'experiments/_sources/family_folds/record_family_folds.csv')
    allf=pd.read_parquet(B/'features_all.parquet')
    contract=json.loads((R/'controls/baseline_contract_family.json').read_text(encoding='utf-8'))
    progress=json.loads((R/'controls/training_progress.json').read_text(encoding='utf-8'))
    assert not f.duplicated(KEY).any()
    assert not fam.duplicated('sample_id').any()
    mapfold=fam.set_index('sample_id').fold
    mapfam=fam.set_index('sample_id').component_id
    assert f.fold.eq(f.sample_id.map(mapfold)).all()
    assert f.component_id.eq(f.sample_id.map(mapfam)).all()
    assert fam.groupby('component_id').fold.nunique().max()==1
    assert sha(ROOT/'experiments/_sources/family_folds/record_family_folds.csv')==contract['family_split_sha256']
    audit=dict(n_rows=len(f),n_records=int(f.sample_id.nunique()),n_pdb=int(f.pdb_id.nunique()),n_family=int(f.component_id.nunique()),duplicated_residue_keys=0,original_vs_family_fold_changed_records=int(f[f.fold.ne(f.original_fold)].sample_id.nunique()),frozen_contract_hash_verified=True,jobs_complete=progress['status']=='complete',job_count=len(progress['jobs']),intermediate_prediction_columns=[x for x in pd.read_parquet(R/'controls/candidate_features56.parquet').columns if x.startswith('pred_')])
    # Trainer rewrites this intermediate table after fitting the first *new* family baseline.
    audit['intermediate_prediction_columns_are_new_family_cache']=audit['intermediate_prediction_columns']==['pred_static47_sasa_rmsf','pred_static47_sasa_hchi']
    assert audit['intermediate_prediction_columns_are_new_family_cache']
    OOF=[]
    for variant,method in [('full','full'),('intrachain','intra')]:
        o=old[old.variant.eq(variant)].copy()
        assert not o.duplicated(KEY).any()
        x=f[KEY+['chain','author_residue_id','amino_acid','fold','valid_rmsf','valid_hchi','md_rmsf','md_hchi']].merge(o,on=KEY,validate='one_to_one',suffixes=('_feature','_oof'))
        assert len(x)==len(f)
        for name in ['chain','amino_acid','fold']:
            assert x[name+'_feature'].eq(x[name+'_oof']).all()
        assert x.author_residue_id_feature.astype(str).eq(x.author_residue_id_oof.astype(str)).all()
        err={}
        for metric,raw in [('rmsf','rmsf'),('hchi','chi_entropy')]:
            valid=x['valid_'+metric+('_feature' if metric=='rmsf' else '')].to_numpy(bool)
            assert np.array_equal(valid,x['valid_'+raw+('_oof' if raw=='rmsf' else '')])
            err[metric]=float(np.max(np.abs(x.loc[valid,'md_'+metric]-x.loc[valid,'target_'+raw])))
            assert err[metric]<1e-5
        f=f.merge(o[KEY+['pred_rmsf','pred_chi_entropy']].rename(columns={'pred_rmsf':'pred_'+method+'_rmsf','pred_chi_entropy':'pred_'+method+'_hchi'}),on=KEY,validate='one_to_one')
        OOF.append(dict(variant=variant,target_max_abs_difference=err,fold_mask_and_chemical_keys_match=True))
    audit['oof_integrity']=OOF
    # Validate row-index alignment and per-fold prediction caches over all 961321 residues.
    rowbase=np.zeros(allf.record_index.max()+1,dtype=int)
    first=allf.reset_index().groupby('record_index')['index'].min()
    rowbase[first.index.to_numpy()]=first.to_numpy()
    rowids=rowbase[allf.record_index.to_numpy()]+allf.residue_index.to_numpy()
    assert np.array_equal(rowids,np.arange(len(allf)))
    fold=allf.sample_id.map(mapfold).to_numpy()
    recordids=fam.groupby('component_id').sample_id.nunique().sort_values(ascending=False)
    global_largest=recordids.index[0]
    audit['global_largest_family']=dict(component_id=global_largest,n_records=int(recordids.iloc[0]),fold=int(fam.loc[fam.component_id.eq(global_largest),'fold'].iloc[0]),candidate_records=int(f[f.component_id.eq(global_largest)].sample_id.nunique()),candidate_pdb=int(f[f.component_id.eq(global_largest)].pdb_id.nunique()))
    jobaudit=[]
    for j in progress['jobs']:
        mi=['rmsf','hchi'].index(j['metric']); metric=j['metric']; ff=j['fold'];variant=j['variant']
        valid=allf['valid_'+metric].to_numpy(bool)
        train=valid&(fold!=ff);test=valid&(fold==ff)
        tr=set(fam.loc[fam.fold.ne(ff),'component_id']);te=set(fam.loc[fam.fold.eq(ff),'component_id'])
        assert not tr&te
        assert j['n_train']==int(train.sum()) and j['n_test']==int(test.sum())
        assert j['n_train_records']==allf.loc[train,'sample_id'].nunique() and j['n_test_records']==allf.loc[test,'sample_id'].nunique()
        cache=np.load(R/f'controls/models/{variant}_familyfold{ff}_{metric}.npz')
        assert np.array_equal(cache['rowids'],rowids[test])
        assert str(cache['contract_sha256'])==sha(R/'controls/baseline_contract_family.json')
        joined=np.load(R/f'controls/{variant}_family_oof.npz')
        assert np.array_equal(joined['rowids'],rowids)
        assert np.array_equal(joined['prediction'][test,mi],cache['prediction'])
        candindex=f.merge(allf.reset_index()[KEY+['index']],on=KEY,validate='one_to_one')['index'].to_numpy()
        vp=f['valid_'+metric].to_numpy(bool)
        assert np.array_equal(f.loc[vp,'pred_'+variant+'_'+metric].to_numpy(),joined['prediction'][candindex[vp],mi])
        jobaudit.append(dict(variant=variant,fold=ff,metric=metric,n_train=int(train.sum()),n_test=int(test.sum()),family_overlap=len(tr&te),cache_rowkeys_and_prediction_values_match=True))
    audit['jobs']=jobaudit
    audit['all_fits_verified']=len(jobaudit)==20
    # Independent same-fold AA means, including source/system/record/residue weighting.
    aa=np.load(R/'controls/aa_mean_family_oof.npz')['prediction']
    aadiff=[]
    for ff in range(5):
        for mi,metric in enumerate(['rmsf','hchi']):
            mask=(fold!=ff)&allf['valid_'+metric].to_numpy(bool)
            means=weighted_means(allf,'md_'+metric,mask)
            test=(fold==ff)&allf['valid_'+metric].to_numpy(bool)
            expected=allf.loc[test,'aatype'].map(means).to_numpy()
            d=float(np.abs(expected-aa[test,mi]).max())
            assert d<1e-6
            aadiff.append(dict(fold=ff,metric=metric,max_abs_difference=d))
    audit['aa_means_independently_verified']=aadiff
    rec=record_maes(f)
    rec.to_csv(OUT/'independent_record_mae.csv',index=False)
    original=comparisons(rec,'original')
    ref=pd.read_csv(R/'controls/paired_full_vs_baselines_family_cluster_ci.csv').query("subset=='all'")
    check=original.merge(ref,on=['role','metric','comparison'],suffixes=('_independent','_source'))
    audit['metric_max_abs_diff_to_main']=float(abs(check.full_mae-check.full_pdb_macro_mae).max())
    audit['gain_max_abs_diff_to_main']=float(abs(check.gain-check.pdb_macro_gain).max())
    audit['interval_max_abs_diff_to_main']=float(max(abs(check.lo-check.pdb_gain_ci_lo).max(),abs(check.hi-check.pdb_gain_ci_hi).max()))
    assert audit['metric_max_abs_diff_to_main']<1e-12 and audit['interval_max_abs_diff_to_main']<1e-12
    original.to_csv(OUT/'independent_main_comparisons.csv',index=False)
    # Correctly join independently recomputed common-duration labels: source fold is not family fold.
    u=pd.read_csv(R/'uniform_md/all_uniform_labels_per_residue.csv.gz',low_memory=False)
    u=u[u.variant.eq('uniform_rep1')].rename(columns={'fold':'original_label_split_fold'})
    assert not u.duplicated(KEY).any()
    z=f.merge(u[KEY+['chain','author_residue_id','amino_acid','family_id','original_label_split_fold','rmsf','chi_entropy','valid_rmsf','valid_chi_entropy']],on=KEY,validate='one_to_one',suffixes=('','_uniform'))
    for n in ['chain','amino_acid']:assert z[n].eq(z[n+'_uniform']).all()
    assert z.author_residue_id.astype(str).eq(z.author_residue_id_uniform.astype(str)).all()
    assert z.component_id.eq(z.family_id).all()
    audit['uniform']=dict(n_records=int(z.sample_id.nunique()),n_pdb=int(z.pdb_id.nunique()),n_family=int(z.component_id.nunique()),n_rows=len(z),family_fold_retained=z.fold.eq(z.sample_id.map(mapfold)).all().item(),source_label_fold_differs_from_family_records=int(z[z.fold.ne(z.original_label_split_fold)].sample_id.nunique()),old_vs_uniform_rmsf_mask_changes=int(z.valid_rmsf.ne(z.valid_rmsf_uniform).sum()),old_vs_uniform_hchi_mask_changes=int(z.valid_hchi.ne(z.valid_chi_entropy).sum()))
    z['valid_rmsf']=z.valid_rmsf_uniform;z['valid_hchi']=z.valid_chi_entropy
    z['md_rmsf']=z.rmsf;z['md_hchi']=z.chi_entropy
    uniform=comparisons(record_maes(z,'uniform_rep1'),'uniform_rep1')
    uniform.to_csv(OUT/'independent_uniform_comparisons.csv',index=False)
    annotation=pd.read_csv(R/'biological_annotation/all56_biological_annotations.csv')
    assert not annotation.sample_id.duplicated().any()
    coreids=annotation.loc[annotation.core_functional_peptide,'sample_id']
    core_original=comparisons(record_maes(f[f.sample_id.isin(coreids)],'core_original'),'core_original')
    core_uniform=comparisons(record_maes(z[z.sample_id.isin(coreids)],'core_uniform_rep1'),'core_uniform_rep1')
    core_original.to_csv(OUT/'independent_core_original_comparisons.csv',index=False)
    core_uniform.to_csv(OUT/'independent_core_uniform_comparisons.csv',index=False)
    audit['uniform']['interpretation']='Frozen models trained on original (often 500 ns) labels evaluated against recomputed 100 ns reference: time-scale robustness, not newly trained model performance.'
    audit['core_uniform']=dict(n_records=int(z[z.sample_id.isin(coreids)].sample_id.nunique()),n_pdb=int(z[z.sample_id.isin(coreids)].pdb_id.nunique()),n_family=int(z[z.sample_id.isin(coreids)].component_id.nunique()))
    # AA domination: variance explained by per-AA means vs within-AA differentiation.
    aarows=[]
    for role,a in f[f.valid_hchi].groupby('role'):
        for method in METHODS:
            target=a.md_hchi.to_numpy();pred=a['pred_'+method+'_hchi'].to_numpy()
            tm=a.groupby('amino_acid').md_hchi.transform('mean').to_numpy()
            pm=a.groupby('amino_acid')['pred_'+method+'_hchi'].transform('mean').to_numpy()
            aarows.append(dict(role=role,method=method,n_residues=len(a),pooled_rho=float(spearmanr(target,pred).statistic),within_AA_centered_rho=float(spearmanr(target-tm,pred-pm).statistic),target_variance_fraction_between_AA=float(np.var(tm)/np.var(target)),prediction_variance_fraction_between_AA=float(np.var(pm)/np.var(pred))))
    pd.DataFrame(aarows).to_csv(OUT/'AA_composition_diagnostic_descriptive_only.csv',index=False)
    # Matched pairs: verify no pseudo-independent inference and recompute from stored pairs.
    pairs=pd.read_parquet(R/'controls/sameAA_backbone_matched_all_pairs.parquet')
    assert not pairs.duplicated(['sample_id','role','rmsf_tolerance_A','residue_i','residue_j']).any()
    for side in ['i','j']:
        pp=f[KEY+['amino_acid','role','md_rmsf','md_hchi']].rename(columns={'residue_index':'residue_'+side,'amino_acid':'aa_'+side,'role':'role_'+side,'md_rmsf':'rmsf_'+side,'md_hchi':'hchi_'+side})
        pairs=pairs.merge(pp,on=['sample_id','residue_'+side],validate='many_to_one')
    assert pairs.aa_i.eq(pairs.aa_j).all() and pairs.role.eq(pairs.role_i).all() and pairs.role.eq(pairs.role_j).all()
    assert (abs(pairs.rmsf_j-pairs.rmsf_i)<=pairs.rmsf_tolerance_A+1e-12).all()
    assert np.allclose(pairs.hchi_j-pairs.hchi_i,pairs.target_hchi_difference,atol=1e-12,rtol=0)
    p25=pairs[pairs.rmsf_tolerance_A.eq(.25)]
    pairrows=[]
    for role,a in p25.groupby('role'):
        for method in METHODS[1:]:
            b=a.assign(gain=a[method+'_absolute_error']-a.full_absolute_error).groupby(['sample_id','pdb_id','component_id']).gain.mean().reset_index()
            p=b.groupby(['pdb_id','component_id']).gain.mean().reset_index()
            pairrows.append(dict(role=role,comparison='full_vs_'+method,**cluster(p)))
    pd.DataFrame(pairrows).to_csv(OUT/'independent_matched_pair_comparisons.csv',index=False)
    audit['matched_pairs']=dict(n_pairs_all_tolerances=len(pairs),unique_keys=True,same_AA_role_and_tolerance_verified=True,no_residue_or_pair_p_values=True,pairs_at_025=int(len(p25)),records_at_025=int(p25.sample_id.nunique()))
    duplicate_pdb=f[['sample_id','pdb_id','component_id']].drop_duplicates().groupby('pdb_id').filter(lambda a:len(a)>1)
    duplicate_pdb.to_csv(OUT/'duplicate_pdb_records.csv',index=False)
    famsizes=f[['sample_id','pdb_id','component_id']].drop_duplicates().groupby('component_id').agg(n_records=('sample_id','nunique'),n_pdb=('pdb_id','nunique')).sort_values(['n_pdb','n_records'],ascending=False)
    famsizes.to_csv(OUT/'candidate_family_sizes.csv')
    audit['candidate_largest_cluster']=dict(component_id=famsizes.index[0],n_records=int(famsizes.iloc[0].n_records),n_pdb=int(famsizes.iloc[0].n_pdb),pdb_weight_fraction=float(famsizes.iloc[0].n_pdb/f.pdb_id.nunique()))
    audit['duplicate_pdb_names']=sorted(duplicate_pdb.pdb_id.unique().tolist())
    js(OUT/'independent_statistics_audit.json',audit)
    print(original.to_string(index=False))
    print('UNIFORM')
    print(uniform.to_string(index=False))
    print('CORE UNIFORM')
    print(core_uniform.to_string(index=False))
    print('AUDIT',json.dumps({k:v for k,v in audit.items() if k not in ['jobs','aa_means_independently_verified']},ensure_ascii=False))

if __name__=='__main__':main()
