"""Family-frozen strong baseline, read-only source inputs; sequential <=2 threads."""
import os
for key in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='2'
from pathlib import Path
import json, time, hashlib, gc
from collections import Counter
import numpy as np
import pandas as pd
import joblib, sklearn, psutil
from sklearn.ensemble import HistGradientBoostingRegressor
from threadpoolctl import threadpool_limits

OUT=Path(__file__).resolve().parent
ROOT=Path(r'E:\Project\FlexPPI')
B=ROOT/'paper/revision_v14_20261010/section33_expansion_20261025/feature_baseline'
R=OUT.parent.parent
FAMILY=ROOT/'experiments/_sources/family_folds'
CONFIG=dict(loss='squared_error',learning_rate=.07,max_iter=120,max_leaf_nodes=31,max_depth=None,min_samples_leaf=50,l2_regularization=1.,max_bins=63,early_stopping=False,random_state=7)
METRICS=['rmsf','hchi']

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(2**20),b''): h.update(block)
    return h.hexdigest()

def write_json(path,value):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,ensure_ascii=False),encoding='utf-8')
    tmp.replace(path)

def balanced_weights(records,record_ids,mask):
    counts=np.bincount(record_ids[mask],minlength=len(records))
    active=[r for r in records if counts[r['record_index']]>0]
    systems=Counter((r['dataset'],r['biological_system']) for r in active)
    sources=Counter(source for source,system in systems)
    weight=np.zeros(len(records),float)
    for r in active:
        i=r['record_index']; weight[i]=1/(len(sources)*sources[r['dataset']]*systems[(r['dataset'],r['biological_system'])]*counts[i])
    w=weight[record_ids[mask]]
    return w/w.mean()

def main():
    tick=time.monotonic()
    OUT.mkdir(exist_ok=True)
    models=OUT/'models'; models.mkdir(exist_ok=True)
    source_contract=json.loads((B/'baseline_contract.json').read_text(encoding='utf-8'))
    sets={name:source_contract['feature_sets'][name] for name in ['static47_sasa','aa_sasa_contacts']}
    records=json.loads((B/'cohort_records.json').read_text(encoding='utf-8'))
    f=pd.read_parquet(B/'features_all.parquet')
    cand=pd.read_csv(R/'inventory_agent/short_peptide_5to40_candidates.csv')
    cand=cand[cand.two_chain_unambiguous].copy()
    assert len(cand)==cand.sample_id.nunique()==56 and cand.pdb_id.nunique()==54 and cand.family_id.nunique()==46
    fam=pd.read_csv(FAMILY/'record_family_folds.csv')
    assert len(fam)==fam.sample_id.nunique()==3063
    assert set(fam.sample_id)==set(f.sample_id)
    assert fam.groupby('component_id')['fold'].nunique().max()==1
    mapfold=fam.set_index('sample_id')['fold']; mapfam=fam.set_index('sample_id')['component_id']
    f['original_fold']=f['fold']; f['fold']=f.sample_id.map(mapfold).astype(np.int8)
    f['component_id']=f.sample_id.map(mapfam)
    assert not f.duplicated(['sample_id','residue_index']).any()
    assert all(f.iloc[r['start']:r['stop']].sample_id.eq(r['sample_id']).all() for r in records)
    assert all(int(fam.set_index('sample_id').loc[r.sample_id,'fold'])==int(r.fold) for r in pd.read_csv(R/'short_peptide_preliminary/all56_records_two_roles.csv').itertuples())
    assert all(mapfam[r.sample_id]==r.family_id for r in cand.itertuples())
    iscandidate=f.sample_id.isin(cand.sample_id).to_numpy()
    rolemap=pd.read_csv(R/'short_peptide_preliminary/all56_records_two_roles.csv')[['sample_id','chain','role','pdb_id']]
    c=f[iscandidate].merge(rolemap,on=['sample_id','chain'],how='left',validate='many_to_one')
    assert c.role.notna().all()
    assert len(c)==iscandidate.sum()
    c=c.drop(columns=[col for col in c if col.startswith('pred_')]) # old-fold predictions prohibited
    c.to_parquet(OUT/'candidate_features56.parquet',index=False)
    fold=f['fold'].to_numpy(); rid=f.record_index.to_numpy()
    original_starts=np.array([r['start'] for r in records])
    rowids=original_starts[rid]+f.residue_index.to_numpy()
    assert np.array_equal(rowids,np.arange(len(f)))
    targets=np.stack([f['md_'+m].to_numpy() for m in METRICS],axis=1)
    valid=np.stack([f['valid_'+m].to_numpy() for m in METRICS],axis=1)
    labels_finite=np.isfinite(targets)
    assert labels_finite[valid].all()
    contract=dict(protocol='frozen family splits; fit on other four folds of all3063 original-label records',configuration=CONFIG,feature_sets=sets,threads=2,no_early_stopping=True,no_row_subsampling=True,max_training_rows=None,all_records=3063,all_rows=len(f),candidate_records=56,candidate_pdb=54,candidate_family=46,valid_counts=valid.sum(axis=0).tolist(),input_information='full protein heavy atoms incl sidechains; stronger-input explanatory comparator',family_split_sha256=sha(FAMILY/'record_family_folds.csv'),features_sha256=sha(B/'features_all.parquet'),old_baseline_contract_sha256=sha(B/'baseline_contract.json'),sklearn_version=sklearn.__version__,feature_matrix_MiB=len(f)*len(sets['static47_sasa'])*4/2**20,weighting='source/system/record/valid-residue balanced; same as original static baseline',no_test_tuning=True,secondary_structure='raw backbone phi/psi sin/cos plus N-CA-C angle, no DSSP calls',subsets='all, interface_ca8, noninterface_ca8, complete_sidechains')
    write_json(OUT/'baseline_contract_family.json',contract)
    progress=[]
    aa=f['aatype'].to_numpy()
    aa_prediction=np.full(targets.shape,np.nan,np.float32)
    for ff in range(5):
        for mi,m in enumerate(METRICS):
            train=(fold!=ff)&valid[:,mi]; test=(fold==ff)&valid[:,mi]
            w=balanced_weights(records,rid,train)
            means=np.bincount(aa[train],weights=w*targets[train,mi],minlength=20)/np.bincount(aa[train],weights=w,minlength=20)
            aa_prediction[test,mi]=means[aa[test]]
    np.savez_compressed(OUT/'aa_mean_family_oof.npz',prediction=aa_prediction,rowids=rowids)
    with threadpool_limits(limits=2):
        for variant,cols in sets.items():
            x=f[cols].to_numpy(dtype=np.float32)
            prediction=np.full(targets.shape,np.nan,np.float32)
            for ff in range(5):
                train_fam=set(fam.loc[fam.fold!=ff,'component_id']); test_fam=set(fam.loc[fam.fold==ff,'component_id'])
                assert not train_fam&test_fam
                for mi,m in enumerate(METRICS):
                    stem=f'{variant}_familyfold{ff}_{m}'
                    cache=models/(stem+'.npz')
                    train=(fold!=ff)&valid[:,mi]; test=(fold==ff)&valid[:,mi]
                    if cache.exists():
                        z=np.load(cache)
                        assert np.array_equal(z['rowids'],rowids[test])
                        assert str(z['contract_sha256'])==sha(OUT/'baseline_contract_family.json')
                        pred=z['prediction']; elapsed=0
                    else:
                        start=time.monotonic()
                        model=HistGradientBoostingRegressor(**CONFIG,categorical_features=[name.startswith('aa_offset_') for name in cols])
                        model.fit(x[train],targets[train,mi],sample_weight=balanced_weights(records,rid,train))
                        pred=model.predict(x[test]).astype(np.float32)
                        assert np.isfinite(pred).all()
                        joblib.dump(model,models/(stem+'.joblib'),compress=3)
                        np.savez_compressed(cache,rowids=rowids[test],prediction=pred,contract_sha256=np.array(sha(OUT/'baseline_contract_family.json')))
                        elapsed=time.monotonic()-start
                        del model; gc.collect()
                    prediction[test,mi]=pred
                    info=dict(variant=variant,fold=ff,metric=m,n_train=int(train.sum()),n_test=int(test.sum()),n_train_records=int(f.loc[train,'sample_id'].nunique()),n_test_records=int(f.loc[test,'sample_id'].nunique()),n_train_family=len(train_fam),n_test_family=len(test_fam),family_overlap=0,seconds=elapsed,process_memory_MiB=psutil.Process().memory_info().rss/2**20,peak_memory_MiB=psutil.Process().memory_info().peak_wset/2**20)
                    progress.append(info)
                    write_json(OUT/'training_progress.json',dict(status='running',jobs=progress,total_seconds=time.monotonic()-tick))
                    print('BASELINE',json.dumps(info),flush=True)
            assert np.isfinite(prediction[valid]).all()
            np.savez_compressed(OUT/(variant+'_family_oof.npz'),prediction=prediction,rowids=rowids)
            c.to_parquet(OUT/'candidate_features56.parquet',index=False)
            for mi,m in enumerate(METRICS): c['pred_'+variant+'_'+m]=prediction[iscandidate,mi]
            c.to_parquet(OUT/'candidate_with_baselines56.parquet',index=False)
            del x; gc.collect()
    for mi,m in enumerate(METRICS): c['pred_aa_mean_'+m]=aa_prediction[iscandidate,mi]
    c.to_parquet(OUT/'candidate_with_baselines56.parquet',index=False)
    write_json(OUT/'training_progress.json',dict(status='complete',jobs=progress,total_seconds=time.monotonic()-tick))
    print('COMPLETE',len(progress),round(time.monotonic()-tick,1),flush=True)
if __name__=='__main__': main()
