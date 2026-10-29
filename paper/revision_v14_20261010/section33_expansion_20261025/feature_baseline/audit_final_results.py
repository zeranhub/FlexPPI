"""Additional frozen-contract and same-amino-acid OOF audits (no fitting)."""
from pathlib import Path
from collections import Counter
import json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import spearmanr
import mdtraj
import sklearn
import scipy
import joblib

OUT=Path(__file__).resolve().parent
METHODS=['full','intra','static','aa_only','aa_sasa','aa_sasa_contacts','static47_sasa']


def rho(a,b):
    return float(spearmanr(a,b).statistic) if len(a)>=4 and np.ptp(a)>0 and np.ptp(b)>0 else np.nan


def main():
    contract=json.loads((OUT/'features_all_ready.json').read_text(encoding='utf-8'))
    assert contract['records']==3063 and not contract['missing']
    columns=['sample_id','record_index','residue_index','dataset','fold','amino_acid','chain',
             'sidechain_complete','chi_input_complete','nonstandard_residue','discarded_altloc_atoms',
             'sidechain_sasa_complex_A2','sidechain_sasa_isolated_A2','delta_sidechain_sasa_A2',
             'residue_sasa_complex_A2','residue_sasa_isolated_A2','delta_residue_sasa_A2']
    for metric in ['rmsf','hchi']:
        columns+=['valid_'+metric,'md_'+metric]+['pred_'+method+'_'+metric for method in METHODS[:4]]
    frame=pq.read_table(OUT/'features_all.parquet',columns=columns).to_pandas()
    assert len(frame)==961321 and frame.sample_id.nunique()==3063
    assert frame.valid_rmsf.sum()==961321 and frame.valid_hchi.sum()==824988
    assert not frame.duplicated(['sample_id','residue_index']).any()
    records=json.loads((OUT/'cohort_records.json').read_text(encoding='utf-8'))
    global_rows=np.array([record['start'] for record in records])[frame.record_index.to_numpy()]+frame.residue_index.to_numpy()
    assert np.array_equal(global_rows,np.arange(len(frame)))
    for method in METHODS[4:]:
        with np.load(OUT/(method+'_oof.npz')) as values:
            assert np.array_equal(values['rowids'],global_rows)
            for i,metric in enumerate(['rmsf','hchi']):
                p=values['prediction'][:,i]
                assert np.isfinite(p[frame['valid_'+metric]]).all()
                frame['pred_'+method+'_'+metric]=p
    quality=pd.read_csv(OUT/'quality_all.csv')
    assert len(quality)==3063 and set(quality.status)=={'ok'}
    single=frame.groupby('sample_id').chain.nunique()
    single_ids=set(single[single==1].index)
    singlemask=frame.sample_id.isin(single_ids)
    assert (frame.loc[singlemask,'delta_sidechain_sasa_A2']==0).all()
    assert (frame.loc[singlemask,'delta_residue_sasa_A2']==0).all()
    for name in ['sidechain','residue']:
        assert (frame['delta_'+name+'_sasa_A2']>=0).all()
        assert np.allclose(frame[name+'_sasa_isolated_A2']-frame[name+'_sasa_complex_A2'],frame['delta_'+name+'_sasa_A2'],atol=.001,rtol=0)
    model_audit=[]
    for method in METHODS[4:]:
        for fold in range(5):
            for metric in ['rmsf','hchi']:
                path=OUT/'models'/f'{method}_fold{fold}_{metric}.joblib'
                model=joblib.load(path)
                assert model.n_iter_==120 and model.get_params()['early_stopping'] is False
                model_audit.append(dict(model=path.name,n_iter=model.n_iter_,n_features=model.n_features_in_,random_state=model.random_state))
    audit=dict(status='passed',records=3063,rows=len(frame),rmsf_valid=int(frame.valid_rmsf.sum()),hchi_valid=int(frame.valid_hchi.sum()),
               singlechain_records=len(single_ids),singlechain_rows=int(singlemask.sum()),singlechain_delta_exactly_zero=True,
               incomplete_sidechains=int((~frame.sidechain_complete).sum()),
               incomplete_valid_hchi=int((~frame.sidechain_complete&frame.valid_hchi).sum()),
               zero_sasa_complete_hchi=int(((frame.sidechain_sasa_complex_A2==0)&frame.sidechain_complete&frame.valid_hchi).sum()),
               negative_delta_before_clipping=int(quality.negative_delta_count.sum()),
               min_delta_before_clipping=float(quality.negative_delta_min.min()),
               nonstandard_residues=int(frame.nonstandard_residue.sum()),discarded_altloc_atoms=int(frame.discarded_altloc_atoms.sum()),
               versions=dict(numpy=np.__version__,pandas=pd.__version__,mdtraj=mdtraj.__version__,sklearn=sklearn.__version__,scipy=scipy.__version__),
               no_test_model_selection=True,no_modification_to_original_inputs=True,trained_models=model_audit)
    (OUT/'final_integrity_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    quality.groupby('dataset')[['n','incomplete_sidechains','incomplete_valid_hchi','chi_incomplete_valid_hchi','nonstandard','discarded_altloc_atoms','heavy_atoms','negative_delta_count']].sum().to_csv(OUT/'quality_by_source.csv')
    rows=[]
    for sid,part in frame.groupby('sample_id',sort=False):
        for metric in ['rmsf','hchi']:
            valid=part['valid_'+metric].to_numpy()
            aa=part.amino_acid.to_numpy()
            counts=Counter(aa[valid])
            for label,keep in [('all_types',valid),('type_n_ge_4',valid&np.array([counts[a]>=4 for a in aa])),
                               ('complete_sidechains',valid&part.sidechain_complete.to_numpy())]:
                if keep.sum()<4:
                    continue
                types=aa[keep]
                y=part['md_'+metric].to_numpy()[keep].astype(np.float64)
                centered_truth=y.copy()
                indices={a:types==a for a in set(types)}
                for a,selected in indices.items():
                    centered_truth[selected]=0 if np.ptp(centered_truth[selected])==0 else centered_truth[selected]-centered_truth[selected].mean()
                for method in METHODS:
                    prediction=part['pred_'+method+'_'+metric].to_numpy()[keep].astype(np.float64)
                    centered_prediction=prediction.copy()
                    for a,selected in indices.items():
                        centered_prediction[selected]=0 if np.ptp(centered_prediction[selected])==0 else centered_prediction[selected]-centered_prediction[selected].mean()
                    train_mean=part['pred_aa_only_'+metric].to_numpy()[keep]
                    rows.append(dict(sample_id=sid,dataset=part.dataset.iloc[0],fold=int(part.fold.iloc[0]),
                                     metric=metric,variant=method,subset=label,n=int(keep.sum()),
                                     rho_raw=rho(y,prediction),rho_train_aa_centered=rho(y-train_mean,prediction-train_mean),
                                     rho_within_record_aa_centered=rho(centered_truth,centered_prediction),mae=float(np.abs(prediction-y).mean())))
    result=pd.DataFrame(rows)
    result.to_csv(OUT/'same_type_per_record_metrics.csv',index=False)
    result.groupby(['variant','metric','subset','dataset'])[['rho_raw','rho_train_aa_centered','rho_within_record_aa_centered','mae']].mean().to_csv(OUT/'same_type_macro_by_source.csv')
    result.groupby(['variant','metric','subset'])[['rho_raw','rho_train_aa_centered','rho_within_record_aa_centered','mae']].mean().to_csv(OUT/'same_type_macro_all.csv')
    print(json.dumps(audit,indent=2))


if __name__=='__main__':main()
