"""Read-only independent audit of saved baseline scores; no completion checker."""
from pathlib import Path
import datetime
import hashlib
import json
import numpy as np
import pandas as pd

TASK = Path(r'E:\Project\FlexPPI\experiments\_sources')
BASE = TASK/'baselines'
OUT = TASK/'completed_results_20261012'
BF = Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\backflip2_corrected_scoring')
METRICS = ['mae','rmse','bias','spearman','pearson','z_mae','z_rmse','top10_precision']
METHODS = ['flexppi','intrachain','planet_full','planet_seq','rmsf_net_pdb','anm','gnm']
RELATIVE = {'anm','gnm','rmsf_net_pdb'}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def close(actual,expected):
    np.testing.assert_allclose(np.asarray(actual,dtype=float),np.asarray(expected,dtype=float),rtol=1e-11,atol=1e-11,equal_nan=True)


def raw_statistics(target,prediction,absolute):
    target=np.asarray(target,dtype=float); prediction=np.asarray(prediction,dtype=float)
    target_z=(target-target.mean())/target.std(ddof=0)
    prediction_z=(prediction-prediction.mean())/prediction.std(ddof=0)
    k=max(1,int(np.ceil(.10*len(target))))
    target_top=set(np.argsort(-target,kind='stable')[:k]); pred_top=set(np.argsort(-prediction,kind='stable')[:k])
    result=dict(spearman=pd.Series(target).rank().corr(pd.Series(prediction).rank()),
                pearson=pd.Series(target).corr(pd.Series(prediction)),
                z_mae=float(np.mean(abs(target_z-prediction_z))),
                z_rmse=float(np.sqrt(np.mean((target_z-prediction_z)**2))),
                top10_precision=len(target_top&pred_top)/k)
    if absolute:
        result.update(mae=float(np.mean(abs(target-prediction))),
                      rmse=float(np.sqrt(np.mean((target-prediction)**2))),bias=float(np.mean(prediction-target)))
    return result


def main():
    hashes={}
    for line in (OUT/'baseline_remote_sha256.txt').read_text(encoding='utf-8-sig').splitlines():
        expected,remote=line.split(None,1)
        relative=remote.split('/baselines/',1)[1]
        local=BASE/relative
        assert sha(local)==expected,(relative,expected,sha(local))
        hashes[str(local)]=expected
    protocols={}; per_by_protocol={}; spot_checks=[]
    for protocol,folder in [('original','additional_original_oof'),('family','additional_family_oof')]:
        root=BASE/folder
        per=pd.read_csv(root/'per_record_metrics.csv',float_precision='round_trip')
        macro=pd.read_csv(root/'record_macro_metrics.csv',float_precision='round_trip').set_index(['method','dataset'])
        fold_macro=pd.read_csv(root/'source_fold_macro_metrics.csv',float_precision='round_trip').set_index(['method','dataset','fold'])
        coverage=pd.read_csv(root/'coverage_audit.csv',keep_default_na=False)
        coverage_summary=pd.read_csv(root/'coverage_summary.csv').set_index('method')
        audit=json.loads((root/'scoring_audit.json').read_text())
        assert audit['all_methods_complete'] and audit['complete_required']
        assert len(per)==3063*7 and not per.duplicated(['method','sample_id']).any()
        assert set(per.method)==set(METHODS) and per.common_fraction.eq(1.0).all()
        assert per['count'].eq(per.expected_valid).all()
        assert coverage.complete.all() and coverage.error.eq('').all()
        assert coverage.expected_valid.eq(coverage.covered_valid).all()
        assert coverage_summary.records.eq(3063).all() and coverage_summary.complete_records.eq(3063).all()
        assert coverage_summary.expected_valid.eq(961321).all() and coverage_summary.covered_valid.eq(961321).all()
        for method in RELATIVE:
            assert per[per.method==method][['mae','rmse','bias']].isna().all().all()
            assert not audit['methods'][method]['absolute_angstrom']
        for method in set(METHODS)-RELATIVE:
            assert audit['methods'][method]['absolute_angstrom']
            assert np.isfinite(per[per.method==method][['mae','rmse','bias']]).all().all()
        assert np.isfinite(per[['spearman','pearson','z_mae','z_rmse','top10_precision']]).all().all()
        source=per.groupby(['method','dataset'])[METRICS].mean()
        overall=per.groupby('method')[METRICS].mean()
        folded=per.groupby(['method','dataset','fold'])[METRICS].mean()
        assert set(source.index)==set(index for index in macro.index if index[1]!='all')
        for index,row in source.iterrows(): close(row,macro.loc[index,METRICS])
        for method,row in overall.iterrows(): close(row,macro.loc[(method,'all'),METRICS])
        assert set(folded.index)==set(fold_macro.index)
        for index,row in folded.iterrows(): close(row,fold_macro.loc[index,METRICS])
        original=pd.read_csv(OUT/f'baseline_residue_spotcheck_{protocol}.csv',float_precision='round_trip')
        for sample,rows in original.groupby('sample_id',sort=False):
            assert len(rows) in (415,130)
            for method in METHODS:
                actual=raw_statistics(rows.target_rmsf,rows[method],method not in RELATIVE)
                expected=per[(per.sample_id==sample)&(per.method==method)].iloc[0]
                for metric,value in actual.items(): close(value,expected[metric])
                spot_checks.append(dict(protocol=protocol,sample_id=sample,method=method,
                                        residues=len(rows),all_reported_metrics_recomputed_and_match=True))
        support=per[per.method=='flexppi']
        source_support=support.groupby('dataset').agg(records=('sample_id','size'),residues=('count','sum'),biological_systems=('biological_system','nunique'))
        fold_support=support.groupby(['dataset','fold']).agg(records=('sample_id','size'),residues=('count','sum'))
        protocols[protocol]=dict(records=3063,common_rmsf_residues=961321,methods=7,
                                 per_record_rows=len(per),common_fraction=1.0,
                                 source_support=source_support.reset_index().to_dict('records'),
                                 source_fold_support=fold_support.reset_index().to_dict('records'),
                                 overall=overall.reset_index().replace({np.nan:None}).to_dict('records'),
                                 source=source.reset_index().replace({np.nan:None}).to_dict('records'),
                                 numeric_aggregate_max_tolerance=1e-11,
                                 common_residues_complete_no_failed_record_exclusion=True,
                                 family_fold_sha256=audit['family_fold_sha256'])
        per_by_protocol[protocol]=per
        source.reset_index().to_csv(OUT/f'baseline_{protocol}_source_values_verified.csv',index=False)
        overall.reset_index().to_csv(OUT/f'baseline_{protocol}_overall_values_verified.csv',index=False)
    # Public weights are inferred once and unchanged across the two FlexPPI splits.
    public=list(RELATIVE)+['planet_full','planet_seq']
    for method in public:
        a=per_by_protocol['original'].query('method==@method').set_index('sample_id').sort_index()
        b=per_by_protocol['family'].query('method==@method').set_index('sample_id').sort_index()
        assert a.index.equals(b.index)
        close(a[METRICS],b[METRICS])
    certificates={}
    for method in ['rmsf_net','planet_full','planet_seq']:
        root=BASE/method
        cert=json.loads((root/'run_audit.json').read_text())
        verification=json.loads((root/'verification_audit.json').read_text())
        metadata=json.loads((root/'metadata.json').read_text())
        assert cert['status']=='complete' and cert['failure_count']==0
        assert cert['verified_records']==3063 and cert['verified_residues']==961321
        assert cert['verified_items']==(3063 if method=='rmsf_net' else 4973)
        assert metadata==cert['provenance']['inference_metadata']
        certificates[method]=dict(status=cert['status'],verified_records=cert['verified_records'],
                                  verified_items=cert['verified_items'],verified_residues=cert['verified_residues'],
                                  failure_count=0,completed_utc=cert['completed_utc'],
                                  existing_certificate_sha256=sha(root/'run_audit.json'),
                                  existing_verification_sha256=sha(root/'verification_audit.json'),
                                  metadata=metadata,completion_checker_rerun=False)
    unit=json.loads((OUT/'planet_formal_unit_spotcheck.json').read_text(encoding='utf-8-sig'))
    assert len(unit)==4 and all(row['conversion_exact'] and row['dtype']=='float32' for row in unit)
    family_path=TASK/'family_folds/record_family_folds.csv'
    assert sha(family_path)==protocols['family']['family_fold_sha256']
    b=pd.read_csv(BF/'per_record_metrics.csv',float_precision='round_trip')
    b=b[(b.role=='training')&(b.variant=='backflip2_official')&(b.metric=='rmsf')]
    assert len(b)==3063 and b.sample_id.is_unique and int(b['count'].sum())==961321
    support=per_by_protocol['original'].query('method=="flexppi"').set_index('sample_id').sort_index()
    b=b.set_index('sample_id').sort_index()
    assert b.index.equals(support.index) and b['count'].equals(support['count'])
    bf_overall=b[['mae','rmse','bias','spearman','pearson']].mean().to_dict()
    bf_sources=b.groupby('dataset')[['mae','rmse','bias','spearman','pearson']].mean().reset_index()
    bf_sources.to_csv(OUT/'backflip2_prior_source_values_verified.csv',index=False)
    hashes[str(BF/'per_record_metrics.csv')]=sha(BF/'per_record_metrics.csv')
    hashes[str(family_path)]=sha(family_path)
    original=pd.read_csv(BASE/'additional_original_oof/record_macro_metrics.csv').set_index(['method','dataset'])
    family=pd.read_csv(BASE/'additional_family_oof/record_macro_metrics.csv').set_index(['method','dataset'])
    differences=dict(planet_full_vs_original_Full_mae_percent_reduction=100*(original.loc[('flexppi','all'),'mae']-original.loc[('planet_full','all'),'mae'])/original.loc[('flexppi','all'),'mae'],
                     Full_vs_Intra_original_mae_percent_reduction=100*(original.loc[('intrachain','all'),'mae']-original.loc[('flexppi','all'),'mae'])/original.loc[('intrachain','all'),'mae'],
                     Full_vs_Intra_family_mae_percent_reduction=100*(family.loc[('intrachain','all'),'mae']-family.loc[('flexppi','all'),'mae'])/family.loc[('intrachain','all'),'mae'],
                     planet_full_original_Full_MAE_difference=original.loc[('planet_full','all'),'mae']-original.loc[('flexppi','all'),'mae'],
                     planet_full_original_Full_RMSE_difference=original.loc[('planet_full','all'),'rmse']-original.loc[('flexppi','all'),'rmse'])
    table=[]
    for method in ['planet_full','planet_seq','rmsf_net_pdb']:
        row=original.loc[(method,'all')]
        table.append(dict(method=method,RMSF_MAE_angstrom=None if method in RELATIVE else row.mae,
                          RMSF_Spearman=row.spearman,Hchi_MAE_nats=None,Hchi_Spearman=None,
                          zMAE_unitless=row.z_mae,scope_records=3063,scope_residues=961321,
                          released_weights=True,same_training_labels_as_FlexPPI=False))
    pd.DataFrame(table).to_csv(OUT/'table1_new_baseline_rows.csv',index=False)
    report=dict(status='passed_independent_saved_score_numeric_audit',
                created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                scope='Aggregate and coverage recomputation from all saved per-record metrics plus raw-residue re-scoring of two fixed records in both protocols; no new inference or official completion checker run.',
                remote_files_copied_and_hash_verified=len(hashes)-2,
                source_sha256=hashes,protocols=protocols,certificates=certificates,
                raw_residue_spotchecks=spot_checks,PLANET_formal_unit_spotchecks=unit,
                public_weight_predictions_identical_between_original_and_family_scores=True,
                backflip2_prior_scoring=dict(records=3063,common_residues=961321,overall=bf_overall,
                    caveat='Existing prior corrected scoring supports the same records and residue counts; BackFlip2 is not a method inside the newly generated seven-way intersection CSV.'),
                comparisons=differences,table1_additions=table,
                family_ITrace_single_fold=dict(records=244,fold=0,
                    caveat='All ITrace records are assigned to one family fold, so source-level across-five-fold uncertainty for ITrace is not estimable from this table.'),
                scientific_limits=['RMSF-net PDB uses a reproduced simulated 4A density grid, not measured cryo-EM density; numerical identity to the original Chimera molmap executable has not been established.',
                                   'RMSF-net, ANM and GNM scores have no fitted Angstrom scale, so their MAE/RMSE/bias fields are not absolute RMSF errors.',
                                   'PLANET full uses released encoded features; seq ignores structure features; defaults use 300K and float32 weights without calibration.',
                                   'Public pretraining data and possible ATLAS/mdCATH overlap have not been audited per item, so this is a released-weight comparison, not same-training-corpus retraining.',
                                   'Released baseline predictions remain unchanged when FlexPPI is retrained under family folds; do not call their full-cohort results new family-held-out training results.',
                                   'SHP is not the paper joint discrete chi-state entropy; no Hchi scores are supplied by these methods.',
                                   'Overall means weight records equally; ATLAS has 1929/3063 records and can dominate the pooled rank and MAE. Source strata are required.'],
                completion_checker_rerun=False,certificates_modified=False,DOCX_modified=False,pipeline_modified=False)
    (OUT/'baseline_numeric_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['status'],'records':3063,'residues':961321,
                      'raw_recomputed_comparisons':len(spot_checks),'differences':differences,
                      'backflip2':bf_overall,'table1_additions':table},ensure_ascii=False,indent=2))


if __name__=='__main__': main()
