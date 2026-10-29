"""Independent local CSV audit and compact source tables, no inference."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

ROOT=Path('E:/Project/FlexPPI/experiments/_sources')
OUT=ROOT/'completed_results_20261012'
SUMMARY=OUT/'summary'
AUDIT=OUT/'audit'
SINGLE=ROOT/'single_target'
FAMILY=ROOT/'family_folds'
NUMERIC=['mae','rmse','bias','spearman','pearson']


def require(condition,message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compare(actual,expected,keys,columns):
    joint=actual.merge(expected,on=keys,suffixes=('_recomputed','_saved'),validate='one_to_one')
    require(len(joint)==len(actual)==len(expected),'Summary key membership differs')
    for column in columns:
        require(np.allclose(joint[column+'_recomputed'],joint[column+'_saved'],atol=1e-12,equal_nan=True),
            'Summary numeric mismatch: '+column)


def overall(data,keys):
    means=data.groupby(keys)[NUMERIC].mean().reset_index()
    numbers=data.groupby(keys).agg(records=('sample_id','size'),valid_residues=('count','sum')).reset_index()
    return means.merge(numbers,on=keys,validate='one_to_one')


def main():
    SUMMARY.mkdir(exist_ok=True)
    AUDIT.mkdir(exist_ok=True)
    single=pd.read_csv(SINGLE/'analysis/per_record_metrics.csv')
    paired=pd.read_csv(SINGLE/'analysis/single_vs_dual_per_record.csv')
    family=pd.read_csv(FAMILY/'evaluation/per_record_metrics.csv')
    external=pd.read_csv(FAMILY/'external_2430_evaluation/per_record_metrics.csv')
    require(not single[['variant','metric','sample_id']].duplicated().any(),'Duplicated single analysis record')
    require(not family[['variant','metric','sample_id']].duplicated().any(),'Duplicated family analysis record')
    require(not external[['variant','metric','sample_id']].duplicated().any(),'Duplicated external analysis record')
    require(len(single)==12252 and len(paired)==6126 and len(family)==12252 and len(external)==256,'Wrong complete row counts')
    ssummary=overall(single,['variant','metric'])
    fsummary=overall(family,['variant','metric'])
    esummary=overall(external,['variant','metric'])
    for data in [ssummary,fsummary]:
        require(set(data.records)=={3063},'Incomplete main evaluation records')
        require(set(data[data.metric.eq('rmsf')].valid_residues)=={961321},'Wrong main RMSF coverage')
        require(set(data[data.metric.eq('chi_entropy')].valid_residues)=={824988},'Wrong main Hchi coverage')
    require(set(esummary.records)=={32},'Incomplete external records')
    require(set(esummary[esummary.metric.eq('rmsf')].valid_residues)=={18748},'Wrong external RMSF coverage')
    require(set(esummary[esummary.metric.eq('chi_entropy')].valid_residues)=={16265},'Wrong external Hchi coverage')
    compare(fsummary,pd.read_csv(FAMILY/'evaluation/overall_record_macro_metrics.csv'),['variant','metric'],NUMERIC)
    compare(esummary,pd.read_csv(FAMILY/'external_2430_evaluation/record_macro_metrics.csv'),['variant','metric'],NUMERIC)
    source=overall(family,['variant','dataset','metric'])
    compare(source,pd.read_csv(FAMILY/'evaluation/record_macro_metrics.csv'),['variant','dataset','metric'],NUMERIC)
    # Distinguish equal-record and equal-fold summaries with explicit column
    # names and output files; never attach fold SD to a record-macro mean.
    sfolds=single.groupby(['variant','metric','fold'])[NUMERIC].mean().reset_index()
    ffolds=family.groupby(['variant','metric','fold'])[NUMERIC].mean().reset_index()
    def foldtable(data,keys):
        table=data.groupby(keys)[NUMERIC].agg(['mean','std'])
        table.columns=[name+'_'+scope for name,scope in table.columns]
        return table.reset_index()
    sft=foldtable(sfolds,['variant','metric'])
    fft=foldtable(ffolds,['variant','metric'])
    oldfold=pd.read_csv(SINGLE/'analysis/metrics_five_fold_mean_sample_sd.csv')
    oldfold=oldfold[oldfold.dataset.eq('all')].drop(columns='dataset')
    newfold=pd.read_csv(FAMILY/'evaluation/five_fold_mean_sample_sd.csv')
    newfold=newfold[newfold.dataset.eq('all')].drop(columns='dataset')
    compare(sft,oldfold,['variant','metric'],[name+'_'+scope for name in NUMERIC for scope in ['mean','std']])
    compare(fft,newfold,['variant','metric'],[name+'_'+scope for name in NUMERIC for scope in ['mean','std']])
    lookup=single.set_index(['variant','metric','sample_id'])
    for row in paired.itertuples(index=False):
        single_row=lookup.loc[(row.variant,row.metric,row.sample_id)]
        dual_row=lookup.loc[('full',row.metric,row.sample_id)]
        require(int(row.fold)==int(single_row.fold)==int(dual_row.fold),'Paired fold differs')
        require(row.dataset==single_row.dataset==dual_row.dataset,'Paired dataset differs')
        require(int(row.count)==int(single_row['count'])==int(dual_row['count']),'Paired valid counts differ')
        for metric in ['mae','rmse','spearman','pearson']:
            require(abs(getattr(row,'single_'+metric)-single_row[metric])<1e-12,'Single paired score differs')
            require(abs(getattr(row,'dual_'+metric)-dual_row[metric])<1e-12,'Dual paired score differs')
    records=pd.read_csv(FAMILY/'record_family_folds.csv').set_index('sample_id')
    require(len(records)==3063 and records.index.is_unique,'Family record mapping incomplete')
    for row in family.itertuples(index=False):
        require(row.fold==records.loc[row.sample_id,'fold'] and row.dataset==records.loc[row.sample_id,'dataset'],'Wrong family test fold mapping')
    itrace=records[records.dataset.eq('itrace')]
    require(len(itrace)==244 and set(itrace.fold)=={0} and itrace.component_id.nunique()==1,'ITrace not one component in fold0')
    # Verify no component crosses any of the five test partitions.
    require(records.groupby('component_id').fold.nunique().max()==1,'A family component crosses folds')
    components=json.loads((FAMILY/'family_components.json').read_text())
    relevant=[row for row in components if row['external_records']]
    require(len(relevant)==1 and relevant[0]['external_records']==32 and relevant[0]['training_records']==633,'Wrong strict external component purge')
    external_audit=json.loads((FAMILY/'external_2430_evaluation/evaluation_audit.json').read_text())
    require(external_audit['status']=='complete' and external_audit['training_records']==2430 and external_audit['purged_records']==633,'External protocol not complete2430')
    require(external_audit['external_records']==32 and external_audit['external_family_components']==1,'External scope changed')
    require(external_audit['record_order_exact'] and external_audit['static_and_neural_labels_masks_identical'],'External input/mask audit failed')
    for variant in external.variant.unique():
        for metric,target_index in [('rmsf',0),('chi_entropy',1)]:
            part=external[external.variant.eq(variant)&external.metric.eq(metric)].set_index('sample_id')
            require(set(part.index)=={row['sample_id'] for row in external_audit['inputs']},'External sample identities differ')
            require(all(int(part.loc[row['sample_id'],'count'])==row['expected_valid_counts'][target_index] for row in external_audit['inputs']),
                'External per-record mask counts differ')
    oldnew=fsummary[fsummary.variant.eq('full')].merge(ssummary[ssummary.variant.eq('full')],on=['variant','metric'],suffixes=('_family','_original'))
    oldnew['family_minus_original_mae']=oldnew.mae_family-oldnew.mae_original
    oldnew['family_minus_original_spearman']=oldnew.spearman_family-oldnew.spearman_original
    extpair=[]
    for metric in ['rmsf','chi_entropy']:
        matrix=external[external.metric.eq(metric)].pivot(index='sample_id',columns='variant',values=['mae','spearman'])
        for control in ['intrachain','static_sequence_complex','aa_only']:
            delta=matrix[('mae',control)]-matrix[('mae','full')]
            extpair.append(dict(metric=metric,control=control,records=32,
                mean_control_minus_full_mae=float(delta.mean()),full_lower_mae_records=int((delta>0).sum()),
                full_spearman_gain=float((matrix[('spearman','full')]-matrix[('spearman',control)]).mean())))
    outputs={'single_target_overall_record_macro.csv':ssummary,'single_target_fold_mean_sd.csv':sft,
        'single_target_source_record_macro.csv':overall(single,['variant','dataset','metric']),
        'family_overall_record_macro.csv':fsummary,'family_fold_mean_sd.csv':fft,
        'family_source_record_macro.csv':source,'family_source_fold_macro.csv':family.groupby(['variant','dataset','metric','fold'])[NUMERIC].mean().reset_index(),
        'itrace_single_family_fold0.csv':source[source.dataset.eq('itrace')].assign(test_fold=0,independent_family_components=1,fold_SD_not_identifiable=True),
        'family_vs_original_full_record_macro.csv':oldnew,'external_2430_record_macro.csv':esummary,'external_2430_paired_descriptive.csv':pd.DataFrame(extpair)}
    for name,data in outputs.items():
        data.to_csv(SUMMARY/name,index=False)
    cache=json.loads((AUDIT/'single_dual_cache_comparison.json').read_text(encoding='utf-8-sig'))
    require(cache['status']=='passed_readonly_comparison','Independent cache comparison did not pass')
    inputs=[SINGLE/'analysis/per_record_metrics.csv',SINGLE/'analysis/single_vs_dual_per_record.csv',SINGLE/'analysis/analysis_audit.json',
        SINGLE/'evaluation/evaluation_audit.json',FAMILY/'evaluation/per_record_metrics.csv',FAMILY/'evaluation/evaluation_audit.json',
        FAMILY/'external_2430_evaluation/per_record_metrics.csv',FAMILY/'external_2430_evaluation/evaluation_audit.json',
        FAMILY/'record_family_folds.csv',FAMILY/'family_components.json']
    audit=dict(status='passed_independent_readonly_review',no_training=True,no_inference=True,no_gpu=True,no_existing_audit_modified=True,
        single_dual_all_cached_masks_targets_exact=True,main_records=3063,valid_rmsf=961321,valid_hchi=824988,
        external_records=32,external_train_records=2430,external_purged_records=633,external_family_components=1,
        external_valid_rmsf=18748,external_valid_hchi=16265,ITrace_records=244,ITrace_components=1,ITrace_fold=0,
        ITrace_fold_SD_identifiable=False,family_component_crossfold_count=0,
        source_hashes={str(path):sha(path) for path in inputs},output_hashes={name:sha(SUMMARY/name) for name in outputs},
        script_sha256=sha(__file__),independent_cache_audit_sha256=sha(AUDIT/'single_dual_cache_comparison.json'),
        statistical_interpretation='Recordmacro mean and equal-fold mean/sample SD are separate estimands; no independent-seed significance or independent-family uncertainty claimed')
    (AUDIT/'single_family_csv_review_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('Single overall\n',ssummary.to_string(index=False))
    print('Family overall\n',fsummary.to_string(index=False))
    print('Old/new Full\n',oldnew[['metric','mae_original','mae_family','spearman_original','spearman_family']].to_string(index=False))
    print('External\n',esummary.to_string(index=False))
    print('External pair\n',pd.DataFrame(extpair).to_string(index=False))
    print('Audit PASSED',AUDIT/'single_family_csv_review_audit.json')


if __name__=='__main__':
    main()
