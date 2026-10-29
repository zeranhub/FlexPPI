"""Local CPU read-only audit of downloaded small completion outputs."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
import numpy as np
import pandas as pd

ROOT = Path('E:/Project/FlexPPI/experiments/_sources')
OUT = ROOT/'progress_readonly_20261012'
SINGLE = ROOT/'single_target'
FAMILY = ROOT/'family_folds'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition,message):
    if not condition:
        raise ValueError(message)


def parse(value):
    return datetime.fromisoformat(value.replace('Z','+00:00'))


def main():
    pipeline = json.loads((OUT/'pipeline_state.json').read_text())
    files = []
    for path in [SINGLE/'evaluation/evaluation_audit.json',SINGLE/'analysis/analysis_audit.json',FAMILY/'evaluation/evaluation_audit.json']:
        report=json.loads(path.read_text())
        require(report['status']=='complete','Incomplete analysis '+str(path))
        corresponding = next(row for row in pipeline['steps'] if any(p.endswith(str(path.relative_to(ROOT)).replace('\\','/')) for p in row.get('completion_files',[])))
        recorded = [value for name,value in corresponding['completion_sha256'].items() if name.endswith(str(path.relative_to(ROOT)).replace('\\','/'))]
        require(recorded==[sha(path)],'Pipeline committed audit hash mismatch '+str(path))
        files.append(dict(path=str(path),sha256=sha(path),status='complete'))
    single=pd.read_csv(SINGLE/'analysis/per_record_metrics.csv')
    paired=pd.read_csv(SINGLE/'analysis/single_vs_dual_per_record.csv')
    family=pd.read_csv(FAMILY/'evaluation/per_record_metrics.csv')
    require(len(single)==12252 and len(paired)==6126 and len(family)==12252,'Unexpected per-record comparison sizes')
    require(not family[['variant','metric','sample_id']].duplicated().any(),'Duplicated family record-target')
    effective=family.groupby(['variant','metric']).agg(records=('sample_id','size'),residues=('count','sum')).reset_index()
    require(set(effective.records)=={3063},'Family arm not complete')
    require(set(effective[effective.metric.eq('rmsf')].residues)=={961321},'RMSF mask count changed')
    require(set(effective[effective.metric.eq('chi_entropy')].residues)=={824988},'Hchi mask count changed')
    effective.to_csv(OUT/'family_effective_counts.csv',index=False)
    numeric=['mae','rmse','bias','spearman','pearson']
    computed=family.groupby(['variant','metric'])[numeric].mean().reset_index()
    frozen=pd.read_csv(FAMILY/'evaluation/overall_record_macro_metrics.csv')
    match=computed.merge(frozen,on=['variant','metric'],suffixes=('_recomputed','_frozen'))
    require(all(np.allclose(match[name+'_recomputed'],match[name+'_frozen'],atol=1e-12) for name in numeric),'Family record macro summary mismatch')
    itrace=family[family.dataset.eq('itrace')]
    require(set(itrace.fold)=={0} and set(itrace.groupby(['variant','metric']).size())=={244},'ITrace is not the frozen single heldout component')
    pairs=family.pivot(index=['sample_id','dataset','metric'],columns='variant',values=['mae','spearman'])
    fp=[]
    for metric in ['rmsf','chi_entropy']:
        subset=pairs[pairs.index.get_level_values('metric')==metric]
        gains=subset[('mae','intrachain')]-subset[('mae','full')]
        fp.append(dict(metric=metric,records=len(subset),mean_intra_minus_full_mae=float(gains.mean()),
            full_lower_mae_records=int((gains>0).sum()),full_lower_mae_fraction=float((gains>0).mean()),
            full_spearman_gain=float((subset[('spearman','full')]-subset[('spearman','intrachain')]).mean())))
    pd.DataFrame(fp).to_csv(OUT/'family_paired_descriptive.csv',index=False)
    sp=[]
    for metric,group in paired.groupby('metric'):
        require(len(group)==3063,'Single/dual comparison incomplete')
        sp.append(dict(metric=metric,records=len(group),record_macro_single_mae=float(group.single_mae.mean()),
            record_macro_dual_mae=float(group.dual_mae.mean()),mean_single_minus_dual_mae=float(group.single_minus_dual_mae.mean()),
            dual_lower_mae_records=int((group.single_minus_dual_mae>0).sum()),
            dual_lower_mae_fraction=float((group.single_minus_dual_mae>0).mean()),
            mean_single_minus_dual_spearman=float(group.single_minus_dual_spearman.mean())))
    pd.DataFrame(sp).to_csv(OUT/'single_dual_paired_descriptive.csv',index=False)
    runtime=[]
    for row in pipeline['steps']:
        if row.get('started_utc') and row.get('finished_utc'):
            runtime.append(dict(step=row['name'],status=row['status'],start_utc=row['started_utc'],finish_utc=row['finished_utc'],
                seconds=(parse(row['finished_utc'])-parse(row['started_utc'])).total_seconds(),source='pipeline state actual UTC timestamps'))
    fq=json.loads((OUT/'queue_all_all.json').read_text())
    sq=json.loads((OUT/'single_target_queue_state.json').read_text())
    for row in fq['jobs']:
        runtime.append(dict(step='family_'+row['variant']+'_fold'+str(row['fold']),status=row['status'],
            seconds=row['finished']-row['started'],source='family queue actual wallclock floats'))
    for row in sq['jobs']:
        runtime.append(dict(step='single_'+row['arm']+'_fold'+str(row['fold']),status=row['status'],
            seconds=(parse(row['finished_utc'])-parse(row['started_utc'])).total_seconds(),source='single queue actual UTC timestamps'))
    family_end=next(row for row in pipeline['steps'] if row['name']=='train_evaluate_family_folds')['finished_utc']
    family_last=max(row['finished'] for row in fq['jobs'])
    runtime.append(dict(step='family_oof_export_after_training',status='complete',seconds=parse(family_end).timestamp()-family_last,
        source='last training job completion to family pipeline export completion; includes orchestration overhead'))
    for variant in ['full','seq']:
        state=json.loads((OUT/('planet_'+variant+'_cpu_queue_state.json')).read_text())
        require(state['status']=='complete','PLANET queue incomplete '+variant)
        require(state['output_dtype_audit']['status']=='passed' and state['output_dtype_audit']['chains']==4973,'PLANET float32 audit incomplete')
        for row in state['children']:
            runtime.append(dict(step='planet_'+variant+'_'+row['phase'],status='complete',seconds=(parse(row['ended_utc'])-parse(row['started_utc'])).total_seconds(),
                source='CPU queue actual phase timestamps; phase polling adds up to30s overhead'))
    pd.DataFrame(runtime).to_csv(OUT/'actual_runtime_steps.csv',index=False)
    coverage=pd.read_csv(SINGLE/'analysis/uncertainty_coverage_five_fold_mean_sample_sd.csv')
    coverage=coverage[coverage.dataset.eq('all')&coverage.variance_kind.eq('model')&coverage.nominal_coverage.eq(.95)]
    coverage.to_csv(OUT/'single_dual_nominal95_coverage.csv',index=False)
    for folder in [SINGLE/'analysis',FAMILY/'evaluation']:
        for path in folder.glob('*.csv'):
            files.append(dict(path=str(path),sha256=sha(path),bytes=path.stat().st_size))
    audit=dict(status='complete_readonly_snapshot',checked_utc=datetime.now(timezone.utc).isoformat(),
        pipeline_snapshot_status=pipeline['status'],pipeline_snapshot_step=pipeline['current_step'],
        trained_fits=20,completed_oof_analyses=['single_target','all_chain_family'],
        family_complete_records_per_arm=3063,valid_rmsf=961321,valid_hchi=824988,
        ITrace_test_fold=0,ITrace_records=244,ITrace_five_fold_SD_identifiable=False,
        PLANET_full_seq_queue_certified=True,RMSF_run_audit_seen=False,
        notes=['RMSF raw prediction count completion is not an independent run_audit certificate.',
            'Primary family macro means weight each record equally; unweighted fold means differ because folds are unequal and source-concentrated.',
            'Single-target comparison uses original five heldout folds; compare original dual, not new-family dual.',
            'Descriptive record wins and fold SD are not independent-seed significance tests.'],sources=files,
        code_sha256=sha(__file__))
    (OUT/'completed_analyses_readonly_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('Single paired',json.dumps(sp,indent=2))
    print('Family paired',json.dumps(fp,indent=2))
    print('Actual runtimes',pd.DataFrame(runtime)[['step','seconds']].to_string(index=False))
    print('Nominal95 coverage',coverage[['variant','metric','observed_coverage_mean','mean_interval_width_mean']].to_string(index=False))


if __name__=='__main__':
    main()
