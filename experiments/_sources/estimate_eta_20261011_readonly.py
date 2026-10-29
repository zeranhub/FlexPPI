"""Read-only revised ETA from committed queue timing; never rewrites certificates."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import re
import statistics

ROOT=Path('/root/flexppi_completion_20261010')
OLD=Path('/root/flexppi_revision_20261002/corrected_training_v1')
TZ=datetime.timezone(datetime.timedelta(hours=8))


def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def date(value):return datetime.datetime.fromisoformat(value.replace('Z','+00:00'))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--local-utc',required=True)
    args=parser.parse_args()
    local=date(args.local_utc)
    server=datetime.datetime.now(datetime.timezone.utc)
    skew=(local-server).total_seconds()
    family=read(ROOT/'family_folds/queue_all_all.json')
    single=read(ROOT/'single_target/single_target_queue_state.json')
    pipeline=read(ROOT/'pipeline_state.json')
    old=read(OLD/'queue_state.json')
    done=[row for row in single['jobs'] if row['status'].startswith('complete')]
    single_durations=[dict(arm=row['arm'],fold=row['fold'],started_utc=row['started_utc'],finished_utc=row['finished_utc'],
        elapsed_seconds=(date(row['finished_utc'])-date(row['started_utc'])).total_seconds()) for row in done]
    full_history=[row['finished']-row['started'] for row in old['jobs'] if row['variant']=='full' and row['status']=='complete']
    intra_history=[row['finished']-row['started'] for row in old['jobs'] if row['variant']=='intrachain' and row['status']=='complete']
    family_done=[row for row in family['jobs'] if row['status'] in ['complete','already_complete']]
    current=next(row for row in family['jobs'] if row['status']=='running')
    family_times=[dict(variant=row['variant'],fold=row['fold'],started_unix=row['started'],finished_unix=row['finished'],
        elapsed_seconds=row['finished']-row['started']) for row in family_done]
    log=ROOT/'family_folds/logs'/f"{current['variant']}_fold{current['fold']}.log"
    lines=log.read_text(errors='replace').splitlines()
    epochs=[tuple(map(int,match.groups())) for line in lines if (match:=re.search(r'^epoch=(\d+) step=(\d+)',line))]
    bad=[line for line in lines if 'Traceback (most recent call last)' in line or 'CUDA out of memory' in line or
         (line.startswith('epoch=') and re.search(r"'loss':\s*(?:nan|[+-]?inf)\b",line,re.I))]
    remaining_full=5-sum(row['variant']=='full' for row in family_done)
    remaining_intra=5-sum(row['variant']=='intrachain' for row in family_done)
    elapsed_current=server.timestamp()-current['started']
    scale=statistics.mean(row['elapsed_seconds'] for row in single_durations)/statistics.mean(full_history)
    point=(remaining_full*statistics.mean(full_history)+remaining_intra*statistics.mean(intra_history))*scale-elapsed_current
    bases={}
    for name,total in [('rmsf_net',3063),('planet_full',4973),('planet_seq',4973)]:
        counts={}
        for path in (ROOT/'baselines'/name).glob('*.json'):
            value=read(path)
            if ('sample_id' in value or 'input_id' in value) and 'status' in value:counts[value['status']]=counts.get(value['status'],0)+1
        audit_path=ROOT/'baselines'/name/'run_audit.json'
        audit=read(audit_path) if audit_path.exists() else None
        bases[name]=dict(observed_state_counts=counts,expected_items=total,
            certificate_exists=audit is not None,certificate_sha256=sha(audit_path) if audit else None,
            certificate_summary={key:audit.get(key) for key in ['status','model','verified_items','verified_records','verified_residues','failure_count','completed_utc']} if audit else None,
            inference_done=counts.get('complete')==total and counts.get('failed',0)==0,
            certificate_read_only_no_reverification=True)
    control=read(ROOT/'pipeline_control.json')
    remaining_family_hours=[14.0,18.0]
    postprocess_hours=[0.5,2.0]
    overall_hours=[14.5,20.0]
    eta=[(local+datetime.timedelta(hours=h)).astimezone(TZ).isoformat() for h in overall_hours]
    report=dict(status='read_only_eta_updated_from_actual_completed_work',local_reference_utc=local.isoformat(),
        local_reference_beijing=local.astimezone(TZ).isoformat(),server_snapshot_utc=server.isoformat(),
        local_ahead_of_server_seconds=skew,clock_policy='Durations use differences within server timing; calendar ETA uses local clock, which is approximately12min ahead. Local reference supplied immediately before SSH execution.',
        evidence=dict(single_queue_sha256=sha(ROOT/'single_target/single_target_queue_state.json'),
            family_queue_sha256=sha(ROOT/'family_folds/queue_all_all.json'),old_queue_sha256=sha(OLD/'queue_state.json'),
            control_sha256=sha(ROOT/'pipeline_control.json'),current_log=str(log)),
        pipeline=dict(status=pipeline['status'],current_step=pipeline['current_step'],pid=pipeline['pid']),
        single=dict(completed_fits=len(done),planned_fits=10,actual_durations=single_durations,
            elapsed_seconds_mean=statistics.mean(row['elapsed_seconds'] for row in single_durations),
            elapsed_seconds_range=[min(row['elapsed_seconds'] for row in single_durations),max(row['elapsed_seconds'] for row in single_durations)],
            oof_export_complete=read(ROOT/'single_target/evaluation/evaluation_audit.json').get('status')=='complete',
            paired_analysis_complete=read(ROOT/'single_target/analysis/analysis_audit.json').get('status')=='complete'),
        family=dict(completed_fits=len(family_done),planned_fits=10,committed_actual_durations=family_times,
            current_job=current,current_epoch_observed=epochs[-1][0] if epochs else None,current_step_observed=epochs[-1][1] if epochs else None,
            current_fit_elapsed_seconds=elapsed_current,current_log_tail=lines[-4:],observed_error_lines=bad,
            remaining_full_fits=remaining_full,remaining_intrachain_fits=remaining_intra,
            historical_full_seconds_range=[min(full_history),max(full_history)],historical_intra_seconds_range=[min(intra_history),max(intra_history)],
            historical_full_seconds_mean=statistics.mean(full_history),historical_intra_seconds_mean=statistics.mean(intra_history),
            current_single_to_historical_full_speed_ratio=scale,empirical_analogue_remaining_seconds_point=point,
            remaining_hours_working_range=remaining_family_hours,
            limitation='New family fold0 heldouts all244 ITrace and462 of870 multichain records; its60.77min fit has a different training source/length mixture. Do not multiply that lone fit time by all remaining9. No completed new intrachain fit is available yet; use prior matching40epoch fits with current10-single timing as reference.'),
        baselines=bases,
        remaining_postprocess=dict(hours_allowance=postprocess_hours,new_walltime_not_yet_measured=True,
            work=['family OOF predictions','strict2430 external prediction','RMSF-net standalone verifier after existing inference completion','CPU original/family matched scoring'],
            no_remaining_planet_or_rmsf_inference_wait=all(row['inference_done'] for row in bases.values())),
        computational_eta=dict(remaining_hours_working_range=overall_hours,eta_beijing_range=eta,
            point_hours_with1h_postprocess=point/3600+1.0,
            point_eta_beijing=(local+datetime.timedelta(seconds=point+3600)).astimezone(TZ).isoformat(),
            interpretation_zh='计算更可能在今晚深夜至明日凌晨完成；不能承诺今天白天全部完成。14.5–20小时是工作估计，家族训练只有首项新耗时，需在首个链内模型完成后更新。',
            conditions=['same isolated GPU queue continues without failures or resource contention','9 remaining40epoch family fits still run serially','12 newKIX production runs remain cancelled','no additional retraining requested']),
        manuscript_eta=dict(additional_active_work_hours=[2.0,4.0],
            immediate_followup_calendar_range=[(local+datetime.timedelta(hours=overall_hours[0]+2)).astimezone(TZ).isoformat(),(local+datetime.timedelta(hours=overall_hours[1]+4)).astimezone(TZ).isoformat()],
            calendar_range_only_if_immediate_active_followup=True,heartbeat_created=False,automatic_unattended_manuscript_update=False,
            interpretation_zh='服务器计算结束后，还需约2–4小时实际整理结果、更新图表正文及重新渲染核验；当前没有heartbeat，不保证无人值守自动回填。'),
        md_cancel_still_in_control=control['md_reuse_policy']['allow_new_kix_production'] is False,
        read_only=True,no_training_queue_threads_processes_changed=True,no_planet_completion_certificate_rewritten=True)
    output=ROOT/'estimate_20261011_0845_readonly.json'
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({key:report[key] for key in ['local_reference_beijing','server_snapshot_utc','local_ahead_of_server_seconds','family','baselines','computational_eta','manuscript_eta']},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
