"""Read progress without loading weights, opening CUDA, or changing queues."""
from pathlib import Path
import collections,datetime,json,re,subprocess
ROOT=Path(__file__).resolve().parent
def read_json(path):
    if not path.exists():return None
    try:return json.loads(path.read_text())
    except Exception as error:return {'read_error':str(error),'path':str(path)}
def tail(path,byte_limit=65536):
    if not path.exists():return []
    with path.open('rb') as handle:
        handle.seek(max(0,path.stat().st_size-byte_limit))
        return handle.read().decode('utf8','replace').splitlines()[-4:]
def predictions(directory):
    counts=collections.Counter()
    seconds=[]
    for path in directory.glob('*.json'):
        item=read_json(path)
        if item and ('sample_id' in item or 'input_id' in item) and 'status' in item:
            counts[item['status']]+=1
            if item['status']=='complete':seconds.append(item.get('seconds',item.get('total_seconds',0)))
    return {'observed_record_or_chain_states':dict(counts),'observed_prediction_seconds_sum':sum(seconds),
            'completion_audit':read_json(directory/'run_audit.json')}
single_queue=read_json(ROOT/'single_target/single_target_queue_state.json')
current_single_job=next((job for job in (single_queue or {}).get('jobs',[]) if job.get('status')=='running'),None)
current_single_log=(ROOT/'single_target/logs'/f"{current_single_job['arm']}_fold{current_single_job['fold']}.log") if current_single_job else None
snapshot={'snapshot_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'pipeline':read_json(ROOT/'pipeline_state.json'),
 'single_target':{'queue':single_queue,'current_job':current_single_job,
                  'final_checkpoint_files_observed':len(list((ROOT/'single_target/runs').glob('*/fold_*/epoch_040.pt'))),
                  'current_log_path':str(current_single_log) if current_single_log else None,
                  'current_log_tail':tail(current_single_log) if current_single_log else []},
 'family_final_checkpoint_files_observed':len(list((ROOT/'family_folds/runs').glob('*/fold_*/epoch_040.pt'))),
 'baselines':{name:predictions(ROOT/'baselines'/name) for name in ['rmsf_net','planet_full','planet_seq']},
 'planet_cpu_states':{name:read_json(ROOT/f'baselines/planet_{name}_cpu_queue_state.json') for name in ['full','seq']},
 'md_done_files_observed':len(list((ROOT/'kix_md_triplicate_rebuild_20261010/runs').glob('*/rep*/DONE.json'))),
 'md_worker':read_json(ROOT/'kix_md_triplicate_rebuild_20261010/extension_worker_apo_myb_mll_both.json'),
 'caution':'Observed file or item counts are progress indicators; only independent completion audits certify full scientific outputs. Manuscript update remains required after computations.'}
try:
    snapshot['gpu_read_only']=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used','--format=csv,noheader'],text=True).strip()
    lines=subprocess.check_output(['ps','-eo','pid,ppid,etimes,args'],text=True).splitlines()
    snapshot['owned_task_processes']=[line.strip() for line in lines if str(ROOT) in line and 'collect_remote_progress' not in line and 'bash -c' not in line]
except Exception as error:snapshot['process_query_error']=str(error)
temporary=ROOT/'completion_progress.latest.tmp.json'
temporary.write_text(json.dumps(snapshot,indent=2),encoding='utf8')
temporary.replace(ROOT/'completion_progress.latest.json')
compact={k:v for k,v in snapshot.items() if k in ['snapshot_utc','family_final_checkpoint_files_observed','md_done_files_observed','gpu_read_only']}
compact['pipeline_status']=snapshot['pipeline']['status'] if snapshot['pipeline'] else None
compact['current_step']=snapshot['pipeline'].get('current_step') if snapshot['pipeline'] else None
compact['single_target_current_log_tail']=snapshot['single_target']['current_log_tail']
compact['single_target_current_job']=snapshot['single_target']['current_job']
compact['single_target_final_checkpoint_files_observed']=snapshot['single_target']['final_checkpoint_files_observed']
compact['baselines']={name:{'observed_states':data['observed_record_or_chain_states'],
    'has_completion_audit':data['completion_audit'] is not None} for name,data in snapshot['baselines'].items()}
compact['planet_cpu_states']={name:({k:state.get(k) for k in ['status','phase','parent_pid','child_pid']} if state else None)
    for name,state in snapshot['planet_cpu_states'].items()}
print(json.dumps(compact,indent=2))
