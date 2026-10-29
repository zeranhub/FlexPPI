"""Cancel only unstarted KIX work; safely replace the monitoring parent.

Never signals the existing single-target queue, its training child, RMSF-net,
or either PLANET queue. Exact PID/start-time/argv identities guard every action.
"""
from pathlib import Path
import datetime
import fcntl
import hashlib
import json
import os
import signal
import subprocess
import time

ROOT=Path('/root/flexppi_completion_20261010')
BEFORE_SHA='6e3c82f5354cb4d0022fbe3ed312c2f6d6ce8e10fcdfdeb4838eef3fa44e07a4'
AFTER_SHA='009ee9c32e73ea65c52a7364c83f3fe26c70bbba1b731cbbbf7fe9bfec481a69'
CANCEL=['md_triplicate_production','analyze_md_triplicates','quality_md_rep01','quality_md_rep02','quality_md_rep03']


def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,value):
    path=Path(path)
    temp=Path(str(path)+'.cancellation.tmp')
    temp.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
    os.replace(temp,path)


def identity(pid):
    path=Path('/proc')/str(pid)
    try:
        value=(path/'stat').read_text();parts=value[value.rfind(')')+2:].split()
        return dict(pid=int(pid),start_ticks=int(parts[19]),state=parts[0],
                    command=(path/'cmdline').read_bytes().decode().rstrip('\0').split('\0'))
    except FileNotFoundError:return None


def same(expected):
    actual=identity(expected['pid'])
    return actual is not None and actual['state']!='Z' and all(actual[key]==expected[key] for key in ['pid','start_ticks','command'])


def processes():
    items=[]
    for path in Path('/proc').iterdir():
        if not path.name.isdigit():continue
        item=identity(int(path.name))
        if item and item['state']!='Z' and str(ROOT) in ' '.join(item['command']):items.append(item)
    return items


def main():
    audit_dir=ROOT/'runtime_cancellation_audit_20261010'
    audit_dir.mkdir(exist_ok=True)
    with (ROOT/'gpu_queue_restart.lock').open('a+') as restart_lock:
        fcntl.flock(restart_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        source=ROOT/'run_completion_pipeline.py'
        assert sha(source)==BEFORE_SHA,'Parent source differs from inspected original'
        ready=ROOT/'run_completion_pipeline.cancel_ready.py'
        assert sha(ready)==AFTER_SHA,'Replacement differs from reviewed source'
        state=read(ROOT/'pipeline_state.json')
        assert state['pid']==71135 and state['status']=='running' and state['current_step']=='train_single_targets'
        parent=identity(71135)
        assert parent and parent['command']==[str(ROOT/'env/bin/python'),'-u',str(source)]
        running=[row for row in state['steps'] if row.get('status')=='running']
        assert len(running)==1 and running[0]['name']=='train_single_targets'
        owned=processes()
        queues=[item for item in owned if 'single_target/run_single_target_queue.py' in ' '.join(item['command'])]
        assert len(queues)==1 and queues[0]['pid']==71188
        queue=queues[0]
        assert queue['command']==running[0]['command']
        trainees=[item for item in owned if '-m flexppi train --config' in ' '.join(item['command'])]
        assert len(trainees)==1
        training=trainees[0]
        assert not any('queue_triplicate.py' in ' '.join(item['command']) and '--verify-only' not in item['command'] for item in owned)
        assert not any(row['name'] in CANCEL and row.get('status') in ['running','complete'] for row in state['steps'])
        assert not list((ROOT/'kix_md_triplicate_rebuild_20261010/runs').glob('*/rep*/DONE.json'))
        protected=[item for item in owned if item['pid'] in [queue['pid'],training['pid'],70967,73982,73983,73988,73989]]
        assert len(protected)==7,'Expected protected CPU/training process set changed'
        (audit_dir/'run_completion_pipeline.before.py').write_bytes(source.read_bytes())
        write(audit_dir/'pipeline_state.before.json',state)
        proof_before={row['name']:row.get('completion_sha256') for row in state['steps'] if row['status']=='complete'}
        control=dict(protocol='user_md_reuse_priority_20261010_v1',updated_utc=now(),
            reason='User requested reuse of existing local/hmu-gpu KIX MD, one trajectory per condition is sufficient, and manuscript priority; cancel all unstarted new triplicate simulations.',
            cancel_steps=CANCEL,preserved_running_steps={'train_single_targets':queue},
            md_reuse_policy=dict(strategy='reuse_existing_audited_single_trajectory',priority='manuscript',
                allow_new_kix_production=False,additional_kix_production_runs=0,cancelled_planned_new_runs=12,
                independent_replicates_required=1,never_claim_triplicate_uncertainty_from_single_trajectory=True,
                existing_md_search_hosts=['local','hmu-gpu','hit-gpu'],
                provenance_requirement='Only actual audited production trajectories; short CPU/CUDA technical smoke does not count as biological MD validation'),
            handoff=dict(old_parent=parent,protected_processes=protected,old_source_sha256=BEFORE_SHA,new_source_sha256=AFTER_SHA,
                         restart_policy='Signal exactly the old monitoring parent PID; preserve children; new parent reacquires the existing gpu_queue.lock path and monitors original queue without starting training'))
        write(ROOT/'pipeline_control.json',control)
        replacement=source.with_suffix('.cancellation.tmp.py')
        replacement.write_bytes(ready.read_bytes());os.replace(replacement,source)
        assert same(parent) and all(same(item) for item in protected),'Process identity changed before handoff'
        os.kill(parent['pid'],signal.SIGTERM)
        deadline=time.monotonic()+10
        while same(parent) and time.monotonic()<deadline:time.sleep(.1)
        assert not same(parent),'Old monitor did not stop; no second monitor started'
        assert all(same(item) for item in protected),'A protected child changed during parent-only restart'
        log=ROOT/'pipeline_md_cancel_resumed_20261010.log'
        command=['nohup',str(ROOT/'env/bin/python'),'-u',str(source)]
        with log.open('ab',buffering=0) as handle:
            new=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=handle,stderr=subprocess.STDOUT,
                                 start_new_session=True,close_fds=True)
        deadline=time.monotonic()+20
        effective=None
        while time.monotonic()<deadline:
            latest=read(ROOT/'pipeline_state.json')
            job=next((row for row in latest['steps'] if row['name']=='train_single_targets'),{})
            if latest.get('pid')==new.pid and job.get('preserved_process',{}).get('start_ticks')==queue['start_ticks']:
                effective=latest;break
            if new.poll() is not None:raise RuntimeError('New monitor exited: '+str(new.returncode))
            time.sleep(.2)
        assert effective is not None,'New monitor did not commit runtime handoff'
        assert effective['status']=='running' and effective['current_step']=='train_single_targets'
        assert effective['runtime_control_sha256']==sha(ROOT/'pipeline_control.json')
        assert effective['cancelled_pending_steps']==CANCEL
        assert all(same(item) for item in protected),'A protected process changed after new monitor started'
        assert proof_before=={row['name']:row.get('completion_sha256') for row in effective['steps'] if row['status']=='complete'}
        with (ROOT/'gpu_queue.lock').open('a+') as test_lock:
            try:
                fcntl.flock(test_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:lock_confirmed=True
            else:
                fcntl.flock(test_lock,fcntl.LOCK_UN);raise RuntimeError('New monitor has not reacquired GPU queue lock')
        report=dict(status='runtime_effective_pending_kix_md_cancelled',completed_utc=now(),
            cancellation_scope=CANCEL,cancelled_unstarted_production_runs=12,production_runs_were_started=False,
            old_parent=parent,new_parent=identity(new.pid),new_parent_command=command,
            original_single_queue=queue,current_training_process=training,protected_processes_before_and_after=protected,
            all_protected_pid_starttime_argv_unchanged=True,only_signal_sent=dict(pid=parent['pid'],signal='SIGTERM'),
            same_gpu_queue_lock_path=str(ROOT/'gpu_queue.lock'),new_parent_exclusive_lock_verified=lock_confirmed,
            completed_step_audits_preserved=True,control_sha256=sha(ROOT/'pipeline_control.json'),
            pipeline_before_sha256=BEFORE_SHA,pipeline_after_sha256=sha(source),restart_log=str(log),
            effective_runtime_state_sha256=sha(ROOT/'pipeline_state.json'),md_reuse_policy=control['md_reuse_policy'],
            no_training_restarted=True,no_cpu_baseline_process_touched=True,
            boundary='The monitoring parent was replaced; existing training/CPU children continue. The persistent lock file was retained, and its exclusive ownership was reacquired by the new monitor after the original monitor released it.')
        write(ROOT/'pipeline_md_cancellation_audit_20261010.json',report)
        print(json.dumps(report,indent=2))


if __name__=='__main__':main()
