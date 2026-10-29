"""Resume only after a guarded evaluation-adapter repair; no training restart."""
from pathlib import Path
import datetime
import fcntl
import hashlib
import json
import os
import subprocess
import time

ROOT=Path('/root/flexppi_completion_20261010')
OLD='03f865c7c477f05a67d747325af22b251965d69141a2cc3c72ed5fb24564a07a'
NEW='0d28fb4f8778f31441ada7b812cab48fe7f2ae0b9845ec4873c451e78eff291b'
CONTROL='837494a7ecbdc8cb87d464b34d1e8bdfdd5bb77cbe8dd525dbc8896f8d2661a5'
PIPELINE='009ee9c32e73ea65c52a7364c83f3fe26c70bbba1b731cbbbf7fe9bfec481a69'


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''):digest.update(block)
    return digest.hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,value):
    path=Path(path);temp=Path(str(path)+'.repair.tmp')
    temp.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n');os.replace(temp,path)


def main():
    with (ROOT/'gpu_queue_restart.lock').open('a+') as restart:
        fcntl.flock(restart,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for folder in Path('/proc').iterdir():
            if not folder.name.isdigit():continue
            try:cmd=(folder/'cmdline').read_bytes().decode().rstrip('\0').split('\0')
            except (FileNotFoundError,PermissionError):continue
            assert str(ROOT/'run_completion_pipeline.py') not in cmd,'A live pipeline already owns the queue'
            assert str(ROOT/'family_folds/evaluate_external_2430.py') not in cmd,'External evaluation already running'
        assert sha(ROOT/'external_adapter_repair_20261012/evaluate_external_2430.before.py')==OLD
        assert sha(ROOT/'family_folds/evaluate_external_2430.py')==NEW
        assert sha(ROOT/'run_completion_pipeline.py')==PIPELINE
        assert sha(ROOT/'pipeline_control.json')==CONTROL
        control=read(ROOT/'pipeline_control.json')
        assert control['md_reuse_policy']['allow_new_kix_production'] is False and len(control['cancel_steps'])==5
        with (ROOT/'gpu_queue.lock').open('r+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            fcntl.flock(lock,fcntl.LOCK_UN)
        state=read(ROOT/'pipeline_state.json')
        assert state['status']=='failed' and state['current_step']=='evaluate_external_2430'
        failed=next(row for row in state['steps'] if row['name']=='evaluate_external_2430')
        assert failed['status']=='failed' and failed['exit_code']==1
        for row in state['steps']:
            if row['status']=='complete':
                for path,expected in row.get('completion_sha256',{}).items():assert sha(path)==expected,'Prior committed completion artifact changed'
        completed_weights=[]
        family=read(ROOT/'family_folds/queue_all_all.json')
        single=read(ROOT/'single_target/single_target_queue_state.json')
        assert family['status']=='complete' and len(family['jobs'])==10
        assert single['status']=='complete' and len(single['jobs'])==10
        for row in family['jobs']:
            assert row['status'] in ['complete','already_complete']
            path=ROOT/f"family_folds/runs/{row['variant']}/fold_{row['fold']}/epoch_040.pt"
            assert sha(path)==row['checkpoint_sha256']
            completed_weights.append(dict(path=str(path),sha256=row['checkpoint_sha256']))
        for row in single['jobs']:
            assert row['status'].startswith('complete') and sha(row['final_checkpoint'])==row['checkpoint_sha256']
            completed_weights.append(dict(path=row['final_checkpoint'],sha256=row['checkpoint_sha256']))
        checkpoints={'full':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',
                     'intrachain':'6b069d9d3521e95f5f2d63a3f19d64fb4071862efa4682473387b0ca45ff73c3'}
        for variant,expected in checkpoints.items():
            assert sha(ROOT/f'family_folds/runs/{variant}/fold_0/epoch_040.pt')==expected
        certs={name:sha(ROOT/f'baselines/planet_{name}/run_audit.json') for name in ['full','seq']}
        history=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            reason='Evaluation-only adapter supplied three source label/mask channels to a two-output count array. Frozen metrics map rmsf→0 and chi_entropy→1 explicitly; SASA excluded by name. No labels, structures, training code/configuration or weights changed.',
            adapter_before_sha256=OLD,adapter_after_sha256=NEW,source_pipeline_sha256=PIPELINE,control_sha256=CONTROL,
            previous_failed_step=failed,completed_checkpoints_verified=completed_weights,
            reuse_existing_purged_family_fold0=True,expected_external_records=32,expected_valid_rmsf=18748,expected_valid_chi_entropy=16265,
            planet_completion_certificates_before=certs,training_not_restarted=True,md_cancellation_preserved=True)
        write(ROOT/'external_adapter_repair_20261012/repair_provenance.json',history)
        state['repair_provenance']=str(ROOT/'external_adapter_repair_20261012/repair_provenance.json')
        state.pop('error',None)
        write(ROOT/'pipeline_state.json',state)
        log=ROOT/'pipeline_external_adapter_repair_20261012.log'
        command=['nohup',str(ROOT/'env/bin/python'),'-u',str(ROOT/'run_completion_pipeline.py')]
        with log.open('ab',buffering=0) as handle:
            child=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=handle,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
        deadline=time.monotonic()+15
        while time.monotonic()<deadline:
            latest=read(ROOT/'pipeline_state.json')
            if latest.get('pid')==child.pid and latest['status'] in ['running','computations_complete_pending_manuscript_update']:
                break
            if child.poll() is not None:raise RuntimeError('Resumed parent exited: '+str(child.returncode))
            time.sleep(.2)
        else:raise RuntimeError('Resumed parent did not commit state')
        history.update(status='resumed_after_adapter_repair',parent_pid=child.pid,parent_command=command,log=str(log),
            resumed_state=latest['status'],current_step=latest.get('current_step'))
        write(ROOT/'external_adapter_repair_20261012/repair_provenance.json',history)
        print(json.dumps({key:history[key] for key in ['status','parent_pid','adapter_before_sha256','adapter_after_sha256','training_not_restarted','md_cancellation_preserved','expected_external_records','expected_valid_rmsf','expected_valid_chi_entropy','current_step']},indent=2))


if __name__=='__main__':main()
