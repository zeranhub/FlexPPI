"""CPU-only negative tests on independent fixtures; no inference is launched."""
import json
from pathlib import Path
import subprocess
import sys
import check_additional_baseline_freeze as guard
from wait_resume_rmsf_net_cpu import active_workers


def main():
    folder = guard.BASE / 'guard_negative_tests_20261010'
    folder.mkdir(exist_ok=True)
    output = folder / 'fake_output'
    output.mkdir(exist_ok=True)
    guard.write_json(output/'run_audit.json',{'status':'complete','note':'Deliberately stale fixture; must never bypass failed/running/missing queue'})
    cases = [('failed_state',{'status':'failed','error':'fixture failure'}),
        ('failed_item',{'status':'running','observed_items':{'failed':1}}),
        ('missing_state',None)]
    results = []
    for name, state in cases:
        state_path = folder / (name+'.json')
        if state is not None:
            guard.write_json(state_path,state)
        elif state_path.exists():
            state_path.unlink()
        command = [sys.executable,str(guard.BASE/'wait_additional_cpu_queue.py'),'--state',str(state_path),
            '--output',str(output),'--model','planet-full','--timeout-seconds','0.2','--poll-seconds','0.1',
            '--wait-audit',str(folder/(name+'_wait_audit.json'))]
        run = subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        guard.require(run.returncode == 1,'Negative queue fixture was accepted: '+name)
        report = guard.read_json(folder/(name+'_wait_audit.json'))
        guard.require(report['status'] == 'failed','Negative fixture wrote complete audit: '+name)
        results.append(dict(case=name,exit_code=run.returncode,error=report['error']))
    key = 'cohort/cohort_freeze.json'
    original = guard.HASHES[key]
    guard.HASHES[key] = '0'*64
    try:
        guard.preflight()
    except ValueError as error:
        results.append(dict(case='frozen_hash_changed',rejected=True,error=str(error)))
    else:
        raise ValueError('Changed dispatched hash was accepted')
    finally:
        guard.HASHES[key] = original
    contract_path = guard.BASE/'additional_cpu_execution_contract.json'
    contract = guard.read_json(contract_path)
    for name, expected in contract['files'].items():
        guard.require(guard.sha(guard.BASE/name) == expected,'Real deployed dispatch code is not frozen: '+name)
    workers = active_workers(guard.BASE/'run_rmsf_net_3063.py',guard.BASE/'rmsf_net')
    result = dict(status='passed',no_gpu_used=True,no_inference_performed=True,negative_checks=results,
        execution_contract_sha256=guard.sha(contract_path),frozen_files_checked=len(contract['files']),
        original_rmsf_workers_detected_read_only=workers)
    guard.write_json(folder/'test_audit.json',result)
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
