"""Avoid duplicate RMSF-net workers, then check cached run and independently verify.

The old wrapper is deliberately untouched. This scheduler takes an inherited
exclusive lock and waits for every original wrapper targeting this output.
No --retry-failed is passed: existing failed records must fail verification.
"""
import argparse
import fcntl
import os
from pathlib import Path
import subprocess
import time
from check_additional_baseline_freeze import BASE, ROOT, COHORT_HASH, require, preflight, write_json


def active_workers(wrapper, output):
    result = []
    for folder in Path('/proc').iterdir():
        if not folder.name.isdigit() or int(folder.name) == os.getpid():
            continue
        try:
            command = (folder / 'cmdline').read_bytes().split(b'\0')
            args = [part.decode(errors='replace') for part in command if part]
            if str(wrapper) not in args or '--output' not in args:
                continue
            target = args[args.index('--output')+1]
            if Path(target) == output:
                result.append({'pid':int(folder.name),'command':args})
        except (OSError,ValueError,IndexError):
            continue
    return sorted(result,key=lambda row:row['pid'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--timeout-seconds',type=float,default=86400)
    parser.add_argument('--poll-seconds',type=float,default=20)
    args = parser.parse_args()
    require(0 < args.timeout_seconds <= 86400, 'Timeout must be positive and <=24h')
    require(0 < args.poll_seconds <= 20, 'Polling interval must be <=20s')
    output, wrapper = BASE / 'rmsf_net', BASE / 'run_rmsf_net_3063.py'
    audit_path, lock_path = BASE / 'rmsf_net_resume_state.json', BASE / 'rmsf_net_resume.lock'
    started, observed = time.monotonic(), []
    env = dict(os.environ,CUDA_VISIBLE_DEVICES='',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
    with lock_path.open('a+') as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit('A guarded RMSF-net resume is already active; no duplicate launched')
        report = dict(status='waiting',pid=os.getpid(),output=str(output),no_gpu_used=True,observed_original_workers=[])
        write_json(audit_path,report)
        try:
            while True:
                workers = active_workers(wrapper,output)
                if not workers:
                    break
                for row in workers:
                    if row['pid'] not in [x['pid'] for x in observed]:
                        observed.append(row)
                        print('Waiting for original RMSF-net CPU PID',row['pid'],flush=True)
                report.update(observed_original_workers=observed,current_worker_pids=[w['pid'] for w in workers])
                write_json(audit_path,report)
                require(time.monotonic()-started < args.timeout_seconds, 'Timed out waiting for original RMSF-net process')
                time.sleep(min(args.poll_seconds,args.timeout_seconds-(time.monotonic()-started)))
            report['preflight'] = preflight()
            require(not active_workers(wrapper,output), 'Another unguarded original RMSF-net worker appeared during preflight')
            report.update(status='cache_resume',observed_original_workers=observed,current_worker_pids=[])
            write_json(audit_path,report)
            command = [str(ROOT/'env/bin/python'),'-u',str(wrapper),'--cohort',str(BASE/'cohort'),
                '--output',str(output),'--source-root',str(BASE/'RMSF-net'),'--weights',str(BASE/'RMSF-net/model_cpu/sim_model.pth'),
                '--threads','2','--batch-size','2']
            report['cache_resume_command'] = command
            write_json(audit_path,report)
            # A surviving child keeps this lock if the parent is interrupted.
            result = subprocess.run(command,env=env,stdin=subprocess.DEVNULL,pass_fds=(lock.fileno(),))
            require(result.returncode == 0, 'RMSF-net cache resume exited with '+str(result.returncode))
            verifier = [str(ROOT/'cpu_env/bin/python'),'-u',str(BASE/'verify_additional_baselines.py'),
                '--cohort',str(BASE/'cohort'),'--output',str(output),'--model','rmsf-net',
                '--expected-cohort-freeze-sha256',COHORT_HASH,'--weights',str(BASE/'RMSF-net/model_cpu/sim_model.pth')]
            report.update(status='verification',verifier_command=verifier)
            write_json(audit_path,report)
            result = subprocess.run(verifier,env=env,stdin=subprocess.DEVNULL,pass_fds=(lock.fileno(),))
            require(result.returncode == 0, 'Independent RMSF-net verification failed')
            report.update(status='complete',seconds=time.monotonic()-started)
            write_json(audit_path,report)
        except Exception as error:
            report.update(status='failed',error=type(error).__name__+': '+str(error),seconds=time.monotonic()-started)
            write_json(audit_path,report)
            print(report['error'],flush=True)
            raise SystemExit(1)


if __name__ == '__main__':
    main()
