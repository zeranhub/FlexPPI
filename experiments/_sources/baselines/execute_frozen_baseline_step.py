"""CPU guard around each dispatched waiting, verification, or scoring step."""
import argparse
import os
from pathlib import Path
import subprocess
from check_additional_baseline_freeze import BASE, COHORT_HASH, HASHES, preflight, require, read_json, sha, write_json
from wait_additional_cpu_queue import completion


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--contract',required=True)
    parser.add_argument('--contract-sha256',required=True)
    parser.add_argument('--step',required=True)
    args = parser.parse_args()
    contract_path = Path(args.contract)
    audit_path = BASE / ('dispatch_' + args.step + '_audit.json')
    report = dict(status='preflight',step=args.step,no_gpu_used=True)
    try:
        require(sha(contract_path) == args.contract_sha256, 'Dispatched execution contract changed')
        contract = read_json(contract_path)
        require(contract['status'] == 'frozen_cpu_execution_contract', 'Execution contract is not frozen')
        require(contract['cohort_freeze_sha256'] == COHORT_HASH, 'Contract cohort hash changed')
        for name, expected in contract['files'].items():
            require(sha(BASE / name) == expected, 'Execution code/weight/manifest changed: ' + name)
        report['frozen_preflight'] = preflight()
        step = contract['steps'][args.step]
        if step.get('requires_all_completion'):
            # Scoring cannot use failed queue state even if an old run_audit
            # exists. Also compare fresh audits to the queue's float32 sweep.
            for variant in ['full','seq']:
                state = read_json(BASE / ('planet_' + variant + '_cpu_queue_state.json'))
                completion(state, BASE / ('planet_' + variant), 'planet-' + variant)
            state = read_json(BASE / 'rmsf_net_resume_state.json')
            require(state['status'] == 'complete', 'Guarded RMSF-net cache verification is not complete')
            audit = read_json(BASE / 'rmsf_net/run_audit.json')
            require(audit['status'] == 'complete' and audit['model'] == 'rmsf-net', 'RMSF-net completion audit missing')
            require(audit['verified_items'] == audit['verified_records'] == 3063 and audit['verified_residues'] == 961321, 'RMSF-net cohort incomplete')
            require(audit['failure_count'] == 0 and audit['checker_sha256'] == HASHES['verify_additional_baselines.py'], 'RMSF-net audit has failures/altered checker')
            require(audit['provenance']['cohort_freeze_sha256'] == COHORT_HASH, 'RMSF-net cohort changed')
        report.update(status='running',command=step['command'],contract_sha256=args.contract_sha256)
        write_json(audit_path,report)
        env = dict(os.environ,CUDA_VISIBLE_DEVICES='',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
        result = subprocess.run(step['command'],env=env,stdin=subprocess.DEVNULL)
        require(result.returncode == 0, 'Guarded child exited with code ' + str(result.returncode))
        for path in step['check']:
            require(Path(path).is_file(), 'Required completion artifact missing: ' + path)
        report.update(status='complete',checks={path:sha(path) for path in step['check']})
        write_json(audit_path,report)
    except Exception as error:
        report.update(status='failed',error=type(error).__name__ + ': ' + str(error))
        write_json(audit_path,report)
        print(report['error'],flush=True)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
