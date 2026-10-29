"""Persistent CPU-only PLANET queue; frozen adapter/model files stay unchanged.

One invocation serves one variant. The variant lock is inherited by the child
so a surviving inference process still prevents a duplicate queue. Completion
requires the separate frozen verifier plus exact float32 output checks.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import runpy
import shlex
import subprocess
import sys
import time
import traceback

ROOT = Path('/root/flexppi_completion_20261010')
BASE = ROOT / 'baselines'
PYTHON = ROOT / 'planet_env' / 'bin' / 'python'
GUARDS = {
    'cohort/cohort_freeze.json': 'b58795275fa2843a502fdbb68f417f9a2b0b80f76e28d226ad5d58c9bf430a91',
    'planet_payload_freeze.json': '8f06e3b0adf617b6d3f846421996e6ff44ce070da1b6f831ddb59c51b5bf1d39',
    'planet_cpu_smoke_audit_20261010/smoke_audit.json': 'b92b7b88afd64f74c8c22d9d77a6594d8c5cc3be6515ee7de27d84d3d07c4e56',
    'planet_cpu_smoke_audit_20261010/SMOKE_PASSED': 'c42616cf3f7224622b6eb5428bbb7bc44b8f9169a6e757748dc91f90dfe06345',
    'run_planet_md_3063.py': 'fae56fc2df6dfa50f8aa1301020d3a07b7bb70ccffc5408d615e3e99c93dc910',
    'verify_additional_baselines.py': '7353a155f4be8a79656cf81309cb2a60fe05329c1bd77abeda17d82fcdd49d18',
}
WEIGHTS = {
    'full': ('rocketshp_pretrained_20250505_v1.ckpt', 'b31a5eb353e773c346a413fbca7e63735a6ab53249a628073a7c01ceee393bf9'),
    'seq': ('rocketshp_pretrained_20250507_v1_seq.ckpt', '54f8ee6df8ed0c45e33995431b720cb7676b38d7521a5101602fae58f574e243'),
}
CPU_ENV = {
    'CUDA_VISIBLE_DEVICES': '', 'OMP_NUM_THREADS': '2', 'MKL_NUM_THREADS': '2',
    'OPENBLAS_NUM_THREADS': '2', 'NUMEXPR_NUM_THREADS': '2', 'VECLIB_MAXIMUM_THREADS': '2',
    'PYTHONDONTWRITEBYTECODE': '1', 'HF_HUB_OFFLINE': '1', 'TRANSFORMERS_OFFLINE': '1',
    'HF_HOME': str(BASE / 'hf_cache'), 'HF_HUB_CACHE': str(BASE / 'hf_cache' / 'hub'),
}


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path, value):
    temporary = Path(str(path) + '.' + str(os.getpid()) + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    os.replace(temporary, path)


def commands(variant):
    output = BASE / ('planet_' + variant)
    adapter_args = [
        '--cohort', str(BASE / 'cohort'), '--source-root', str(BASE / 'PLANET-MD'),
        '--weights', str(BASE / 'weights' / WEIGHTS[variant][0]), '--variant', variant,
        '--output', str(output), '--device', 'cpu', '--retry-failed',
    ]
    worker = [str(PYTHON), '-u', str(Path(__file__).resolve()), '--variant', variant, '--worker']
    verifier = [
        str(PYTHON), '-u', str(BASE / 'verify_additional_baselines.py'),
        '--cohort', str(BASE / 'cohort'), '--output', str(output), '--model', 'planet-' + variant,
        '--expected-cohort-freeze-sha256', GUARDS['cohort/cohort_freeze.json'],
        '--inference-script', str(BASE / 'run_planet_md_3063.py'),
        '--weights', str(BASE / 'weights' / WEIGHTS[variant][0]),
    ]
    return adapter_args, worker, verifier


def preflight(variant):
    require(Path(sys.executable).resolve() == PYTHON.resolve(), 'Wrong isolated PLANET environment')
    marker = ROOT / 'PLANET_ENV_READY'
    require(marker.is_file() and marker.stat().st_size > 0, 'PLANET_ENV_READY missing or empty')
    hashes = {}
    for name, expected in GUARDS.items():
        actual = sha(BASE / name)
        require(actual == expected, 'Frozen contract hash mismatch: ' + name)
        hashes[name] = actual
    smoke = read_json(BASE / 'planet_cpu_smoke_audit_20261010' / 'smoke_audit.json')
    require(smoke['status'] == 'full_and_seq_cpu_smoke_passed_2_chains_each', 'CPU smoke did not pass')
    passed = read_json(BASE / 'planet_cpu_smoke_audit_20261010' / 'SMOKE_PASSED')
    require(passed['status'] == smoke['status'] and passed['audit_sha256'] == hashes['planet_cpu_smoke_audit_20261010/smoke_audit.json'], 'CPU smoke PASSED marker mismatch')
    require(all(item['unchanged'] for item in smoke['frozen_source_cache_and_weights_unchanged_after_smoke']), 'Smoke observed source/cache/weight change')
    require(sha(BASE / 'planet_cpu_smoke_audit_20261010' / 'payload_audit.json') == smoke['payload_audit_sha256'], 'Smoke payload audit changed')
    require(smoke['environment'] == str(PYTHON), 'CPU smoke environment mismatch')
    require(smoke['released_float32_weights_preserved'] and not smoke['dtype_conversion'], 'Smoke dtype policy mismatch')
    require(not smoke['gpu_used'] and not smoke['trained'] and not smoke['calibrated'], 'Smoke protocol mismatch')
    smoke_variant = smoke['outputs'][variant]
    require(sha(Path(smoke_variant['folder']) / 'metadata.json') == smoke_variant['metadata_sha256'], 'Smoke metadata changed')
    require(sha(smoke_variant['log']) == smoke_variant['log_sha256'], 'Smoke log changed')
    require(len(smoke_variant['samples']) == 2, 'Smoke chain count mismatch')
    for item in smoke_variant['samples']:
        require(item['status'] == 'complete' and item['dtype'] == 'float32' and item['nm_to_angstrom_exact'], 'Smoke output was not valid float32')
        require(sha(Path(smoke_variant['folder']) / (item['input_id'] + '.npz')) == item['output_sha256'], 'Smoke prediction changed')
    frozen = read_json(BASE / 'planet_payload_freeze.json')
    require(frozen['status'] == 'frozen_payload', 'Payload is not frozen')
    verified = []
    for entry in frozen['files']:
        path = BASE / entry['path']
        actual = sha(path)
        require(path.stat().st_size == entry['bytes'] and actual == entry['sha256'], 'Payload changed: ' + entry['path'])
        verified.append(dict(path=entry['path'], sha256=actual, bytes=path.stat().st_size))
    actual_tar = sha(BASE / 'planet_payload.tar')
    require(actual_tar == frozen['payload_sha256'], 'Frozen payload archive changed')
    require(sha(BASE / 'weights' / WEIGHTS[variant][0]) == WEIGHTS[variant][1], 'Released checkpoint changed')
    freeze = read_json(BASE / 'cohort' / 'cohort_freeze.json')
    require(freeze['status'] == 'frozen' and freeze['records'] == 3063 and freeze['chains'] == 4973, 'Cohort size/status mismatch')
    require(freeze['total_residues'] == freeze['valid_rmsf'] == 961321, 'RMSF residue count mismatch')
    require(freeze['valid_chi_entropy'] == 824988, 'Entropy residue count mismatch')
    for name, expected in freeze['files'].items():
        require(sha(BASE / 'cohort' / name) == expected, 'Frozen cohort file changed: ' + name)
    with (BASE / 'cohort' / 'chain_manifest.csv').open(encoding='utf-8', newline='') as handle:
        chains = list(csv.DictReader(handle))
    require(len(chains) == len({c['input_id'] for c in chains}) == 4973, 'Duplicate/missing frozen chain identity')
    require(len({c['sample_id'] for c in chains}) == 3063, 'Missing source record identity')
    for chain in chains:
        require(len(chain['sequence']) == int(chain['n_residues']), 'Sequence length mismatch: ' + chain['input_id'])
        require(sha(Path(chain['pdb_file'])) == chain['input_sha256'], 'Frozen structure changed: ' + chain['input_id'])
    require(max(int(c['n_residues']) for c in chains) == 1024, 'Frozen maximum chain length changed')
    adapter = (BASE / 'run_planet_md_3063.py').read_text(encoding='utf-8')
    require("get_model(device='cpu').to(device).eval()" in adapter, 'ESM CPU loading policy changed')
    require(not any(term in adapter for term in ('.half(', '.bfloat16(', 'autocast(', 'set_default_dtype(')), 'Adapter converts dtype')
    report = dict(status='passed', variant=variant, checked_utc=utc(), cpu_environment=CPU_ENV,
        guard_hashes=hashes, environment_ready_sha256=sha(marker), environment_ready_text=marker.read_text().strip(),
        payload_archive_sha256=actual_tar, verified_payload_files=verified,
        cohort_files=freeze['files'], chains_verified=4973, records_verified=3063, residues_verified=961321,
        max_chain_residues=1024, weights_sha256=WEIGHTS[variant][1], released_float32_weights=True,
        no_calibration=True, no_training=True, no_gpu=True, wrapper_sha256=sha(__file__))
    write_json(BASE / ('planet_' + variant + '_cpu_queue_preflight.json'), report)
    return report


def worker(variant):
    os.environ.update(CPU_ENV)
    import torch
    torch.set_num_threads(2)
    torch.set_num_interop_threads(2)
    require(torch.get_num_threads() == torch.get_num_interop_threads() == 2, 'Torch thread limit not applied')
    require(torch.get_default_dtype() == torch.float32, 'Torch default dtype is not float32')
    require(os.environ['CUDA_VISIBLE_DEVICES'] == '', 'CUDA was not hidden')
    adapter_args, _, _ = commands(variant)
    worker_path = BASE / ('planet_' + variant + '_cpu_worker_audit.json')
    report = dict(
        status='running', started_utc=utc(), pid=os.getpid(), parent_pid=os.getppid(), variant=variant,
        torch=str(torch.__version__), intra_op_threads=2, inter_op_threads=2, default_dtype='torch.float32',
        environment=CPU_ENV, adapter_sha256=sha(BASE / 'run_planet_md_3063.py'),
        adapter_arguments=adapter_args, trained=False, calibrated=False, device='cpu')
    write_json(worker_path, report)
    print('CPU WORKER', variant, 'pid', os.getpid(), 'Torch threads=2/2; float32; device=cpu', flush=True)
    sys.argv = [str(BASE / 'run_planet_md_3063.py')] + adapter_args
    try:
        runpy.run_path(str(BASE / 'run_planet_md_3063.py'), run_name='__main__')
        report.update(status='inference_exited', ended_utc=utc(), independently_verified=False)
        write_json(worker_path, report)
    except Exception as error:
        report.update(status='failed', ended_utc=utc(), error=type(error).__name__ + ': ' + str(error))
        write_json(worker_path, report)
        raise


def exact_outputs(variant):
    import numpy as np
    output = BASE / ('planet_' + variant)
    audit_path = output / 'run_audit.json'
    require(audit_path.is_file(), 'Independent run_audit.json was not produced')
    audit = read_json(audit_path)
    require(audit['status'] == 'complete' and audit['model'] == 'planet-' + variant, 'Independent verifier did not complete')
    require(audit['verified_items'] == audit['expected_items'] == 4973, 'Not all 4973 chains passed')
    require(audit['verified_records'] == 3063 and audit['verified_residues'] == 961321, 'Complete cohort totals mismatch')
    require(audit['failure_count'] == 0 and audit['checker_sha256'] == GUARDS['verify_additional_baselines.py'], 'Verifier failures or altered checker')
    metadata = read_json(output / 'metadata.json')
    require(metadata['device'] == 'cpu', 'Inference device was not fixed CPU')
    count = 0
    for row in audit['verified_outputs']:
        with np.load(output / (row['input_id'] + '.npz'), allow_pickle=False) as saved:
            native, converted = saved['rmsf_native_nm'], saved['rmsf_angstrom']
            require(native.dtype == converted.dtype == np.dtype('float32'), 'Prediction is not exact float32: ' + row['input_id'])
            require(np.array_equal(converted, 10.0 * native), 'Unit conversion differs: ' + row['input_id'])
        count += 1
    require(count == 4973, 'Float32 sweep incomplete')
    return dict(status='passed', checked_utc=utc(), chains=count, exact_float32=True,
                exact_nm_to_angstrom=True, run_audit_sha256=sha(audit_path))


def wait_child(command, phase, state, state_path, lock):
    print('START', phase, shlex.join(command), flush=True)
    child = subprocess.Popen(command, env=dict(os.environ, **CPU_ENV), stdin=subprocess.DEVNULL,
                             close_fds=True, pass_fds=(lock.fileno(),))
    event = dict(phase=phase, pid=child.pid, parent_pid=os.getpid(), command=command,
                 command_shell=shlex.join(command), started_utc=utc())
    state['children'].append(event)
    state.update(phase=phase, child_pid=child.pid, heartbeat_utc=utc())
    write_json(state_path, state)
    while child.poll() is None:
        time.sleep(30)
        state['heartbeat_utc'] = utc()
        if phase == 'inference':
            completed, failed = 0, 0
            for path in Path(state['output']).glob('chain_*.json'):
                try:
                    item = read_json(path)
                    completed += item.get('status') == 'complete'
                    failed += item.get('status') == 'failed'
                except (OSError, ValueError):
                    continue
            state['observed_items'] = dict(complete=completed, failed=failed, expected=4973,
                                           independently_verified=False)
        write_json(state_path, state)
    event.update(ended_utc=utc(), exit_code=child.returncode)
    state.update(child_pid=None, heartbeat_utc=utc())
    write_json(state_path, state)
    require(child.returncode == 0, phase + ' exited with code ' + str(child.returncode))


def queue(variant):
    os.environ.update(CPU_ENV)
    state_path = BASE / ('planet_' + variant + '_cpu_queue_state.json')
    lock_path = BASE / ('planet_' + variant + '_cpu_queue.lock')
    with lock_path.open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print('Another queue/child already holds', lock_path, flush=True)
            return 2
        adapter_args, inference, verifier = commands(variant)
        state = dict(status='running', phase='preflight', variant=variant, device='cpu',
            started_utc=utc(), ended_utc=None, parent_pid=os.getpid(), launcher_parent_pid=os.getppid(),
            child_pid=None, children=[], command=sys.argv, command_shell=shlex.join([str(PYTHON)] + sys.argv),
            output=str(BASE / ('planet_' + variant)), log=str(BASE / ('planet_' + variant + '_cpu_queue.log')),
            lock=str(lock_path), environment=CPU_ENV, expected_chains=4973, expected_records=3063,
            wrapper_sha256=sha(__file__), inference_command=inference, verifier_command=verifier,
            frozen_contract=GUARDS, weights_sha256=WEIGHTS[variant][1], trained=False, calibrated=False,
            completion_rule='Frozen standalone verifier exit0 AND run_audit complete 4973 chains/3063 records/961321 residues AND exact float32 sweep AND final frozen-input recheck')
        write_json(state_path, state)
        try:
            state['preflight'] = preflight(variant)
            write_json(state_path, state)
            wait_child(inference, 'inference', state, state_path, lock)
            wait_child(verifier, 'verification', state, state_path, lock)
            state['output_dtype_audit'] = exact_outputs(variant)
            state['final_frozen_input_recheck'] = preflight(variant)
            state.update(status='complete', phase='finished', ended_utc=utc(), heartbeat_utc=utc())
            write_json(state_path, state)
            print('COMPLETE', variant, 'independently verified 4973 chains', flush=True)
            return 0
        except Exception as error:
            state.update(status='failed', phase='failed', ended_utc=utc(), heartbeat_utc=utc(),
                         error=type(error).__name__ + ': ' + str(error), traceback=traceback.format_exc())
            write_json(state_path, state)
            traceback.print_exc()
            return 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--variant', choices=['full', 'seq'], required=True)
    parser.add_argument('--worker', action='store_true', help='Internal CPU inference child only')
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args()
    if args.worker:
        worker(args.variant)
    elif args.preflight_only:
        print(json.dumps(preflight(args.variant), indent=2))
    else:
        raise SystemExit(queue(args.variant))


if __name__ == '__main__':
    main()
