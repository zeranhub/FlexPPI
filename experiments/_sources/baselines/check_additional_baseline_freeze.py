"""CPU-only, independent execution-contract guard; never starts inference."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path('/root/flexppi_completion_20261010')
BASE = ROOT / 'baselines'
COHORT_HASH = 'b58795275fa2843a502fdbb68f417f9a2b0b80f76e28d226ad5d58c9bf430a91'
SMOKE_HASH = 'b92b7b88afd64f74c8c22d9d77a6594d8c5cc3be6515ee7de27d84d3d07c4e56'
HASHES = {
    'cohort/cohort_freeze.json': COHORT_HASH,
    'planet_payload_freeze.json': '8f06e3b0adf617b6d3f846421996e6ff44ce070da1b6f831ddb59c51b5bf1d39',
    'run_rmsf_net_3063.py': 'e7fc76342048cb2827142f8e1a22e0a069bbcfc7b1db5c1d1c7a71113bc5ebf2',
    'rmsf_net_pdb_adapter.py': 'd4e3606d5d7678251e7615e87987dadfb65a2162b9dd0ab5427368fe13e5c6e1',
    'RMSF-net/model.py': '5a59cfde470c3b1ec18397517dfbf1f6d9bb1e1284beb69c2dffc308e1ea22e8',
    'RMSF-net/model_cpu/sim_model.pth': 'cd12dce9489ed57c5b0b8d5400161e1b4e165963622dc78323fdeb07651d8c9b',
    'run_planet_md_3063.py': 'fae56fc2df6dfa50f8aa1301020d3a07b7bb70ccffc5408d615e3e99c93dc910',
    'verify_additional_baselines.py': '7353a155f4be8a79656cf81309cb2a60fe05329c1bd77abeda17d82fcdd49d18',
    'score_common_residues.py': '2823f0a6cb72da47b48144a175d5eb4011811160d7a1c6a744a074c0d71bb619',
    'planet_cpu_smoke_audit_20261010/smoke_audit.json': SMOKE_HASH,
    'planet_cpu_smoke_audit_20261010/SMOKE_PASSED': 'c42616cf3f7224622b6eb5428bbb7bc44b8f9169a6e757748dc91f90dfe06345',
    'run_planet_cpu_queue.py': '034a7fdfc5b6bc39164bb43f23518d09586863af186ece495275965c3910e602',
}


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


def smoke_check():
    folder = BASE / 'planet_cpu_smoke_audit_20261010'
    marker = read_json(folder / 'SMOKE_PASSED')
    require(marker['status'] == 'full_and_seq_cpu_smoke_passed_2_chains_each', 'CPU SMOKE_PASSED status mismatch')
    require(marker['audit_sha256'] == SMOKE_HASH, 'CPU SMOKE_PASSED audit hash mismatch')
    smoke = read_json(folder / 'smoke_audit.json')
    require(smoke['status'] == marker['status'], 'Smoke audit did not pass')
    require(smoke['environment'] == str(ROOT / 'planet_env/bin/python'), 'Smoke environment mismatch')
    require(smoke['released_float32_weights_preserved'] and not smoke['dtype_conversion'], 'Smoke dtype policy mismatch')
    require(not smoke['gpu_used'] and not smoke['trained'] and not smoke['calibrated'], 'Smoke protocol mismatch')
    require(all(row['unchanged'] for row in smoke['frozen_source_cache_and_weights_unchanged_after_smoke']), 'Smoke changed frozen payload')
    rows = []
    for variant in ['full', 'seq']:
        block = smoke['outputs'][variant]
        output = Path(block['folder'])
        require(sha(output / 'metadata.json') == block['metadata_sha256'], 'Smoke metadata changed: ' + variant)
        require(sha(block['log']) == block['log_sha256'], 'Smoke log changed: ' + variant)
        require(len(block['samples']) == 2, 'Smoke must contain two chains per model')
        metadata = read_json(output / 'metadata.json')
        require(metadata['device'] == 'cpu' and metadata['variant'] == variant, 'Smoke model/device mismatch')
        require(metadata['code_sha256'] == HASHES['run_planet_md_3063.py'], 'Smoke adapter mismatch')
        for row in block['samples']:
            require(row['status'] == 'complete' and row['dtype'] == 'float32' and row['nm_to_angstrom_exact'], 'Smoke prediction did not pass')
            path = output / (row['input_id'] + '.npz')
            require(sha(path) == row['output_sha256'], 'Smoke output changed: ' + row['input_id'])
            status = read_json(output / (row['input_id'] + '.json'))
            require(status['status'] == 'complete' and status['output_sha256'] == row['output_sha256'], 'Smoke item status mismatch')
            with np.load(path, allow_pickle=False) as saved:
                native, converted = saved['rmsf_native_nm'], saved['rmsf_angstrom']
                require(native.shape == converted.shape == (int(status['residues']),), 'Smoke residue shape mismatch')
                require(native.dtype == converted.dtype == np.dtype('float32'), 'Smoke predictions are not float32')
                require(np.isfinite(native).all() and np.isfinite(converted).all(), 'Smoke nonfinite prediction')
                require(np.array_equal(converted, 10.0 * native), 'Smoke unit conversion changed')
            rows.append({'variant': variant, 'input_id': row['input_id'], 'output_sha256': row['output_sha256']})
    return rows


def preflight():
    verified = {}
    for name, expected in HASHES.items():
        actual = sha(BASE / name)
        require(actual == expected, 'Frozen execution file mismatch: ' + name)
        verified[name] = actual
    marker = ROOT / 'PLANET_ENV_READY'
    require(marker.is_file() and marker.stat().st_size > 0, 'PLANET_ENV_READY missing or empty')
    freeze = read_json(BASE / 'cohort/cohort_freeze.json')
    require(freeze['status'] == 'frozen' and freeze['records'] == 3063 and freeze['chains'] == 4973, 'Frozen cohort membership mismatch')
    require(freeze['total_residues'] == freeze['valid_rmsf'] == 961321 and freeze['valid_chi_entropy'] == 824988, 'Frozen effective residues mismatch')
    for name, expected in freeze['files'].items():
        require(sha(BASE / 'cohort' / name) == expected, 'Frozen cohort data mismatch: ' + name)
    for name, delimiter, path_field, expected_n in [('structure_manifest.tsv','\t','container_pdb_path',3063), ('chain_manifest.csv',',','pdb_file',4973)]:
        with (BASE / 'cohort' / name).open(encoding='utf-8', newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter=delimiter))
        require(len(rows) == expected_n, 'Wrong manifest size: ' + name)
        for row in rows:
            require(sha(row[path_field]) == row['input_sha256'], 'Frozen PDB input mismatch: ' + row.get('input_id',row['sample_id']))
    payload = read_json(BASE / 'planet_payload_freeze.json')
    require(payload['status'] == 'frozen_payload', 'PLANET payload status mismatch')
    for row in payload['files']:
        path = BASE / row['path']
        require(path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], 'Frozen PLANET payload mismatch: ' + row['path'])
    smoke = smoke_check()
    return dict(status='passed', no_gpu_used=True, no_inference_performed=True,
        checked_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        cohort_freeze_sha256=COHORT_HASH, records=3063, chains=4973, valid_rmsf=961321,
        verified_execution_files=verified, frozen_payload_files=len(payload['files']),
        frozen_payload_bytes=payload['bytes'], cpu_smoke_passed=smoke,
        environment_ready_sha256=sha(marker), checker_sha256=sha(__file__))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--audit', default=str(BASE / 'execution_freeze_preflight.json'))
    args = parser.parse_args()
    path = Path(args.audit)
    try:
        report = preflight()
    except Exception as error:
        report = dict(status='failed', error=type(error).__name__ + ': ' + str(error), no_gpu_used=True, no_inference_performed=True)
        write_json(path, report)
        print(json.dumps(report, indent=2))
        raise SystemExit(1)
    write_json(path, report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['verified_execution_files','cpu_smoke_passed']}, indent=2))


if __name__ == '__main__':
    main()
