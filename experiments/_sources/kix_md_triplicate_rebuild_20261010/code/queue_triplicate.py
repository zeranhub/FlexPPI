"""Schedule three fresh repeats in a rebuilt batch; historical pilot excluded.

Remote Linux / OpenMM environment required. Per-job flock allows disjoint
state shards on separate GPUs while preventing a duplicate job writer.
Run only after newly built systems pass verification and CUDA smoke.
"""
import argparse
import csv
import fcntl
import json
import math
import os
from pathlib import Path
import subprocess
import sys

from run_md import PARAMS, STATES, now, run_job, sha, write_json


def verify(root):
    policy = json.loads((root / 'SCHEDULE_POLICY.json').read_text())
    if policy['allowed_repeats'] != [1, 2, 3]:
        raise RuntimeError('This new batch is restricted to repeats 1, 2, and 3')
    if json.loads((root / 'protocol.json').read_text()) != PARAMS:
        raise RuntimeError('Physical protocol differs from the historical runner')
    ready_path = root / 'prepared' / 'READY.json'
    ready = json.loads(ready_path.read_text())
    if ready['parameters'] != PARAMS:
        raise RuntimeError('New prepared protocol differs from historical protocol')
    input_audit = json.loads((root / 'prepared' / 'input_audit.json').read_text())
    if input_audit['status'] != 'all_four_rebuilt_inputs_passed' or not all(v['passed'] for v in input_audit['states'].values()):
        raise RuntimeError('New heavy-coordinate/force-field input audit must pass before production')
    for state in STATES:
        folder = root / 'prepared' / state
        metadata = json.loads((folder / 'metadata.json').read_text())
        for name, key in [('system.xml', 'system_sha256'), ('protein.pdb', 'protein_pdb_sha256')]:
            if sha(folder / name) != metadata[key] or metadata[key] != ready['states'][state][key]:
                raise RuntimeError(f'Prepared input hash mismatch: {state}/{name}')
        if sha(folder / 'solvated.pdb') != input_audit['states'][state]['solvated_pdb_sha256']:
            raise RuntimeError('Audited full solvated coordinates changed')
        if not (folder / 'solvated.pdb').is_file():
            raise FileNotFoundError(folder / 'solvated.pdb')
        if sha(root / 'inputs' / f'{state}.pdb') != ready['inputs_sha256'][f'{state}.pdb']:
            raise RuntimeError(f'Original construct hash mismatch: {state}')
    frozen = json.loads((root / 'frozen_predictions' / 'FROZEN.json').read_text())
    if sha(root / 'frozen_predictions' / 'residue_predictions.csv') != frozen['prediction_csv_sha256']:
        raise RuntimeError('Frozen model/Static predictions are missing or changed')
    for state in STATES:
        if sha(root / 'inputs' / f'{state}.pdb') != frozen['input_sha256'][state]:
            raise RuntimeError('Prediction construct differs from MD construct')
    with (root / 'frozen_predictions' / 'residue_predictions.csv').open(encoding='utf-8-sig', newline='') as stream:
        prediction_rows = [row for row in csv.DictReader(stream) if row['chain'] == 'B']
    methods = set(row['method'] for row in prediction_rows)
    if not {'full', 'static'}.issubset(methods):
        raise RuntimeError('Both Full and Static must be frozen before production')
    for state in STATES:
        for method in methods:
            rows = [row for row in prediction_rows if row['state'] == state and row['method'] == method]
            if sorted(int(row['author_residue_number']) for row in rows) != list(range(586, 673)):
                raise RuntimeError('Frozen method must cover all KIX sites exactly once')
            for row in rows:
                if not math.isfinite(float(row['predicted_rmsf_angstrom'])):
                    raise RuntimeError('Non-finite frozen RMSF')
                if row['amino_acid'] not in ('A', 'G') and not math.isfinite(float(row['predicted_hchi_nats'])):
                    raise RuntimeError('Non-finite frozen Hchi at a valid side-chain site')
    if not (root / 'smoke' / 'PASSED.json').is_file():
        raise RuntimeError('Four new-platform CUDA smoke tests must pass first')
    return policy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--states', nargs='+', choices=STATES, default=list(STATES))
    parser.add_argument('--repeats', nargs='+', choices=[1, 2, 3], type=int, default=[1, 2, 3])
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    verify(root)
    if len(set(args.states)) != len(args.states) or len(set(args.repeats)) != len(args.repeats):
        raise ValueError('Duplicate job requested')
    if args.verify_only:
        print('NEW BATCH INPUTS AND CUDA SMOKE VERIFIED', flush=True)
        return
    worker_key = '_'.join(args.states)
    worker_status = root / f'extension_worker_{worker_key}.json'
    completed = []
    for repeat in args.repeats:
        for state in args.states:
            if (root / 'STOP_AFTER_SEGMENT').exists():
                raise SystemExit(75)
            folder = root / 'runs' / state / f'rep{repeat:02d}'
            folder.mkdir(parents=True, exist_ok=True)
            with (folder / 'writer.lock').open('a+') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                write_json(worker_status, dict(pid=os.getpid(), state=state, repeat=repeat,
                           status='running', cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES'),
                           updated_utc=now(), completed=completed))
                run_job(root, state, repeat)
                done = json.loads((folder / 'DONE.json').read_text())
                if done['state'] != state or done['repeat'] != repeat or done['production_ns'] != 100 or done['equilibration_ns'] != 10 or done['final_step'] != 55000000:
                    raise RuntimeError(f'Invalid completion marker: {state}/rep{repeat:02d}')
                summary_path = folder / 'analysis' / 'summary.json'
                if not summary_path.exists():
                    subprocess.run([sys.executable, str(root / 'code' / 'analyze_job.py'),
                                    '--protein-pdb', str(root / 'prepared' / state / 'protein.pdb'),
                                    '--job-dir', str(folder)], check=True)
                summary = json.loads(summary_path.read_text())
                if summary['state'] != state or summary['replicate'] != f'rep{repeat:02d}' or summary['n_frames'] != 100000 or len(summary['segments']) != 100:
                    raise RuntimeError('Incomplete or wrong-repeat analysis')
                if sha(folder / 'analysis' / 'residue_metrics.csv') != summary['output_csv_sha256']:
                    raise RuntimeError('Analysis output hash mismatch')
                completed.append(f'{state}/rep{repeat:02d}')
                write_json(worker_status, dict(pid=os.getpid(), status='job_completed',
                           updated_utc=now(), completed=completed))
    write_json(worker_status, dict(status='worker_completed', finished_utc=now(), completed=completed,
               caution='MD production and per-repeat metrics complete; replicate and delta QA still required'))


if __name__ == '__main__':
    main()
