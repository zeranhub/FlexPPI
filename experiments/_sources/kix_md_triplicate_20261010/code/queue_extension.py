"""Schedule only prespecified rep02/03; never rebuild or touch historical MD.

Remote Linux / OpenMM environment required. Per-job flock allows disjoint
state shards on separate GPUs while preventing a duplicate job writer.
Run only after restored prepared systems pass verification and new GPU smoke.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys

from run_md import PARAMS, STATES, now, run_job, sha, write_json


def verify(root):
    policy = json.loads((root / 'SCHEDULE_POLICY.json').read_text())
    if policy['allowed_repeats'] != [2, 3]:
        raise RuntimeError('This extension is restricted to repeats 2 and 3')
    if json.loads((root / 'protocol.json').read_text()) != PARAMS:
        raise RuntimeError('Physical protocol differs from the historical runner')
    ready_path = root / 'prepared' / 'READY.json'
    ready = json.loads(ready_path.read_text())
    if ready['parameters'] != PARAMS:
        raise RuntimeError('Restored prepared protocol differs from historical protocol')
    for state in STATES:
        folder = root / 'prepared' / state
        metadata = json.loads((folder / 'metadata.json').read_text())
        for name, key in [('system.xml', 'system_sha256'), ('protein.pdb', 'protein_pdb_sha256')]:
            if sha(folder / name) != metadata[key] or metadata[key] != ready['states'][state][key]:
                raise RuntimeError(f'Restored input hash mismatch: {state}/{name}')
        if not (folder / 'solvated.pdb').is_file():
            raise FileNotFoundError(folder / 'solvated.pdb')
        if sha(root / 'inputs' / f'{state}.pdb') != ready['inputs_sha256'][f'{state}.pdb']:
            raise RuntimeError(f'Original construct hash mismatch: {state}')
    if not (root / 'smoke' / 'PASSED.json').is_file():
        raise RuntimeError('Four new-platform CUDA smoke tests must pass first')
    return policy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--states', nargs='+', choices=STATES, default=list(STATES))
    parser.add_argument('--repeats', nargs='+', choices=[2, 3], type=int, default=[2, 3])
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    verify(root)
    if len(set(args.states)) != len(args.states) or len(set(args.repeats)) != len(args.repeats):
        raise ValueError('Duplicate job requested')
    if args.verify_only:
        print('RESTORED INPUTS AND NEW GPU SMOKE VERIFIED', flush=True)
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
