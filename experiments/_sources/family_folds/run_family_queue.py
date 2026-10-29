"""Execute the frozen ten-family-fold fits only when explicitly dispatched."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-root', required=True)
    parser.add_argument('--variant', choices=['full', 'intrachain', 'all'], default='all')
    parser.add_argument('--fold', type=int, choices=range(5))
    parser.add_argument('--evaluate', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    freeze = json.loads((root / 'training_freeze.json').read_text())
    assert freeze['seed'] == 7 and freeze['epochs'] == 40 and freeze['fits'] == 10
    for name, expected in freeze['files'].items():
        assert sha(root / name) == expected, 'Frozen asset changed: ' + name
    os.chdir(root)
    env = dict(os.environ, PYTHONPATH=args.source_root, OMP_NUM_THREADS='2', MKL_NUM_THREADS='2')
    variants = ['full', 'intrachain'] if args.variant == 'all' else [args.variant]
    folds = list(range(5)) if args.fold is None else [args.fold]
    state_path = root / f'queue_{args.variant}_{args.fold if args.fold is not None else "all"}.json'
    state = dict(status='running', freeze_sha256=sha(root / 'training_freeze.json'), jobs=[], started=time.time())
    for variant in variants:
        for fold in folds:
            checkpoint = root / f'runs/{variant}/fold_{fold}/epoch_040.pt'
            item = dict(variant=variant, fold=fold, status='running', started=time.time())
            state['jobs'].append(item)
            state_path.write_text(json.dumps(state, indent=2))
            if checkpoint.exists():
                import torch
                saved = torch.load(checkpoint, map_location='cpu', weights_only=False)
                assert saved['epoch'] == 40 and saved['config']['seed'] == 7
                assert saved['config']['model']['information_ablation'] == variant
                item.update(status='already_complete', checkpoint_sha256=sha(checkpoint), finished=time.time())
            else:
                log = root / f'logs/{variant}_fold{fold}.log'
                log.parent.mkdir(exist_ok=True)
                with log.open('a') as handle:
                    result = subprocess.run([sys.executable, '-u', '-m', 'flexppi', 'train', '--config', f'configs/{variant}_fold{fold}.yaml'], env=env, stdout=handle, stderr=subprocess.STDOUT)
                item.update(exit_code=result.returncode, finished=time.time(), status='complete' if result.returncode == 0 and checkpoint.exists() else 'failed')
                if item['status'] == 'complete':
                    item['checkpoint_sha256'] = sha(checkpoint)
                else:
                    state['status'] = 'failed'
                    state_path.write_text(json.dumps(state, indent=2))
                    raise SystemExit('Training failed; inspect ' + str(log))
            state_path.write_text(json.dumps(state, indent=2))
    state['status'] = 'training_complete'
    if args.evaluate:
        assert all((root / f'runs/{v}/fold_{f}/epoch_040.pt').exists() for v in ['full', 'intrachain'] for f in range(5))
        result = subprocess.run([sys.executable, '-u', 'evaluate_family_models.py'], env=env)
        state.update(status='complete' if result.returncode == 0 else 'evaluation_failed', evaluation_exit_code=result.returncode)
    state_path.write_text(json.dumps(state, indent=2))


if __name__ == '__main__':
    main()
