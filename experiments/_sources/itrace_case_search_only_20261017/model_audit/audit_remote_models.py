"""Read-only inspection of saved checkpoints and their actual source manifests."""
import collections
import hashlib
import json
from pathlib import Path
import torch


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def manifest(path):
    p = Path(path)
    rows = [json.loads(x) for x in p.read_text().splitlines() if x.strip()]
    return dict(path=str(p), sha256=sha(p), records=len(rows),
                sources=dict(collections.Counter(r.get('dataset') for r in rows)),
                itrace_sample_ids=[r['sample_id'] for r in rows if str(r.get('dataset', '')).lower() == 'itrace'],
                pdb_3kxf_records=[r for r in rows if '3kxf' in json.dumps(r).lower()])


roots = [Path('/root/flexppi_completion_20261010/family_folds'),
         Path('/root/flexppi_revision_20261002/corrected_training_v1')]
out = {'checkpoints': [], 'source_trainer_sha256': {}, 'note': 'No models trained; read-only checkpoint and manifest inspection.'}
for root in roots:
    paths = sorted((root / 'runs').glob('*/fold_*/epoch_040.pt'))
    for p in paths:
        saved = torch.load(p, map_location='cpu', weights_only=False)
        c = saved['config']
        item = dict(path=str(p), sha256=sha(p), epoch=saved.get('epoch'), keys=list(saved), config=c, manifests={})
        for key in ['state_manifest', 'validation_state_manifest', 'test_state_manifest']:
            value = c.get('data', {}).get(key)
            if value:
                mp = Path(value)
                if not mp.is_absolute():
                    mp = root / mp
                item['manifests'][key] = manifest(mp)
        out['checkpoints'].append(item)
trainer = roots[1] / 'src/flexppi/trainer.py'
out['source_trainer_sha256'][str(trainer)] = sha(trainer)
print(json.dumps(out, ensure_ascii=False, indent=2))
