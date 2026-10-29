"""Read frozen cohort without changing a frozen byte or residue order."""
from pathlib import Path
import csv
import hashlib
import json
import os
import numpy as np

ROOT = Path(__file__).resolve().parent
COHORT = Path('E:/Project/FlexPPI/experiments/_sources/baselines/cohort')
OLD = COHORT.parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    tmp = path.with_suffix('.tmp.json')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
    os.replace(tmp, path)


def load_records(smoke=False):
    freeze = json.loads((COHORT / 'cohort_freeze.json').read_text())
    assert freeze['records'] == 3063 and freeze['valid_rmsf'] == 961321
    for filename in ['structure_manifest.tsv', 'record_index.json']:
        assert sha(COHORT / filename) == freeze['files'][filename], filename
    rows = list(csv.DictReader((COHORT / 'structure_manifest.tsv').open(encoding='utf-8', newline=''), delimiter='\t'))
    assert len(rows) == 3063
    assert sum(int(r['n_residues']) for r in rows) == 961321
    if smoke:
        selected = {}
        for row in rows:
            selected.setdefault(row['cohort'], row)
        rows = list(selected.values())
        assert len(rows) == 3, list(selected)
    return rows, freeze


def local_input(record):
    # Filename relocation only; the canonical bytes must remain exact.
    filename = record['container_pdb_path'].replace('\\', '/').rsplit('/', 1)[-1]
    path = COHORT / 'structures' / filename
    assert sha(path) == record['input_sha256'], record['sample_id']
    return path


def record_key(sample):
    return 'record_' + hashlib.sha256(sample.encode()).hexdigest()[:24]


def parse_ca(path):
    coords, chain, author, name = [], [], [], []
    seen = set()
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        if line.startswith('ENDMDL') and coords:
            break
        if not line.startswith('ATOM  ') or line[12:16].strip() != 'CA':
            continue
        if line[16:17] not in (' ', 'A'):
            continue
        ch = line[21:22].strip() or '_'
        rid = str(int(line[22:26])) + line[26:27].strip()
        key = (ch, rid)
        assert key not in seen, (path, key)
        seen.add(key)
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        chain.append(ch)
        author.append(rid)
        name.append(line[17:20].strip())
    xyz = np.asarray(coords, dtype=np.float64)
    assert xyz.ndim == 2 and xyz.shape[1] == 3 and np.isfinite(xyz).all()
    return xyz, np.asarray(chain), np.asarray(author), np.asarray(name)


def valid_cached(path, status_path, record, value_fields):
    if not path.exists() or not status_path.exists():
        return False
    status = json.loads(status_path.read_text())
    if status.get('status') != 'complete':
        return False
    assert status['input_sha256'] == record['input_sha256']
    assert status['output_sha256'] == sha(path)
    with np.load(path, allow_pickle=False) as saved:
        for field in value_fields:
            assert saved[field].shape == (int(record['n_residues']),)
            assert np.isfinite(saved[field]).all()
    return True
