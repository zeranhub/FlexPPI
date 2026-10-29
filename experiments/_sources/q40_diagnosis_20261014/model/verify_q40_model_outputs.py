"""Independently verify the small Q40 diagnostic tables against frozen sources."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
result = json.loads((ROOT / 'diagnosis_summary.json').read_text(encoding='utf8'))
paths = list(result['sources_sha256'])
find = lambda suffix: Path(next(p for p in paths if p.replace('\\', '/').endswith(suffix)))
checks = []

def check(name, passed):
    checks.append({'name': name, 'passed': bool(passed)})
    assert passed, name

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

for path, expected in result['sources_sha256'].items():
    check('source unchanged: ' + Path(path).name, sha(path) == expected)
full = np.load(find('evaluation/full/A003R.npz'), allow_pickle=False)
intra = np.load(find('evaluation/intrachain/A003R.npz'), allow_pickle=False)
all66 = pd.read_csv(ROOT / 'all66_frozen_model_comparison.csv')
q = pd.read_csv(ROOT / 'all6_Gln_frozen_model_comparison.csv')
indices = np.flatnonzero((full['chain'] == 'A') & full['valid'][:, 1])
check('all66 positions exact and unique', len(all66) == 66 and all66.residue_id.nunique() == 66
      and np.array_equal(all66.residue_index, indices))
q_expected = all66[all66.amino_acid == 'Q'].reset_index(drop=True)
check('all six GLN in source order', np.array_equal(q.residue_id, [2, 31, 40, 41, 49, 62]))
for col in all66:
    if col == 'amino_acid':
        check('Q subset ' + col, np.array_equal(q[col], q_expected[col]))
    else:
        check('Q subset ' + col, np.allclose(q[col], q_expected[col], atol=1e-13, rtol=0))
for label, source, key in [('MD', full, 'target'), ('full', full, 'mean'), ('intra', intra, 'mean')]:
    column = label + '_Hchi_nats'
    check(column + ' exact frozen values', np.allclose(all66[column], source[key][indices, 1], atol=1e-13, rtol=0))
for label in ['full', 'intra', 'static', 'AA_balanced']:
    difference = all66[label + '_Hchi_nats'] - all66.MD_Hchi_nats
    check(label + ' signed error', np.allclose(difference, all66[label + '_signed_error_nats'], atol=1e-13, rtol=0))
    check(label + ' absolute error', np.allclose(np.abs(difference), all66[label + '_absolute_error_nats'], atol=1e-13, rtol=0))
schema = json.loads(find('static_corrected/training_schema.json').read_text(encoding='utf8'))
entry = next(r for r in schema['records'] if r['sample_id'] == 'A003R')
local = np.arange(entry['start'], entry['stop'])[indices]
for label, suffix in [('static', 'static_sequence_complex_oof.npz'), ('AA_balanced', 'aa_only_oof.npz')]:
    values = np.load(find('static_corrected/' + suffix), allow_pickle=False)['prediction'][local, 1]
    check(label + ' exact source-array values', np.allclose(values, all66[label + '_Hchi_nats'], atol=1e-13, rtol=0))
windows = pd.read_csv(find('joint_entropy_original_five_windows.csv'))
for _, row in q.iterrows():
    values = windows[windows.residue == 'Q' + str(row.residue_id)].Hchi_nats.to_numpy()
    check('Q' + str(row.residue_id) + ' five windows', len(values) == 5)
    check('Q' + str(row.residue_id) + ' target mean reproduced', abs(values.mean() - row.MD_Hchi_nats) < 2e-7)
    check('Q' + str(row.residue_id) + ' sample variance reproduced', abs(values.var(ddof=1) - row.label_window_sample_variance_nats2) < 2e-8)
check('Q40 largest Full error among 66', int(all66.loc[all66.full_absolute_error_nats.idxmax(), 'residue_id']) == 40)
check('Q40 largest GLN label variance', int(q.loc[q.label_window_sample_variance_nats2.idxmax(), 'residue_id']) == 40)
check('Q40 Full above both type means', result['Q40']['full_Hchi_nats'] > max(
      result['training_Q_means']['unweighted_residue_mean_nats'], result['training_Q_means']['source_system_record_balanced_Q_mean_nats']))
payload = {'status': 'passed', 'checks_passed': len(checks), 'checks': checks,
           'scope': 'Separate reload and arithmetic verification of all frozen model values, complete site coverage, errors, five-window summaries, source hashes, and core one-site assertions.',
           'outputs_sha256': {p.name: sha(p) for p in sorted(ROOT.iterdir())
                              if p.is_file() and p.suffix in ['.csv', '.json', '.py'] and p.name != 'verification_audit.json'}}
(ROOT / 'verification_audit.json').write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf8')
print(json.dumps({'status': payload['status'], 'checks_passed': len(checks)}, ensure_ascii=False))
