from pathlib import Path
import hashlib,json
import numpy as np
root=Path(__file__).resolve().parent
source=root/'snapshots/initial_7_entries'
out=root/'evaluation'
names=['Full','Intra','Static','AA','DynaProt-M']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
audit=json.loads((source/'scoring_audit.json').read_text())
assert set(names)<=set(audit['complete_methods'])
with np.load(source/'aligned_predictions.npz') as f:
    np.savez_compressed(out/'frozen_predictions.npz',**{name:f[name] for name in names})
metadata={name:audit['methods'][name] for name in names}
proof=dict(cache_sha256=sha(out/'frozen_predictions.npz'),truth_sha256=audit['truth_sha256'],
    original_scoring_audit_sha256=sha(source/'scoring_audit.json'),
    original_scoring_script_sha256=audit['script_sha256'],
    original_prediction_manifest_sha256=sha(root/'cohort/existing_prediction_manifest.csv'),methods=metadata)
(out/'frozen_predictions_audit.json').write_text(json.dumps(proof,indent=2),encoding='utf-8')
print('Prepared identity-verified frozen cache:',names)
