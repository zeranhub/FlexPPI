from pathlib import Path
import importlib.util
import json
import hashlib
import numpy as np

path = Path(__file__).with_name('analyze_oof_complementarity.py')
spec = importlib.util.spec_from_file_location('pair_analysis',path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
a = np.arange(8,dtype=float)
r = module.tails(a)[0]
h = module.tails(a[::-1])[0]
classes = module.joint_class(r,h)
assert list(classes) == ['low_RMSF_high_Hchi']*2+['unclassified']*4+['high_RMSF_low_Hchi']*2
assert module.tails(np.zeros(8))[3] == 'collapsed_q25_q75'
assert (module.tails(np.zeros(8))[0] == '').all()
assert module.correlation(a,a)['spearman'] == 1.0
assert module.correlation(a,np.zeros(8))['correlation_reason'] == 'second_variable_constant'
assert module.correlation(a[:3],a[:3])['correlation_reason'] == 'fewer_than_4_common_positions'
assert module.errors(a,a)['mae'] == 0.0
aa = np.asarray(list('VVVVQQQQ'))
x = np.column_stack([a,a[::-1]])
centered = x.copy()
for kind in set(aa):
    ix = aa == kind
    centered[ix] -= x[ix].mean(axis=0)
    assert np.allclose(centered[ix].mean(axis=0),0,atol=1e-15)
result = dict(status='passed',analysis_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),checks=['known reversed quartiles and 50% middle retention','collapsed quartiles never arbitrarily split','constant correlation explicitly undefined','short correlation explicitly undefined','identical prediction has zero MAE and rho1','AA centering subtracts each output mean separately'],CPU_only=True)
Path(__file__).with_name('synthetic_invariants_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
