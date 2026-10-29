"""Independently verify the protected 4g9f example data and stored smoothing."""
from pathlib import Path
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent.parent
ARC = ROOT / 'figure1/source_archive'
OLD = Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\figures\reference_rebuild_20261005\flexppi_main\revision_v3')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
copies = []
for name in ['curve_data_all_840.csv','curve_local_means.csv','curve_display_receipt.json']:
    assert (OLD/name).read_bytes() == (ARC/name).read_bytes()
    copies.append({'file':str(ARC/name),'original':str(OLD/name),'sha256':sha(ARC/name)})
with (ARC/'curve_data_all_840.csv').open(newline='',encoding='utf-8-sig') as f:
    raw = list(csv.DictReader(f))
with (ARC/'curve_local_means.csv').open(newline='',encoding='utf-8-sig') as f:
    smooth = list(csv.DictReader(f))
assert len(raw) == 1680 and len(smooth) == 840
positions = {(int(r['global_index']),r['chain']) for r in raw}
assert len(positions) == 840
preds = {(int(r['global_index']),r['metric']):r for r in raw}
chain_counts = {}
maxerr = 0.0
valid_counts = {'rmsf':0,'hchi':0}
for r in smooth:
    i = int(r['global_index'])
    chain = r['chain']
    chain_counts[chain] = chain_counts.get(chain,0)+1
    for metric, prefix in [('rmsf','rmsf'),('hchi','hchi')]:
        source = preds[(i,metric)]
        valid = int(source['valid']) == 1
        valid_counts[metric] += valid
        if not valid:
            assert r[prefix+'_raw'] == '' and r[prefix+'_local_mean'] == ''
            continue
        assert float(r[prefix+'_raw']) == float(source['prediction'])
        values=[]
        for j in range(max(0,i-5),min(840,i+6)):
            other=preds[(j,metric)]
            if other['chain']==chain and int(other['valid'])==1:
                values.append(float(other['prediction']))
        expected=sum(values)/len(values)
        actual=float(r[prefix+'_local_mean'])
        assert len(values)==int(r[prefix+'_window_valid_count'])
        error=abs(expected-actual)
        maxerr=max(error,maxerr)
        assert error<2e-14
assert list(chain_counts.values()) == [276,100,10,204,250]
assert valid_counts == {'rmsf':840,'hchi':734}
report = {'status':'passed','source_files':copies,'raw_rows':len(raw),'unique_positions':len(positions),
          'valid_counts':valid_counts,'chain_lengths':chain_counts,'all_smoothing_values_recomputed':True,
          'mean_method':'Centered ±5 within same chain, valid targets only; invalid central Hχ positions remain absent.',
          'maximum_smoothing_error':maxerr,'all_raw_values_and_valid_masks_unchanged':True}
(ROOT/'audit/protected_840_data.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','positions':840,'raw_rows':1680,'max_smoothing_error':maxerr}))
