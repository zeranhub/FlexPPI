from pathlib import Path
import json,sys
import compute_metrics as cm
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
df=pd.read_csv(OUT/'per_record_metrics.csv')
records=pd.read_csv(OUT/'record_codes.csv').sort_values('record_code')
draws=cm.stratified_draws(records.dataset.to_numpy())[:1000]
rows=[];counts={}
for target in cm.TARGETS:
 for method in cm.METHODS:
  g=df[(df.target==target)&(df.method==method)].sort_values('record_code')
  with np.load(OUT/f'raw_curves_{target}_{method}.npz') as z:
   code=z['record_code'];end=np.r_[np.flatnonzero(np.diff(code)!=0),len(code)-1]
   assert np.array_equal(code[end],g.record_code)
   assert np.array_equal(z['recall'][end],np.ones(len(g)))
   assert np.array_equal(z['precision'][end],g.prevalence.to_numpy()) or np.allclose(z['precision'][end],g.prevalence.to_numpy(),rtol=0,atol=2e-16)
   pp=z['precision'][end].copy()
  valid=g.eligible.to_numpy();pp[~valid]=np.nan
  interval=np.quantile(np.nanmean(pp[draws],axis=1),[.025,.975])
  rows.append(dict(target=target,method=method,recall=1.,precision=float(np.nanmean(pp)),precision_ci_low=float(interval[0]),precision_ci_high=float(interval[1]),n_records=int(valid.sum()),endpoint='all_selected_original_last_threshold'))
  counts[f'{target}/{method}']=len(g)
pd.DataFrame(rows).to_csv(OUT/'all_selected_endpoint.csv',index=False)
protocol_sha=cm.sha(OUT/'PRESPECIFIED_PROTOCOL.md')
summary=json.loads((OUT/'summary.json').read_text(encoding='utf8'))
summary['protocol_sha256']=protocol_sha
summary['display_full_selection_endpoint']={'file':'all_selected_endpoint.csv','source':'last original score-threshold point per record','added_after_grid':True,'changes_exact_AP_or_AUROC':False}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf8')
audit=json.loads((OUT/'input_alignment_audit.json').read_text(encoding='utf8'))
audit['protocol_sha256']=protocol_sha
audit['endpoint_exporter_sha256']=cm.sha(__file__)
(OUT/'input_alignment_audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf8')
out={'status':'passed','raw_curve_all_selected_endpoints_checked_per_method_target':counts,'record_total':3063,'all_selected_precision_matches_record_prevalence':True,'all_selected_recall_equals_one':True,'unchanged_exact_AP_and_AUROC':True,'exporter_sha256':cm.sha(__file__),'protocol_sha256':protocol_sha,'summary_sha256':cm.sha(OUT/'summary.json')}
(OUT/'curve_endpoint_audit.json').write_text(json.dumps(out,indent=2),encoding='utf8')
print(json.dumps(out,indent=2))
