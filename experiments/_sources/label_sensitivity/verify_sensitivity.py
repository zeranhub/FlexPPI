from pathlib import Path
import hashlib,json
import pandas as pd
root=Path(__file__).resolve().parent
plan=json.loads((root/'sensitivity_plan.json').read_text(encoding='utf8'))
reports=[json.loads(p.read_text(encoding='utf8')) for p in (root/'record_audit').glob('*.json')]
assert len(reports)==83 and all(r['status']=='complete' for r in reports)
masks=pd.read_csv(root/'valid_mask_audit.csv')
geo=[(r['sample_id'],a['replica_index'],a['chi_bond_instances_gt5A']) for r in reports for a in r['replica_audits'] if a['chi_bond_instances_gt5A']]
mapping=[(r['sample_id'],a['replica_index'],a['residue_coverage']) for r in reports for a in r['replica_audits'] if a['residue_coverage']<1]
audit=dict(selected_records=plan['selected_records'],eligible_records=plan['eligible_records'],computed_records=len(reports),
   exclusions=plan['timestamp_exclusions'],common_window_ns=plan['common_end_ps']/1000,step_ps=plan['common_step_ps'],frames=plan['common_frames'],
   nonzero_long_chi_bond_diagnostics=geo,partial_residue_mapping=mapping,
   masks=masks.groupby(['variant','dataset','metric'])[['old_valid','new_valid','common','old_only','new_only']].sum().reset_index().to_dict(orient='records'),
   source_code_sha256=hashlib.sha256((root/'run_label_sensitivity.py').read_bytes()).hexdigest(),plan_sha256=hashlib.sha256((root/'sensitivity_plan.json').read_bytes()).hexdigest())
(root/'verification_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
print(json.dumps(audit,indent=2))
