from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
ROOT=Path(__file__).resolve().parent
read=lambda p:json.loads(p.read_text(encoding='utf8'))
write=lambda p,v:p.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf8')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
p=ROOT/'fixed_analysis_protocol.json';b=ROOT/'protocol_freeze_binding.json'
protocol=read(p);binding=read(b)
write(ROOT/'fixed_analysis_protocol_initial_v1.json',protocol)
write(ROOT/'protocol_freeze_binding_initial_v1.json',binding)
time=datetime.now(timezone.utc).isoformat()
source=ROOT.parent/'experimental_labels/raw/pettmann2021_fullText.xml'
protocol['label_contract']['Pettmann_primary_fit']='Original-study reported hybrid: use constrained KD if constrained KD>20 microM, otherwise Bmax-fitted KD. This is the published Methods rule, not a model-performance choice.'
protocol['label_contract']['Pettmann_fit_sensitivities']=['all-constrained (HERMES compatibility)','all-Bmax-fitted']
protocol['label_contract']['measurement_batch_status']='Study/protocol/37C equivalence confirmed in source curation, but per-candidate shared-run/batch IDs are not supplied. Do not describe these as known same-batch measurements.'
protocol['label_contract']['fit_not_independent_samples']=True
protocol['label_contract']['study_sensitivity']='Report each 1G4 study separately in addition to candidate-level median within-study delta. Do not use study identity as a predictive feature or increase independent n for multiple study rows.'
protocol['amendments']=[dict(timestamp_utc=time,kind='complete_predeclared_primary_source_label_rule',new_model_outputs_seen=False,
    reason='Primary source retrieved and Methods explicitly selects hybrid KD. Initial protocol required primary source governs and had not selected a fit version.',
    source_path=str(source),source_sha256=sha(source),rule_quote='For constrained KD above 20 μM, we reported the constrained KD, otherwise we use the Bmax fitted KD.',
    initial_protocol_sha256=binding['files']['fixed_analysis_protocol.json'])]
write(p,protocol)
binding['amendments']=protocol['amendments']
binding['files']={x.name:sha(x) for x in [p,ROOT/'candidate_splits.json']}
binding['label_rule_completed_at_utc']=time
write(b,binding)
print(json.dumps(dict(status='completed_before_new_panel_model_outputs',protocol_sha256=sha(p),split_sha256=sha(ROOT/'candidate_splits.json'),amended_at_utc=time),ensure_ascii=False))
