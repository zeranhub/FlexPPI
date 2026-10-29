"""Complete source-based prespecification before seeing new model outcomes."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib
HERE=Path(__file__).resolve().parent
p=HERE/'multimode_protocol.json'
a=json.loads(p.read_text(encoding='utf8'))
archive=HERE/'multimode_protocol_initial_contract.json'
assert not archive.exists(),'Retain immutable history; version a later amendment explicitly.'
archive.write_bytes(p.read_bytes())
a['status']='analysis_rules_frozen_pending_label_only_eligibility_binding'
a['frozen_at_utc']=datetime.now(timezone.utc).isoformat()
a['source_based_completion']={
 'new_predictions_seen':False,'new_functional_label_values_seen_by_evaluator':False,
 'reason':'Seven HERMES/MSKCC complete peptide saturation panels include published HERMES-fixed and ProteinMPNN scores. These are stronger, feasible static comparators than rebuilding all1200 candidate sidechains for two geometric summaries.',
 'not_changed_from_old_pilot':'The old20261021 study remains unchanged and negative. This contract defines a separately labelled exploratory study.'}
a['candidate_features']['static']=['HERMES_fixed_score_candidate_minus_parent','ProteinMPNN_score_candidate_minus_parent']
a['candidate_features']['static_source_contract']='Select the original published fixed-backbone score column for each method from documented source code, not by correlation with EC50. Preserve raw score, score definition, orientation and parent subtraction; no AF3-relaxed/closest-template or target-outcome-selected column. Both external methods also get standalone zero-shot ranking comparisons using their published direction.'
a['candidate_features']['static_unavailable']='If either predeclared published comparator is absent for candidates, record coverage failure; do not silently replace with a weaker method. Any restricted cohort is defined before new prediction/label evaluation and uses the same candidate IDs for all policies.'
a['mode1_site_nomination']['strategies']=['random_positions','static_contact_degree','static_relative_SASA','parent_R_percentile','parent_H_percentile','parent_joint_percentile','parent_joint_low','parent_partner_context_response']
a['mode1_site_nomination']['primary_strategy']='parent_joint_percentile'
a['mode1_site_nomination']['companion_low']='parent_joint_low=1-parent_joint_percentile. It is a predeclared competing hypothesis, never substituted for the primary on the basis of performance.'
a['mode1_site_nomination']['companion_partner_context']='If feasible, remove the binding partner chains while keeping the parent target chain identities/backbone/residue order and using the same model. Site score is mean of abs(R_full-R_partner_removed)/training_R_sameAA_RMSscale and abs(H_full-H_partner_removed)/training_H_sameAA_RMSscale over defined metrics; G/A use R only. Higher indicates computational partner-context sensitivity, not real apo dynamics or a binding experiment.'
a['mode1_site_nomination']['parent_H_only_Gly_Ala']='JointHchi undefined at Gly/Ala. In H-only ranking assign the predeclared neutral percentile0.5, flag undefined, and retain equal candidate/site scope; do not invent physical Hchi. Joint-high/low use defined R alone for these sites.'
a['mode1_site_nomination']['site_companions_role']='Report every companion and both budgets separately; a companion win cannot establish primary joint-high success or be pooled with a different mode.'
a['data_contract']['new_antigen_source']='Seven distinct TCR saturation-peptide EC50 panels from pinned HERMES mutation_effects/mskcc and original avidity_curves.xlsx/infer_ec50s source. Verify receptor relationships, assay construct and all candidate coverage before treating as independent families.'
a['data_contract']['EC50']='Orient as pEC50=-log10(EC50 in verified common molar unit) or parent-relative pEC50 if a measured parent is supplied. State T-cell activity endpoint; do not call itKD. No-response or outside-tested-dose estimates require censoring review. Do not remove all no-response candidates then claim complete-pool functional-discovery validity.'
a['data_contract']['eligibility_before_evaluation']='Curator supplies immutable per-panel candidate scope/endpoint/censoring/threshold contract and hash before evaluation. Label-derived distribution gates assess task suitability only; they do not select a better-performing model, threshold or subset.'
p.write_text(json.dumps(a,ensure_ascii=False,indent=2),encoding='utf8')
sha=lambda q:hashlib.sha256(Path(q).read_bytes()).hexdigest()
b=dict(protocol_id=a['protocol_id'],status=a['status'],frozen_at_utc=a['frozen_at_utc'],
    initial_contract_sha256=sha(archive),current_protocol_sha256=sha(p),
    no_new_model_predictions_seen=True,no_new_functional_label_values_seen=True,
    data_eligibility_binding_pending=True,feature_implementation_freeze_pending=True)
(HERE/'protocol_freeze_binding.json').write_text(json.dumps(b,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(b,ensure_ascii=False))
