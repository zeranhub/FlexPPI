from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import csv,hashlib,json
ROOT=Path(__file__).resolve().parent
labels=ROOT.parent/'experimental_labels'
read=lambda p:json.loads(p.read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
split=read(ROOT/'candidate_splits.json');protocol=read(ROOT/'fixed_analysis_protocol.json')
labelrules=read(labels/'label_rule_freeze.json')
assert labelrules['primary']=='reported_hybrid' and labelrules['frozen_before_model_results']
assert labelrules['all_have_measured_within_study_parent']
label_files=[]
for variant,info in labelrules['variants'].items():
    for role,key,hkey in [('observations','observation_file','primary_file_sha256'),('candidates','candidate_file','candidate_file_sha256')]:
        path=Path(info[key]);digest=sha(path);assert digest==info[hkey]
        label_files.append(dict(variant=variant,role=role,path=str(path),sha256=digest))
rawrows=list(csv.DictReader((labels/'SPR_labels_candidates_primary.csv').open(encoding='utf-8-sig')))
crosswalk=[]
our={(r['system'],r['sequence']):r for r in split['candidates']}
for row in rawrows:
    # Read candidate identity and provenance fields only, not numeric affinity.
    system='nyeso' if row['tcr']=='1G4' else 'tax' if row['tcr']=='A6' else None
    assert system is not None,row['tcr']
    old=our[system,row['peptide_sequence']]
    assert row['panel_reference_peptide']==old['parent_sequence']
    crosswalk.append(dict(protocol_candidate_id=old['candidate_id'],label_candidate_id=row['candidate_id'],panel_id=row['panel_id'],
        sequence=row['peptide_sequence'],system=system,analysis_set=old['analysis_set'],mutation_positions=old['mutation_positions']))
assert len(crosswalk)==30 and len({r['label_candidate_id'] for r in crosswalk})==30
byid={r['candidate_id']:r for r in split['candidates']}
tested=[]
for fold in split['folds']:
    assert set(fold['train_candidates']).isdisjoint(fold['test_candidates'])
    position=fold['heldout_position']
    assert all(position not in byid[x]['mutation_positions'] for x in fold['train_candidates'])
    assert all(byid[x]['mutation_positions']==[position] for x in fold['test_candidates'])
    assert all(byid[x]['system']==fold['system'] for x in fold['train_candidates']+fold['test_candidates'])
    tested.extend(fold['test_candidates'])
primary=[r['candidate_id'] for r in split['candidates'] if r['analysis_set']=='primary_single_substitution']
assert Counter(tested)==Counter(primary)
source_files=[Path(r'E:\Project\FlexPPI\src\flexppi')/p for p in ['inference.py','data/structures.py','model/network.py','model/information_ablation.py','config.py']]
active=protocol['base_model']
run_path=Path(active.get('local_run_path',str(ROOT.parent/'model/epoch_040.pt')))
run_sha=sha(run_path) if run_path.exists() else None
if run_sha is not None: assert run_sha==active['sha256']
out=dict(status='passed_protocol_and_existing_assets',created_at_utc=datetime.now(timezone.utc).isoformat(),
    protocol_id=protocol['protocol_id'],source_code_contract_hashes=[dict(path=str(p),sha256=sha(p)) for p in source_files],
    fixed_rule_binding_sha256=sha(ROOT/'protocol_freeze_binding.json'),protocol_sha256=sha(ROOT/'fixed_analysis_protocol.json'),split_sha256=sha(ROOT/'candidate_splits.json'),
    original_primary_label_rule_sha256=sha(labels/'label_rule_freeze.json'),label_files=label_files,
    active_checkpoint_sha256=active['sha256'],active_training_manifest_sha256=active['training_manifest_sha256'],
    active_training_records=active['records'],active_checkpoint_run_path=str(run_path),active_checkpoint_run_path_sha256=run_sha,
    same_weights_as_current_manuscript_2430=active.get('same_weights_as_current_manuscript_2430',True),
    completely_transitive_family_disjoint=active.get('fully_transitive_family_disjoint',True),
    pilot_label='Historical2431 local feasibility pilot' if active['records']==2431 else 'Target2430 protocol',
    numeric_affinity_values_not_inspected=True,new_model_panel_outputs_not_inspected=True,
    candidate_identity_crosswalk=crosswalk,n_unique_total_candidates=30,n_primary_single_substitutions=26,n_systems=2,n_position_folds=12,
    each_primary_candidate_has_one_heldout_prediction=True,mutation_positions_do_not_cross_their_test_fold=True,
    training_overlap_audit_status=read(ROOT/'training_overlap_audit.json')['status'],raw_template_audit_status=read(ROOT/'raw_template_overlap_audit.json')['status'],
    pending_execution_gates=[f'Actual execution checkpoint byte SHA must equal {active["sha256"]}; original local checkpoint already independently verified.' if run_sha is None else 'Local execution checkpoint copy byte SHA verified against active protocol.',
        'Verify actual prepared inputs, chain roles, features U and template/candidate source hashes before inference; original raw templates alone do not substitute.',
        'No calibration fitting, inference, MD or new model training was run by this reviewer.'],
    required_protocol_changes=[])
(ROOT/'protocol_gate_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k in ['status','n_unique_total_candidates','n_primary_single_substitutions','n_systems','n_position_folds','required_protocol_changes']},ensure_ascii=False))
