"""Append-only technical interpretation; no feature/statistical rule changes."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
protocol=ROOT/'structure_preparation_protocol.json'
script=ROOT/'prepare_complete_structures.py'
fixer=ROOT/'tool_runtime/pdbfixer/pdbfixer.py'
feature=ROOT/'feature_extraction/extract_features.py'
text=fixer.read_text(encoding='utf8')
assert 'system.setParticleMass(atom.index, 0.0)' in text
assert 'newPositions2[a2.index] = state.getPositions()[a1.index]' in text
result=dict(status='reviewed_technical_details_without_rule_changes',reviewed_at_utc=datetime.now(timezone.utc).isoformat(),
    inputs={str(p):sha(p) for p in [protocol,script,fixer]},
    algorithm='PDBFixer1.12.0 addMissingAtoms seed17 and OpenMM8.6.1 CPU. Missing-residue list empty; no added loops, no added hydrogen step, no global relaxation.',
    original_atom_preservation='PDBFixer assigns existing protein atoms mass0 in its optimization system and writes back only the newly added atom coordinates. Actual source protein heavy coordinates are independently compared at the stored PDB0.001Angstrom precision.',
    internal_clash_repair='PDBFixer may use at most10x200 Langevin steps at5fs for newly added atoms when contacts are below0.13nm, followed by local minimization. Original atoms remain fixed. This is an atom-completion procedure, not production protein MD or an additional trajectory experiment.',
    U_reference_clarification='For each of5 templates, the panel parent sequence goes through the same preparation as all candidates; its completed heavy-atom structure defines the unchanged8Angstrom HLA/TCR environment. Under2BNR, the panel parent remains9V, not the template-native9C peptide. No candidate-specific U is selected.',
    static_source_requirement='Use the explicit static_pdb/output path tied to preparation_audit.output_sha256. The inherited prepared_input_pdb points to the earlier deleted-sidechain neural-input file and must not be used for SASA or U when a complete static file is present.',
    input_water_policy='PDBFixer removes heterogens/water before atom completion. The invariant and static calculations concern canonical protein heavy atoms; source crystallographic waters are not counted as model or SASA protein atoms.',
    dynamic_model_input_invariant='Residue identities, author IDs, chain ordering, N/CA/C/O coordinates and masks must equal the sequence-synchronized fixed-backbone input. Sidechain optimization does not add model input channels.',
    feature_or_statistical_rule_changed=False,analysis_protocol_sha256=sha(HERE/'fixed_analysis_protocol.json'),
    candidate_split_sha256=sha(HERE/'candidate_splits.json'),
    candidate_SPR_values_read=False,model_predictions_read=False,neural_inference_or_training_performed=False,
    all76_final_static_QC_pending=True,frozen_U_independent_review_pending=True)
output=HERE/'structure_preparation_technical_binding.json'
if output.exists():raise RuntimeError('Preserve prior technical review; version a subsequent amendment explicitly')
output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(dict(status=result['status'],analysis_protocol_sha256=result['analysis_protocol_sha256'],statistical_rules_changed=False),ensure_ascii=False))
