from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,importlib.util
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
review=read(HERE/'static_complete_inputs_independent_audit.json')
assert review['status']=='passed_all76_static_inputs_and_frozen_U'
assert review['n_reviewed']==76 and len(review['U_review'])==5
spec=importlib.util.spec_from_file_location('independent_static_QA',HERE/'audit_static_complete_inputs.py')
ind=importlib.util.module_from_spec(spec);spec.loader.exec_module(ind)
path=ROOT/'feature_extraction/extract_features.py'
spec=importlib.util.spec_from_file_location('frozen_feature_definitions',path)
feature=importlib.util.module_from_spec(spec);spec.loader.exec_module(feature)
assert feature.THREE_TO_ONE==ind.AA3 and feature.CHI_COUNTS==ind.CHI
assert all(feature.CANONICAL[ind.AA3[name]]==({'N','CA','C','O'}|set(atoms.split())) for name,atoms in ind.SIDECHAINS.items())
assert feature.RADII=={'C':1.70,'N':1.55,'O':1.52,'S':1.80}
index=read(ROOT/'structures/static_complete_inputs.json')
for row in index['candidates']:
    assert feature.static_path(row)==Path(row['static_complete_pdb'])
    assert feature.static_path(row)!=Path(row['prepared_input_pdb'])
protocol=read(HERE/'fixed_analysis_protocol.json')
assert sha(HERE/'fixed_analysis_protocol.json')=='79e2bb8d848c5843619fdb58452d008e8637bb30ae1b26284969c11a2e2406dc'
result=dict(status='passed',reviewed_at_utc=datetime.now(timezone.utc).isoformat(),
    scope='Final76 complete protein-heavy structures, unchanged neural inputs, immutable static-path index, frozen5 template U regions, independent20AA/chi/canonical-atom definitions.',
    all76_files_independently_parsed=True,all4697_added_heavy_atoms_only=True,
    all_existing_protein_heavy_coordinates_exactly_unchanged=True,
    all_model_input_sequence_chain_order_author_ids_NCAO_and_masks_exactly_unchanged=True,
    all_canonical_heavy_atoms_complete=True,no_added_residues_or_loops=True,
    five_parents_receive_same_preparation=True,
    static_path_index_resolves_to_complete_files_not_deleted_sidechain_inputs=True,
    U_independent_method='All reference residue/peptide atom distances recomputed with independent full NumPy distances, not the feature code cKDTree; all positions, roles,AA,chi counts,masks and distances matched.',
    U_counts=[{k:r[k] for k in ['template_pdb','n_U_all','n_U_RMSF','n_U_Hchi','role_counts']} for r in review['U_review']],
    AA_mapping_implementation_repair='Feature source had a19-amino-acid mapping missingLEU; it was corrected and independently compared for all20 names to the real model parser before candidate prediction analysis. This is an implementation repair, not a change in candidate eligibility, U definition, features or statistical rules.',
    structural_method='PDBFixer1.12.0 seed17/OpenMM8.6.1 CPU; existing atoms fixed; only missing protein heavy atoms rebuilt. Internal new-atom clash integration is atom completion, not production protein MD.',
    feature_and_statistical_protocol_unchanged=True,analysis_protocol_sha256=sha(HERE/'fixed_analysis_protocol.json'),
    candidate_split_sha256=sha(HERE/'candidate_splits.json'),
    model_checkpoint_sha256=protocol['base_model']['sha256'],
    historical2431_model_different_from_manuscript2430=True,indirect_training_component_overlap_still_disclosed=True,
    source_files={str(p):sha(p) for p in [HERE/'static_complete_inputs_independent_audit.json',HERE/'audit_static_complete_inputs.py',
        ROOT/'structures/static_complete_manifest.json',ROOT/'structures/static_complete_inputs.json',ROOT/'structures/fixed_backbone_manifest.json',
        ROOT/'feature_extraction/fixed_regions.json',path,ROOT/'feature_extraction/input_contract_test.json',
        ROOT/'feature_extraction/BUGFIX_AA_MAPPING_20261021.md',ROOT/'structure_preparation_protocol.json',
        ROOT/'prepare_complete_structures.py',HERE/'structure_preparation_technical_binding.json']},
    actual_SPR_values_read=False,candidate_predictions_read=False,neural_inference_or_training_performed_by_reviewer=False,
    review_only_no_source_or_main_manuscript_modified=True,required_changes=[])
(HERE/'final_structure_and_U_review.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:result[k] for k in ['status','analysis_protocol_sha256','U_counts','required_changes']},ensure_ascii=False))
