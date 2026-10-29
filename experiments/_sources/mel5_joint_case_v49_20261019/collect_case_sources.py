"""Archive the small analysis sources; keep original trajectories at their audited paths."""
from pathlib import Path
import hashlib,json,shutil
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parent
out=ROOT/'case_sources';out.mkdir(exist_ok=True)
selection={
 'mel5_validation':('tcr_expanded_screen_v41_20261017/mel5_validation',[
   'corrected_pair_all_aligned_positions.csv','corrected_alignment_and_CDR_audit.json',
   'plot_sources/mel5_all6_CDR_per_position_source.csv','plot_sources/mel5_all6_CDR_summary.csv',
   'plot_sources/mel5_natural_pMHC_interface_source.csv','plot_sources/mel5_CDR_definitions.csv']),
 'model_context':('tcr_expanded_screen_v41_20261017/mel5_context',[
   'all_natural_pMHC_context_comparison.csv','mel5_context_inference_audit.json']),
 'primary_literature':('tcr_expanded_screen_v41_20261017/biology_sources',[
   'mel5_thermodynamics_2013.md']),
 'model_training_audit':('itrace_case_search_only_20261017/model_audit',['verified_summary.json']),
 'original_predictions':('tcr_case_v40_20261017/family_predictions/full',['itrace_3hg1.npz','itrace_4jff.npz']),
 'residue_story':('tcr_joint_case_search_20261019/residue_story',[
   'paired_Q72_beta_CDR2_static_and_model.json','root_geometry.py','completion_status.json',
   'independent_review/final_evidence_chain_review.json',
   'trajectory_link/analysis_audit.json','trajectory_link/key_findings.json',
   'trajectory_link/per_replica_Q72_beta_contact_and_RMSF.csv','trajectory_link/three_replica_summary.csv',
   'trajectory_link/Q72_joint_chi_state_occupancies.csv',
   'trajectory_link/Q72_two_end_backbone_Hbond_same_frame_per_replica.csv',
   'trajectory_link/Q72_two_end_backbone_Hbond_by_joint_state.csv',
   'trajectory_link/per_frame_Q72_beta_contacts_and_states.csv',
   'trajectory_link/nonoverlapping_20ns_window_descriptors.csv',
   'trajectory_link/within_replica_20ns_descriptive_correlations.csv',
   'trajectory_link/analyze_Q72_beta_link.py','trajectory_link/validate_outputs.py',
   'trajectory_link/validation_and_inventory.json'])
}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rows=[]
for label,(folder,files) in selection.items():
    for rel in files:
        source=EXP/folder/rel;assert source.is_file(),source
        dest=out/label/rel;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,dest)
        assert sha(source)==sha(dest)
        rows.append({'original_path':str(source),'archive_path':dest.relative_to(ROOT).as_posix(),
                     'sha256':sha(source),'bytes':source.stat().st_size})
(ROOT/'case_source_inventory.json').write_text(json.dumps({'status':'passed','files':rows,
 'raw_trajectories':'Kept in the original locations recorded and hashed in trajectory_link/analysis_audit.json; no new MD production.',
 'model_checkpoint_SHA256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5'},
 ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'status':'passed','copied_files':len(rows),'total_bytes':sum(r['bytes'] for r in rows)}))
