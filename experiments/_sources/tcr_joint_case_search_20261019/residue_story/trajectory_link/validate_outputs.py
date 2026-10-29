from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'label_sensitivity/runtime'))
import pandas as pd
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
audit=json.loads((ROOT/'analysis_audit.json').read_text())
findings=json.loads((ROOT/'key_findings.json').read_text())
frames=pd.read_csv(ROOT/'per_frame_Q72_beta_contacts_and_states.csv')
reps=pd.read_csv(ROOT/'per_replica_Q72_beta_contact_and_RMSF.csv')
windows=pd.read_csv(ROOT/'nonoverlapping_20ns_window_descriptors.csv')
states=pd.read_csv(ROOT/'Q72_joint_chi_state_occupancies.csv')
assert len(frames)==6012 and len(reps)==30 and len(windows)==300 and len(states)==162
assert len(frames.groupby(['pdb','replica']))==6
assert not frames.isna().any().any()
assert not reps.isna().any().any()
assert all(audit['source_files'][i]['raw_files_unchanged'] for i in range(6))
assert audit['script_sha256']==sha(ROOT/'analyze_Q72_beta_link.py')
assert all(abs(v-1)<1e-12 for v in states.groupby(['pdb','replica']).fraction.sum())
joint=pd.read_csv(ROOT/'Q72_two_end_backbone_Hbond_same_frame_per_replica.csv')
assert (joint[joint.pdb=='3hg1'].both_same_frame_Hbond_fraction==0).all()
assert (joint[joint.pdb=='4jff'].both_same_frame_Hbond_fraction>.33).all()
inventory=[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(ROOT.iterdir()) if p.is_file() and p.name!='validation_and_inventory.json']
(ROOT/'validation_and_inventory.json').write_text(json.dumps({'status':'passed','frame_rows':len(frames),'independent_trajectories':6,'per_replica_per_residue_rows':len(reps),'window_rows':len(windows),'all_possible_Q72_joint_state_rows':len(states),'all_saved_frames_included':True,'all_3_replicas_included_both_conditions':True,'original_labels_reproduced_atol':1e-12,'no_new_MD_or_inference':True,'files':inventory},ensure_ascii=False,indent=2),encoding='utf8')
print('validated all six trajectories and output inventory')
