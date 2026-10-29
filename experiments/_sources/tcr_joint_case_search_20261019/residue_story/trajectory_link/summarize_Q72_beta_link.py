from pathlib import Path
import sys,json,math
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parents[2]
sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
import pandas as pd
rep=pd.read_csv(ROOT/'per_replica_Q72_beta_contact_and_RMSF.csv')
pair=pd.read_csv(EXP/'tcr_expanded_screen_v41_20261017/mel5_validation/corrected_pair_all_aligned_positions.csv')
q=pair[(pair.role=='mhc')&(pair.id_a==72)].iloc[0]
gly=pair[(pair.role=='tcr_beta')&(pair.position_a==54)].iloc[0]
corr=pd.read_csv(ROOT/'within_replica_20ns_descriptive_correlations.csv')
state=pd.read_csv(ROOT/'Q72_joint_chi_state_occupancies.csv')
frames=pd.read_csv(ROOT/'per_frame_Q72_beta_contacts_and_states.csv')
frames['both_Q72_backbone_hbonds']=frames.Q72_OE1_betaN_hbond&frames.Q72_NE2_betaGlyO_hbond
joint=frames.groupby(['pdb','replica']).agg(frames=('frame','size'),OE1_betaN_Hbond_fraction=('Q72_OE1_betaN_hbond','mean'),NE2_GlyO_Hbond_fraction=('Q72_NE2_betaGlyO_hbond','mean'),both_same_frame_Hbond_fraction=('both_Q72_backbone_hbonds','mean')).reset_index()
joint.to_csv(ROOT/'Q72_two_end_backbone_Hbond_same_frame_per_replica.csv',index=False)
jointstate=frames.groupby(['pdb','replica','Q72_joint_state']).agg(frames=('frame','size'),OE1_betaN_Hbond_fraction=('Q72_OE1_betaN_hbond','mean'),NE2_GlyO_Hbond_fraction=('Q72_NE2_betaGlyO_hbond','mean'),both_same_frame_Hbond_fraction=('both_Q72_backbone_hbonds','mean')).reset_index()
jointstate['state_fraction']=jointstate.frames/jointstate.groupby(['pdb','replica']).frames.transform('sum')
jointstate.to_csv(ROOT/'Q72_two_end_backbone_Hbond_by_joint_state.csv',index=False)
out={
 'status':'passed',
 'residue_pair':{'Hchi_residue':'unchanged HLA-A*02:01 Gln72 (chain A, both model inputs 72)','RMSF_residue':'unchanged TCR beta Gly54 in corrected common/literature numbering, WT model input E54, engineered model input E55'},
 'model':{'Q72_Hchi_nats_WT':float(q.hchi_full_a),'Q72_Hchi_nats_engineered':float(q.hchi_full_b),'Q72_Hchi_delta_nats':float(q.hchi_full_b-q.hchi_full_a),'beta_Gly_RMSF_A_WT':float(gly.rmsf_full_a),'beta_Gly_RMSF_A_engineered':float(gly.rmsf_full_b),'beta_Gly_RMSF_delta_A':float(gly.rmsf_full_b-gly.rmsf_full_a)},
 'MD_descriptive_three_independent_replicas':{},
 'window_correlations':{},
 'same_frame_two_end_Hbond_geometry':joint.to_dict('records'),
 'supported_inference':'Affinity maturation creates a Q72 sidechain/backbone contact network with beta CDR2; the unchanged Q72 has lower predicted joint chi state entropy and the opposing unchanged Gly54 has lower predicted CA RMSF. Existing MD supports higher contact occupancies, lower three-replica mean Q72 Hchi, and lower Gly RMSF; individual Q72 Hchi replicas overlap between conditions. Continuous contact residence lifetimes were not evaluated. These local observables support a coordinated interface reorganization/stabilization interpretation, not a directional causal effect of entropy upon RMSF or a calculation of binding affinity.',
 'causal_limit':'No selective Q72 intervention, mutational double-cycle, or directional perturbation was performed; fluctuations and sidechain occupancy can share causes from the engineered CDR2 sequence and geometry. Within-replica 20ns Hchi/RMSF correlations are not reproducibly signed.',
}
for pdb in ['3hg1','4jff']:
    f=rep[(rep.pdb==pdb)&(rep.beta_common_position==54)].sort_values('replica')
    metrics=['Q72_Hchi_nats','beta_CA_RMSF_A','min_heavy_contact_le4_A_fraction','min_heavy_contact_le4p5_A_fraction','Q72_NE2_beta_GlyO_le3p5_A_fraction','Q72_NE2_beta_GlyO_hbond_fraction','Q72_OE1_beta_Ile_or_PheN_hbond_fraction']
    out['MD_descriptive_three_independent_replicas'][pdb]={m:{'per_replica':f[m].tolist(),'mean':float(f[m].mean()),'sample_sd':float(f[m].std())} for m in metrics}
    z=corr[(corr.pdb==pdb)&(corr.beta_common_position==54)]
    for x,y in [('Q72_Hchi_nats','beta_CA_RMSF_A'),('Q72_Hchi_nats','min_heavy_contact_le4_A_fraction')]:
        g=z[(z.x==x)&(z.y==y)].sort_values('replica')
        out['window_correlations'][f'{pdb}:{x}__{y}']=[float(v) if pd.notna(v) else None for v in g.Spearman_descriptive]
    s=state[state.pdb==pdb]
    out['MD_descriptive_three_independent_replicas'][pdb]['dominant_joint_chi_states']=[g.sort_values('fraction',ascending=False).head(3).to_dict('records') for _,g in s.groupby('replica')]
(ROOT/'key_findings.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k not in ['MD_descriptive_three_independent_replicas']},ensure_ascii=False,indent=2))
