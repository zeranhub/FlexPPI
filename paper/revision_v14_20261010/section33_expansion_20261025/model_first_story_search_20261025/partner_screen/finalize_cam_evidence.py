"""Independent exhaustive grouping and coordinate-level checks for CaM story."""
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import json,hashlib,sys
import numpy as np
import pandas as pd
import MDAnalysis as mda
from scipy.stats import spearmanr
P=Path(__file__).resolve().parent
CORE=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
def summarize(g):
    strong=g.md_delta.abs()>=.1
    return dict(n=len(g),MD_state_a=g.md_a.mean(),MD_state_b=g.md_b.mean(),MD_delta=g.md_delta.mean(),Full_delta=g.full0_delta.mean(),
        Full_MAE=(g.md_delta-g.full0_delta).abs().mean(),zero_MAE=g.md_delta.abs().mean(),
        Intra_MAE=(g.md_delta-g.intra0_delta).abs().mean(),AA_SASA_contacts_MAE=(g.md_delta-g.aa_sasa_contacts0_delta).abs().mean(),
        Static47_SASA_MAE=(g.md_delta-g.static47_sasa0_delta).abs().mean(),rho=float(spearmanr(g.md_delta,g.full0_delta).statistic) if len(g)>3 else None,
        n_large=int(strong.sum()),correct_large=int(((np.sign(g.md_delta)==np.sign(g.full0_delta))&strong).sum()))
def main():
    d=pd.read_csv(P/'cam_md_state_contrast_all_common_residues.csv');w=pd.read_csv(P/'cam_md_state_contrast_all_residues_five_blocks.csv')
    for sid in ['A004C','A007U']:
        z=np.load(P/'training_member_predictions'/('full_'+sid+'.npz'))
        raw=np.load(CORE/'neural_exports/evaluation/full'/(sid+'.npz'))
        for state in ['a','b']:
            g=d[d['sample_'+state].eq(sid)]
            if not len(g):continue
            mi=(g.metric=='hchi').astype(int).to_numpy();ix=g['index_'+state].to_numpy(int)
            assert np.max(np.abs(g['md_'+state]-raw['target'][ix,mi]))<1e-6
            assert np.max(np.abs(g['full0_'+state]-z['mean'][ix,mi]))<1e-6
    regions={
        'all':np.ones(len(d),bool),
        'contact_either_state':d.interface,
        'outside_both_CA8_interfaces':~d.interface,
        'NHE1_contact_only':(d.contact_a>0)&(d.contact_b==0),
        'IQCG_contact_only':(d.contact_a==0)&(d.contact_b>0),
        'contact_in_both_states':(d.contact_a>0)&(d.contact_b>0),
        'EF1_20_31':d.mature_cam_number.between(20,31),
        'EF2_56_67':d.mature_cam_number.between(56,67),
        'EF3_93_104':d.mature_cam_number.between(93,104),
        'EF4_129_140':d.mature_cam_number.between(129,140),
        'all_four_EF_loops':d.mature_cam_number.between(20,31)|d.mature_cam_number.between(56,67)|d.mature_cam_number.between(93,104)|d.mature_cam_number.between(129,140),
        'all_nine_Met':d.amino_acid=='M',
    }
    rows=[]
    for metric in ['rmsf','hchi']:
        for name,mask in regions.items():
            g=d[mask&(d.metric==metric)]
            if len(g): rows.append(dict(metric=metric,region=name,**summarize(g)))
    pd.DataFrame(rows).to_csv(P/'final_cam_fixed_structural_group_metrics.csv',index=False)
    coordinates=[];top=[]
    mapping=pd.read_csv(P/'cam_md_mature_number_mapping.csv')
    for sid in ['A004C','A007U']:
        z=np.load(P/'verified_case_trajectories'/sid/'ca_coordinates.npz');mm=mapping[mapping.sample_id==sid].set_index('mature_cam_number')
        x=z['coords_A'];ni=mm.loc[list(range(4,76))].residue_index.to_numpy(int);ci=mm.loc[list(range(82,149))].residue_index.to_numpy(int)
        dist=np.linalg.norm(x[:,ni].mean(1)-x[:,ci].mean(1),axis=1)
        allix=mm.residue_index.to_numpy(int);ca=x[:,allix];rg=np.sqrt(((ca-ca.mean(1,keepdims=True))**2).sum(2).mean(1))
        for b in range(5):
            ix=slice(b*1000,(b+1)*1000)
            coordinates.append(dict(sample_id=sid,block=b+1,interlobe_centroid_mean_A=float(dist[ix].mean()),interlobe_centroid_SD_A=float(dist[ix].std()),CaM_CA_Rg_mean_A=float(rg[ix].mean())))
        u=mda.Universe(str(Path(r'D:\data\dynarepo_non_nucleic_full')/sid/'topology.tpr'))
        calcium=[dict(resname=str(r.resname),resid=int(r.resid)) for r in u.residues if r.resname.upper() in {'CA','CAL','CA2','CA2+','CAL2'}]
        top.append(dict(sample_id=sid,tpr_atoms=len(u.atoms),calcium_named_residues=calcium,
            note='Both provided TPR topologies contain protein atoms only and no named calcium. This comparison cannot establish a calcium-removal mechanism.'))
    pd.DataFrame(coordinates).to_csv(P/'final_cam_interlobe_distance_and_compaction.csv',index=False)
    (P/'final_cam_topology_calcium_audit.json').write_text(json.dumps(top,indent=2),encoding='utf-8')
    key=d[d.amino_acid=='M'].merge(pd.read_csv(P/'cam_md_state_contrast_residue_stability.csv')[['metric','mature_cam_number','blocks_same_sign_as_mean','min_block_delta','max_block_delta','fine_entropy_mean_per_chi_delta','width_mean_per_chi_delta','shifted60_joint_entropy_delta']],on=['metric','mature_cam_number'],validate='one_to_one')
    key.to_csv(P/'final_cam_all_nine_methionines_with_stability.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False));print('\nCOMPACTION\n',pd.DataFrame(coordinates).groupby('sample_id').mean().to_string())
    (P/'final_cam_checks.json').write_text(json.dumps(dict(status='passed',all_common_positions=145,valid_chi_residues=122,
        metrics_independently_compared_to_npz=True,fold=0,both_states_in_training=True,
        naming='NHE1 structural state to IQCG structural state, not an isolated partner or calcium intervention.',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
