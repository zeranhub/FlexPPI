"""Audit sidechain contact geometry of shortlisted true-environment matched pairs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from run_feature_baselines import heavy_structure

BASE=Path(__file__).resolve().parent
OUT=BASE/'case_candidates'
PREFERRED=[('A00QD','B'),('A00LH','C'),('A00LZ','B'),('A00L0','C'),('A00BS','B')]


def main():
    records=pd.read_csv(BASE/'cohort_records.csv').set_index('sample_id')
    ranking=pd.read_csv(OUT/'candidate_chain_ranking.csv')
    tcr=ranking.description.str.contains(r'T[ -]?cell receptor|TCR',case=False,na=False)
    excluded=ranking.loc[tcr,['sample_id','pdb_id','description']].drop_duplicates()
    excluded.to_csv(OUT/'excluded_tcr_records.csv',index=False)
    ranking=ranking.loc[~tcr].copy()
    ranking.to_csv(OUT/'non_tcr_candidate_chain_ranking.csv',index=False)
    allpairs=pd.read_csv(OUT/'successful_pairs_against_both_strong_baselines.csv')
    contexts=[];chosen=[]
    for sid,chain in PREFERRED:
        pairs=allpairs[(allpairs.sample_id==sid)&(allpairs.chain==chain)].sort_values('true_gap',ascending=False)
        used=set();selected=[]
        for _,pair in pairs.iterrows():
            ids={int(pair.low_residue_index),int(pair.high_residue_index)}
            if not ids&used:
                used|=ids;selected.append(pair.to_dict())
        chosen.extend(selected)
        trajectory,keys,atom_residue,sidechain,quality=heavy_structure(records.loc[sid,'structure_path'])
        atoms=list(trajectory.topology.atoms)
        coords=trajectory.xyz[0]*10
        for residue_index in sorted(used):
            target=(atom_residue==residue_index)&sidechain
            assert target.any()
            target_atoms=np.flatnonzero(target)
            distance=cdist(coords[target],coords)
            nearest_atom=distance.argmin(axis=0)
            atomdist=distance.min(axis=0)
            for other in range(len(keys)):
                if other==residue_index:
                    continue
                otheratoms=np.flatnonzero(atom_residue==other)
                closest=int(otheratoms[np.argmin(atomdist[otheratoms])])
                minimum=float(atomdist[closest])
                if minimum>5:
                    continue
                source_atom=int(target_atoms[nearest_atom[closest]])
                contexts.append(dict(sample_id=sid,pdb_id=records.loc[sid,'pdb_id'],target_chain=keys[residue_index][0],
                                     target_residue_id=keys[residue_index][1]+keys[residue_index][2],
                                     target_residue_index=residue_index,target_residue_name=atoms[target_atoms[0]].residue.name,
                                     neighbor_chain=keys[other][0],neighbor_residue_id=keys[other][1]+keys[other][2],
                                     neighbor_residue_index=other,neighbor_residue_name=atoms[closest].residue.name,
                                     cross_chain=keys[other][0]!=keys[residue_index][0],
                                     source_atom=atoms[source_atom].name,neighbor_atom=atoms[closest].name,
                                     minimum_sidechain_heavy_distance_A=minimum,
                                     n_neighbor_heavy_atoms_within_4A=int(np.sum(atomdist[otheratoms]<=4)),
                                     n_neighbor_heavy_atoms_within_5A=int(np.sum(atomdist[otheratoms]<=5))))
    pd.DataFrame(contexts).sort_values(['sample_id','target_chain','target_residue_index','minimum_sidechain_heavy_distance_A']).to_csv(OUT/'shortlist_sidechain_neighbor_contexts.csv',index=False)
    pd.DataFrame(chosen).to_csv(OUT/'shortlist_disjoint_residue_pairs.csv',index=False)
    shortlist=pd.concat([ranking[(ranking.sample_id==sid)&(ranking.chain==chain)] for sid,chain in PREFERRED],ignore_index=True)
    shortlist.to_csv(OUT/'shortlist_five_case_metrics.csv',index=False)
    (OUT/'shortlist_interpretation.json').write_text(json.dumps(dict(
        primary_strict_candidates=['3AAA / A00QD chain B','4DW2 / A00LH chain C'],
        secondary_two_pair_candidates=['6A77 / A00LZ chain B','2FGJ(raw metadata; chain annotation says 2FJG) / A00L0 chain C'],
        genuine_paired_but_no_chain_level_strong_baseline_advantage=['1I2M / A00BS chain B'],
        context_measurement='Minimum distance from observed target sidechain heavy atoms to every other observed residue heavy atom. Geometry only; not a verified hydrogen bond or a dynamics measurement.',
        caution='Cross-chain contact may connect antibody light/heavy chains or cap alpha/beta rather than the functional ligand. Inspect neighbor_chain labels before attributing response to antigen/myotrophin.',
        excluded_tcr_records=excluded.to_dict('records'),
        all_shortlisted_mainfold_family_clean=False,
        selection_scope='Exploratory post-training screen; not new confirmatory evidence, no significance testing.'),indent=2),encoding='utf-8')
    print(shortlist[['sample_id','pdb_id','chain','full_mae','aa_sasa_contacts_mae','static47_sasa_mae','full_typed_rho','aa_sasa_contacts_typed_rho','static47_sasa_typed_rho','interface_full_mae','interface_intra_mae','n_disjoint_superior_pairs']].to_string(index=False))


if __name__=='__main__':main()
