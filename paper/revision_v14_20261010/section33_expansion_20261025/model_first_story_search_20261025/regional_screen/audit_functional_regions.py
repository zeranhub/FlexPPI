"""Functional intervals were taken from the source papers after model screening.

All residues of each declared interval remain in the saved denominator. Point
examples are descriptions, not selected tests of general performance.
"""
from pathlib import Path
from collections import defaultdict
from difflib import SequenceMatcher
import json
import numpy as np
import pandas as pd
from screen_regions import summarize, METHODS, rho

OUT=Path(__file__).resolve().parent
AA=dict(zip(['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL'],'ARNDCQEGHILKKMMFPSTWYV'))

def original_residues(pid):
    rows=[];seen=set()
    for line in (OUT/(pid+'_original.pdb')).read_text().splitlines():
        if line[:6] not in ['ATOM  ','HETATM'] or line[12:16].strip()!='CA' or line[16:17] not in [' ','A'] or line[17:20] not in AA:continue
        key=(line[21],line[22:27].strip())
        if key in seen:continue
        seen.add(key)
        rows.append(dict(original_chain=key[0],original_author=key[1],amino_acid=AA[line[17:20]]))
    return pd.DataFrame(rows)

def main():
    d=pd.read_csv(OUT/'seven_candidate_residue_predictions.csv')
    mapping=[]
    for sid,pid in [('A004F','2ZXX'),('A007U','4LZX')]:
        original=original_residues(pid)
        for chain,q in d[d.sample_id==sid].groupby('chain'):
            q=q.sort_values('residue_index')
            raw=original[original.original_chain==chain]
            if ''.join(q.amino_acid)==''.join(raw.amino_acid):
                pairs=zip(range(len(q)),range(len(raw)))
            else:
                blocks=SequenceMatcher(None,''.join(q.amino_acid),''.join(raw.amino_acid),autojunk=False).get_matching_blocks()
                pairs=[(i,j) for z in blocks for i,j in zip(range(z.a,z.a+z.size),range(z.b,z.b+z.size))]
                assert len(pairs)==len(raw),(pid,chain,len(pairs),len(raw))
            for i,j in pairs:
                r=q.iloc[i];orig=raw.iloc[j]
                mapping.append(dict(sample_id=sid,pdb_id=pid,residue_index=r.residue_index,chain=chain,input_author=r.author_residue_id,original_author=orig.original_author,amino_acid=r.amino_acid))
    m=pd.DataFrame(mapping)
    m.to_csv(OUT/'original_structure_residue_mapping.csv',index=False)
    d=d.merge(m[['sample_id','residue_index','original_author']],on=['sample_id','residue_index'],how='left',validate='one_to_one')
    d['original_number']=pd.to_numeric(d.original_author,errors='coerce')
    iq=d[(d.sample_id=='A007U')&(d.chain=='B')&(d.original_number.between(400,410))].copy()
    gem=d[(d.sample_id=='A004F')&d.chain.isin(['A','B'])&d.original_number.between(94,113)].copy()
    iq.to_csv(OUT/'iqcg_IQ_motif_400_410_all_residues.csv',index=False)
    d[(d.sample_id=='A007U')&(d.chain=='B')].to_csv(OUT/'iqcg_full_37_residues_original_mapping.csv',index=False)
    reg=[summarize(iq,'A007U','B','literature_IQ_motif_400_410')]
    for chain,q in gem.groupby('chain'):
        reg.append(summarize(q,'A004F',chain,'geminin_N_terminal_binding_segment_94_113'))
    pd.DataFrame(reg).to_csv(OUT/'functional_region_metrics.csv',index=False)
    # Each 100 ns time segment belongs to one long trajectory, not an independent replicate.
    winrows=[]
    for sid,select in [('A007U',iq),('A004F',gem)]:
        w=pd.read_csv(OUT/f'verified_case_trajectories/{sid}/all_sidechains_per_window.csv')
        r=pd.read_csv(OUT/f'verified_case_trajectories/{sid}/all_residues_rmsf_per_window.csv')
        h=select.merge(w[['residue_index','block','hchi','fine_entropy_mean_per_chi','width_mean_per_chi']],on='residue_index')
        h=h.merge(r[['residue_index','block','rmsf']],on=['residue_index','block'],validate='one_to_one')
        h.to_csv(OUT/f'{sid}_functional_region_all_windows.csv',index=False)
        for (chain,block),q in h.groupby(['chain','block']):
            for metric,yfield in [('rmsf','rmsf'),('hchi','hchi')]:
                good=q['valid_'+metric]
                if metric=='hchi':good &= q.sidechain_complete & (q.amino_acid!='P')
                z=q[good]
                for method in METHODS:
                    winrows.append(dict(sample_id=sid,chain=chain,block=block,metric=metric,method=method,n=len(z),mae=float(abs(z[yfield]-z['pred_'+method+'_'+metric]).mean())))
    pd.DataFrame(winrows).to_csv(OUT/'functional_region_per_window_metrics.csv',index=False)
    print(pd.DataFrame(reg).to_string(index=False))
    print(iq[['original_number','amino_acid','md_rmsf','pred_full_rmsf','md_hchi','pred_full_hchi','pred_intra_hchi','pred_aa_sasa_contacts_hchi','pred_static47_sasa_hchi']].to_string(index=False))
    print(pd.DataFrame(winrows).query("sample_id=='A007U' and metric=='hchi'").pivot(index='block',columns='method',values='mae'))

if __name__=='__main__':main()
