"""A fixed all-four-state same-chemistry Hchi subset, in addition to pairwise sets."""
from pathlib import Path
import pandas as pd,numpy as np
from Bio.Align import PairwiseAligner
OUT=Path(__file__).resolve().parent
d=pd.read_csv(OUT/'MEL5_four_states_all_residues.csv',dtype={'id':str})
f={p:g for p,g in d.groupby('pdb',sort=False)};ref=f['3hg1'];maps={}
for p,g in f.items():
    maps[p]={}
    for role in ['mhc','peptide']:
        a=ref[ref.role.eq(role)].reset_index(drop=True);b=g[g.role.eq(role)].reset_index(drop=True)
        al=PairwiseAligner();al.mode='global';al.match_score=2;al.mismatch_score=-1;al.open_gap_score=-6;al.extend_gap_score=-.5
        hit=al.align(''.join(a.aa),''.join(b.aa))[0]
        for i,j in zip(*hit.indices):
            if i>=0 and j>=0:maps[p][(role,str(a.iloc[i].id))]=b.iloc[j]
rows=[]
for key in maps['3hg1']:
    rr=[maps[p].get(key) for p in ['3hg1','4jff','4jfd','4jfe']]
    interface=any(bool(r.interface) for r in rr if r is not None)
    common=all(r is not None and bool(r.hchi_valid) for r in rr) and len({r.aa for r in rr if r is not None})==1
    row={'reference_role':key[0],'reference_id':key[1],'four_state_interface_union':interface,'same_AA_valid_Hchi_all_four':common}
    for p,r in zip(['3hg1','4jff','4jfd','4jfe'],rr):
        for m in ['chain','id','aa','hchi_valid','hchi_full','hchi_md','hchi_var']:
            row[p+'_'+m]=r[m] if r is not None else None
    rows.append(row)
q=pd.DataFrame(rows)
q[q.four_state_interface_union].to_csv(OUT/'MEL5_four_state_union_all30_Hchi_eligibility.csv',index=False,encoding='utf-8-sig',float_format='%.17g')
v=q[q.four_state_interface_union&q.same_AA_valid_Hchi_all_four]
v.to_csv(OUT/'MEL5_fixed_sameAA_Hchi_interface_all_four_states_22positions.csv',index=False,encoding='utf-8-sig',float_format='%.17g')
print('Fixed four-state common Hchi',len(v),'of union',q.four_state_interface_union.sum())
