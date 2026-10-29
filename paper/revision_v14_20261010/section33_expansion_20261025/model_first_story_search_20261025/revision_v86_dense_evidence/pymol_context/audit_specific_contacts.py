from pathlib import Path
import json,numpy as np
out=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v86_dense_evidence\pymol_context')
a=[]
for s in (out/'inputs/1A06_observed_native.pdb').read_text().splitlines():
 if s.startswith('ATOM'):a.append((s[21],int(s[22:26]),s[17:20],s[12:16].strip(),np.array([float(s[k:k+8]) for k in [30,38,46]])))
f=[q for q in a if q[1]==316 and q[3] not in ['N','CA','C','O']]
rows=[]
for site in [45,99]:
 t=[q for q in a if q[1]==site];d=np.array([[np.linalg.norm(q[4]-z[4]) for z in t] for q in f]);i,j=np.unravel_index(d.argmin(),d.shape)
 rows.append({'site316_sidechain_atom':f[i][3],'contact_residue':site,'contact_residue_name':t[j][2],'closest_atom':t[j][3],'min_distance_A':float(d[i,j])})
(out/'M316_specific_contacts.json').write_text(json.dumps(rows,indent=2),encoding='utf-8');print(json.dumps(rows,indent=2))