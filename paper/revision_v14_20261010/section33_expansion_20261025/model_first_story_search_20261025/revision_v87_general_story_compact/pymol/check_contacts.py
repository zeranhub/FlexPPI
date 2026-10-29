from pathlib import Path
import numpy as np
p=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v85_1MXE\pymol\inputs\1MXE_CaM_native4_147_CaMKI294_318.pdb')
a=[]
for s in p.read_text().splitlines():
 if s.startswith('ATOM'):a.append((s[21],int(s[22:26]),s[17:20],s[12:16].strip(),np.array([float(s[k:k+8]) for k in [30,38,46]])))
for site in [303,316]:
 t=np.array([r[4] for r in a if r[0]=='B' and r[1]==site and r[3] not in ['N','CA','C','O']]);d={}
 for r in a:
  if r[0]=='A':
   dd=np.linalg.norm(t-r[4],axis=1).min();k=(r[1],r[2]);d[k]=min(d.get(k,100),dd)
 print(site,[(k,round(v,2)) for k,v in sorted(d.items()) if v<=4.5])