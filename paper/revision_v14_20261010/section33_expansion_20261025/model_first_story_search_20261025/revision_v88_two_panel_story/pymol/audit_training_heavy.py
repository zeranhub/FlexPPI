from pathlib import Path
import json, numpy as np
out=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v85_1MXE\pymol')
j=json.loads((out/'coordinate_audit.json').read_text())
def heavy(path):
 d={}
 for s in Path(path).read_text().splitlines():
  if s.startswith(('ATOM  ','HETATM')) and s[76:78].strip() not in ['H','D']:
   d[(s[12:16].strip(),s[17:20],s[21],int(s[22:26]),s[26])]=[float(s[i:i+8]) for i in [30,38,46]]
 return d
a=heavy(j['source']);b=heavy(j['training_payload_structure']);common=set(a)&set(b)
j['source_training_heavy_atom_audit']={'source_heavy_atom_count':len(a),'training_heavy_atom_count':len(b),'common_heavy_atom_keys':len(common),'only_source_keys':list(set(a)-set(b)),'only_training_keys':list(set(b)-set(a)),'maximum_common_coordinate_change_A':max(float(np.max(np.abs(np.array(a[k])-np.array(b[k])))) for k in common),'common_coordinates_exactly_equal':all(a[k]==b[k] for k in common)}
(out/'coordinate_audit.json').write_text(json.dumps(j,indent=2),encoding='utf-8');print(json.dumps(j['source_training_heavy_atom_audit'],indent=2))