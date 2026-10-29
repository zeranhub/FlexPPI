from pathlib import Path
import json,hashlib,shutil,collections
import numpy as np
out=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v85_1MXE\pymol')
out.mkdir(parents=True,exist_ok=True);(out/'inputs').mkdir(exist_ok=True)
src=Path(r'E:\Project\FlexPPI\experiments\_sources\baselines\cohort\structures\record_aba5e1f3db118b0aba8895de.pdb')
train=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\dataset\structures\A001K.pdb')
shutil.copy2(src,out/'inputs'/'family_fold1_input.pdb')
def parse(p):
 rows=[]
 for s in p.read_text().splitlines():
  if s.startswith(('ATOM  ','HETATM')):
   rows.append((s[12:16].strip(),s[17:20],s[21],int(s[22:26]),s[26],*[float(s[i:i+8]) for i in [30,38,46]]))
 return rows
rows=parse(src); rr=parse(train)
a=np.array([r[5:] for r in rows]);b=np.array([r[5:] for r in rr])
relabels=[]
for s in src.read_text().splitlines():
 if s.startswith(('ATOM  ','HETATM','TER   ')) and s[21:22] in ['A','B']:
  chain=s[21];old=int(s[22:26]);new=old+3 if chain=='A' else old+149
  s=s[:22]+f'{new:4d}'+s[26:]
 relabels.append(s)
path=out/'inputs'/'1MXE_CaM_native4_147_CaMKI294_318.pdb';path.write_text('\n'.join(relabels)+'\n')
after=parse(path)
audit={
 'source':str(src),'training_payload_structure':str(train),
 'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
 'training_sha256':hashlib.sha256(train.read_bytes()).hexdigest(),
 'source_training_same_atom_records':bool([r[:5] for r in rows]==[r[:5] for r in rr]),
 'source_training_coordinate_max_change_A':float(np.max(np.abs(a-b))) if a.shape==b.shape else None,
 'source_atom_count':len(rows),'native_atom_count':len(after),
 'native_coordinate_max_change_A':float(np.max(np.abs(a-np.array([r[5:] for r in after])))),
 'renumbering':{'CaM_chain_A':'source residue 1–144 → native 4–147 (+3)','CaMKI_chain_B':'source residue 145–169 → native 294–318 (+149)'},
 'atoms_added':0,'atoms_removed':0,
 'native_chain_residues':{c:sorted(set(r[3] for r in after if r[2]==c)) for c in sorted(set(r[2] for r in after))},
 'native_CA_sequence_3letter':{c:[r[1] for r in after if r[2]==c and r[0]=='CA'] for c in sorted(set(r[2] for r in after))}
}
(out/'coordinate_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8');print(json.dumps(audit,indent=2))