from pathlib import Path
import json,sys,csv,hashlib
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028')
sys.path.insert(0,r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\numbering\deps')
from anarcii import Anarcii
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
states=json.loads((R/'states.json').read_text(encoding='utf-8'))
chains={};unique={};records={}
for s in states:
 cr={}
 for l in Path(s['path']).read_text().splitlines():
  if l.startswith('ENDMDL'):break
  if l.startswith('ATOM  ') and l[12:16].strip()=='CA' and l[17:20] in AA:cr.setdefault(l[21],{})[l[22:27].strip()]=AA[l[17:20]]
 for ch,rs in cr.items():
  seq=''.join(rs.values());key=s['id']+'_'+ch;records[key]=(s['id'],ch,list(rs),seq)
  # Canonical role or existing MD antibody chains A/B; antigen chains excluded.
  if ch in ['H','L'] or (s.get('sample_id')=='A00KW' and ch in ['A','B']):
   chains[key]=seq;unique.setdefault(seq,key)
numberer=Anarcii(seq_type='antibody',mode='accuracy',cpu=True,ncpu=2)
raw=numberer.number({name:seq for seq,name in unique.items()})
(R/'reference/anarcii_antibody_numbering_raw.json').write_text(json.dumps(raw,indent=2),encoding='utf-8')
byseq={seq:raw[name] for seq,name in unique.items()}
rows=[]
for key,seq in chains.items():
 sid,ch,authors,_=records[key];v=byseq[seq]
 assert v['error'] is None,(key,v)
 cursor=v['query_start'];num_map={}
 for (position,insertion),aa in v['numbering']:
  if aa=='-':continue
  assert seq[cursor]==aa,(key,cursor,aa,seq[cursor])
  num_map[authors[cursor]]=(int(position),insertion,aa);cursor+=1
 for i,author in enumerate(authors):
  pos,ins,aa=num_map.get(author,(None,'',seq[i]));region='framework'
  if pos is None:region='constant_or_unnumbered'
  elif 27<=pos<=38:region='CDR1'
  elif 56<=pos<=65:region='CDR2'
  elif 105<=pos<=117:region='CDR3'
  role='H' if v['chain_type']=='H' else 'L'
  rows.append(dict(state_id=sid,chain=ch,author_residue_id=author,aa=seq[i],antibody_role=role,imgt_number=pos,imgt_insertion=ins,imgt_region=region,anarcii_score=v['score']))
with (R/'reference/antibody_residue_numbering.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps({'unique_sequences_numbered':len(unique),'state_chains':len(chains),'mapped_rows':len(rows),'all_numbered_residue_identities_match':True}))
