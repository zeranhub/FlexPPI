from pathlib import Path
import csv,json,sys
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
sys.path.insert(0,r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\numbering\deps')
from anarcii import Anarcii
states=json.loads((R/'pair_screen_states.json').read_text())
unique={}
for s in states:
 for c,r in s['chain_records'].items():
  if c in ['H','L']:unique.setdefault(r['sequence'],s['id']+'_'+c)
n=Anarcii(seq_type='antibody',mode='accuracy',cpu=True,ncpu=2).number({name:seq for seq,name in unique.items()})
(R/'audits/paired_anarcii_raw.json').write_text(json.dumps(n,indent=2),encoding='utf-8')
byseq={s:n[name] for s,name in unique.items()};rows=[]
for s in states:
 for c,r in s['chain_records'].items():
  if c not in ['H','L']:continue
  seq=r['sequence'];v=byseq[seq];assert v['error'] is None
  assert v['chain_type']==('H' if c=='H' else 'K')
  cursor=v['query_start'];mapping={}
  for (p,ins),aa in v['numbering']:
   if aa=='-':continue
   assert seq[cursor]==aa;mapping[cursor]=(p,ins);cursor+=1
  for i,(author,aa) in enumerate(zip(r['authors'],seq)):
   pos,ins=mapping.get(i,(None,''));region='constant_or_unnumbered' if pos is None else 'framework'
   if pos is not None:
    if 27<=pos<=38:region='CDR1'
    elif 56<=pos<=65:region='CDR2'
    elif 105<=pos<=117:region='CDR3'
   rows.append({'state_id':s['id'],'chain':c,'author':author,'aa':aa,'imgt_number':pos,'imgt_insertion':ins,'region':region})
with (R/'paired_numbering.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('COMPLETE',len(unique),len(rows))
