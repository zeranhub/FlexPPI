from pathlib import Path
import csv,json,sys
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
sys.path.insert(0,r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\numbering\deps')
from anarcii import Anarcii
audit=json.loads((R/'audits/current_inference.json').read_text())
seq={ch:audit['chains'][ch] for ch in ['B','C']}
n=Anarcii(seq_type='antibody',mode='accuracy',cpu=True,ncpu=2).number(seq)
(R/'audits/6a77_anarcii_raw.json').write_text(json.dumps(n,indent=2),encoding='utf-8')
with (R/'6A77_all_positions_predictions_and_MD.csv').open(encoding='utf-8-sig') as f: rows=[r for r in csv.DictReader(f) if r['variant']=='full']
out=[]
for ch,s in seq.items():
 num=n[ch];assert num['error'] is None
 assert num['chain_type']==('K' if ch=='B' else 'H'),(ch,num['chain_type'])
 rr=[r for r in rows if r['chain']==ch];cursor=num['query_start'];positions={}
 for (p,ins),aa in num['numbering']:
  if aa=='-':continue
  assert s[cursor]==aa;positions[cursor]=(p,ins);cursor+=1
 for i,r in enumerate(rr):
  p,ins=positions.get(i,(None,''));region='framework' if p is not None else 'constant_or_unnumbered'
  if p is not None:
   if 27<=p<=38:region='CDR1'
   elif 56<=p<=65:region='CDR2'
   elif 105<=p<=117:region='CDR3'
  out.append({'chain':ch,'author':r['author'],'aa':r['aa'],'role':'light' if ch=='B' else 'heavy','imgt_number':p,'imgt_insertion':ins,'region':region})
with (R/'6A77_antibody_numbering.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
print('COMPLETE',len(out))
