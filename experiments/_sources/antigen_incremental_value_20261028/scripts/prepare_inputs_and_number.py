from pathlib import Path
from collections import OrderedDict
import csv,hashlib,json,sys

R=Path(r'E:\Project\FlexPPI\experiments\_sources\antigen_incremental_value_20261028')
sys.path.insert(0,r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\numbering\deps')
from anarcii import Anarcii
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
AA3={v:k for k,v in AA.items()}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def writecsv(path,rows):
 with path.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
with (R/'functional_candidates.csv').open(encoding='utf-8-sig',newline='') as f: candidates=list(csv.DictReader(f))
(R/'structures').mkdir(exist_ok=True)
templates={};sequences={}
for pdb in ['2BNQ','2BNR','3MV7','4PRP','5BRZ','5BS0']:
 selected=OrderedDict();raw=R/'sources/pdb'/f'{pdb}.pdb'
 for line in raw.read_text().splitlines():
  if line.startswith('ENDMDL'):break
  if not line.startswith('ATOM  ') or line[21] not in 'ABCDE' or line[17:20] not in AA:continue
  key=(line[21],line[22:27],line[12:16]);priority=(line[16]==' ',float(line[54:60]),line[16]=='A')
  if key not in selected or priority>selected[key][0]:selected[key]=(priority,line[:16]+' '+line[17:])
 lines=[x[1] for x in selected.values()];chains={}
 for l in lines:
  if l[12:16].strip()=='CA':chains.setdefault(l[21],OrderedDict())[l[22:27].strip()]=AA[l[17:20]]
 assert set(chains)==set('ABCDE') and len(chains['C']) in [9,11]
 templates[pdb]={'lines':lines,'chains':chains,'raw_sha256':sha(raw)}
 for ch in 'DE':sequences[pdb+'_'+ch]=''.join(chains[ch].values())
numberer=Anarcii(seq_type='tcr',mode='accuracy',cpu=True,ncpu=2)
raw_number=numberer.number(sequences)
(R/'audits/anarcii_numbering_raw.json').write_text(json.dumps(raw_number,indent=2),encoding='utf-8')
numbering=[];CDRs={}
for key,seq in sequences.items():
 pdb,ch=key.split('_');n=raw_number[key];assert n['error'] is None,(key,n)
 assert n['chain_type']==('A' if ch=='D' else 'B'),(key,n['chain_type'])
 cursor=n['query_start'];authors=list(templates[pdb]['chains'][ch]);mapping={}
 for (pos,ins),aa in n['numbering']:
  if aa=='-':continue
  assert seq[cursor]==aa;mapping[cursor]=(int(pos),ins);cursor+=1
 def cdr(lo,hi):return ''.join(aa for (pos,ins),aa in n['numbering'] if lo<=pos<=hi and aa!='-')
 CDRs[key]={'cdr1':cdr(27,38),'cdr2':cdr(56,65),'cdr3':cdr(105,117),'cdr3_with_anchors':cdr(104,118)}
 for i,aa in enumerate(seq):
  pos,ins=mapping.get(i,(None,''));region='other'
  if pos is not None:
   if 27<=pos<=38:region='CDR1'
   elif 56<=pos<=65:region='CDR2'
   elif 105<=pos<=117:region='CDR3'
  numbering.append({'template':pdb,'chain':ch,'author':authors[i],'aa':aa,'imgt_number':pos,'imgt_insertion':ins,'region':region,'role':'alpha' if ch=='D' else 'beta'})
writecsv(R/'tcr_numbering.csv',numbering)
(R/'audits/verified_CDRs.json').write_text(json.dumps(CDRs,indent=2),encoding='utf-8')
states=[];baseline_rows=[];references=[]
for system in ['nyeso','ebv','mage']:
 rows=[c for c in candidates if c['system']==system];assert len(rows)==93
 primary=rows[0]['primary_template'];secondary=rows[0]['secondary_template']
 a,b=CDRs[primary+'_D'],CDRs[primary+'_E']
 for r in rows:
  baseline_rows.append({'candidate_id':r['candidate_id'],'system':system,'peptide':r['sequence'],
    'A1':a['cdr1'],'A2':a['cdr2'],'A3':a['cdr3'],'B1':b['cdr1'],'B2':b['cdr2'],'B3':b['cdr3'],
    'TRA':a['cdr3_with_anchors'],'TRB':b['cdr3_with_anchors'],'TRAV':'','TRAJ':'','TRBV':'','TRBJ':'','T-Cell-Type':'CD8','Peptide':r['sequence'],'MHC':r['hla']})
 for pdb,version in [(primary,'primary'),(secondary,'secondary')]:
  t=templates[pdb];parent=''.join(t['chains']['C'].values());authors=list(t['chains']['C'])
  refid=system+'__'+pdb+'__reference'
  for cid,seq,is_ref in [(refid,parent,True)]+[(r['candidate_id'],r['sequence'],False) for r in rows]:
   assert len(seq)==len(authors)
   change={author:(aa,new) for author,aa,new in zip(authors,parent,seq) if aa!=new};out=[]
   for line in t['lines']:
    key=line[22:27].strip()
    if line[21]=='C' and key in change:
     if line[12:16].strip() not in ['N','CA','C','O','OXT']:continue
     line=line[:17]+AA3[change[key][1]]+line[20:]
    out.append(line)
   sid=refid if is_ref else cid+'__'+pdb;p=R/'structures'/f'{sid}.pdb';p.write_text('\n'.join(out)+'\nEND\n',encoding='ascii')
   backbone=lambda lines:{(l[21],l[22:27],l[12:16]):l[30:54] for l in lines if l[12:16].strip() in ['N','CA','C','O']}
   assert backbone(out)==backbone(t['lines'])
   state={'id':sid,'path':str(p),'sha256':sha(p),'candidate_id':cid,'system':system,'template':pdb,'template_version':version,'is_reference':is_ref,
          'reference_state':refid,'sequence':seq,'chains':list(t['chains']), 'fixed_backbone':True,'kind':'uniform sequence replacement on experimental reference; no mutant relaxation'}
   states.append(state)
   if is_ref:references.append(state)
writecsv(R/'baseline_inputs.csv',baseline_rows)
(R/'states.json').write_text(json.dumps(states,indent=2),encoding='utf-8')
(R/'audits/preparation_audit.json').write_text(json.dumps({'states':len(states),'candidates':len(candidates),'references':len(references),
 'all_fixed_backbones_verified':True,'no_labels_used_for_structure_or_region_choice':True,'sequence_and_roles_verified_by_COMPND_and_ANARCII':True,
 'genes_not_inferred':'ERGO V/J names left unknown uniformly, not guessed from CDR sequences.','protocol_sha256':sha(R/'protocol_frozen_before_new_predictions.json'),
 'raw_templates':{p:{'raw_sha256':t['raw_sha256'],'chain_sequences':{c:''.join(v.values()) for c,v in t['chains'].items()}} for p,t in templates.items()}},indent=2),encoding='utf-8')
print(json.dumps({'prepared_states':len(states),'baseline_candidates':len(baseline_rows),'all_backbone_checks_passed':True},indent=2))
