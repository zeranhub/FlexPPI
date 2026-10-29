from pathlib import Path
from collections import OrderedDict
import hashlib,json,urllib.request
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
D=R/'structures';D.mkdir(exist_ok=True)
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'));AA3={v:k for k,v in AA.items()}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(pdb):
 p=D/f'{pdb}_raw.pdb'
 if not p.exists():
  with urllib.request.urlopen(f'https://files.rcsb.org/download/{pdb}.pdb',timeout=35) as f:p.write_bytes(f.read())
 selected=OrderedDict()
 for l in p.read_text().splitlines():
  if l.startswith('ENDMDL'):break
  if not l.startswith('ATOM  ') or l[17:20] not in AA:continue
  key=(l[21],l[22:27],l[12:16]);priority=(l[16]==' ',float(l[54:60]),l[16]=='A')
  if key not in selected or priority>selected[key][0]:selected[key]=(priority,l[:16]+' '+l[17:])
 return [x[1] for x in selected.values()]
raw={p:load(p) for p in ['6A76','6A77','6A78','6A79']}
states=[];cons=[]
def make(sid,pdb,mapping,kind,mutation=None,keep=None):
 out=[];changed=False
 for l in raw[pdb]:
  if l[21] not in mapping:continue
  role=mapping[l[21]];author=l[22:27].strip()
  if keep is not None and (role,author) not in keep:continue
  if mutation and (role,author)==('H','103'):
   old,new=mutation;assert AA[l[17:20]]==old
   if l[12:16].strip() not in ['N','CA','C','O','OXT']:continue
   l=l[:17]+AA3[new]+l[20:];changed=True
  out.append(l[:21]+role+l[22:])
 assert mutation is None or changed
 p=D/f'{sid}.pdb';p.write_text('\n'.join(out)+'\nEND\n',encoding='ascii')
 chains={}
 for l in out:
  if l[12:16].strip()=='CA':chains.setdefault(l[21],OrderedDict())[l[22:27].strip()]=AA[l[17:20]]
 state={'id':sid,'pdb':pdb,'path':str(p),'sha256':sha(p),'kind':kind,'chains':list(chains),'sequence_changes':[] if not mutation else [{'chain':'H','author':'103','from':mutation[0],'to':mutation[1]}],
        'chain_records':{c:{'sequence':''.join(v.values()),'authors':list(v)} for c,v in chains.items()}}
 states.append(state);return state
for pdb,label in [('6A78','WT'),('6A79','P103A')]:
 for copy,mapping in [(1,{'A':'A','L':'L','H':'H'}),(2,{'B':'A','M':'L','I':'H'})]:
  make(f'scfv_{label}_copy{copy}',pdb,mapping,'observed scFv-Ig5 biological assembly')
for copy,mapping in [(1,{'A':'A','L':'L','H':'H'}),(2,{'B':'A','M':'L','I':'H'})]:
 cons.append({'id':f'P103A_minus_WT_copy{copy}','before':f'scfv_WT_copy{copy}','after':f'scfv_P103A_copy{copy}','type':'matched construct, independent crystals','counterfactual':False})
 for pdb,label,mut in [('6A78','WT_P103A_cf',('P','A')),('6A79','P103A_A103P_cf',('A','P'))]:
  s=make(f'scfv_{label}_copy{copy}',pdb,mapping,'fixed-backbone sequence counterfactual',mutation=mut)
  before=f'scfv_WT_copy{copy}' if pdb=='6A78' else f'scfv_P103A_copy{copy}'
  cons.append({'id':s['id']+'_minus_parent','before':before,'after':s['id'],'type':'fixed-backbone identity sensitivity','counterfactual':True})
# Explicit shared-coverage sensitivity; no coordinate reconstruction.
observed=states[:4]
shared=set.intersection(*[{(c,a) for c,v in s['chain_records'].items() for a in v['authors']} for s in observed])
for s in observed:
 mapping={'A':'A','L':'L','H':'H'} if s['id'].endswith('copy1') else {'B':'A','M':'L','I':'H'}
 make(s['id']+'_shared',s['pdb'],mapping,'observed coordinates cropped to all-four common residue coverage',keep=shared)
for copy in [1,2]:cons.append({'id':f'P103A_minus_WT_copy{copy}_shared','before':f'scfv_WT_copy{copy}_shared','after':f'scfv_P103A_copy{copy}_shared','type':'matched construct common-coverage sensitivity','counterfactual':False})
for label in ['WT','P103A']:
 cons.append({'id':f'{label}_copy2_minus_copy1','before':f'scfv_{label}_copy1','after':f'scfv_{label}_copy2','type':'within-crystal copy sensitivity','counterfactual':False})
# Fab comparison, chain roles verified in source COMPND.
for pdb in ['6A76','6A77']:
 mapping={'L':'L','H':'H'}
 if pdb=='6A77':mapping['A']='A'
 make('fab_'+('apo' if pdb=='6A76' else 'bound'),pdb,mapping,'observed WT Fab crystal')
make('fab_bound_antigen_stripped','6A77',{'L':'L','H':'H'},'antigen removed from bound coordinates; not observed apo')
cons.extend([{'id':'Fab_bound_minus_apo','before':'fab_apo','after':'fab_bound','type':'same WT Fab separate apo/holo crystals','counterfactual':False},
 {'id':'Fab_same_backbone_antigen_addition','before':'fab_bound_antigen_stripped','after':'fab_bound','type':'synthetic context perturbation','counterfactual':True}])
contract={'status':'frozen_before_new_pair_predictions','primary':'Observed first biological assembly: 6A78 WT vs6A79 P103A. Second copies and common-coverage inputs are sensitivity, not independent repetitions.',
 'scope':'All common positions, both variable domains, all six IMGT CDRs, fixed network H31D/H97E/Ag67R and complete sequence neighbors +/-2. Ag67 is the only antigen Arg contacting the two acidic sites in these structures; literatureR477 vs current UniProt numbering is unresolved. Ag31 is Val and must not be used as Arg. No pre-imposed RMSF/Hchi direction from salt-bridge stability.',
 'Hchi':'Exclude any changed AA identity and chemically undefined A/G. Do not interpret Hchi as water entropy, salt-bridge occupancy, binding entropy or KD.',
 'reference':'Only WT Fab6A77 has matching local frozen MD labels. This reference does not quantitatively validate WT/mutant scFv differences. One local XTC only.',
 'coverage':'Use author residue ID joins and retain actual sequence checks. Extra termini and missing atoms recorded, not row subtraction. scFv unresolved linker segments not built.',
 'model':'Current2430 all-chain-family fold0 full and intrachain epoch40; strict SHA verified.',
 'no_retraining':True,'no_new_MD':True,'do_not_modify_manuscript':True}
(R/'pair_screen_states.json').write_text(json.dumps(states,indent=2),encoding='utf-8')
(R/'pair_screen_contrasts.json').write_text(json.dumps(cons,indent=2),encoding='utf-8')
(R/'pair_screen_contract.json').write_text(json.dumps(contract,indent=2),encoding='utf-8')
(R/'audits/pair_preparation.json').write_text(json.dumps({'state_count':len(states),'contrasts':len(cons),'shared_all_four_residues':len(shared),'raw_sources_sha256':{p:sha(D/f'{p}_raw.pdb') for p in raw},'frozen_states_sha256':sha(R/'pair_screen_states.json'),'contract_sha256':sha(R/'pair_screen_contract.json')},indent=2),encoding='utf-8')
print('PREPARED',len(states),len(cons),len(shared),flush=True)
