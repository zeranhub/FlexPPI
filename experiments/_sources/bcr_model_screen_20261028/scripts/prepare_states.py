from pathlib import Path
from collections import OrderedDict
import hashlib,json
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028')
D=R/'structures';S=R/'states';S.mkdir(exist_ok=True)
AA3='ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split();AA1='ARNDCQEGHILKMFPSTWYV';TO1=dict(zip(AA3,AA1));TO3=dict(zip(AA1,AA3))
states=[];contrasts=[]
def load(pdb):
 data=(D/(pdb.lower()+'_raw.pdb')).read_text().splitlines();selected=OrderedDict()
 for line in data:
  if line.startswith('ENDMDL'):break
  if not line.startswith('ATOM  ') or line[17:20] not in TO1:continue
  key=(line[21],line[22:27],line[12:16]);alt=line[16]
  priority=(1 if alt==' ' else 0,float(line[54:60] or 0),1 if alt=='A' else 0)
  if key not in selected or priority>selected[key][0]:selected[key]=(priority,line[:16]+' '+line[17:])
 return [v[1] for v in selected.values()]
raw={p:load(p) for p in ['1c08','5ayu','6p4a','6p4b','6p4c','6p4d']}
def make(sid,pdb,mapping,kind,changes=None,variable=False,note=''):
 changes=changes or {};out=[];audit=[];seen={}
 for line in raw[pdb]:
  if line[21] not in mapping:continue
  role=mapping[line[21]];num=int(line[22:26]);author=line[22:27].strip()
  if variable and role in ['H','L'] and num>({'H':113,'L':107}[role]):continue
  key=(role,author)
  if key in changes:
   expected,replacement=changes[key]
   assert TO1[line[17:20]]==expected,(sid,key,line[17:20],expected)
   if line[12:16].strip() not in ['N','CA','C','O','OXT']:continue
   seen[key]=(expected,replacement)
   line=line[:17]+TO3[replacement]+line[20:]
  line=line[:21]+role+line[22:]
  out.append(line)
 assert set(seen)==set(changes),(sid,set(changes)-set(seen))
 p=S/(sid+'.pdb');p.write_text('\n'.join(out)+'\nEND\n',encoding='ascii')
 chains={}
 for line in out:chains.setdefault(line[21],OrderedDict())[line[22:27].strip()]=TO1[line[17:20]]
 spec={'id':sid,'pdb_id':pdb.upper(),'path':str(p),'kind':kind,'chains':list(chains),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'note':note,'sequence_only_changes':[{'chain':k[0],'author':k[1],'before':v[0],'after':v[1]} for k,v in seen.items()],'chain_records':{c:{'length':len(v),'sequence':''.join(v.values()),'author_residue_ids':list(v)} for c,v in chains.items()}}
 states.append(spec);return sid
mutations={'HR71A':{('H','71'):('R','A')},'LR45A':{('L','45'):('R','A')},'HQ3A_control':{('H','3'):('Q','A')},'LH34A_control':{('L','34'):('H','A')}}
# Mutants are sequence-only sensitivity inputs: N/CA/C/O are not moved or relaxed.
for context,pdb,mapping in [('bound','1c08',{'A':'L','B':'H','C':'A'}),('apo','5ayu',{'L':'L','H':'H'}),('stripped','1c08',{'A':'L','B':'H'})]:
 wt=make('hh10_'+context+'_WT',pdb,mapping,'observed crystal' if context!='stripped' else 'partner stripped from bound crystal',note='Fv construct; stripped structure is not an independent apo experiment.' if context=='stripped' else 'Deposited Fv coordinates; no missing residues fabricated.')
 for name,change in mutations.items():
  sid=make('hh10_'+context+'_'+name,pdb,mapping,'fixed-backbone sequence counterfactual',changes=change,note='Only residue identity changes. Not an experimentally determined mutant backbone or equilibrated MD structure; side-chain atoms at changed site removed because model reads backbone plus identity.')
  contrasts.append({'id':context+'_'+name,'case':'HyHEL10_2017','before':wt,'after':sid,'type':'fixed-backbone sequence sensitivity','same_context':True,'reference':'>10-fold affinity increase; altered fluctuations' if not name.endswith('control') else 'little changed experimental affinity; not an asserted zero-dynamic-change control','counterfactual':True})
contrasts.append({'id':'hh10_observed_apo_to_bound','case':'HyHEL10_2017','before':'hh10_apo_WT','after':'hh10_bound_WT','type':'independent crystals plus partner context','same_context':False,'counterfactual':False})
contrasts.append({'id':'hh10_same_backbone_partner_context','case':'HyHEL10_2017','before':'hh10_stripped_WT','after':'hh10_bound_WT','type':'fixed-backbone antigen addition','same_context':False,'counterfactual':True})
wt4={('H','4'):('L','F'),('H','33'):('Y','H'),('H','56'):('S','N'),('H','58'):('Y','F')};rev4={k:(b,a) for k,(a,b) in wt4.items()}
for suffix,var in [('full',False),('fv',True)]:
 make('flex_rigid_WT_'+suffix,'6p4a',{'L':'L','H':'H','C':'A'},'observed complex crystal',variable=var)
 make('flex_flexible_4X_'+suffix,'6p4b',{'L':'L','H':'H','C':'A'},'observed complex crystal',variable=var)
 make('flex_apo_4X_'+suffix,'6p4c',{'L':'L','H':'H'},'observed apo Fab crystal',variable=var)
 contrasts.extend([{'id':'flex_joint_state_change_'+suffix,'case':'FlexHEL_2020','before':'flex_rigid_WT_'+suffix,'after':'flex_flexible_4X_'+suffix,'type':'independent crystals; antibody and antigen both changed','same_context':False,'counterfactual':False},{'id':'flex_4X_apo_to_bound_'+suffix,'case':'FlexHEL_2020','before':'flex_apo_4X_'+suffix,'after':'flex_flexible_4X_'+suffix,'type':'same antibody sequence apo/bound crystals plus antigen context','same_context':False,'counterfactual':False}])
make('flex_flexible_4X_copy2_full','6p4b',{'A':'L','B':'H','D':'A'},'observed second crystal copy',note='Crystal copy sensitivity, not independent biological repeat.')
make('flex_rigid_antigen_stripped','6p4a',{'C':'A'},'partner stripped antigen')
make('flex_flexible_antigen_stripped','6p4b',{'C':'A'},'partner stripped antigen')
make('flex_free_self_antigen_6p4d','6p4d',{'A':'A'},'observed free HEL3x crystal',note='Author101 is Arg; not sequence-identical to the rigid/flex D101 foreign antigens.')
make('flex_rigid_backbone_4X_cf','6p4a',{'L':'L','H':'H','C':'A'},'fixed-backbone sequence counterfactual',changes=wt4,variable=True)
make('flex_flexible_backbone_WT_cf','6p4b',{'L':'L','H':'H','C':'A'},'fixed-backbone sequence counterfactual',changes=rev4,variable=True)
contrasts.extend([{'id':'flex_antibody_only_forward_cf','case':'FlexHEL_2020','before':'flex_rigid_WT_fv','after':'flex_rigid_backbone_4X_cf','type':'fixed-backbone four-antibody-substitution sensitivity','same_context':True,'counterfactual':True},{'id':'flex_antibody_only_reverse_cf','case':'FlexHEL_2020','before':'flex_flexible_4X_fv','after':'flex_flexible_backbone_WT_cf','type':'fixed-backbone four-antibody-reversion sensitivity','same_context':True,'counterfactual':True},{'id':'flex_crystal_copy_sensitivity','case':'FlexHEL_2020','before':'flex_flexible_4X_full','after':'flex_flexible_4X_copy2_full','type':'two crystal copies of same complex, not biological repeats','same_context':True,'counterfactual':False}])
for sample in ['A00EJ','A00HE','A00KW']:
 p=Path(r'D:\data\dynarepo_non_nucleic_full')/sample/'structure.pdb'
 assert p.is_file(),p
 states.append({'id':'mdref_'+sample,'pdb_id':{'A00EJ':'6P4B','A00HE':'1VFB','A00KW':'1MLC'}[sample],'path':str(p),'kind':'original saved MD input for exact-label evaluation','sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'sample_id':sample,'note':'Same structure as frozen label source. One local trajectory; time windows are not independent simulation repeats.'})
(R/'states.json').write_text(json.dumps(states,ensure_ascii=False,indent=2),encoding='utf-8')
(R/'contrasts.json').write_text(json.dumps(contrasts,ensure_ascii=False,indent=2),encoding='utf-8')
(R/'screen_design_frozen_before_inference.json').write_text(json.dumps({'states':len(states),'contrasts':len(contrasts),'selection':'Candidate identities and states fixed from literature before predictions. Full residue/region output retained, no successful site filtering.','primary_checkpoint':'Current all-chain-family fold0 full model, pending exact checkpoint identity verification.','method':'Independent absolute-state predictions; deltas computed only afterwards. No retraining or new MD.','Hchi_comparison':'Exclude changed residue identity from Hchi deltas; A/G have no sidechain chi and proline eligibility reported. All positions retained in raw CSV.','evidence_limits':['Counterfactual amino-acid replacement on a fixed backbone is only model sensitivity, not experimental mutant coordinates.','NMR millisecond exchange is qualitative external context, not the same numerical observable as RMSF/Hchi.','Two crystal copies are not simulation replicates.','Existing MD labels provide absolute-state checks only for their exact recorded input context.']},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'states':len(states),'contrasts':len(contrasts),'output':str(R),'counterfactual_states':sum(s['kind']=='fixed-backbone sequence counterfactual' for s in states)},ensure_ascii=False))
