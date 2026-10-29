from pathlib import Path
import argparse,csv,gzip,hashlib,json,sys,time
import numpy as np
import torch
R=Path(r'E:\Project\FlexPPI\experiments\_sources\antigen_incremental_value_20261028')
SRC=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\src')
sys.path.insert(0,str(SRC))
from flexppi.inference import _structure_example
from flexppi.data.dataset import collate_states
from flexppi.data.structures import parse_pdb
from flexppi.config import ModelConfig
from flexppi.model.information_ablation import build_information_ablation
p=argparse.ArgumentParser();p.add_argument('--variant',choices=['full','intrachain'],required=True);p.add_argument('--threads',type=int,default=2);p.add_argument('--resume',action='store_true');args=p.parse_args()
v=args.variant;torch.set_num_threads(args.threads);torch.manual_seed(7)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
states=json.loads((R/'states.json').read_text());assert len(states)==564
cp_path=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028\runtime\family_folds\runs')/v/'fold_0/epoch_040.pt'
expected={'full':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5','intrachain':'6b069d9d3521e95f5f2d63a3f19d64fb4071862efa4682473387b0ca45ff73c3'}
assert sha(cp_path)==expected[v]
cp=torch.load(cp_path,map_location='cpu',weights_only=False)
config=ModelConfig(**cp['config']['model']);assert list(config.metrics)==['rmsf','chi_entropy']
assert cp['epoch']==40 and not config.use_condition_features and not config.enable_response_head
model=build_information_ablation(config,v).eval();model.load_state_dict(cp['model'],strict=True)
number=list(csv.DictReader((R/'tcr_numbering.csv').open(encoding='utf-8-sig')))
P=R/'predictions'/v;P.mkdir(parents=True,exist_ok=True)
regions=[];entries=[];start=time.monotonic();total=0;cached_count=0
if args.resume:states=sorted(states,key=lambda s:(s['template_version']!='primary',s['system'],not s['is_reference'],s['id']))
with gzip.open(P/'all_residue_predictions.csv.gz','wt',encoding='utf-8',newline='') as out:
 fields=['state_id','candidate_id','system','template','model_index','chain','author','aa','model_valid','geometry_valid','hchi_valid','rmsf','hchi']
 writer=csv.DictWriter(out,fieldnames=fields);writer.writeheader()
 for j,s in enumerate(states):
  assert sha(s['path'])==s['sha256'];st=parse_pdb(s['path']);assert st.num_residues<=1024
  ex=_structure_example(s['path'],st.unique_chains(),300.,100.,0.)
  cache_path=P/(s['id']+'.npz');cached=False
  if args.resume and cache_path.is_file():
   try:
    z=np.load(cache_path)
    assert str(z['input_sha256'])==s['sha256'] and str(z['checkpoint_sha256'])==expected[v]
    assert np.array_equal(z['chain'],st.chain_ids) and np.array_equal(z['author'],st.author_residue_ids) and ''.join(z['aa'])==st.sequence
    mean=z['mean'];lv=z['log_variance'];mask=z['model_mask'];geo=z['geometry_mask'];assert mean.shape==(st.num_residues,2) and np.isfinite(mean[mask]).all()
    cached=True;cached_count+=1
   except (AssertionError,ValueError,EOFError,OSError):cached=False
  if not cached:
   with torch.no_grad():pred=model.forward_state(collate_states([ex]))
   mean=pred['mean'][0].cpu().numpy();lv=pred['log_variance'][0].cpu().numpy();mask=pred['mask'][0].cpu().numpy();geo=ex['geometry_mask'].numpy()
  aa=np.array(list(st.sequence));valid_chi=mask&geo&~np.isin(aa,['A','G']);assert np.isfinite(mean[mask]).all()
  if not cached:np.savez_compressed(cache_path,mean=mean,log_variance=lv,model_mask=mask,geometry_mask=geo,hchi_valid=valid_chi,chain=np.array(st.chain_ids),author=np.array(st.author_residue_ids),aa=aa,input_sha256=s['sha256'],checkpoint_sha256=expected[v])
  index={(c,a):i for i,(c,a) in enumerate(zip(st.chain_ids,st.author_residue_ids))}
  for role,ch in [('alpha','D'),('beta','E')]:
   for region,allowed in [('CDR3',['CDR3']),('CDR12',['CDR1','CDR2'])]:
    indices=[index[(n['chain'],n['author'])] for n in number if n['template']==s['template'] and n['chain']==ch and n['region'] in allowed]
    assert indices
    for k,metric in enumerate(['rmsf','hchi']):
     eligible=[i for i in indices if ((mask&geo)[i] if k==0 else valid_chi[i])]
     assert eligible
     regions.append({'state_id':s['id'],'candidate_id':s['candidate_id'],'system':s['system'],'template':s['template'],'template_version':s['template_version'],'is_reference':s['is_reference'],'reference_state':s['reference_state'],'variant':v,'region':role+'_'+region,'metric':metric,'n':len(eligible),'mean':float(mean[eligible,k].mean()),'author_ids':'|'.join(st.author_residue_ids[i] for i in eligible)})
  for i in range(st.num_residues):
   writer.writerow(dict(state_id=s['id'],candidate_id=s['candidate_id'],system=s['system'],template=s['template'],model_index=i,chain=st.chain_ids[i],author=st.author_residue_ids[i],aa=st.sequence[i],model_valid=bool(mask[i]),geometry_valid=bool(geo[i]),hchi_valid=bool(valid_chi[i]),rmsf=float(mean[i,0]),hchi=float(mean[i,1])))
  total+=st.num_residues
  entries.append({'state_id':s['id'],'input_sha256':s['sha256'],'residues':st.num_residues,'chains':st.unique_chains(),'finite':True,'reused_verified_prediction':cached})
  if j%20==0 or j==len(states)-1:print(v,j+1,'/',len(states),'elapsed',round(time.monotonic()-start,1),flush=True)
with (P/'regional_absolute_predictions.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(regions[0]));w.writeheader();w.writerows(regions)
(R/'audits'/f'inference_{v}.json').write_text(json.dumps({'status':'complete','distinct_model_inputs_completed':len(states),'new_forward_passes_this_execution':len(states)-cached_count,'verified_cached_predictions':cached_count,'residue_rows':total,'strict_load':True,'checkpoint':str(cp_path),'checkpoint_sha256':expected[v],'source':str(SRC),'source_files_sha256':{str(p.relative_to(SRC)):sha(p) for p in SRC.rglob('*.py')},'states_sha256':sha(R/'states.json'),'protocol_sha256':sha(R/'protocol_frozen_before_new_predictions.json'),'runtime':{'torch':torch.__version__,'numpy':np.__version__,'python':sys.version,'threads':args.threads},'seconds':time.monotonic()-start,'entries':entries},indent=2),encoding='utf-8')
print('COMPLETE',v,total,flush=True)
