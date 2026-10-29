from pathlib import Path
import json,sys,time,hashlib,csv
import numpy as np
import torch
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028')
SRC=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\src')
sys.path.insert(0,str(SRC))
from flexppi.inference import _structure_example
from flexppi.data.dataset import collate_states
from flexppi.data.structures import parse_pdb
from flexppi.config import ModelConfig
from flexppi.model.information_ablation import build_information_ablation
P=R/'predictions';P.mkdir(exist_ok=True)
states=json.loads((R/'states.json').read_text(encoding='utf-8'))
torch.set_num_threads(2);torch.manual_seed(7)
prepared={}
for s in states:
 assert hashlib.sha256(Path(s['path']).read_bytes()).hexdigest()==s['sha256'],s['id']
 st=parse_pdb(s['path']);chs=st.unique_chains()
 assert 2<=st.num_residues<=1024,(s['id'],st.num_residues)
 if 'chains' in s:assert set(chs)==set(s['chains'])
 ex=_structure_example(s['path'],chs,300.,100.,0.)
 prepared[s['id']]=(st,ex)
allrows=[];audit=[]
for variant,expected in [('full','ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5'),('intrachain','6b069d9d3521e95f5f2d63a3f19d64fb4071862efa4682473387b0ca45ff73c3')]:
 checkpoint=R/f'runtime/family_folds/runs/{variant}/fold_0/epoch_040.pt'
 digest=hashlib.sha256(checkpoint.read_bytes()).hexdigest();assert digest==expected
 cp=torch.load(checkpoint,map_location='cpu',weights_only=False)
 conf=ModelConfig(**cp['config']['model']);assert list(conf.metrics)==['rmsf','chi_entropy']
 assert cp['epoch']==40 and not conf.use_condition_features and not conf.enable_response_head
 model=build_information_ablation(conf,variant);model.load_state_dict(cp['model'],strict=True);model.eval()
 (P/variant).mkdir(exist_ok=True)
 for j,s in enumerate(states):
  st,ex=prepared[s['id']];chs=st.unique_chains();batch=collate_states([ex]);start=time.perf_counter()
  with torch.no_grad():out=model.forward_state(batch)
  mean=out['mean'][0].cpu().numpy();logvar=out['log_variance'][0].cpu().numpy();mask=out['mask'][0].cpu().numpy();duration=time.perf_counter()-start
  assert mean.shape==(st.num_residues,2) and np.isfinite(mean[mask]).all()
  geometry=ex['geometry_mask'].numpy();aas=np.array(list(st.sequence));chi_valid=mask & geometry & ~np.isin(aas,['A','G'])
  np.savez_compressed(P/variant/(s['id']+'.npz'),mean=mean,log_variance=logvar,model_mask=mask,geometry_mask=geometry,hchi_chemical_valid=chi_valid,chain=np.array(st.chain_ids),author=np.array(st.author_residue_ids),aa=aas,input_sha256=s['sha256'],checkpoint_sha256=digest)
  for i in range(st.num_residues):
   allrows.append({'state_id':s['id'],'case':'HyHEL10_2017' if s['id'].startswith('hh10') else ('existing_MD_label_check' if s['id'].startswith('mdref') else 'FlexHEL_2020'),'variant':variant,'pdb_id':s['pdb_id'],'input_kind':s['kind'],'model_index':i,'chain':st.chain_ids[i],'author_residue_id':st.author_residue_ids[i],'aa':st.sequence[i],'geometry_valid':bool(geometry[i]),'model_mask':bool(mask[i]),'rmsf':float(mean[i,0]),'hchi':float(mean[i,1]),'hchi_chemical_valid':bool(chi_valid[i]),'rmsf_log_variance':float(logvar[i,0]),'hchi_log_variance':float(logvar[i,1])})
  audit.append({'state_id':s['id'],'variant':variant,'residues':st.num_residues,'chains':chs,'input_sha256':s['sha256'],'checkpoint_sha256':digest,'seconds':duration,'finite_masked_outputs':True})
  print(variant,j+1,len(states),s['id'],st.num_residues,round(duration,3),flush=True)
 del model
with (R/'all_states_per_residue_predictions.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
impl={str(p.relative_to(SRC)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(SRC.rglob('*.py'))}
(R/'inference_audit.json').write_text(json.dumps({'status':'completed','states':len(states),'model_variants':['full','intrachain'],'forward_passes':len(audit),'rows':len(allrows),'runtime':{'torch':torch.__version__,'numpy':np.__version__,'python':sys.version,'cpu_threads':2},'source':str(SRC),'implementation_files_sha256':impl,'frozen_states_sha256':hashlib.sha256((R/'states.json').read_bytes()).hexdigest(),'frozen_contrasts_sha256':hashlib.sha256((R/'contrasts.json').read_bytes()).hexdigest(),'predictions':audit,'no_retraining':True,'no_new_MD':True,'conditions_not_model_inputs':True,'caveat':'Whole-family input audit is stored separately; counterfactual sequence-only states are not independent experimental mutant structures.'},ensure_ascii=False,indent=2),encoding='utf-8')
print('COMPLETE',len(allrows),flush=True)
