from pathlib import Path
import csv,hashlib,json,sys,time
import numpy as np,torch
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
S=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\src')
sys.path.insert(0,str(S))
from flexppi.inference import _structure_example
from flexppi.data.structures import parse_pdb
from flexppi.data.dataset import collate_states
from flexppi.config import ModelConfig
from flexppi.model.information_ablation import build_information_ablation
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
torch.set_num_threads(2);torch.manual_seed(7)
states=json.loads((R/'pair_screen_states.json').read_text());assert len(states)==15
rows=[];aud=[];start=time.monotonic()
for v,expected in [('full','ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5'),('intrachain','6b069d9d3521e95f5f2d63a3f19d64fb4071862efa4682473387b0ca45ff73c3')]:
 cp_path=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028\runtime\family_folds\runs')/v/'fold_0/epoch_040.pt';assert sha(cp_path)==expected
 cp=torch.load(cp_path,map_location='cpu',weights_only=False);cfg=ModelConfig(**cp['config']['model']);assert list(cfg.metrics)==['rmsf','chi_entropy']
 m=build_information_ablation(cfg,v).eval();m.load_state_dict(cp['model'],strict=True)
 P=R/'paired_predictions'/v;P.mkdir(parents=True,exist_ok=True)
 for s in states:
  assert sha(s['path'])==s['sha256'];st=parse_pdb(s['path']);assert st.num_residues<=1024
  assert set(st.unique_chains())==set(s['chains'])
  for ch,record in s['chain_records'].items():assert st.chain_sequence(ch)==record['sequence']
  ex=_structure_example(s['path'],st.unique_chains(),300.,100.,0.)
  with torch.no_grad():o=m.forward_state(collate_states([ex]))
  mean=o['mean'][0].numpy();mask=o['mask'][0].numpy();geo=ex['geometry_mask'].numpy();chi=mask&geo&~np.isin(list(st.sequence),['A','G'])
  assert np.isfinite(mean[mask]).all()
  np.savez_compressed(P/(s['id']+'.npz'),mean=mean,log_variance=o['log_variance'][0].numpy(),model_mask=mask,geometry_mask=geo,hchi_valid=chi,chain=st.chain_ids,author=st.author_residue_ids,aa=list(st.sequence),input_sha256=s['sha256'],checkpoint_sha256=expected)
  for i in range(st.num_residues):rows.append({'state_id':s['id'],'variant':v,'pdb':s['pdb'],'input_kind':s['kind'],'model_index':i,'chain':st.chain_ids[i],'author':st.author_residue_ids[i],'aa':st.sequence[i],'geometry_valid':bool(geo[i]),'model_valid':bool(mask[i]),'hchi_valid':bool(chi[i]),'rmsf':float(mean[i,0]),'hchi':float(mean[i,1])})
  aud.append({'state_id':s['id'],'variant':v,'residues':st.num_residues,'chains':st.unique_chains(),'input_sha256':s['sha256'],'checkpoint_sha256':expected,'finite':True})
  print(v,s['id'],st.num_residues,flush=True)
 del m
with (R/'paired_all_residue_predictions.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(R/'audits/paired_inference.json').write_text(json.dumps({'status':'complete','states':len(states),'forward_passes':len(aud),'residue_rows':len(rows),'strict_load':True,'contract_sha256':sha(R/'pair_screen_contract.json'),'states_sha256':sha(R/'pair_screen_states.json'),'source':str(S),'torch':torch.__version__,'numpy':np.__version__,'seconds':time.monotonic()-start,'entries':aud},indent=2),encoding='utf-8')
print('COMPLETE',len(rows),flush=True)
