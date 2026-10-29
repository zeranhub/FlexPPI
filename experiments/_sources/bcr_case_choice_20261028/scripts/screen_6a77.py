from pathlib import Path
import csv,hashlib,json,sys
import numpy as np,torch
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
S=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\src')
sys.path.insert(0,str(S))
from flexppi.inference import _structure_example
from flexppi.data.structures import parse_pdb
from flexppi.data.dataset import collate_states
from flexppi.config import ModelConfig
from flexppi.model.information_ablation import build_information_ablation
torch.set_num_threads(2)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
p=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\dataset\structures\A00LZ.pdb')
lp=p.parent.parent/'labels/A00LZ.npz';st=parse_pdb(str(p));lab=np.load(lp,allow_pickle=False)
names=lab['metric_names'].tolist();assert st.num_residues==len(lab['values'])==519
rows=[];aud=[]
for v,expected in [('full','ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5'),('intrachain','6b069d9d3521e95f5f2d63a3f19d64fb4071862efa4682473387b0ca45ff73c3')]:
 cp_path=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028\runtime\family_folds\runs')/v/'fold_0/epoch_040.pt'
 assert sha(cp_path)==expected
 cp=torch.load(cp_path,map_location='cpu',weights_only=False);cfg=ModelConfig(**cp['config']['model'])
 m=build_information_ablation(cfg,v).eval();m.load_state_dict(cp['model'],strict=True)
 ex=_structure_example(str(p),st.unique_chains(),300.,100.,0.)
 with torch.no_grad():o=m.forward_state(collate_states([ex]))
 mean=o['mean'][0].numpy();mask=o['mask'][0].numpy();geo=ex['geometry_mask'].numpy()
 np.savez_compressed(R/f'6A77_{v}_predictions.npz',mean=mean,log_variance=o['log_variance'][0].numpy(),model_mask=mask,geometry_mask=geo,chain=st.chain_ids,author=st.author_residue_ids,aa=list(st.sequence),input_sha256=sha(p),checkpoint_sha256=expected)
 for i in range(st.num_residues):
  row={'variant':v,'model_index':i,'chain':st.chain_ids[i],'author':st.author_residue_ids[i],'aa':st.sequence[i],'geometry_valid':bool(geo[i]),'model_valid':bool(mask[i])}
  for j,metric in enumerate(['rmsf','chi_entropy']):
   k=names.index(metric);row[metric+'_pred']=float(mean[i,j]);row[metric+'_MD']=float(lab['values'][i,k]);row[metric+'_valid']=bool(lab['mask'][i,k]&mask[i]&geo[i])
  rows.append(row)
 aud.append({'variant':v,'checkpoint':str(cp_path),'checkpoint_sha256':expected,'strict_load':True})
 del m
with (R/'6A77_all_positions_predictions_and_MD.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(R/'audits/current_inference.json').write_text(json.dumps({'status':'completed','sample':'A00LZ','residues':519,'input':str(p),'input_sha256':sha(p),'label':str(lp),'label_sha256':sha(lp),'label_names':names,'models':aud,'chains':{c:st.chain_sequence(c) for c in st.unique_chains()},'important':'Frozen519-residue version, not the legacy514-residue local cache. No rowwise mixing of different structure versions.'},indent=2),encoding='utf-8')
print('COMPLETE 6A77',st.num_residues,st.unique_chains(),flush=True)
