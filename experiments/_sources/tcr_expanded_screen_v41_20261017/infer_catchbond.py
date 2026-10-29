from pathlib import Path
import sys,json,hashlib,csv
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,'/root/flexppi_revision_20261002/corrected_training_v1/src')
import numpy as np
import torch
from flexppi.config import ModelConfig
from flexppi.data.dataset import collate_states
from flexppi.data.structures import parse_pdb
from flexppi.inference import _structure_example
from flexppi.model.information_ablation import build_information_ablation
from flexppi.trainer import move_to_device
torch.set_num_threads(2)
check=Path('/root/flexppi_completion_20261010/family_folds/runs/full/fold_0/epoch_040.pt')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(check)=='ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5'
c=torch.load(check,map_location='cpu',weights_only=False)
assert c['config']['model']['metrics']==['rmsf','chi_entropy']
model=build_information_ablation(ModelConfig(**c['config']['model']),'full').eval();model.load_state_dict(c['model'],strict=True)
rows=[];audit=[]
for cond in ['WT','S32M']:
    for kind in ['prepared','experimental_matched']:
        path=ROOT/f'catchbond_{cond}_{kind}.pdb'
        ex=_structure_example(str(path),None,300.,1000.,2.)
        batch=move_to_device(collate_states([ex]),torch.device('cpu'))
        with torch.inference_mode():mean=model.forward_state(batch)['mean'][0].numpy()
        st=parse_pdb(str(path))
        np.savez_compressed(ROOT/f'catchbond_{cond}_{kind}_prediction.npz',mean=mean,chain=st.chain_ids,author_residue_id=st.author_residue_ids,amino_acid=np.asarray(list(st.sequence)))
        for i in range(len(mean)):rows.append({'condition':cond,'input_kind':kind,'chain':st.chain_ids[i],'residue_id':st.author_residue_ids[i],'aa':st.sequence[i],'rmsf_pred':float(mean[i,0]),'Hchi_pred':float(mean[i,1])})
        audit.append({'condition':cond,'kind':kind,'input_sha256':sha(path),'checkpoint_sha256':sha(check),'residues':len(mean)})
with (ROOT/'catchbond_predictions.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(ROOT/'inference_audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit,indent=2))
