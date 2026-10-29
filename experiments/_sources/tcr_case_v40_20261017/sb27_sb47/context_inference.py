"""Same joint model, fixed input coordinates, with/without TCR chains."""
from pathlib import Path
import sys, json, hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,'/root/flexppi_revision_20261002/corrected_training_v1/src')
import torch
from flexppi.config import ModelConfig
from flexppi.data.dataset import collate_states
from flexppi.data.structures import parse_pdb
from flexppi.inference import _structure_example
from flexppi.model.information_ablation import build_information_ablation
from flexppi.trainer import move_to_device
torch.set_num_threads(2)
checkpoint=Path('/root/flexppi_completion_20261010/family_folds/runs/full/fold_0/epoch_040.pt')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(checkpoint)=='ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5'
try:c=torch.load(checkpoint,map_location='cpu',weights_only=False)
except TypeError:c=torch.load(checkpoint,map_location='cpu')
assert c['epoch']==40 and c['config']['model']['metrics']==['rmsf','chi_entropy']
model=build_information_ablation(ModelConfig(**c['config']['model']),'full').eval()
model.load_state_dict(c['model'],strict=True)
rows=[]; audit={'checkpoint_sha256':sha(checkpoint),'full_reproduction':{},'intervention':'Remove D/E TCR atoms; identical A/B/C coordinates and chain order retained; same model.'}
for pdb in ['2ak4','4jry']:
    full=ROOT/(pdb+'_model_input.pdb')
    pmhc=ROOT/(pdb+'_pMHC_only.pdb')
    pmhc.write_text(''.join(l for l in full.read_text().splitlines(True) if l[:6].strip() in ['ATOM','HETATM','TER'] and l[21:22] in 'ABC')+'END\n')
    arrays={}
    for name,path in [('full',full),('pMHC_only',pmhc)]:
        ex=_structure_example(str(path),None,300.,200.,2.)
        batch=move_to_device(collate_states([ex]),torch.device('cpu'))
        with torch.inference_mode():pred=model.forward_state(batch)['mean'][0].numpy()
        st=parse_pdb(str(path)); mask=np.asarray(st.chain_ids)=='C'
        arrays[name]=pred[mask]
        np.savez_compressed(ROOT/(pdb+'_'+name+'_inference.npz'),mean=pred,chain=st.chain_ids,aa=np.asarray(list(st.sequence)),author_residue_id=st.author_residue_ids)
        for i,(aa,v) in enumerate(zip(np.asarray(list(st.sequence))[mask],pred[mask]),1):
            rows.append(dict(pdb=pdb,complex={'2ak4':'SB27','4jry':'SB47'}[pdb],input=name,position=i,aa=str(aa),rmsf=float(v[0]),Hchi=float(v[1])))
        if name=='full':
            old=np.load(ROOT/('itrace_'+pdb+'.npz'))
            error=float(np.max(np.abs(pred-old['mean'])))
            assert error<2e-5,(pdb,error)
            audit['full_reproduction'][pdb]=dict(max_absolute_error=error,input_sha256=sha(path))
    audit[pdb]={'Q7_full':arrays['full'][6].tolist(),'Q7_pMHC_only':arrays['pMHC_only'][6].tolist()}
import csv
with (ROOT/'context_inference.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
audit['status']='passed'
(ROOT/'context_inference_audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit,indent=2))
