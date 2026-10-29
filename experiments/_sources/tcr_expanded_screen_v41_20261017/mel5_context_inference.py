"""Known MEL5 complexes: fixed pMHC coordinates, remove TCR inputs only."""
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
model=build_information_ablation(ModelConfig(**c['config']['model']),'full').eval()
model.load_state_dict(c['model'],strict=True)
rows=[];audit={'checkpoint_sha256':sha(check),'intervention':'Remove D/E TCR from input; preserve every A/B/C pMHC atom coordinate, sequence, and chain order. This probes model input dependence; not a physical unbinding experiment.','full_reproduction':{}}
for pdb in ['3hg1','4jff']:
    full=ROOT/f'{pdb}_model_input.pdb'
    pmhc=ROOT/f'{pdb}_pMHC_only.pdb'
    lines=full.read_text().splitlines(True)
    keep=[l for l in lines if l[:6].strip() in ['ATOM','HETATM','TER'] and l[21:22] in 'ABC']
    pmhc.write_text(''.join(keep)+'END\n')
    assert [l for l in pmhc.read_text().splitlines(True) if l.startswith(('ATOM  ','HETATM'))]==[l for l in lines if l.startswith(('ATOM  ','HETATM')) and l[21:22] in 'ABC']
    for kind,path in [('full',full),('pMHC_only',pmhc)]:
        ex=_structure_example(str(path),None,300.,200.,2.)
        batch=move_to_device(collate_states([ex]),torch.device('cpu'))
        with torch.inference_mode():mean=model.forward_state(batch)['mean'][0].numpy()
        st=parse_pdb(str(path))
        np.savez_compressed(ROOT/f'{pdb}_{kind}_context_prediction.npz',mean=mean,chain=st.chain_ids,author_residue_id=st.author_residue_ids,amino_acid=np.asarray(list(st.sequence)))
        for i in range(len(mean)):
            rows.append({'pdb':pdb,'input_kind':kind,'chain':st.chain_ids[i],'residue_id':st.author_residue_ids[i],'aa':st.sequence[i],'rmsf_pred':float(mean[i,0]),'Hchi_pred':float(mean[i,1])})
        if kind=='full':
            old=np.load(ROOT/f'itrace_{pdb}.npz')
            assert mean.shape==old['mean'].shape
            err=float(np.max(abs(mean-old['mean'])))
            assert err<2e-5,(pdb,err)
            audit['full_reproduction'][pdb]={'max_absolute_error':err,'input_sha256':sha(path),'residues':len(mean)}
        else:audit[pdb]={'pMHC_sha256':sha(path),'pMHC_residues':len(mean)}
        print(pdb,kind,len(mean),flush=True)
with (ROOT/'mel5_context_predictions.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
audit['status']='passed'
(ROOT/'mel5_context_inference_audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit,indent=2),flush=True)
