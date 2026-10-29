"""Read-only application of two existing models with no ITrace training records."""
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
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
plan=json.loads((ROOT/'analysis_plan.json').read_text())
pdbs=list(dict.fromkeys(p for sections in [plan['series'],plan.get('additional_biology_verified_comparisons',{})] for group in sections.values() for p in group))
device=torch.device('cuda')
models={
 'family_full_no_ITrace':(Path('/root/flexppi_completion_20261010/family_folds/runs/full/fold_0/epoch_040.pt'),'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5'),
 'old_strict_full_no_ITrace':(Path('/root/flexppi_revision_20261002/corrected_training_v1/runs/strict_full/fold_0/epoch_040.pt'),'b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f')}
audits=[]
for name,(check,expectedsha) in models.items():
    assert sha(check)==expectedsha
    c=torch.load(check,map_location='cpu',weights_only=False)
    assert c['config']['model']['metrics']==['rmsf','chi_entropy']
    model=build_information_ablation(ModelConfig(**c['config']['model']),'full').to(device).eval()
    model.load_state_dict(c['model'],strict=True)
    out=ROOT/name;out.mkdir(exist_ok=True)
    for pdb in pdbs:
        path=Path('/root/flexppi_revision_20261002/corrected_training_v1/dataset/structures')/f'itrace_{pdb}.pdb'
        ex=_structure_example(str(path),None,310.,200.,2.)
        batch=move_to_device(collate_states([ex]),device)
        with torch.inference_mode():mean=model.forward_state(batch)['mean'][0].cpu().numpy()
        st=parse_pdb(str(path))
        assert len(mean)==len(st.sequence)
        np.savez_compressed(out/f'itrace_{pdb}.npz',mean=mean,chain=np.asarray(st.chain_ids),author_residue_id=np.asarray(st.author_residue_ids),amino_acid=np.asarray(list(st.sequence)))
        rec={'model':name,'pdb':pdb,'checkpoint_sha256':expectedsha,'input_sha256':sha(path),'residues':len(mean)}
        if name=='family_full_no_ITrace':
            old=np.load(Path('/root/flexppi_completion_20261010/family_folds/evaluation/full')/f'itrace_{pdb}.npz')
            assert np.array_equal(st.chain_ids,old['chain'])
            assert np.array_equal(st.author_residue_ids,old['author_residue_id'])
            err=float(np.max(abs(mean-old['mean'])))
            rec['max_error_vs_existing_cache']=err
            assert err<5e-5,(pdb,err)
        audits.append(rec)
        print(name,pdb,len(mean),flush=True)
    del model,c
(ROOT/'fresh_inference_audit.json').write_text(json.dumps({'status':'passed','locked_series':plan['series'],'forward_passes':len(audits),'new_training':False,'new_MD':False,'predictions':audits},indent=2))
print('All inference complete:',len(audits),flush=True)
