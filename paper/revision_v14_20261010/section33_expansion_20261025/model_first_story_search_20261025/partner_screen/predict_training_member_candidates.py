"""Fixed fold 0 current-model predictions for exploratory training-member cases.

No retraining. Fold 0 chosen once for all states, no per-case best-fold search.
Models see no label data during this program. Evaluation is separate.
"""
import os
for key in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']: os.environ[key]='2'
os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
import sys,json,hashlib,time
import numpy as np
import pandas as pd
import torch
CORE=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
OUT=Path(__file__).resolve().parent
R=OUT.parent.parent
sys.path.insert(0,str(CORE/'training_payload/src'))
from flexppi.config import ModelConfig
from flexppi.constants import METRIC_NAMES
from flexppi.data.dataset import collate_states
from flexppi.data.structures import parse_pdb,target_mask
from flexppi.model.information_ablation import build_information_ablation
torch.set_num_threads(2);torch.set_num_interop_threads(1)

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    t=time.time(); dest=OUT/'training_member_predictions';dest.mkdir(exist_ok=True)
    inventory=pd.read_csv(R/'paired_and_experimental/paired_inventory.csv')
    selected=inventory[inventory.pdb_id.isin(['2Y5B','4ZRJ','2XWT','2I3T','3QML'])]
    ids=sorted(set(selected.apo_sample)|set(selected.holo_sample)|{'A001K','A004C','A007U'})
    recs={r['sample_id']:r for r in map(json.loads,(CORE/'training_payload/splits/all_ready.jsonl').read_text().splitlines())}
    inputs=[];prepared=[]
    for sid in ids:
        p=CORE/'training_payload/dataset/structures'/(sid+'.pdb')
        s=parse_pdb(str(p));n=s.num_residues
        b=dict(sample_id=sid,dataset='exploratory_training_member',source_group='exploratory_training_member',metadata={},
            aatype=torch.from_numpy(s.aatype),residue_index=torch.from_numpy(s.residue_index),chain_index=torch.from_numpy(s.chain_index),
            backbone=torch.from_numpy(s.backbone),backbone_mask=torch.from_numpy(s.backbone_mask),geometry_mask=torch.from_numpy(s.frame_mask),
            target_mask=torch.from_numpy(target_mask(s,[])),condition=torch.zeros(3),
            label_values=torch.full((n,len(METRIC_NAMES)),float('nan')),label_variances=torch.zeros(n,len(METRIC_NAMES)),label_mask=torch.zeros(n,len(METRIC_NAMES),dtype=torch.bool))
        prepared.append((sid,s,collate_states([b])))
        fold=recs[sid]['metadata']['fold']
        inputs.append(dict(sample_id=sid,path=str(p),sha256=sha(p),n=n,original_fold=fold,present_in_fold0_training=fold!=0))
    (OUT/'training_member_input_audit.json').write_text(json.dumps(inputs,indent=2),encoding='utf-8')
    print('Prepared',ids,flush=True)
    manifest=[]
    for variant in ['full','intrachain']:
        ck=CORE/'trained_checkpoints'/variant/'fold_0/epoch_040.pt'
        checkpoint=torch.load(ck,map_location='cpu',weights_only=False);mc=checkpoint['config']['model']
        assert checkpoint['epoch']==40 and mc['information_ablation']==variant
        model=build_information_ablation(ModelConfig(**mc),variant).cpu().eval();model.load_state_dict(checkpoint['model'],strict=True);model.requires_grad_(False)
        manifest.append(dict(variant=variant,fold=0,path=str(ck),sha256=sha(ck)))
        for sid,s,b in prepared:
            with torch.inference_mode(): result=model.forward_state(b)
            v=result['mean'][0].numpy();mask=result['mask'][0].numpy()
            assert np.isfinite(v).all() and np.array_equal(mask,s.frame_mask)
            np.savez_compressed(dest/(variant+'_'+sid+'.npz'),mean=v,mask=mask,chain=np.array(s.chain_ids),author_residue_id=np.array(s.author_residue_ids),amino_acid=np.array(list(s.sequence)))
            print('DONE',variant,sid,s.num_residues,flush=True)
        del model,checkpoint
    (OUT/'training_member_completion.json').write_text(json.dumps(dict(status='complete',elapsed_seconds=time.time()-t,inputs=inputs,checkpoints=manifest,script_sha256=sha(Path(__file__)),scope='Fixed fold 0; most states occur in its training set. No independent validation claim.'),indent=2),encoding='utf-8')
    print('COMPLETE',time.time()-t,flush=True)
if __name__=='__main__': main()
