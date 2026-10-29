"""Freeze corrected fixed-epoch predictions before collecting the new MD batch."""
from pathlib import Path
import argparse,sys,json,hashlib
from collections import Counter
import numpy as np
import pandas as pd
import torch

parser=argparse.ArgumentParser()
parser.add_argument('--training-root',type=Path,required=True)
parser.add_argument('--static-root',type=Path,required=True)
parser.add_argument('--md-root',type=Path,required=True)
args=parser.parse_args()
sys.path.insert(0,str(args.training_root/'src'))
sys.path.insert(0,str(args.static_root))
sys.path.insert(0,str(args.static_root/'runtime_dependencies'))
from flexppi.config import ModelConfig
from flexppi.model.information_ablation import build_information_ablation
from flexppi.inference import _structure_example
from flexppi.data.dataset import collate_states
from flexppi.data.structures import parse_pdb,target_mask
from static_feature_definitions import features
import joblib

torch.set_num_threads(2)
out=args.md_root/'frozen_predictions'
out.mkdir(exist_ok=True)
if (out/'FROZEN.json').exists():
    raise RuntimeError('Predictions already frozen; do not overwrite')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rows=[]
models=[]
inputs={state:args.md_root/'inputs'/f'{state}.pdb' for state in ['apo','myb','mll','both']}
for fold in range(5):
    for arm in ['full','intrachain']:
        checkpoint_path=args.training_root/'runs'/arm/f'fold_{fold}'/'epoch_040.pt'
        checkpoint=torch.load(checkpoint_path,map_location='cpu',weights_only=False)
        assert checkpoint['epoch']==40 and checkpoint['config']['seed']==7
        config=ModelConfig(**checkpoint['config']['model'])
        model=build_information_ablation(config,arm).eval()
        model.load_state_dict(checkpoint['model'],strict=True)
        models.append({'method':arm,'fold':fold,'checkpoint':str(checkpoint_path),'sha256':sha(checkpoint_path)})
        for state,path in inputs.items():
            structure=parse_pdb(str(path))
            example=_structure_example(str(path),structure.unique_chains(),298,100,0)
            batch=collate_states([example])
            with torch.inference_mode(): prediction=model.forward_state(batch)
            mean=prediction['mean'][0].numpy()
            mask=prediction['mask'][0].numpy()
            for i in np.flatnonzero(mask):
                rows.append(dict(state=state,method=arm,fold=fold,chain=structure.chain_ids[i],
                    author_residue_number=int(structure.author_residue_ids[i]),amino_acid=structure.sequence[i],
                    predicted_rmsf_angstrom=float(mean[i,config.metrics.index('rmsf')]),
                    predicted_hchi_nats=float(mean[i,config.metrics.index('chi_entropy')])))
        del model
    static_models={metric:joblib.load(args.static_root/'static_corrected'/f'static_fold{fold}_{metric}.joblib')
                   for metric in ['rmsf','chi_entropy']}
    for metric in static_models:
        path=args.static_root/'static_corrected'/f'static_fold{fold}_{metric}.joblib'
        models.append(dict(method='static',metric=metric,fold=fold,path=str(path),sha256=sha(path)))
    manifest=[json.loads(s) for s in (args.training_root/'splits'/f'fold_{fold}_train.jsonl').read_text().splitlines()]
    # Match the existing source/system/record/residue weighting exactly.
    totals=np.zeros((21,2)); mass=np.zeros((21,2))
    prior_records=[]
    for record in manifest:
        with np.load(record['label_path'],allow_pickle=False) as data:
            assert list(data['metric_names'])[:2]==['rmsf','chi_entropy']
            values=data['values'][:,:2]
            mask=data['mask'][:,:2]
        structure=parse_pdb(record['structure_path'])
        valid=mask & structure.frame_mask[:,None] & target_mask(structure,record['target_chains'])[:,None] & np.isfinite(values)
        prior_records.append((record,structure.aatype,values,valid))
    for m in range(2):
        active=[item for item in prior_records if item[3][:,m].any()]
        system_records=Counter((item[0]['dataset'],item[0]['metadata'].get('biological_system',item[0]['sample_id'])) for item in active)
        source_systems=Counter(source for source,system in system_records)
        for record,aatype,values,valid in active:
            mask=valid[:,m]
            key=(record['dataset'],record['metadata'].get('biological_system',record['sample_id']))
            weight=1.0/(len(source_systems)*source_systems[key[0]]*system_records[key]*mask.sum())
            np.add.at(totals[:,m],aatype[mask],values[mask,m]*weight)
            np.add.at(mass[:,m],aatype[mask],weight)
    priors=totals/np.maximum(mass,1e-12)
    for state,path in inputs.items():
        structure=parse_pdb(str(path))
        x,names,sets,categorical=features(structure)
        ids=[names.index(name) for name in sets['sequence_complex']]
        static=np.stack([static_models[metric].predict(x[:,ids]) for metric in ['rmsf','chi_entropy']],axis=1)
        for method,pred in [('static',static),('aa',priors[structure.aatype])]:
            for i in np.flatnonzero(structure.frame_mask):
                rows.append(dict(state=state,method=method,fold=fold,chain=structure.chain_ids[i],
                    author_residue_number=int(structure.author_residue_ids[i]),amino_acid=structure.sequence[i],
                    predicted_rmsf_angstrom=float(pred[i,0]),predicted_hchi_nats=float(pred[i,1])))
frame=pd.DataFrame(rows)
assert not frame.duplicated(['state','method','fold','chain','author_residue_number']).any()
frame.to_csv(out/'all_folds_predictions.csv',index=False)
primary=frame[frame.fold==1].drop(columns='fold')
primary.to_csv(out/'residue_predictions.csv',index=False)
freeze=dict(protocol='corrected_pipeline_epoch040_fold1_primary_all_folds_sensitivity_v1',
    primary_fold=1,selection='Fold1 was already used by the v13 manuscript; frozen before the new MD batch; retain every fold in sensitivity.',
    prediction_csv_sha256=sha(out/'residue_predictions.csv'),all_folds_csv_sha256=sha(out/'all_folds_predictions.csv'),
    input_sha256={state:sha(path) for state,path in inputs.items()},models=models,
    training_freeze_sha256=sha(args.training_root/'training_freeze.json'),code_sha256=sha(__file__))
(out/'FROZEN.json').write_text(json.dumps(freeze,indent=2),encoding='utf-8')
print(json.dumps(dict(rows=len(frame),primary_rows=len(primary),output=str(out)),indent=2))
