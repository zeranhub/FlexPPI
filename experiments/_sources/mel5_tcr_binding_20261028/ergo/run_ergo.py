"""Inference-only adapter for official ERGO-II models; no weight/architecture edits.
Run with --input official-schema.csv --output-dir result_dir.
Old Lightning 0.7.1 is replaced only by nn.Module base and no-op training decorators.
"""
from pathlib import Path
import sys, os, json, csv, hashlib, types, argparse, pickle, shutil
import torch

ROOT=Path(__file__).resolve().parent
UPSTREAM=ROOT/'download'/'ERGO-II-master'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
# The upstream model loader points to TCR_Autoencoder, while released weights live in Models/AE.
for p in (UPSTREAM/'Models'/'AE').glob('*.pt'):
    dest=UPSTREAM/'TCR_Autoencoder'/p.name
    dest.parent.mkdir(exist_ok=True)
    if not dest.exists():shutil.copyfile(p,dest)
    assert sha(p)==sha(dest)
# ERGO-II only uses the Lightning superclass as torch.nn.Module during inference.
pl=types.ModuleType('pytorch_lightning')
pl.LightningModule=torch.nn.Module
pl.Trainer=type('TrainingUnavailable',(),{})
pl.data_loader=lambda function:function
logging=types.ModuleType('pytorch_lightning.logging')
logging.TensorBoardLogger=type('TrainingUnavailable',(),{})
callbacks=types.ModuleType('pytorch_lightning.callbacks')
callbacks.EarlyStopping=type('TrainingUnavailable',(),{})
sys.modules['pytorch_lightning']=pl
sys.modules['pytorch_lightning.logging']=logging
sys.modules['pytorch_lightning.callbacks']=callbacks
# Checkpoint tensors were saved on GPU; only placement is changed to CPU.
original_load=torch.load
def cpu_load(*args,**kwargs):
    kwargs.setdefault('map_location','cpu')
    return original_load(*args,**kwargs)
torch.load=cpu_load
os.chdir(UPSTREAM)
sys.path.insert(0,str(UPSTREAM))
import Predict
from Loader import SignedPairsDataset
from torch.utils.data import DataLoader

def run(input_path,out_dir):
    import pandas as pd, numpy as np, sklearn
    torch.set_num_threads(2)
    torch.manual_seed(0)
    out_dir.mkdir(parents=True,exist_ok=True)
    incoming=pd.read_csv(input_path,keep_default_na=False)
    expected=['TRA','TRB','TRAV','TRAJ','TRBV','TRBJ','T-Cell-Type','Peptide','MHC']
    assert all(c in incoming.columns for c in expected)
    audit={'official_repository':'https://github.com/IdoSpringer/ERGO-II','upstream_commit':'85d320ab03ade33460be9a81ea3a51b8d37cd998',
      'inference_modifications':['LightningModule replaced by torch.nn.Module for inference; no training operations invoked.','Lightning training-only decorator and import classes replaced by inert import placeholders.','torch.load map_location=cpu; weight state_dict loaded by untouched upstream Predict.load_model with strict=True default.','Official Models/AE weight files copied to upstream hardcoded TCR_Autoencoder path without byte changes.'],
      'runtime':{'torch':torch.__version__,'numpy':np.__version__,'pandas':pd.__version__,'sklearn':sklearn.__version__},
      'input':str(input_path),'input_sha256':sha(input_path),'score_interpretation':'Official 0-1 binding classifier score; not a calibrated KD or measured probability. Higher means greater binding propensity under model.',
      'threshold':'No fixed binary cutoff is endorsed by upstream prediction help; 0.5 is recorded only as an untuned reference, not selected on this panel.',
      'models':{}}
    all_rows=[]
    for dataset in ['mcpas','vdjdb']:
        model, train_file=Predict.get_model(dataset)
        model.eval()
        dictionaries=Predict.get_train_dicts(train_file)
        frame=incoming.copy()
        vocab_audit=[]
        for col,d in zip(['TRAV','TRBV','TRAJ','TRBJ','MHC'],dictionaries):
            for val in sorted(set(frame[col])):
                normalized=val
                if col=='MHC' and val not in d and val.startswith('HLA-A*02:') and 'HLA-A*02' in d:
                    normalized='HLA-A*02'
                    frame.loc[frame[col]==val,col]=normalized
                vocab_audit.append({'column':col,'original':val,'normalized':normalized,'known_in_training_vocabulary':normalized in d,'index':d.get(normalized,0)})
        model_input=out_dir/f'input_{dataset}.csv'
        frame[expected].to_csv(model_input,index=False)
        samples, official_df=Predict.read_input_file(model_input)
        assert len(samples)==len(incoming),'Invalid rows must not be silently discarded'
        ds=SignedPairsDataset(samples,dictionaries)
        loader=DataLoader(ds,batch_size=1000,shuffle=False,collate_fn=lambda b:ds.collate(b,tcr_encoding=model.tcr_encoding_model,cat_encoding=model.cat_encoding))
        with torch.no_grad():
            scores=[]
            for i,batch in enumerate(loader): scores.extend(model.validation_step(batch,i)['y_hat'].reshape(-1).tolist())
            singleton=[]
            for i in range(len(ds)):
                batch=ds.collate([ds[i]],tcr_encoding=model.tcr_encoding_model,cat_encoding=model.cat_encoding)
                singleton.extend(model.validation_step(batch,i)['y_hat'].reshape(-1).tolist())
        maxerr=max(abs(a-b) for a,b in zip(scores,singleton))
        assert maxerr<1e-6,(dataset,maxerr)
        train=pickle.load(open(train_file,'rb'))
        overlaps=[]
        for i,s in enumerate(samples):
            exact=[t for t in train if t['tcra']==s['tcra'] and t['tcrb']==s['tcrb'] and t['peptide']==s['peptide']]
            beta=[t for t in train if t['tcrb']==s['tcrb'] and t['peptide']==s['peptide']]
            antigen=[t for t in train if t['peptide']==s['peptide']]
            overlaps.append({'row':i,'peptide':s['peptide'],'same_alpha_beta_peptide_train_rows':len(exact),'same_beta_peptide_train_rows':len(beta),'same_peptide_train_rows':len(antigen)})
        output=incoming.copy()
        output['model']='ERGO-II_'+dataset
        output['binding_score']=scores
        output['predicted_binding_at_reference_0_5']=[s>=.5 for s in scores]
        output.to_csv(out_dir/f'predictions_{dataset}.csv',index=False,float_format='%.12g')
        all_rows.extend(output.to_dict('records'))
        version='1meajht' if dataset=='mcpas' else '1veajht'
        weights=list((UPSTREAM/'Models'/f'version_{version}'/'checkpoints').iterdir())
        audit['models'][dataset]={
          'version':version,'hparams':vars(model.hparams),'checkpoint':str(weights[0]),'checkpoint_sha256':sha(weights[0]),
          'ae_weights':{p.name:sha(p) for p in (UPSTREAM/'Models'/'AE').glob('*.pt')},
          'train_samples_sha256':sha(train_file),'upstream_source_sha256':{name:sha(UPSTREAM/name) for name in ['Predict.py','Trainer.py','Models.py','Loader.py']},
          'vocabulary_audit':vocab_audit,'train_overlap':overlaps,'all_rows_retained':True,'singleton_batch_max_abs_diff':maxerr,
          'strict_checkpoint_loading':True,'state_parameter_count':sum(p.numel() for p in model.parameters())}
    pd.DataFrame(all_rows).to_csv(out_dir/'all_ergo_scores.csv',index=False,float_format='%.12g')
    (out_dir/'ergo_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(pd.DataFrame(all_rows)[['model','TRA','TRB','Peptide','binding_score']].to_string(index=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    run(args.input.resolve(),args.output_dir.resolve())
