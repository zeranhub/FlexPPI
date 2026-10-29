from pathlib import Path
import hashlib,json,os,sys
os.environ['TF_CPP_MIN_LOG_LEVEL']='2';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd,tensorflow as tf
R=Path(r'E:\Project\FlexPPI\experiments\_sources\antigen_incremental_value_20261028')
U=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\nettcr\upstream')
sys.path.insert(0,str(U/'src'));import keras_utils
out=R/'baselines/nettcr';out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
df=pd.read_csv(R/'baseline_inputs.csv').fillna('');cols=['peptide','A1','A2','A3','B1','B2','B3'];maxima=[12,7,8,22,6,7,23]
train=pd.read_csv(U/'data/nettcr_2_2_limited_dataset.csv')
selections={};overlaps=[]
for system,g in df.groupby('system'):
 pair=(g.A3.iloc[0],g.B3.iloc[0]);m=train[(train.A3==pair[0])&(train.B3==pair[1])].copy()
 m['query_system']=system;overlaps.append(m)
 counts=m.partition.value_counts().to_dict()
 model_candidates=[]
 for t in range(5):
  for val in range(5):
   if t==val:continue
   train_overlap=sum(int(n) for part,n in counts.items() if part not in [t,val]);val_overlap=int(counts.get(val,0))
   model_candidates.append((train_overlap,val_overlap,t,val))
 best=min((x[0],x[1]) for x in model_candidates)
 chosen=[f't.{t}.v.{val}' for tr,va,t,val in model_candidates if (tr,va)==best]
 selections[system]={'known_pair_partition_counts':{str(k):int(v) for k,v in counts.items()},'selected_models':chosen,'minimum_training_pair_overlap':best[0],'minimum_validation_pair_overlap':best[1], 'rule':'Minimize number of exact paired-CDR3 rows in training, then validation, before predictions; average all tied models.'}
pd.concat(overlaps,ignore_index=True).to_csv(out/'training_pair_overlap.csv',index=False)
(out/'selection_before_predictions.json').write_text(json.dumps(selections,indent=2),encoding='utf-8')
for c,n in zip(cols,maxima):assert df[c].str.len().between(1,n).all(),c
data={key:np.float32(keras_utils.enc_list_bl_max_len(df[col],keras_utils.blosum50_20aa,n)/5) for key,col,n in zip(['pep','a1','a2','a3','b1','b2','b3'],cols,maxima)}
tf.config.threading.set_inter_op_parallelism_threads(1);tf.config.threading.set_intra_op_parallelism_threads(1)
predictions={};aud=[]
for path in sorted((U/'models/nettcr_2_2_pan/checkpoint').glob('*.tflite')):
 interpreter=tf.lite.Interpreter(model_path=str(path),num_threads=1);inputs=interpreter.get_input_details();outputs=interpreter.get_output_details()
 interpreter.resize_tensor_input(outputs[0]['index'],[len(df),outputs[0]['shape'][1]])
 for item in inputs:interpreter.resize_tensor_input(item['index'],[len(df),item['shape'][1],item['shape'][2]])
 interpreter.allocate_tensors()
 for item in inputs:interpreter.set_tensor(item['index'],data[item['name'].split(':')[0].split('_')[-1]])
 interpreter.invoke();p=interpreter.get_tensor(outputs[0]['index']).reshape(-1)
 assert np.isfinite(p).all() and ((p>=0)&(p<=1)).all();predictions[path.stem]=p;aud.append({'model':path.stem,'sha256':sha(path)})
 assert len(p)==279
raw=df[['candidate_id','system','peptide']].copy()
for name,p in predictions.items():raw[name]=p
raw.to_csv(out/'all_model_scores.csv',index=False)
result=df[['candidate_id','system','peptide']].copy()
result['NetTCR22_selected']=np.nan;result['NetTCR22_all20']=np.mean(list(predictions.values()),axis=0)
for system,selection in selections.items():
 mask=df.system==system;p=np.mean([predictions[n] for n in selection['selected_models']],axis=0);result.loc[mask,'NetTCR22_selected']=p[mask]
result.to_csv(out/'scores.csv',index=False)
(out/'audit.json').write_text(json.dumps({'rows':len(df),'official_model_count':len(aud),'all_rows_retained':True,'weights':aud,'input_sha256':sha(R/'baseline_inputs.csv'),'training_csv_sha256':sha(U/'data/nettcr_2_2_limited_dataset.csv'),'selections':selections,'score_is_not_affinity_or_function_probability':True,'training_overlap_limit':'Exact paired CDR3 list audit, not full homologous-family isolation. Functional design labels not used in checkpoint selection.','runtime':{'python':sys.version,'tensorflow':tf.__version__,'numpy':np.__version__}},indent=2),encoding='utf-8')
print(json.dumps({'rows':len(result),'selection':selections},indent=2))
