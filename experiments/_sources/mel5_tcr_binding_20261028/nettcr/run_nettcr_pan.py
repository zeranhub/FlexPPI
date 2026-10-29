from pathlib import Path
import os,sys,json,hashlib,argparse,platform,subprocess
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL','2')
os.environ.setdefault('OMP_NUM_THREADS','1')
p=argparse.ArgumentParser();p.add_argument('--input',required=True);a=p.parse_args()
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'upstream/src'))
import numpy as np
import pandas as pd
import tensorflow as tf
import keras_utils

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,value):Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
input_path=Path(a.input);df=pd.read_csv(input_path).fillna('')
cols=['peptide','A1','A2','A3','B1','B2','B3'];maxima=[12,7,8,22,6,7,23]
assert set(cols)<=set(df)
for col,maxlen in zip(cols,maxima):
 assert df[col].str.len().between(1,maxlen).all(),(col,df[col].str.len().tolist())
 assert df[col].str.fullmatch('[ARNDCQEGHILKMFPSTWYV]+').all(),col
# Keep common input and exact model-specific input separately.
df.to_csv(root/'input_panel_with_metadata.csv',index=False)
df[cols].to_csv(root/'input_nettcr.csv',index=False)
# External target-identity match was determined before inference.
train=pd.read_csv(root/'upstream/data/nettcr_2_2_limited_dataset.csv')
cv=pd.read_csv(root/'upstream/models/nettcr_2_2_pan/cv_pred_df.csv')
known_pairs={(r.A3,r.B3) for _,r in df.iterrows()}
identity=train.apply(lambda r:(r.A3,r.B3) in known_pairs,axis=1)
matched=train[identity].copy();matched.to_csv(root/'known_CDR3_pair_training_overlap.csv',index=False)
assert set(matched.partition)=={0},set(matched.partition)
cvsel=cv[cv.apply(lambda r:(r.A3,r.B3) in known_pairs,axis=1)].copy()
assert set(cvsel.partition)=={0}
cvsel.to_csv(root/'reference_CV_rows_for_validation.csv',index=False)
combined=pd.concat([df[cols],cvsel[cols]],ignore_index=True)
data={key:np.float32(keras_utils.enc_list_bl_max_len(combined[col],keras_utils.blosum50_20aa,n)/5) for key,col,n in zip(['pep','a1','a2','a3','b1','b2','b3'],cols,maxima)}
tf.config.threading.set_inter_op_parallelism_threads(1)
tf.config.threading.set_intra_op_parallelism_threads(1)
model_dir=root/'upstream/models/nettcr_2_2_pan/checkpoint'
raw=[];models=[];fold_predictions=[];all_predictions=[]
for path in sorted(model_dir.glob('*.tflite')):
 name=path.stem;parts=name.split('.');t=int(parts[1]);v=int(parts[3])
 interpreter=tf.lite.Interpreter(model_path=str(path),num_threads=1)
 inputs=interpreter.get_input_details();outputs=interpreter.get_output_details()
 # Same dynamic batching and input-name mapping as official src/predict.py.
 interpreter.resize_tensor_input(outputs[0]['index'],[len(combined),outputs[0]['shape'][1]])
 for item in inputs:
  interpreter.resize_tensor_input(item['index'],[len(combined),item['shape'][1],item['shape'][2]])
 interpreter.allocate_tensors()
 for item in inputs:
  key=item['name'].split(':')[0].split('_')[-1]
  interpreter.set_tensor(item['index'],data[key])
 interpreter.invoke();pred=interpreter.get_tensor(outputs[0]['index']).reshape(-1)
 assert np.isfinite(pred).all() and ((pred>=0)&(pred<=1)).all()
 all_predictions.append(pred)
 if t==0:fold_predictions.append(pred)
 models.append({'model':name,'sha256':sha(path),'test_partition':t,'validation_partition':v,'training_partitions':[x for x in range(5) if x not in [t,v]],'primary_excludes_known_CDR3_pair_partition':t==0})
 for i,row in df.iterrows():
  r=row.to_dict();r.update({'model':name,'prediction':float(pred[i]),'test_partition':t,'validation_partition':v,'primary_model':t==0});raw.append(r)
npall=np.stack(all_predictions);npfold=np.stack(fold_predictions)
assert len(all_predictions)==20 and len(fold_predictions)==4
# Native float32 summation reproduces the official prediction script's averaging.
mean20=sum(all_predictions)/20
mean4=sum(fold_predictions)/4
summary=df.copy();summary['NetTCR22_pan_known_overlap_excluded_4_mean']=mean4[:len(df)];summary['NetTCR22_pan_all20_mean']=mean20[:len(df)]
summary['primary_four_models_min']=npfold[:,:len(df)].min(axis=0);summary['primary_four_models_max']=npfold[:,:len(df)].max(axis=0)
summary['exploratory_class_at_0p5']=summary['NetTCR22_pan_known_overlap_excluded_4_mean']>=0.5
pd.DataFrame(raw).to_csv(root/'predictions_per_model.csv',index=False)
summary.to_csv(root/'nettcr22_pan_predictions.csv',index=False)
cvout=cvsel.copy();cvout['recomputed_t0_four_model_mean']=mean4[len(df):];cvout['abs_difference']=np.abs(cvout.prediction-cvout.recomputed_t0_four_model_mean)
cvout.to_csv(root/'CV_reproduction_validation.csv',index=False)
maxdiff=float(cvout.abs_difference.max());assert maxdiff<2e-6,maxdiff
# Confirm exact single-model predictions against the unmodified official CLI.
validation_out=root/'official_cli_validation';(validation_out/'checkpoint').mkdir(parents=True,exist_ok=True)
import shutil
shutil.copy2(model_dir/'t.0.v.1.tflite',validation_out/'checkpoint/t.0.v.1.tflite')
command=[sys.executable,str(root/'upstream/src/predict.py'),'--test_data',str(root/'input_nettcr.csv'),'--outdir',str(validation_out),'--model_name','t.0.v.1','--model_type','pan']
run=subprocess.run(command,text=True,capture_output=True,encoding='utf-8',errors='replace')
(validation_out/'run.log').write_text(run.stdout+'\n'+run.stderr,encoding='utf-8');assert run.returncode==0,run.stderr
cli=pd.read_csv(validation_out/'t.0.v.1_prediction.csv').prediction.to_numpy()
wrapper=pd.DataFrame(raw).query("model == 't.0.v.1'").prediction.to_numpy()
cli_diff=float(np.max(np.abs(cli-wrapper)));assert cli_diff<1e-7,cli_diff
save(root/'inference_audit.json',{
 'input_file':str(input_path),'input_sha256':sha(input_path),'input_rows':len(df),'official_repository':'https://github.com/mnielLab/NetTCR-2.2','official_commit':json.loads((root/'download_audit.json').read_text())['commit'],
 'CDR_definition':'ANARCII IMGT CDR1 27-38, CDR2 56-65, CDR3 105-117, excluding C104/F118; matches official training schema.',
 'model_mode':'Uniform pan-specific NetTCR2.2 for every peptide; no per-peptide model switching and no TCRbase scaling.',
 'primary_selection_before_inference':'t.0.v.1/2/3/4, because all detected exact CDR3 alpha/beta pairs in the final training list lie in partition0; exclude known training and validation overlap.',
 'primary_selection_scope':'Only verified exact CDR3-pair overlap exclusion. Does not prove all near-homologous TCR families excluded.',
 'model_members_are_not_biological_repeats':True,
 'threshold':'No universal author binary cutoff used. A score>=0.5 column is exploratory and explicitly labelled; no threshold tuned on this panel.',
 'score_semantics':'Binding specificity sigmoid score, not KD or calibrated affinity; rank differences are descriptive.',
 'training_overlap_rows':len(matched),'training_overlap_partitions':sorted(int(x) for x in matched.partition.unique()),
 'reference_CV_rows':len(cvout),'max_reference_CV_absolute_difference':maxdiff,'unmodified_official_CLI_max_difference':cli_diff,
 'models':models,
 'runtime':{'python':sys.version,'platform':platform.platform(),'tensorflow':tf.__version__,'numpy':np.__version__,'pandas':pd.__version__},
 'limits':['ELAGIGILTV is a training peptide; Ala substitutions are unseen peptide sequences.','Uniform pan model is a declared model component; official server would use epitope-specific/TCRbase models for some known peptides.','12 related conditions are a case-specific panel, not a general benchmark.','Score 0.5 is not asserted to be an official decision threshold.']
})
(root/'requirements_runtime_frozen.txt').write_bytes(subprocess.check_output([sys.executable,'-m','pip','freeze']))
print(summary[['peptide','NetTCR22_pan_known_overlap_excluded_4_mean','NetTCR22_pan_all20_mean']].to_string(index=False))
print(json.dumps({'CV_reproduction_rows':len(cvout),'CV_max_error':maxdiff,'official_CLI_max_error':cli_diff},indent=2))
