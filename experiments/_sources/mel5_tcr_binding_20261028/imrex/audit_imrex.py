from pathlib import Path
import os,sys,json,hashlib,subprocess
os.environ['TF_CPP_MIN_LOG_LEVEL']='2';os.environ['OMP_NUM_THREADS']='2';os.environ['TF_NUM_INTRAOP_THREADS']='2';os.environ['TF_NUM_INTEROP_THREADS']='2'
import numpy as np,pandas as pd,tensorflow as tf,h5py
r=Path(__file__).resolve().parent;up=r/'upstream';sys.path.insert(0,str(up))
from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features,parse_operator
from src.processing.data_stream import DataStream
from src.processing.image_generator import ImageGenerator
from src.processing.image_padding import ImagePadding

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
panel=pd.read_csv(r.parent/'binding_panel_inputs.csv')
inp=pd.read_csv(r/'input_imrex.csv',sep=';')
assert inp.case_id.tolist()==panel.case_id.tolist()
assert inp.cdr3.tolist()==panel.CDR3_beta.tolist() and inp['antigen.epitope'].tolist()==panel.peptide.tolist()
assert inp.cdr3.str.len().between(10,20).all() and inp['antigen.epitope'].str.len().between(8,11).all()
fb=CombinedPeptideFeatureBuilder(parse_features('hydrophob,isoelectric,mass,hydrophil'),parse_operator('absdiff'))
gen=ImageGenerator(DataStream(inp[['cdr3','antigen.epitope']].to_numpy()),fb,has_label=False)
x=np.array(list(ImagePadding(gen,20,11,pad_value=0,has_label=False)))
assert x.shape==(12,20,11,4) and np.isfinite(x).all()
np.savez_compressed(r/'official_encoded_inputs.npz',x=x)
audit={'model':'ImRex','repository':'https://github.com/pmoris/ImRex','commit':'d3975b744c030c4bc67c942c82a32b1c83857c4e','paper_doi':'10.1093/bib/bbaa318','default_checkpoint':'2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001','sensitivity_checkpoint':'2020-07-30_11-30-27_trbmhci-shuffle-padded-b32-lre4-reg001','primary_selection':'Default checkpoint named in official README, not selected using this panel results. Both official full-data checkpoints retained.','input':'CDR3 beta (with anchors C104/F118) plus epitope only; no alpha chain, CDR1/2, HLA or structure.','input_limits':{'CDR3_beta':[10,20],'epitope':[8,11]},'features':'Official CombinedPeptideFeatureBuilder: hydrophob,isoelectric,mass,hydrophil; absdiff; centered zero padding to 20x11x4.','input_panel_sha256':sha(r.parent/'binding_panel_inputs.csv'),'official_input_sha256':sha(r/'input_imrex.csv'),'encoded_shape':list(x.shape),'inference_code_sha256':sha(up/'src/scripts/predict/predict.py'),'inference_code_changes':'None. Original official CLI used for all 12 rows, both checkpoints.','classification_reference':{'cutoff':0.5,'rule':'score > 0.5','basis':'Fixed generic sigmoid reference consistent with original trainer BinaryAccuracy()/Precision()/Recall() defaults; official prediction CLI itself exports scores only.','calibrated_for_assay':False},'score_interpretation':'Recognition score, not KD, physical binding affinity, or a calibrated probability at SPR detection threshold.','runtime':{'python':sys.version,'tensorflow':tf.__version__,'numpy':np.__version__,'pandas':pd.__version__,'keras':'2.15.0','biopython':'1.79','pyteomics':'4.7.5','platform':'Windows CPU'},'runtime_note':'Original repository pinned TF2.1/Python-era packages; compatible TF2.15 used. pyteomics4.2 had pandas metadata API incompatibility before inference; upgraded to4.7.5 without source/model changes. Deployment failure saved separately and never counted as prediction failure.','checkpoints':{}}
allrows=[]
for tag in ['downsample','full']:
 train=pd.read_csv(r/'training'/('train_'+tag+'.csv'),sep=';')
 pred=pd.read_csv(r/('predictions_'+tag+'.csv'))
 assert pred.case_id.tolist()==panel.case_id.tolist()
 model=tf.keras.models.load_model(r/'weights'/(tag+'.h5'))
 api=model.predict(tf.data.Dataset.from_tensor_slices(x).batch(128),verbose=0).reshape(-1)
 maxdiff=float(np.max(np.abs(api-pred.prediction_score.to_numpy())))
 assert maxdiff<1e-7
 out=panel.copy();out['model']='ImRex';out['checkpoint']=tag;out['binding_score']=[float(v) for v in api];out['predicted_binding_at_0p5']=(api>0.5).astype(int);out['official_cli_score']=pred.prediction_score
 covers=[]
 for i,row in out.iterrows():
  exact=train[(train.cdr3==row.CDR3_beta)&(train['antigen.epitope']==row.peptide)]
  out.loc[i,'exact_pair_train_rows']=len(exact)
  out.loc[i,'CDR3_train_rows']=int((train.cdr3==row.CDR3_beta).sum())
  out.loc[i,'epitope_train_rows']=int((train['antigen.epitope']==row.peptide).sum())
  covers.append({'case_id':row.case_id,'exact_rows':exact.to_dict('records'),'cdr3_train_rows':int((train.cdr3==row.CDR3_beta).sum()),'epitope_train_rows':int((train['antigen.epitope']==row.peptide).sum())})
 out.to_csv(r/('audited_scores_'+tag+'.csv'),index=False,float_format='%.17g')
 allrows+=out.to_dict('records')
 truth=out.experimental_binding.to_numpy().astype(int);labels=out.predicted_binding_at_0p5.to_numpy().astype(int)
 cm={'TP':int(((truth==1)&(labels==1)).sum()),'FN':int(((truth==1)&(labels==0)).sum()),'TN':int(((truth==0)&(labels==0)).sum()),'FP':int(((truth==0)&(labels==1)).sum())}
 state=out.set_index('case_id').binding_score.to_dict();baseline=state['4jff']
 with h5py.File(r/'weights'/(tag+'.h5'),'r') as f:config=json.loads(f.attrs['model_config'])
 audit['checkpoints'][tag]={'weight_sha256':sha(r/'weights'/(tag+'.h5')),'weight_bytes':(r/'weights'/(tag+'.h5')).stat().st_size,'train_csv_sha256':sha(r/'training'/('train_'+tag+'.csv')),'training_rows':len(train),'training_label_counts':{str(k):int(v) for k,v in train.y.value_counts().items()},'official_cli_api_max_absolute_difference':maxdiff,'confusion_at_reference_cutoff':cm,'all_case_scores':state,'G4A_minus_mature_ELA':state['4jfd']-baseline,'I7A_minus_mature_ELA':state['4jfe']-baseline,'training_overlap':covers,'serialized_input_shape':model.input_shape,'serialized_final_activation':config['config']['layers'][-1]['config'].get('activation')}
pd.concat([pd.read_csv(r/'audited_scores_downsample.csv'),pd.read_csv(r/'audited_scores_full.csv')]).to_csv(r/'all_imrex_scores.csv',index=False,float_format='%.17g')
(r/'inference_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
(r/'requirements_frozen.txt').write_bytes(subprocess.check_output([sys.executable,'-m','pip','freeze']))
print(json.dumps({k:{'scores':v['all_case_scores'],'confusion':v['confusion_at_reference_cutoff'],'api_cli_delta':v['official_cli_api_max_absolute_difference']} for k,v in audit['checkpoints'].items()},ensure_ascii=False,indent=2))
