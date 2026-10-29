from pathlib import Path
import subprocess,sys,os,json,hashlib
r=Path(__file__).resolve().parent
up=r/'upstream'
env=os.environ.copy();env['PYTHONPATH']=str(up);env['TF_CPP_MIN_LOG_LEVEL']='2';env['OMP_NUM_THREADS']='2';env['TF_NUM_INTRAOP_THREADS']='2';env['TF_NUM_INTEROP_THREADS']='2'
runs=[]
for tag in ['downsample','full']:
 cmd=[sys.executable,'-X','utf8',str(up/'src/scripts/predict/predict.py'),'--model',str(r/'weights'/(tag+'.h5')),'--input',str(r/'input_imrex.csv'),'--output',str(r/('predictions_'+tag+'.csv'))]
 p=subprocess.run(cmd,cwd=up,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace')
 (r/('run_'+tag+'.log')).write_text(p.stdout,encoding='utf-8')
 runs.append({'tag':tag,'command':cmd,'exit_code':p.returncode,'log':str(r/('run_'+tag+'.log'))})
 print(tag,p.returncode,p.stdout[-3500:])
 if p.returncode:break
(r/'inference_runs.json').write_text(json.dumps(runs,indent=2),encoding='utf-8')
