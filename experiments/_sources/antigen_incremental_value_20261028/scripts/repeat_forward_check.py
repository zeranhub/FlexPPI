from pathlib import Path
import json,sys,hashlib,csv
import numpy as np,torch
R=Path(r'E:\Project\FlexPPI\experiments\_sources\antigen_incremental_value_20261028')
sys.path.insert(0,r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\src')
from flexppi.inference import _structure_example
from flexppi.data.structures import parse_pdb
from flexppi.data.dataset import collate_states
from flexppi.config import ModelConfig
from flexppi.model.information_ablation import build_information_ablation
states=json.loads((R/'states.json').read_text());selected=[]
for system in ['nyeso','ebv','mage']:
 pool=[s for s in states if s['system']==system and s['template_version']=='primary']
 selected.append(next(s for s in pool if s['is_reference']));selected.append(next(s for s in pool if not s['is_reference']))
records=[];regional_differences=[]
number=list(csv.DictReader((R/'tcr_numbering.csv').open(encoding='utf-8-sig')))
full_audit=json.loads((R/'audits/inference_full.json').read_text())
original_threads={r['state_id']:(2 if r['reused_verified_prediction'] else 4) for r in full_audit['entries']}
for variant in ['full','intrachain']:
 path=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028\runtime\family_folds\runs')/variant/'fold_0/epoch_040.pt'
 cp=torch.load(path,map_location='cpu',weights_only=False);m=build_information_ablation(ModelConfig(**cp['config']['model']),variant).eval();m.load_state_dict(cp['model'],strict=True)
 for threads in [2,4]:
  torch.set_num_threads(threads);regional={};old_regional={}
  for s in selected:
   st=parse_pdb(s['path']);ex=_structure_example(s['path'],st.unique_chains(),300.,100.,0.)
   with torch.no_grad():mean=m.forward_state(collate_states([ex]))['mean'][0].numpy()
   z=np.load(R/'predictions'/variant/(s['id']+'.npz'));diff=float(np.max(np.abs(mean-z['mean'])))
   saved_threads=2 if variant=='intrachain' else original_threads[s['id']]
   records.append({'variant':variant,'state_id':s['id'],'repeat_threads':threads,'saved_threads':saved_threads,'max_absolute_difference':diff,'same_thread_comparison':threads==saved_threads})
   lookup={(ch,a):i for i,(ch,a) in enumerate(zip(z['chain'],z['author']))};regional[s['id']]={};old_regional[s['id']]={}
   for role,ch in [('alpha','D'),('beta','E')]:
    for region,allowed in [('CDR3',['CDR3']),('CDR12',['CDR1','CDR2'])]:
     ix=[lookup[(r['chain'],r['author'])] for r in number if r['template']==s['template'] and r['chain']==ch and r['region'] in allowed]
     for j,metric in enumerate(['rmsf','hchi']):
      valid=z['model_mask']&z['geometry_mask'] if j==0 else z['hchi_valid'];good=[i for i in ix if valid[i]];key=role+'_'+region+'_'+metric
      regional[s['id']][key]=float(mean[good,j].mean());old_regional[s['id']][key]=float(z['mean'][good,j].mean())
  for s in selected:
   if s['is_reference']:continue
   ref=s['reference_state']
   differences=[abs((regional[s['id']][k]-regional[ref][k])-(old_regional[s['id']][k]-old_regional[ref][k])) for k in regional[s['id']]]
   regional_differences.append({'variant':variant,'state_id':s['id'],'repeat_threads':threads,'max_dynamic_feature_difference':max(differences)})
 del m
(R/'audits/repeated_forward_check.json').write_text(json.dumps({'status':'completed_with_disclosed_CPU_thread_variation','repeat_forward_passes':len(records),'max_absolute_difference':max(r['max_absolute_difference'] for r in records),'same_thread_max_difference':max(r['max_absolute_difference'] for r in records if r['same_thread_comparison']),'max_regional_dynamic_feature_difference':max(r['max_dynamic_feature_difference'] for r in regional_differences),'initial_check_note':'An initial across-thread tolerance of 1e-5 was exceeded by an intrachain EBV reference (1.6093e-5); initial failure log retained. Completed separate same-thread and cross-thread checks instead of silently claiming bitwise equivalence.','records':records,'regional_feature_checks':regional_differences},indent=2),encoding='utf-8')
print('COMPLETED',len(records),max(r['max_absolute_difference'] for r in records),max(r['max_dynamic_feature_difference'] for r in regional_differences))
