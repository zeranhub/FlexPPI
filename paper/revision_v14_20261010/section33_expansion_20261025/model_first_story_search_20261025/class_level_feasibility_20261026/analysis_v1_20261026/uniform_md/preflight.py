import os,sys,json,importlib.util,hashlib,collections
from pathlib import Path
os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
ROOT=Path(r'E:\Project\FlexPPI');OUT=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\class_level_feasibility_20261026\analysis_v1_20261026\uniform_md')
spec=importlib.util.spec_from_file_location('sens',ROOT/'experiments/_sources/label_sensitivity/run_label_sensitivity.py');sens=importlib.util.module_from_spec(spec);spec.loader.exec_module(sens)
import pandas as pd,numpy as np
cand=pd.read_csv(OUT.parent.parent/'inventory_agent/short_peptide_5to40_candidates.csv');cand=cand[cand.two_chain_unambiguous].to_dict('records')
records={r['sample_id']:r for r in sens.jsonl(ROOT/'experiments/core_corrected_pipeline/training_payload/splits/all_ready.jsonl')}
orig={r['sample_id']:r for r in sens.jsonl(ROOT/'data/splits/state_rmsf_hchi_v2_performance_first/all_ready.jsonl')}
rows=[]; summary=[]
for c in cand:
 sid=c['sample_id'];r=records[sid];o=orig[sid];dst=OUT/'preflight'/f'{sid}.json';dst.parent.mkdir(exist_ok=True)
 if dst.exists(): row=json.loads(dst.read_text(encoding='utf-8'))
 else:
  row=dict(sample_id=sid,pdb_id=c['pdb_id'],family_id=c['family_id'],dataset=r['dataset'],short_chain=c['short_chain'],short_length=c['short_length'],receptor_chains=c['receptor_chains'],fold=r['metadata']['fold'],temperature_declared_K=r['temperature'],temperature_provenance='frozen dataset metadata, not newly recovered thermostat file',target_chains=r['target_chains'],structure=str(ROOT/'experiments/core_corrected_pipeline/training_payload/dataset/structures'/f'{sid}.pdb'),old_label=str(ROOT/'experiments/core_corrected_pipeline/training_payload/dataset/labels'/f'{sid}.npz'),trajectories=[])
  if r['dataset']=='dynarepo': pairs=[(sens.legacy.resolve(o['structure_path']),sens.legacy.resolve(t)) for t in o['trajectory_paths']]
  elif r['dataset']=='itrace': pairs=[(Path(f'E:/data/ITrace/trajectories/{c["pdb_id"].lower()}_run{k}_topology.pdb'),Path(f'E:/data/ITrace/trajectories/{c["pdb_id"].lower()}_run{k}_md_processed.xtc')) for k in [1,2,3]]
  elif r['dataset']=='atlas':
   f=Path('E:/data/protein_md_datasets/ATLAS_extracted')/sid;pairs=[(f/f'{sid}.pdb',f/f'{sid}_R{k}.xtc') for k in [1,2,3]]
  for k,(top,xtc) in enumerate(pairs,1):
   z=dict(replica_index=k,topology=str(top),trajectory=str(xtc),topology_exists=top.exists(),trajectory_exists=xtc.exists())
   try:
    z.update(sens.scan_xtc_headers(xtc));z['trajectory_bytes']=xtc.stat().st_size;z['trajectory_mtime_ns']=xtc.stat().st_mtime_ns;z['topology_sha256']=sens.sha256(top)
    z['uniform_0_100ns_eligible']=z['first_ps']==0 and z['last_ps']>=100000 and not z['backwards_time_indices']
   except Exception as ex:z.update(error=repr(ex),uniform_0_100ns_eligible=False)
   row['trajectories'].append(z)
  row['status']='preflight_complete';dst.write_text(json.dumps(row,ensure_ascii=False,indent=2),encoding='utf-8')
 rows.append(row)
 for t in row['trajectories']: summary.append({k:row.get(k) for k in ['sample_id','pdb_id','family_id','dataset','fold','temperature_declared_K']}|{k:t.get(k) for k in ['replica_index','topology','trajectory','topology_exists','trajectory_exists','frames','atom_count','first_ps','last_ps','median_step_ps','trajectory_bytes','uniform_0_100ns_eligible','error']}|dict(duplicate_times=len(t.get('duplicate_time_indices',[])),backwards_times=len(t.get('backwards_time_indices',[]))))
 print('preflight',len(rows),len(cand),sid,r['dataset'],[(t.get('frames'),t.get('last_ps'),t.get('trajectory_bytes')) for t in row['trajectories']],flush=True)
pd.DataFrame(summary).to_csv(OUT/'trajectory_inventory.csv',index=False)
(OUT/'preflight_summary.json').write_text(json.dumps(dict(records=len(rows),PDBs=len({r['pdb_id'] for r in rows}),datasets=dict(collections.Counter(r['dataset'] for r in rows)),trajectories=len(summary),eligible_trajectories=sum(t['uniform_0_100ns_eligible'] for r in rows for t in r['trajectories']),total_GB=sum(t.get('trajectory_bytes',0) for r in rows for t in r['trajectories'])/1e9,records_with_all_trajectories_eligible=sum(all(t['uniform_0_100ns_eligible'] for t in r['trajectories']) for r in rows)),indent=2),encoding='utf-8')
print((OUT/'preflight_summary.json').read_text(),flush=True)
