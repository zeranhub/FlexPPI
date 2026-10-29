"""Version 1: score-blind 0-100 ns labels from existing local MD, source files read-only.
Derived from the original corrected chain-fit and joint-three-bin chi estimators.
Every local trajectory is retained. Continuous 100-ns blocks are sensitivity windows,
not independent replicas. Selection is determined by saved timestamps, never prediction.
"""
from __future__ import annotations
import os
for variable in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']: os.environ[variable]='1'
import sys,json,hashlib,importlib.util,time,tempfile,shutil,argparse,traceback,warnings
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
ROOT=Path(r'E:\Project\FlexPPI');OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('sens',ROOT/'experiments/_sources/label_sensitivity/run_label_sensitivity.py');sens=importlib.util.module_from_spec(spec);spec.loader.exec_module(sens)
import numpy as np,pandas as pd

def write(path,obj):
 path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')

def freeze():
 rows=[json.loads(p.read_text(encoding='utf-8')) for p in sorted((OUT/'preflight').glob('*.json'))]
 assert len(rows)==56
 for r in rows:
  for t in r['trajectories']:
   t['windows']=[]
   if not t['uniform_0_100ns_eligible']:continue
   times=np.asarray(t['times_ps']); count=int(np.floor(t['last_ps']/100000.0))
   for k in range(count):
    grid=np.arange(k*100000,(k+1)*100000+.1,200.0)
    ids=np.searchsorted(times,grid);ids=np.clip(ids,0,len(times)-1);prior=np.maximum(ids-1,0);use=abs(times[prior]-grid)<abs(times[ids]-grid);ids[use]=prior[use]
    error=abs(times[ids]-grid)
    if error.max()>.05 or len(np.unique(ids))!=len(grid):raise ValueError(f'{r["sample_id"]}: unsupported exact grid')
    t['windows'].append(dict(window_index=k,start_ps=float(grid[0]),end_ps=float(grid[-1]),frame_count=len(grid),source_indices=ids.tolist(),max_grid_time_error_ps=float(error.max())))
  r['eligible']=bool(r['trajectories']) and all(t['uniform_0_100ns_eligible'] for t in r['trajectories'])
  r['exclusion_reason']=None if r['eligible'] else 'Saved trajectory does not cover monotone 0-100 ns. Membership retained; uniform labels unavailable.'
 plan=dict(protocol='short_peptide_common_0_100ns_200ps_chain_Kabsch_joint3chi_v1',members=56,eligible=sum(r['eligible'] for r in rows),records=rows,main_window=dict(start_ps=0,end_ps=100000,step_ps=200,frames=501,endpoints_inclusive=True),replicas='all available local trajectories; currently one per record',blocks='All complete consecutive 100-ns blocks; shared endpoint allowed. Dependent sensitivity windows, not independent replicas.',geometry='No PBC or coordinate reconstruction introduced; covalent-geometry diagnostics retained, not silently filtered.',rmsf='Per-chain all-C-alpha Kabsch to each window first frame; RMS fluctuation around aligned mean; Angstrom.',hchi='All defined chi atoms required; joint 3^q state from fixed [-180,-60),[-60,60),[60,180] bins, Shannon natural-log entropy in nats. PRO has two chi angles and is valid when complete. A/G undefined.',mask='Finite metric, canonical backbone frame mask, target-chain mask; compare original mask per position.',source_code_sha256={str(p):sens.sha256(p) for p in [Path(__file__),Path(sens.__file__),Path(sens.legacy.__file__),ROOT/'experiments/core_corrected_pipeline/training_payload/src/flexppi/constants.py']})
 path=OUT/'uniform_protocol_frozen.json'
 if path.exists():
  old=json.loads(path.read_text(encoding='utf-8'));assert old==plan,'Frozen protocol changed; create an explicit new version.'
 else:write(path,plan)
 return plan

def staged_sparse_xtc(item,dst):
 source=Path(item['trajectory']);before=source.stat();assert before.st_size==item['trajectory_bytes'] and before.st_mtime_ns==item['trajectory_mtime_ns']
 indices=sorted({i for w in item['windows'] for i in w['source_indices']});chosen=set(indices);offsets=item['frame_offsets'];ends=offsets[1:]+[item['trajectory_bytes']]
 digest=hashlib.sha256();subset=hashlib.sha256();copied=0
 with source.open('rb') as src,dst.open('wb') as target:
  for i,(start,end) in enumerate(zip(offsets,ends)):
   assert src.tell()==start;data=src.read(end-start);assert len(data)==end-start;digest.update(data)
   if i in chosen:target.write(data);subset.update(data);copied+=len(data)
 after=source.stat();assert before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns
 return indices,dict(source_sha256=digest.hexdigest(),selected_xtc_sha256=subset.hexdigest(),selected_bytes=copied,source_bytes=before.st_size,source_size_mtime_unchanged=True)

def extract(row,item,folder):
 import MDAnalysis as mda
 top=OUT/'topologies'/f'{row["sample_id"]}_rep{item["replica_index"]}.pdb';changes=sens.legacy.canonical_pdb(item['topology'],top)
 assert sens.sha256(item['topology'])==item['topology_sha256']
 xtc=folder/'selected.xtc';selected,stageaudit=staged_sparse_xtc(item,xtc);order={src:i for i,src in enumerate(selected)}
 reference=sens.parse_pdb(row['structure']);replica=sens.parse_pdb(str(top));mapping,mapstats=sens.legacy.replica_map(reference,replica)
 univ=mda.Universe(str(top),str(xtc))
 try:
  assert len(univ.trajectory)==len(selected)
  residues=[r for r in univ.residues if r.resname in sens.THREE_TO_ONE and 'CA' in r.atoms.names]
  keys=[(str(r.atoms.chainIDs[0]).strip() or '_',str(r.resid)+str(r.icode).strip(),str(r.resname)) for r in residues]
  assert keys==list(zip(replica.chain_ids,replica.author_residue_ids,replica.residue_names))
  ca=[];quads=[];chis={};residue_audit=[]
  for i,r in enumerate(residues):
   lookup={a.name:a.index for a in r.atoms};assert len(lookup)==len(r.atoms),'Duplicate atom names in residue'
   ca.append(lookup['CA']);definitions=sens.CHI_ATOMS.get(r.resname,[]);missing=sorted({a for quad in definitions for a in quad if a not in lookup})
   if definitions and not missing:
    chis[i]=list(range(len(quads),len(quads)+len(definitions)));quads += [[lookup[a] for a in q] for q in definitions]
   residue_audit.append(dict(replica_residue_index=i,replica_chain=keys[i][0],replica_author_residue_id=keys[i][1],replica_resname=r.resname,amino_acid=sens.THREE_TO_ONE[r.resname],chi_count_defined=len(definitions),chi_valid_topology=bool(definitions and not missing),chi_missing_atoms='|'.join(missing),chi_reason='valid' if definitions and not missing else ('undefined_no_chi' if not definitions else 'missing_chi_atom')))
  qarr=np.asarray(quads,dtype=int);coords=[];states=[];badchi=0;maxbond=0.;max_step=0.;prevca=None
  for n,frame in enumerate(univ.trajectory):
   expected=item['times_ps'][selected[n]];assert abs(frame.time-expected)<=.05
   capos=frame.positions[ca].copy();assert np.isfinite(capos).all()
   if prevca is not None:max_step=max(max_step,float(np.linalg.norm(capos-prevca,axis=-1).max()))
   prevca=capos;coords.append(capos)
   if len(quads):
    points=frame.positions[qarr];assert np.isfinite(points).all();bonds=np.linalg.norm(np.diff(points,axis=1),axis=-1);badchi+=int((bonds>5).sum());maxbond=max(maxbond,float(bonds.max()))
    angles=sens._dihedral_degrees(points);assert np.isfinite(angles).all();states.append(np.clip(np.floor((angles+180)/120).astype(np.int8),0,2))
  coords=np.asarray(coords);states=np.asarray(states);windows=[]
  for w in item['windows']:
   take=np.asarray([order[i] for i in w['source_indices']]);sample=np.full((replica.num_residues,len(sens.METRIC_NAMES)),np.nan);sample[:,0]=sens.legacy.chain_fit(coords[take],replica.chain_ids)
   for i,columns in chis.items():sample[i,1]=sens._entropy(sum(states[take,col].astype(np.int32)*3**p for p,col in enumerate(columns)))
   mapped=np.full((reference.num_residues,len(sens.METRIC_NAMES)),np.nan)
   for target,source in mapping:mapped[target]=sample[source]
   windows.append(mapped)
  rows=[];sourcemap={target:source for target,source in mapping}
  for i in range(reference.num_residues):
   if i in sourcemap:extra=residue_audit[sourcemap[i]]
   else:extra=dict(chi_reason='unmapped_canonical_residue')
   rows.append(dict(sample_id=row['sample_id'],replica_index=item['replica_index'],residue_index=i,chain=reference.chain_ids[i],author_residue_id=reference.author_residue_ids[i],amino_acid=reference.sequence[i],canonical_resname=reference.residue_names[i],frame_valid=bool(reference.frame_mask[i]))|extra)
  return windows,dict(replica_index=item['replica_index'],source_frames=item['frames'],selected_frames=len(selected),windows=len(windows),canonical_name_changes=changes,mapping=mapstats,residue_coverage=len(mapping)/reference.num_residues,chi_residues=len(chis),chi_bond_instances_gt5A=badchi,chi_max_bond_A=maxbond,max_CA_step_A=max_step,topology_sha256=item['topology_sha256'],canonical_topology_sha256=sens.sha256(top),**stageaudit),rows
 finally:univ.trajectory.close()

def process(args):
 row,planhash=args;sid=row['sample_id'];auditpath=OUT/'record_audit'/f'{sid}.json';started=time.monotonic()
 if auditpath.exists():
  cached=json.loads(auditpath.read_text(encoding='utf-8'))
  if cached.get('status')=='complete' and cached.get('plan_sha256')==planhash:return cached
 report=dict(sample_id=sid,pdb_id=row['pdb_id'],family_id=row['family_id'],dataset=row['dataset'],fold=row['fold'],plan_sha256=planhash)
 if not row['eligible']:
  report.update(status='unavailable_uniform_window',reason=row['exclusion_reason']);write(auditpath,report);return report
 try:
  reference=sens.parse_pdb(row['structure']);assert sens.sha256(row['structure'])
  olds=np.load(row['old_label'],allow_pickle=False);oldv=olds['values'];oldm=olds['mask'];assert len(oldv)==reference.num_residues
  replica_windows=[];audits=[];mapping_rows=[]
  staging=OUT/'_staging';staging.mkdir(exist_ok=True)
  for item in row['trajectories']:
   with tempfile.TemporaryDirectory(prefix=f'{sid}_',dir=staging) as td:
    samples,audit,maps=extract(row,item,Path(td));replica_windows.append(samples);audits.append(audit);mapping_rows+=maps
  maskbase=reference.frame_mask[:,None]&sens.target_mask(reference,row['target_chains'])[:,None]
  labels={};outrows=[];maskrows=[]
  variants=[('uniform_rep1',[replica_windows[0][0]]),('uniform_available_replicas',[windows[0] for windows in replica_windows])]
  max_windows=min(len(w) for w in replica_windows)
  variants += [(f'block_{k*100:03d}_{(k+1)*100:03d}ns',[windows[k] for windows in replica_windows]) for k in range(max_windows)]
  for name,samples in variants:
   values,variances,mask=sens._summarise_samples(samples,reference.num_residues);mask &= maskbase;path=OUT/'labels'/name/f'{sid}.npz';sens.save_label_file(str(path),values,variances,mask)
   labels[name]=dict(path=str(path),sha256=sens.sha256(path),valid_counts=mask[:,:2].sum(0).tolist(),replicas=len(samples),variance_interpretation='zero when one local replica; does not estimate sampling uncertainty')
   for i in range(reference.num_residues):
    outrows.append(dict(sample_id=sid,pdb_id=row['pdb_id'],family_id=row['family_id'],fold=row['fold'],dataset=row['dataset'],variant=name,residue_index=i,chain=reference.chain_ids[i],author_residue_id=reference.author_residue_ids[i],amino_acid=reference.sequence[i],rmsf=float(values[i,0]),chi_entropy=float(values[i,1]),valid_rmsf=bool(mask[i,0]),valid_chi_entropy=bool(mask[i,1])))
    if name=='uniform_rep1':
     for j,metric in enumerate(sens.METRIC_NAMES[:2]):
      if bool(mask[i,j])!=bool(oldm[i,j]):maskrows.append(dict(sample_id=sid,residue_index=i,chain=reference.chain_ids[i],author_residue_id=reference.author_residue_ids[i],amino_acid=reference.sequence[i],metric=metric,old_valid=bool(oldm[i,j]),new_valid=bool(mask[i,j]),reason='see per-residue chemical topology/mapping audit and canonical backbone/target mask'))
  dest=OUT/'per_record';dest.mkdir(exist_ok=True);pd.DataFrame(outrows).to_csv(dest/f'{sid}_labels.csv.gz',index=False,compression='gzip');pd.DataFrame(mapping_rows).to_csv(dest/f'{sid}_residue_mapping.csv',index=False)
  write(dest/f'{sid}_mask_changes.json',maskrows)
  report.update(status='complete',replicas=len(replica_windows),replica_audits=audits,labels=labels,old_mask_counts=oldm[:,:2].sum(0).tolist(),mask_differences=maskrows,structure_sha256=sens.sha256(row['structure']),old_label_sha256=sens.sha256(row['old_label']),elapsed_seconds=time.monotonic()-started)
 except Exception as ex:report.update(status='failed',error=repr(ex),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-started)
 write(auditpath,report);return report

def aggregate(plan):
 reports=[json.loads(p.read_text(encoding='utf-8')) for p in (OUT/'record_audit').glob('*.json')];frames=[pd.read_csv(p) for p in sorted((OUT/'per_record').glob('*_labels.csv.gz'))]
 if frames:pd.concat(frames,ignore_index=True).to_csv(OUT/'all_uniform_labels_per_residue.csv.gz',index=False,compression='gzip')
 maps=[pd.read_csv(p) for p in sorted((OUT/'per_record').glob('*_residue_mapping.csv'))]
 if maps:pd.concat(maps,ignore_index=True).to_csv(OUT/'residue_chemical_mapping_audit.csv.gz',index=False,compression='gzip')
 maskchanges=[m for r in reports for m in r.get('mask_differences',[])];write(OUT/'old_new_mask_changes.json',maskchanges)
 summary=dict(status='complete' if len(reports)==56 and all(r['status'] in ['complete','unavailable_uniform_window'] for r in reports) else 'partial',members=56,reports=len(reports),computed=sum(r['status']=='complete' for r in reports),unavailable=[r for r in reports if r['status']=='unavailable_uniform_window'],failed=[r for r in reports if r['status']=='failed'],all_one_local_replica=all(r.get('replicas',1)==1 for r in reports),mask_change_count=len(maskchanges),geometry_flags=[dict(sample_id=r['sample_id'],replica_index=a['replica_index'],chi_bond_instances_gt5A=a['chi_bond_instances_gt5A'],max_CA_step_A=a['max_CA_step_A']) for r in reports for a in r.get('replica_audits',[]) if a['chi_bond_instances_gt5A'] or a['residue_coverage']<1],plan_sha256=sens.sha256(OUT/'uniform_protocol_frozen.json'))
 write(OUT/'compute_summary.json',summary);print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--samples',nargs='*');parser.add_argument('--workers',type=int,default=2);parser.add_argument('--aggregate-only',action='store_true');args=parser.parse_args();assert 1<=args.workers<=2
 plan=freeze()
 if not args.aggregate_only:
  todo=[r for r in plan['records'] if not args.samples or r['sample_id'] in args.samples];ph=sens.sha256(OUT/'uniform_protocol_frozen.json')
  with ProcessPoolExecutor(max_workers=args.workers) as pool:
   futs=[pool.submit(process,(r,ph)) for r in todo]
   for n,f in enumerate(as_completed(futs),1):
    a=f.result();print('compute',n,len(todo),a['sample_id'],a['status'],round(a.get('elapsed_seconds',0),2),flush=True)
 aggregate(plan)
if __name__=='__main__':main()
