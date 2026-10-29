"""Compute only complete-structure static caches; never read NPZ or SPR labels."""
from pathlib import Path
import argparse,json,datetime,time
from concurrent.futures import ProcessPoolExecutor,as_completed
import extract_features as f

def initialize(paths):f.dependencies(paths)
def compute_one(payload):
    r,s=payload;iid=f.input_id(r);roles=f.role_map(r);sp=f.static_path(s)
    fixed,_=f.read_structure(f.input_path(r));st,extra=f.read_structure(sp)
    q=f.structure_qc(st,roles,r['peptide_sequence']);bb=f.backbone_comparison(fixed,st)
    if not (q['peptide_sequence_matches'] and q['chain_roles_present'] and q['all_canonical_heavy_atoms_complete'] and not q['unsupported_atom_elements'] and not extra and bb['passed']):
        return {'input_id':iid,'status':'QC_failed','missing_count':len(q['canonical_missing']),'missing_first3':q['canonical_missing'][:3],'backbone':bb}
    settings={'static_pdb_sha256':f.digest(sp),'radii_A':f.RADII,'SASA_probe_A':1.4,'SASA_n_points':960,'roles':roles,'contact_cutoff_A':4.5,'biopython_version':f.Bio.__version__}
    out=f.HERE/'static_cache'/(iid+'.json')
    if out.exists() and f.read_json(out).get('settings')==settings:return {'input_id':iid,'status':'cached'}
    start=time.perf_counter();values=f.calculate_static(st,roles)
    f.dump(out,{'input_id':iid,'settings':settings,'values':values})
    return {'input_id':iid,'status':'passed','seconds':time.perf_counter()-start}

def main():
    p=argparse.ArgumentParser();p.add_argument('--runtime-path',action='append',default=[])
    p.add_argument('--static-manifest',default=str(f.ROOT/'structures/static_complete_inputs.json'))
    p.add_argument('--fixed-manifest',default=str(f.ROOT/'structures/fixed_backbone_manifest.json'))
    p.add_argument('--workers',type=int,default=2)
    a=p.parse_args();fr=f.records(f.read_json(a.fixed_manifest));sr={f.input_id(x):x for x in f.records(f.read_json(a.static_manifest))}
    f.HERE.joinpath('static_cache').mkdir(exist_ok=True)
    start=time.perf_counter();results=[]
    with ProcessPoolExecutor(max_workers=a.workers,initializer=initialize,initargs=(a.runtime_path,)) as pool:
        jobs=[pool.submit(compute_one,(r,sr[f.input_id(r)])) for r in fr]
        for job in as_completed(jobs):
            rec=job.result();results.append(rec)
            print(json.dumps({'done':len(results),'total':len(fr),**rec}),flush=True)
    report={'status':'passed' if all(x['status'] in {'passed','cached'} for x in results) else 'complete_with_QC_failures','created_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'static_manifest_sha256':f.digest(a.static_manifest),'fixed_manifest_sha256':f.digest(a.fixed_manifest),'n_inputs':len(results),'results':results,'elapsed_seconds':time.perf_counter()-start,'workers':a.workers,'model_predictions_or_SPR_read':False,'feature_manifest_generated':False}
    f.dump(f.HERE/'static_precomputation_audit.json',report)
    print(json.dumps({k:report[k] for k in ['status','n_inputs','elapsed_seconds','model_predictions_or_SPR_read']}),flush=True)
if __name__=='__main__':main()
