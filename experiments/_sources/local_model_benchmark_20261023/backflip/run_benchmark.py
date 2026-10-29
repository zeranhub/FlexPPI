"""Fresh CPU inference with official frozen weights, paired input contexts.

The optional torch_scatter guard only permits unused MD-label helper imports;
it fails loudly if invoked.  No numerical network operation is substituted.
"""
from pathlib import Path
import argparse
import csv
import faulthandler
import hashlib
import importlib.util
import json
import os
import shutil
import sys
import time
import traceback
import types

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dependencies(source, threads):
    os.environ['OMP_NUM_THREADS'] = str(threads)
    os.environ['MKL_NUM_THREADS'] = str(threads)
    os.environ['OPENBLAS_NUM_THREADS'] = str(threads)
    os.environ['WANDB_MODE'] = 'disabled'
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(source))
    for dependency in (ROOT/'dependencies').glob('*'):
        if dependency.is_dir():
            sys.path.append(str(dependency))
    import torch
    torch.set_num_threads(threads)
    torch.set_num_interop_threads(1)
    guarded = False
    if importlib.util.find_spec('torch_scatter') is None:
        module = types.ModuleType('torch_scatter')
        def unavailable(*args, **kwargs):
            raise RuntimeError('torch_scatter is unavailable; an MD-label helper was unexpectedly invoked during inference')
        module.scatter = unavailable
        module.scatter_add = unavailable
        sys.modules['torch_scatter'] = module
        guarded = True
    from backflip.deployment.inference_class import BackFlip
    from backflip.deployment import utils as upstream_utils
    frames = getattr(upstream_utils,'frames_from_pdb_tite',upstream_utils.frames_from_pdb)
    return torch, BackFlip, frames, guarded


def metrics(y, p):
    import numpy as np
    from scipy.stats import spearmanr, pearsonr
    valid = np.isfinite(y) & np.isfinite(p)
    y,p = y[valid],p[valid]
    out = dict(count=len(y), mae=float(np.mean(abs(y-p))),
        rmse=float(np.sqrt(np.mean((y-p)**2))), bias=float(np.mean(p-y)),
        spearman=float('nan'), pearson=float('nan'))
    if len(y)>=4 and np.ptp(y)>0 and np.ptp(p)>0:
        out['spearman'] = float(spearmanr(y,p).statistic)
        out['pearson'] = float(pearsonr(y,p).statistic)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source', required=True)
    ap.add_argument('--checkpoint', required=True)
    ap.add_argument('--tag', default='backflip-2.1')
    ap.add_argument('--manifest', default=str(ROOT/'manifest.json'))
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--context', choices=['whole_complex','independent_chain','both'], default='both')
    ap.add_argument('--output-subdir', default='results')
    ap.add_argument('--limit', type=int)
    args = ap.parse_args()
    assert 1 <= args.threads <= 4
    source = Path(args.source)
    checkpoint = Path(args.checkpoint)
    print('IMPORTING',args.tag,flush=True)
    faulthandler.dump_traceback_later(90, repeat=False)
    torch,BackFlip,frames,guarded = dependencies(source,args.threads)
    faulthandler.cancel_dump_traceback_later()
    import numpy as np
    rows = json.loads(Path(args.manifest).read_text(encoding='utf-8'))
    if args.limit:
        rows = rows[:args.limit]
    output = ROOT/args.output_subdir/args.tag
    output.mkdir(parents=True, exist_ok=True)
    audit = dict(tag=args.tag, device='cpu', threads=args.threads, batch_size=1,
        torch_version=torch.__version__, checkpoint_sha256=sha(checkpoint),
        checkpoint_config_sha256=sha(checkpoint.parent/'config.yaml'),
        source_path=str(source), source_commit=('887c6a73b216e3fe2ff9f551ae6e2c983c21b292' if args.tag.startswith('backflip-2') else '60b8c83511b965925bbe990bb2be578bf3529460'),
        runner_sha256=sha(__file__), input_manifest_sha256=sha(args.manifest),
        torch_scatter_unreachable_import_guard=guarded,
        guard_definition='Missing optional torch_scatter training/MD-label functions raise if called. Network forward and official frame conversion remain unchanged.',
        output_origin='Fresh inference on this local workstation; no reused predictions.',
        primary_reference='Existing per-chain Kabsch RMSF labels. Official globally aligned covariance target differs for multichain input; no claim of training target equivalence.',
        positional_encoding='Official default concatenated arange across all valid input residues; no chain identity supplied to model.',
        coverage='Only complete N/CA/C geometry predicted; original residue identity retained. Failures disclosed.',
        timings_scope='Parse canonical PDB, official frame conversion, official CPU forward, output conversion; excludes checkpoint import/load and input preparation.')
    (output/'run_metadata.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    begin = time.perf_counter()
    model = BackFlip(checkpoint, device='cpu', progress_bar=False, confidence_intervals=False)
    audit['model_load_seconds'] = time.perf_counter()-begin
    audit['parameters'] = sum(p.numel() for p in model.model.parameters())
    print('LOADED',args.tag,audit['parameters'],audit['model_load_seconds'],flush=True)
    per_record = []
    for count,row in enumerate(rows,1):
        for context in (['whole_complex','independent_chain'] if args.context=='both' else [args.context]):
            path = output/(row['sample_id']+'__'+context+'.npz')
            if path.exists():
                continue
            status = dict(sample_id=row['sample_id'],dataset=row['dataset'],role=row['role'],
                tag=args.tag, context=context, n_residues=row['n_residues'])
            started = time.perf_counter()
            try:
                source_output = output/(row['sample_id']+'__whole_complex.npz')
                if context=='independent_chain' and row['n_chains']==1 and source_output.exists():
                    shutil.copyfile(source_output,path)
                    saved = dict(np.load(path))
                    status.update(status='complete',seconds=0.,output_sha256=sha(path),
                        predicted=int(np.isfinite(saved['prediction']).sum()),
                        reused_fresh_context='whole_complex (single-chain contexts are identical)',
                        metrics=metrics(saved['target'][saved['valid']],saved['prediction'][saved['valid']]))
                    with (output/'status.jsonl').open('a',encoding='utf-8') as f:
                        f.write(json.dumps(status)+'\n')
                    print(count,len(rows),row['sample_id'],context,'shared_current_single_chain_forward',flush=True)
                    continue
                pred = np.full(row['n_residues'],np.nan,np.float64)
                visited = np.zeros(row['n_residues'],np.int32)
                chunks = [x for x in row['contexts'] if (x['context']=='whole_complex') == (context=='whole_complex')]
                chunk_times = []
                for chunk in chunks:
                    assert sha(chunk['pdb_path']) == chunk['sha256']
                    tick = time.perf_counter()
                    frame_result = frames(Path(chunk['pdb_path']))
                    features = frame_result[0] if isinstance(frame_result,tuple) else frame_result
                    ids = np.array(chunk['record_indices'],dtype=int)
                    assert len(features['trans_1']) == len(ids)
                    with torch.inference_mode():
                        results = model.predict_from_frames(
                            translations=[features['trans_1']], rotations=[features['rotmats_1']],
                            batch_size=1, stop_grad=True,
                            aatypes=[features['aatype']] if model.use_aatype else None)
                        result = results[0]
                    if 'per_res_covariance' in result:
                        covariance = result['per_res_covariance'].detach().cpu().numpy()
                        trace = np.trace(covariance,axis1=-2,axis2=-1)
                        assert (trace>=0).all()
                        values = 10*np.sqrt(trace)
                    else:
                        values = result['global_rmsf'].detach().cpu().numpy().reshape(-1)
                    assert values.shape == (len(ids),) and np.isfinite(values).all()
                    pred[ids] = values
                    visited[ids] += 1
                    chunk_times.append(dict(context=chunk['context'], n=len(ids), seconds=time.perf_counter()-tick))
                    del result,results,features
                reference = dict(np.load(row['reference_path']))
                assert (visited[reference['frame_mask']]==1).all()
                assert (visited[~reference['frame_mask']]==0).all()
                assert np.isfinite(pred[reference['valid']]).all()
                elapsed = time.perf_counter()-started
                np.savez_compressed(path,prediction=pred,valid=reference['valid'],target=reference['rmsf'],
                    sequence=reference['sequence'],chain=reference['chain'],author_residue_id=reference['author_residue_id'])
                status.update(status='complete',seconds=elapsed,output_sha256=sha(path),
                    predicted=int(np.isfinite(pred).sum()),chunk_times=chunk_times,
                    metrics=metrics(reference['rmsf'][reference['valid']],pred[reference['valid']]))
            except Exception as exc:
                status.update(status='failed',seconds=time.perf_counter()-started,error=str(exc),traceback=traceback.format_exc())
            with (output/'status.jsonl').open('a',encoding='utf-8') as f:
                f.write(json.dumps(status)+'\n')
            print(count,len(rows),row['sample_id'],context,status['status'],round(status['seconds'],3),flush=True)
            if status['status']=='failed':
                print(status['traceback'],flush=True)
    audit['completed_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    (output/'run_metadata.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')


if __name__ == '__main__':
    main()
