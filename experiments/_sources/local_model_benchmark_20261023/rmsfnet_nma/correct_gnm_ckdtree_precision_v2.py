"""Versioned GNM adjacency correction to the native ProDy CKDTree float32 path.

Scan EVERY record for graph changes before recomputing only changed graphs.
No original prediction is overwritten, and no reference performance selects a value.
"""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import time
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[name]='2'
sys.dont_write_bytecode=True
import numpy as np
from scipy.linalg import eigvalsh
from scipy.stats import spearmanr
from cohort_common import ROOT, OLD, sha, load_records, record_key, write_json
from run_nma_scipy import fluctuations


def build_gnm_pair_lists(coords):
    xyz32=coords.astype(np.float32)
    edges64=[];edges32=[];changed=[]
    for i in range(len(coords)):
        delta64=coords[i+1:]-coords[i]
        squared64=(delta64**2).sum(1)
        delta32=xyz32[i+1:]-xyz32[i]
        squared32=np.zeros(len(delta32),dtype=np.float32)
        # The C KDTree_dist loop explicitly rounds a float accumulator at each axis.
        for axis in range(3):
            squared32 += delta32[:,axis]*delta32[:,axis]
        mask64=squared64<=225.0
        mask32=squared32<=np.float32(225.0)
        edges64.extend((i,i+1+int(off)) for off in np.flatnonzero(mask64))
        edges32.extend((i,i+1+int(off)) for off in np.flatnonzero(mask32))
        for off in np.flatnonzero(mask64!=mask32):
            j=i+1+int(off)
            changed.append(dict(i=i,j=j,squared_distance_float64=float(squared64[off]),
                squared_distance_ckdtree_float32=float(squared32[off]),included_float64=bool(mask64[off]),included_ckdtree_float32=bool(mask32[off])))
    return np.asarray(edges64,dtype=np.int32),np.asarray(edges32,dtype=np.int32),changed


def kirchhoff_from_edges(count,edges):
    matrix=np.zeros((count,count),dtype=np.float64)
    for i,j in edges:
        matrix[i,j]=-1.0;matrix[j,i]=-1.0
        matrix[i,i]+=1.0;matrix[j,j]+=1.0
    return matrix


def ap90(y,p):
    positives=y>=np.quantile(y,.9)
    order=np.argsort(-p,kind='stable');truth=positives[order].astype(int);scores=p[order]
    boundaries=np.r_[np.flatnonzero(scores[1:]!=scores[:-1]),len(scores)-1]
    tp=np.cumsum(truth)[boundaries];recall=tp/truth.sum();precision=tp/(boundaries+1)
    return float(np.sum(np.diff(np.r_[0.,recall])*precision))


def metrics(y,p):
    return dict(rho=float(spearmanr(y,p).statistic),ap90=ap90(y,p),
        standardized_mae=float(np.mean(np.abs((y-y.mean())/y.std()-(p-p.mean())/p.std()))))


def main():
    tick=time.perf_counter()
    rows,_=load_records(False)
    v1=ROOT/'nma_scipy/predictions'
    output=ROOT/'gnm_ckdtree_v2'
    predictions=output/'predictions'
    predictions.mkdir(parents=True,exist_ok=True)
    original_audit=json.loads((ROOT/'nma_scipy/full_reference_audit.json').read_text())
    failed_samples={r['sample_id'] for r in original_audit['failed_records']}
    scan=[];diagnostics={};overlays=[]
    for index,record in enumerate(rows,1):
        sample=record['sample_id'];key=record_key(sample)
        source=v1/(key+'.npz')
        with np.load(source,allow_pickle=False) as saved:
            coords=saved['ca_coords'].copy();chain=saved['chain'].copy();author=saved['author_residue_id'].copy()
        old_status=json.loads((v1/(key+'.json')).read_text())
        edges64,edges32,changed=build_gnm_pair_lists(coords)
        item=dict(sample_id=sample,dataset=record['cohort'],residues=len(coords),
            edges_float64=len(edges64),edges_ckdtree_float32=len(edges32),changed_edge_count=len(changed),
            edges_float64_sha256=hashlib.sha256(edges64.tobytes()).hexdigest(),
            edges_ckdtree_float32_sha256=hashlib.sha256(edges32.tobytes()).hexdigest())
        scan.append(item)
        if changed or sample in failed_samples:
            matrix=kirchhoff_from_edges(len(coords),edges32)
            values=eigvalsh(matrix)
            nonzero=values[values>=1e-6]
            retained=min(20,len(coords)-1,3*len(coords)-6)
            eig20=float(nonzero[retained-1]);eig21=float(nonzero[retained])
            for pair in changed:
                pair.update(chain_i=str(chain[pair['i']]),chain_j=str(chain[pair['j']]),
                    author_i=str(author[pair['i']]),author_j=str(author[pair['j']]))
            diagnostics[sample]=dict(**item,changed_pairs=changed,zero_modes_v1=old_status['gnm_zero_modes'],
                zero_modes_ckdtree_float32=int(np.sum(values<1e-6)),eigenvalue20=eig20,eigenvalue21=eig21,
                spectrum_gap=eig21-eig20,cutoff_eigenvalue_multiplicity=int(np.sum(np.isclose(nonzero,eig20,atol=1e-10,rtol=1e-12))),
                retained_modes=retained,degenerate_cutoff=bool(np.isclose(eig20,eig21,atol=1e-10,rtol=1e-12)))
        if changed:
            forward_tick=time.perf_counter()
            gnm,variance,info=fluctuations(matrix,retained,1,False)
            path=predictions/(key+'.npz')
            np.savez_compressed(path,gnm=gnm,gnm_sq_fluct=variance,chain=chain,author_residue_id=author,
                                residue_index=np.arange(len(coords)),ca_coords=coords)
            status=dict(sample_id=sample,dataset=record['cohort'],status='complete',residues=len(coords),
                input_sha256=record['input_sha256'],source_v1_path=str(source),source_v1_sha256=sha(source),
                output_sha256=sha(path),new_local_recalculation=True,gnm_modes=info['modes'],
                gnm_zero_modes=info['zero_modes'],recalculation_seconds=time.perf_counter()-forward_tick,
                reason='Match official default CKDTree float32 coordinate and distance accumulator precision; not selected by benchmark score')
            write_json(predictions/(key+'.json'),status)
            overlays.append(dict(sample_id=sample,path=str(path),status_path=str(predictions/(key+'.json'))))
        if index%250==0:
            print(index,len(rows),'changed_graph_records',len(overlays),flush=True)
    with (output/'all_record_adjacency_audit.csv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(scan[0]));writer.writeheader();writer.writerows(scan)
    # Read reference only AFTER determining the implementation-driven correction set.
    reference_rows={sample:[] for sample in diagnostics}
    with gzip.open(OLD/'prody_common_evaluation/common_residue_predictions.csv.gz','rt',encoding='utf-8',newline='') as handle:
        for row in csv.DictReader(handle):
            if row['sample_id'] in reference_rows:
                reference_rows[row['sample_id']].append(row)
    for sample,refs in reference_rows.items():
        key=record_key(sample)
        with np.load(v1/(key+'.npz'),allow_pickle=False) as saved:
            before=saved['gnm'].copy()
        path=predictions/(key+'.npz')
        if path.exists():
            with np.load(path,allow_pickle=False) as saved:
                after=saved['gnm'].copy()
        else:
            after=before
        old=np.asarray([float(r['gnm']) for r in refs]);target=np.asarray([float(r['target_rmsf']) for r in refs])
        old_metrics=metrics(target,old);before_metrics=metrics(target,before);after_metrics=metrics(target,after)
        diag=diagnostics[sample]
        difference=np.abs(after-old)
        diag.update(old_metrics=old_metrics,v1_metrics=before_metrics,v2_metrics=after_metrics,
            v2_minus_old_metrics={m:after_metrics[m]-old_metrics[m] for m in old_metrics},
            v2_minus_v1_metrics={m:after_metrics[m]-before_metrics[m] for m in old_metrics},
            v2_reference_maxabs=float(difference.max()),v2_reference_positions_above_1e_6=int(np.sum(difference>1e-6)),
            old_total_sq_fluct=float(np.sum(old**2)),v2_total_sq_fluct=float(np.sum(after**2)))
    write_json(output/'difference_diagnosis_and_metric_impact.json',dict(records=list(diagnostics.values()),
        interpretation='Nine graph-identical degenerate-cutoff records are retained with disclosed SciPy basis dependence. Changed adjacency graphs alone receive versioned CKDTree-precision recalculation.',
        script_sha256=sha(__file__)))
    metadata=dict(version='gnm_ckdtree_float32_v2',source_commit='08e4e0a5168551a51ab3a9bb2fb97d462236454a',
        source_code_hashes={p.name:sha(p) for p in (ROOT/'upstream_prody/kdtree').glob('*')},
        official_path='gnm.GNMBase.buildKirchhoff default kdtree=True -> ProDy _CKDTree: coordinate and sum float32; Kirchhoff float64',
        geometry_rule='Each coordinate cast float32, sequential float32 squared-distance accumulation, <=float32(15)^2; gamma1; same20-mode solver',
        scanned_records=len(rows),changed_graph_records=len(overlays),overlays=overlays,
        base_predictions=str(v1),unchanged_graphs_keep_fresh_local_v1=True,original_v1_untouched=True,
        benchmark_performance_used_for_selection=False,new_local_recalculations=True,
        wall_seconds=time.perf_counter()-tick,script_sha256=sha(__file__))
    write_json(output/'overlay_manifest.json',metadata)
    print(json.dumps(metadata,indent=2),flush=True)


if __name__=='__main__':
    main()
