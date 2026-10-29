"""Exact selected-box count, without densities or neural inference."""
import csv
import json
import os
from pathlib import Path
import time
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import numpy as np
from cohort_common import ROOT, COHORT, sha, load_records, local_input, record_key, write_json


def main():
    rows, _ = load_records(False)
    output = ROOT/'rmsfnet/profile'
    output.mkdir(parents=True, exist_ok=True)
    geometries=[]
    started=time.perf_counter()
    for index, record in enumerate(rows,1):
        coordinates=[]
        for line in local_input(record).read_text(encoding='utf-8').splitlines():
            if line.startswith('ENDMDL') and coordinates:
                break
            if not line.startswith('ATOM  ') or len(line)<54 or line[16:17] not in (' ','A'):
                continue
            coordinates.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
        coordinates=np.asarray(coordinates,dtype=np.float64)
        origin=np.floor((coordinates.min(0)-12.0)/1.5)*1.5
        upper=np.ceil((coordinates.max(0)+12.0)/1.5)*1.5
        shape=np.rint((upper-origin)/1.5).astype(int)+1
        indices=np.floor((coordinates-origin)/1.5).astype(int)[:,::-1]
        indices=np.clip(indices,0,shape[::-1]-1)
        cell_set={tuple(int(v) for v in row//10) for row in indices}
        result=dict(sample_id=record['sample_id'],dataset=record['cohort'],residues=int(record['n_residues']),
                    selected_box_count=len(cell_set),grid_voxel_count=int(np.prod(shape)),input_sha256=record['input_sha256'])
        smoke_status=ROOT/'rmsfnet/smoke'/(record_key(record['sample_id'])+'.json')
        if smoke_status.exists():
            smoke=json.loads(smoke_status.read_text())
            assert smoke['selected_box_count']==result['selected_box_count']
            assert smoke['grid_voxel_count']==result['grid_voxel_count']
        geometries.append(result)
        if index%250==0:
            print(index,len(rows),flush=True)
    with (output/'complete_queue_geometry.csv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(geometries[0]))
        writer.writeheader()
        writer.writerows(geometries)
    baseline=json.loads((ROOT/'rmsfnet/smoke/summary.json').read_text())
    statuses=[json.loads(p.read_text()) for p in (ROOT/'rmsfnet/smoke').glob('record_*.json')]
    seconds_per_box=sum(s['inference_seconds'] for s in statuses)/sum(s['selected_box_count'] for s in statuses)
    seconds_preprocess_per_residue=sum(s['preprocessing_seconds'] for s in statuses)/baseline['total_residues']
    total_boxes=sum(g['selected_box_count'] for g in geometries)
    total_residues=sum(g['residues'] for g in geometries)
    summary=dict(records=len(geometries),total_residues=total_residues,total_boxes=total_boxes,
        largest_record_boxes=max(g['selected_box_count'] for g in geometries),
        largest_grid_voxels=max(g['grid_voxel_count'] for g in geometries),
        smoke_inference_seconds_per_box=seconds_per_box,smoke_preprocess_seconds_per_residue=seconds_preprocess_per_residue,
        estimated_original_batch2_total_hours=(seconds_per_box*total_boxes+seconds_preprocess_per_residue*total_residues)/3600,
        estimation='Exact complete-assembly selected-box count multiplied by 3-source CPU2-thread smoke time per box; not a finished full run',
        geometry_seconds=time.perf_counter()-started,manifest_sha256=sha(COHORT/'structure_manifest.tsv'),script_sha256=sha(__file__))
    write_json(output/'queue_estimate.json',summary)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
