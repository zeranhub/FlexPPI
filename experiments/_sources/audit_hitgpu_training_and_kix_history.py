"""CPU/read-only checkpoint and historical trajectory audit; no inference/MD."""
from pathlib import Path
import datetime
import hashlib
import json
import math
import os
import re
import subprocess

ROOT=Path('/root/flexppi_completion_20261010')
TRAIN=Path('/root/flexppi_revision_20261002/corrected_training_v1')
SINGLE=ROOT/'single_target'
MD=ROOT/'kix_md_triplicate_rebuild_20261010'
os.environ['CUDA_VISIBLE_DEVICES']=''
os.environ['OMP_NUM_THREADS']=os.environ['MKL_NUM_THREADS']=os.environ['OPENBLAS_NUM_THREADS']='2'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def checkpoint_audit(path):
    import torch
    torch.set_num_threads(2)
    checkpoint=torch.load(path,map_location='cpu',weights_only=False)
    weights=checkpoint['model']
    bad=[key for key,tensor in weights.items() if tensor.is_floating_point() and not bool(torch.isfinite(tensor).all())]
    optimizer=[value for row in checkpoint['optimizer']['state'].values() for value in row.values() if torch.is_tensor(value) and value.is_floating_point()]
    cfg=checkpoint['config']
    return dict(path=str(path),sha256=sha(path),epoch=checkpoint['epoch'],global_step=checkpoint['global_step'],
        train_loss=checkpoint['train_loss'],train_loss_finite=math.isfinite(checkpoint['train_loss']),
        val_loss_is_nan=math.isnan(checkpoint['val_loss']),val_fraction=cfg['train']['val_fraction'],
        validation_state_manifest=cfg['data']['validation_state_manifest'],seed=cfg['seed'],metrics=cfg['model']['metrics'],
        from_scratch=cfg['init_checkpoint'] is None,model_tensor_count=len(weights),nonfinite_model_tensor_names=bad,
        optimizer_float_tensors=len(optimizer),optimizer_all_finite=all(bool(torch.isfinite(value).all()) for value in optimizer),
        all_history_train_losses_finite=all(math.isfinite(row['train_loss']) for row in checkpoint['history']))


def main():
    queue=read(SINGLE/'single_target_queue_state.json')
    completed=[row for row in queue['jobs'] if row['status'].startswith('complete')]
    current=next((row for row in queue['jobs'] if row['status']=='running'),None)
    audits=[checkpoint_audit(Path(row['final_checkpoint'])) for row in completed]
    log=str(SINGLE/'logs'/f"{current['arm']}_fold{current['fold']}.log") if current else None
    lines=Path(log).read_text(errors='replace').splitlines() if log else []
    epochs=[int(match.group(1)) for line in lines if (match:=re.search(r'^epoch=(\d+)',line))]
    batch_bad=[line for line in lines if line.startswith('epoch=') and re.search(r"'loss':\s*(?:nan|[+-]?inf)\b",line,re.I)]
    traces=[line for line in lines if 'Traceback (most recent call last)' in line or 'CUDA out of memory' in line]
    if current:
        paths=sorted(Path(current['output_dir']).glob('epoch_*.pt'))
        if paths:audits.append(checkpoint_audit(paths[-1]))
    guard=read(SINGLE/'frozen_source_guard.json')
    source_mismatch=[name for name,expected in guard['source_hashes'].items() if sha(TRAIN/name)!=expected]
    plan=read(SINGLE/'single_target_plan.json')
    input_mismatch=[name for name,expected in plan['input_file_hashes'].items() if sha(TRAIN/name)!=expected]
    unhealthy=bool(batch_bad or traces or source_mismatch or input_mismatch or any(
        not row['train_loss_finite'] or row['nonfinite_model_tensor_names'] or not row['optimizer_all_finite'] or not row['all_history_train_losses_finite'] for row in audits))
    # Search actual filesystem trajectories; plans or transferred audit notes do
    # not certify a historical production run. Exclude only dependency/cache trees.
    search_command=['find','/root','/workspace','-path','*/site-packages','-prune','-o','-path','/root/.cache','-prune','-o',
        '-path','*/hf_cache','-prune','-o','-type','f','(','-iname','*kix*','-o','-iname','*.xtc','-o','-iname','*.dcd','-o','-iname','*.trr',')','-print']
    found=subprocess.run(search_command,capture_output=True,text=True,check=True).stdout.strip().splitlines()
    trajectory_paths=[value for value in found if Path(value).suffix.lower() in ['.xtc','.dcd','.trr']]
    dcd_code="""
from pathlib import Path
import json
from MDAnalysis.coordinates.DCD import DCDReader
root=Path('/root/flexppi_completion_20261010/kix_md_triplicate_rebuild_20261010/smoke')
rows=[]
for state in ['apo','myb','mll','both']:
    path=root/(state+'_protein.dcd')
    with DCDReader(str(path)) as reader:
        first=float(reader[0].time);last=float(reader[-1].time)
        rows.append(dict(state=state,path=str(path),frames=len(reader),atoms=reader.n_atoms,
            dt_ps=float(reader.dt),first_saved_time_ps=first,last_saved_time_ps=last,saved_span_ps=last-first))
print(json.dumps(rows))
"""
    dcd=subprocess.run([str(ROOT/'md_conda/bin/python'),'-'],input=dcd_code,capture_output=True,text=True,
                       env=dict(os.environ,CUDA_VISIBLE_DEVICES='',PYTHONDONTWRITEBYTECODE='1'),check=True)
    smoke_rows=json.loads(dcd.stdout)
    for row in smoke_rows:
        metadata=read(MD/'smoke'/(row['state']+'.json'))
        row.update(trajectory_sha256=sha(row['path']),metadata_path=str(MD/'smoke'/(row['state']+'.json')),
            metadata_sha256=sha(MD/'smoke'/(row['state']+'.json')),
            simulated_endpoint_ps=metadata['samples'][-1]['time_ps'],classification='40ps_technical_smoke_only',
            energy_temperature_samples_finite=all(math.isfinite(sample[key]) for sample in metadata['samples'] for key in ['potential_kj_mol','kinetic_kj_mol','temperature_k','time_ps']))
    frozen=read(MD/'frozen_predictions/FROZEN.json')
    primary=next(row for row in frozen['models'] if row['method']=='full' and row['fold']==frozen['primary_fold'])
    assert sha(primary['checkpoint'])==primary['sha256']
    assert sha(MD/'frozen_predictions/residue_predictions.csv')==frozen['prediction_csv_sha256']
    report=dict(status='current_training_no_numerical_failure_observed' if not unhealthy else 'review_required',
        checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),audit_code_sha256=sha(__file__),
        training=dict(single_target_completed_jobs=len(completed),single_target_planned_jobs=10,
            family_completed_jobs=len(list((ROOT/'family_folds/runs').glob('*/fold_*/epoch_040.pt'))),family_planned_jobs=10,
            current_job=current,current_log_path=log,current_epoch_observed=max(epochs) if epochs else None,
            current_log_tail=lines[-4:],nonfinite_batch_loss_lines=batch_bad,traceback_or_cuda_oom_lines=traces,
            checkpoint_audits=audits,source_files_checked=len(guard['source_hashes']),source_mismatches=source_mismatch,
            frozen_manifest_mismatches=input_mismatch,
            val_nan_explanation='val_fraction=0 and no validation manifest produce an empty validation set; val_loss NaN is the expected placeholder. Fixed epoch040 defines all evaluations.',
            no_inference_no_optimizer_step_no_gpu=True),
        cancellation=read(ROOT/'pipeline_md_cancellation_audit_20261010.json'),
        historical_kix=dict(search_roots=['/root','/workspace'],search_command=search_command,all_matching_files=found,
            actual_trajectory_files=trajectory_paths,actual_rebuilt_smoke=smoke_rows,
            new_production_done_files=[str(path) for path in (MD/'runs').glob('*/rep*/DONE.json')],
            new_production_summary_files=[str(path) for path in (MD/'runs').glob('*/rep*/analysis/summary.json')],
            long_historical_kix_production_trajectories_found_on_searched_hit_gpu_roots=False,
            conclusion_zh='已搜索hit-gpu的/root和/workspace；仅发现当前重建批次四个40ps技术烟雾DCD，每个40帧、1ps间隔，没有发现可复用的旧100ns KIX生产轨迹。原始本地及hmu-gpu历史单轨迹由另行审计负责；本服务器烟雾不能替代生产MD或用于宣称收敛。',
            primary_prediction_fold=frozen['primary_fold'],primary_full_checkpoint=primary,
            primary_full_checkpoint_content_hash_verified=True,primary_prediction_csv_sha256=frozen['prediction_csv_sha256'],
            primary_prediction_csv_content_hash_verified=True),
        pipeline_source_sha256=sha(ROOT/'run_completion_pipeline.py'),collector_source_sha256=sha(ROOT/'collect_remote_progress.py'),
        no_docx_modified=True)
    (ROOT/'hitgpu_training_kix_status_audit_20261010.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    compact={key:report[key] for key in ['status','checked_utc','pipeline_source_sha256','collector_source_sha256']}
    compact['training']={key:report['training'][key] for key in ['single_target_completed_jobs','single_target_planned_jobs','family_completed_jobs','current_job','current_epoch_observed','source_mismatches','frozen_manifest_mismatches']}
    compact['checkpoints']=audits
    compact['kix_trajectories']=smoke_rows
    print(json.dumps(compact,ensure_ascii=False,indent=2,allow_nan=False))


if __name__=='__main__':main()
