"""Read-only startup acceptance of both independent, already running CPU queues."""
import csv
import datetime
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

BASE = Path('/root/flexppi_completion_20261010/baselines')
WRAPPER_SHA = '034a7fdfc5b6bc39164bb43f23518d09586863af186ece495275965c3910e602'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def process(pid):
    folder = Path('/proc') / str(pid)
    assert folder.is_dir(), 'Required queue process is no longer live: ' + str(pid)
    environment = dict(item.split('=', 1) for item in (folder/'environ').read_bytes().decode().split('\0') if '=' in item)
    assert environment.get('CUDA_VISIBLE_DEVICES') == ''
    assert all(environment.get(name) == '2' for name in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'])
    device_fds = []
    for path in (folder/'fd').iterdir():
        try:
            target = str(path.resolve(strict=True))
        except (FileNotFoundError, OSError):
            continue
        if target.startswith('/dev/nvidia'):
            device_fds.append(target)
    assert not device_fds, 'CPU process has NVIDIA device descriptors'
    return dict(pid=pid, command=(folder/'cmdline').read_bytes().decode().replace('\0',' ').strip(),
                environment={name:environment.get(name) for name in ['CUDA_VISIBLE_DEVICES','OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']},
                nvidia_device_fds=device_fds,
                namespace_pid_status=[line for line in (folder/'status').read_text().splitlines() if line.startswith(('NSpid:', 'PPid:'))])


def main():
    with (BASE/'cohort'/'chain_manifest.csv').open(newline='') as handle:
        first_two = list(csv.DictReader(handle))[:2]
    assert sha(BASE/'run_planet_cpu_queue.py') == WRAPPER_SHA
    gpu = subprocess.run(['nvidia-smi','--query-compute-apps=pid,process_name,used_gpu_memory','--format=csv,noheader,nounits'],
                         check=True,capture_output=True,text=True)
    gpu_pids = [int(line.split(',')[0]) for line in gpu.stdout.strip().splitlines() if line.strip()]
    variants = {}
    for variant in ['full','seq']:
        state = read(BASE/('planet_'+variant+'_cpu_queue_state.json'))
        assert state['status'] == 'running' and state['phase'] == 'inference'
        assert state['wrapper_sha256'] == WRAPPER_SHA
        parent_pid, child_pid = state['parent_pid'], state['child_pid']
        parent, child = process(parent_pid), process(child_pid)
        assert parent_pid not in gpu_pids and child_pid not in gpu_pids
        worker = read(BASE/('planet_'+variant+'_cpu_worker_audit.json'))
        assert worker['device'] == 'cpu' and worker['intra_op_threads'] == worker['inter_op_threads'] == 2
        assert worker['default_dtype'] == 'torch.float32' and not worker['trained'] and not worker['calibrated']
        output = BASE/('planet_'+variant)
        metadata = read(output/'metadata.json')
        assert metadata['device'] == 'cpu' and metadata['variant'] == variant
        assert metadata['code_sha256'] == state['frozen_contract']['run_planet_md_3063.py']
        assert metadata['weights_sha256'] == state['weights_sha256']
        samples = []
        for entry in first_two:
            key = entry['input_id']
            item = read(output/(key+'.json'))
            assert item['status'] == 'complete' and item['input_id'] == key
            assert item['sample_id'] == entry['sample_id'] and item['chain'] == entry['chain']
            assert item['input_sha256'] == entry['input_sha256'] == sha(entry['pdb_file'])
            assert item['residues'] == int(entry['n_residues'])
            path = output/(key+'.npz')
            assert sha(path) == item['output_sha256']
            with np.load(path,allow_pickle=False) as saved:
                native, converted = saved['rmsf_native_nm'], saved['rmsf_angstrom']
                assert native.dtype == converted.dtype == np.dtype('float32')
                assert native.shape == converted.shape == (int(entry['n_residues']),)
                assert np.isfinite(native).all() and np.isfinite(converted).all()
                assert np.array_equal(converted,10.0*native)
                assert str(saved['sequence'].item()) == entry['sequence']
            samples.append(dict(input_id=key,sample_id=entry['sample_id'],residues=item['residues'],
                                output_sha256=item['output_sha256'],prediction_seconds=item['seconds'],
                                actual_prediction=True,float32=True,finite=True,exact_unit_conversion=True))
        statuses = [read(path)['status'] for path in output.glob('chain_*.json')]
        variants[variant] = dict(status='startup_passed_full_cohort_still_running',parent=parent,child=child,
            complete_items_observed=statuses.count('complete'),failed_items_observed=statuses.count('failed'),
            actual_first_two_predictions=samples,metadata_sha256=sha(output/'metadata.json'),
            worker_audit_sha256=sha(BASE/('planet_'+variant+'_cpu_worker_audit.json')),
            preflight_sha256=sha(BASE/('planet_'+variant+'_cpu_queue_preflight.json')),
            queue_log_snapshot_sha256=sha(BASE/('planet_'+variant+'_cpu_queue.log')),
            full_run_independently_verified=False,run_audit_exists=(output/'run_audit.json').exists())
    report = dict(status='both_cpu_queues_startup_passed',checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        wrapper_sha256=WRAPPER_SHA,auditor_sha256=sha(__file__),variants=variants,
        nvidia_smi_compute_snapshot=gpu.stdout.strip(),gpu_process_pids=gpu_pids,
        no_planet_gpu_context_observed=True,
        gpu_context_evidence='CUDA_VISIBLE_DEVICES empty in both parent/child environments; no NVIDIA device file descriptor in any of those processes; none appears in nvidia-smi compute-app snapshot.',
        no_model_source_or_running_process_changed=True,completion_pending=True)
    (BASE/'planet_cpu_queue_startup_audit_20261010.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
