"""Strict verification of full and seq small-sample actual CPU predictions."""
from pathlib import Path
import datetime
import hashlib
import json
import numpy as np

ROOT = Path('/root/flexppi_completion_20261010/baselines')
OUT = ROOT / 'planet_cpu_smoke_audit_20261010'


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    audit = json.loads((OUT / 'payload_audit.json').read_text())
    outputs, post_files = {}, []
    freeze = json.loads((ROOT / 'planet_payload_freeze.json').read_text())
    for item in freeze['files']:
        if item['path'].startswith(('weights/', 'hf_cache/')) or item['path'].startswith('PLANET-MD/'):
            path = ROOT / item['path']
            actual = sha(path)
            assert actual == item['sha256']
            post_files.append(dict(path=item['path'], sha256=actual, unchanged=True))
    assert sha(ROOT / 'run_planet_md_3063.py') == audit['adapter_sha256']
    for variant in ('full', 'seq'):
        folder = ROOT / f'planet_{variant}_cpu_smoke_20261010'
        metadata_path = folder / 'metadata.json'
        metadata = json.loads(metadata_path.read_text())
        assert metadata['device'] == 'cpu' and metadata['variant'] == variant
        assert metadata['weights_sha256'] == audit['checkpoints'][variant]['sha256']
        assert metadata['structure_feature_stage'] == 'encoded'
        assert metadata['code_sha256'] == audit['adapter_sha256']
        rows = []
        assert len(list(folder.glob('*.npz'))) == 2
        for entry in audit['smoke_inputs']:
            name = entry['input_id']
            npz = folder / (name + '.npz')
            status_path = folder / (name + '.json')
            status = json.loads(status_path.read_text())
            assert status['status'] == 'complete' and status['input_sha256'] == entry['input_sha256']
            assert sha(npz) == status['output_sha256']
            with np.load(npz, allow_pickle=False) as prediction:
                native = prediction['rmsf_native_nm']
                angstrom = prediction['rmsf_angstrom']
                assert native.shape == angstrom.shape == (int(entry['n_residues']),)
                assert np.isfinite(native).all() and np.isfinite(angstrom).all()
                assert np.array_equal(angstrom, 10.0 * native)
                assert str(prediction['sequence'].item()) == entry['sequence']
                assert native.dtype == angstrom.dtype == np.float32
                rows.append(dict(input_id=name, sample_id=entry['sample_id'],
                                 chain=entry['chain'], residues=len(native),
                                 status='complete', actual_prediction=True,
                                 input_sha256=entry['input_sha256'], output_sha256=sha(npz),
                                 dtype=str(native.dtype), nm_to_angstrom_exact=True,
                                 rmsf_native_nm_range=[float(native.min()),float(native.max())],
                                 rmsf_angstrom_range=[float(angstrom.min()),float(angstrom.max())],
                                 prediction_seconds=status['seconds']))
        outputs[variant] = dict(folder=str(folder), metadata=metadata,
                                metadata_sha256=sha(metadata_path), samples=rows,
                                log=str(ROOT/f'planet_{variant}_cpu_smoke_20261010.log'),
                                log_sha256=sha(ROOT/f'planet_{variant}_cpu_smoke_20261010.log'))
    report = dict(status='full_and_seq_cpu_smoke_passed_2_chains_each',
                  created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  environment='/root/flexppi_completion_20261010/planet_env/bin/python',
                  offline_cache_only=True, gpu_used=False, trained=False, calibrated=False,
                  OMP_NUM_THREADS=2,MKL_NUM_THREADS=2,OPENBLAS_NUM_THREADS=2,torch_num_threads=2,
                  released_float32_weights_preserved=True, dtype_conversion=False,
                  frozen_source_cache_and_weights_unchanged_after_smoke=post_files,
                  payload_audit_sha256=sha(OUT/'payload_audit.json'), outputs=outputs,
                  conclusion='Public full and seq checkpoints load and predict all residues of the two frozen input chains under offline CPU inference; full struct_features=True and encoded stage, seq=False. This is a technical smoke only, not the 3063-record evaluation.',
                  limitations=['Independent chain input does not encode partners.',
                               'Released weights may overlap ATLAS/mdCATH; not a clean family retraining baseline.',
                               'Only RMSF is scored here; released SHP is not joint discrete chi-state entropy.',
                               'No GPU context or GPU smoke was executed by this agent.'])
    (OUT / 'smoke_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    (OUT / 'SMOKE_PASSED').write_text(json.dumps({'status':report['status'],'audit_sha256':sha(OUT/'smoke_audit.json')}, indent=2) + '\n')
    print(json.dumps({'status':report['status'], 'gpu_used':False,
                      'variants': {variant: data['samples'] for variant,data in outputs.items()}}, indent=2))


if __name__ == '__main__':
    main()
