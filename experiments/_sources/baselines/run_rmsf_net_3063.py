"""Run existing audited PDB density adapter, retaining complete-cohort coverage."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--source-root', required=True, help='Official RMSF-net directory containing model.py')
    parser.add_argument('--weights', required=True)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--batch-size', type=int, default=2)
    parser.add_argument('--limit', type=int)
    parser.add_argument('--retry-failed', action='store_true')
    args = parser.parse_args()
    sys.path.insert(0, args.source_root)
    adapter_path = Path(__file__).with_name('rmsf_net_pdb_adapter.py')
    spec = importlib.util.spec_from_file_location('rmsf_net_adapter', adapter_path)
    adapter = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = adapter
    spec.loader.exec_module(adapter)
    adapter.torch.set_num_threads(args.threads)
    cohort, output = Path(args.cohort), Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    freeze = json.loads((cohort / 'cohort_freeze.json').read_text())
    manifest_path = cohort / 'structure_manifest.tsv'
    assert sha(manifest_path) == freeze['files']['structure_manifest.tsv']
    metadata = dict(manifest_sha256=sha(manifest_path), cohort_freeze_sha256=sha(cohort / 'cohort_freeze.json'),
        weights_sha256=sha(args.weights), adapter_sha256=sha(adapter_path), code_sha256=sha(__file__),
        source_commit=adapter.SOURCE_COMMIT, density='Chimera molmap reproduction: resolution4A, spacing1.5A, atomic-number Gaussian, edge12A',
        device='cpu', threads=args.threads, batch_size=args.batch_size,
        output='Relative pre-normalization score expm1(log(10)*network_output), mean over protein heavy atoms; no fitted Angstrom calibration')
    metadata_path = output / 'metadata.json'
    if metadata_path.exists():
        assert json.loads(metadata_path.read_text()) == metadata
    else:
        metadata_path.write_text(json.dumps(metadata, indent=2))
    model = adapter.load_model(Path(args.weights))
    with manifest_path.open(encoding='utf-8') as handle:
        records = list(csv.DictReader(handle, delimiter='\t'))
    if args.limit is not None:
        records = records[:args.limit]
    for position, record in enumerate(records):
        sample = record['sample_id']
        key = 'record_' + hashlib.sha256(sample.encode()).hexdigest()[:24]
        path, status_path = output / (key + '.npz'), output / (key + '.json')
        if path.exists() and status_path.exists():
            prior = json.loads(status_path.read_text())
            if prior.get('status') == 'complete':
                assert prior['input_sha256'] == record['input_sha256'] and sha(path) == prior['output_sha256']
                with np.load(path) as cached:
                    assert cached['score'].shape == (int(record['n_residues']),) and np.isfinite(cached['score']).all()
                continue
        if status_path.exists() and json.loads(status_path.read_text()).get('status') == 'failed' and not args.retry_failed:
            continue
        started = time.perf_counter()
        status = dict(sample_id=sample, input_sha256=record['input_sha256'], status='running')
        try:
            assert sha(record['container_pdb_path']) == record['input_sha256']
            atoms, residues = adapter.parse_pdb(record['container_pdb_path'])
            assert len(residues) == int(record['n_residues'])
            origin, shape = adapter.grid_geometry(atoms)
            density = adapter.simulated_density(atoms, origin, shape)
            atom_indices = adapter.atom_indices_zyx(atoms, origin, density.shape)
            cells = adapter.selected_cells(atom_indices)
            preprocessing = time.perf_counter() - started
            grid_prediction, inference_seconds = adapter.predict_grid(model, density, cells, args.batch_size)
            values = np.asarray(adapter.residue_predictions(atoms, residues, atom_indices, grid_prediction), dtype=np.float64)
            assert np.isfinite(values).all()
            temporary = output / (key + '.tmp.npz')
            np.savez_compressed(temporary, score=values, chain=np.asarray([r.chain for r in residues]), author_residue_id=np.asarray([r.number+r.insertion for r in residues]))
            os.replace(temporary, path)
            status.update(status='complete', residues=len(residues), preprocessing_seconds=preprocessing,
                inference_seconds=inference_seconds, boxes=len(cells), voxels=int(density.size), output_sha256=sha(path))
        except Exception as error:
            status.update(status='failed', error=type(error).__name__ + ': ' + str(error))
        status['total_seconds'] = time.perf_counter() - started
        status_path.write_text(json.dumps(status, indent=2))
        print(position+1, len(records), sample, status['status'], f"{status['total_seconds']:.3f}s", flush=True)


if __name__ == '__main__':
    main()
