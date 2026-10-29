"""ANM/GNM over every frozen complete structure; atomic per-record outputs."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import time

os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('MKL_NUM_THREADS', '2')
import numpy as np
from prody import ANM, GNM, calcSqFlucts, confProDy, parsePDB


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--modes', type=int, default=20)
    parser.add_argument('--cutoff', type=float, default=15.0)
    parser.add_argument('--gamma', type=float, default=1.0)
    parser.add_argument('--limit', type=int)
    parser.add_argument('--shard-index', type=int, default=0)
    parser.add_argument('--shard-count', type=int, default=1)
    parser.add_argument('--retry-failed', action='store_true')
    args = parser.parse_args()
    cohort, output = Path(args.cohort), Path(args.output)
    freeze = json.loads((cohort / 'cohort_freeze.json').read_text())
    manifest_path = cohort / 'structure_manifest.tsv'
    assert sha(manifest_path) == freeze['files']['structure_manifest.tsv']
    assert freeze['records'] == 3063 and freeze['valid_rmsf'] == 961321
    assert 0 <= args.shard_index < args.shard_count and args.modes > 0
    output.mkdir(parents=True, exist_ok=True)
    metadata = dict(manifest_sha256=sha(manifest_path), cohort_freeze_sha256=sha(cohort / 'cohort_freeze.json'),
        code_sha256=sha(__file__), modes=args.modes, cutoff=args.cutoff, gamma=args.gamma,
        input='Complete original protein assembly; all frozen CA residues in original order',
        output_units='Relative sqrt(calcSqFlucts); no fitted Angstrom calibration',
        device='cpu', shard_count=args.shard_count, shard_index=args.shard_index)
    metadata_path = output / f'metadata_shard{args.shard_index}.json'
    if metadata_path.exists():
        assert json.loads(metadata_path.read_text()) == metadata, 'Run parameters changed'
    else:
        metadata_path.write_text(json.dumps(metadata, indent=2))
    with manifest_path.open(encoding='utf-8', newline='') as handle:
        records = list(csv.DictReader(handle, delimiter='\t'))
    if args.limit is not None:
        records = records[:args.limit]
    records = records[args.shard_index::args.shard_count]
    confProDy(verbosity='none')
    for position, record in enumerate(records):
        sample = record['sample_id']
        key = 'record_' + hashlib.sha256(sample.encode()).hexdigest()[:24]
        result_path, status_path = output / (key + '.npz'), output / (key + '.json')
        if result_path.exists() and status_path.exists():
            status = json.loads(status_path.read_text())
            if status.get('status') == 'complete':
                assert status['input_sha256'] == record['input_sha256']
                assert sha(result_path) == status['output_sha256']
                with np.load(result_path) as saved:
                    assert saved['anm'].shape == saved['gnm'].shape == (int(record['n_residues']),)
                    assert np.isfinite(saved['anm']).all() and np.isfinite(saved['gnm']).all()
                continue
        if status_path.exists() and json.loads(status_path.read_text()).get('status') == 'failed' and not args.retry_failed:
            continue
        status = dict(sample_id=sample, dataset=record['cohort'], status='running', input_sha256=record['input_sha256'])
        started = time.perf_counter()
        try:
            pdb = Path(record['container_pdb_path'])
            assert sha(pdb) == record['input_sha256']
            atoms = parsePDB(str(pdb), subset='ca', model=1)
            assert atoms is not None and atoms.numAtoms() == int(record['n_residues'])
            count = atoms.numAtoms()
            assert count >= 4
            n_modes = min(args.modes, count - 1, 3 * count - 6)
            tick = time.perf_counter()
            anm_model = ANM(sample + '_anm')
            anm_model.buildHessian(atoms, cutoff=args.cutoff, gamma=args.gamma)
            anm_model.calcModes(n_modes=n_modes)
            anm_variance = np.asarray(calcSqFlucts(anm_model), dtype=np.float64)
            anm_seconds = time.perf_counter() - tick
            tick = time.perf_counter()
            gnm_model = GNM(sample + '_gnm')
            gnm_model.buildKirchhoff(atoms, cutoff=args.cutoff, gamma=args.gamma)
            gnm_model.calcModes(n_modes=n_modes)
            gnm_variance = np.asarray(calcSqFlucts(gnm_model), dtype=np.float64)
            gnm_seconds = time.perf_counter() - tick
            assert np.isfinite(anm_variance).all() and np.isfinite(gnm_variance).all()
            assert min(anm_variance.min(), gnm_variance.min()) >= -1e-12
            chain = np.asarray([str(c).strip() or '_' for c in atoms.getChids()])
            author = np.asarray([str(int(n)) + str(i).strip() for n, i in zip(atoms.getResnums(), atoms.getIcodes())])
            temporary = output / (key + '.tmp.npz')
            np.savez_compressed(temporary, anm=np.sqrt(np.maximum(anm_variance, 0)), gnm=np.sqrt(np.maximum(gnm_variance, 0)),
                anm_sq_fluct=anm_variance, gnm_sq_fluct=gnm_variance, chain=chain, author_residue_id=author)
            os.replace(temporary, result_path)
            status.update(status='complete', residues=count, anm_modes=anm_model.numModes(), gnm_modes=gnm_model.numModes(),
                anm_seconds=anm_seconds, gnm_seconds=gnm_seconds, total_seconds=time.perf_counter()-started, output_sha256=sha(result_path))
        except Exception as error:
            status.update(status='failed', error=type(error).__name__ + ': ' + str(error), total_seconds=time.perf_counter()-started)
        status_path.write_text(json.dumps(status, indent=2), encoding='utf-8')
        if (position + 1) % 25 == 0 or status['status'] == 'failed':
            print(position+1, len(records), sample, status['status'], f"{status['total_seconds']:.3f}s", flush=True)
    statuses = [json.loads(p.read_text()) for p in output.glob('record_*.json')]
    summary = dict(total_frozen_records=3063, observed_statuses=len(statuses),
        complete=sum(s['status'] == 'complete' for s in statuses), failed=sum(s['status'] == 'failed' for s in statuses),
        timed_seconds_sum=sum(s.get('total_seconds', 0) for s in statuses), run_metadata=metadata)
    (output / f'summary_shard{args.shard_index}.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
