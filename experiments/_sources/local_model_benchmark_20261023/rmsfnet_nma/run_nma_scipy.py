"""Local ANM/GNM reproduction of ProDy equations, validated against frozen ProDy output.

No compiled ProDy extension is required. The original ProDy source snapshots are
saved in upstream_prody/ and are the algorithm reference, not an installed model.
Outputs are relative sqrt(sum(eigenvector**2/eigenvalue)), not calibrated Angstroms.
"""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import sys
import time

for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[key] = '2'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.dont_write_bytecode = True

import numpy as np
import scipy
from scipy.linalg import eigh, eigvalsh
from cohort_common import ROOT, COHORT, OLD, sha, load_records, local_input, record_key, parse_ca, write_json, valid_cached

ZERO = 1e-6  # ProDy utilities/eigtools.py constant; do not change to machine epsilon.


def build_matrices(coords, cutoff=15.0, gamma=1.0):
    """Exact pair traversal and accumulation order of ProDy's default dense ANM."""
    count = len(coords)
    hessian = np.zeros((3*count, 3*count), dtype=np.float64)
    kirchhoff = np.zeros((count, count), dtype=np.float64)
    cutoff2 = cutoff * cutoff
    for i in range(count):
        i3 = 3*i
        deltas = coords[i+1:] - coords[i]
        for offset, distance2 in enumerate((deltas**2).sum(1)):
            if distance2 > cutoff2:
                continue
            if distance2 == 0:
                raise ValueError('Coincident CA coordinates have undefined ANM spring direction')
            j = i+1+offset
            j3 = 3*j
            block = np.outer(deltas[offset], deltas[offset]) * (-gamma / distance2)
            hessian[i3:i3+3, j3:j3+3] = block
            hessian[j3:j3+3, i3:i3+3] = block
            hessian[i3:i3+3, i3:i3+3] -= block
            hessian[j3:j3+3, j3:j3+3] -= block
            kirchhoff[i, j] = -gamma
            kirchhoff[j, i] = -gamma
            kirchhoff[i, i] += gamma
            kirchhoff[j, j] += gamma
    return hessian, kirchhoff


def lowest_nonzero_modes(matrix, n_modes, expected_zeros):
    """Reproduce ProDy solveEig(..., zeros=False, turbo=True) dense SciPy branch."""
    dof = len(matrix)
    n_modes = min(n_modes, dof)
    subset = None if n_modes >= dof else (0, n_modes + expected_zeros - 1)
    if subset is None:
        values, vectors = eigh(matrix, driver='evd')
    else:
        values, vectors = eigh(matrix, subset_by_index=subset)
    n_zeros = int(np.sum(values < ZERO))
    if n_zeros > expected_zeros:
        if n_zeros == n_modes + expected_zeros and n_modes < dof:
            n_zeros = int(np.sum(eigvalsh(matrix) < ZERO))
        if n_modes < dof:
            first = min(n_modes + expected_zeros, dof-1)
            last = min(n_modes + n_zeros - 1, dof-1)
            extra_values, extra_vectors = eigh(matrix, subset_by_index=(first, last))
            values = np.concatenate((values, extra_values))
            vectors = np.hstack((vectors, extra_vectors))
    values = values[n_zeros:n_zeros+n_modes]
    vectors = vectors[:, n_zeros:n_zeros+n_modes]
    if len(values) == 0 or np.any(values < ZERO):
        raise ValueError('No stable nonzero modes')
    return values, vectors, n_zeros


def fluctuations(matrix, n_modes, expected_zeros, is_3d):
    values, vectors, n_zeros = lowest_nonzero_modes(matrix, n_modes, expected_zeros)
    variance = np.sum(vectors**2 / values[None, :], axis=1)
    if is_3d:
        variance = variance.reshape((-1, 3)).sum(1)
    assert np.isfinite(variance).all() and variance.min() >= 0
    return np.sqrt(variance), variance, dict(modes=len(values), zero_modes=n_zeros)


def reference_comparison(rows, output):
    reference_path = OLD / 'prody_common_evaluation/common_residue_predictions.csv.gz'
    wanted = {r['sample_id']: r for r in rows}
    reference = {sample: [] for sample in wanted}
    with gzip.open(reference_path, 'rt', encoding='utf-8', newline='') as handle:
        for row in csv.DictReader(handle):
            if row['sample_id'] in wanted:
                reference[row['sample_id']].append(row)
    audit = dict(reference_path=str(reference_path), reference_sha256=sha(reference_path),
                 tolerance_maxabs=1e-6, comparison='Frozen original ProDy CA-level ANM/GNM values, no scaling or fitting', records=[])
    for sample, record in wanted.items():
        refs = reference[sample]
        assert len(refs) == int(record['n_residues'])
        assert [int(r['residue_index']) for r in refs] == list(range(len(refs)))
        with np.load(output / (record_key(sample)+'.npz'), allow_pickle=False) as saved:
            assert np.array_equal(saved['chain'], np.asarray([r['chain'] for r in refs]))
            assert np.array_equal(saved['author_residue_id'], np.asarray([r['author_residue_id'] for r in refs]))
            result = dict(sample_id=sample, dataset=record['cohort'], residues=len(refs))
            for method in ['anm', 'gnm']:
                old_values = np.asarray([float(r[method]) for r in refs], dtype=np.float64)
                difference = np.abs(saved[method] - old_values)
                result[method+'_maxabs'] = float(difference.max())
                result[method+'_meanabs'] = float(difference.mean())
                assert difference.max() < audit['tolerance_maxabs'], result
            audit['records'].append(result)
    audit['pass'] = True
    write_json(output / 'reference_validation.json', audit)
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--verify-reference', action='store_true')
    args = parser.parse_args()
    rows, freeze = load_records(args.smoke)
    output = ROOT / 'nma_scipy' / ('smoke' if args.smoke else 'predictions')
    output.mkdir(parents=True, exist_ok=True)
    if not args.smoke:
        validation = json.loads((ROOT / 'nma_scipy/smoke/reference_validation.json').read_text())
        assert validation['pass'] and len(validation['records']) == 3
        assert {r['dataset'] for r in validation['records']} == {'atlas','dynarepo','itrace'}
        assert validation['tolerance_maxabs'] == 1e-6
    metadata = dict(methods=['ANM','GNM'], implementation='NumPy/SciPy reproduction of original ProDy dense equations',
        upstream_code={p.name:sha(p) for p in sorted((ROOT/'upstream_prody').glob('*.py'))},
        script_sha256=sha(__file__), common_sha256=sha(ROOT/'cohort_common.py'),
        manifest_sha256=sha(COHORT/'structure_manifest.tsv'), freeze_sha256=sha(COHORT/'cohort_freeze.json'),
        output_units='Relative sqrt(calcSqFlucts), no Angstrom calibration', cutoff_angstrom=15.0, gamma=1.0,
        requested_modes=20, zero_threshold=ZERO, numpy_version=np.__version__, scipy_version=scipy.__version__,
        device='cpu', threads=2, input_protocol='Complete original assembly, frozen bytes and original CA identity/order',
        new_local_calculation=True, posthoc_fit=False)
    metadata_path = output/'metadata.json'
    if metadata_path.exists():
        assert json.loads(metadata_path.read_text()) == metadata
    write_json(metadata_path, metadata)
    statuses = []
    for index, record in enumerate(rows, 1):
        key = record_key(record['sample_id'])
        path, status_path = output/(key+'.npz'), output/(key+'.json')
        if valid_cached(path, status_path, record, ['anm','gnm']):
            statuses.append(json.loads(status_path.read_text()))
            continue
        status = dict(sample_id=record['sample_id'], dataset=record['cohort'], input_sha256=record['input_sha256'],
                      status='running', expected_residues=int(record['n_residues']))
        write_json(status_path, status)
        tick = time.perf_counter()
        try:
            xyz, chain, author, amino = parse_ca(local_input(record))
            count = len(xyz)
            assert count == int(record['n_residues'])
            n_modes = min(20, count-1, 3*count-6)
            hessian, kirchhoff = build_matrices(xyz)
            tick_anm = time.perf_counter()
            anm, anm_variance, anm_info = fluctuations(hessian, n_modes, 6, True)
            anm_seconds = time.perf_counter()-tick_anm
            tick_gnm = time.perf_counter()
            gnm, gnm_variance, gnm_info = fluctuations(kirchhoff, n_modes, 1, False)
            gnm_seconds = time.perf_counter()-tick_gnm
            temp = output/(key+'.tmp.npz')
            np.savez_compressed(temp, anm=anm, gnm=gnm, anm_sq_fluct=anm_variance, gnm_sq_fluct=gnm_variance,
                chain=chain, author_residue_id=author, amino_acid_3=amino, residue_index=np.arange(count), ca_coords=xyz)
            os.replace(temp, path)
            status.update(status='complete', residues=count, output_sha256=sha(path),
                anm_modes=anm_info['modes'], gnm_modes=gnm_info['modes'],
                anm_zero_modes=anm_info['zero_modes'], gnm_zero_modes=gnm_info['zero_modes'],
                anm_seconds=anm_seconds, gnm_seconds=gnm_seconds)
        except Exception as exc:
            status.update(status='failed', error=type(exc).__name__+': '+str(exc))
        status['total_seconds'] = time.perf_counter()-tick
        write_json(status_path, status)
        statuses.append(status)
        if args.smoke or index % 25 == 0 or status['status'] == 'failed':
            print(index, len(rows), record['sample_id'], status['status'], round(status['total_seconds'],3), status.get('error',''), flush=True)
        if args.smoke and status['status'] == 'failed':
            raise RuntimeError(status['error'])
    summary = dict(expected=len(rows), complete=sum(s['status']=='complete' for s in statuses),
        failed=sum(s['status']=='failed' for s in statuses), timed_seconds_sum=sum(s['total_seconds'] for s in statuses), metadata=metadata)
    write_json(output/'summary.json', summary)
    print(json.dumps(summary, indent=2), flush=True)
    if args.verify_reference:
        print(json.dumps(reference_comparison(rows, output), indent=2), flush=True)
    assert summary['complete'] == len(rows) and summary['failed'] == 0


if __name__ == '__main__':
    main()
