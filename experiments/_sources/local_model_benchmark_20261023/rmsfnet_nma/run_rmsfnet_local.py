"""New local inference using the exact frozen PDB-density adapter and weights."""
import argparse
import importlib.util
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
from cohort_common import ROOT, COHORT, OLD, sha, load_records, local_input, record_key, parse_ca, write_json, valid_cached

SOURCE = Path('E:/Project/FlexPPI/tmp/third_party/RMSF-net')
ADAPTER = OLD / 'rmsf_net_pdb_adapter.py'
WEIGHTS = SOURCE / 'model_cpu/sim_model.pth'
EXPECTED_ADAPTER = 'd4e3606d5d7678251e7615e87987dadfb65a2162b9dd0ab5427368fe13e5c6e1'
EXPECTED_MODEL = '5a59cfde470c3b1ec18397517dfbf1f6d9bb1e1284beb69c2dffc308e1ea22e8'
EXPECTED_WEIGHTS = 'cd12dce9489ed57c5b0b8d5400161e1b4e165963622dc78323fdeb07651d8c9b'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--smoke', action='store_true')
    ap.add_argument('--batch-size', type=int, choices=[1, 2], default=2)
    args = ap.parse_args()
    rows, freeze = load_records(args.smoke)
    assert sha(ADAPTER) == EXPECTED_ADAPTER
    assert sha(SOURCE / 'model.py') == EXPECTED_MODEL
    assert sha(WEIGHTS) == EXPECTED_WEIGHTS
    old_meta = json.loads((OLD / 'rmsf_net/metadata.json').read_text())
    assert old_meta['weights_sha256'] == EXPECTED_WEIGHTS
    assert old_meta['adapter_sha256'] == EXPECTED_ADAPTER
    sys.path.insert(0, str(SOURCE))
    spec = importlib.util.spec_from_file_location('frozen_rmsfnet_adapter', ADAPTER)
    adapter = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = adapter
    start_import = time.perf_counter()
    spec.loader.exec_module(adapter)
    adapter.torch.set_num_threads(2)
    adapter.torch.set_num_interop_threads(1)
    adapter.torch.manual_seed(20261023)
    output = ROOT / 'rmsfnet' / ('smoke' if args.smoke else 'predictions')
    output.mkdir(parents=True, exist_ok=True)
    metadata = dict(model='RMSF-net PDB-only', device='cpu', threads=2, interop_threads=1,
        batch_size=args.batch_size, source_commit=adapter.SOURCE_COMMIT,
        adapter_sha256=sha(ADAPTER), model_code_sha256=sha(SOURCE / 'model.py'), weights_sha256=sha(WEIGHTS),
        script_sha256=sha(__file__), common_sha256=sha(Path(__file__).with_name('cohort_common.py')),
        manifest_sha256=sha(COHORT / 'structure_manifest.tsv'), freeze_sha256=sha(COHORT / 'cohort_freeze.json'),
        old_output_manifest_sha256=old_meta['manifest_sha256'],
        density='Frozen adapter:4A molmap Gaussian,1.5A spacing,12A padding,40^3 boxes,10^3 cores',
        output_units='Relative pre-normalization expm1(log(10)*network_output) averaged over protein heavy atoms; no Angstrom calibration',
        input_protocol='Complete original assembly, canonical bytes verified, original CA residue order and author identities retained',
        torch_version=adapter.torch.__version__, numpy_version=np.__version__, python=sys.version,
        new_local_inference=True, posthoc_fit=False)
    metadata_path = output / 'metadata.json'
    if metadata_path.exists():
        assert json.loads(metadata_path.read_text()) == metadata
    write_json(metadata_path, metadata)
    tick = time.perf_counter()
    model = adapter.load_model(WEIGHTS)
    print('Load model %.3fs; adapter import %.3fs' % (time.perf_counter()-tick, tick-start_import), flush=True)
    new_statuses = []
    for i, record in enumerate(rows, 1):
        key = record_key(record['sample_id'])
        path, status_path = output / (key + '.npz'), output / (key + '.json')
        if valid_cached(path, status_path, record, ['score']):
            new_statuses.append(json.loads(status_path.read_text()))
            continue
        status = dict(sample_id=record['sample_id'], dataset=record['cohort'], input_sha256=record['input_sha256'],
                      status='running', expected_residues=int(record['n_residues']))
        write_json(status_path, status)
        tick = time.perf_counter()
        try:
            pdb = local_input(record)
            xyz, ca_chain, ca_author, ca_name = parse_ca(pdb)
            assert len(xyz) == int(record['n_residues'])
            residues, values, timing = adapter.run_sample(model, str(pdb), args.batch_size)
            values = np.asarray(values, dtype=np.float64)
            chain = np.asarray([r.chain for r in residues])
            author = np.asarray([str(int(r.number)) + r.insertion for r in residues])
            amino = np.asarray([r.name for r in residues])
            assert np.array_equal(chain, ca_chain) and np.array_equal(author, ca_author) and np.array_equal(amino, ca_name)
            assert values.shape == (len(xyz),) and np.isfinite(values).all()
            temp = output / (key + '.tmp.npz')
            np.savez_compressed(temp, score=values, chain=chain, author_residue_id=author,
                                amino_acid_3=amino, residue_index=np.arange(len(xyz)), ca_coords=xyz)
            os.replace(temp, path)
            status.update(status='complete', residues=len(xyz), output_sha256=sha(path), **timing)
        except Exception as exc:
            status.update(status='failed', error=type(exc).__name__ + ': ' + str(exc))
        status['total_seconds'] = time.perf_counter() - tick
        write_json(status_path, status)
        new_statuses.append(status)
        if args.smoke or i % 25 == 0 or status['status'] == 'failed':
            print(i, len(rows), record['sample_id'], status['status'], round(status['total_seconds'],3), status.get('error',''), flush=True)
        if args.smoke and status['status'] == 'failed':
            raise RuntimeError(status['error'])
    summary = dict(expected=len(rows), complete=sum(s['status']=='complete' for s in new_statuses),
                   failed=sum(s['status']=='failed' for s in new_statuses),
                   timed_seconds_sum=sum(s.get('total_seconds',0) for s in new_statuses),
                   total_boxes=sum(s.get('selected_box_count',0) for s in new_statuses),
                   total_residues=sum(s.get('residues',0) for s in new_statuses), metadata=metadata)
    write_json(output / 'summary.json', summary)
    print(json.dumps(summary, indent=2), flush=True)
    assert summary['complete'] == len(rows) and summary['failed'] == 0


if __name__ == '__main__':
    main()
