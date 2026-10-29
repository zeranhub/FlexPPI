"""Local CPU DynaProt-M, using the published checkpoint on every frozen chain.

Input PDB bytes and sequence must match the original frozen cohort. Outputs are
atomic and resumable. No fitting, calibration, residue selection, or ensembles.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent
COHORT = Path('E:/Project/FlexPPI/experiments/_sources/baselines/cohort')
for key in ['OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[key] = '2'
os.environ['HF_HOME'] = str(ROOT / 'hf_cache')
os.environ['HF_HUB_DISABLE_XET'] = '1'

import numpy as np
import torch
from dynaprot.inference import DynaProtPredictor


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()
    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    torch.manual_seed(20261023)
    manifest = COHORT / 'chain_manifest.csv'
    records = list(csv.DictReader(manifest.open(encoding='utf-8', newline='')))
    assert len(records) == 4973
    if args.smoke:
        systems = {r['sample_id']: r for r in json.loads((COHORT / 'record_index.json').read_text())}
        # One first chain from each source; smoke data are never used as benchmark summaries.
        source = {}
        for row in records:
            entry = systems[row['sample_id']]
            dataset = entry.get('dataset', entry.get('cohort', ''))
            source.setdefault(dataset, row)
        records = list(source.values())
        assert len(records) == 3, list(source)
    output = ROOT / 'dynaprot' / ('smoke' if args.smoke else 'predictions')
    output.mkdir(parents=True, exist_ok=True)
    checkpoint = ROOT / 'weights/dynaprot/dynaprot-M'
    metadata = dict(model='DynaProt-M', upstream='https://github.com/MihirBafna/dynaprot',
        upstream_commit='e1bc8e9f5d80afd91f79ea6e8def8029121a8e6e',
        checkpoint_sha256=sha(checkpoint / 'model.ckpt'),
        config_sha256=sha(checkpoint / 'model_config.yaml'),
        source_manifest_sha256=sha(manifest), script_sha256=sha(__file__),
        input_protocol='Published independent-chain inference; every frozen input chain; canonical chain PDB bytes verified before prediction',
        units='Angstrom, official sqrt(trace(predicted marginal covariance))',
        device='cpu', threads=2, torch_version=torch.__version__, numpy_version=np.__version__,
        fit=False, calibration=False, only_marginal_checkpoint=True)
    metadata_path = output / 'metadata.json'
    if metadata_path.exists():
        assert json.loads(metadata_path.read_text()) == metadata, 'Parameters or source changed'
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    start_load = time.perf_counter()
    predictor = DynaProtPredictor.from_pretrained(marginal=checkpoint, joint=None, device='cpu')
    print('Loaded DynaProt-M in %.3fs, parameters=%d' %
          (time.perf_counter() - start_load, sum(p.numel() for p in predictor.marginal_model.parameters())), flush=True)
    for i, row in enumerate(records):
        identifier = row['input_id']
        path = output / (identifier + '.npz')
        status_path = output / (identifier + '.json')
        if status_path.exists() and path.exists():
            old = json.loads(status_path.read_text())
            if old.get('status') == 'complete':
                assert sha(path) == old['output_sha256']
                continue
        pdb = COHORT / 'chains' / Path(row['pdb_file']).name
        status = dict(input_id=identifier, sample_id=row['sample_id'], chain=row['chain'], status='running',
                      expected_residues=int(row['n_residues']), input_sha256=row['input_sha256'])
        tick = time.perf_counter()
        try:
            assert sha(pdb) == row['input_sha256']
            prediction = predictor.predict(pdb, chain_id=None, outputs='rmsf')
            assert prediction.sequence == row['sequence']
            values = prediction.rmsf.detach().cpu().numpy()
            assert values.shape == (int(row['n_residues']),) and np.isfinite(values).all()
            temporary = output / (identifier + '.tmp.npz')
            np.savez_compressed(temporary, rmsf_angstrom=values,
                                sequence=np.asarray(prediction.sequence),
                                residue_index=prediction.residue_index)
            os.replace(temporary, path)
            status.update(status='complete', output_sha256=sha(path),
                          model_seconds=prediction.timings['marginal'])
        except Exception as error:
            status.update(status='failed', error=type(error).__name__ + ': ' + str(error))
        status['wall_seconds'] = time.perf_counter() - tick
        status_path.write_text(json.dumps(status, indent=2), encoding='utf-8')
        if args.smoke or (i + 1) % 25 == 0 or status['status'] == 'failed':
            print(i + 1, len(records), row['sample_id'], identifier, status['status'],
                  round(status['wall_seconds'], 3), status.get('error', ''), flush=True)
        if args.smoke and status['status'] == 'failed':
            raise RuntimeError(status['error'])
    statuses = [json.loads(p.read_text()) for p in output.glob('chain_*.json')]
    summary = dict(expected=len(records), complete=sum(s['status'] == 'complete' for s in statuses),
                   failed=sum(s['status'] == 'failed' for s in statuses),
                   wall_seconds=sum(s.get('wall_seconds', 0) for s in statuses), metadata=metadata)
    (output / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)
    assert summary['complete'] == len(records) and summary['failed'] == 0


if __name__ == '__main__':
    main()
