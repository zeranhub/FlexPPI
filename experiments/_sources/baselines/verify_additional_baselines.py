"""Independently verify complete RMSF-net or PLANET-MD baseline outputs.

This checker imports no PyTorch or model code and never performs inference.
Only a fully validated3063-record/4973-chain run receives run_audit.json and
exit0. Any missing/failed/corrupt item is retained in verification_failures.csv
and verification_audit.json; an old completion certificate is invalidated.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import time
import numpy as np

EXPECTED_WEIGHTS = {
    'rmsf-net': 'cd12dce9489ed57c5b0b8d5400161e1b4e165963622dc78323fdeb07651d8c9b',
    'planet-full': 'b31a5eb353e773c346a413fbca7e63735a6ab53249a628073a7c01ceee393bf9',
    'planet-seq': '54f8ee6df8ed0c45e33995431b720cb7676b38d7521a5101602fae58f574e243',
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_rows(path, delimiter=','):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def atomic_json(path, value):
    temporary = path.with_suffix('.verification.tmp.json')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
    os.replace(temporary, path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', required=True)
    parser.add_argument('--output', required=True, help='Existing prediction directory, never a new inference destination')
    parser.add_argument('--model', required=True, choices=EXPECTED_WEIGHTS)
    parser.add_argument('--expected-cohort-freeze-sha256', help='Independent hash from dispatched experiment contract')
    parser.add_argument('--inference-script', help='Defaults to corresponding wrapper next to this checker')
    parser.add_argument('--weights', help='Optional deployed checkpoint; if supplied, content hash must match released frozen weights')
    args = parser.parse_args()
    cohort, output = Path(args.cohort), Path(args.output)
    if not output.is_dir():
        raise SystemExit('Prediction directory does not exist: ' + str(output))
    failures, checked, provenance = [], [], {}
    started = time.time()
    audit_path = output/'run_audit.json'
    if audit_path.exists():
        # A completion certificate cannot survive a failed fresh verification.
        invalid = output/f'run_audit.previous_{time.time_ns()}.json'
        os.replace(audit_path, invalid)
        provenance['previous_completion_certificate'] = str(invalid)

    def fail(scope, identifier, error):
        failures.append(dict(scope=scope, identifier=identifier, error=str(error)))

    try:
        freeze_path = cohort/'cohort_freeze.json'
        freeze_sha = sha(freeze_path)
        if args.expected_cohort_freeze_sha256:
            require(freeze_sha == args.expected_cohort_freeze_sha256, 'Dispatched cohort freeze hash mismatch')
        freeze = json.loads(freeze_path.read_text(encoding='utf-8'))
        require(freeze['status'] == 'frozen', 'Cohort is not frozen')
        require(freeze['records'] == 3063 and freeze['chains'] == 4973, 'Wrong complete cohort size')
        require(freeze['valid_rmsf'] == 961321 and freeze['valid_chi_entropy'] == 824988, 'Wrong effective residue counts')
        required_files = ['structure_manifest.tsv','chain_manifest.csv','truth_residues.csv.gz','record_index.json']
        for name in required_files:
            require(name in freeze['files'], 'Missing frozen file hash: ' + name)
            require(sha(cohort/name) == freeze['files'][name], 'Frozen file hash mismatch: ' + name)
        records = read_rows(cohort/'structure_manifest.tsv', '\t')
        chains = read_rows(cohort/'chain_manifest.csv')
        require(len(records) == 3063 and len({r['sample_id'] for r in records}) == 3063, 'Record manifest size or duplicate identity')
        require(len(chains) == 4973 and len({c['input_id'] for c in chains}) == 4973, 'Chain manifest size or duplicate input identity')
        record_lookup = {r['sample_id']:r for r in records}
        chain_lookup = {c['input_id']:c for c in chains}
        require({c['sample_id'] for c in chains} == set(record_lookup), 'Record/chain manifest membership mismatch')
        require(len({(c['sample_id'],c['chain']) for c in chains}) == 4973, 'Duplicate source chain identity')
        for chain in chains:
            require(len(chain['sequence']) == int(chain['n_residues']), 'Sequence length mismatch: ' + chain['input_id'])
        # Use the already frozen truth table to check original residue identity;
        # prediction values and errors never determine coverage or membership.
        record_identity = {sample:[] for sample in record_lookup}
        chain_identity = {key:[] for key in chain_lookup}
        valid_counts, truth_rows = [0,0], 0
        with gzip.open(cohort/'truth_residues.csv.gz','rt',encoding='utf-8',newline='') as handle:
            for row in csv.DictReader(handle):
                sample, key = row['sample_id'], row['input_id']
                require(sample in record_lookup and key in chain_lookup, 'Unknown identity in frozen truth')
                require(chain_lookup[key]['sample_id'] == sample and chain_lookup[key]['chain'] == row['chain'], 'Truth chain source mismatch')
                index = int(row['residue_index'])
                require(index == len(record_identity[sample]), 'Truth record residue order mismatch: ' + sample)
                require(int(row['input_index']) == len(chain_identity[key]), 'Truth chain residue order mismatch: ' + key)
                record_identity[sample].append((row['chain'],row['author_residue_id']))
                chain_identity[key].append(row['amino_acid'])
                valid_counts[0] += int(row['valid_rmsf'])
                valid_counts[1] += int(row['valid_chi_entropy'])
                truth_rows += 1
        require(valid_counts == [961321,824988], 'Truth effective residue counts changed')
        require(truth_rows == freeze['total_residues'], 'Truth total residue count changed')
        for record in records:
            require(len(record_identity[record['sample_id']]) == int(record['n_residues']), 'Truth record shape mismatch')
        for chain in chains:
            require(''.join(chain_identity[chain['input_id']]) == chain['sequence'], 'Truth chain sequence mismatch')
        metadata_path = output/'metadata.json'
        metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
        manifest_name = 'structure_manifest.tsv' if args.model == 'rmsf-net' else 'chain_manifest.csv'
        require(metadata['manifest_sha256'] == freeze['files'][manifest_name], 'Inference manifest hash mismatch')
        require(metadata['weights_sha256'] == EXPECTED_WEIGHTS[args.model], 'Released checkpoint hash mismatch')
        inference_name = 'run_rmsf_net_3063.py' if args.model == 'rmsf-net' else 'run_planet_md_3063.py'
        inference_script = Path(args.inference_script) if args.inference_script else Path(__file__).with_name(inference_name)
        require(metadata['code_sha256'] == sha(inference_script), 'Inference wrapper changed since output metadata was frozen')
        if args.weights:
            require(sha(args.weights) == EXPECTED_WEIGHTS[args.model], 'Deployed checkpoint content mismatch')
        if args.model == 'rmsf-net':
            require(metadata['cohort_freeze_sha256'] == freeze_sha, 'Inference cohort freeze mismatch')
            require(metadata['source_commit'] == 'eb353b0d9d909a6c46d2808467a9d0f16e2bbda7', 'Wrong RMSF-net source commit')
            require(metadata['adapter_sha256'] == sha(inference_script.with_name('rmsf_net_pdb_adapter.py')), 'Density adapter changed since inference')
        else:
            require(metadata['variant'] == ('full' if args.model == 'planet-full' else 'seq'), 'PLANET model variant mismatch')
            require(metadata['rmsf_conversion'] == '10*native_nm to Angstrom, no fitted calibration', 'Unexpected PLANET unit conversion')
        provenance.update(cohort_freeze_sha256=freeze_sha, frozen_files=freeze['files'], inference_metadata_sha256=sha(metadata_path),
            inference_metadata=metadata, inference_script_sha256=sha(inference_script), weights_content_checked=bool(args.weights))
    except Exception as error:
        fail('contract','cohort_or_metadata',type(error).__name__ + ': ' + str(error))

    if not failures:
        entries = records if args.model == 'rmsf-net' else chains
        expected_keys = set()
        for entry in entries:
            if args.model == 'rmsf-net':
                key = 'record_' + hashlib.sha256(entry['sample_id'].encode()).hexdigest()[:24]
                input_path = Path(entry['container_pdb_path'])
            else:
                key = entry['input_id']
                input_path = Path(entry['pdb_file'])
            expected_keys.add(key)
            status_path, prediction_path = output/(key+'.json'), output/(key+'.npz')
            try:
                require(status_path.is_file(), 'Missing per-item status JSON')
                require(prediction_path.is_file(), 'Missing per-item NPZ')
                status = json.loads(status_path.read_text(encoding='utf-8'))
                require(status.get('status') == 'complete', 'Status is ' + str(status.get('status')) + ': ' + str(status.get('error','')))
                require(status['sample_id'] == entry['sample_id'], 'Status sample identity mismatch')
                require(status['input_sha256'] == entry['input_sha256'], 'Status input hash mismatch')
                require(sha(input_path) == entry['input_sha256'], 'Actual frozen structure content mismatch')
                n = int(entry['n_residues'])
                require(int(status['residues']) == n, 'Status residue count mismatch')
                output_hash = sha(prediction_path)
                require(status['output_sha256'] == output_hash, 'Prediction file hash mismatch')
                elapsed = float(status['total_seconds' if args.model == 'rmsf-net' else 'seconds'])
                require(math.isfinite(elapsed) and elapsed >= 0, 'Invalid timing in status')
                with np.load(prediction_path, allow_pickle=False) as saved:
                    if args.model == 'rmsf-net':
                        require(set(['score','chain','author_residue_id']) <= set(saved.files), 'Missing RMSF-net NPZ variables')
                        score = saved['score']
                        require(score.shape == (n,) and score.dtype.kind in 'fiu' and np.isfinite(score).all(), 'Invalid RMSF-net score shape/type/finiteness')
                        identity = record_identity[entry['sample_id']]
                        require(saved['chain'].shape == saved['author_residue_id'].shape == (n,), 'RMSF-net identity array shape mismatch')
                        require(list(zip(saved['chain'].astype(str), saved['author_residue_id'].astype(str))) == identity, 'RMSF-net residue identity/order mismatch')
                    else:
                        require(status['input_id'] == key and status['chain'] == entry['chain'], 'Status chain identity mismatch')
                        require(set(['rmsf_native_nm','rmsf_angstrom','sequence']) <= set(saved.files), 'Missing PLANET NPZ variables')
                        native, converted = saved['rmsf_native_nm'], saved['rmsf_angstrom']
                        require(native.shape == converted.shape == (n,), 'PLANET residue shape mismatch')
                        require(native.dtype.kind in 'fiu' and converted.dtype.kind in 'fiu', 'PLANET nonnumeric prediction array')
                        require(np.isfinite(native).all() and np.isfinite(converted).all(), 'PLANET nonfinite predictions')
                        require(np.allclose(converted,10.0*native,atol=2e-6,rtol=2e-6), 'PLANET nm-to-Angstrom conversion mismatch')
                        require(saved['sequence'].shape == () and str(saved['sequence'].item()) == entry['sequence'], 'PLANET sequence identity mismatch')
                checked.append(dict(input_id=key,sample_id=entry['sample_id'],chain=entry.get('chain'),residues=n,
                    status_sha256=sha(status_path),input_sha256=entry['input_sha256'],output_sha256=output_hash,seconds=elapsed))
            except Exception as error:
                fail('prediction',key,type(error).__name__ + ': ' + str(error))
        prefix = 'record_' if args.model == 'rmsf-net' else 'chain_'
        for extension in ['json','npz']:
            observed = {path.stem for path in output.glob(prefix+'*.'+extension)}
            for key in sorted(observed-expected_keys):
                fail('unexpected_file',key+'.'+extension,'Output is outside the frozen cohort or an unfinished temporary file')

    with (output/'verification_failures.csv').open('w',encoding='utf-8',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=['scope','identifier','error'])
        writer.writeheader()
        writer.writerows(failures)
    required_items = 3063 if args.model == 'rmsf-net' else 4973
    report = dict(status='complete' if not failures and len(checked) == required_items else 'failed',
        model=args.model,expected_records=3063,expected_chains=4973,expected_items=required_items,
        verified_items=len(checked),verified_records=len({row['sample_id'] for row in checked}),
        verified_residues=sum(row['residues'] for row in checked),failure_count=len(failures),
        valid_rmsf=961321,valid_chi_entropy=824988,seconds=time.time()-started,
        checker_sha256=sha(__file__),completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        checks='Frozen manifests/truth hashes, metadata code/weight hashes, every status complete, exact input/output hashes, finite per-residue arrays, frozen residue identity, and PLANET unit conversion.',
        failure_details_file='verification_failures.csv',no_inference_performed=True,no_gpu_used=True,
        provenance=provenance,verified_outputs=checked)
    atomic_json(output/'verification_audit.json',report)
    if report['status'] != 'complete':
        print(json.dumps({key:value for key,value in report.items() if key not in ['provenance','verified_outputs']},indent=2))
        raise SystemExit(1)
    require(report['verified_records'] == 3063 and report['verified_residues'] == 961321,'Final complete-cohort totals mismatch')
    atomic_json(audit_path,report)
    print(json.dumps({key:value for key,value in report.items() if key not in ['provenance','verified_outputs']},indent=2))


if __name__ == '__main__':
    main()
