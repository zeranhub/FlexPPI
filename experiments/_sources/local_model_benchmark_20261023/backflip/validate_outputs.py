"""Audit fresh BackFlip outputs against the complete, frozen residue map."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tag', required=True)
    args = parser.parse_args()
    rows = json.loads((ROOT/'manifest.json').read_text(encoding='utf-8'))
    output = ROOT/'results'/args.tag
    statuses = {}
    for line in (output/'status.jsonl').read_text(encoding='utf-8').splitlines():
        status = json.loads(line)
        statuses[status['sample_id'], status['context']] = status
    audit = dict(tag=args.tag, expected_records=3063, expected_valid_residues=961321,
                 contexts={}, missing=[], failed=[], checked_output_hashes=0,
                 singleton_shared_prediction_arrays_equal=True,
                 audit_script_sha256=sha(Path(__file__)))
    for context in ['whole_complex', 'independent_chain']:
        counts = dict(completed_records=0, valid_residues=0, groups={},
                      actual_inference_records=0, shared_single_chain_records=0)
        for row in rows:
            key = row['sample_id'], context
            path = output/(row['sample_id']+'__'+context+'.npz')
            status = statuses.get(key)
            if status is None or not path.exists():
                audit['missing'].append(dict(sample_id=key[0], context=context))
                continue
            if status['status'] != 'complete':
                audit['failed'].append(status)
                continue
            assert sha(path) == status['output_sha256'], path
            audit['checked_output_hashes'] += 1
            with np.load(row['reference_path']) as truth, np.load(path) as saved:
                assert len(saved['prediction']) == row['n_residues']
                assert np.array_equal(truth['rmsf'], saved['target'], equal_nan=True)
                assert np.array_equal(truth['valid'], saved['valid'])
                for field in ['sequence', 'chain', 'author_residue_id']:
                    assert np.array_equal(truth[field], saved[field]), (path, field)
                prediction = saved['prediction']
                assert np.isfinite(prediction[truth['frame_mask']]).all(), path
                assert np.isnan(prediction[~truth['frame_mask']]).all(), path
                counts['valid_residues'] += int(saved['valid'].sum())
                if context == 'independent_chain' and row['n_chains'] == 1:
                    with np.load(output/(row['sample_id']+'__whole_complex.npz')) as whole:
                        assert np.array_equal(whole['prediction'], prediction, equal_nan=True)
                    assert status['seconds'] == 0 and 'reused_fresh_context' in status
            group = ('tcr_pmhc_itrace' if row['dataset'] == 'itrace' else
                     'other_multichain' if row['n_chains'] > 1 else 'single_chain')
            counts['groups'][group] = counts['groups'].get(group, 0) + 1
            counts['completed_records'] += 1
            counts['shared_single_chain_records'] += int('reused_fresh_context' in status)
            counts['actual_inference_records'] += int('reused_fresh_context' not in status)
        audit['contexts'][context] = counts
    audit['complete'] = (not audit['missing'] and not audit['failed'] and
                         all(x['completed_records'] == 3063 and x['valid_residues'] == 961321
                             for x in audit['contexts'].values()))
    (output/'completion_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == '__main__':
    main()
