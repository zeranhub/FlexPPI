"""Wait <=24h for a CPU queue; failed state always overrides old run_audit."""
import argparse
import json
from pathlib import Path
import time
from check_additional_baseline_freeze import COHORT_HASH, HASHES, require, sha, read_json, write_json


def completion(state, output, model):
    require(state['status'] == 'complete', 'Queue is not complete')
    require(state['device'] == 'cpu', 'Queue device changed')
    require(state['variant'] == model.removeprefix('planet-'), 'Queue model mismatch')
    require(Path(state['output']) == output, 'Queue output directory mismatch')
    require(state['expected_chains'] == 4973 and state['expected_records'] == 3063, 'Queue cohort size mismatch')
    audit_path = output / 'run_audit.json'
    audit = read_json(audit_path)
    require(audit['status'] == 'complete' and audit['model'] == model, 'Independent completion audit mismatch')
    require(audit['expected_items'] == audit['verified_items'] == 4973 and audit['verified_records'] == 3063, 'Not all chains and records verified')
    require(audit['verified_residues'] == 961321 and audit['failure_count'] == 0, 'Coverage/failure audit mismatch')
    require(audit['checker_sha256'] == HASHES['verify_additional_baselines.py'], 'Independent verifier changed')
    require(audit['provenance']['cohort_freeze_sha256'] == COHORT_HASH, 'Completion used a different frozen cohort')
    require(audit['provenance']['inference_metadata']['device'] == 'cpu', 'Completed inference was not CPU')
    require(state['output_dtype_audit']['status'] == 'passed', 'Final float32 output sweep failed')
    require(state['output_dtype_audit']['run_audit_sha256'] == sha(audit_path), 'Queue completion certificate was replaced/changed')
    require(sha(output / 'metadata.json') == audit['provenance']['inference_metadata_sha256'], 'Inference metadata changed after verification')
    rows = audit['verified_outputs']
    require(len(rows) == len({row['input_id'] for row in rows}) == 4973, 'Certificate output list is incomplete/duplicated')
    for row in rows:
        key = row['input_id']
        require(sha(output / (key+'.npz')) == row['output_sha256'], 'Prediction changed after full verification: ' + key)
        require(sha(output / (key+'.json')) == row['status_sha256'], 'Per-chain status changed after full verification: ' + key)
    require({path.stem for path in output.glob('chain_*.npz')} == {row['input_id'] for row in rows}, 'Unexpected/missing prediction outputs')
    require({path.stem for path in output.glob('chain_*.json')} == {row['input_id'] for row in rows}, 'Unexpected/missing item statuses')
    require(state['final_frozen_input_recheck']['status'] == 'passed', 'Final execution freeze check failed')
    return sha(audit_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--state', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--model', required=True, choices=['planet-full','planet-seq'])
    parser.add_argument('--timeout-seconds', type=float, default=86400)
    parser.add_argument('--poll-seconds', type=float, default=20)
    parser.add_argument('--wait-audit', required=True)
    args = parser.parse_args()
    require(0 < args.timeout_seconds <= 86400, 'Timeout must be positive and <=24h')
    require(0 < args.poll_seconds <= 20, 'Polling interval must be <=20s')
    state_path, output = Path(args.state), Path(args.output)
    started, last_phase = time.monotonic(), None
    try:
        while True:
            state = read_json(state_path) if state_path.is_file() else {'status':'not_created'}
            # Never read a prior run_audit as sufficient evidence while state is
            # failed/running/missing. Any failed item is actionable immediately.
            require(state.get('status') != 'failed', 'CPU queue failed: ' + str(state.get('error','')))
            require(state.get('observed_items',{}).get('failed',0) == 0, 'CPU queue has failed per-item predictions')
            phase = (state.get('status'),state.get('phase'))
            if phase != last_phase:
                print(json.dumps({'model':args.model,'queue_status':phase[0],'phase':phase[1]},ensure_ascii=False), flush=True)
                last_phase = phase
            if state.get('status') == 'complete':
                audit_hash = completion(state, output, args.model)
                # Re-read state after audit validation so a newer failed state
                # cannot be masked by an older complete snapshot.
                current = read_json(state_path)
                require(current['status'] == 'complete', 'Queue state changed during completion validation')
                report = dict(status='complete', model=args.model, seconds=time.monotonic()-started,
                    state_sha256=sha(state_path), run_audit_sha256=audit_hash, no_gpu_used=True, no_inference_performed=True)
                write_json(Path(args.wait_audit),report)
                print(json.dumps(report,indent=2))
                return
            require(state.get('status') in ['running','not_created'], 'Unknown CPU queue status: ' + str(state.get('status')))
            elapsed = time.monotonic()-started
            require(elapsed < args.timeout_seconds, 'Timed out waiting for independently verified CPU queue')
            time.sleep(min(args.poll_seconds,args.timeout_seconds-elapsed))
    except Exception as error:
        report = dict(status='failed', model=args.model, error=type(error).__name__+': '+str(error),
            seconds=time.monotonic()-started, no_gpu_used=True, no_inference_performed=True)
        write_json(Path(args.wait_audit),report)
        print(json.dumps(report,indent=2),flush=True)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
