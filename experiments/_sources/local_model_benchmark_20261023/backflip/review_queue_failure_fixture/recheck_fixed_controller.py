"""Isolated controller review: no native model processes are started."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import datetime
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
EXPECTED = ['Full','Intra','Static','AA','DynaProt-M','ANM','GNM','FastProtFlex',
    'backflip-2.1 (assembly)','backflip-2.1 (chains)','backflip-1.0 (assembly)',
    'backflip-1.0 (chains)','backflip-1.0-seq (assembly)','backflip-1.0-seq (chains)']


def case(name, expected_phase):
    specification = importlib.util.spec_from_file_location('review_'+name, ROOT/'continue_local_queue.py')
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    fixture = HERE/name
    (fixture/'evaluation').mkdir(parents=True, exist_ok=True)
    audit_path = fixture/'evaluation/scoring_audit.json'
    audit_path.write_text(json.dumps(dict(generated_utc='2026-01-01T00:00:00+00:00',
                                          complete_methods=['STALE_OLD_METHOD'])))
    module.ROOT = fixture
    module.STATE = fixture/'controller_state.json'
    module.LEDGER = fixture/'controller_phase_results.json'
    module.active_previous = lambda: []
    module.config = lambda: dict(run_planet=False, run_rmsfnet=False)
    calls = []

    def fake_run(script, args=(), label=None):
        assert script not in ('planet/run_local.py', 'rmsfnet_nma/run_rmsfnet_local.py')
        calls.append(label)
        ok = name != 'all_failed'
        if label.startswith('score_') and ok and name != 'stale_audit':
            entries = EXPECTED[:-1] if name == 'missing_coverage' else EXPECTED
            audit_path.write_text(json.dumps(dict(
                generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                complete_methods=entries)))
        module.record_result(label, ok, return_code=0 if ok else 1)
        return ok

    module.run = fake_run
    module.main()
    state = json.loads(module.STATE.read_text())
    assert state['phase'] == expected_phase, state
    return dict(case=name, actual_phase=state['phase'], expected_phase=expected_phase,
                calls=calls, native_processes_started=0,
                missing=state.get('missing_selected_entries', []))


results = [case('all_failed', 'attention_final_publication_failed'),
           case('stale_audit', 'attention_final_publication_failed'),
           case('fresh_success', 'finished_selected_cohort'),
           case('missing_coverage', 'finished_with_failures_or_missing_coverage')]
comparisons = 0
for microseconds in [0, 1, 999999]:
    before = datetime.datetime(2026, 10, 23, 23, 59, 59, microseconds,
                               tzinfo=datetime.timezone.utc)
    for delta in [datetime.timedelta(0), datetime.timedelta(microseconds=1),
                  datetime.timedelta(milliseconds=1), datetime.timedelta(seconds=1),
                  datetime.timedelta(days=1)]:
        after = before + delta
        assert (after.isoformat() >= before.isoformat()) == (after >= before)
        comparisons += 1
report = dict(cases=results, utc_isoformat_comparison_checks=comparisons,
              production_config_modified=False, native_processes_started=0)
(HERE/'fixed_controller_recheck.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
