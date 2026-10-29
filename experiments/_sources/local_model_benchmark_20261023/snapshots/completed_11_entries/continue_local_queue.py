"""Finite local job continuation, not a scheduler or recurring automation."""
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
import psutil

ROOT=Path(__file__).resolve().parent
PYTHON=ROOT/'runtime_py312/Scripts/python.exe'
QA_PYTHON=Path('C:/Users/DELL/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe')
STATE=ROOT/'controller_state.json'
FLAGS=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0
LEDGER=ROOT/'controller_phase_results.json'
phase_results=[]

def record_result(label,ok,**extra):
    phase_results.append(dict(label=label,ok=bool(ok),updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**extra))
    LEDGER.write_text(json.dumps(phase_results,indent=2),encoding='utf-8')

def save(phase,**extra):
    state=dict(phase=phase,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),pid=os.getpid(),**extra)
    temp=STATE.with_suffix('.tmp');temp.write_text(json.dumps(state,indent=2),encoding='utf-8');temp.replace(STATE)
    print(json.dumps(state),flush=True)

def config():return json.loads((ROOT/'long_queue_config.json').read_text())

def active_previous():
    found=[]
    signatures=['backflip/run_benchmark.py','backflip/run_remaining_variants.ps1',
                'rmsfnet_nma/run_nma_scipy.py','fastprotflex/run_local.py']
    for p in psutil.process_iter(['pid','cmdline']):
        try:
            command=' '.join(p.info['cmdline'] or []).replace('\\','/').lower()
            if str(ROOT).replace('\\','/').lower() in command and any(s in command for s in signatures):
                found.append(p.pid)
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    return found

def run(script,args=(),label=None,interpreter=None):
    label=label or Path(script).stem
    save('running_'+label,script=str(script),arguments=list(args))
    with (ROOT/('controller_'+label+'.log')).open('a',encoding='utf-8') as log:
        proc=subprocess.Popen([str(interpreter or PYTHON),'-u',str(ROOT/script),*args],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,creationflags=FLAGS)
        while proc.poll() is None:
            save('running_'+label,child_pid=proc.pid,free_ram_gib=round(psutil.virtual_memory().available/1024**3,2))
            time.sleep(30)
        save('finished_'+label,return_code=proc.returncode)
        record_result(label,proc.returncode==0,return_code=proc.returncode)
        return proc.returncode==0

def publish(label):
    started=datetime.datetime.now(datetime.timezone.utc).isoformat()
    ok=run('evaluate_benchmark.py',label='score_'+label)
    if ok:
        fresh=json.loads((ROOT/'evaluation/scoring_audit.json').read_text())
        ok=fresh.get('generated_utc','')>=started
    if ok:ok=run('plot_benchmark.py',label='plot_'+label)
    if ok:ok=run('build_benchmark_report.py',label='report_'+label)
    if ok:ok=run('qa_current.py',label='qa_'+label,interpreter=QA_PYTHON)
    run('rmsfnet_nma/summarize_execution_inventory.py',label='inventory_'+label)
    record_result('publish_'+label,ok)
    save('published_'+label if ok else 'publication_failed_'+label)
    return ok

def main():
    save('waiting_existing_local_jobs')
    start=time.monotonic()
    while active_previous():
        if time.monotonic()-start>172800:
            save('attention_existing_jobs_exceeded_48_hours',active_pids=active_previous());return
        save('waiting_existing_local_jobs',active_pids=active_previous(),free_ram_gib=round(psutil.virtual_memory().available/1024**3,2))
        time.sleep(30)
    last_publish_ok=publish('light_models')
    expected=['Full','Intra','Static','AA','DynaProt-M','ANM','GNM','FastProtFlex',
        'backflip-2.1 (assembly)','backflip-2.1 (chains)','backflip-1.0 (assembly)',
        'backflip-1.0 (chains)','backflip-1.0-seq (assembly)','backflip-1.0-seq (chains)']
    if config()['run_planet']:
        expected+=['PLANET-MD','PLANET-MD-seq','PLANET-MD-mini']
        if run('planet/run_local.py',['--mode','smoke','--minimum-free-gib','5','--memory-wait-seconds','600'],label='planet_smoke'):
            run('planet/run_local.py',['--mode','full'],label='planet_full')
        last_publish_ok=publish('planet')
    if config()['run_rmsfnet']:
        expected+=['RMSF-net']
        run('rmsfnet_nma/run_rmsfnet_local.py',label='rmsfnet_full')
        last_publish_ok=publish('rmsfnet')
    if not last_publish_ok:
        save('attention_final_publication_failed',phase_results=phase_results,old_outputs_are_not_current_success=True)
        return
    audit=json.loads((ROOT/'evaluation/scoring_audit.json').read_text())
    missing=[name for name in expected if name not in audit['complete_methods']]
    failed_phases=[row['label'] for row in phase_results if not row['ok']]
    save('finished_selected_cohort' if not missing and not failed_phases else 'finished_with_failures_or_missing_coverage',
         complete_entries=audit['complete_methods'],missing_selected_entries=missing,failed_phases=failed_phases,
         visual_review_status=json.loads((ROOT/'QA_REPORT.json').read_text()).get('status','unknown'),
         report=str(ROOT/'RESULTS_zh.md'),figure=str(ROOT/'figures/benchmark_current.pdf'),deferred_heavy_models=config())

if __name__=='__main__':
    try:main()
    except Exception:
        save('controller_failed',traceback=traceback.format_exc());raise
