"""Own one GPU queue, preserve logs, and stop at any failed completion gate."""
from pathlib import Path
import fcntl,hashlib,json,os,subprocess,sys,time

ROOT=Path(__file__).resolve().parent
TRAIN=Path('/root/flexppi_revision_20261002/corrected_training_v1')
STATIC=Path('/root/flexppi_revision_20261002/deployment_cost_v1')
PYTHON=ROOT/'env/bin/python'
MDPY=ROOT/'md_conda/bin/python'
MD=ROOT/'kix_md_triplicate_rebuild_20261010'
SINGLE=ROOT/'single_target'
FAMILY=ROOT/'family_folds'
state_path=ROOT/'pipeline_state.json'
control_path=ROOT/'pipeline_control.json'
ROOT.joinpath('pipeline_logs').mkdir(exist_ok=True)
now=lambda:time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(value):
    temporary=state_path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2),encoding='utf-8')
    temporary.replace(state_path)


def controls():
    return json.loads(control_path.read_text()) if control_path.exists() else {}


def process_identity(pid):
    folder=Path('/proc')/str(pid)
    try:
        raw=(folder/'stat').read_text()
        fields=raw[raw.rfind(')')+2:].split()
        return {'pid':int(pid),'start_ticks':int(fields[19]),'state':fields[0],
                'command':(folder/'cmdline').read_bytes().decode().rstrip('\0').split('\0')}
    except FileNotFoundError:return None


def preserved_process_is_alive(identity):
    actual=process_identity(identity['pid'])
    return actual is not None and actual['state']!='Z' and all(actual[key]==identity[key] for key in ['pid','start_ticks','command'])


def certify_preserved_single_queue():
    # The surviving original child commits this only after all ten original
    # strict config/epoch checks. Recheck the actual checkpoints on CPU.
    queued=json.loads((SINGLE/'single_target_queue_state.json').read_text())
    plan=json.loads((SINGLE/'single_target_plan.json').read_text())
    if queued.get('status')!='complete' or len(queued['jobs'])!=10:raise RuntimeError('Preserved single-target queue did not commit all ten jobs')
    lookup={(row['arm'],row['fold']):row for row in queued['jobs']}
    if set(lookup)!={(row['arm'],row['fold']) for row in plan['jobs']}:raise RuntimeError('Preserved queue job identities differ')
    sys.path.insert(0,str(SINGLE))
    from run_single_target_queue import valid_final
    for job in plan['jobs']:
        row=lookup[(job['arm'],job['fold'])]
        if not row['status'].startswith('complete'):raise RuntimeError('Preserved job incomplete')
        actual=valid_final(Path(job['final_checkpoint']),job)
        if actual!=row['checkpoint_sha256']:raise RuntimeError('Preserved final checkpoint hash changed')
    return {'status':'complete','jobs':10,'queue_state_sha256':sha(SINGLE/'single_target_queue_state.json'),
            'all_epoch040_strictly_validated':True,'no_training_restarted':True}


def run(name,command,check,cwd=ROOT):
    previous=next((item for item in state['steps'] if item['name']==name),None)
    if previous and previous['status']=='complete':
        if not all(Path(path).exists() for path in check):raise RuntimeError('Missing previously committed output: '+name)
        return
    control=controls()
    if name in control.get('cancel_steps',[]):
        item={'name':name,'status':'cancelled_per_user','cancelled_utc':now(),'never_started':True,
              'reason':control.get('reason'),'control_file':str(control_path),'control_sha256':sha(control_path),
              'command_not_executed':[str(part) for part in command]}
        if previous:state['steps'].remove(previous)
        state['steps'].append(item);state.update(updated_utc=now());save(state)
        print(now(),name,'cancelled_per_user; not launched',flush=True)
        return
    handoff=control.get('preserved_running_steps',{}).get(name)
    if handoff and not (previous and previous['status']=='complete'):
        if handoff['command']!=[str(part) for part in command]:raise RuntimeError('Preserved process command differs: '+name)
        item=dict(previous or {},name=name,status='running',command=handoff['command'],
                  completion_files=[str(path) for path in check],preserved_process=handoff,
                  resumed_monitor_utc=now(),parent_pid=os.getpid(),training_not_restarted=True)
        if previous:state['steps'].remove(previous)
        state['steps'].append(item);state.update(status='running',current_step=name,updated_utc=now());save(state)
        print(now(),name,'waiting for preserved child',handoff['pid'],flush=True)
        while preserved_process_is_alive(handoff):
            state.update(updated_utc=now(),runtime_control_sha256=sha(control_path));save(state);time.sleep(20)
        if not all(Path(path).exists() for path in check):raise RuntimeError('Preserved child exited without final outputs: '+name)
        if name!='train_single_targets':raise RuntimeError('Unsupported preserved completion gate: '+name)
        item.update(status='complete',finished_utc=now(),exit_code=None,
                    completion_method='Original surviving child committed all ten jobs; final checkpoint configs/hashes independently rechecked',
                    preserved_completion_audit=certify_preserved_single_queue(),
                    completion_sha256={str(path):sha(path) for path in check if Path(path).is_file()})
        state['updated_utc']=now();save(state)
        return
    item={'name':name,'status':'running','started_utc':now(),'command':[str(part) for part in command],
          'completion_files':[str(path) for path in check]}
    if previous:state['steps'].remove(previous)
    state['steps'].append(item);state.update(status='running',current_step=name,updated_utc=now());save(state)
    log=ROOT/'pipeline_logs'/f'{name}.log'
    env=dict(os.environ,PYTHONPATH=str(TRAIN/'src'),OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',CUDA_VISIBLE_DEVICES='0',PYTHONUNBUFFERED='1')
    if str(command[0])==str(MDPY):
        # Conda's CUDA plugin RPATH selects NVRTC 12.6. The host's driver
        # supports 12.4 PTX; an explicit 12.4 preload passed all four real smokes.
        nvrtc=ROOT/'env/lib/python3.11/site-packages/nvidia/cuda_nvrtc/lib'
        env['LD_PRELOAD']=str(nvrtc/'libnvrtc.so.12')
        env['LD_LIBRARY_PATH']=str(nvrtc)+':'+str(ROOT/'md_conda/lib')
    print(now(),name,flush=True)
    with log.open('ab') as handle:
        process=subprocess.run([str(part) for part in command],cwd=cwd,env=env,stdout=handle,stderr=subprocess.STDOUT)
    item.update(finished_utc=now(),exit_code=process.returncode,log=str(log))
    if process.returncode or not all(Path(path).exists() for path in check):
        item['status']='failed';state.update(status='failed',updated_utc=now());save(state)
        raise RuntimeError('Completion step failed: '+name+'; inspect '+str(log))
    item.update(status='complete',completion_sha256={str(path):sha(path) for path in check if Path(path).is_file()})
    state['updated_utc']=now();save(state)


with (ROOT/'gpu_queue.lock').open('a+') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    state=json.loads(state_path.read_text()) if state_path.exists() else {'started_utc':now(),'steps':[],'status':'initializing','pid':os.getpid()}
    state['pid']=os.getpid()
    control=controls()
    state.update(runtime_control_file=str(control_path),runtime_control_sha256=sha(control_path) if control_path.exists() else None,
                 md_reuse_policy=control.get('md_reuse_policy'),cancelled_pending_steps=control.get('cancel_steps',[]))
    save(state)
    # Installer processes are independent CPU/network work. Neither ready flag
    # nor mere import success is considered a scientific experiment result.
    try:
        deadline=time.monotonic()+7200
        while not (ROOT/'ENV_READY').exists() or not (ROOT/'MD_ENV_READY').exists():
            if time.monotonic()>deadline:raise RuntimeError('Environment initialization did not finish within two hours')
            state.update(status='waiting_for_environments',updated_utc=now());save(state);time.sleep(20)
        run('freeze_kix_predictions',[PYTHON,'-u',ROOT/'freeze_kix_predictions.py','--training-root',TRAIN,'--static-root',STATIC,'--md-root',MD],
            [MD/'frozen_predictions/FROZEN.json',MD/'frozen_predictions/residue_predictions.csv'])
        run('verify_prior_kix_predictions',[PYTHON,'-u',ROOT/'verify_kix_predictions.py','--md-root',MD,'--prior-csv',ROOT/'kix_fold1_all_residues.csv'],
            [MD/'frozen_predictions/prior_comparison_audit.json'])
        run('prepare_md',[MDPY,'-u',MD/'code/run_md.py','prepare','--root',MD],[MD/'prepared/READY.json'])
        run('audit_md_inputs',[MDPY,'-u',MD/'code/audit_rebuilt_prepared.py','--root',MD],[MD/'prepared/input_audit.json'])
        run('smoke_md_cuda',[MDPY,'-u',MD/'code/run_md.py','smoke','--root',MD],[MD/'smoke/PASSED.json'])
        run('verify_md_queue',[MDPY,'-u',MD/'code/queue_triplicate.py','--root',MD,'--verify-only'],[])
        run('smoke_single_targets',[PYTHON,'-u',SINGLE/'smoke_single_target.py','--plan',SINGLE/'single_target_plan.json','--device','cuda'],[SINGLE/'smoke_audit.json'])
        run('train_single_targets',[PYTHON,'-u',SINGLE/'run_single_target_queue.py','--plan',SINGLE/'single_target_plan.json'],
            [SINGLE/f'runs/{arm}/fold_{fold}/epoch_040.pt' for fold in range(5) for arm in ['rmsf_only','chi_entropy_only']])
        run('export_single_targets',[PYTHON,'-u',SINGLE/'export_completion_oof.py','--plan',SINGLE/'single_target_plan.json'],[SINGLE/'evaluation/evaluation_audit.json'])
        run('analyze_single_targets',[PYTHON,'-u',SINGLE/'analyze_completion_oof.py','--dual-evaluation',TRAIN/'evaluation',
            '--single-evaluation',SINGLE/'evaluation','--output',SINGLE/'analysis'],[SINGLE/'analysis/analysis_audit.json',SINGLE/'analysis/single_vs_dual_per_record.csv'])
        run('train_evaluate_family_folds',[PYTHON,'-u',FAMILY/'run_family_queue.py','--source-root',TRAIN/'src','--evaluate'],[FAMILY/'evaluation/evaluation_audit.json'])
        run('evaluate_external_2430',[PYTHON,'-u',FAMILY/'evaluate_external_2430.py','--source-root',TRAIN/'src',
            '--external-root',TRAIN/'external','--static-output',ROOT/'baselines/static_external_2430',
            '--static-reference',ROOT/'baselines/static_reference'],[FAMILY/'external_2430_evaluation/evaluation_audit.json'])
        additions=ROOT/'baselines/gpu_additional_steps.json'
        if not additions.exists():raise RuntimeError('Additional baseline protocol not frozen yet: '+str(additions))
        for job in json.loads(additions.read_text()):
            run(job['name'],job['command'],job['check'])
        run('md_triplicate_production',[MDPY,'-u',MD/'code/queue_triplicate.py','--root',MD],
            [MD/f'runs/{state}/rep{repeat:02d}/analysis/summary.json' for repeat in [1,2,3] for state in ['apo','myb','mll','both']])
        run('analyze_md_triplicates',[MDPY,'-u',MD/'code/analyze_triplicate_delta.py','--root',MD],[MD/'triplicate_delta_analysis/summary.json'])
        for repeat in ['rep01','rep02','rep03']:
            run('quality_md_'+repeat,[MDPY,'-u',MD/'code/quality_one_repeat.py','--root',MD,'--repeat',repeat],
                [MD/'quality'/repeat/'summary.json'])
        state.update(status='computations_complete_pending_manuscript_update',current_step=None,finished_utc=now());save(state)
    except Exception as error:
        state.update(status='failed',error=repr(error),updated_utc=now());save(state)
        raise
