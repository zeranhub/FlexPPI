"""Build local reviewable CPU-only protocol; does not run any remote command."""
import hashlib
import json
from pathlib import Path
from check_additional_baseline_freeze import HASHES, COHORT_HASH


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    local = Path(__file__).resolve().parent
    root = '/root/flexppi_completion_20261010'
    base = root+'/baselines'
    python = root+'/cpu_env/bin/python'
    steps = {}
    steps['wait_resume_verify_rmsf_net'] = dict(
        command=[python,'-u',base+'/wait_resume_rmsf_net_cpu.py','--timeout-seconds','86400','--poll-seconds','20'],
        check=[base+'/rmsf_net/run_audit.json',base+'/rmsf_net_resume_state.json'])
    for variant in ['full','seq']:
        steps['wait_verify_planet_'+variant+'_cpu'] = dict(
            command=[python,'-u',base+'/wait_additional_cpu_queue.py','--state',base+'/planet_'+variant+'_cpu_queue_state.json',
                '--output',base+'/planet_'+variant,'--model','planet-'+variant,'--timeout-seconds','86400','--poll-seconds','20',
                '--wait-audit',base+'/planet_'+variant+'_cpu_wait_audit.json'],
            check=[base+'/planet_'+variant+'/run_audit.json',base+'/planet_'+variant+'_cpu_wait_audit.json'])
    common = [python,'-u',base+'/score_common_residues.py','--cohort',base+'/cohort',
        '--prody',base+'/prody','--rmsf-net',base+'/rmsf_net','--planet-full',base+'/planet_full','--planet-seq',base+'/planet_seq']
    old = '/root/flexppi_revision_20261002/corrected_training_v1/evaluation'
    for mode, evaluation in [('original',old),('family',root+'/family_folds/evaluation')]:
        output = base+'/additional_'+mode+'_oof'
        command = common+['--flexppi',evaluation+'/full','--intrachain',evaluation+'/intrachain']
        if mode == 'family':
            command += ['--family-folds',root+'/family_folds/record_family_folds.csv']
        command += ['--output',output,'--require-complete']
        steps['score_additional_'+mode+'_oof'] = dict(command=command,
            check=[output+'/scoring_audit.json',output+'/record_macro_metrics.csv'],requires_all_completion=True)
    files = dict(HASHES)
    for name in ['check_additional_baseline_freeze.py','wait_additional_cpu_queue.py','wait_resume_rmsf_net_cpu.py','execute_frozen_baseline_step.py']:
        files[name] = sha(local/name)
    contract = dict(status='frozen_cpu_execution_contract',cohort_freeze_sha256=COHORT_HASH,
        no_gpu_tasks=True,policy='Wait existing CPU queues; PLANET standalone verifier is owned by its queue. Never replace its run_audit certificate. RMSF waits original processes and then independently verifies cached outputs.',
        files=files,steps=steps)
    path = local/'additional_cpu_execution_contract.json'
    path.write_text(json.dumps(contract,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
    contract_hash = sha(path)
    dispatch = []
    for name, step in steps.items():
        dispatch.append(dict(name=name,command=[python,'-u',base+'/execute_frozen_baseline_step.py',
            '--contract',base+'/additional_cpu_execution_contract.json','--contract-sha256',contract_hash,'--step',name],
            check=step['check']+[base+'/dispatch_'+name+'_audit.json']))
    output = local/'gpu_additional_steps.json'
    output.write_text(json.dumps(dispatch,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(contract_sha256=contract_hash,dispatch_sha256=sha(output),steps=len(dispatch),gpu_tasks=0),indent=2))


if __name__ == '__main__':
    main()
