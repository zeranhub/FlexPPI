"""Read-only snapshot of new local execution, partial smoke, and frozen OOF provenance.

This is a completion/timing inventory, never a metric leaderboard. Read status
files only; model output identity verification belongs to the scoring pipeline.
"""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BENCH = ROOT.parent
EXPECTED_RECORDS = 3063
EXPECTED_CHAINS = 4973
EXPECTED_RESIDUES = 961321


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_file(path, default=None):
    return json.loads(Path(path).read_text(encoding='utf-8')) if Path(path).exists() else default


def status_records(directory, pattern='record_*.json'):
    return [json_file(p) for p in Path(directory).glob(pattern) if not p.name.endswith('.tmp.json')]


def csv_records(path):
    if not Path(path).exists():
        return []
    with Path(path).open(encoding='utf-8',newline='') as handle:
        return list(csv.DictReader(handle))


def main():
    generated = datetime.now(timezone.utc).isoformat()
    runtime_choice=json_file(BENCH/'long_queue_config.json',{})
    inventory=[]
    def add(method,provenance,context,units,complete=0,failed=0,residues=0,
            chain_inputs='',seconds='',threads='',timing_scope='',source='',note='',expected=EXPECTED_RECORDS):
        status='complete_local' if provenance=='new_local_cpu' and complete==expected else 'frozen_oof_reuse' if provenance=='frozen_oof' else 'partial_local' if complete else 'not_run_local'
        inventory.append(dict(method=method,execution_status=status,provenance=provenance,
            input_context=context,output_units=units,expected_records=expected,complete_records=complete,
            failed_records=failed,covered_residues=residues,expected_residues=EXPECTED_RESIDUES,
            complete_chain_inputs=chain_inputs,threads=threads,measured_seconds_sum=seconds,
            timing_scope=timing_scope,status_source=source,notes=note,snapshot_utc=generated))
    for name,context in [('FlexPPI Full','complete_structure'),('FlexPPI Intra','within_chain'),('Static','complete_structure'),('AA','residue_identity')]:
        add(name,'frozen_oof',context,'Angstrom',EXPECTED_RECORDS,0,EXPECTED_RESIDUES,
            source=str(BENCH/'cohort/existing_prediction_manifest.csv'),
            note='Reused frozen five-fold OOF predictions; input/output identity/hash inventory already audited. This run contains no new model training or inference timing.')
    dyna_summary=json_file(BENCH/'dynaprot/predictions/summary.json',{})
    statuses=[] if dyna_summary.get('complete')==EXPECTED_CHAINS else status_records(BENCH/'dynaprot/predictions','chain_*.json')
    complete_chains=sum(s.get('status')=='complete' for s in statuses)
    if dyna_summary.get('complete')==EXPECTED_CHAINS:
        complete_chains=EXPECTED_CHAINS
    add('DynaProt-M','new_local_cpu','independent_chain','Angstrom',
        EXPECTED_RECORDS if complete_chains==EXPECTED_CHAINS else 0,
        dyna_summary.get('failed',sum(s.get('status')=='failed' for s in statuses)),
        EXPECTED_RESIDUES if complete_chains==EXPECTED_CHAINS else sum(s.get('expected_residues',0) for s in statuses if s.get('status')=='complete'),
        chain_inputs=complete_chains,seconds=dyna_summary.get('wall_seconds',sum(s.get('wall_seconds',0) for s in statuses if s.get('status')=='complete')),threads=2,
        timing_scope='Per-chain parsing, preprocessing, forward and save; model load/download excluded. See summary wall_seconds separately.',
        source=str(BENCH/'dynaprot/predictions'),note='Chain outputs must all be present before a record is eligible; no multichain joint inference.')
    nma_statuses=status_records(ROOT/'nma_scipy/predictions')
    nma_v2=json_file(ROOT/'gnm_ckdtree_v2/overlay_manifest.json',{})
    for method in ['ANM','GNM']:
        statuses=nma_statuses
        key=method.lower()+'_seconds'
        complete=[s for s in statuses if s.get('status')=='complete']
        add(method,'new_local_cpu','complete_structure','relative',len(complete),sum(s.get('status')=='failed' for s in statuses),
            sum(s.get('residues',0) for s in complete),seconds=sum(s.get(key,0) for s in complete),threads=2,
            timing_scope='Solver only for this model; shared CA parse and Hessian/Kirchhoff construction omitted here. Pair total measured separately.',
            source=str(ROOT/'gnm_ckdtree_v2/resolved_prediction_index.csv') if method=='GNM' and nma_v2 else str(ROOT/'nma_scipy/predictions'),
            note=('NumPy/SciPy reproduction, cutoff15A/gamma1/20 nonzero modes; default CKDTree float32 geometry with one versioned correction. Nine retained degenerate-mode records are not pointwise equivalent to old cache; see final_reference_audit_v2.json. Not calibrated Angstroms.' if method=='GNM' and nma_v2 else
                  'NumPy/SciPy reproduction of ProDy dense equations, cutoff15A/gamma1/20 nonzero modes. Full961321-position ANM maximum difference1.4071e-7, no positions above1e-6; not calibrated Angstroms.'))
    for tag in ['backflip-2.1','backflip-1.0','backflip-1.0-seq','backflip-1.0-aa']:
        for context in ['whole_complex','independent_chain']:
            directory=BENCH/'backflip/results'/tag
            run_source='full_queue'
            if not directory.exists():
                directory=BENCH/'backflip/smoke_results'/tag
                run_source='smoke_only'
            if not directory.exists() and tag=='backflip-1.0-aa':
                continue
            latest={}
            if (directory/'status.jsonl').exists():
                for line in (directory/'status.jsonl').read_text(encoding='utf-8').splitlines():
                    if line.strip():
                        item=json.loads(line)
                        if item['context']==context:
                            latest[item['sample_id']]=item
            complete=[s for s in latest.values() if s.get('status')=='complete']
            failed=[s for s in latest.values() if s.get('status')=='failed']
            add(tag+' ('+context+')','new_local_cpu',context,'Angstrom',len(complete),len(failed),
                sum(s.get('predicted',0) for s in complete),seconds=sum(s.get('seconds',0) for s in complete),threads=2,
                timing_scope='Canonical PDB parse/frame conversion/official forward/output conversion; excludes model import/load/input preparation.',
                source=str(directory),note=run_source+'; independent-chain single-chain records share their own fresh assembly forward and contribute zero additional seconds, not zero model cost.')
    fast_path=BENCH/'fastprotflex/full_results.csv'
    fast=csv_records(fast_path)
    fast_complete=[s for s in fast if s.get('status')=='complete']
    add('FastProtFlex','new_local_cpu','complete_structure','relative',len(fast_complete),sum(str(s.get('status','')).startswith('failed') for s in fast),
        sum(int(s.get('ca_positions_verified',0)) for s in fast_complete),seconds=sum(float(s.get('wall_seconds',0)) for s in fast_complete),threads=1,
        timing_scope='Per-record R subprocess including launch, graph construction, official function, save; memory waits are separate.',
        source=str(fast_path),note='One-thread CPU affinity; this runtime is not comparable to two-thread deep model forward-only measurements. Full-precision values, not rounded PDB B factors.')
    statuses=status_records(ROOT/'rmsfnet/smoke')
    complete=[s for s in statuses if s.get('status')=='complete']
    add('RMSF-net PDB-only','new_local_cpu','complete_structure','relative',len(complete),sum(s.get('status')=='failed' for s in statuses),
        sum(s.get('residues',0) for s in complete),seconds=sum(s.get('total_seconds',0) for s in complete),threads=2,
        timing_scope='Three-source smoke, full-assembly Gaussian density + official forward; batch2. No full queue was started.',
        source=str(ROOT/'rmsfnet/smoke'),note='Only 3/3063 records, 1433/961321 residues. Full queue exact153284 blocks; original batch2 estimate31.77h, not measured completion. Layout microtest does not establish full equivalence.')
    for name,context,units,reason in [
        ('PLANET-MD full','independent_chain','Angstrom (nm x10)','Historical per-record results exist, but this inventory counts only new local inference and none has been launched here.'),
        ('PLANET-MD sequence','sequence_chain','Angstrom (nm x10)','Historical per-record results exist, no new local whole-cohort per-residue output yet.'),
        ('PLANET-MD mini','sequence_chain','Angstrom (nm x10)','Sequence model with shared ESM3 backbone; checkpoint preparation only, no new local model forward yet.'),
        ('Flexpert-3D','independent_chain','Angstrom (nm x10)','Deployment incomplete: partial weight is not loadable. Official range probes remain slow; no active full download or local forward. See references/flexpert_weight_recovery_audit_20261024.json.'),
        ('Flexpert-Seq','sequence_chain','Angstrom (nm x10)','Deployment incomplete: partial weight is not loadable. Official range probes remain slow; no active full download or local forward. See references/flexpert_weight_recovery_audit_20261024.json.'),
        ('Pegasus','sequence_chain','method-dependent','Resource assessment only; official large protein language-model inputs have not run on the shared 16GB CPU machine.')]:
        if name.startswith('PLANET-MD') and not runtime_choice.get('run_planet',False):
            reason+=' Full inference deferred by user: complete quicker models first.'
        add(name,'candidate_not_run',context,units,note=reason)
    output=ROOT/'execution_inventory'
    output.mkdir(exist_ok=True)
    nma_complete=[s for s in nma_statuses if s.get('status')=='complete']
    nma_total=sum(s.get('total_seconds',0) for s in nma_complete)
    nma_anm=sum(s.get('anm_seconds',0) for s in nma_complete)
    nma_gnm=sum(s.get('gnm_seconds',0) for s in nma_complete)
    (output/'nma_timing_breakdown.json').write_text(json.dumps(dict(
        complete_records=len(nma_complete),expected_records=EXPECTED_RECORDS,
        complete_residues=sum(s.get('residues',0) for s in nma_complete),
        cpu_threads=2,pair_total_seconds_sum=nma_total,
        anm_solver_seconds_sum=nma_anm,gnm_solver_seconds_sum=nma_gnm,
        shared_parse_build_save_and_other_seconds_sum=nma_total-nma_anm-nma_gnm,
        gnm_v2_precision_scan_and_diagnostic_seconds=nma_v2.get('wall_seconds'),
        caution='ANM and GNM share this matrix-construction pass; solver columns alone omit shared work. Different thread counts and timing scopes prevent cross-model speed ranking.',
        snapshot_utc=generated),indent=2),encoding='utf-8')
    for filename,rows in [('model_execution_inventory.csv',inventory),('timing_table.csv',[r for r in inventory if r['measured_seconds_sum']!=''])]:
        with (output/filename).open('w',encoding='utf-8-sig',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
    (output/'inventory_audit.json').write_text(json.dumps(dict(snapshot_utc=generated,script_sha256=sha(__file__),
        user_runtime_config=runtime_choice,
        records=EXPECTED_RECORDS,chains=EXPECTED_CHAINS,residues=EXPECTED_RESIDUES,
        purpose='Operational completion/timing snapshot, not performance scoring or output-identity verification',
        no_unrun_models_as_complete=True,no_partial_smoke_as_full_cohort=True,
        timing_caution='Different thread counts, input contexts, and timing scopes prevent interpreting these sums as a speed ranking. Running queues create partial snapshots.'),indent=2),encoding='utf-8')
    print(json.dumps(inventory,indent=2),flush=True)


if __name__=='__main__':
    main()
