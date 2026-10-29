"""Bind recovered frozen split metadata; no torch, scoring or sequence search."""
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
import csv,hashlib,json,sys
sys.dont_write_bytecode=True
from audit_local_thermopcd import topology

ROOT=Path(r'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023')
CORE=Path(r'E:/Project/FlexPPI/experiments/core_corrected_pipeline')
OLD=Path(r'C:/Users/DELL/.codex/worktrees/4666/FlexPPI/paper/revision_20261002/experiments/corrected_pipeline')
OUT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):return [json.loads(s) for s in Path(p).read_text(encoding='utf-8-sig').splitlines() if s.strip()]
def csvrows(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def main():
    audit=read(CORE/'neural_exports/strict_external_evaluation/evaluation_audit.json')
    freeze_path=OLD/'strict_external_payload/strict_external_freeze.json'
    freeze=read(freeze_path)
    trainpath=OLD/'strict_external_payload/splits/strict_external_train.jsonl'
    assert sha(freeze_path)==audit['training_freeze_sha256']
    assert sha(trainpath)==freeze['execution_hashes']['splits\\strict_external_train.jsonl']
    train=rows(trainpath); trainids={r['sample_id'] for r in train}
    primary=rows(CORE/'training_payload/splits/all_ready.jsonl');primaryids={r['sample_id'] for r in primary}
    assert len(train)==len(trainids)==2431 and trainids<=primaryids and len(primaryids)==3063
    inventory_path=OLD/'sequence_audit/chain_inventory.json';chains=read(inventory_path)
    exact_primary=defaultdict(set);exact_strict=defaultdict(set)
    for r in chains:
        if r['role']=='training':
            exact_primary[r['sequence']].add(r['sample_id'])
            if r['sample_id'] in trainids:exact_strict[r['sequence']].add(r['sample_id'])
    external=OLD/'external_payload/external'
    assert sha(external/'external_payload_freeze.json')==audit['external_payload_sha256']
    assert sha(external/'prespecified_protocol.json')==audit['prespecified_protocol_sha256']
    ready=rows(external/'ready_test.jsonl');frozen=read(external/'external_payload_freeze.json')
    manifests={r['sample_id']:r for r in rows(CORE/'neural_exports/strict_external_evaluation/prediction_manifest.jsonl') if r['variant']=='full'}
    prior=[]
    for r in ready:
        sid=r['sample_id'];pdb=external/'structures'/(sid+'.pdb');label=external/'labels'/(sid+'.npz')
        assert sha(pdb)==r['metadata']['input_sha256']==manifests[sid]['input_sha256']==frozen['data_sha256']['structures\\'+sid+'.pdb']
        assert sha(label)==r['metadata']['label_sha256']==manifests[sid]['label_sha256']==frozen['data_sha256']['labels\\'+sid+'.npz']
        t=topology(pdb)
        prior.append(dict(sample_id=sid,source_pdb_sha256=sha(pdb),source_label_sha256=sha(label),actual_chain_count=t['chain_count'],chain_lengths=t['chain_lengths'],same_sequence_multichain=t['exact_homomer_topology_candidate'],strict_training_sample_id_overlap=sid in trainids,
            exact_primary_training_records=sorted(set().union(*(exact_primary[s] for s in t['chain_sequences'].values()))),
            exact_strict_training_records=sorted(set().union(*(exact_strict[s] for s in t['chain_sequences'].values())))))
    production=[]
    alltop={r['path']:r for r in csvrows(OUT/'thermopcd_topology_inventory.csv')}
    for r in csvrows(OUT/'thermopcd_directory_readiness.csv'):
        if r['production_candidate']!='True':continue
        t=alltop.get(r['selected_initial_topology'],{})
        seqs=json.loads(t.get('chain_sequences','{}'))
        production.append(dict(directory=r['directory'],pdb_id=r['pdb_id'],previous_strict32_pdb=r['previous_strict32_pdb'],
            selected_topology=r['selected_initial_topology'],chain_count=r['selected_topology_chain_count'],classification=r['selected_topology_classification'],
            exact_primary_training_records=sorted(set().union(*(exact_primary[s] for s in seqs.values()))),
            exact_strict_training_records=sorted(set().union(*(exact_strict[s] for s in seqs.values()))),
            new_candidate_homology_status='not_newly_searched; original32 validated under frozen MMseqs rule, new sequences need declared search before family-heldout claim'))
    membership=csvrows(ROOT/'protein_complex_benchmark/classification/membership_frozen.csv')
    original=csvrows(ROOT/'cohort/records.csv');origdict={r['sample_id']:r for r in original}
    primary_multi=[r['sample_id'] for r in membership if r['benchmark_group']=='protein' and int(origdict[r['sample_id']]['n_chains'])>=2]
    strict_overlap=read(OLD/'sequence_audit/strict_overlap_validation.json')
    assert strict_overlap['strict_training_records']==2431 and strict_overlap['maximum_external_chain_hits_after_purge']==0
    output=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='PASS_RECOVERED_FROZEN_STRICT_SOURCES',
        recovered_base=str(OLD),training_records=2431,training_source_counts=dict(Counter(r['dataset'] for r in train)),primary_training_universe=3063,
        existing_external_records=len(ready),existing_external_original_input_and_label_hashes_passed=len(prior),
        existing_external_actual_chain_count_distribution=dict(Counter(r['actual_chain_count'] for r in prior)),existing_external_same_sequence_multichain=sum(r['same_sequence_multichain'] for r in prior),
        existing_external_sample_id_overlap_with_primary=sorted({r['sample_id'] for r in ready}&primaryids),
        existing_external_sample_id_overlap_with_strict=sorted({r['sample_id'] for r in ready}&trainids),
        existing_external_exact_chain_overlap_strict_records=sum(bool(r['exact_strict_training_records']) for r in prior),
        existing_frozen_homology_validation=strict_overlap['scope'],existing_frozen_homology_external_chains=strict_overlap['external_chain_count'],
        existing_primary10_multichain_protein_ids=primary_multi,primary10_still_in_strict_training=sorted(set(primary_multi)&trainids),
        primary10_not_in_strict_training=sorted(set(primary_multi)-trainids),
        strict_checkpoint_selection='A distinct predeclared purged experiment; do not relabel as arbitrary original OOF fold. Original3063 use corresponding frozen OOF fold. New test records require explicit overlap audit of fixed checkpoint training membership before inference.',
        new_external_test_selection='For candidates demonstrably absent from the complete3063 training universe, a predeclared equal five-fold Full ensemble is an available source-out protocol, with family exposure separately audited. A single strict_full alternative must be declared separately and checked against its2431 records; never choose weights by test score.',
        local_thermopcd_production_directory_count=len(production),no_inference=True,no_new_sequence_search=True,no_scores_read=True)
    for name,rs in [('recovered_strict32_source_binding.csv',prior),('thermopcd_production_exact_training_exposure.csv',production)]:
        with (OUT/name).open('x',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader()
            for r in rs:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})
    files=[freeze_path,trainpath,inventory_path,OLD/'sequence_audit/sequence_audit_summary.json',OLD/'sequence_audit/strict_overlap_validation.json',external/'ready_test.jsonl',external/'external_payload_freeze.json',external/'prespecified_protocol.json',Path(__file__),OUT/'readiness_audit.json']
    output['source_sha256']={str(p):sha(p) for p in files}
    with (OUT/'recovered_strict_source_audit.json').open('x',encoding='utf-8') as f:json.dump(output,f,indent=2,ensure_ascii=False)
    print(json.dumps({k:v for k,v in output.items() if k!='source_sha256'},ensure_ascii=False),flush=True)

if __name__=='__main__':main()
