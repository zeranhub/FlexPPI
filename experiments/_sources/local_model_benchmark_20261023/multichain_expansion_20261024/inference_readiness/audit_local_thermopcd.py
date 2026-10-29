"""Read small local topology/protocol metadata only; no MD or model execution."""
from pathlib import Path
from collections import Counter, defaultdict, OrderedDict
from datetime import datetime, timezone
import csv
import hashlib
import json
import pickletools
import re
import zipfile
import numpy as np

ROOT=Path(r'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023')
CORE=Path(r'E:/Project/FlexPPI/experiments/core_corrected_pipeline')
THERMO=Path(r'D:/data/ThermoPCD_trajectories')
OUT=Path(__file__).resolve().parent
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
ALIASES={'HID':'HIS','HIE':'HIS','HIP':'HIS','HSD':'HIS','HSE':'HIS','HSP':'HIS','MSE':'MET','CYX':'CYS','CYM':'CYS','ASH':'ASP','GLH':'GLU','LYN':'LYS'}

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))

def rows(p):
    with Path(p).open(encoding='utf-8-sig') as f:
        return [json.loads(x) for x in f if x.strip()]

def topology(p):
    chains=OrderedDict(); seen=set(); unknown=Counter(); header=[]; ter=0
    with p.open(encoding='utf-8',errors='replace') as f:
        for line in f:
            if line.startswith('ENDMDL'):break
            if line.startswith(('HEADER','TITLE ','COMPND','SOURCE')):header.append(line.rstrip())
            if line.startswith('TER'):ter+=1
            if not line.startswith(('ATOM  ','HETATM')) or line[12:16].strip()!='CA' or line[16:17] not in (' ','A','1'):continue
            rn=line[17:20].strip(); rn=ALIASES.get(rn,rn)
            if rn not in AA:
                unknown[rn]+=1;continue
            chain=line[21:22]; key=(chain,line[22:27])
            if key in seen:continue
            seen.add(key); chains.setdefault(chain,[]).append(AA[rn])
    seqs={c:''.join(s) for c,s in chains.items()}
    n=len(seqs); exact=n>=2 and len(set(seqs.values()))==1
    return dict(chain_count=n,chain_ids=list(seqs),chain_lengths={k:len(v) for k,v in seqs.items()},
                chain_sequences=seqs,residues=sum(map(len,seqs.values())),distinct_sequences=len(set(seqs.values())),
                all_sequences_identical=exact,blank_chain_present=' ' in seqs,ter_count=ter,
                exact_homomer_topology_candidate=exact and ' ' not in seqs,
                classification=('exact_same_sequence_multichain_candidate' if exact and ' ' not in seqs else 'single_chain' if n==1 else 'multichain_distinct_sequences_identity_unresolved' if n>=2 else 'no_parseable_protein_CA'),
                unknown_ca_residue_names=dict(unknown),pdb_header='\n'.join(header))

def pdb_id(text,names):
    patterns=[r'(?i)\bPDB\s*[-_:,]?\s*([0-9][A-Z0-9]{3})\b',r'(?i)^(?:DB\s+)?([0-9][A-Z0-9]{3})(?:\b|[_ ,.-])']
    for pat in patterns:
        m=re.search(pat,text)
        if m:return m.group(1).upper()
    for n in names:
        m=re.search(r'(?i)^([0-9][A-Z0-9]{3})(?:\b|[-_.])',n)
        if m and not m.group(1).upper().endswith('NS'):return m.group(1).upper()
    return ''

def write_csv(p,rs):
    keys=list(dict.fromkeys(k for r in rs for k in r))
    with p.open('x',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for r in rs:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})

def main():
    primary=rows(CORE/'training_payload/splits/all_ready.jsonl')
    trained_pdb=defaultdict(set)
    for r in primary:
        for p in r.get('metadata',{}).get('pdb_ids',[]):trained_pdb[p.upper()].add(r['sample_id'])
        if r.get('metadata',{}).get('pdb_id'):trained_pdb[r['metadata']['pdb_id'].upper()].add(r['sample_id'])
    candidates=rows(Path(r'E:/Project/FlexPPI/artifacts/data/thermopcd_performance_first_candidates.jsonl'))
    expected={r['sample_id']:r for r in candidates}
    previous_ids={r['sample_id'] for r in rows(CORE/'neural_exports/strict_external_evaluation/prediction_manifest.jsonl')}
    previous_pdbs={expected[s]['metadata']['pdb_id'].upper() for s in previous_ids if s in expected}
    inventory=[]; topology_rows=[]; topology_cache={}
    for collection in ('data','ThermoPCD_labels'):
        for d in sorted((THERMO/collection).iterdir()):
            if not d.is_dir():continue
            files=sorted(p for p in d.iterdir() if p.is_file())
            prod=[p for p in files if p.stat().st_size>0 and (p.name.lower() in ('dynamic-nopbc.xtc','system_reduced.xtc') or re.search(r'(?i)(?:production|dynamic|md)[-_].*\.xtc$',p.name))]
            pdbs=[p for p in files if p.suffix.lower()=='.pdb' and p.stat().st_size>0]
            initial=sorted(pdbs,key=lambda p:(0 if p.name.lower().endswith('-m.pdb') else 1 if '0ns' in p.name.lower() else 2,p.name))
            # Parse every nonzero topology; snapshot-only files cannot become production trajectories.
            prs=[]
            for p in initial:
                h=sha(p)
                if h not in topology_cache:topology_cache[h]=topology(p)
                pr=dict(path=str(p),sha256=h,**topology_cache[h]);prs.append(pr)
                topology_rows.append(dict(collection=collection,directory=str(d),pdb_id=pdb_id(d.name,[x.name for x in files]),**pr))
            selected=prs[0] if prs else {}
            pid=pdb_id(d.name,[p.name for p in files])
            inventory.append(dict(collection=collection,directory=str(d),pdb_id=pid,production_named_nonzero_files=[dict(path=str(p),bytes=p.stat().st_size) for p in prod],
                production_candidate=bool(prod),nonzero_pdb_count=len(pdbs),selected_initial_topology=selected.get('path',''),
                selected_topology_chain_count=selected.get('chain_count',0),selected_topology_classification=selected.get('classification','topology_missing'),
                exact_homomer_topology_candidate=selected.get('exact_homomer_topology_candidate',False),
                any_topology_exact_homomer_candidate=any(p['exact_homomer_topology_candidate'] for p in prs),
                selected_chain_lengths=selected.get('chain_lengths',{}),previous_strict32_pdb=pid in previous_pdbs,
                same_pdb_in_primary3063=sorted(trained_pdb.get(pid,[])),
                label_extraction_readiness='candidate_requires_uniform_extraction_and_topology_trajectory_order_QC' if prod and prs else 'no_continuous_MD_or_initial_topology_confirmed'))
    prior=[]
    for sid in sorted(previous_ids):
        p=CORE/'neural_exports/strict_external_evaluation/strict_full/fold_0'/(sid+'.npz')
        with np.load(p,allow_pickle=False) as z:
            chain=z['chain']; amino=z['amino_acid']; seqs=OrderedDict()
            for c,a in zip(chain,amino):seqs.setdefault(str(c),[]).append(str(a))
            seqs={c:''.join(s) for c,s in seqs.items()}
        r=expected.get(sid,{})
        source=Path(r.get('structure_path',''))
        prior.append(dict(sample_id=sid,pdb_id=r.get('metadata',{}).get('pdb_id',''),chain_count=len(seqs),
            chain_lengths={k:len(v) for k,v in seqs.items()},distinct_sequences=len(set(seqs.values())),
            exact_homomer_topology_candidate=len(seqs)>=2 and len(set(seqs.values()))==1,
            source_title=source.parent.name,original_topology_exists=source.is_file(),
            original_topology_path=str(source),strict_prediction_identity_source=str(p),
            original_primary_same_pdb_records=sorted(trained_pdb.get(r.get('metadata',{}).get('pdb_id','').upper(),[]))))
    checkpoint=CORE/'trained_checkpoints/strict_full/fold_0/epoch_040.pt'
    collection_audit=read_json(CORE/'trained_checkpoints/completed_strict_members_audit.json')
    declared=next(r for r in collection_audit['checkpoints'] if r['member']=='strict_full')
    checkpoint_sha=sha(checkpoint)
    assert checkpoint_sha==declared['checkpoint_sha256']
    # Decode only primitive pickle opcode arguments; never unpickle or load tensor storage.
    with zipfile.ZipFile(checkpoint) as z:
        pickle_bytes=z.read(next(n for n in z.namelist() if n.endswith('/data.pkl')))
    ops=list(pickletools.genops(pickle_bytes));start=next(i for i,(_,a,_) in enumerate(ops) if a=='config')
    end=next(i for i,(_,a,_) in enumerate(ops[start:],start) if a=='train_loss')
    primitives=[dict(offset=p,opcode=o.name,value=(True if o.name=='NEWTRUE' else False if o.name=='NEWFALSE' else a)) for o,a,p in ops[start:end] if o.name in ('BINUNICODE','SHORT_BINUNICODE','BININT','BININT1','BININT2','BINFLOAT','NEWTRUE','NEWFALSE','NONE')]
    strict_source=CORE/'final_training_run_records_20261004/strict_external_freeze.json'
    missing=[str(CORE/p) for p in ('sequence_audit/chain_inventory.json','sequence_audit/purged_training_local.jsonl','sequence_audit/strict_overlap_validation.json','strict_external_payload/splits/strict_external_train.jsonl','external_payload/external/ready_test.jsonl') if not (CORE/p).exists()]
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='metadata_topology_readiness_only_no_inference',
        scope='All nonzero local ThermoPCD PDB topologies; directory/file stat only for trajectories; existing strict32 NPZ identity arrays only; no prediction scores or tensor storage loaded.',
        directories=len(inventory),production_named_nonzero_directories=sum(r['production_candidate'] for r in inventory),
        topology_files=len(topology_rows),unique_topology_sha256=len(topology_cache),
        selected_topology_class_counts=dict(Counter(r['selected_topology_classification'] for r in inventory)),
        selected_exact_homomer_directories=sum(r['exact_homomer_topology_candidate'] for r in inventory),
        selected_exact_homomer_with_production_directories=sum(r['exact_homomer_topology_candidate'] and r['production_candidate'] for r in inventory),
        any_topology_exact_homomer_directories=sum(r['any_topology_exact_homomer_candidate'] for r in inventory),
        prior32=dict(records=len(prior),chain_count_distribution=dict(Counter(r['chain_count'] for r in prior)),exact_same_sequence_multichain=sum(r['exact_homomer_topology_candidate'] for r in prior)),
        strict_full=dict(checkpoint_path=str(checkpoint),checkpoint_sha256=checkpoint_sha,collection_manifest_sha_binding=True,primitive_config_opcodes=primitives,training_freeze=read_json(strict_source),
            applicability='The existing source/family-purged strict experiment may be used only after new-record membership/chain-homology checks against its2431 training records. It is not an alternate OOF fold for primary records.',
            missing_local_evidence=missing),
        no_inference=True,no_training=True,no_trajectory_content_read=True)
    for name,rs in [('thermopcd_directory_readiness.csv',inventory),('thermopcd_topology_inventory.csv',topology_rows),('existing_strict32_identity.csv',prior)]:write_csv(OUT/name,rs)
    sources=[Path(__file__),CORE/'training_payload/splits/all_ready.jsonl',Path(r'E:/Project/FlexPPI/artifacts/data/thermopcd_performance_first_candidates.jsonl'),THERMO/'manifest.jsonl',CORE/'neural_exports/strict_external_evaluation/evaluation_audit.json',CORE/'neural_exports/strict_external_evaluation/prediction_manifest.jsonl',CORE/'trained_checkpoints/completed_strict_members_audit.json',strict_source,checkpoint]
    report['source_sha256']=[dict(path=str(p),sha256=sha(p)) for p in sources]
    report['output_sha256']={p.name:sha(p) for p in OUT.glob('*.csv')}
    with (OUT/'readiness_audit.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2,ensure_ascii=False)
    print(json.dumps({k:v for k,v in report.items() if k not in ('strict_full','source_sha256','output_sha256')},ensure_ascii=False),flush=True)

if __name__=='__main__':main()
