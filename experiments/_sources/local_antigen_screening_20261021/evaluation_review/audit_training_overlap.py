"""Local read-only sequence/provenance audit; never loads or runs a model."""
from pathlib import Path
from collections import defaultdict, Counter
import csv, hashlib, json, sys

ROOT = Path(__file__).resolve().parent
EXP = ROOT.parents[1]
SOURCE = ROOT.parent / 'experimental_labels/raw/github/mutation_effects'
REF = Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
INVENTORY = REF / 'sequence_audit/chain_inventory.json'
TRAIN = EXP / 'family_folds/splits/fold_0_train.jsonl'
TEST = EXP / 'family_folds/splits/fold_0_test.jsonl'
FREEZE = EXP / 'family_folds/training_freeze.json'
GROUPS = EXP / 'family_folds/family_components.json'
AUDIT = EXP / 'itrace_case_search_only_20261017/model_audit/remote_checkpoint_manifest_audit.json'
seq_inventory_path = EXP / 'tcr_case_v40_20261017/local_pair_screen/sequence_inventory.json'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
read = lambda p: json.loads(Path(p).read_text(encoding='utf-8-sig'))
write = lambda p,x: p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
records = [json.loads(l) for l in TRAIN.read_text().splitlines() if l.strip()]
tests = [json.loads(l) for l in TEST.read_text().splitlines() if l.strip()]
inventory = read(INVENTORY)
freeze = read(FREEZE)
assert sha(TRAIN) == freeze['files']['splits/fold_0_train.jsonl']
assert sha(INVENTORY) == freeze['chain_inventory_sha256']
trainids = {r['sample_id'] for r in records}
testids = {r['sample_id'] for r in tests}
assert not trainids.intersection(testids)
trainchains = [c for c in inventory if c['sample_id'] in trainids]
seq_to_inventory = defaultdict(list)
trainseqs = defaultdict(list)
for c in inventory:
    assert c['sequence_id'] == 'seq_'+hashlib.sha256(c['sequence'].encode()).hexdigest()
    seq_to_inventory[c['sequence_id']].append(c)
for c in trainchains:
    trainseqs[c['sequence']].append(c)
queries = read(seq_inventory_path)
pdblist = ['2bnq','2bnr','1ao7','1qrn','1qse','1qsf']
queryrows = []
queryids = set()
for pdb in pdblist:
    q = queries[pdb]
    for role,c in q['roles'].items():
        seq = c['sequence']
        sid = 'seq_'+hashlib.sha256(seq.encode()).hexdigest()
        queryids.add(sid)
        queryrows.append(dict(pdb=pdb, role=role, chain=c['chain'], length=len(seq), sequence=seq, sequence_id=sid,
            inventory_sequence_present=sid in seq_to_inventory,
            exact_training_chains=[{k:x[k] for k in ['sample_id','chain','dataset','biological_system','length']} for x in trainseqs.get(seq,[])],
            structure=q['structure'],structure_sha256=sha(q['structure'])))

# Re-audit the original recorded MMseqs alignments. Integer span coverages
# deliberately reproduce prepare_family_folds.py, including its 3-decimal identity.
paths = [REF/'sequence_audit/all_chain_similarity.tsv',
         EXP/'all_chain_unmasked_similarity.tsv',
         EXP/'family_folds/all_chain_unmasked_similarity.tsv']
edges = defaultdict(dict)
query_search_seen = set()
searchfiles = []
for path in paths:
    searchfiles.append(dict(path=str(path),sha256=sha(path)))
    with path.open() as f:
        for line in f:
            a = line.rstrip().split('\t')
            if a[0] not in queryids and a[1] not in queryids: continue
            qid,tid=a[:2]
            if qid in queryids: query_search_seen.add(qid)
            if tid in queryids: query_search_seen.add(tid)
            if qid==tid: continue
            identity=float(a[2]); qc=(int(a[5])-int(a[4])+1)/int(a[6]);tc=(int(a[8])-int(a[7])+1)/int(a[9])
            if identity < .30 or min(qc,tc) < .70: continue
            for source_id,other_id,covq,covt in [(qid,tid,qc,tc),(tid,qid,tc,qc)]:
                if source_id not in queryids: continue
                for hit in seq_to_inventory[other_id]:
                    if hit['sample_id'] not in trainids: continue
                    key=(hit['sample_id'],hit['chain'])
                    edges[source_id][key]=dict(sample_id=hit['sample_id'],chain=hit['chain'],biological_system=hit['biological_system'],identity=identity,query_coverage=covq,target_coverage=covt)

groups=read(GROUPS)
recordgroup = {s:g for g in groups for s in g['members']}
record_summaries=[]
for pdb in pdblist:
    sid='itrace_'+pdb
    group=recordgroup[sid]
    assert sid in testids and sid not in trainids
    assert not trainids.intersection(group['members'])
    record_summaries.append(dict(pdb=pdb,sample_id=sid,in_train=False,in_fold0_test=True,
        component_id=group['component_id'],component_fold=group['fold'],component_training_pool_records=group['training_records'],
        component_external_records=group['external_records'],component_members_present_in_actual_fold0_train=0,
        exact_pdb_train_hits=[r['sample_id'] for r in records if pdb.upper() in str(r.get('metadata',{}).get('pdb_ids',[])).upper() or r.get('metadata',{}).get('pdb_id','').upper()==pdb.upper()]))
for q in queryrows:
    q['recorded_mmseqs_query_seen']=q['sequence_id'] in query_search_seen
    q['recorded_mmseqs_30pct_70pct_both_training_hits']=list(edges.get(q['sequence_id'],{}).values())

# Local candidate peptide exact-identity checks; do not read affinity values.
candidate_audit=[]
for system,reference in [('nyeso','SLLMWITQC'),('tax','LLFGYPVYV')]:
    path=SOURCE/system/f'{system}_peptide_kd.csv'
    sequences=sorted({r['sequence'].strip() for r in csv.DictReader(path.open(encoding='utf-8-sig'))})
    for seq in sequences:
        positions=[i+1 for i,(a,b) in enumerate(zip(reference,seq)) if a!=b]
        exact=trainseqs.get(seq,[])
        candidate_audit.append(dict(system=system,reference_peptide=reference,sequence=seq,length=len(seq),
            same_length_as_reference=len(seq)==len(reference),mutation_positions=positions,
            training_exact_peptide_hits=[{k:h[k] for k in ['sample_id','chain','dataset','biological_system']} for h in exact]))
write(ROOT/'candidate_sequence_overlap.json',candidate_audit)

# Supplemental local alignment is an independent descriptive nearest-neighbour
# screen. It is not falsely relabelled as the original MMseqs family search.
sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
from Bio.Align import PairwiseAligner
aligner=PairwiseAligner(mode='local',match_score=2,mismatch_score=-1,open_gap_score=-3,extend_gap_score=-0.5)
nearest={}
unique_train={s:c for s,c in trainseqs.items()}
for seq in sorted({q['sequence'] for q in queryrows if q['length']>=20}):
    hits=[]
    for target,chainrecords in unique_train.items():
        # This prefilter cannot discard an alignment with both span coverages>=0.7.
        if min(len(seq),len(target))/max(len(seq),len(target))<.70: continue
        a=aligner.align(seq,target)[0]
        qb,tb=a.aligned
        paired=sum(int(e-s) for s,e in qb)
        if not paired: continue
        qspan=int(qb[-1,1]-qb[0,0]);tspan=int(tb[-1,1]-tb[0,0])
        qc=qspan/len(seq);tc=tspan/len(target)
        if min(qc,tc)<.7: continue
        identical=sum(sum(aa==bb for aa,bb in zip(seq[int(qs):int(qe)],target[int(ts):int(te)])) for (qs,qe),(ts,te) in zip(qb,tb))
        columns=max(qspan,tspan)  # also report stricter paired+gap column identity below
        qcoord,tcoord=a.coordinates
        ncolumns=sum(max(int(qcoord[i+1]-qcoord[i]),int(tcoord[i+1]-tcoord[i])) for i in range(len(qcoord)-1))
        ident=identical/ncolumns
        hit=dict(identity_including_gap_columns=ident,identity_paired_positions=identical/paired,query_span_coverage=qc,target_span_coverage=tc,
            score=float(a.score),alignment_columns=ncolumns,paired_residues=paired,identical_residues=identical,
            training_chain_records=[{k:c[k] for k in ['sample_id','chain','biological_system','length']} for c in chainrecords])
        hits.append(hit)
    hits.sort(key=lambda h:(h['identity_including_gap_columns'],h['score']),reverse=True)
    nearest['seq_'+hashlib.sha256(seq.encode()).hexdigest()]={
        'top10':hits[:10],
        'hits_identity_ge_0_30_count':sum(h['identity_including_gap_columns']>=.3 for h in hits),
        'hits_identity_ge_0_30':[h for h in hits if h['identity_including_gap_columns']>=.3],
    }
    print('audited query length',len(seq),'nearest',hits[0]['identity_including_gap_columns'] if hits else None,flush=True)
write(ROOT/'independent_local_alignment_neighbours.json',nearest)
audit=read(AUDIT)
checkpoint=audit['checkpoints'][0]
assert checkpoint['manifests']['state_manifest']['sha256']==sha(TRAIN)
summary=dict(status='complete_local_audit',no_inference_or_training_run=True,no_remote_access=True,
    candidate_numeric_affinities_not_read=True,
    checkpoint_provenance=dict(reported_checkpoint_sha256=checkpoint['sha256'],reported_checkpoint_path=checkpoint['path'],
        local_checkpoint_byte_hash_pending=True,sequence_backend=checkpoint['config']['model']['sequence_backend'],epoch=checkpoint['epoch'],seed=checkpoint['config']['seed']),
    training_manifest=dict(path=str(TRAIN),sha256=sha(TRAIN),records=len(records),sources=dict(Counter(r['dataset'] for r in records)),itrace_records=sum(r['dataset']=='itrace' for r in records),chain_count=len(trainchains),unique_chain_sequences=len(trainseqs)),
    source_hashes=[dict(path=str(p),sha256=sha(p)) for p in [INVENTORY,TEST,FREEZE,GROUPS,AUDIT,seq_inventory_path]]+searchfiles,
    record_membership=record_summaries,query_chain_audit=queryrows,
    exact_chain_overlap_count=sum(len(q['exact_training_chains']) for q in queryrows),
    accepted_recorded_mmseqs_overlap_count=sum(len(q['recorded_mmseqs_30pct_70pct_both_training_hits']) for q in queryrows),
    supplemental_alignment_methods=dict(mode='local',match=2,mismatch=-1,open_gap=-3,extend_gap=-.5,min_both_integer_span_coverage=.7,identity='identical / all alignment columns including gap columns',not_equivalent_to_original_mmseqs=True),
    candidate_unique_sequences=dict(Counter(r['system'] for r in candidate_audit)),
    candidate_exact_training_peptide_sequences=sum(bool(r['training_exact_peptide_hits']) for r in candidate_audit),
    conclusions_not_supported=['No-ITrace alone does not establish independence; actual training membership and related chains were audited here.','This audit excludes direct or recorded threshold-family overlap, not all remote structural or functional relationship.','Candidate sequence equality in an affinity panel is not itself a train/test partition; downstream calibration must remain grouped.','Reference PDBs and published SPR outcomes were available during earlier case selection; this is retrospective testing, not prospective blind design.'])
write(ROOT/'training_overlap_audit.json',summary)
print(json.dumps({k:summary[k] for k in ['status','training_manifest','exact_chain_overlap_count','accepted_recorded_mmseqs_overlap_count','candidate_unique_sequences','candidate_exact_training_peptide_sequences']},ensure_ascii=False))
