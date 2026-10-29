"""Freeze a performance-blind protein / protein-complex membership audit.

Only cohort identity metadata, MD input sequences, and official RCSB annotations
are read. No model arrays, predictions, scores, or evaluation results are read.
MD chain letters are recorded but never used to select a source polymer entity.
"""
from pathlib import Path
import argparse, collections, concurrent.futures, datetime, hashlib, json, re, sys, time
import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE / 'runtime'))
from Bio.Align import PairwiseAligner

POLICY = {
    'version': 'protein_identity_v1_20261024',
    'human_confirmation_client_date': '2026-10-24',
    'human_confirmation': '同源多聚体归蛋白（推荐）',
    'groups': {'protein': 'One biological protein identity, including one chain, homomultimers, truncations, and split chains of one mature functional protein.',
               'protein_complex': 'Two or more different biological proteins represented in the actual MD input.'},
    'unknown': 'Evidence insufficient or ambiguous; excluded from frozen binary membership until explicitly resolved.',
    'single_input_chain': 'protein, as directly authorized; classification is about the simulated input, not all entities in the source PDB.',
    'tcr_pmhc_itrace': 'protein_complex, as directly authorized for all 244 source-confirmed records.',
    'matching': 'Match each actual MD chain sequence against every official canonical polymer sequence for its source PDB. Never infer correspondence from renamed author-chain letters.',
    'automatic_sequence_threshold': {'known_residue_match_fraction': 0.97, 'maximum_unknown_query_fraction': 0.05, 'short_chain_under_20_aa': 'exact substring only'},
    'identity_rule': 'An entity ID or a distinct sequence alone is not a distinct biological protein. Same UniProt is a candidate identity link, not proof when distinct mature products are annotated.',
    'processed_products': 'Split chains still forming one mature protein (for example insulin A/B) are protein; explicitly different mature products (for example HCV NS3 plus NS4A) are protein_complex even if their precursor UniProt accession is the same. Unclear cases remain unknown.',
    'performance_blinding': 'Policy and source hashes frozen before reading any binary-group performance; no performance input is permitted.'
}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

def inputs():
    records = pd.read_csv(ROOT/'cohort/records.csv', keep_default_na=False)
    chains = pd.read_csv(ROOT/'cohort/chains.csv', keep_default_na=False)
    annotations = pd.read_csv(ROOT/'cohort/annotation_architecture_candidates.csv', keep_default_na=False)
    assert len(records) == 3063 and records.sample_id.is_unique
    counts = chains.groupby('sample_id').size()
    assert all(int(r.n_chains) == counts[r.sample_id] for r in records.itertuples())
    return records, chains, annotations

def freeze_policy():
    HERE.mkdir(parents=True, exist_ok=True)
    policy_path = HERE/'policy_frozen.json'
    if policy_path.exists():
        assert json.loads(policy_path.read_text(encoding='utf-8')) == POLICY
    else:
        dump(policy_path, POLICY)
    hashes = {str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'cohort/records.csv', ROOT/'cohort/chains.csv', ROOT/'cohort/annotation_architecture_candidates.csv']}
    cache = Path('E:/Project/FlexPPI/experiments/_sources/functional_annotations_complete')
    hashes.update({str(p):sha(p) for p in sorted(cache.glob('batch_*.json'))})
    manifest = {'policy_sha256':sha(policy_path), 'inputs_sha256':hashes,
                'read_scope':'identity metadata and official annotations only; model performance is not read'}
    target = HERE/'source_hashes_frozen.json'
    if target.exists():
        assert json.loads(target.read_text(encoding='utf-8')) == manifest
    else:
        dump(target, manifest)
    return manifest

QUERY_BODY = '''rcsb_id struct {title} polymer_entities {
 rcsb_id entity_poly {pdbx_seq_one_letter_code_can type}
 rcsb_polymer_entity {pdbx_description pdbx_fragment}
 rcsb_polymer_entity_container_identifiers {auth_asym_ids uniprot_ids}
 rcsb_polymer_entity_align {reference_database_accession reference_database_name aligned_regions {entity_beg_seq_id ref_beg_seq_id length}}
}'''

def fetch():
    records, chains, annotations = inputs()
    ids = sorted(set(annotations.set_index('sample_id').loc[records[records.main_group_candidate=='other_multichain'].sample_id, 'pdb_ids']))
    assert len(ids) == 581 and all(re.fullmatch('[0-9A-Z]{4}', p) for p in ids)
    folder = HERE/'official_rcsb_payloads'
    folder.mkdir(exist_ok=True)
    batches = [ids[i:i+25] for i in range(0,len(ids),25)]
    def work(item):
        i, batch = item
        path=folder/('batch_%03d.json'%i)
        if path.exists():
            cached=json.loads(path.read_text(encoding='utf-8'))
            if cached['requested_ids']==batch and cached.get('response',{}).get('data',{}).get('entries'):
                return i, len(cached['response']['data']['entries']), 'cached'
        query = '{entries(entry_ids:' + json.dumps(batch) + '){' + QUERY_BODY + '}}'
        for attempt in range(3):
            try:
                result=requests.get('https://data.rcsb.org/graphql',params={'query':query},timeout=60)
                result.raise_for_status(); response=result.json()
                if response.get('errors'): raise ValueError(response['errors'])
                dump(path, {'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(), 'source_url':'https://data.rcsb.org/graphql', 'method':'GET', 'requested_ids':batch, 'query':query, 'response':response})
                return i,len(response['data']['entries']), 'downloaded'
            except Exception as exc:
                if attempt==2:
                    dump(folder/('failure_%03d.json'%i),{'requested_ids':batch,'error':repr(exc)})
                    return i,0,repr(exc)
                time.sleep(1+attempt)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for answer in pool.map(work, enumerate(batches)):
            print(json.dumps({'fetch_batch':answer[0],'entries':answer[1],'status':answer[2]}),flush=True)

def load_entities():
    output={}; sources={}
    for p in sorted((HERE/'official_rcsb_payloads').glob('batch_*.json')):
        value=json.loads(p.read_text(encoding='utf-8'))
        for entry in value.get('response',{}).get('data',{}).get('entries',[]) or []:
            if not entry: continue
            for entity in entry.get('polymer_entities',[]) or []:
                poly=entity.get('entity_poly') or {}
                seq=re.sub('[^A-Z]','',poly.get('pdbx_seq_one_letter_code_can',''))
                if not seq or not poly.get('type','').startswith('polypeptide'): continue
                ids=entity.get('rcsb_polymer_entity_container_identifiers') or {}
                desc=entity.get('rcsb_polymer_entity') or {}
                output.setdefault(entry['rcsb_id'],[]).append({
                    'pdb_id':entry['rcsb_id'],'entity_id':entity['rcsb_id'], 'sequence':seq,
                    'description':desc.get('pdbx_description',''),'fragment':desc.get('pdbx_fragment',''),
                    'uniprot_ids':sorted(set(ids.get('uniprot_ids') or [])),
                    'auth_asym_ids':ids.get('auth_asym_ids') or [],
                    'reference_alignments':entity.get('rcsb_polymer_entity_align') or [],
                    'payload_path':str(p.relative_to(HERE)), 'payload_sha256':sha(p)})
    return output

def sequence_match(query, target):
    if query in target:
        start=target.index(query)
        return {'match_fraction':1.,'identity':1.,'matches':len(query),'mismatches':0,'query_gaps':0,'target_gaps':0,'target_start':start+1,'target_stop':start+len(query),'method':'exact_substring','alignment_score':2*len(query)}
    if len(query)<20: return None
    aligner=PairwiseAligner(mode='global',match_score=2,mismatch_score=-2,open_gap_score=-5,extend_gap_score=-.1)
    # target is the official sequence, query is the observed MD input sequence.
    # Free terminal gaps in query permit MD truncations of the official entity.
    aligner.query_left_gap_score=0; aligner.query_right_gap_score=0
    alignment=aligner.align(target,query)[0]
    coords=alignment.coordinates
    matches=mismatches=qgap=tgap=0; starts=[]; stops=[]
    for j in range(coords.shape[1]-1):
        t0,t1=int(coords[0,j]),int(coords[0,j+1]); q0,q1=int(coords[1,j]),int(coords[1,j+1])
        if t1>t0 and q1>q0:
            for x,y in zip(target[t0:t1],query[q0:q1]):
                if y=='X': continue
                if x==y: matches+=1
                else: mismatches+=1
            starts.append(t0);stops.append(t1)
        elif q1>q0: tgap+=q1-q0
        elif q0>0 and q0<len(query): qgap+=t1-t0
    known=len(query)-query.count('X')
    return {'match_fraction':matches/max(known,1), 'identity':matches/max(matches+mismatches,1),
            'matches':matches,'mismatches':mismatches,'query_gaps':qgap,'target_gaps':tgap,
            'target_start':min(starts)+1 if starts else None,'target_stop':max(stops) if stops else None,
            'method':'semiglobal_sequence_alignment','alignment_score':float(alignment.score)}

def classify():
    records, chains, annotations=inputs(); entities=load_entities()
    source_overrides=json.loads((HERE/'source_identity_overrides.json').read_text(encoding='utf-8')) if (HERE/'source_identity_overrides.json').exists() else {}
    reviewed_matches=json.loads((HERE/'reviewed_sequence_exceptions.json').read_text(encoding='utf-8')) if (HERE/'reviewed_sequence_exceptions.json').exists() else {}
    ann=annotations.set_index('sample_id'); all_matches=[]; membership=[]
    chain_groups={sid:sub for sid,sub in chains.groupby('sample_id')}
    for ri,record in enumerate(records.itertuples()):
        chain_table=chain_groups[record.sample_id]
        base={'sample_id':record.sample_id,'benchmark_group':'unknown','classification_status':'unknown',
              'confidence':'unresolved','actual_input_n_chains':len(chain_table),'pdb_id':ann.loc[record.sample_id,'pdb_ids'],
              'evidence_rule':'','evidence':'','policy_version':POLICY['version']}
        base['mapping_pdb_id']=source_overrides.get(record.sample_id,{}).get('mapping_pdb_id',base['pdb_id'])
        base['source_identity_conflict_note']=source_overrides.get(record.sample_id,{}).get('reason','')
        if len(chain_table)==1:
            base.update(benchmark_group='protein',classification_status='confirmed_single_input_chain',confidence='high',evidence_rule='authorized_single_chain',evidence='Actual MD input contains exactly one chain; source assembly composition is not used.')
        elif record.main_group_candidate=='tcr_pmhc_itrace':
            base.update(benchmark_group='protein_complex',classification_status='confirmed_tcr_pmhc_source',confidence='high',evidence_rule='authorized_tcr_pmhc',evidence='Source-confirmed ITrace TCR-pMHC; explicit human-authorized fixed membership.')
        else:
            source=entities.get(base['mapping_pdb_id'],[]); accepted=[]
            for chain in chain_table.itertuples():
                candidates=[]
                if chain.sequence.count('X')/len(chain.sequence)<=.05:
                    for entity in source:
                        m=sequence_match(chain.sequence,entity['sequence'])
                        exception=reviewed_matches.get(record.sample_id+':'+chain.chain)
                        reviewed=bool(m and exception and entity['entity_id']==exception['entity_id'] and m['match_fraction']>=exception['minimum_match_fraction'])
                        if m and (m['match_fraction']>=.97 or reviewed):
                            m['reviewed_sequence_exception']=exception['evidence'] if reviewed else ''
                            candidates.append((entity,m))
                # Keep all close candidates, rather than choosing by MD chain letter.
                if candidates:
                    best=max(m['alignment_score'] for e,m in candidates)
                    candidates=[(e,m) for e,m in candidates if m['alignment_score']>=best-max(1.,.02*len(chain.sequence))]
                accepted.append((chain,candidates))
                for entity,m in candidates:
                    all_matches.append({'sample_id':record.sample_id,'input_chain':chain.chain,'input_id':chain.input_id,
                        'input_sequence_sha256':chain.chain_sequence_sha256,'input_length':len(chain.sequence),
                        **{k:v for k,v in entity.items() if k not in ['sequence','reference_alignments']},**m})
                if not candidates:
                    all_matches.append({'sample_id':record.sample_id,'input_chain':chain.chain,'input_id':chain.input_id,
                        'input_sequence_sha256':chain.chain_sequence_sha256,'input_length':len(chain.sequence),
                        'pdb_id':base['pdb_id'],'entity_id':'','description':'','method':'unmatched'})
            base['evidence']=json.dumps([{'chain':c.chain,'sequence_sha256':c.chain_sequence_sha256,
                'candidate_entities':[e['entity_id'] for e,m in cand],
                'candidate_uniprot_ids':[e['uniprot_ids'] for e,m in cand],
                'descriptions':[e['description'] for e,m in cand]} for c,cand in accepted],ensure_ascii=False)
            # First pass intentionally conservative. Identity refinement is explicit
            # and audited separately, never inferred from entity count alone.
            if accepted and all(cand for c,cand in accepted):
                unique_entities={e['entity_id'] for c,cand in accepted for e,m in cand}
                if len(unique_entities)==1:
                    base.update(benchmark_group='protein',classification_status='confirmed_same_official_entity_by_sequence',confidence='high',evidence_rule='all_MD_chains_sequence_map_to_one_official_polymer_identity')
                else:
                    # A disjoint pair of nonempty UniProt identity sets proves
                    # multiple proteins without depending on unresolved peptides.
                    forced=[]
                    for c,cand in accepted:
                        common=set(cand[0][0]['uniprot_ids'])
                        for e,m in cand[1:]: common &= set(e['uniprot_ids'])
                        forced.append(common)
                    distinct=any(a and b and a.isdisjoint(b) for i,a in enumerate(forced) for b in forced[i+1:])
                    if distinct:
                        base.update(benchmark_group='protein_complex',classification_status='confirmed_distinct_uniprot_by_sequence',confidence='high',evidence_rule='sequence_mapped_input_chains_have_disjoint_official_UniProt_identities')
                    else:
                        base.update(classification_status='unknown_shared_accession_or_unannotated_identity',evidence_rule='needs_mature_product_or_identity_review')
            else:
                base.update(classification_status='unknown_input_sequence_mapping_incomplete',evidence_rule='one_or_more_MD_chains_lack_a_confident_official_sequence_match')
        membership.append(base)
        if ri%200==0: print(json.dumps({'classified_records':ri,'total':3063}),flush=True)
    match_df=pd.DataFrame(all_matches)
    for col in ['uniprot_ids','auth_asym_ids']:
        if col in match_df: match_df[col]=match_df[col].apply(lambda v:json.dumps(v) if isinstance(v,list) else '')
    match_df.to_csv(HERE/'chain_entity_matches.csv',index=False)
    members=pd.DataFrame(membership)
    members.to_csv(HERE/'membership_automatic_before_reviews.csv',index=False)
    # Decisions are written with their evidence before they can affect membership.
    # This is an explicit record list, never a default for unresolved records.
    for decision_file in ['reviewed_identity_decisions.json','boundary_decisions_approved.json']:
        path=HERE/decision_file
        if not path.exists():continue
        for decision in json.loads(path.read_text(encoding='utf-8')):
            selected=members.index[members.sample_id==decision['sample_id']]
            assert len(selected)==1,decision['sample_id']
            ix=selected[0]
            assert decision['benchmark_group'] in ['protein','protein_complex']
            members.loc[ix,'benchmark_group']=decision['benchmark_group']
            members.loc[ix,'classification_status']=decision['classification_status']
            members.loc[ix,'confidence']=decision['confidence']
            members.loc[ix,'evidence_rule']='explicit_review_with_saved_official_annotation_evidence'
            members.loc[ix,'review_decision_file']=decision_file
            members.loc[ix,'review_decision_sha256']=sha(path)
            members.loc[ix,'review_reason']=decision['reason']
    members.to_csv(HERE/'membership_frozen.csv',index=False)
    members[members.benchmark_group=='unknown'].to_csv(HERE/'unresolved_records.csv',index=False)
    members[members.source_identity_conflict_note!=''][['sample_id','pdb_id','mapping_pdb_id','source_identity_conflict_note','benchmark_group']].to_csv(HERE/'source_identity_conflicts.csv',index=False)
    dump(HERE/'official_entities_normalized.json',entities)
    summary={'records':len(members),'counts':members.benchmark_group.value_counts().to_dict(),
             'classification_status_counts':members.classification_status.value_counts().to_dict(),
             'source_PDB_with_sequences':len(entities),'other_multichain_records':626,
             'policy_sha256':sha(HERE/'policy_frozen.json'),'membership_sha256':sha(HERE/'membership_frozen.csv'),
             'model_performance_read':False,'model_inference_run':False}
    dump(HERE/'classification_summary.json',summary)
    print(json.dumps(summary,ensure_ascii=False),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--fetch',action='store_true');p.add_argument('--classify',action='store_true')
    args=p.parse_args(); freeze_policy()
    if args.fetch: fetch()
    if args.classify: classify()

if __name__=='__main__': main()
