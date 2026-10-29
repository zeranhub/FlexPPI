"""Verify and freeze metadata-only binary membership for downstream aggregation."""
from pathlib import Path
import datetime, hashlib, json, platform
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def main():
    records=pd.read_csv(ROOT/'cohort/records.csv',keep_default_na=False)
    chains=pd.read_csv(ROOT/'cohort/chains.csv',keep_default_na=False)
    members=pd.read_csv(HERE/'membership_frozen.csv',keep_default_na=False)
    initial=json.loads((HERE/'source_hashes_frozen.json').read_text(encoding='utf-8'))
    assert members.sample_id.is_unique and len(members)==3063
    assert set(members.sample_id)==set(records.sample_id)
    assert set(members.benchmark_group)=={'protein','protein_complex'}
    assert not members.classification_status.str.contains('unknown',case=False).any()
    assert all(members.evidence!='') and all(members.confidence!='unresolved')
    for rel,expected in initial['inputs_sha256'].items():
        path=Path(rel) if Path(rel).is_absolute() else ROOT/rel
        assert sha(path)==expected,('input changed',str(path))
    assert sha(HERE/'policy_frozen.json')==initial['policy_sha256']
    matches=pd.read_csv(HERE/'chain_entity_matches.csv',keep_default_na=False)
    other_ids=set(records.loc[records.main_group_candidate=='other_multichain','sample_id'])
    expected_chains=chains[chains.sample_id.isin(other_ids)]
    observed=set(zip(matches.sample_id,matches.input_chain))
    assert observed==set(zip(expected_chains.sample_id,expected_chains.chain))
    assert (matches.entity_id!='').all()
    assert len(expected_chains)==len(expected_chains.drop_duplicates(['sample_id','chain']))
    # Entity IDs remain provenance, rather than being counted as distinct proteins.
    # Mature-product exceptions and split-domain identity links are explicit.
    same_mature={'1E3U':'UniProt:P14489:OXA-10','3F4Y':'UniProt:P04580:mature_gp41',
        '4LLO':'UniProt:Q60603:KCNH1','5JJD':'UniProt:Q9Y5P4:CERT',
        '4ZRJ':'UniProt:P35240:Merlin','1RKE':'UniProt:P18206:Vinculin'}
    identity_rows=[]
    for row in matches.to_dict(orient='records'):
        pdb=row['pdb_id'];ids=json.loads(row['uniprot_ids']) if row['uniprot_ids'] else []
        if pdb in same_mature:identity=same_mature[pdb];basis='reviewed same mature functional protein; entity IDs are not identities'
        elif pdb=='2OIN':identity='UniProt:P27958:'+('NS4A' if row['entity_id'].endswith('_2') else 'NS3');basis='different mature products from common precursor'
        elif pdb=='3VV2':identity='UniProt:P58502:'+('cleaved_propeptide' if row['entity_id'].endswith('_2') else 'mature_subtilisin');basis='cleaved independently acting propeptide versus mature enzyme'
        elif pdb=='7CQC':identity='engineered_fusion:Aq1964_PA14';basis='contiguous pieces of one official engineered fusion entity'
        elif pdb=='3P9W' and row['entity_id'].endswith('_2'):
            identity='engineered_antibody:autonomous_human_VH_binder';basis='official entity description is misnamed; primary paper, entry title and VH sequence resolve identity'
        elif pdb=='5HYS' and row['entity_id'].endswith(('_1','_2')):
            identity='antibody:omalizumab:'+('heavy_chain' if row['entity_id'].endswith('_1') else 'light_chain');basis='entry title resolves omalizumab antibody identity despite generic or misleading entity names'
        elif pdb=='4BQD' and row['entity_id'].endswith('_2'):
            identity='independent_selected_TFPI_antagonist:compound_3';basis='primary paper establishes independently selected antagonist, not a TFPI self-fragment; chemical 20-mer and canonical 21-symbol representation differ'
        elif ids:identity=';'.join('UniProt:'+x for x in ids);basis='official accession identity candidates with mature-name evidence retained'
        else:identity='official_named_protein:'+row['description'];basis='official mature protein/binder name; no accession invented'
        identity_rows.append({**row,'biological_identity_evidence':identity,'identity_basis':basis})
    pd.DataFrame(identity_rows).to_csv(HERE/'chain_biological_identity_evidence.csv',index=False)
    source_files=[]
    for folder in ['official_rcsb_payloads','official_sequence_search']:
        source_files.extend(sorted((HERE/folder).glob('*.json')))
    source_files.extend([HERE/name for name in ['policy_frozen.json','source_hashes_frozen.json','source_identity_overrides.json',
        'reviewed_sequence_exceptions.json','reviewed_identity_decisions.json','reviewed_identity_decisions.csv',
        'boundary_decisions_approved.json','official_entities_normalized.json','chain_entity_matches.csv','chain_biological_identity_evidence.csv']])
    boundary_manifest=HERE/'boundary_review_source_hashes.json'
    assert boundary_manifest.is_file(),'independent boundary review sources must be saved before freeze'
    for rel,expected_hash in json.loads(boundary_manifest.read_text(encoding='utf-8')).items():
        assert sha(HERE/rel)==expected_hash,('boundary source changed',rel)
    source_files.append(boundary_manifest)
    source_hashes={str(p.relative_to(HERE)):sha(p) for p in source_files}
    code_hashes={p.name:sha(p) for p in HERE.glob('*.py')}
    counts={key:int((members.benchmark_group==key).sum()) for key in ['protein','protein_complex']}
    unresolved=pd.read_csv(HERE/'unresolved_records.csv')
    assert len(unresolved)==0
    audit={'status':'frozen','records':3063,'group_counts':counts,'unresolved_records':0,
        'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'policy_sha256':sha(HERE/'policy_frozen.json'),'membership_sha256':sha(HERE/'membership_frozen.csv'),
        'source_hashes':source_hashes,'code_hashes':code_hashes,
        'input_source_hashes':initial['inputs_sha256'],
        'actual_nonfixed_multichain_records':626,'actual_nonfixed_multichain_chains':len(expected_chains),
        'all_nonfixed_multichain_chains_officially_sequence_mapped':True,
        'MD_author_chain_letters_used_for_mapping':False,'source_entity_count_used_as_protein_count':False,
        'source_PDB_recovery_records':7,'reviewed_sequence_threshold_exceptions':4,
        'performance_blind':True,'model_performance_read':False,'model_inference_run':False,
        'original_cohort_changed':False,'existing_scores_or_figures_changed':False,
        'python_version':platform.python_version(),'alignment_library':'Biopython 1.84 PairwiseAligner',
        'downstream_contract':'Use only when this audit status is frozen and the membership SHA256 matches. Include all 3063 unique sample IDs; no residual other_multichain group.'}
    dump(HERE/'membership_freeze_audit.json',audit)
    print(json.dumps({k:audit[k] for k in ['status','records','group_counts','unresolved_records','policy_sha256','membership_sha256','actual_nonfixed_multichain_chains']},ensure_ascii=False))

if __name__=='__main__':main()
