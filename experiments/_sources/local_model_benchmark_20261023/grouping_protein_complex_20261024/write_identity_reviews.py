"""Record explicit per-record metadata review before applying final membership.

This script reads only chain identity evidence and official entry annotations.
The list is explicit so a newly unresolved record can never be silently defaulted
to protein_complex by a catch-all rule.
"""
from pathlib import Path
import hashlib, json
import pandas as pd

HERE=Path(__file__).resolve().parent
COMPLEX_PDBS=set('''
1KXT 1LP1 1OP9 2AEQ 2BDN 2IFF 2OIN 2V5Q 2YBR
3HG0 3K2M 3OGO 3QHT 3QWQ 3RZW 3UYO 3ZU7
4EEF 4EIZ 4GN4 4HEP 4I0C 4JW3 4K5A 4LGP 4UU9 4XL5 4XVP
5DJU 5DMJ 5E7F 5ELU 5IP4 5KY5 5L21 5M2M 5ML9 5N88 5NH3 5OMN 5OTJ
6BB4 6BX5 6CWG 6H3T 6H71 6IVZ 4OYD 2XWT 4LVN 5J57 5TIH
1WEJ 2VXT 3EOA 3G6D 3HI6 3L5W 3MJ9 3MXW 3RJQ 3SE8 3V6Z 3WD5
4DN4 4DW2 4ETQ 4FP8 4FQI 4G6J 4M5Z 5C7X 5GRJ 5HYS 5O14 5WK3 5WUX
5X0T 5Y9J 6A0Z 6A77 6AL0 6BPC 4M3K 4POU 4Y7M 5E5M 5HGG 5JMO 5SV3
5VNW 6DBG 6EY6 3BX7 2FJG 2W9E 5VPG
'''.split())
PENDING_BOUNDARIES=set('1E3U 3F4Y 3P9W 4LLO 5JJD 3VV2 4ZRJ 1RKE 4BQD 4JO6'.split())

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    members=pd.read_csv(HERE/'membership_automatic_before_reviews.csv',keep_default_na=False)
    entities=json.loads((HERE/'official_entities_normalized.json').read_text(encoding='utf-8'))
    reviews=[]; pending=[]
    for row in members[members.benchmark_group=='unknown'].itertuples():
        pdb=row.mapping_pdb_id
        assert pdb in COMPLEX_PDBS|PENDING_BOUNDARIES,(row.sample_id,pdb)
        if pdb in PENDING_BOUNDARIES:
            pending.append(row._asdict());continue
        evidence=json.loads(row.evidence)
        ids={eid for ch in evidence for eid in ch['candidate_entities']}
        selected=[e for e in entities[pdb] if e['entity_id'] in ids]
        payload=HERE/selected[0]['payload_path']
        source=json.loads(payload.read_text(encoding='utf-8'))
        entry=next(e for e in source['response']['data']['entries'] if e and e['rcsb_id']==pdb)
        title=entry['struct']['title']
        names=[{'entity_id':e['entity_id'],'name':e['description'],'fragment':e['fragment'],'uniprot_ids':e['uniprot_ids']} for e in selected]
        reason='The actual MD sequences match the listed official proteins; their mature names and the official entry title identify a target protein and a separate antibody, designed binder, or independently named protein partner. Missing UniProt accessions on engineered binders do not make them the target protein.'
        if pdb=='2OIN':
            reason='Official fragment field identifies HCV protease domain NS3 and the second entity is NS4A peptide. They are separately functioning mature products, despite sharing precursor accession P27958; the frozen policy explicitly classifies NS3 plus NS4A as protein_complex.'
        if pdb=='5HYS':
            reason='Official entity names for two chains are generic, but the official entry title explicitly identifies IgE complexed with omalizumab. Actual MD chains sequence-map to those entities and Ig epsilon C region (P01854), establishing distinct antibody and IgE proteins.'
        reviews.append({'sample_id':row.sample_id,'original_pdb_id':row.pdb_id,'mapping_pdb_id':pdb,
            'benchmark_group':'protein_complex','classification_status':'reviewed_official_protein_names_and_MD_sequence',
            'confidence':'high','entry_title':title,'matched_mature_identity_annotations':names,
            'reason':reason,'actual_chain_sequence_evidence':evidence,
            'source_payload':str(payload.relative_to(HERE)),'source_payload_sha256':sha(payload),
            'review_scope':'Pure metadata curation; no performance or model output read.'})
    target=HERE/'reviewed_identity_decisions.json'
    target.write_text(json.dumps(reviews,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    pd.DataFrame([{k:json.dumps(v,ensure_ascii=False) if isinstance(v,list) else v for k,v in x.items()} for x in reviews]).to_csv(HERE/'reviewed_identity_decisions.csv',index=False)
    pd.DataFrame(pending).to_csv(HERE/'pending_boundary_records.csv',index=False)
    print(json.dumps({'explicit_reviewed_complex_records':len(reviews),'boundary_records_pending':len(pending),'decision_file_sha256':sha(target)}))

if __name__=='__main__':main()
