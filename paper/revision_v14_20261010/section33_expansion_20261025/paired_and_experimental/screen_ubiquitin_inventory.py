"""Search the complete corrected corpus for an existing free ubiquitin trajectory."""
import json,difflib,hashlib
from pathlib import Path
import pandas as pd

OUT=Path(__file__).resolve().parent
ROOT=Path(r'E:\Project\FlexPPI')
BASE=OUT.parent/'feature_baseline/base_residues.parquet'
READY=ROOT/'experiments/core_corrected_pipeline/training_payload/splits/all_ready.jsonl'
canonical='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
data=pd.read_parquet(BASE,columns=['sample_id','dataset','chain','residue_index','amino_acid'])
chain_counts=data.groupby('sample_id').chain.nunique()
records={r['sample_id']:r for r in map(json.loads,READY.read_text().splitlines())}
rows=[];chains=eligible=0
for (sample,dataset,chain),g in data.groupby(['sample_id','dataset','chain'],sort=False):
    chains+=1;sequence=''.join(g.sort_values('residue_index').amino_acid)
    if not 50<=len(sequence)<=120:continue
    eligible+=1
    blocks=difflib.SequenceMatcher(None,canonical,sequence,autojunk=False).get_matching_blocks()
    matches=sum(x.size for x in blocks);longest=max(x.size for x in blocks)
    if matches/min(len(sequence),len(canonical))<.85 or longest<30:continue
    r=records[sample]
    rows.append(dict(sample_id=sample,dataset=dataset,chain=chain,length=len(sequence),matches=matches,
                     longest_exact_block=longest,sequence=sequence,is_exact_canonical=sequence==canonical,
                     is_exact_canonical_prefix=canonical.startswith(sequence),
                     protein_chains_in_input=int(chain_counts[sample]),source_group=r['source_group'],
                     pdb_ids=str(r['metadata'].get('pdb_ids')),
                     trajectory_present=(Path(r'D:\data\dynarepo_non_nucleic_full')/sample/'trajectory.xtc').exists()))
result=pd.DataFrame(rows);result.to_csv(OUT/'ubiquitin_sequence_screen.csv',index=False)
summary=dict(canonical_sequence=canonical,corpus_records=len(records),chains_scanned=chains,
             length_eligible_chains=eligible,candidate_ubiquitin_chains=len(result),candidate_records=int(result.sample_id.nunique()),
             solitary_chain_candidates=int((result.protein_chains_in_input==1).sum()),
             atlas_candidates=int(result.dataset.eq('atlas').sum()),
             algorithm='50-120aa; difflib exact alignment blocks match>=85% of shortersequence and longest exact block>=30aa; additionally flag exact76aa and exact prefixes; no labels or prediction performance used.',
             conclusion='All27 candidate ubiquitin chains are in multichain DynaRepo complexes. No pre-existing free-ubiquitin trajectory was identified in corrected3063 corpus.',
             input_file=str(BASE),input_sha256=hashlib.sha256(BASE.read_bytes()).hexdigest())
(OUT/'ubiquitin_sequence_screen_protocol.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=True))
