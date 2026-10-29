"""Model-fold sensitivity and conservative *non-exhaustive* homology flags.

The k-mer prefilter finds obvious relatives; it cannot certify sequence novelty.
Exported queries/training FASTA permit a later exhaustive MMseqs2/BLAST audit.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict, Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/apo_holo_exploration_20260913'
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from flexppi.data.structures import parse_pdb
from flexppi.data.pairing import global_sequence_alignment
from analyze_apo_holo_exploration import corr,partial_rank_correlation


def main():
    index=json.loads((OUT/'training_sequence_index.json').read_text())
    sequences=defaultdict(list)
    for r in index:
        for seq in r['sequences']:
            if len(seq)>=30:
                sequences[seq].append(dict(sample_id=r['sample_id'],held_out_fold=r['held_out_fold']))
    seqs=list(sequences)
    inverted=defaultdict(set)
    for i,s in enumerate(seqs):
        for kmer in set(s[j:j+4] for j in range(len(s)-3)):
            inverted[kmer].add(i)
    query_cases=defaultdict(set)
    # Examine whole target components, not just residues observed in both states.
    for path in (OUT/'mapping').glob('*.json'):
        record=json.loads(path.read_text())
        case=record['case']
        for spec,chains in [(record['selected_apo'],record['selected_apo']['chains']),
                            (case['holo'],case['target_chains'])]:
            s=parse_pdb(spec['path'])
            for chain in chains:
                seq=s.chain_sequence(chain)
                if len(seq)>=30:
                    query_cases[seq].add(case['case_id'])
    (OUT/'analysis/training_chains.fasta').write_text(''.join(f'>train_{i}\n{s}\n' for i,s in enumerate(seqs)))
    (OUT/'analysis/query_chains.fasta').write_text(''.join(f'>query_{i}\n{s}\n' for i,s in enumerate(query_cases)))
    hits=[]
    for qi,(query,cases) in enumerate(query_cases.items()):
        kmers=set(query[j:j+4] for j in range(len(query)-3))
        counts=Counter(i for k in kmers for i in inverted[k])
        for idx,count in counts.items():
            seq=seqs[idx]
            if count<max(2,.04*len(kmers)) or min(len(query),len(seq))/max(len(query),len(seq))<.6:
                continue
            alignment=global_sequence_alignment(query,seq)
            # Coverage here counts exact matches, not aligned mismatches; record
            # explicitly to avoid confusing this heuristic with BLAST coverage.
            exact_coverage=min(len(alignment.indices_a)/len(query),len(alignment.indices_b)/len(seq))
            if alignment.identity>=.3 and exact_coverage>=.3:
                hits.append(dict(query_id=f'query_{qi}',case_ids=sorted(cases),query_length=len(query),train_id=f'train_{idx}',
                    train_length=len(seq),identity=alignment.identity,bidirectional_exact_match_coverage=exact_coverage,
                    training_records=sequences[seq]))
        print('homology screen',qi+1,'/',len(query_cases),flush=True)
    (OUT/'analysis/detected_homology.json').write_text(json.dumps(dict(
        limitation='Non-exhaustive 4-mer prefilter (>=2 and >=4% query kmers); length ratio >=0.6; global alignment. Absence of hits is NOT proof of independence.',hits=hits),indent=2))
    overlap=json.loads((OUT/'training_overlap.json').read_text())
    clean={r['case_id']:r['folds_without_detected_exact_overlap'] for r in overlap}
    for h in hits:
        folds={r['held_out_fold'] for r in h['training_records']}
        for cid in h['case_ids']:
            if cid in clean:
                clean[cid]=[k for k in clean[cid] if folds=={k}]
    (OUT/'analysis/folds_without_detected_overlap.json').write_text(json.dumps(clean,indent=2))
    raw=pd.read_csv(OUT/'residue_predictions.csv',dtype={'apo_author_id':str,'author_id':str})
    nmr=pd.read_csv(OUT/'analysis/calmodulin_prediction_nmr_paired.csv',dtype={'apo_author_id':str})
    rows=[]
    for (cid,fold),g in raw[raw.direction=='backbone_sidechain'].groupby(['case_id','fold']):
        m=g.merge(nmr,on=['case_id','apo_author_id','amino_acid'])
        for observable,metric in [('backbone','rmsf'),('methyl','chi_entropy')]:
            r=m[m.observable==observable]
            rows.append(dict(case_id=cid,fold=fold,observable=observable,n=len(r),
                spearman=corr(r['delta_'+metric],r.minus_delta_s2),
                partial_aa_distance=partial_rank_correlation(r['delta_'+metric],r.minus_delta_s2,r.partner_distance,r.amino_acid),
                detected_overlap_absent=fold in clean[cid]))
    pd.DataFrame(rows).to_csv(OUT/'analysis/calmodulin_fold_sensitivity.csv',index=False)
    rows=[]
    for (cid,fold),g in raw[raw.direction=='md_delta'].groupby(['case_id','fold']):
        for metric in ('rmsf','chi_entropy'):
            ok=g['true_delta_'+metric].notna()&g['delta_'+metric].notna()
            p,y=g.loc[ok,'delta_'+metric],g.loc[ok,'true_delta_'+metric]
            rows.append(dict(case_id=cid,group=g.group.iloc[0],fold=fold,metric=metric,spearman=corr(p,y),
                mae=np.mean(np.abs(p-y)),zero_mae=np.mean(np.abs(y)),
                detected_overlap_absent=fold in clean[cid]))
    pd.DataFrame(rows).to_csv(OUT/'analysis/md_fold_sensitivity.csv',index=False)
    print('Candidate homology hits:',len(hits))
    print('CaM/PDZ folds without detected overlap:',{k:v for k,v in clean.items() if k.startswith(('cam_','pdz'))})


if __name__=='__main__':
    main()
