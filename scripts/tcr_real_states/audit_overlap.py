"""Conservative homology flags, separately for TCR targets and protein partners."""
from pathlib import Path
import json
import hashlib
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'artifacts/tcr_real_states_20260915'


def main():
    path=ROOT/'artifacts/apo_holo_exploration_20260913/training_sequence_index.json'
    records=json.loads(path.read_text())
    train={}
    for r in records:
        for s in r['sequences']:
            train.setdefault(s,[]).append(dict(sample_id=r['sample_id'],held_out_fold=r['held_out_fold'],pdb_ids=r.get('pdb_ids',[])))
    queries={}
    for s in json.loads((BASE/'states.json').read_text()):
        if s['coverage']!='full' or s['context']=='tcr_only':
            continue
        for ch,r in s['residues'].items():
            queries.setdefault(r['sequence'],[]).append(dict(state=s['id'],case=s['case'],chain=ch,
                role='target_TCR' if ch in s['target_chains'] else 'protein_partner'))
    al=PairwiseAligner(mode='local',open_gap_score=-10,extend_gap_score=-.5)
    al.substitution_matrix=substitution_matrices.load('BLOSUM62')
    hits=[]
    for qi,(q,origins) in enumerate(queries.items()):
        for seq,recs in train.items():
            alignment=next(iter(al.align(q,seq)),None)
            if alignment is None:  # No positive-score local alignment, e.g. a short peptide.
                continue
            count=exact=0
            for (b,e),(bb,ee) in zip(*alignment.aligned):
                count+=e-b
                exact+=sum(x==y for x,y in zip(q[b:e],seq[bb:ee]))
            identity=exact/max(count,1);coverage=count/len(q)
            if identity>=.3 and coverage>=.8:
                hits.append(dict(query_id=qi,origins=origins,query_length=len(q),training_length=len(seq),identity=identity,
                    query_aligned_coverage=coverage,exact_sequence=(q==seq),records=recs))
        print('Overlap query',qi+1,'/',len(queries),flush=True)
    summaries=[]
    for case in ['B4.2.3','1G4-CD3']:
        for role in ['target_TCR','protein_partner']:
            hs=[h for h in hits if any(o['case']==case and o['role']==role for o in h['origins'])]
            held={r['held_out_fold'] for h in hs for r in h['records']}
            summaries.append(dict(case=case,role=role,n_alignment_hits=len(hs),n_unique_training_samples=len({r['sample_id'] for h in hs for r in h['records']}),
                maximum_identity=max([h['identity'] for h in hs],default=None),n_exact_sequence_hits=sum(h['exact_sequence'] for h in hs),
                training_samples_held_out_folds=sorted(held),
                checkpoints_without_detected_homolog=[f for f in range(5) if not held or held=={f}]))
    report=dict(method='Exhaustive local BLOSUM62, gap open 10 / extend 0.5; identity >= 0.30 and query aligned coverage >= 0.80.',
        training_index_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),n_unique_training_sequences=len(train),n_query_sequences=len(queries),
        caveat='Flags candidate homologs, not a family-independent evaluation certificate. No BLAST E-values or functional conclusions.',summary=summaries,hits=hits)
    (BASE/'training_overlap_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(summaries,indent=2),flush=True)


if __name__=='__main__':
    main()
