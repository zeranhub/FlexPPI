"""Exhaustive local sequence-alignment screen of target chains vs training chains.

This detects candidate homologs, not a certified family-independent benchmark.
No k-mer or length prefilter: every unique target sequence is compared against
all3879 unique training-chain sequences. Thresholds are conservative flags.
"""
from pathlib import Path
import json
import hashlib
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/mechanistic_case_studies_20260913'

def main():
    idx=ROOT/'artifacts/apo_holo_exploration_20260913/training_sequence_index.json'
    records=json.loads(idx.read_text())
    train={}
    for r in records:
        for s in r['sequences']:
            train.setdefault(s,[]).append(dict(sample_id=r['sample_id'],held_out_fold=r['held_out_fold']))
    specs={s['id']:s for s in json.loads((BASE/'states.json').read_text())}
    cons=json.loads((BASE/'contrasts.json').read_text())
    queries={}
    for c in cons:
        for sid,chains in [(c['before'],c['chainmap'].keys()),(c['after'],c['chainmap'].values())]:
            for ch in chains:
                seq=specs[sid]['residues'][ch]['sequence']
                queries.setdefault(seq,set()).add(c['case'])
    al=PairwiseAligner(mode='local',open_gap_score=-10,extend_gap_score=-.5)
    al.substitution_matrix=substitution_matrices.load('BLOSUM62')
    hits=[]
    for qi,(q,cases) in enumerate(queries.items()):
        for s,recs in train.items():
            a=al.align(q,s)[0]
            count=exact=0
            for (b,e),(bb,ee) in zip(*a.aligned):
                count+=e-b
                exact+=sum(x==y for x,y in zip(q[b:e],s[bb:ee]))
            cov=count/len(q);identity=exact/max(1,count)
            if cov>=.8 and identity>=.3:
                hits.append(dict(query_id=qi,cases=sorted(cases),query_length=len(q),training_length=len(s),
                    identity=identity,query_aligned_coverage=cov,exact_sequence=q==s,records=recs))
        print('Query',qi+1,'/',len(queries),flush=True)
    clean={}
    for case in sorted({v for cases in queries.values() for v in cases}):
        hs=[h for h in hits if case in h['cases']]
        held={r['held_out_fold'] for h in hs for r in h['records']}
        clean[case]=[f for f in range(5) if not held or held=={f}]
    report=dict(training_index_sha256=hashlib.sha256(idx.read_bytes()).hexdigest(),unique_training_sequences=len(train),
        unique_query_sequences=len(queries),method='Exhaustive pairwise local BLOSUM62; gap open10/extend0.5; identity>=0.3 and aligned query coverage>=0.8',
        caveat='Descriptive homology flag, not calibrated BLAST E-value or proof of family independence. Partner sequences not included.',
        folds_without_detected_target_homolog=clean,hits=hits)
    (BASE/'training_overlap_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(clean,flush=True)

if __name__=='__main__':
    main()
