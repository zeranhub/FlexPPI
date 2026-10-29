"""Conservative sequence novelty screen; not a formal IMGT annotation."""
from pathlib import Path
import json
import re
import hashlib
from collections import defaultdict

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/tcr_external_validation_search_20260917'


def main():
    source = ROOT / 'artifacts/apo_holo_exploration_20260913/training_sequence_index.json'
    training = defaultdict(list)
    for record in json.loads(source.read_text()):
        for seq in record['sequences']:
            training[seq].append({k: record[k] for k in ['sample_id', 'held_out_fold', 'pdb_ids']})
    static = json.loads((OUT / 'static_agent/static_structure_audit.json').read_text())
    aligner = PairwiseAligner(mode='local', open_gap_score=-10, extend_gap_score=-0.5)
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    rows, detailed = [], []
    for pdb in ['9YW4', '9ZCL', '9PBG', '9WC2']:
        for entity in static['entries'][pdb]['entities']:
            if not entity['is_tcr']:
                continue
            seq = entity['sequence']
            # Sequence-anchor window is used for novelty screening ONLY.
            # IMGT assignment and exact CDR endpoints require independent annotation.
            anchor = re.search(r'[FW]G.G', seq[75:135])
            assert anchor is not None
            fpos = 75 + anchor.start()
            cpos = seq.rfind('C', 65, fpos)
            assert cpos >= 0
            loop = seq[cpos+1:fpos]
            v_window = seq[:fpos+10]
            exact_full, exact_v, exact_loop, hits = [], [], [], []
            for trained, records in training.items():
                if seq == trained:
                    exact_full.extend(records)
                if v_window in trained:
                    exact_v.extend(records)
                if loop in trained:
                    exact_loop.extend(records)
                alignments = aligner.align(v_window, trained)
                alignment = next(iter(alignments), None)
                if alignment is None:
                    continue
                aligned = identical = 0
                for (a, b), (c, d) in zip(*alignment.aligned):
                    aligned += b-a
                    identical += sum(x == y for x, y in zip(v_window[a:b], trained[c:d]))
                coverage = aligned / len(v_window)
                identity = identical / max(aligned, 1)
                if coverage >= 0.8 and identity >= 0.7:
                    hits.append(dict(identity=identity, query_coverage=coverage,
                                     training_sequence=trained, records=records,
                                     alignment=str(alignment)))
            hits.sort(key=lambda h: (h['identity'], h['query_coverage']), reverse=True)
            samples = lambda records: sorted({r['sample_id'] for r in records})
            nearest = hits[0] if hits else None
            row = dict(pdb=pdb, entity_id=entity['entity_id'], description=entity['description'],
                provisional_v_window_length=len(v_window), provisional_cdr3_sequence=loop,
                loop_label_seq_start=cpos+2, loop_label_seq_end=fpos,
                exact_full_sequence_samples=';'.join(samples(exact_full)),
                exact_v_window_samples=';'.join(samples(exact_v)),
                exact_cdr3_fragment_samples=';'.join(samples(exact_loop)),
                nearest_v_identity=nearest['identity'] if nearest else None,
                nearest_v_coverage=nearest['query_coverage'] if nearest else None,
                nearest_training_samples=';'.join(samples(nearest['records'])) if nearest else '',
                homologous_training_samples=';'.join(samples([r for h in hits for r in h['records']])))
            rows.append(row)
            detailed.append(dict(summary=row, query_full_sequence=seq, query_v_window=v_window,
                exact_full_records=exact_full, exact_v_records=exact_v, exact_loop_records=exact_loop,
                high_identity_hits=hits))
            print(json.dumps({k:v for k,v in row.items() if k!='homologous_training_samples'}), flush=True)
    pd.DataFrame(rows).to_csv(OUT / 'candidate_sequence_screen.csv', index=False)
    result = dict(training_index_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        unique_training_sequences=len(training), queries=detailed,
        method='Exhaustive local BLOSUM62 with gap-open 10, gap-extend 0.5; reported hits identity >=0.70 and query aligned coverage >=0.80.',
        caveats=['V-window and CDR3 fragments are conservatively located using Cys and [FW]GxG anchors, NOT formal IMGT numbering.',
          'Exact fragment presence is a warning; a paired clonotype requires both chains and proper CDR assignment.',
          'Absence of a full-sequence exact hit is not proof of independent clonotype because training inputs may omit coordinates.',
          'The nearest hit is maximum identity among alignments covering at least 80% of this provisional V window.',
          'No MD quality or experimental outcome is inferred from this sequence audit.'])
    (OUT / 'candidate_sequence_screen.json').write_text(json.dumps(result, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
