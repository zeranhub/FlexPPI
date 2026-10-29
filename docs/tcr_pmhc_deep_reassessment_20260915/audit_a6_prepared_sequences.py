"""Map complete first-frame PDB members to the observed 3QH3/4FTV sequences.

Run with D:/Anaconda/envs/pvb_cpu/python.exe; no trajectory data used.
"""
from pathlib import Path
import json
import math
import re
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'artifacts/tcr_pmhc_reassessment_20260915'
inventory = json.loads((BASE / 'knapp2019/prefix_inventory.json').read_text())
sources = json.loads((BASE / 'structure_sequence_audit.json').read_text())
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -.5
selected = []
for token in ['A6_TCR_run6/', 'A6_TCRpMHC_run5/']:
    selected.append(next(p for p in inventory['complete_PDB_members'] if token in p['member']))
comparisons = []
backbone_gaps = []
for entry in selected:
    atom_map = {}
    for line in (ROOT / entry['local']).read_text().splitlines():
        if line.startswith('ATOM  '):
            key = (line[21], line[22:27].strip())
            atom_map.setdefault(key, {})[line[12:16].strip()] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    prepared_chains = ['D', 'E'] if 'TCRpMHC' in entry['member'] else ['A', 'B']
    for prepared_chain, reference_chain in zip(prepared_chains, ['A', 'B']):
        prepared = [r for r in entry['CA_residues'] if r['chain'] == prepared_chain]
        for previous, current in zip(prepared, prepared[1:]):
            current_number = int(re.match(r'-?\d+', current['resid']).group())
            previous_number = int(re.match(r'-?\d+', previous['resid']).group())
            if current_number - previous_number > 1:
                c = atom_map[(prepared_chain, previous['resid'])].get('C')
                n = atom_map[(prepared_chain, current['resid'])].get('N')
                distance = math.sqrt(sum((a-b)**2 for a,b in zip(c,n))) if c and n else None
                backbone_gaps.append({'member': entry['member'], 'chain': prepared_chain,
                                      'before': previous['resid'], 'after': current['resid'],
                                      'C_to_N_distance_A': distance})
        residues = sources['3QH3']['chains'][reference_chain]
        reference = [{'resid': key, 'resname': residues[key]} for key in sorted(residues, key=int)]
        left = ''.join(seq1(r['resname']) for r in reference)
        right = ''.join(seq1(r['resname']) for r in prepared)
        alignment = aligner.align(left, right)[0]
        mapped = []
        absent = []
        for i, j in alignment.indices.T:
            if i < 0:
                continue
            row = {'reference_3QH3_chain': reference_chain, 'reference_3QH3_resid': reference[i]['resid'],
                   'reference_aa': left[i]}
            if j < 0:
                absent.append(row)
            else:
                row.update({'prepared_chain': prepared_chain, 'prepared_resid': prepared[j]['resid'],
                            'prepared_index_1based': int(j)+1, 'prepared_aa': right[j],
                            'identical': left[i] == right[j]})
                mapped.append(row)
        comparisons.append({'member': entry['member'], 'prepared_chain': prepared_chain,
                            'reference_3QH3_chain': reference_chain,
                            'prepared_length': len(prepared), 'reference_length': len(reference),
                            'reference_observed_residues_absent_in_prepared': absent,
                            'mismatches': [r for r in mapped if not r['identical']], 'mapped': mapped})
report = {'description': 'Exploratory global observed-sequence alignment; first optimal alignment only; '
                         'not a complete IMGT or experimental-region mapping. PDB members read completely '
                         'from incomplete archive, full-archive MD5 not yet checked at extraction time.',
          'comparisons': comparisons, 'prepared_author_number_gaps': backbone_gaps}
(BASE / 'knapp2019/prepared_sequence_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
for row in comparisons:
    print(json.dumps({k:v for k,v in row.items() if k != 'mapped'}))
print(json.dumps({'prepared_author_number_gaps': backbone_gaps}, indent=2))
