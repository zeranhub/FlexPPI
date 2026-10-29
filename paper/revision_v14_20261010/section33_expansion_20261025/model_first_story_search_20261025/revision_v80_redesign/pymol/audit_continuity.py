"""Audit visible backbone continuity without rebuilding or rebonding input atoms."""
from pathlib import Path
import json, math, hashlib
from collections import defaultdict

OUT=Path(__file__).resolve().parent
def distance(a,b):return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))

def main():
    report={}
    for p in (OUT/'inputs').glob('*.pdb'):
        chains=defaultdict(dict)
        for line in p.read_text().splitlines():
            if not line.startswith('ATOM'):continue
            atom=line[12:16].strip()
            if atom not in ('N','CA','C'):continue
            chain=line[21];res=int(line[22:26]);xyz=tuple(float(line[i:i+8]) for i in (30,38,46))
            chains[chain].setdefault(res,{})[atom]=xyz
        entry={}
        for ch,residues in chains.items():
            nums=sorted(residues);ca=[];peptide=[];missing=[]
            for n in nums:
                if not {'N','CA','C'} <= set(residues[n]):missing.append(n)
            for a,b in zip(nums,nums[1:]):
                if 'CA' in residues[a] and 'CA' in residues[b]:ca.append((a,b,distance(residues[a]['CA'],residues[b]['CA'])))
                if 'C' in residues[a] and 'N' in residues[b]:peptide.append((a,b,distance(residues[a]['C'],residues[b]['N'])))
            entry[ch]={'n_residues':len(nums),'first':nums[0],'last':nums[-1],
                'missing_numbers':[n for n in range(nums[0],nums[-1]+1) if n not in residues],
                'missing_backbone_atoms_at_residues':missing,
                'max_sequential_CA_distance_A':max(v[2] for v in ca),
                'max_sequential_C_N_distance_A':max(v[2] for v in peptide),
                'CA_gaps_over_4_5_A':[v for v in ca if v[2]>4.5],
                'peptide_gaps_over_1_8_A':[v for v in peptide if v[2]>1.8]}
        report[p.name]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'chains':entry}
    (OUT/'backbone_continuity_QA.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
