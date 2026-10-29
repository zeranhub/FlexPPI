from pathlib import Path
import json,csv,itertools,sys
from collections import defaultdict,Counter
from difflib import SequenceMatcher
ROOT=Path(__file__).parent
BASE=ROOT.parents[1]
FAMILY=BASE/'family_folds/splits/fold_0_test.jsonl'
rows=[json.loads(x) for x in FAMILY.read_text().splitlines() if x.strip()]
ids=sorted(r['sample_id'][7:] for r in rows if r['dataset']=='itrace')
meta={r['pdb_id']:r for r in json.loads(Path('E:/data/ITrace/tables/manifest.json').read_text())['trajectories'] if r['traj_id'].endswith('run1')}
roles=json.loads(Path('E:/Project/FlexPPI/output/itrace_benchmark/raw/itrace_chain_roles.json').read_text())
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
AA.update(HID='H',HIE='H',HIP='H',HSD='H',HSE='H',HSP='H',CYX='C',MSE='M')
inventory={}
for pdb in ids:
    path=Path(f'E:/data/ITrace/trajectories/{pdb}_run1_topology.pdb')
    chains=defaultdict(list);seen=set()
    for line in path.read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')) and line[12:16].strip()=='CA':
            chain,resid,resn=line[21:22],line[22:27].strip(),line[17:20].strip()
            if (chain,resid) in seen:continue
            seen.add((chain,resid));chains[chain].append(dict(id=resid,aa=AA.get(resn,'X')))
    byrole={roles[pdb].get(chain,'unknown'):dict(chain=chain,sequence=''.join(r['aa'] for r in residues),residues=residues) for chain,residues in chains.items()}
    inventory[pdb]=dict(pdb_id=pdb,roles=byrole,metadata=meta[pdb],structure=str(path))
(ROOT/'sequence_inventory.json').write_text(json.dumps(inventory,ensure_ascii=False,indent=2),encoding='utf-8')
pairrows=[]
def compare(a,b):
    if a==b:return 1.0,0,len(a)
    if len(a)==len(b):
        d=sum(x!=y for x,y in zip(a,b))
        return 1-d/len(a),d,len(a)
    m=SequenceMatcher(None,a,b,autojunk=False)
    matches=sum(x.size for x in m.get_matching_blocks())
    return matches/min(len(a),len(b)),max(len(a),len(b))-matches,matches
for aid,bid in itertools.combinations(ids,2):
    a,b=inventory[aid]['roles'],inventory[bid]['roles']
    common=sorted(set(a)&set(b))
    if not all(r in common for r in ['mhc','peptide','tcr_alpha','tcr_beta']):continue
    scores={r:compare(a[r]['sequence'],b[r]['sequence']) for r in common}
    mhc_match=scores['mhc'][0]>=.995 and scores['mhc'][1]<=3
    tcr_match=all(scores[r][0]>=.995 and scores[r][1]<=3 for r in ['tcr_alpha','tcr_beta'])
    pep_a,pep_b=a['peptide']['sequence'],b['peptide']['sequence']
    same_peptide=pep_a==pep_b
    tcr_d=sum(scores[r][1] for r in ['tcr_alpha','tcr_beta'])
    category=None
    if mhc_match and tcr_match and not same_peptide:category='same_TCR_MHC_peptide_change'
    elif mhc_match and same_peptide and 0<tcr_d<=6:category='same_pMHC_small_TCR_change'
    elif tcr_match and same_peptide and 0<scores['mhc'][1]<=5:category='same_TCR_peptide_small_MHC_change'
    elif mhc_match and same_peptide and tcr_d>6:category='same_pMHC_different_TCR'
    if not category:continue
    row=dict(a=aid,b=bid,category=category,peptide_a=pep_a,peptide_b=pep_b,peptide_diff=scores['peptide'][1],mhc_diff=scores['mhc'][1],tcr_alpha_diff=scores['tcr_alpha'][1],tcr_beta_diff=scores['tcr_beta'][1],hla_a=meta[aid]['hla_allele'],hla_b=meta[bid]['hla_allele'],title_a=meta[aid]['title'],title_b=meta[bid]['title'],doi_a=meta[aid].get('citation',{}).get('doi',''),doi_b=meta[bid].get('citation',{}).get('doi',''))
    pairrows.append(row)
with (ROOT/'biological_candidate_pairs.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(pairrows[0]));w.writeheader();w.writerows(pairrows)
print('inventory',len(inventory),'pairs',len(pairrows),Counter(r['category'] for r in pairrows))
for row in pairrows:
    if row['category']=='same_TCR_MHC_peptide_change' or {row['a'],row['b']} & {'1ao7','3d39','3d3v','2bnq','2bnr'}:print(json.dumps(row,ensure_ascii=False))
