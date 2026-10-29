"""Read RCSB mmCIFs in memory and audit auth-numbered TCR chains."""
import io
import json
import urllib.request
import hashlib
from pathlib import Path
from collections import defaultdict
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1
import numpy as np

PDBS = ['8TW4', '8TW6', '9C3E', '9BBC', '3QH3', '1AO7', '1QRN', '4FTV', '3QEU', '3QDG', '5IW1', '5IVX', '1KGC', '1MI5', '1TCR', '2CKB', '2BNU', '2BNR', '2BNQ']
OUT = Path('artifacts/tcr_pmhc_reassessment_20260915')
(OUT / 'sources').mkdir(parents=True, exist_ok=True)
structures = {}
for pdb in PDBS:
    url = 'https://files.rcsb.org/download/' + pdb + '.cif'
    source_path = OUT / 'sources' / (pdb + '.cif')
    data = source_path.read_bytes() if source_path.exists() else urllib.request.urlopen(url, timeout=60).read()
    source_path.write_bytes(data)
    text = data.decode()
    cif = MMCIF2Dict(io.StringIO(text))
    entities = dict(zip(cif['_entity_poly.entity_id'], cif['_entity_poly.pdbx_seq_one_letter_code_can']))
    chains = defaultdict(dict)
    backbone = defaultdict(lambda: defaultdict(dict))
    chain_entity = {}
    cols = ['group_PDB', 'label_atom_id', 'auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code', 'label_comp_id', 'label_entity_id', 'pdbx_PDB_model_num', 'Cartn_x', 'Cartn_y', 'Cartn_z']
    for row in zip(*[cif['_atom_site.' + col] for col in cols]):
        group, atom, chain, resi, ins, aa, entity, model, x, y, z = row
        if group != 'ATOM' or atom not in ['N','CA','C','O'] or model != '1':
            continue
        key = resi + ('' if ins in ['?', '.'] else ins)
        backbone[chain][key][atom] = [float(x),float(y),float(z)]
        if atom != 'CA':
            continue
        chains[chain][key] = aa
        chain_entity[chain] = entity
    structures[pdb] = {'chains': dict(chains), 'backbone': {c:dict(v) for c,v in backbone.items()}, 'entity_sequences': entities, 'chain_entity': chain_entity, 'source_url':url, 'sha256':hashlib.sha256(data).hexdigest(), 'missing_backbone_atoms':{c:{k:[a for a in ['N','CA','C','O'] if a not in atoms] for k,atoms in v.items() if len(atoms)<4} for c,v in backbone.items()}}
    print(json.dumps({'pdb': pdb, 'method': cif.get('_exptl.method'), 'resolution': cif.get('_em_3d_reconstruction.resolution', cif.get('_refine.ls_d_res_high')), 'chains': {c:{'modeled_n':len(v), 'first':list(v)[0], 'last':list(v)[-1], 'seqres_n':len(entities[chain_entity[c]].replace('\n',''))} for c,v in chains.items()}}))

comparisons = []
for other in ['8TW4', '8TW6', '9BBC']:
    for chain in ['A', 'B', 'D', 'E', 'F', 'G', 'X', 'Y']:
        left = structures[other]['chains'].get(chain, {})
        right = structures['9C3E']['chains'].get(chain, {})
        common = sorted(set(left) & set(right), key=lambda x: int(x))
        mismatch = [(k,left[k],right[k]) for k in common if left[k] != right[k]]
        comparisons.append({'comparison':other+' vs 9C3E', 'chain':chain, 'n_common_auth':len(common), 'common_auth_ids':common, 'mismatches':mismatch, 'left_only':sorted(set(left)-set(right),key=int), 'right_only':sorted(set(right)-set(left),key=int)})

def fit_rmsd(left_xyz, right_xyz):
    left = np.array(left_xyz)
    right = np.array(right_xyz)
    left0 = left - left.mean(axis=0)
    right0 = right - right.mean(axis=0)
    u, s, vh = np.linalg.svd(left0.T @ right0)
    correction = np.eye(3)
    correction[2,2] = np.linalg.det(u @ vh)
    rotation = u @ correction @ vh
    moved = left0 @ rotation + right.mean(axis=0)
    distances = np.linalg.norm(moved-right,axis=1)
    return float(np.sqrt(np.mean(distances**2))), distances.tolist()

geometry = []
for other in ['8TW4','8TW6','9BBC']:
    for label, selector in [
        ('TCR_alpha_beta_all_observed_common', lambda c,k: c in ['A','B']),
        ('TCR_alpha_all',lambda c,k:c=='A'),
        ('TCR_beta_all',lambda c,k:c=='B'),
        ('alpha_variable_21_130_exploratory',lambda c,k:c=='A' and 21<=int(k)<=130),
        ('alpha_constant_140_230_exploratory',lambda c,k:c=='A' and 140<=int(k)<=230),
        ('beta_variable_22_140_exploratory',lambda c,k:c=='B' and 22<=int(k)<=140),
        ('beta_constant_150_270_exploratory',lambda c,k:c=='B' and 150<=int(k)<=270),
    ]:
        pairs=[]
        for chain in ['A','B']:
            left=structures[other]['chains'][chain]
            right=structures['9C3E']['chains'][chain]
            for key in sorted(set(left)&set(right),key=int):
                if left[key]==right[key] and selector(chain,key):
                    pairs.append((chain,key))
        left_xyz=[structures[other]['backbone'][c][k]['CA'] for c,k in pairs]
        right_xyz=[structures['9C3E']['backbone'][c][k]['CA'] for c,k in pairs]
        rmsd, distances = fit_rmsd(left_xyz,right_xyz)
        row={'comparison':other+' vs 9C3E','selection':label,'n_common':len(pairs),'fitted_CA_RMSD_A':rmsd,'residue_pairs':pairs,'per_residue_distance_after_own_fit_A':distances}
        geometry.append(row)
        print(json.dumps({k:v for k,v in row.items() if k not in ['residue_pairs','per_residue_distance_after_own_fit_A']}))

(OUT / 'notti_matched_residue_comparisons.json').write_text(json.dumps(comparisons, indent=2), encoding='utf-8')
(OUT / 'notti_geometry_audit.json').write_text(json.dumps(geometry, indent=2), encoding='utf-8')
(OUT / 'structure_sequence_audit.json').write_text(json.dumps(structures, indent=2), encoding='utf-8')

aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -2
aligner.open_gap_score = -5
aligner.extend_gap_score = -.1
classical = []
for label, apo, holo, chain_pairs in [
    ('A6_nativeTax','3QH3','1AO7',[('A','D'),('B','E')]),
    ('A6_P6A','3QH3','1QRN',[('A','D'),('B','E')]),
    ('A6_c134','3QH3','4FTV',[('A','D'),('B','E')]),
    ('DMF5','3QEU','3QDG',[('D','D'),('E','E')]),
    ('B4.2.3','5IW1','5IVX',[('A','E'),('B','F')]),
    ('LC13','1KGC','1MI5',[('D','D'),('E','E')]),
    ('2C','1TCR','2CKB',[('A','A'),('B','B')]),
    ('1G4_native_C9','2BNU','2BNR',[('A','D'),('B','E')]),
    ('1G4_C9V','2BNU','2BNQ',[('A','D'),('B','E')]),
]:
    for apo_chain, holo_chain in chain_pairs:
        if apo_chain not in structures[apo]['chains'] or holo_chain not in structures[holo]['chains']:
            print(json.dumps({'pair_chain_missing':label,'apo_chains':list(structures[apo]['chains']),'holo_chains':list(structures[holo]['chains'])}))
            continue
        left = list(structures[apo]['chains'][apo_chain].items())
        right = list(structures[holo]['chains'][holo_chain].items())
        left_seq = ''.join(seq1(aa) for _,aa in left)
        right_seq = ''.join(seq1(aa) for _,aa in right)
        alignment = aligner.align(left_seq,right_seq)[0]
        matched=[]
        mismatches=[]
        for leftblock,rightblock in zip(*alignment.aligned):
            for li,ri in zip(range(*leftblock),range(*rightblock)):
                lk,la=left[li]
                rk,ra=right[ri]
                row={'apo_resi':lk,'holo_resi':rk,'apo_aa':la,'holo_aa':ra}
                if la == ra:
                    matched.append(row)
                else:
                    mismatches.append(row)
        left_xyz=[structures[apo]['backbone'][apo_chain][r['apo_resi']]['CA'] for r in matched]
        right_xyz=[structures[holo]['backbone'][holo_chain][r['holo_resi']]['CA'] for r in matched]
        rmsd, distances = fit_rmsd(left_xyz,right_xyz)
        row={'case':label,'apo':apo,'holo':holo,'apo_chain':apo_chain,'holo_chain':holo_chain,'apo_observed_n':len(left),'holo_observed_n':len(right),'n_matched_identical':len(matched),'mismatches':mismatches,'mapping':matched,'CA_RMSD_over_identical_residues_A':rmsd,'alignment':str(alignment),'warning':'Global observed-sequence alignment; missing segments are not completed. Review mapping before modeling.'}
        classical.append(row)
        print(json.dumps({k:v for k,v in row.items() if k not in ['mapping','alignment','warning']}))
(OUT / 'classical_apo_holo_sequence_mapping.json').write_text(json.dumps(classical, indent=2), encoding='utf-8')
