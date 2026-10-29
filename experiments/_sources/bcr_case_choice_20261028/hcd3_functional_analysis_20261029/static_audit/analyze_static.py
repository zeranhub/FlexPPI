"""Exploratory complete-region structural audit of frozen A00LZ (6A77).

No mutation effects, affinity estimates, or residue importance labels are inferred.
Distances are in Angstrom and SASA is in square Angstrom.
"""
from pathlib import Path
import json, hashlib, copy
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

ROOT = Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
OUT = ROOT / 'hcd3_functional_analysis_20261029' / 'static_audit'
PDB = Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\dataset\structures\A00LZ.pdb')
OUT.mkdir(parents=True, exist_ok=True)
structure = PDBParser(QUIET=True).get_structure('6A77_frozen', str(PDB))
for residue in structure.get_residues():
    for atom in list(residue):
        if atom.element.upper() in {'H', 'D'}:
            residue.detach_child(atom.id)
residues = list(structure[0].get_residues())
keys = [(r.parent.id, str(r.id[1]) + r.id[2].strip()) for r in residues]
coords = [np.asarray([a.coord for a in r], dtype=float) for r in residues]
chains = [r.parent.id for r in residues]
chain_index = {}
for c in structure[0]:
    for i, r in enumerate(c):
        chain_index[(c.id, str(r.id[1]) + r.id[2].strip())] = i
numbering = pd.read_csv(ROOT / '6A77_antibody_numbering.csv', dtype={'author': str})
number_map = {(r.chain, r.author): r._asdict() for r in numbering.itertuples(index=False)}
pred = pd.read_csv(ROOT / '6A77_all_positions_predictions_and_MD.csv', dtype={'author': str})
pred_map = {(r.variant, r.chain, r.author): r._asdict() for r in pred.itertuples(index=False)}

def sasa_for(keep):
    st = copy.deepcopy(structure)
    for c in list(st[0]):
        if c.id not in keep:
            st[0].detach_child(c.id)
    ShrakeRupley(probe_radius=1.4, n_points=240).compute(st, level='R')
    return {(r.parent.id, str(r.id[1])+r.id[2].strip()): float(r.sasa) for r in st.get_residues()}

sasa_complex = sasa_for({'A', 'B', 'C'})
sasa_fab = sasa_for({'B', 'C'})
sasa_single = {}
for c in ['A', 'B', 'C']:
    sasa_single.update(sasa_for({c}))
rows = []
pairs = []
for i, (r, key) in enumerate(zip(residues, keys)):
    rec = {'chain': key[0], 'author': key[1], 'aa': seq1(r.resname), 'resname': r.resname,
           'chain_residue_index': chain_index[key], 'role': 'antigen' if key[0]=='A' else ('light' if key[0]=='B' else 'heavy')}
    num = number_map.get(key)
    rec['region'] = num['region'] if num else ('antigen' if key[0]=='A' else 'constant_or_unnumbered')
    rec['imgt_number'] = num['imgt_number'] if num else None
    rec['imgt_insertion'] = str(num['imgt_insertion']).strip() if num else None
    rec['scope'] = rec['role'] + '_' + rec['region']
    rec['sasa_complex'] = sasa_complex[key]
    rec['sasa_single_chain'] = sasa_single[key]
    rec['sasa_fab_no_antigen'] = sasa_fab.get(key, np.nan)
    rec['sasa_lost_on_add_antigen'] = sasa_fab.get(key, np.nan) - sasa_complex[key]
    min_d = {'antigen': np.inf, 'antibody_other_chain': np.inf, 'nonlocal_same_chain': np.inf}
    buckets = {t: {s: [] for s in min_d} for t in [4.0, 4.5, 5.0]}
    for j, key2 in enumerate(keys):
        if i == j:
            continue
        category = None
        if key[0] != 'A' and key2[0] == 'A':
            category = 'antigen'
        elif key[0] in {'B','C'} and key2[0] in {'B','C'} and key[0] != key2[0]:
            category = 'antibody_other_chain'
        elif key[0] == key2[0] and abs(chain_index[key] - chain_index[key2]) > 2:
            category = 'nonlocal_same_chain'
        if category is None:
            continue
        ds = cdist(coords[i], coords[j])
        d = float(ds.min())
        min_d[category] = min(min_d[category], d)
        for t in buckets:
            if d <= t:
                buckets[t][category].append(f'{key2[0]}{key2[1]}{seq1(residues[j].resname)}')
        if key[0]=='C' and 95 <= int(key[1]) <= 109 and d <= 5.0:
            k1,k2=np.unravel_index(np.argmin(ds), ds.shape)
            pairs.append({'chain':key[0], 'author':key[1], 'aa':seq1(r.resname),
                          'partner_chain':key2[0], 'partner_author':key2[1], 'partner_aa':seq1(residues[j].resname),
                          'category':category,'min_heavy_distance':d,
                          'closest_atom':list(r)[k1].name,'closest_partner_atom':list(residues[j])[k2].name,
                          **{f'atom_pairs_le_{t:g}':int(np.count_nonzero(ds<=t)) for t in buckets}})
    for category,d in min_d.items():
        rec[f'min_heavy_distance_{category}'] = d if np.isfinite(d) else np.nan
    for t, dct in buckets.items():
        for category, hits in dct.items():
            rec[f'n_{category}_residues_le_{t:g}'] = len(hits)
            rec[f'{category}_residues_le_{t:g}'] = ';'.join(hits)
    for variant in ['full', 'intrachain']:
        pr = pred_map[(variant, *key)]
        assert pr['aa'] == rec['aa'], (key, pr['aa'], rec['aa'])
        rec[f'{variant}_rmsf'] = pr['rmsf_pred']
        rec[f'{variant}_chi_entropy'] = pr['chi_entropy_pred'] if pr['chi_entropy_valid'] else np.nan
        if variant == 'full':
            rec['MD_rmsf'] = pr['rmsf_MD']
            rec['MD_chi_entropy'] = pr['chi_entropy_MD'] if pr['chi_entropy_valid'] else np.nan
            rec['rmsf_valid'] = pr['rmsf_valid']
            rec['chi_entropy_valid'] = pr['chi_entropy_valid']
    rows.append(rec)

all_table = pd.DataFrame(rows)
all_table.to_csv(OUT/'all_positions_structure_and_dynamics.csv',index=False,encoding='utf-8-sig')
h3=all_table[(all_table.chain=='C') & (all_table.region=='CDR3')].copy()
assert h3.author.tolist() == [str(n) for n in range(95,110)]
h3.to_csv(OUT/'complete_HCDR3_structure_and_dynamics.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(pairs).to_csv(OUT/'complete_HCDR3_contact_pairs_le_5A.csv',index=False,encoding='utf-8-sig')

summaries=[]
for scope, part in all_table.groupby('scope',sort=False):
    rec={'scope':scope,'n':len(part)}
    for c in ['MD_rmsf','full_rmsf','intrachain_rmsf','MD_chi_entropy','full_chi_entropy','intrachain_chi_entropy','sasa_complex','sasa_lost_on_add_antigen']:
        rec[f'mean_{c}']=float(part[c].mean())
    for t in [4.0,4.5,5.0]:
        rec[f'n_antigen_contact_positions_le_{t:g}']=int((part[f'n_antigen_residues_le_{t:g}']>0).sum())
    summaries.append(rec)
pd.DataFrame(summaries).to_csv(OUT/'all_regions_structure_dynamics_summary.csv',index=False,encoding='utf-8-sig')

descriptive=[]
for scope,part in [('all_six_CDR',all_table[all_table.region.str.startswith('CDR')]),('HCDR3',h3)]:
    for t in [4.0,4.5,5.0]:
        contact = part[f'n_antigen_residues_le_{t:g}']>0
        for state,mask in [('antigen_contact',contact),('no_antigen_contact',~contact)]:
            p=part[mask]
            descriptive.append({'scope':scope,'threshold':t,'contact_group':state,'n':len(p),
                                **{f'mean_{c}':float(p[c].mean()) for c in ['MD_rmsf','full_rmsf','intrachain_rmsf','MD_chi_entropy','full_chi_entropy','intrachain_chi_entropy']}})
pd.DataFrame(descriptive).to_csv(OUT/'exploratory_contact_group_descriptives.csv',index=False,encoding='utf-8-sig')
audit={'input_pdb':str(PDB),'input_pdb_sha256':hashlib.sha256(PDB.read_bytes()).hexdigest(),
       'n_residues':len(rows),'n_heavy_atoms':sum(len(c) for c in coords),
       'HCDR3_author_positions':h3.author.tolist(),'heavy_distance_thresholds_A':[4.0,4.5,5.0],
       'same_chain_exclusion':'exclude sequence index separations <=2; number by coordinate order within chain',
       'sasa_method':'Bio.PDB.SASA.ShrakeRupley, 240 sphere points, probe radius1.4A, H/D removed, default atomic radii',
       'sasa_antigen_burial':'Fab (chainsB+C) minus complex (A+B+C), identical coordinates; per-residue values, no /2 convention',
       'reference_scope':'Single bound-state MD reference, no WT-mutant or apo-holo validation',
       'analysis_status':'Exploratory complete-region descriptive structural analysis; no energetic hotspot or affinity attribution'}
(OUT/'audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
cols=['author','aa','n_antigen_residues_le_4.5','antigen_residues_le_4.5','antibody_other_chain_residues_le_4.5','sasa_complex','sasa_lost_on_add_antigen','MD_rmsf','full_rmsf','intrachain_rmsf','MD_chi_entropy','full_chi_entropy','intrachain_chi_entropy']
print(h3[cols].to_string(index=False))
print('\nALL REGIONS\n'+pd.DataFrame(summaries).to_string(index=False))
print('\nCONTACT GROUPS\n'+pd.DataFrame(descriptive).to_string(index=False))
