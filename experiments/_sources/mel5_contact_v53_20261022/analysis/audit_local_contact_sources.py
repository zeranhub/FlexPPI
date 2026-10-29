"""Recompute descriptive input-structure contact geometry; no model/MD runs."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd

BASE = Path(r'E:\Project\FlexPPI\experiments\_sources')
ROOT = BASE / 'mel5_contact_v53_20261022' / 'analysis'
OLD = BASE / 'mel5_joint_case_v49_20261019' / 'figure4'
V51 = BASE / 'mel5_interface_v51_20261022' / 'analysis'
SOURCES = {
    'WT_input_PDB': OLD / 'source_archive/3hg1_model_input.pdb',
    'mature_input_PDB': OLD / 'source_archive/4jff_model_input.pdb',
    'aligned_model_outputs': OLD / 'source_archive/corrected_pair_all_aligned_positions.csv',
    'all23_plotted_Hchi': OLD / 'source_plotted_interface_Hchi.csv',
    'all45_plotted_CDR_RMSF': OLD / 'source_plotted_CDR_RMSF.csv',
    'v51_local_network': V51 / 'local_Q72_betaCDR2_neighbor_network.csv',
}
ROOT.mkdir(parents=True, exist_ok=True)

def parse(path):
    atoms = {}
    for line in path.read_text().splitlines():
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        element = line[76:78].strip() or line[12:16].strip()[0]
        if element in ('H', 'D') or line[16:17] not in (' ', 'A'):
            continue
        key = (line[21], int(line[22:26]))
        atoms.setdefault(key, {})[line[12:16].strip()] = np.array([
            float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return atoms

atoms = [parse(SOURCES['WT_input_PDB']), parse(SOURCES['mature_input_PDB'])]
aligned = pd.read_csv(SOURCES['aligned_model_outputs'])

def calculate_pair(p, t):
    rec = {'HLA_common_position': int(p.position_a), 'HLA_WT_aa': p.aa_a,
           'HLA_mature_aa': p.aa_b, 'TCR_common_position': int(t.position_a),
           'TCR_WT_aa': t.aa_a, 'TCR_mature_aa': t.aa_b,
           'TCR_region': t.region, 'HLA_Hchi_WT_nats': p.hchi_full_a,
           'HLA_Hchi_mature_nats': p.hchi_full_b,
           'TCR_RMSF_WT_A': t.rmsf_full_a, 'TCR_RMSF_mature_A': t.rmsf_full_b}
    for i, (state, suffix) in enumerate([('WT', 'a'), ('mature', 'b')]):
        p_atoms = atoms[i][(p['chain_'+suffix], int(p['id_'+suffix]))]
        t_atoms = atoms[i][(t['chain_'+suffix], int(t['id_'+suffix]))]
        distances = [(float(np.linalg.norm(x-y)), pn, tn)
                     for pn,x in p_atoms.items() for tn,y in t_atoms.items()]
        d,pn,tn = min(distances)
        n = sum(datum[0] <= 4.0 for datum in distances)
        rec.update({state+'_minimum_heavy_distance_A': d,
                    state+'_minimum_HLA_atom': pn, state+'_minimum_TCR_atom': tn,
                    state+'_heavy_atom_pair_count_le4A': n,
                    state+'_residue_pair_contact_le4A': bool(n)})
    rec['delta_minimum_heavy_distance_A'] = rec['mature_minimum_heavy_distance_A']-rec['WT_minimum_heavy_distance_A']
    rec['delta_HLA_Hchi_nats'] = rec['HLA_Hchi_mature_nats']-rec['HLA_Hchi_WT_nats']
    rec['delta_TCR_RMSF_A'] = rec['TCR_RMSF_mature_A']-rec['TCR_RMSF_WT_A']
    a,b = rec['WT_residue_pair_contact_le4A'], rec['mature_residue_pair_contact_le4A']
    rec['contact_status'] = 'retained' if a and b else ('lost' if a else ('gained' if b else 'absent_both'))
    return rec

hla = aligned[(aligned.role=='mhc') & aligned.position_a.isin([65,66,68,72,75,76])]
beta = aligned[(aligned.chain_a=='E') & aligned.position_a.between(48,56)]
local = pd.DataFrame([calculate_pair(p,t) for _,p in hla.iterrows() for _,t in beta.iterrows()])
local.to_csv(ROOT/'local_HLA6_beta48_56_all54_pairs.csv', index=False)
local[local.contact_status!='absent_both'].to_csv(ROOT/'local_HLA6_beta48_56_contact_edges.csv', index=False)
q72 = local[(local.HLA_common_position==72) & local.TCR_common_position.between(50,54)]
assert len(q72)==5
q72.to_csv(ROOT/'Q72_betaCDR2_all5_source.csv', index=False)

# Compare every existing v51 entry in this complete window, not just selected edges.
v51 = pd.read_csv(SOURCES['v51_local_network'])
checked = 0
for _, r in local.iterrows():
    matches = v51[(v51.pMHC_common_position==r.HLA_common_position) &
                  (v51.TCR_common_position==r.TCR_common_position)]
    if matches.empty:
        assert min(r.WT_minimum_heavy_distance_A,r.mature_minimum_heavy_distance_A)>6
        continue
    s=matches.iloc[0]
    for st in ['WT','mature']:
        assert abs(r[st+'_minimum_heavy_distance_A']-s[st+'_distance_A'])<1e-12
        assert r[st+'_heavy_atom_pair_count_le4A']==s[st+'_n_pairs_le4A']
    checked += 1

hchi = pd.read_csv(SOURCES['all23_plotted_Hchi'])
assert len(hchi)==23
hchi['delta_Hchi_nats']=hchi.hchi_full_b-hchi.hchi_full_a
hchi['rank_increasing_delta']=hchi.delta_Hchi_nats.rank(method='min').astype(int)
hchi.to_csv(ROOT/'all23_interface_Hchi_source.csv',index=False)
hchi[(hchi.chain_a=='A') & hchi.position_a.isin([65,66,68,72,75,76])].to_csv(ROOT/'HLA6_Hchi_source.csv',index=False)
cdr=pd.read_csv(SOURCES['all45_plotted_CDR_RMSF'])
assert len(cdr)==45
cdr[cdr.region=='CDR2β'].to_csv(ROOT/'betaCDR2_RMSF_all5_source.csv',index=False)

definitions = {
    'status':'passed',
    'input_scope':'Archived simulation-prepared starting/model-input structures from 3HG1 and 4JFF; not newly retrieved raw crystals, MD means, or predicted coordinates.',
    'WT_label':'MEL5; mature label alpha24beta17.',
    'numbering':'Common MEL5 numbering throughout. Mature TCR chain E input residue is common+1; HLA A numbering unchanged.',
    'distance_definition':'Minimum Euclidean distance across every nonhydrogen atom in the two residues, including backbone and sidechain. The minimizing atom identities may differ across states. Do not label this as a fixed OE1-N or NE2-O distance.',
    'contact_definition':'A residue pair is a contact when at least one pair of its heavy atoms is at distance <=4.0 A in that one input structure. Gained/lost/retained refer only to this binary geometrical definition, not dynamics or binding energy.',
    'atom_pair_count_definition':'Number of distinct inter-residue heavy-atom pairs <=4.0 A in the one input structure. Not independent contacts, bond number, hydrogen-bond occupancy, or trajectory frame fraction. Mutations alter atom counts, so contact-count magnitudes are not comparable as an affinity score.',
    'curve_definition':'RMSF and Hchi are the existing per-residue Full joint model output of each complete bound complex. Delta = mature minus MEL5. Hchi nats; RMSF A. No new inference was performed.',
    'all23_scope':'Exactly all23 valid same-amino-acid pMHC interface residues in the archived plotted source:18 HLA +5 peptide; not all HLA residues or all interface atom pairs.',
    'local_scope':'Explicit local window of six HLA sites65,66,68,72,75,76 and aligned TCR beta48-56; complete54pair grid. Q72/CDR2beta primary table contains all50-54 five sites without a proximity filter. This is a local mechanism view, not a claim of prospective residue selection.',
    'causal_scope':'Diagram edges show geometric proximity and node values show model response. They do not establish a directional Hchi-to-RMSF influence, energy transfer, or predict the measured affinity ratio.',
    'source_sha256':{k:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in SOURCES.items()},
    'v51_pairs_recomputed_exactly':checked,
    'local_contact_edge_count':int(sum(local.contact_status!='absent_both')),
    'all23_decreases':int(sum(hchi.delta_Hchi_nats<0)),
}
# Reproducible check of the transparent local display-selection rule.
all12 = aligned[(aligned.role=='mhc') & aligned.position_a.between(65,76)]
selection=[]
for _,p in all12.iterrows():
    z={'common_position':int(p.position_a),'aa':p.aa_a,'hchi_valid':bool(p.hchi_valid)}
    for i,(state,suffix) in enumerate([('WT','a'),('mature','b')]):
        pp=np.array(list(atoms[i][(p['chain_'+suffix],int(p['id_'+suffix]))].values()))
        distances=[]
        for _,t in beta.iterrows():
            tt=np.array(list(atoms[i][(t['chain_'+suffix],int(t['id_'+suffix]))].values()))
            d=float(np.linalg.norm(pp[:,None]-tt[None],axis=2).min())
            distances.append((d,int(t.position_a)))
        d,pos=min(distances)
        z[state+'_closest_beta_residue']=pos
        z[state+'_minimum_to_beta48_56_A']=d
    z['eligible']=z['hchi_valid'] and min(z['WT_minimum_to_beta48_56_A'],z['mature_minimum_to_beta48_56_A'])<=6
    selection.append(z)
selection=pd.DataFrame(selection)
assert selection[selection.eligible].common_position.tolist()==[65,66,68,72,75,76]
selection.to_csv(ROOT/'HLA65_76_selection_rule_audit_all12.csv',index=False)
definitions['selection_rule']={
    'status':'verified_exactly_six',
    'region':'All aligned HLA chain A residues65-76 and aligned TCRbeta48-56',
    'eligibility':'Hchi valid and minimum heavy-atom distance to any beta48-56 residue <=6.0A in either archived structure',
    'selected_positions':selection[selection.eligible].common_position.tolist(),
    'excluded_positions':selection[~selection.eligible].to_dict('records'),
    'n_HLA_candidates':len(selection),'n_beta_candidates':len(beta),
    'full_rule_table':str(ROOT/'HLA65_76_selection_rule_audit_all12.csv'),
    'scope':'Rule is a transparent geometric definition of a local display window, not a preregistered/prospective selection criterion.'}
(ROOT/'data_scope_audit.json').write_text(json.dumps(definitions,indent=2,ensure_ascii=False),encoding='utf8')
print(q72[['TCR_common_position','TCR_WT_aa','TCR_mature_aa','WT_minimum_heavy_distance_A','mature_minimum_heavy_distance_A','WT_heavy_atom_pair_count_le4A','mature_heavy_atom_pair_count_le4A','contact_status','TCR_RMSF_WT_A','TCR_RMSF_mature_A']].to_string(index=False))
print(local[local.contact_status!='absent_both'][['HLA_common_position','TCR_common_position','contact_status','WT_minimum_heavy_distance_A','mature_minimum_heavy_distance_A']].to_string(index=False))
