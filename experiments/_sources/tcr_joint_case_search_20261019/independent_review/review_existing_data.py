from pathlib import Path
import csv
import json
import math
import statistics

ROOT = Path(r'E:\Project\FlexPPI\experiments\_sources')
OLD = ROOT / 'tcr_expanded_screen_v41_20261017'
OUT = ROOT / 'tcr_joint_case_search_20261019' / 'independent_review'
OUT.mkdir(parents=True, exist_ok=True)

def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

def yes(value):
    return str(value).lower() == 'true'

def ranks(values):
    result = [0.0] * len(values)
    order = sorted(range(len(values)), key=lambda index: values[index])
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2
        for index in order[start:end]:
            result[index] = rank
        start = end
    return result

def metrics(rows, kind, prediction='full'):
    chosen = [row for row in rows if yes(row[kind + '_valid'])]
    reference = [float(row[kind + '_md_delta']) for row in chosen]
    model = [float(row[kind + '_' + prediction + '_delta']) for row in chosen]
    mae = statistics.mean(abs(a - b) for a, b in zip(reference, model))
    zero = statistics.mean(abs(a) for a in reference)
    return {
        'n': len(chosen),
        'rho': statistics.correlation(ranks(reference), ranks(model)),
        'MAE': mae, 'zero_change_MAE': zero, 'relative_MAE_gain': 1 - mae / zero,
        'same_direction_n': sum(a * b > 0 for a, b in zip(reference, model)),
        'MD_mean_a': statistics.mean(float(row[kind + '_md_a']) for row in chosen),
        'MD_mean_b': statistics.mean(float(row[kind + '_md_b']) for row in chosen),
        'model_mean_a': statistics.mean(float(row[kind + '_' + prediction + '_a']) for row in chosen),
        'model_mean_b': statistics.mean(float(row[kind + '_' + prediction + '_b']) for row in chosen)
    }

data = read_csv(OLD / 'mel5_validation' / 'corrected_pair_all_aligned_positions.csv')
context = read_csv(OLD / 'mel5_context' / 'all_natural_pMHC_context_comparison.csv')
cdrs = [row for row in data if row['role'].startswith('tcr') and row['region'].startswith('CDR')]
cdr3 = [row for row in cdrs if row['region'].startswith('CDR3')]
pmhc = [row for row in context if row['role'] in ('mhc', 'peptide') and yes(row['interface'])]
summary = {
    'method': 'Independent arithmetic using only corrected aligned position tables; no new prediction, training, or MD.',
    'CDR_all45_rmsf': metrics(cdrs, 'rmsf'),
    'CDR3_all22_rmsf': metrics(cdr3, 'rmsf'),
    'CDR3_sameAA14_rmsf': metrics([row for row in cdr3 if yes(row['same_aa'])], 'rmsf'),
    'CDR_sameAA18_hchi': metrics([row for row in cdrs if yes(row['same_aa'])], 'hchi'),
    'natural_pMHC_interface_full_hchi': metrics(pmhc, 'hchi'),
    'natural_pMHC_interface_without_TCR_hchi': metrics(pmhc, 'hchi', 'pMHC_only'),
    'natural_pMHC_interface_full_rmsf': metrics(pmhc, 'rmsf'),
    'natural_pMHC_all_hchi': metrics(context, 'hchi'),
    'MHC_key_positions': {},
    'six_CDR_model_rmsf_state_means': {}
}
sb_root = ROOT / 'tcr_case_v40_20261017' / 'sb27_sb47'
sb_replicas = read_csv(sb_root / 'all13_peptide_per_replica_labels.csv')
sb_predictions = read_csv(sb_root / 'context_inference.csv')
sb_rows = []
for position in range(1, 14):
    row = {'rmsf_valid': 'True'}
    for state, pdb in [('a', '2ak4'), ('b', '4jry')]:
        predicted = next(item for item in sb_predictions if item['pdb'] == pdb and item['input'] == 'full' and int(item['position']) == position)
        replicates = [item for item in sb_replicas if item['pdb'] == pdb and int(item['position']) == position]
        row['hchi_valid'] = 'False' if predicted['aa'] in ('G', 'A') else 'True'
        row['rmsf_md_' + state] = statistics.mean(float(item['rmsf_A']) for item in replicates)
        row['hchi_md_' + state] = statistics.mean(float(item['Hchi_nats']) for item in replicates) if yes(row['hchi_valid']) else float('nan')
        row['rmsf_full_' + state] = float(predicted['rmsf'])
        row['hchi_full_' + state] = float(predicted['Hchi'])
    for kind in ('rmsf', 'hchi'):
        row[kind + '_md_delta'] = row[kind + '_md_b'] - row[kind + '_md_a']
        row[kind + '_full_delta'] = row[kind + '_full_b'] - row[kind + '_full_a']
    sb_rows.append(row)
summary['SB27_SB47_all13_rmsf'] = metrics(sb_rows, 'rmsf')
summary['SB27_SB47_all11_hchi'] = metrics(sb_rows, 'hchi')
for region in sorted(set(row['region'] for row in cdrs)):
    region_rows = [row for row in cdrs if row['region'] == region]
    summary['six_CDR_model_rmsf_state_means'][region] = metrics(region_rows, 'rmsf')
for position in (58, 65, 66, 68, 72, 75, 76, 159):
    row = next(row for row in context if row['role'] == 'mhc' and int(row['id_a']) == position)
    summary['MHC_key_positions'][str(position)] = {
        key: row[key] for key in ('aa_a', 'interface', 'hchi_md_a', 'hchi_md_b', 'hchi_full_a', 'hchi_full_b', 'hchi_md_delta', 'hchi_full_delta', 'hchi_pMHC_only_delta', 'hchi_rep1_a', 'hchi_rep2_a', 'hchi_rep3_a', 'hchi_rep1_b', 'hchi_rep2_b', 'hchi_rep3_b')
    }

def atoms(path):
    result = []
    seen = set()
    for line in path.read_text().splitlines():
        if not line.startswith('ATOM  '):
            continue
        name = line[12:16].strip()
        if name.startswith('H') or (len(name) > 1 and name[0].isdigit() and name[1] == 'H'):
            continue
        if line[16] not in (' ', 'A'):
            continue
        key = (line[21], int(line[22:26]), line[26], name)
        if key in seen:
            continue
        seen.add(key)
        result.append({'chain': line[21], 'id': int(line[22:26]), 'res': line[17:20], 'atom': name, 'xyz': tuple(float(line[start:start + 8]) for start in (30, 38, 46))})
    return result

contacts = {'definition': 'Geometric heavy-atom distances in existing model input snapshots; <=4.0 A reports proximity, not a chemical H-bond assignment or trajectory occupancy.', 'states': {}}
beta2_contacts = {}
for pdb in ('3hg1', '4jff'):
    existing = atoms(OLD / 'mel5_validation' / (pdb + '_model_input.pdb'))
    chain_res = {}
    for atom in existing:
        chain_res.setdefault(atom['chain'], set()).add((atom['id'], atom['res']))
    state = {'chain_lengths': {key: len(value) for key, value in chain_res.items()}, 'MHC_contacts': {}}
    # The corrected table supplies role-to-chain mapping rather than guessing from letter order.
    state_field = 'a' if pdb == '3hg1' else 'b'
    role_chain = {row['role']: row['chain_' + state_field] for row in data}
    mhc_chain = role_chain['mhc']
    tcr_chains = {role_chain['tcr_alpha'], role_chain['tcr_beta']}
    beta2_ids = {int(row['id_' + state_field]) for row in cdrs if row['region'] == 'CDR2β'}
    beta2_atoms = [atom for atom in existing if atom['chain'] == role_chain['tcr_beta'] and atom['id'] in beta2_ids]
    beta2_mhc = {}
    for a in existing:
        if a['chain'] != mhc_chain:
            continue
        for b in beta2_atoms:
            distance = math.dist(a['xyz'], b['xyz'])
            if distance > 4.0:
                continue
            info = beta2_mhc.setdefault(str(a['id']), {'MHC_resname': a['res'], 'atom_pair_n_le4A': 0, 'min_distance_A': distance, 'nearest_MHC_atom': a['atom'], 'nearest_TCR_atom': b['atom'], 'nearest_TCR_input_residue': b['id'], 'beta_residues': set()})
            info['atom_pair_n_le4A'] += 1
            if distance < info['min_distance_A']:
                info['min_distance_A'] = distance
                info['nearest_MHC_atom'] = a['atom']
                info['nearest_TCR_atom'] = b['atom']
                info['nearest_TCR_input_residue'] = b['id']
            info['beta_residues'].add((b['id'], b['res']))
    for info in beta2_mhc.values():
        info['beta_residues'] = sorted(info['beta_residues'])
    beta2_contacts[pdb] = beta2_mhc
    for position in (58, 65, 68):
        first = [atom for atom in existing if atom['chain'] == mhc_chain and atom['id'] == position]
        partner_res = {}
        for second in existing:
            if second['chain'] not in tcr_chains:
                continue
            distance, nearest = min((math.dist(a['xyz'], second['xyz']), a['atom']) for a in first)
            key = (second['chain'], second['id'], second['res'])
            old = partner_res.get(key)
            record = {'distance_A': distance, 'MHC_atom': nearest, 'TCR_atom': second['atom']}
            if old is None or distance < old['distance_A']:
                partner_res[key] = record
        listed = []
        for (chain, residue_id, res), record in sorted(partner_res.items(), key=lambda item: item[1]['distance_A']):
            if record['distance_A'] > 4.0:
                continue
            paired = next((row for row in data if row['chain_' + state_field] == chain and int(row['id_' + state_field]) == residue_id), None)
            listed.append({'TCR_chain': chain, 'TCR_input_residue': residue_id, 'TCR_resname': res, 'role': paired['role'] if paired else None, 'corrected_region': paired['region'] if paired else None, **record})
        state['MHC_contacts'][str(position)] = listed
    contacts['states'][pdb] = state

beta2_union = sorted(set(beta2_contacts['3hg1']) | set(beta2_contacts['4jff']), key=int)
contacts['all_MHC_positions_contacting_CDR2beta_either_state'] = []
for residue_id in beta2_union:
    row = next(row for row in context if row['role'] == 'mhc' and row['id_a'] == residue_id)
    contacts['all_MHC_positions_contacting_CDR2beta_either_state'].append({
        'MHC_position': int(residue_id), 'aa': row['aa_a'],
        'WT_contacts': beta2_contacts['3hg1'].get(residue_id),
        'engineered_contacts': beta2_contacts['4jff'].get(residue_id),
        **{key: (float(row[key]) if row[key] else None) for key in ('hchi_md_a', 'hchi_md_b', 'hchi_full_a', 'hchi_full_b', 'hchi_md_delta', 'hchi_full_delta', 'hchi_pMHC_only_delta')}
    })

(OUT / 'independent_metric_recheck.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
(OUT / 'existing_input_contact_check.json').write_text(json.dumps(contacts, ensure_ascii=False, indent=2), encoding='utf-8')
verdict = {
    'date': '2026-10-19',
    'status': 'independent_review_completed',
    'recommended_case': {'name': 'MEL5 to alpha24beta17', 'parent_PDB': '3HG1', 'engineered_PDB': '4JFF'},
    'recommendation': 'Promote for a joint RMSF/Hchi mechanistic case; do not change manuscript yet.',
    'main_model_claim': 'One joint prediction localizes constrained TCR beta recognition-loop backbone motion and context-dependent side-chain state redistribution on the sequence-identical pMHC interface.',
    'direct_spatial_bridge': {
        'region': 'TCR CDR2beta and HLA Q72',
        'CDR2beta_model_mean_RMSF_A': [0.8961336731910705, 0.6356326818466187],
        'CDR2beta_all5_model_RMSF_decrease': True,
        'Q72_model_Hchi_nats': [1.792449712753296, 1.3750494718551636],
        'Q72_model_Hchi_change_without_TCR_nats': 0.10356616973876953,
        'Q72_MD_Hchi_nats': [1.6835638613759356, 1.2186422462639088],
        'Q72_MD_triplicate_ranges_overlap': True,
        'geometric_contact_source': 'Existing GROMACS initial model-input snapshots, not re-refined crystals or trajectory occupancy.',
        'CDR2beta_MHC_contact_union_unfiltered': ['K68', 'Q72', 'R75', 'V76'],
        'all_contact_union_observations': contacts['all_MHC_positions_contacting_CDR2beta_either_state']
    },
    'whole_region_support': {
        'all_six_CDR_RMSF': summary['CDR_all45_rmsf'],
        'complete_CDR3_RMSF': summary['CDR3_all22_rmsf'],
        'unchanged_CDR3_RMSF': summary['CDR3_sameAA14_rmsf'],
        'sequence_identical_pMHC_interface_Hchi_full': summary['natural_pMHC_interface_full_hchi'],
        'same_interface_Hchi_without_TCR': summary['natural_pMHC_interface_without_TCR_hchi']
    },
    'primary_study': {
        'title': 'T-cell Receptor Specificity Maintained by Altered Thermodynamics',
        'year': 2013, 'journal': 'Journal of Biological Chemistry',
        'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3696650/',
        'experiments': ['phage display engineering', 'SPR equilibrium and kinetic binding', 'X-ray crystallography', 'peptide alanine substitution', 'temperature-dependent SPR binding thermodynamics'],
        'SPR_KD_parent_uM': 18, 'SPR_KD_engineered_nM': 0.6,
        'MHC_contacts_parent_engineered': [56, 97],
        'peptide_contacts_parent_engineered': [39, 42],
        'CDR2beta_total_pMHC_contacts_parent_engineered': [6, 38],
        'reported_MD': False,
        'binding_entropy_interpretation': 'Favorable total binding entropy linked to solvent reorganization/expulsion; cannot be substituted with bound-state sidechain Hchi.'
    },
    'material_limits_for_review': [
        'Use corrected alignment only; old 14-point rho=0.938 is invalid.',
        'Alpha CDR2 RMSF direction differs from MD; do not claim all loops reproduced.',
        'TCR CDR same-amino-acid Hchi is weak; use sequence-identical pMHC interface for the Hchi claim.',
        'All-pMHC Hchi response is weak; interface claim must retain its complete 23-position scope.',
        'R75 has increased static contacts but MD Hchi increases and model is approximately unchanged; extra contacts do not imply lower Hchi at every site.',
        'K68 has good numerical Hchi change but should not be called a newly locked CDR2beta sidechain contact.',
        'Context removal diagnoses model input dependence and is not an apo-pMHC biological experiment.',
        'Do not equate Hchi with total binding entropy, combine RMSF and Hchi into affinity, or claim prediction of the 30000-fold SPR affinity change.'
    ],
    'figure_concept': ['structural affinity-maturation context and Q72 contact patch', 'absolute predicted RMSF across all six CDR loops, 45 positions', 'absolute predicted Hchi across all 23 valid sequence-identical pMHC interface positions', 'Q72 spatial bridge plus TCR-context removal with region-level corroboration'],
    'new_inference': False, 'new_training': False, 'new_MD': False,
    'manuscript_modified': False,
    'files': ['case_argument_review.md', 'independent_metric_recheck.json', 'existing_input_contact_check.json', 'review_existing_data.py']
}
(OUT / 'case_argument_review.json').write_text(json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2))
print(json.dumps(contacts, ensure_ascii=False, indent=2))
