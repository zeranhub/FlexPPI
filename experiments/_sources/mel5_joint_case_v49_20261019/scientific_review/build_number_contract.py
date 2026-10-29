from pathlib import Path
import csv
import json
import statistics
import hashlib

ROOT = Path(r'E:\Project\FlexPPI\experiments\_sources')
OLD = ROOT / 'tcr_expanded_screen_v41_20261017'
SEARCH = ROOT / 'tcr_joint_case_search_20261019'
OUT = ROOT / 'mel5_joint_case_v49_20261019' / 'scientific_review'
OUT.mkdir(parents=True, exist_ok=True)
PAIR = OLD / 'mel5_validation' / 'corrected_pair_all_aligned_positions.csv'
CONTEXT = OLD / 'mel5_context' / 'all_natural_pMHC_context_comparison.csv'
TRAJ = SEARCH / 'residue_story' / 'trajectory_link'

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def yes(value):
    return value == 'True'

def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    result = [0.] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        mean = (start + 1 + end) / 2
        for index in order[start:end]:
            result[index] = mean
        start = end
    return result

def metrics(rows, kind, pred='full'):
    rows = [r for r in rows if yes(r[kind+'_valid'])]
    x = [float(r[kind+'_md_delta']) for r in rows]
    y = [float(r[kind+'_'+pred+'_delta']) for r in rows]
    return {'n': len(rows), 'Spearman_change': statistics.correlation(ranks(x), ranks(y)),
            'same_direction_n': sum(a*b > 0 for a,b in zip(x,y)),
            'MAE_change': statistics.mean(abs(a-b) for a,b in zip(x,y)),
            'model_mean_WT': statistics.mean(float(r[kind+'_'+pred+'_a']) for r in rows),
            'model_mean_engineered': statistics.mean(float(r[kind+'_'+pred+'_b']) for r in rows),
            'MD_mean_WT': statistics.mean(float(r[kind+'_md_a']) for r in rows),
            'MD_mean_engineered': statistics.mean(float(r[kind+'_md_b']) for r in rows)}

rows, context = read(PAIR), read(CONTEXT)
cdr = [r for r in rows if r['role'].startswith('tcr') and r['region'].startswith('CDR')]
beta2 = [r for r in cdr if r['region'] == 'CDR2β']
interface = [r for r in rows if r['role'] in ('mhc','peptide') and yes(r['interface']) and yes(r['hchi_valid'])]
interface_context = [r for r in context if r['role'] in ('mhc','peptide') and yes(r['interface']) and yes(r['hchi_valid'])]
assert len(cdr) == 45 and len(beta2) == 5 and len(interface) == 23 and len(interface_context) == 23
assert all(yes(r['same_aa']) for r in interface)
assert all(float(r['rmsf_full_delta']) < 0 for r in beta2)
q72 = next(r for r in rows if r['role'] == 'mhc' and r['id_a'] == '72')
gly54 = next(r for r in beta2 if r['id_a'] == '54')
assert gly54['aa_a'] == gly54['aa_b'] == 'G' and gly54['id_b'] == '55'

rep = read(TRAJ/'per_replica_Q72_beta_contact_and_RMSF.csv')
dual = read(TRAJ/'Q72_two_end_backbone_Hbond_same_frame_per_replica.csv')
windows = read(TRAJ/'within_replica_20ns_descriptive_correlations.csv')
traj = {}
for pdb in ('3hg1','4jff'):
    selected = [r for r in rep if r['pdb']==pdb and r['beta_common_position']=='54']
    assert len(selected)==3
    traj[pdb] = {key: [float(r[key]) for r in selected] for key in (
        'Q72_Hchi_nats','beta_CA_RMSF_A','min_heavy_contact_le4_A_fraction',
        'Q72_NE2_beta_GlyO_hbond_fraction','Q72_OE1_beta_Ile_or_PheN_hbond_fraction')}
    traj[pdb]['dual_same_frame_hbond_geometry_fraction'] = [float(r['both_same_frame_Hbond_fraction']) for r in dual if r['pdb']==pdb]
    traj[pdb]['three_replica_arithmetic_means'] = {k: statistics.mean(v) for k,v in traj[pdb].items()}
    traj[pdb]['within_replica_20ns_Hchi_Gly_RMSF_Spearman'] = [float(r['Spearman_descriptive']) for r in windows if r['pdb']==pdb and r['beta_common_position']=='54' and r['x']=='Q72_Hchi_nats' and r['y']=='beta_CA_RMSF_A']
assert max(traj['4jff']['beta_CA_RMSF_A']) < min(traj['3hg1']['beta_CA_RMSF_A'])

fieldnames = list(rows[0])
for name, chosen in [('all45_CDR_source_for_review.csv',cdr),('all23_pMHC_Hchi_interface_source_for_review.csv',interface)]:
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(chosen)

files = [PAIR, CONTEXT, OLD/'mel5_context/mel5_context_inference_audit.json',
         OLD/'mel5_validation/corrected_alignment_and_CDR_audit.json',
         OLD/'biology_sources/mel5_thermodynamics_2013.md',
         TRAJ/'analysis_audit.json', TRAJ/'per_replica_Q72_beta_contact_and_RMSF.csv',
         TRAJ/'three_replica_summary.csv', TRAJ/'within_replica_20ns_descriptive_correlations.csv',
         TRAJ/'Q72_two_end_backbone_Hbond_same_frame_per_replica.csv',
         SEARCH/'residue_story/independent_review/final_evidence_chain_review.json',
         ROOT/'tcr_case_v40_20261017/family_predictions/full/itrace_3hg1.npz',
         ROOT/'tcr_case_v40_20261017/family_predictions/full/itrace_4jff.npz']

contract = {
    'status':'independent_number_and_mechanism_contract_ready', 'date':'2026-10-19',
    'scope':'Existing model predictions and already available ITrace trajectories; no new inference/training/MD. Contract is not final manuscript or figure approval.',
    'state_order':['MEL5 / 3HG1','alpha24beta17 / 4JFF'],
    'unit_contract':{'RMSF':'angstrom; absolute output, not change, unless delta is explicitly labelled','Hchi':'nats; entropy of fixed sidechain-chi joint states, not thermodynamic binding entropy','delta':'engineered minus WT'},
    'prediction_source_of_record':str(PAIR),
    'prediction_roundoff_note':'The original corrected table is authoritative for figure absolute predictions. The context experiment rerun reproduces full outputs within 2e-5 and can differ in last decimals; it is valid only for context sensitivity analyses. Do not silently replace the authoritative original curve with the rerun.',
    'key_model_values':{
        'HLA_Q72_Hchi_WT_engineered_nats':[float(q72['hchi_full_a']),float(q72['hchi_full_b'])],
        'unchanged_beta_Gly54_RMSF_WT_engineered_A':[float(gly54['rmsf_full_a']),float(gly54['rmsf_full_b'])],
        'all5_CDR2beta_RMSF_WT_engineered_mean_A':[statistics.mean(float(r['rmsf_full_a']) for r in beta2),statistics.mean(float(r['rmsf_full_b']) for r in beta2)],
        'all5_CDR2beta_model_RMSF_decreases':True,
        'Q72_model_change_without_TCR_nats':float(next(r for r in context if r['role']=='mhc' and r['id_a']=='72')['hchi_pMHC_only_delta'])},
    'scope_contract':{
        'RMSF':'All six projected CDR loops, all 45 aligned positions, including substitutions. No best-loop/same-AA filtering for the overall plot.',
        'Hchi':'All 23 valid chi-bearing positions in the sequence-identical pMHC interface, from the pre-existing geometric union interface definition (8 Å in either input). This is 18 HLA plus 5 peptide residues, not the entire pMHC.',
        'CDR2_focus':'All five aligned beta positions are retained; Q72/Gly54 is a structural focus within complete regional data.',
        'sidechain_comparability':'TCR substitutions between different residue types do not permit a simple same-chi-state Hchi comparison; use unchanged pMHC for the Hchi response.'},
    'all45_CDR_RMSF':metrics(cdr,'rmsf'),
    'all23_pMHC_interface_Hchi_original_predictions':metrics(interface,'hchi'),
    'all23_pMHC_interface_Hchi_context_full_reproduction':metrics(interface_context,'hchi'),
    'all23_pMHC_interface_Hchi_without_TCR':metrics(interface_context,'hchi','pMHC_only'),
    'all6_CDR_means':{region:metrics([r for r in cdr if r['region']==region],'rmsf') for region in sorted(set(r['region'] for r in cdr))},
    'all23_site_identities':[{'role':r['role'],'WT_input_id':r['id_a'],'engineered_input_id':r['id_b'],'aa':r['aa_a']} for r in interface],
    'all5_beta_CDR2_residue_contract':[{'paper_compatible_beta_position':int(r['id_a']),'WT_input_id':r['id_a'],'engineered_input_id':r['id_b'],'AA_WT':r['aa_a'],'AA_engineered':r['aa_b'],'RMSF_WT_A':float(r['rmsf_full_a']),'RMSF_engineered_A':float(r['rmsf_full_b'])} for r in beta2],
    'numbering_contract':'Primary article explicitly names beta Ile53 -> Phe53. Five-position publication-compatible beta50–54 is the continuous aligned SVGIG -> WGPFG patch. Engineered model input is beta51–55 due an extra N-terminal residue; common beta Gly54 is WT input E54 and engineered E55. Never conflate engineered input E54 Phe with paper/common beta54 Gly.',
    'trajectory_independent_support':traj,
    'trajectory_protocol':{'sampling':'All saved frames from 0–200 ns, 1001–1003 per trajectory, all three independent replicas per state. Replicas, not frames, are replication units.','RMSF':'Entire beta chain C-alpha aligned per replica to its first frame; per-residue C-alpha RMSF.','Hchi':'3 fixed bins per chi with natural log; Q72 has 3 chis; all original labels reproduced.','contact':'Q72 sidechain heavy atoms to target beta residue heavy atoms <=4 Å, periodic minimum image.','Hbond_geometry':'D-A <=3.5 Å and D-H-A >=150 degrees, explicit donor hydrogens; name as occupancy of hydrogen-bond geometry.','residence':'Occupancy is not a continuous lifetime or dissociation rate.','twenty_ns_windows':'Descriptive only; mixed correlation signs across replicas do not support Hchi -> RMSF direction.'},
    'context_control_contract':'TCR coordinates/chains D/E are removed while all pMHC coordinates, sequence, chain order remain. This demonstrates the model output depends on binding-partner context; it is not a simulated biological apo state, causal residue perturbation, or experimental unbinding.',
    'primary_article_contract':{'title':'T-cell Receptor Specificity Maintained by Altered Thermodynamics','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC3696650/','experiments':['phage-display-derived receptor comparison','SPR affinity and kinetics','X-ray crystallography','peptide Ala substitution','temperature-dependent SPR thermodynamics'],'reported_MD':False,'KD_WT_uM':18,'KD_engineered_nM':0.6,'affinity_fold':30000,'engineered_koff_per_second':1.09e-4,'WT_koff_per_second':'>0.1; too fast to measure','total_MHC_contacts_WT_engineered':[56,97],'CDR2beta_contacts_WT_engineered':[6,38],'binding_entropy':'The article discusses total thermodynamic entropy and solvent reorganization; bound-state sidechain Hchi cannot substitute for that quantity.'},
    'argument_contract':{
        'supported':'The previously published affinity-matured receptor has a rearranged MHC contact network. Model outputs locate Q72 sidechain state restriction and lower motion of unchanged adjacent Gly54 in this patch. Existing trajectories independently support greater local contact occupancy and reduced Gly54 motion. This complements the primary structural/SPR result at residue resolution.',
        'not_supported':['Hchi directly causes Gly54 RMSF','Q72 alone explains the 30000-fold affinity change or a measured off-rate contribution','Lower bound-state Hchi is favorable total binding entropy or implies larger affinity','The original paper performed MD','All loops or all sidechains are predicted accurately or uniformly become rigid','All Q72 MD replica Hchi values are separated by state','Contact occupancy establishes longer residence times','Removing the TCR from the model is a physical apo experiment'],
        'figure_arrow_rule':'Contact or common-residue association may be shown with a line/brace; never draw a one-way causal arrow Hchi -> RMSF or Hchi -> affinity.'},
    'source_SHA256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
}
assert contract['all23_pMHC_interface_Hchi_original_predictions']['same_direction_n']==19
(OUT/'number_contract.json').write_text(json.dumps(contract,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'status':contract['status'],'key_model_values':contract['key_model_values'],'all23_original':contract['all23_pMHC_interface_Hchi_original_predictions'],'all23_context':contract['all23_pMHC_interface_Hchi_without_TCR'],'contract':str(OUT/'number_contract.json')},ensure_ascii=False,indent=2))
