"""Read existing Hchi predictions and existing replica labels, without inference."""
from pathlib import Path
import csv
import hashlib
import json
from collections import defaultdict
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
EXP=ROOT.parent
OLD=EXP/'kras_case_v41_20261017'
PRED=EXP/'tcr_case_v40_20261017/family_predictions/full'
REPS=EXP/'itrace_case_search_only_20261017/replica_recheck'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def readatoms(p):
    out=[]
    for l in p.read_text().splitlines():
        if l.startswith('ATOM  ') and l[76:78].strip() not in ('H','D'):
            out.append({'chain':l[21],'resid':l[22:27].strip(),'atom':l[12:16].strip(),'aa':l[17:20].strip(),
             'xyz':np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])})
    return out

sources=[];files={};atoms={};replicas={}
for state,pdb in [('WT','7ow5'),('G12D','7ow6')]:
    p=PRED/f'itrace_{pdb}.npz';files[state]=np.load(p);sources.append({'path':str(p),'sha256':sha(p)})
    p=OLD/f'case_sources/{pdb}_model_input.pdb';atoms[state]=readatoms(p);sources.append({'path':str(p),'sha256':sha(p)})
    p=REPS/pdb/f'{pdb}_all_residue_per_replica_labels.csv';rr=list(csv.DictReader(p.open(encoding='utf8',newline='')))
    sources.append({'path':str(p),'sha256':sha(p)})
    replicas[state]=defaultdict(list)
    for row in rr:
        replicas[state][int(row['index'])].append(row)
    assert len(rr)==835*3

a=files['WT'];b=files['G12D']
assert a['chain'].tolist()==b['chain'].tolist() and a['author_residue_id'].tolist()==b['author_residue_id'].tolist()
assert len(a['mean'])==len(b['mean'])==835
geodist={};pepdist={}
for state in ['WT','G12D']:
    aa=atoms[state]
    d6=np.array([r['xyz'] for r in aa if r['chain']=='C' and r['resid']=='6'])
    peptide=np.array([r['xyz'] for r in aa if r['chain']=='C'])
    grouped=defaultdict(list)
    for r in aa: grouped[(r['chain'],r['resid'])].append(r['xyz'])
    geodist[state]={};pepdist[state]={}
    for key,xyz in grouped.items():
        xyz=np.array(xyz)
        geodist[state][key]=float(np.sqrt(((xyz[:,None,:]-d6[None,:,:])**2).sum(2)).min())
        pepdist[state][key]=float(np.sqrt(((xyz[:,None,:]-peptide[None,:,:])**2).sum(2)).min())

allrows=[]
for i in range(835):
    chain=str(a['chain'][i]);resid=str(a['author_residue_id'][i]);aa=str(a['amino_acid'][i]);aab=str(b['amino_acid'][i])
    key=(chain,resid)
    rg={}
    for state in ['WT','G12D']:
        reps=sorted(replicas[state][i],key=lambda r:int(r['replica']))
        assert len(reps)==3
        assert all(r['chain']==chain and r['author_residue_id']==resid for r in reps)
        valid=[float(r['Hchi_nats']) for r in reps if r['Hchi_nats']]
        rg[state]={'nchi':int(reps[0]['nchi']),'values_nats':valid,
            'mean_nats':float(np.mean(valid)) if valid else None,
            'variance_nats2':float(np.var(valid,ddof=1)) if valid else None,
            'n_replicas':len(valid)}
    common=aa==aab and bool(a['valid'][i,1]) and bool(b['valid'][i,1]) and rg['WT']['nchi']>0 and rg['G12D']['nchi']>0
    reasons=[]
    if chain=='A' and min(geodist['WT'][key],geodist['G12D'][key])<=6.0:reasons.append('HLA_within_6A_of_peptide_P6_in_either_state')
    if chain in ('D','E') and min(pepdist['WT'][key],pepdist['G12D'][key])<=4.5:reasons.append('TCR_heavy_atom_contact_within_4.5A_of_any_peptide_atom')
    if chain=='A' and resid in ('70','114','116'):reasons.append('preidentified_HLA_D6_network_residue')
    if chain=='D' and resid in ('50','70'):reasons.append('original_article_named_indirect_electrostatic_TCR_alpha_residue')
    if chain=='C':reasons.append('all_peptide_positions')
    row={'index':i,'chain':chain,'input_residue_id':resid,'WT_aa':aa,'G12D_aa':aab,
       'WT_Hchi_valid':bool(a['valid'][i,1]),'G12D_Hchi_valid':bool(b['valid'][i,1]),
       'common_same_type_Hchi_valid':common,'selection_reasons':reasons,
       'WT_predicted_Hchi_nats':float(a['mean'][i,1]) if a['valid'][i,1] else None,
       'G12D_predicted_Hchi_nats':float(b['mean'][i,1]) if b['valid'][i,1] else None,
       'prediction_delta_nats':float(b['mean'][i,1]-a['mean'][i,1]) if common else None,
       'WT_nchi':rg['WT']['nchi'],'G12D_nchi':rg['G12D']['nchi'],
       'WT_existing_MD_Hchi_nats':rg['WT']['mean_nats'],'G12D_existing_MD_Hchi_nats':rg['G12D']['mean_nats'],
       'existing_MD_delta_nats':rg['G12D']['mean_nats']-rg['WT']['mean_nats'] if common else None,
       'WT_replica_Hchi_nats':rg['WT']['values_nats'],'G12D_replica_Hchi_nats':rg['G12D']['values_nats'],
       'WT_min_heavy_distance_to_P6_A':geodist['WT'][key], 'G12D_min_heavy_distance_to_P6_A':geodist['G12D'][key],
       'WT_min_heavy_distance_to_peptide_A':pepdist['WT'][key],'G12D_min_heavy_distance_to_peptide_A':pepdist['G12D'][key]}
    if common:
        assert np.isfinite(row['prediction_delta_nats'])
        assert abs(float(a['target'][i,1])-rg['WT']['mean_nats'])<1.3e-7
        assert abs(float(b['target'][i,1])-rg['G12D']['mean_nats'])<1.3e-7
        row['prediction_MD_change_direction_agree']=row['prediction_delta_nats']*row['existing_MD_delta_nats']>0
        row['replica_ranges_nonoverlap']=(max(rg['WT']['values_nats'])<min(rg['G12D']['values_nats']) or max(rg['G12D']['values_nats'])<min(rg['WT']['values_nats']))
    allrows.append(row)

def summarize(rows):
    rows=[r for r in rows if r['common_same_type_Hchi_valid']]
    ds=np.array([r['prediction_delta_nats'] for r in rows])
    return {'n_common_valid':len(rows),'prediction_decreased':int((ds<0).sum()),'prediction_increased':int((ds>0).sum()),
      'mean_prediction_delta_nats':float(ds.mean()), 'mean_absolute_prediction_delta_nats':float(abs(ds).mean()),
      'min_prediction_delta_nats':float(ds.min()),'max_prediction_delta_nats':float(ds.max()),
      'existing_MD_direction_agreement_count':sum(r['prediction_MD_change_direction_agree'] for r in rows),
      'existing_MD_replica_range_nonoverlap_count':sum(r['replica_ranges_nonoverlap'] for r in rows)}

selected=[r for r in allrows if r['selection_reasons']]
network=[r for r in allrows if r['chain']=='A' and r['input_residue_id'] in ('70','114','116')]
out={'status':'passed','scope':'Read-only second-output analysis using existing Full0 NPZ and previously computed three-replica labels; no new trajectories, inference, training, figure or manuscript changes.',
 'units':{'Hchi':'nats, natural-log joint side-chain categorical state entropy','delta':'G12D minus WT; compare only unchanged same-amino-acid positions with valid chi labels in both states'},
 'selection_policy':{'HLA':'All HLA residues within 6.0 A heavy-atom distance of peptide P6 in either static input; Q70/R114/D116 always explicitly included.',
  'TCR':'All TCR residues within 4.5 A heavy-atom distance of any peptide atom in either input; additionally alpha R50/K70 named in original study are included.',
  'peptide':'All 10 positions shown, invalid and mutation-type-confounded comparisons explicitly null.',
  'selection_not_based_on_Hchi_prediction_or_MD_change':True,'all_835_aligned_positions_retained_in_audit':True},
 'source_files':sources,'counts':{'all_aligned':len(allrows),'same_type_common_valid':sum(r['common_same_type_Hchi_valid'] for r in allrows),
  'selected_total_including_invalid':len(selected),'selected_common_valid':sum(r['common_same_type_Hchi_valid'] for r in selected)},
 'summaries':{'entire_complex':summarize(allrows),'selected_network_and_contact_residues':summarize(selected),
   'HLA_near_P6':summarize([r for r in selected if r['chain']=='A']),
   'TCR_contact_and_literature_named':summarize([r for r in selected if r['chain'] in ('D','E')]),
   'all_peptide_common_Hchi':summarize([r for r in selected if r['chain']=='C']),
   'Q70_R114_D116_network':summarize(network)},
 'predefined_HLA_network_all_three_positions':network,
 'all_selected_positions':selected,'all_835_aligned_positions':allrows,
 'assessment':{'worth_adding_as_main_affinity_mechanism_panel':False,
   'reason':'Hchi adds location-specific side-chain information, but does not yet give a sufficiently clear additional affinity-mechanism conclusion to justify a main panel. Mixed side-chain responses or differences from RMSF are not failures: the outputs measure different properties and stronger affinity does not require all side chains to order.',
   'positive_information':'Q70/D116 lower Hchi is consistent with altered pocket-side-chain constraints; K10 shows a decrease. Different side-chain responses could complement backbone suppression as evidence of local dynamical redistribution.',
   'negative_information':'R114 does not decrease, whereas the existing MD labels decrease with non-overlapping replica ranges. TCR beta H101/N102 also increase in the prediction but decrease in existing MD. V8 slightly increases despite reduced backbone RMSF. Most central peptide positions are Gly/Ala without chi, and G6 vs D6 cannot be compared at fixed amino-acid type.',
   'independent_information_relative_to_RMSF':'Side-chain joint-state diversity is different from backbone positional fluctuations. Q70/D116 change Hchi while their backbone RMSF changes are small; V8 shows lower backbone RMSF and nearly unchanged/slightly higher Hchi. This is independent dynamical content, but it does not by itself make the affinity-mechanism claim stronger.',
   'R114_interpretation':'Nearly unchanged coarse joint-state entropy does not imply that the D6-R114 salt bridge is absent. Specific state identity or contact geometry may change without changing the entropy summary. The crystallographic contact is independently established by the original article.',
   'MD_role':'The existing MD labels are supplementary evidence with their own replica variability, not an absolute truth used alone to accept or reject the model story.',
   'limitations':['Joint Hchi does not identify a specific chi angle, salt-bridge occupancy, electrostatic energy or affinity.',
    'Gly/Ala model slots are masked placeholders; they must never be graphed as physically meaningful side-chain entropy.',
    'Bound-state Hchi differences do not measure the entropy change on binding; no unbound state is analyzed.',
    'Reported MD labels are existing independent ITrace analyses, with per-replica variability shown rather than treated as exact truth.']}}
(ROOT/'audit/Hchi_case_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
with (ROOT/'audit/Hchi_all_selected_positions.csv').open('w',encoding='utf8',newline='') as f:
    fields=['chain','input_residue_id','WT_aa','G12D_aa','WT_nchi','G12D_nchi','common_same_type_Hchi_valid',
     'WT_predicted_Hchi_nats','G12D_predicted_Hchi_nats','prediction_delta_nats','WT_existing_MD_Hchi_nats',
     'G12D_existing_MD_Hchi_nats','existing_MD_delta_nats','WT_min_heavy_distance_to_P6_A','G12D_min_heavy_distance_to_P6_A','selection_reasons']
    writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(selected)
print(json.dumps({'counts':out['counts'],'summaries':out['summaries'],
 'selected':[{'site':r['chain']+r['WT_aa']+r['input_residue_id'],'same_type_valid':r['common_same_type_Hchi_valid'],
 'WT':r['WT_predicted_Hchi_nats'],'G12D':r['G12D_predicted_Hchi_nats'],'delta':r['prediction_delta_nats'],'MD_delta':r['existing_MD_delta_nats']} for r in selected]},ensure_ascii=False,indent=2))
