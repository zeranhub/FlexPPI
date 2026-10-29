"""Freeze original all-residue prediction and descriptive environment comparisons."""
from pathlib import Path
import json, hashlib, itertools, shutil
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

O=Path(__file__).resolve().parent
S=O.parents[1]
OLD=S/'revision_v80_redesign'
U=OLD/'supplement_ubiquitin'/'source_data'
F=OLD/'figure5'/'source'/'original_fig7'
source=F/'experiment'/'all66_joint_predictions.csv'
r=pd.read_csv(source,float_precision='round_trip')
env=pd.read_csv(O/'all119_input_environment.csv')
r=r.merge(env,left_on=['chain','author_residue_id'],right_on=['chain','input_residue_id'],validate='one_to_one')
r['sample_id']='A003R';r['pdb_id']='2QHO'
r['residue_id']=r['author_residue_id']
r['MD_Hchi_nats']=r.frozen_MD_joint_Hchi_nats
r['Full_Hchi_nats']=r.predicted_Hchi_nats
r['absolute_error_nats']=abs(r.Full_Hchi_nats-r.MD_Hchi_nats)
r['model_protocol']='chemical_identity_corrected_v1 Full original fold0 epoch40 seed7'
r['reference_protocol']='mean of five consecutive 100ns windows, 1000 frames/window, first 5000 frames of existing single 500ns MD'
r['environment_protocol']='single exact starting complex, PyMOL3.0.0, solvent probe1.4A, dot_density4'
r.to_csv(O/'ubiquitin_all66_Hchi_environment_plot.csv',index=False,encoding='utf-8-sig')
r[~r.is_proline].to_csv(O/'ubiquitin_all63_nonPro_Hchi_sensitivity.csv',index=False,encoding='utf-8-sig')
gln=r[r.amino_acid=='Q'].copy()
gln.to_csv(O/'ubiquitin_all6_Gln_Hchi_SASA_plot.csv',index=False,encoding='utf-8-sig')
profile=pd.read_csv(U/'ubiquitin_fold0_all_residues.csv',float_precision='round_trip')
profile=profile[profile.valid_hchi & (profile.chain=='A')].copy()
assert len(profile)==66 and len(r)==66
np.testing.assert_allclose(profile.pred_hchi,r.Full_Hchi_nats,rtol=0,atol=6e-7)
np.testing.assert_allclose(profile.md_hchi,r.MD_Hchi_nats,rtol=0,atol=6e-7)
key=r[r.author_residue_id.isin([2,41,62])].copy()
key.to_csv(O/'key_Q2_Q41_Q62_exact_values.csv',index=False,encoding='utf-8-sig')

def rho(a,b):
    return None if len(a)<3 or len(set(a))<2 or len(set(b))<2 else float(spearmanr(a,b).statistic)
def metric(df):
    return {'n':len(df),'MAE_nats':float(df.absolute_error_nats.mean()),
        'Spearman_rho':rho(df.MD_Hchi_nats,df.Full_Hchi_nats)}
summary={'all66':metric(r),'all63_nonPro_sensitivity':metric(r[~r.is_proline]),
    'all6_Gln':metric(gln),'Gln_SASA_vs_MD_Hchi_rho':rho(gln.sidechain_SASA_A2,gln.MD_Hchi_nats),
    'Gln_SASA_vs_Full_Hchi_rho':rho(gln.sidechain_SASA_A2,gln.Full_Hchi_nats)}
center=r[~r.is_proline].copy()
for col in ['MD_Hchi_nats','Full_Hchi_nats','sidechain_SASA_A2','sc_residue_contacts_4p5A_intra_nonlocal']:
    center[col+'_within_AA_centered']=center[col]-center.groupby('amino_acid')[col].transform('mean')
summary['within_AA_centered_excluding_Pro']={
    'n_residues':len(center),'SASA_MD_rho':rho(center.sidechain_SASA_A2_within_AA_centered,center.MD_Hchi_nats_within_AA_centered),
    'SASA_Full_rho':rho(center.sidechain_SASA_A2_within_AA_centered,center.Full_Hchi_nats_within_AA_centered),
    'contact_MD_rho':rho(center.sc_residue_contacts_4p5A_intra_nonlocal_within_AA_centered,center.MD_Hchi_nats_within_AA_centered),
    'contact_Full_rho':rho(center.sc_residue_contacts_4p5A_intra_nonlocal_within_AA_centered,center.Full_Hchi_nats_within_AA_centered),
    'warning':'single protein descriptive; within-AA-centered ranks do not create independent biological samples; no inferential p value'}
pairs=[]
for aa,d in center.groupby('amino_acid'):
    for (_,a),(_,b) in itertools.combinations(d.iterrows(),2):
        lo,hi=(a,b) if a.sidechain_SASA_A2<=b.sidechain_SASA_A2 else (b,a)
        pairs.append({'amino_acid':aa,'low_SASA_residue':lo.residue,'high_SASA_residue':hi.residue,
          'low_SASA_A2':lo.sidechain_SASA_A2,'high_SASA_A2':hi.sidechain_SASA_A2,
          'SASA_difference_A2':hi.sidechain_SASA_A2-lo.sidechain_SASA_A2,
          'MD_Hchi_difference_high_minus_low':hi.MD_Hchi_nats-lo.MD_Hchi_nats,
          'Full_Hchi_difference_high_minus_low':hi.Full_Hchi_nats-lo.Full_Hchi_nats})
p=pd.DataFrame(pairs)
p['Full_agrees_MD_direction']=p.MD_Hchi_difference_high_minus_low*p.Full_Hchi_difference_high_minus_low>0
p['MD_high_SASA_higher_entropy']=p.MD_Hchi_difference_high_minus_low>0
p['Full_high_SASA_higher_entropy']=p.Full_Hchi_difference_high_minus_low>0
p.to_csv(O/'all_same_AA_unordered_pairs.csv',index=False,encoding='utf-8-sig')
groups=[]
for aa,d in center.groupby('amino_acid'):
    buried=d[d.sidechain_SASA_A2<=5]
    exposed=d[d.sidechain_SASA_A2>=40]
    if len(buried) and len(exposed):
        groups.append({'AA':aa,'n_buried_SASA_le5':len(buried),'n_exposed_SASA_ge40':len(exposed),
          'buried_residues':';'.join(buried.residue),'exposed_residues':';'.join(exposed.residue),
          'buried_MD_mean_Hchi':float(buried.MD_Hchi_nats.mean()),'exposed_MD_mean_Hchi':float(exposed.MD_Hchi_nats.mean()),
          'MD_exposed_minus_buried':float(exposed.MD_Hchi_nats.mean()-buried.MD_Hchi_nats.mean()),
          'buried_Full_mean_Hchi':float(buried.Full_Hchi_nats.mean()),'exposed_Full_mean_Hchi':float(exposed.Full_Hchi_nats.mean()),
          'Full_exposed_minus_buried':float(exposed.Full_Hchi_nats.mean()-buried.Full_Hchi_nats.mean())})
pd.DataFrame(groups).to_csv(O/'all_same_AA_buried_exposed_groups_descriptive.csv',index=False,encoding='utf-8-sig')
summary['absolute_SASA_groups_descriptive']={'definition':'buried sidechain SASA<=5A2 versus exposed>=40A2, same AA only, excluding Pro; post-hoc illustrative thresholds, not a universal burial criterion',
  'n_AA_types_with_both_groups':len(groups),'groups':groups}
summary['same_AA_pairs_descriptive']={}
for cutoff in [0,20,40]:
    d=p[p.SASA_difference_A2>cutoff]
    summary['same_AA_pairs_descriptive'][str(cutoff)]={
      'n_pairs':len(d),'n_amino_acid_types':int(d.amino_acid.nunique()),
      'n_MD_high_SASA_higher_entropy':int(d.MD_high_SASA_higher_entropy.sum()),
      'n_Full_high_SASA_higher_entropy':int(d.Full_high_SASA_higher_entropy.sum()),
      'n_Full_agrees_MD_direction':int(d.Full_agrees_MD_direction.sum()),
      'warning':'Pairs share residues and are not independent; this is a complete single-case descriptive diagnostic, not sample size for a statistical test'}
cor=[]
for aa,d in center.groupby('amino_acid'):
    cor.append({'AA':aa,'n':len(d),'SASA_MD_rho':rho(d.sidechain_SASA_A2,d.MD_Hchi_nats),
        'SASA_Full_rho':rho(d.sidechain_SASA_A2,d.Full_Hchi_nats),'Full_MD_rho':rho(d.MD_Hchi_nats,d.Full_Hchi_nats)})
pd.DataFrame(cor).to_csv(O/'per_AA_environment_correlations_descriptive.csv',index=False,encoding='utf-8-sig')
for fn in ['ubiquitin_Q41_Q62_all27_joint_frequencies.csv','ubiquitin_Q41_Q62_validation_summary.csv','representative_frames.csv','provenance.json']:
    shutil.copy2(U/fn,O/fn)
shutil.copy2(Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\mechanistic_section33_20261024\source_data\A003R_Q41_Q62_label_windows.csv'),O/'Q41_Q62_five_windows.csv')
summary['Q41_Q62_differences']={
 'MD_Q62_minus_Q41':float(gln.loc[gln.residue=='Q62','MD_Hchi_nats'].iloc[0]-gln.loc[gln.residue=='Q41','MD_Hchi_nats'].iloc[0]),
 'Full_Q62_minus_Q41':float(gln.loc[gln.residue=='Q62','Full_Hchi_nats'].iloc[0]-gln.loc[gln.residue=='Q41','Full_Hchi_nats'].iloc[0])}
summary['source_manifest']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,U/'ubiquitin_fold0_all_residues.csv',O/'all119_input_environment.csv']}
raw=Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\experiments\corrected_pipeline\neural_exports\evaluation\full\A003R.npz')
labels=Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\experiments\corrected_pipeline\training_payload\dataset\labels\A003R.npz')
z=np.load(raw,allow_pickle=False);lb=np.load(labels,allow_pickle=False)
idx=np.flatnonzero((z['chain']=='A') & z['valid'][:,1])
assert len(idx)==66
np.testing.assert_array_equal(z['author_residue_id'][idx].astype(int),r.author_residue_id)
np.testing.assert_array_equal(z['amino_acid'][idx],r.amino_acid)
np.testing.assert_allclose(z['mean'][idx,1],r.Full_Hchi_nats,rtol=0,atol=1e-14)
np.testing.assert_allclose(z['target'][idx,1],r.MD_Hchi_nats,rtol=0,atol=1e-14)
np.testing.assert_array_equal(z['target'],lb['values'])
np.testing.assert_array_equal(z['valid'],lb['mask'])
summary['direct_raw_cache_verification']={'all66_exact_identity_and_values':True,
 'cache_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),
 'label_sha256':hashlib.sha256(labels.read_bytes()).hexdigest(),
 'label_metric_names':lb['metric_names'].tolist(),
 'Pro_is_chemically_defined_and_valid':True,
 'original_valid_mask_changed':False,
 'Pro_zero_interpretation':'true zero of the frozen 120-degree joint-bin reference; not a forced missing-angle zero and not absence of motion'}
shutil.copy2(raw,O/'original_frozen_Full_A003R.npz')
shutil.copy2(labels,O/'original_frozen_labels_A003R.npz')
padding=profile[['sample_id','pdb_id','fold','chain','residue_id','amino_acid']].copy()
all_a=pd.read_csv(U/'ubiquitin_fold0_all_residues.csv',float_precision='round_trip')
all_a=all_a[all_a.chain=='A'].copy()
padding=all_a[['sample_id','pdb_id','fold','chain','residue_id','amino_acid','valid_hchi']].copy()
all_idx=np.flatnonzero(z['chain']=='A')
np.testing.assert_array_equal(all_a.residue_id,z['author_residue_id'][all_idx].astype(int))
padding['MD_Hchi_nats']=np.where(z['valid'][all_idx,1],z['target'][all_idx,1].astype(float),np.nan)
padding['Full_Hchi_nats']=np.where(z['valid'][all_idx,1],z['mean'][all_idx,1].astype(float),np.nan)
padding.to_csv(O/'ubiquitin_72_position_profile_with_original_mask.csv',index=False,encoding='utf-8-sig')
(O/'ubiquitin_case_audit.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
print(gln[['residue','MD_Hchi_nats','Full_Hchi_nats','absolute_error_nats','sidechain_SASA_A2','sc_residue_contacts_4p5A_intra_nonlocal','sc_residue_contacts_5p0A_partner']].to_string(index=False))
