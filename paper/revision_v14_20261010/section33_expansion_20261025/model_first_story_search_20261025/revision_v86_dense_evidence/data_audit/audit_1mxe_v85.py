from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import mdtraj as md
from scipy.spatial.distance import cdist
from Bio.SeqUtils import seq1
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025')
OUT=BASE/'revision_v85_1MXE/data_audit';OUT.mkdir(parents=True,exist_ok=True)
SOURCE=BASE/'case_selection_audit_20261026/CaMKI_294_318_family_fold1.csv'
STRUCTURE=Path(r'E:\Project\FlexPPI\experiments\_sources\baselines\cohort\structures\record_aba5e1f3db118b0aba8895de.pdb')
WINDOW=BASE/'partner_screen/verified_case_trajectories/A001K'
d=pd.read_csv(SOURCE)
assert len(d)==25 and set(d.fold)=={1} and set(d.variant)=={'full'}
assert ''.join(d.amino_acid)=='IKKNFAKSKWKQAFNATAVVRHMRK'
assert d.valid_rmsf.sum()==25 and d.valid_chi_entropy.sum()==21
assert np.allclose(d.target_rmsf,d.intra_target_rmsf)
assert np.allclose(d.target_chi_entropy,d.intra_target_chi_entropy,equal_nan=True)
d['anchor_span_303_316']=d.CaMKI_native_position.between(303,316)
d['region']=np.where(d.anchor_span_303_316,'anchor_span303_316','remaining11')
d.to_csv(OUT/'peptide_25_family_fold1.csv',index=False)
# Complete peptide and structural-first region summaries, never postselected by predictions.
summary=[]
regions={'all25':d,'anchor_span303_316':d[d.anchor_span_303_316],'remaining11':d[~d.anchor_span_303_316],'continuous299_317':d[d.CaMKI_native_position.between(299,317)]}
for name,g in regions.items():
 for metric in ['rmsf','chi_entropy']:
  q=g[g['valid_'+metric]]
  r=dict(region=name,metric=metric,n_residues=len(q),MD_mean=q['target_'+metric].mean(),Full_mean=q['pred_'+metric].mean(),Intra_mean=q['intra_'+metric].mean(),Full_MAE=(q['target_'+metric]-q['pred_'+metric]).abs().mean(),Intra_MAE=(q['target_'+metric]-q['intra_'+metric]).abs().mean())
  r['Full_minus_Intra']=r['Full_mean']-r['Intra_mean'];r['Full_relative_MAE_reduction_percent']=100*(1-r['Full_MAE']/r['Intra_MAE'])
  summary.append(r)
pd.DataFrame(summary).to_csv(OUT/'region_absolute_means_and_errors.csv',index=False)
# Original windows are reblocked summaries of a single 500 ns trajectory. model_rmsf in old source ignored.
rw=pd.read_csv(WINDOW/'all_residues_rmsf_per_window.csv');hw=pd.read_csv(WINDOW/'all_sidechains_per_window.csv')
mapcols=['residue_index','CaMKI_native_position','peptide_position','anchor_span_303_316','region']
rw=rw.merge(d[mapcols],on='residue_index',how='inner');hw=hw.merge(d[mapcols],on='residue_index',how='inner')
assert len(rw)==125 and len(hw)==105 and set(hw.block)=={1,2,3,4,5}
rw.drop(columns=['model_rmsf'],errors='ignore').to_csv(OUT/'MD_rmsf_five_windows_25_residues.csv',index=False)
hw.to_csv(OUT/'MD_hchi_five_windows_21_sidechains.csv',index=False)
win=[]
for block in range(1,6):
 for name,q in regions.items():
  pos=q.CaMKI_native_position
  for metric,w,col in [('rmsf',rw,'rmsf'),('chi_entropy',hw,'hchi')]:
   v=w[(w.block==block)&w.CaMKI_native_position.isin(pos)]
   win.append(dict(block=block,start_ns=(block-1)*100,end_ns=block*100,region=name,metric=metric,n_residues=len(v),MD_mean=v[col].mean()))
pd.DataFrame(win).to_csv(OUT/'MD_region_five_windows.csv',index=False)
# Use the exact cohort/model coordinates and heavy atoms. Matched deletion is static context removal, NOT free-peptide MD.
t=md.load_pdb(str(STRUCTURE))
heavy=t.topology.select('not element H');t=t.atom_slice(heavy)
assert t.n_residues==169
assert ''.join(seq1(r.name) for r in t.topology.residues if r.chain.index==1)=='IKKNFAKSKWKQAFNATAVVRHMRK'
cam=np.array([a.index for a in t.topology.atoms if a.residue.chain.index==0])
pep=np.array([a.index for a in t.topology.atoms if a.residue.chain.index==1])
xyz=t.xyz[0]*10
sasa=md.shrake_rupley(t,probe_radius=0.14,n_sphere_points=960,mode='atom')[0].astype(float)*100
pept=t.atom_slice(pep)
sasa_pep=md.shrake_rupley(pept,probe_radius=0.14,n_sphere_points=960,mode='atom')[0].astype(float)*100
pep_to_local={int(a):i for i,a in enumerate(pep)}
backbone={'N','CA','C','O','OXT'}
rows=[];pairs=[]
cam_res=list(t.topology.chain(0).residues)
for j,r in enumerate(t.topology.chain(1).residues):
 native=294+j; inds=np.array([a.index for a in r.atoms]); sides=np.array([a.index for a in r.atoms if a.name not in backbone]);bb=np.array([a.index for a in r.atoms if a.name in backbone]);atoms=list(t.topology.atoms)
 dst=cdist(xyz[inds],xyz[cam]);mini=np.unravel_index(np.argmin(dst),dst.shape)
 pa=atoms[int(inds[mini[0]])];ca=atoms[int(cam[mini[1]])]
 in_sasa=sasa[inds].sum();out_sasa=sum(sasa_pep[pep_to_local[int(i)]] for i in inds)
 in_side=sasa[sides].sum();out_side=sum(sasa_pep[pep_to_local[int(i)]] for i in sides)
 neighbors=[];neighbors45=[];neighbors5=[]
 for cr in cam_res:
  cinds=[a.index for a in cr.atoms];c_min=cdist(xyz[inds],xyz[cinds]).min()
  s_min=cdist(xyz[sides],xyz[cinds]).min() if len(sides) else float('nan')
  ca_idx=next(a.index for a in cr.atoms if a.name=='CA');pr_ca=next(a.index for a in r.atoms if a.name=='CA')
  c_alpha=np.linalg.norm(xyz[ca_idx]-xyz[pr_ca])
  if c_min<=4.0: neighbors.append(f'{seq1(cr.name)}{cr.resSeq+3}')
  if c_min<=4.5: neighbors45.append(f'{seq1(cr.name)}{cr.resSeq+3}')
  if c_min<=5.0: neighbors5.append(f'{seq1(cr.name)}{cr.resSeq+3}')
  pairs.append(dict(CaMKI_native_position=native,peptide_aa=seq1(r.name),CaM_local_id=cr.resSeq,CaM_native_id=cr.resSeq+3,CaM_aa=seq1(cr.name),CaM_lobe='N_lobe' if cr.resSeq+3<=75 else ('C_lobe' if cr.resSeq+3>=82 else 'linker'),minimum_heavy_atom_distance_A=c_min,minimum_peptide_sidechain_to_CaM_heavy_distance_A=s_min,CA_distance_A=c_alpha,heavy_atom_contact4=c_min<=4.0,heavy_atom_contact4_5=c_min<=4.5,heavy_atom_contact5=c_min<=5.0,CA8_contact=c_alpha<=8.0))
 row=dict(CaMKI_native_position=native,peptide_position=j+1,cohort_chain='B',cohort_residue_id=r.resSeq,amino_acid=seq1(r.name),anchor_span_303_316=303<=native<=316,all_atom_SASA_complex_A2=in_sasa,all_atom_SASA_without_CaM_same_coordinates_A2=out_sasa,all_atom_SASA_shielded_by_CaM_A2=out_sasa-in_sasa,sidechain_SASA_complex_A2=in_side,sidechain_SASA_without_CaM_same_coordinates_A2=out_side,sidechain_SASA_shielded_by_CaM_A2=out_side-in_side,sidechain_shielding_fraction=(out_side-in_side)/out_side if out_side>0 else np.nan,minimum_heavy_atom_distance_A=dst.min(),minimum_sidechain_heavy_atom_distance_A=cdist(xyz[sides],xyz[cam]).min() if len(sides) else np.nan,minimum_backbone_heavy_atom_distance_A=cdist(xyz[bb],xyz[cam]).min(),closest_peptide_atom=pa.name,closest_CaM_atom=ca.name,closest_CaM_residue=f'{seq1(ca.residue.name)}{ca.residue.resSeq+3}',CaM_heavy_contact_residues4=';'.join(neighbors),CaM_heavy_contact_residues4_5=';'.join(neighbors45),CaM_heavy_contact_residues5=';'.join(neighbors5),n_CaM_contact_residues4=len(neighbors),n_CaM_contact_residues4_5=len(neighbors45),n_CaM_contact_residues5=len(neighbors5),n_heavy_atom_pairs4=int((dst<=4.0).sum()))
 rows.append(row)
s=pd.DataFrame(rows);p=pd.DataFrame(pairs)
assert (s.all_atom_SASA_shielded_by_CaM_A2>=-1e-3).all()
s.to_csv(OUT/'static_SASA_and_heavy_contacts_all25.csv',index=False);p.to_csv(OUT/'all25_by_all144_CaM_residue_distances.csv',index=False)
joined=d.merge(s,on=['CaMKI_native_position','peptide_position','amino_acid','anchor_span_303_316'],validate='one_to_one')
joined.to_csv(OUT/'peptide25_model_MD_structure_all_values.csv',index=False)
# All recurring amino acid types; amino acid count and valid sidechains preserved.
joined['same_AA_count']=joined.groupby('amino_acid').amino_acid.transform('size')
same=joined[joined.same_AA_count>=2].sort_values(['amino_acid','CaMKI_native_position'])
same.to_csv(OUT/'all_repeated_amino_acid_residues.csv',index=False)
ag=[]
for aa,g in same.groupby('amino_acid',sort=True):
 q=g[g.valid_chi_entropy]
 ag.append(dict(amino_acid=aa,n_residues=len(g),n_effective_sidechains=len(q),positions=';'.join(map(str,g.CaMKI_native_position)),MD_Hchi_min=q.target_chi_entropy.min(),MD_Hchi_max=q.target_chi_entropy.max(),MD_Hchi_range=q.target_chi_entropy.max()-q.target_chi_entropy.min(),Full_Hchi_min=q.pred_chi_entropy.min(),Full_Hchi_max=q.pred_chi_entropy.max(),Intra_Hchi_min=q.intra_chi_entropy.min(),Intra_Hchi_max=q.intra_chi_entropy.max(),sidechain_SASA_complex_min_A2=g.sidechain_SASA_complex_A2.min(),sidechain_SASA_complex_max_A2=g.sidechain_SASA_complex_A2.max()))
pd.DataFrame(ag).to_csv(OUT/'all_repeated_amino_acid_groups.csv',index=False)
# Match N297/N308 for optional biological illustration, not separate selected evidence for full-peptide accuracy.
joined[joined.CaMKI_native_position.isin([296,297,298,303,308,316])].to_csv(OUT/'selected_labels_all_values.csv',index=False)
for n in [297,308]:
 print('\nN',n,joined.loc[joined.CaMKI_native_position==n,['sidechain_SASA_complex_A2','sidechain_SASA_without_CaM_same_coordinates_A2','sidechain_shielding_fraction','minimum_heavy_atom_distance_A','minimum_sidechain_heavy_atom_distance_A','CaM_heavy_contact_residues4','target_chi_entropy','pred_chi_entropy','intra_chi_entropy']].to_dict('records'))
# Mapping native CaM numbering follows mature CaM sequence ADQL... and the conserved M36/M51/M71/M72/M76/M109/M124/M144/M145 positions.
cammap=[]
for r in cam_res:
 cammap.append(dict(cohort_chain='A',cohort_residue_id=r.resSeq,CaM_mature_native_position=r.resSeq+3,amino_acid=seq1(r.name)))
pd.DataFrame(cammap).to_csv(OUT/'CaM_native_numbering_map.csv',index=False)
metadata=dict(status='passed',sample_id='A001K',pdb_id='1MXE',prediction_protocol='CaM family held out, family_fold1; no fold0 training-member prediction or old record-split fold4 prediction used',valid_RMSF_count=25,valid_Hchi_count=21,anchor_definition='W303 to M316 inclusive, 14 residues, chosen by structural double anchors before consulting prediction/MD errors',structural_SASA=dict(method='MDTraj Shrake-Rupley, heavy atoms, probe radius 1.4 Angstrom, 960 sphere points, identical coordinates with vs without CaM',static_only=True,not_unbound_simulation=True),contacts=dict(method='all heavy-atom minimum distances and 4.0 Angstrom contact cutoff; all 25 by 144 residue distances exported',CA8='record source CA<=8 Angstrom preserved only as approximate graph-neighborhood indicator; absence is not absence of atom-level contact'),trajectory=dict(duration_ns=500,n_contiguous_windows=5,window_duration_ns=100,independent_replicates=False,MD_targets_recomputed_at_original_bins=True,window_hchi_target_error=float(hw.groupby('residue_index').hchi.mean().reindex(d[d.valid_chi_entropy].residue_index).sub(d[d.valid_chi_entropy].set_index('residue_index').target_chi_entropy).abs().max()),window_rmsf_target_error=float(rw.groupby('residue_index').rmsf.mean().reindex(d.residue_index).sub(d.set_index('residue_index').target_rmsf).abs().max())),model_context_comparison='Full and Intra are alternative models on the bound complex target; Intra is not a true apo trajectory and differences do not measure a binding-induced change',supported=['Bound peptide 303-316 anchor-spanning region has lower mean backbone RMSF and sidechain Hchi than the remaining 11 residues','Full recovers whole-peptide bound-state RMSF and Hchi more accurately than Intra in this held-out family','Static structural environment distinguishes same-type N297/N308 and anchors, but all recurring AA groups are exported for audit'],unsupported=['No direct apo-to-holo dynamics transition, no activation free energy, no proof of N-lobe causal requirement from static SASA/model alone','Lower Hchi does not directly measure rotation speed or intrinsic continuous angular width','Cross-amino-acid Hchi magnitudes cannot be interpreted solely as environment effects because chi counts differ'],source_files=[dict(path=str(x),sha256=hashlib.sha256(x.read_bytes()).hexdigest()) for x in [SOURCE,STRUCTURE,WINDOW/'all_residues_rmsf_per_window.csv',WINDOW/'all_sidechains_per_window.csv',WINDOW/'audit.json']])
(OUT/'audit.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'README.md').write_text('''# 1MXE / A001K data audit\n\nAll model values are from the CaM-family-held-out fold 1. The 25-residue CaMKI peptide is 294–318; the structure-first anchor span is W303–M316, 14 residues, with the remaining 11 as the comparison region. All valid sidechains and all repeated amino-acid groups are retained. Full versus Intra compares two models against the same bound-state MD target; Intra is not an experimental unbound state.\n\nSASA compares identical model/cohort coordinates with and without CaM, using 1.4 Angstrom probe and 960 points. The difference measures static solvent shielding by CaM, not a free-peptide ensemble. Heavy-atom minimum distances and all contacting CaM residues are exported. CA8 absence must not be called no contact. CaM labels use mature native numbering (local +3), peptide labels native CaMKI numbering (local +149).\n\nMD five 100 ns windows come from a single 500 ns trajectory. They describe temporal consistency, not independent experimental or simulation repeats. Hchi is fixed-bin joint torsional entropy, not rotation speed. Cross-amino-acid magnitudes mix residue type and environment, so same-type examples and all groups are provided.\n\nRegenerate with the standalone audit_1mxe_v85.py script in this folder using the project numerical runtime.\n''',encoding='utf-8')
print(pd.DataFrame(summary).to_string(index=False));print('\nSAVED',OUT)
# Sensitivity of quoted static shielding to surface point density, same probe/coordinates.
sasa2880=md.shrake_rupley(t,probe_radius=0.14,n_sphere_points=2880,mode='atom')[0].astype(float)*100
pep2880=md.shrake_rupley(pept,probe_radius=0.14,n_sphere_points=2880,mode='atom')[0].astype(float)*100
sens=[]
for j,r in enumerate(t.topology.chain(1).residues):
 sides=[a.index for a in r.atoms if a.name not in backbone]
 complex2880=float(sasa2880[sides].sum());without2880=float(sum(pep2880[pep_to_local[int(i)]] for i in sides))
 old=s[s.CaMKI_native_position==294+j].iloc[0]
 sens.append(dict(CaMKI_native_position=294+j,sidechain_SASA_complex_960_A2=old.sidechain_SASA_complex_A2,sidechain_SASA_complex_2880_A2=complex2880,sidechain_SASA_shielded_960_A2=old.sidechain_SASA_shielded_by_CaM_A2,sidechain_SASA_shielded_2880_A2=without2880-complex2880,sidechain_shielding_fraction_960=old.sidechain_shielding_fraction,sidechain_shielding_fraction_2880=(without2880-complex2880)/without2880 if without2880 else np.nan))
pd.DataFrame(sens).to_csv(OUT/'static_SASA_point_density_sensitivity.csv',index=False)
# Five-window consistency is descriptive; no independent-replicate p values are computed.
wdf=pd.DataFrame(win);contr=[]
for metric in ['rmsf','chi_entropy']:
 q=wdf[wdf.metric==metric].pivot(index='block',columns='region',values='MD_mean')
 for block,r in q.iterrows():
  contr.append(dict(contrast='anchor_span303_316_minus_remaining11',metric=metric,block=int(block),MD_difference=float(r.anchor_span303_316-r.remaining11)))
q=hw[hw.CaMKI_native_position.isin([297,308])].pivot(index='block',columns='CaMKI_native_position',values='hchi')
for block,r in q.iterrows():contr.append(dict(contrast='N297_minus_N308',metric='chi_entropy',block=int(block),MD_difference=float(r[297]-r[308])))
pd.DataFrame(contr).to_csv(OUT/'MD_five_window_contrasts.csv',index=False)
qc=[]
for metric in ['rmsf','chi_entropy']:
 g=d[d['valid_'+metric]]
 fe=(g['pred_'+metric]-g['target_'+metric]).abs();ie=(g['intra_'+metric]-g['target_'+metric]).abs()
 qc.append(dict(metric=metric,n=len(g),Full_lower_point_error=int((fe<ie).sum()),Intra_lower_point_error=int((fe>ie).sum()),Full_Intra_equal_point_error=int((fe==ie).sum())))
pd.DataFrame(qc).to_csv(OUT/'per_residue_error_direction_counts.csv',index=False)
metadata['structure_preprocessing']=dict(n_protein_residues=169,calcium_atoms_in_cohort_structure=0,caution='Processed model coordinates contain protein only; model outputs and static SASA are not a measurement of calcium effects.')
metadata['static_SASA_point_density_sensitivity']='All 25 residues checked again at 2880 points; separate CSV provided. Main quoted values use 960 points.'
metadata['structural_SASA']['floating_point_note']='Atom areas accumulated in float64; deletion differences that are mathematically zero remain exactly zero for exposed residues.'
metadata['trajectory']['five_window_contrasts']=pd.DataFrame(contr).groupby(['contrast','metric']).MD_difference.agg(['min','max']).reset_index().to_dict('records')
(OUT/'audit.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
manifest=[dict(file=x.name,bytes=x.stat().st_size,sha256=hashlib.sha256(x.read_bytes()).hexdigest()) for x in sorted(OUT.iterdir()) if x.is_file() and x.name!='manifest.json']
(OUT/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')