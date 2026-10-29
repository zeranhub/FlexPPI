from pathlib import Path
import pandas as pd, numpy as np, json, hashlib
ROOT=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025')
SRC=ROOT/'revision_v88_two_panel_story'
OUT=ROOT/'revision_v89_mechanism_story/data_audit_only'
OUT.mkdir(parents=True,exist_ok=True)
source=SRC/'figure3/source_data/all25_dynamics_and_static_lobe_contacts.csv'
d=pd.read_csv(source)
rw=pd.read_csv(SRC/'data_audit/MD_rmsf_five_windows_25_residues.csv')
hw=pd.read_csv(SRC/'data_audit/MD_hchi_five_windows_21_sidechains.csv')
assert len(d)==25 and d.valid_rmsf.sum()==25 and d.valid_chi_entropy.sum()==21
assert ''.join(d.amino_acid)=='IKKNFAKSKWKQAFNATAVVRHMRK'
assert d.CaMKI_native_position.tolist()==list(range(294,319))
assert set(d.fold)=={1} and set(d.variant)=={'full'}
assert np.allclose(d.target_rmsf,d.intra_target_rmsf)
assert np.allclose(d.target_chi_entropy,d.intra_target_chi_entropy,equal_nan=True)
assert set(d.loc[~d.valid_chi_entropy,'CaMKI_native_position'])=={299,306,309,311}
d['native_label']=d.amino_acid+d.CaMKI_native_position.astype(str)
# Export only effective sidechain predictions; four alanine model outputs are masked, not observations.
cols=['native_label','peptide_position','CaMKI_native_position','amino_acid','valid_rmsf','target_rmsf','pred_rmsf','intra_rmsf','valid_chi_entropy','target_chi_entropy','pred_chi_entropy','intra_chi_entropy','sidechain_shielding_fraction','minimum_heavy_atom_distance_A','CaM_heavy_contact_residues4_5','C_lobe_count_4p5A','N_lobe_count_4p5A']
q=d[cols].copy()
q.loc[~q.valid_chi_entropy,['target_chi_entropy','pred_chi_entropy','intra_chi_entropy']]=np.nan
q.to_csv(OUT/'all25_masked_model_MD_structure.csv',index=False)
q[q.valid_chi_entropy].to_csv(OUT/'all21_valid_Hchi.csv',index=False)
regions={
 'all25':(294,318),
 'native_N_flank_294_302':(294,302),
 'structural_anchor_span_303_316':(303,316),
 'native_tail_317_318':(317,318),
 'screenshot_N145_147_native294_296':(294,296),
 'screenshot_mid150_166_native299_315':(299,315),
 'screenshot_C167_169_native316_318':(316,318),
}
summary=[];win=[];mapping={}
for name,(lo,hi) in regions.items():
 g=d[d.CaMKI_native_position.between(lo,hi)]
 mapping[name]={'native_start':lo,'native_end':hi,'native_labels':g.native_label.tolist(),'n_RMSF':len(g),'n_effective_Hchi':int(g.valid_chi_entropy.sum()),'Ala_Hchi_excluded':g.loc[~g.valid_chi_entropy,'native_label'].tolist(),'selection_basis':'complete peptide' if name=='all25' else ('static W303/M316 anchor positions, independent of prediction' if name=='structural_anchor_span_303_316' else ('screenshot definition only; not a structural core boundary' if name.startswith('screenshot') else 'flanks defined relative to the W303/M316 static anchor span'))}
 for metric,cols0,mask in [('RMSF_A',['target_rmsf','pred_rmsf','intra_rmsf'],g.valid_rmsf),('Hchi_nats',['target_chi_entropy','pred_chi_entropy','intra_chi_entropy'],g.valid_chi_entropy)]:
  valid=g[mask]
  for model,col in zip(['MD','Full','Intra'],cols0):
   vals=valid[col]
   summary.append({'region':name,'metric':metric,'source':model,'n':len(vals),'mean':float(vals.mean()),'median':float(vals.median()),'minimum':float(vals.min()),'maximum':float(vals.max())})
  w,col=(rw,'rmsf') if metric=='RMSF_A' else (hw,'hchi')
  for block,sg in w[w.CaMKI_native_position.between(lo,hi)].groupby('block'):
   win.append({'region':name,'metric':metric,'block':int(block),'start_ns':(int(block)-1)*100,'end_ns':int(block)*100,'n':len(sg),'mean':float(sg[col].mean()),'median':float(sg[col].median())})
sumdf=pd.DataFrame(summary);windf=pd.DataFrame(win)
sumdf.to_csv(OUT/'region_means_medians_Full_Intra_MD.csv',index=False)
windf.to_csv(OUT/'MD_window_region_means_medians.csv',index=False)
# Matched-type F298/F307 and N297/N308; K296 vs central Lys K304 avoids cross-type comparison.
pairs=[(298,307),(296,304),(297,308),(317,314),(318,304)]
contrasts=[]
for left,right in pairs:
 aa=d.set_index('CaMKI_native_position').loc[left,'amino_acid'];bb=d.set_index('CaMKI_native_position').loc[right,'amino_acid']
 assert aa==bb
 for metric,src,col in [('RMSF_A',rw,'rmsf'),('Hchi_nats',hw,'hchi')]:
  piv=src[src.CaMKI_native_position.isin([left,right])].pivot(index='block',columns='CaMKI_native_position',values=col)
  vals=(piv[left]-piv[right])
  contrasts.append({'left':f'{aa}{left}','right':f'{bb}{right}','metric':metric,'window_differences':{int(k):float(v) for k,v in vals.items()},'n_positive_windows':int((vals>0).sum()),'n_negative_windows':int((vals<0).sum()),'minimum':float(vals.min()),'maximum':float(vals.max())})
window_directions=[]
for region in ['native_N_flank_294_302','native_tail_317_318','screenshot_N145_147_native294_296','screenshot_C167_169_native316_318']:
 for stat in ['mean','median']:
  piv=windf[windf.metric.eq('RMSF_A')].pivot(index='block',columns='region',values=stat)
  diff=piv[region]-piv.structural_anchor_span_303_316
  window_directions.append({'region':region,'versus':'structural_anchor_span_303_316','statistic':stat,'window_differences':{int(k):float(v) for k,v in diff.items()},'positive_windows':int((diff>0).sum()),'negative_windows':int((diff<0).sum())})
site=[]
for pos in [294,295,296,298,303,304,307,308,316,317,318]:
 r=d.set_index('CaMKI_native_position').loc[pos]
 rr=rw[rw.CaMKI_native_position.eq(pos)]
 hh=hw[hw.CaMKI_native_position.eq(pos)]
 site.append({'label':r.native_label,'rmsf_Full':float(r.pred_rmsf),'rmsf_Intra':float(r.intra_rmsf),'rmsf_MD':float(r.target_rmsf),'Hchi_Full':float(r.pred_chi_entropy),'Hchi_Intra':float(r.intra_chi_entropy),'Hchi_MD':float(r.target_chi_entropy),'MD_window_RMSF_min':float(rr.rmsf.min()),'MD_window_RMSF_max':float(rr.rmsf.max()),'MD_window_Hchi_min':float(hh.hchi.min()),'MD_window_Hchi_max':float(hh.hchi.max()),'shielding_fraction':float(r.sidechain_shielding_fraction),'closest_CaM_distance_A':float(r.minimum_heavy_atom_distance_A),'CaM_contacts_4p5A':(r.CaM_heavy_contact_residues4_5 if isinstance(r.CaM_heavy_contact_residues4_5,str) else '')})
claims=[
 {'claim':'W303–M316 anchor span has continuous low bound-state backbone fluctuation','judgement':'supported','basis':'Full 14/14 RMSF are 0.447–0.689 A; MD 14/14 are 0.453–0.624 A. These values are lower than native N-flank and two-residue tail means in Full and full-trajectory MD.','limitation':'Low RMSF measures positional fluctuations under the used alignment, not alpha-helix content or binding kinetic lifetime.'},
 {'claim':'K296/F298 are more mobile than the central anchor span','judgement':'qualified support','basis':'K296 Full/MD 1.010/1.216 A; F298 0.626/0.874 A; central mean 0.513/0.553 A. F298 has only modest Full RMSF elevation and lies within the whole central Full range.','limitation':'Do not label the whole K296/F298 pair as a sharply distinct flexible latch. Their involvement in autoinhibition requires literature or separate 1A06 structure, not their RMSF alone.'},
 {'claim':'The C terminus retains relative flexibility in the bound complex','judgement':'qualified support','basis':'R317/K318 Full RMSF 0.778/1.064 A, MD 0.768/1.771 A; the increase is strongest at K318. Tail-versus-core mean direction is positive in all five MD windows.','limitation':'Only two tail residues, substantial time-window variation, and no measured terminal–domain coupling. A functional flexible leash or domain engagement mechanism remains a hypothesis.'},
 {'claim':'F298 and F307 illustrate distinct backbone/sidechain readouts','judgement':'supported descriptive contrast','basis':'Same Phe type: F298/F307 Full RMSF 0.626/0.447 A, MD 0.874/0.511 A, whereas Hchi Full 1.033/0.996 nats, MD 1.406/1.418 nats. Backbone difference does not require greater sidechain entropy.','limitation':'Five-window backbone difference direction and Hchi direction are exported separately; no independent-repeat p-values are claimed.'},
 {'claim':'Partner context constrains bound peptide compared with a real free peptide','judgement':'not measured by Full/Intra','basis':'Both separately trained models used the same bound structure and bound MD target; Intra is the within-chain model. Full/Intra values cannot be substituted for two physical binding states.','limitation':'Use partner-aware versus within-chain prediction labels, not bound/free or before/after CaM.'},
 {'claim':'Residual helix in free peptide, conformational selection plus folding funnel, or a C-first N-second temporal pathway','judgement':'hypotheses only; unsupported by these readouts','basis':'No free peptide trajectory, secondary-structure time series, transitions, ensemble reweighting, or time-resolved domain engagement data are provided.','limitation':'Intra predictions of low central RMSF are not proof of residual helix. Static lobe contacts cannot establish order of engagement.'},
 {'claim':'CaM restrains the regulatory segment and reduces return to autoinhibitory contacts','judgement':'mechanistic hypothesis consistent with bound-state spatial fluctuation pattern and prior static/literature context','basis':'A low-fluctuation anchor span and actual lobe contacts provide a plausible retained-segment geometry.','limitation':'1MXE lacks catalytic domain, and no binding/unbinding or regulatory-segment reinsertion event is measured. Do not claim measured reduced return probability or activity.'},
 {'claim':'RMSF/Hchi establish binding free energy, rotational speed, residence time or activation rate','judgement':'not inferable','basis':'RMSF is positional amplitude and Hchi is fixed-bin joint torsional occupancy entropy, each lacking transition timing and energetic reference.','limitation':'MD window variability is not confidence interval from independent replicas.'},
]
audit={'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'sequence':'IKKNFAKSKWKQAFNATAVVRHMRK','numbering':'CaMKI native 294–318; cohort peptide 145–169; first residue I294, not K294','model_scope':'CaM family-held-out fold 1; Full and Intra independently trained; same bound input structure and same bound MD target. No new model inference was run.','valid_counts':{'RMSF':25,'Hchi':21,'Hchi_excluded_Ala':[299,306,309,311]},'regions':mapping,'region_statistics':summary,'MD_window_region_statistics':win,'MD_window_region_direction_checks':window_directions,'same_type_MD_window_contrasts':contrasts,'selected_sites':site,'claim_audit':claims,'window_note':'Five consecutive 100 ns windows from one 500 ns trajectory, not independent repeats. Full-trajectory RMSF target is exported as frozen_rmsf in window source and must not be replaced with arithmetic average of separately calculated window RMSF values.','sidechain_note':'Compare same amino-acid types when discussing environment. Four Ala model outputs are masked because no effective sidechain chi exists. Different residue types and chi counts confound absolute Hchi across W/M/K/F.'}
(OUT/'mechanism_data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
print(sumdf.to_string(index=False))
print('\nWINDOW DIRECTIONS\n',json.dumps(window_directions,ensure_ascii=False))
print('\nSAME TYPE\n',json.dumps(contrasts,ensure_ascii=False))