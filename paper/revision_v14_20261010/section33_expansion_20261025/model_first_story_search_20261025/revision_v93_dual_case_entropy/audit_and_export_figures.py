"""Verify complete source-to-plot correspondence, actual intervals and asset outputs."""
from pathlib import Path
import json,hashlib,inspect,importlib.util,sys,subprocess,dataclasses,shutil
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from PIL import Image
ROOT=Path(__file__).resolve().parent;F=ROOT/'figure3';P=F/'source_data';BASE=ROOT.parents[2]
O=ROOT/'audit';O.mkdir(exist_ok=True)
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(path,data):path.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
plot=pd.read_csv(P/'all_plotted_values.csv');cam=pd.read_csv(P/'cam_all25.csv');ub=pd.read_csv(P/'ubiquitin_72_position_profile_with_original_mask.csv');gln=pd.read_csv(P/'ubiquitin_all6_Gln_Hchi_SASA_plot.csv')
assert plot.groupby('panel').size().to_dict()=={'b':100,'c':132,'e':12,'f':332,'g':332}
assert not plot.fillna('').duplicated(['panel','kind','identity','series','metric']).any()
errors=[]
for row in plot.itertuples():
    if row.panel=='b':
        col='n_CaM_contact_residues4_5' if row.kind.startswith('CaM_contact') else {'Full':'pred_rmsf','Intra':'intra_rmsf','MD':'target_rmsf'}[row.series]
        v=cam.loc[cam.CaMKI_native_position.eq(int(row.identity)),col].item()
    elif row.panel=='c':
        v=ub.loc[ub.residue_id.eq(int(row.identity)),row.series+'_Hchi_nats'].item()
    elif row.panel=='e':
        q=gln[gln.residue.eq(row.identity)].iloc[0];v=q[row.series+'_Hchi_nats'];errors.append(abs(q.sidechain_SASA_A2-row.SASA_A2))
    elif row.panel=='f':
        z=pd.read_csv(P/'retained_absolute_correlations_all498.csv')
        v=z.loc[z.sample_id.eq(row.identity)&z.metric_key.eq(row.metric)&z.predictor.eq(row.series),'rho_type_centered'].item()
    else:
        z=pd.read_csv(P/'retained_prediction_groups_all498.csv');metric=row.kind.removeprefix('MD_')
        v=z.loc[z.sample_id.eq(row.identity)&z.metric_key.eq(metric)&z.predictor.eq('full_hchi'),row.series].item()
    errors.append(abs(v-row.value))
assert max(errors)<1e-12
assert ub.valid_hchi.sum()==66 and ub.loc[ub.amino_acid.eq('P'),'valid_hchi'].all()
assert (ub.loc[ub.amino_acid.eq('P'),'MD_Hchi_nats']==0).all()
assert ub.loc[~ub.valid_hchi,['MD_Hchi_nats','Full_Hchi_nats']].isna().all().all()

# Independently regenerate all source-stratified paired interval indices and means.
corr=pd.read_csv(P/'retained_absolute_correlations_all498.csv');groups=pd.read_csv(P/'retained_prediction_groups_all498.csv');stats=pd.read_csv(P/'retained_panel_statistics.csv')
meta=corr[['sample_id','dataset']].drop_duplicates().sort_values('sample_id').reset_index(drop=True)
assert meta.dataset.value_counts().to_dict()=={'atlas':30,'itrace':30,'dynarepo':23}
rng=np.random.default_rng(20261026);parts=[]
for source in ['atlas','dynarepo','itrace']:
    loc=np.flatnonzero(meta.dataset.to_numpy()==source);parts.append(loc[rng.integers(0,len(loc),size=(10000,len(loc)))])
indices=np.concatenate(parts,axis=1)
stored=np.load(ROOT/'entropy_audit/paired_source_stratified_bootstrap_indices.npz')
assert np.array_equal(indices,stored['record_indices']) and np.array_equal(meta.sample_id.to_numpy(str),stored['record_ids'])
ci_errors=[]
for row in stats.itertuples():
    source=corr if row.panel_kind=='correlation' else groups
    z=source[source.metric_key.eq(row.metric_key)&source.predictor.eq(row.predictor)].set_index('sample_id').loc[meta.sample_id]
    values=z[row.value_key].to_numpy();means=values[indices].mean(1);lo,hi=np.quantile(means,[.025,.975])
    ci_errors.extend([abs(values.mean()-row.mean),abs(lo-row.bootstrap_95_low),abs(hi-row.bootstrap_95_high)])
assert max(ci_errors)<1e-12
q=json.loads((F/'figure_QA.json').read_text(encoding='utf-8'))
assert q['min_font']>=9 and q['print_font_at_164mm']>=8 and not q['outside_text']
assets={ext:sha(F/f'fig3_dual_case_entropy.{ext}') for ext in ['png','pdf','svg']};assert assets==q['sha256']
with Image.open(F/'fig3_dual_case_entropy.png') as im:assert im.size==(3242,3986)
geometry=json.loads((ROOT/'pymol/saved_session_geometry_audit.json').read_text(encoding='utf-8'))
assert all(x['maximum_coordinate_difference_A']==0 for x in geometry.values() if isinstance(x,dict))
assert geometry['two_local_views_share_rotation_and_optical_scale']
q.update(visual_review='Root and independent case reviewer opened final exact PNG; no substantive layout issues.',reviewed_png_sha256=assets['png'],all908_values_reverified=True,all_24_summary_intervals_reverified=True)
write(F/'figure_QA.json',q)

fonts=subprocess.run([r'C:\texlive\2026\bin\windows\pdffonts.exe',str(F/'fig3_dual_case_entropy.pdf')],capture_output=True,text=True,check=True).stdout
(O/'figure3_pdffonts.txt').write_text(fonts,encoding='utf-8');assert 'ArialMT' in fonts and 'Arial-BoldMT' in fonts
assert all('yes yes yes' in ' '.join(line.split()) for line in fonts.splitlines()[2:] if line.strip())
validator=Path(r'C:\Users\DELL\.codex\skills\academic-figure-skill\scripts\qa_validator.py')
spec=importlib.util.spec_from_file_location('style_validator',validator);module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
findings=[]
for name,fn in inspect.getmembers(module,inspect.isfunction):
    if name.startswith('check_'):
        result=fn((F/'compose_fig3_v93.py').read_text(encoding='utf-8'));findings.extend(result if isinstance(result,list) else [result])
write(O/'automated_style_QA_raw.json',[dataclasses.asdict(x) for x in findings])

metrics=dict(CaM25_Full_RMSF_MAE_A=float(abs(cam.pred_rmsf-cam.target_rmsf).mean()),CaM25_Intra_RMSF_MAE_A=float(abs(cam.intra_rmsf-cam.target_rmsf).mean()),Ub66_Hchi_MAE_nats=float(abs(ub.loc[ub.valid_hchi,'Full_Hchi_nats']-ub.loc[ub.valid_hchi,'MD_Hchi_nats']).mean()),Gln6_Hchi_MAE_nats=float(abs(gln.Full_Hchi_nats-gln.MD_Hchi_nats).mean()),Gln6_rho=float(spearmanr(gln.Full_Hchi_nats,gln.MD_Hchi_nats).statistic))
assert abs(metrics['CaM25_Full_RMSF_MAE_A']-.1148592012)<1e-9
assert abs(metrics['Ub66_Hchi_MAE_nats']-.079531)<1e-6
result=dict(status='passed',all908_values_and_SASA_checked=True,source_value_max_abs_difference=max(errors),bootstrap_summary_rows=len(stats),CI_max_abs_difference=max(ci_errors),record_stratified_not_family_intervals=True,metrics=metrics,native_session_geometry_audit='passed',final_main_asset_sha256=assets,actual_font_embedding='passed',Word_managed_by_separate_build_audit=True,style_raw_counts={c:sum(x.category==c for x in findings) for c in ['PASS','WARN','FAIL']})
write(O/'final_figure_source_audit.json',result)
for ext in ['png','pdf','svg']:
    shutil.copy2(F/f'fig3_dual_case_entropy.{ext}',BASE/f'figure3_v93_20261026_双案例与构象占据.{ext}')
    shutil.copy2(ROOT/f'entropy_audit/figS4_joint_readouts_angular_occupancy.{ext}',BASE/f'figureS4_v93_20261026_构象占据补充核查.{ext}')
print(json.dumps(result,ensure_ascii=False))
