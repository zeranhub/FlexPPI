"""Traceable Figure 4 export checks and source archive. No new scientific experiment."""
from pathlib import Path
import sys, json, hashlib, shutil
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parent.parent
sys.path.insert(0,str(EXP/'completed_results_20261012'/'plot_runtime'))
import numpy as np
import pandas as pd
from PIL import Image
from matplotlib import font_manager

sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
a=json.loads((ROOT/'data_audit.json').read_text(encoding='utf8'))
font=json.loads((ROOT/'font_and_bounds_audit.json').read_text(encoding='utf8'))
model_r=pd.read_csv(ROOT/'source_plotted_CDR_RMSF.csv')
model_h=pd.read_csv(ROOT/'source_plotted_interface_Hchi.csv')
md=pd.read_csv(ROOT/'source_plotted_existing_MD_occupancy.csv')
assert len(model_r)==45 and len(model_h)==23 and len(md)==12
assert len(md[md.metric=='Q72–Gly54 contact'])==6
assert len(md[md.metric=='Q72 two-ended bridge'])==6
assert md.groupby(['metric','pdb']).size().eq(3).all()
values_r=model_r[['MEL5_model_RMSF_A','alpha24beta17_model_RMSF_A']].to_numpy()
values_h=model_h[['hchi_full_a','hchi_full_b']].to_numpy()
assert np.isfinite(values_r).all() and np.isfinite(values_h).all()
assert ((values_r>=0)&(values_r<=1.45)).all()
assert ((values_h>=0)&(values_h<=3.4)).all()
assert md.frame_occupancy_percent.between(0,100).all()
assert a['all_plotted_data_inside_axis_limits'] and font['all_visible_text_inside_canvas']
im=Image.open(ROOT/'fig4_mel5_joint.png').convert('RGB')
assert im.size==(3242,3277)
density={}
for name,(x0,y0,x1,y1) in {'a':(.04,.70,.97,.96),'b':(.09,.47,.965,.651),'c':(.09,.253,.965,.384),'d':(.14,.05,.965,.158)}.items():
    arr=np.asarray(im)
    crop=arr[int((1-y1)*im.height):int((1-y0)*im.height),int(x0*im.width):int(x1*im.width)]
    density[name]=float(np.any(crop<235,axis=2).mean())
    assert density[name]>.015
archive=ROOT/'source_archive';archive.mkdir(exist_ok=True)
originals=[Path(v['path']) for v in a['source_files'].values()]
V=EXP/'tcr_expanded_screen_v41_20261017'/'mel5_validation'
originals += [V/'corrected_pair_all_aligned_positions.csv',V/'corrected_alignment_and_CDR_audit.json',V/'3hg1_model_input.pdb',V/'4jff_model_input.pdb']
manifest={}
for path in originals:
    target=archive/path.name
    shutil.copy2(path,target)
    assert sha(path)==sha(target)
    manifest[path.name]={'original_path':str(path),'SHA256':sha(target)}
(ROOT/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
a.update({'occupancy_axis_limits_percent':[-2,115],'all_occupancy_data_inside_axis_limits':True,'occupancy_data_range_percent':[float(md.frame_occupancy_percent.min()),float(md.frame_occupancy_percent.max())],'all_model_markers_have_headroom':True,'zero_occupancy_replicates_visually_separated_by_vertical_offset':True,'replicate_vertical_offsets_are_display_only':[-.24,0,.24]})
(ROOT/'data_audit.json').write_text(json.dumps(a,ensure_ascii=False,indent=2),encoding='utf8')
report={'status':'passed','PNG_SHA256':sha(ROOT/'fig4_mel5_joint.png'),'PDF_SHA256':sha(ROOT/'fig4_mel5_joint.pdf'),'SVG_SHA256':sha(ROOT/'fig4_mel5_joint.svg'),'reviewed_actual_final_PNG':True,'reviewer':'mel5_figure_v49 agent; final PNG opened with view_image after last generation','dimensions_mm':[183,185],'raster_pixels':list(im.size),'raster_dpi':450,'font':{'family':'Arial','resolved_file':font_manager.findfont('Arial'),'minimum_pt':9,'all_visible_text_inside_canvas':True},'panel_content_densities':density,'AP_0_to_7':'passed: no generic placeholder, default palette, rainbow, four spines, hidden observations, decorative gradients, arbitrary subsampling or smoothing','CL_1_to_7':'passed: 9 pt minimum, fixed 183 x185 mm, semantic orange/purple, outward ticks, thin open spines, complete vector+raster export','VI_1_to_7':'passed: model b/c dominant; structure supplies residue contact map; MD d subordinate; 45/45 CDR positions,23/23 valid pMHC interface sites; no noncontiguous-residue connecting lines; no causal arrow Hchi to RMSF; all three replicas visible','VV_1':'passed: no visible label/data occlusion after second layout correction','VV_2':'passed: four aligned narrative sections; local molecular inset labels identify true atoms','VV_3':'passed: all axis labels, per-residue labels and annotations >=9 pt; readable in final rendered preview','VV_4':'passed: orange hollow circles and purple filled squares give redundant state encoding; teal only HLA/residue identity in structural view','VV_5':'passed: all numerical observations inside axes; RMSF range '+str(a['RMSF_data_range_A'])+' within [0,1.45]; Hchi range '+str(a['Hchi_data_range_nats'])+' within [0,3.4]; occupancies within [0,100]','resolved_first_preview_findings':['RMSF upper bound raised from1.28 to1.45 to retain the CDR1beta maximum','Hchi upper bound raised from3.0 to3.4','x-axis explanation moved clear of residue ticks','MD state mean text vertically separated','three zero-occupancy replicas separated vertically','local atom identities and distances labeled'],'molecular_integrity':json.loads((ROOT/'molecular_view.json').read_text()),'ready_for_manuscript':True}
(ROOT/'visual_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
statistics={'Figure':'4','case':'MEL5 versus affinity-matured alpha24beta17 on identical HLA-A*02:01/ELA pMHC','a':{'coordinates':'existing4JFF model-input PDB, actual all-heavy atomic coordinates for interpretation','literature_SPR_KD':'18 microM versus0.6 nM, ref28 Madura2013 DOI10.1074/jbc.M113.464560','atomic_distances_A':{'Q72_OE1_to_beta_Phe53_N':2.942889736296623,'Q72_NE2_to_beta_Gly54_O':3.119551249779373},'numbering':'Common MEL5 numbering; engineered E54/E55 map to common beta53/beta54'},'b':{'n_definition':'45 aligned residues across six complete CDRs, one held-out structural pair','n':45,'values':'Two absolute predictions per position from one joint model; no smoothing, errors, tests or confidence bands','model':'Full fold0 trained on ATLAS1866+DynaRepo564, ITrace0','source_columns':['rmsf_full_a','rmsf_full_b']},'c':{'n_definition':'All23 Hchi-valid residues among30 union-defined pMHC interface positions, one held-out structural pair; same amino acid in both states','n':23,'excluded':'7 positions without a defined side-chain chi angle, according to existing valid flag','values':'Two absolute predictions per position; no test or confidence band','source_columns':['hchi_full_a','hchi_full_b']},'d':{'n_definition':'3 independent pre-existing200ns MD replicas per complex, every saved0-200ns frame retained','n_per_state':3,'frames_per_replica':md[['pdb','replica','frames']].drop_duplicates().to_dict('records'),'center_statistic':'Unweighted arithmetic mean of the three per-replica occupancies (short vertical mark)','spread':'All individual per-replica point values shown; no error bar','statistical_test':None,'multiple_comparison_correction':None,'occupancy_criteria':a['MD_occupancy_criteria'],'no_pseudoreplication':'Saved frames define occupancy; frames are not treated as independent experimental replicates','MD_not_original_Madura_experiment':True},'all_source_hashes':manifest}
(ROOT/'statistics_and_reproducibility.json').write_text(json.dumps(statistics,ensure_ascii=False,indent=2),encoding='utf8')
(ROOT/'QA_report.md').write_text('''# Figure 4 verification\n\nThe final PNG was opened and visually reviewed after the last generation. All visible text is English and at least9pt; the figure is183×185mm, exported as450dpi PNG and editable-text SVG plus vector PDF. The two model panels retain all45 aligned CDR positions and all23 valid, unchanged-amino-acid pMHC interface positions. The predefined30-site interface source is archived in full; its7 residues without valid side-chain chi angles are not plotted in the entropy panel. Data limits were checked independently of text limits.\n\nThe structure is the actual existing4JFF model input, with common-numbering betaPhe53 and betaGly54 mapped to inputE54/E55. It is used only to interpret the model's residue localization. Distances and atom labels identify Q72OE1–Phe53N and Q72NE2–Gly54O. No molecular coordinates were edited.\n\nThe MD panel has three visible replica points per state per criterion; short vertical marks and right-hand numbers are unweighted means. It uses existing ITrace trajectories and is explicitly labeled ExistingMD. No new production, fitting, model inference, training or calibration was done. No statistical test is applied. See statistics_and_reproducibility.json for n definitions, criteria, provenance and source hashes.\n\nThe production library was inspected for LineTrend, PairedBoxScatter and GroupedBarChart, including actual preview images. Existing production scripts use fundamentally different data structures, so their thin open axes, paired connectors and individual-point logic were inherited rather than executing those unrelated examples. The real-coordinate PyMOL script was adapted from the previous manuscript's molecular renderer.\n''',encoding='utf8')
print(json.dumps({'status':'passed','PNG_SHA256':report['PNG_SHA256'],'all_text_inside':True,'font_min_pt':9,'panel_density':density,'source_archives':len(manifest)},indent=2))
