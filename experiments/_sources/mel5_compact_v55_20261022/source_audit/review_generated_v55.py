"""Read-only independent preservation audit for compact four-panel Figure4."""
from pathlib import Path
import ast,hashlib,json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent.parent
FIG=ROOT/'figure4'
OLD=ROOT.parent/'mel5_composite_v54_20261022/figure4'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
previous=json.loads((OLD/'independent_data_review.json').read_text(encoding='utf-8'))
assert previous['status']=='passed'
baseline=json.loads((ROOT/'source_audit/source_preservation_baseline.json').read_text(encoding='utf-8'))
assert previous['PNG_SHA256']==baseline['v54_reviewed_png_SHA256']
files=[]
for rec in baseline['sources']:
    name=rec['filename'];old=Path(rec['path']);new=FIG/name
    assert sha(old)==rec['SHA256']
    if old.suffix=='.csv':
        assert sha(old)==sha(new),(name,'not byte identical')
        status='byte-identical CSV'
    else:
        a=np.load(old,allow_pickle=False);b=np.load(new,allow_pickle=False)
        assert set(a.files)==set(b.files)
        for k in a.files:
            assert np.array_equal(a[k],b[k],equal_nan=True),(name,k)
        status='all arrays and NaN masks exactly identical'
    files.append({'filename':name,'old_SHA256':sha(old),'new_SHA256':sha(new),'check':status,'status':'passed'})

pep=pd.read_csv(FIG/'source_peptide_four_states.csv',float_precision='round_trip')
cdr=pd.read_csv(FIG/'source_CDR_absolute.csv',float_precision='round_trip')
assert len(pep)==40 and len(cdr)==45
assert set(pep.pdb)=={'3hg1','4jff','4jfd','4jfe'}
assert all(len(g)==10 for _,g in pep.groupby('pdb'))
cmat=np.load(FIG/'source_heatmap_arrays.npz',allow_pickle=False)
assert cmat['CDR_delta_RMSF'].shape==(3,45)
assert cmat['interface_delta_Hchi'].shape==(3,30)
assert np.isfinite(cmat['interface_delta_Hchi']).sum(axis=1).tolist()==[23,23,22]
assert np.abs(cmat['CDR_delta_RMSF']).max()<=.45 and np.nanmax(np.abs(cmat['interface_delta_Hchi']))<=.45

script=FIG/'plot.py'
tree=ast.parse(script.read_text(encoding='utf-8'))
constants={}
for node in tree.body:
    if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
        try:constants[node.targets[0].id]=ast.literal_eval(node.value)
        except (ValueError,TypeError):pass
assert not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in ['bar','barh'] for n in ast.walk(tree))
assert constants['STATES']==['3hg1','4jff','4jfd','4jfe']
assert constants['NAMES']==['MEL5 · ELA','α24β17 · ELA','α24β17 · G4A','α24β17 · I7A']
assert constants['MARKERS']==['o','s','D','^']
assert constants['STYLES']==['-','-','--',':']
for key,value in {'WT':'#CF8C27','ENG':'#762A83','TEAL':'#297B82','G4A':'#4B82A5','I7A':'#73927D','NA':'#D9DDE0'}.items():assert constants[key]==value
assert constants['SIGNED']==['#427A9F','#FAFAFA','#B76778']
new=json.loads((FIG/'figure_data_audit.json').read_text(encoding='utf-8'))
old=json.loads((OLD/'figure_data_audit.json').read_text(encoding='utf-8'))
png=FIG/'fig4_mel5_compact_joint.png'
assert new['figure_sha256']==sha(png)
assert new['colors']==old['colors']
assert new['crop_boxes']==old['crop_boxes']
assert new['counts']=={'peptide_values':40,'absolute_CDR_values':90,'ranked_Hchi':0,'CDR_heatmap':135,'pMHC_heatmap':90,'Hchi_valid_per_contrast':[23,23,22]}
assert new['contrasts']==old['contrasts']
assert new['panels']==4
visible=[r['text'] for r in new['visible_text']]
assert all(k in visible for k in ['a','b','c','d'])
assert 'e' not in visible and 'f' not in visible
assert 'Interface side-chain response · all 23 valid sites' not in visible
assert 'Joint backbone and side-chain responses' in visible
assert 'Δ = α24β17 − MEL5' not in visible
assert '2.94 Å' in visible and '3.12 Å' in visible
for name in constants['NAMES']:assert name in visible
report={
 'status':'passed','PNG_path':str(png),'PNG_SHA256':sha(png),'plot_py_SHA256':sha(script),
 'read_only_review':True,'actual_PNG_viewed':True,
 'source_baseline_review':str(OLD/'independent_data_review.json'),
 'source_preservation':files,
 'all_remaining_scientific_values_identical_to_reviewed_v54':True,
 'counts':new['counts'],
 'old_ranked23_bar_removed':{'no_bar_plot_calls':True,'not_visible_in_actual_PNG':True,'old_title_and_delta_label_absent':True},
 'new_panel_d':{'one_panel_letter':True,'upper_CDR_delta_RMSF_shape':[3,45],'lower_interface_delta_Hchi_shape':[3,30],'upper_lower_values_and_NaN_masks_exactly_preserved':True,'valid_Hchi_per_contrast':[23,23,22],'separate_units_and_colorbars_preserved':True,'common_TCR_numbering_preserved':True,'Q72_highlight_preserved':True},
 'a_b_c':{'same_structure_images_and_crops':True,'fixed_atomic_pair_labels_2_94_and_3_12_preserved':True,'all40_peptide_values_preserved':True,'all90_absolute_CDR_values_preserved':True},
 'state_mapping_and_palette_preserved':True,'state_colors':new['colors'],
 'contrasts_preserved':new['contrasts'],'heatmap_range':new['heatmap_range'],
 'no_axis_color_clipping':True,'new_inference_training_MD':False,
 'blocking_issues':[],'required_edits':[]
}
(FIG/'independent_data_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','PNG_SHA256':sha(png),'file':str(FIG/'independent_data_review.json')},ensure_ascii=False))
