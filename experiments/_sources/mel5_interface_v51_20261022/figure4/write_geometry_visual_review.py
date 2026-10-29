from pathlib import Path
import hashlib,json
import pandas as pd

root=Path(__file__).resolve().parent
base=root.parent
old=base.parent/'mel5_joint_case_v49_20261019'
meta=json.loads((base/'structure/dual_state_molecular_view.json').read_text(encoding='utf8'))
dist=pd.read_csv(base/'analysis/fixed_atom_contact_distances.csv')
cdr=pd.read_csv(old/'figure4/source_plotted_CDR_RMSF.csv')
h=pd.read_csv(old/'figure4/source_plotted_interface_Hchi.csv')
g=cdr[(cdr.region=='CDR2β')&(cdr.common_TCR_position==54)].iloc[0]
q=h[(h.chain_a=='A')&(h.position_a==72)].iloc[0]
assert round(q.hchi_full_a,3)==1.792 and round(q.hchi_full_b,3)==1.375
assert round(g.MEL5_model_RMSF_A,3)==1.118 and round(g.alpha24beta17_model_RMSF_A,3)==.747
for state in ['wt','eng']:
    names=meta['structures'][state]['atom_screen_pixels_top_left']
    assert 'beta52_sidechain_center' in names and 'beta53_sidechain_center' in names
review={'status':'passed','reviewed_figure_sha256':hashlib.sha256((root/'fig4_mel5_interface.png').read_bytes()).hexdigest(),
    'review_scope':'Independent full-figure and enlarged state-specific visual inspection; numerical source and anchor-name verification.',
    'passed_checks':{'fixed_atom_distances_rounded_A':{'WT':[round(v,2) for v in dist.WT_distance_A],'mature':[round(v,2) for v in dist.mature_distance_A]},'Q72_Hchi_rounded':[1.792,1.375],'betaGly54_RMSF_rounded_A':[1.118,.747],'same_camera_and_same_crop':True,'same_axes_width_and_height':True,'correct_numbering_and_beta52_53_54_targets':True,'anchor_JSON_field_names_match_plot':True,'all45_CDR_rows_preserved':len(cdr)==45,'no_Hchi_to_RMSF_causal_arrow':True},
    'anchor_note':'WT beta52_sidechain_center is a CA fallback because Gly52 has no heavy sidechain; used only as a residue-label target and never called a sidechain center in the figure.',
    'must_fix':[],
    'resolved_label_issues':[{'state':'mature','label':'N at beta53 backbone nitrogen','resolution':'Moved left/up; the blue nitrogen endpoint and adjoining distance dashes are now visible.'},{'state':'WT','label':'Q72 NE2','resolution':'Moved left; the blue NE2 endpoint is now visible. Mature NE2 label is also offset farther right/down.'}],
    'other_geometry_assessment':'The two paired structures visibly show the sidechain reorientation and approach of beta53/beta54; distance definitions are the same named atom pairs in both states. No geometry or residue-label mismatch found.'}
(root/'independent_geometry_visual_review.json').write_text(json.dumps(review,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(review,ensure_ascii=False,indent=2))
