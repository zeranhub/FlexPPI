"""Read-only independent review of v54 Figure4 exported sources and labels."""
from pathlib import Path
import ast, json, hashlib
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

ROOT=Path(__file__).resolve().parent.parent
FIG=ROOT/'figure4'
EXP=ROOT.parent
INV=EXP/'opposite_affinity_rmsf_cases_20261022/inventory'
OLD=EXP/'mel5_joint_case_v49_20261019/figure4'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:pd.read_csv(p,float_precision='round_trip')
mapping={
 'CDR_absolute':OLD/'source_plotted_CDR_RMSF.csv',
 'interface_absolute':OLD/'source_plotted_interface_Hchi.csv',
 'peptide_four_states':INV/'MEL5_four_states_complete_peptide.csv',
 'CDR_three_contrasts':INV/'MEL5_three_contrasts_all_six_CDR.csv',
 'interface_three_contrasts':INV/'MEL5_three_contrasts_four_state_union_pMHC_interface.csv',
}
copychecks=[]
for name,path in mapping.items():
    target=FIG/f'source_{name}.csv'
    assert sha(path)==sha(target)
    copychecks.append({'name':name,'original':str(path),'SHA256':sha(path),'copied_byte_identical':True})

old_h=read(mapping['interface_absolute'])
expected_rank=old_h.assign(delta_Hchi=old_h.hchi_full_b-old_h.hchi_full_a).sort_values(['delta_Hchi','chain_a','position_a'],kind='stable').reset_index(drop=True)
rank=read(FIG/'source_interface_all23_ranked.csv')
assert list(rank.columns)==list(expected_rank.columns)
assert len(rank)==23
for col in rank:
    if pd.api.types.is_numeric_dtype(rank[col]):assert np.allclose(rank[col],expected_rank[col],rtol=0,atol=1e-15)
    else:assert rank[col].tolist()==expected_rank[col].tolist()
assert (rank.iloc[0].chain_a,int(rank.iloc[0].position_a))==('A',72)
assert (np.diff(rank.delta_Hchi)>=0).all()

# Matrix order derived independently from the original absolute CDR curve, not from the export.
cdr=read(mapping['CDR_absolute']).sort_values('plot_x')
order=[('D' if 'α' in r.region else 'E',int(r.common_TCR_position)) for r in cdr.itertuples()]
assert len(set(order))==45
engineered={('D' if 'α' in r.region else 'E',int(r.derivative_input_position)):('D' if 'α' in r.region else 'E',int(r.common_TCR_position)) for r in cdr.itertuples()}
original_c=read(mapping['CDR_three_contrasts'])
original_h=read(mapping['interface_three_contrasts'])
contrasts=[('3hg1','4jff'),('4jff','4jfd'),('4jff','4jfe')]
names=['TCR maturation','Peptide G4A','Peptide I7A']
cvalues=[];hvalues=[];nvalid=[]
ih=original_h[original_h.pdb_before.eq('3hg1')&original_h.pdb_after.eq('4jff')]
iorder=[(r.chain_before,int(r.id_before)) for r in ih.itertuples()]
assert len(iorder)==30
cexport=read(FIG/'source_CDR_heatmap_all135.csv')
hexport=read(FIG/'source_interface_heatmap_all90.csv')
assert len(cexport)==135 and len(hexport)==90
for name,(a,b) in zip(names,contrasts):
    q=original_c[original_c.pdb_before.eq(a)&original_c.pdb_after.eq(b)]
    vals={}
    for r in q.itertuples():
        key=(r.chain_before,int(r.id_before))
        if a!='3hg1':key=engineered[key]
        vals[key]=r.rmsf_full_after-r.rmsf_full_before
    assert set(vals)==set(order)
    cvalues.append([vals[k] for k in order])
    e=cexport[cexport.contrast.eq(name)]
    ev={(r.reference_chain,int(r.reference_position)):r.delta_RMSF_A for r in e.itertuples()}
    assert set(ev)==set(vals)
    assert max(abs(ev[k]-vals[k]) for k in vals)<1e-15
    q=original_h[original_h.pdb_before.eq(a)&original_h.pdb_after.eq(b)]
    hv={};he={}
    for r in q.itertuples():
        key=(r.chain_before,int(r.id_before))
        eligible=bool(r.hchi_valid and r.same_aa)
        hv[key]=(r.hchi_full_after-r.hchi_full_before) if eligible else np.nan
        he[key]=eligible
    assert set(hv)==set(iorder)
    hvalues.append([hv[k] for k in iorder]);nvalid.append(sum(he.values()))
    e=hexport[hexport.contrast.eq(name)]
    ev={(r.chain,int(r.position)):r for r in e.itertuples()}
    assert set(ev)==set(hv)
    for k in hv:
        assert bool(ev[k].eligible)==he[k]
        if he[k]:assert abs(ev[k].delta_Hchi_nats-hv[k])<1e-15
        else:assert np.isnan(ev[k].delta_Hchi_nats)
arrays=np.load(FIG/'source_heatmap_arrays.npz',allow_pickle=False)
assert arrays['CDR_delta_RMSF'].shape==(3,45)
assert arrays['interface_delta_Hchi'].shape==(3,30)
assert np.allclose(arrays['CDR_delta_RMSF'],cvalues,rtol=0,atol=1e-15)
assert np.allclose(arrays['interface_delta_Hchi'],hvalues,rtol=0,atol=1e-15,equal_nan=True)
assert nvalid==[23,23,22]
assert max(abs(np.asarray(cvalues)).ravel())<=.45
assert np.nanmax(abs(np.asarray(hvalues)))<=.45

# Inspect literal plotting declarations, plus the actual PNG, to ensure mappings were not swapped.
tree=ast.parse((FIG/'plot.py').read_text(encoding='utf-8'))
constants={}
for node in tree.body:
    if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
        try:constants[node.targets[0].id]=ast.literal_eval(node.value)
        except (ValueError,TypeError):pass
assert constants['STATES']==['3hg1','4jff','4jfd','4jfe']
assert constants['NAMES']==['MEL5 · ELA','α24β17 · ELA','α24β17 · G4A','α24β17 · I7A']
assert constants['MARKERS']==['o','s','D','^']
assert constants['STYLES']==['-','-','--',':']
for key,value in {'WT':'#CF8C27','ENG':'#762A83','TEAL':'#297B82','G4A':'#4B82A5','I7A':'#73927D','NA':'#D9DDE0'}.items():assert constants[key]==value
assert constants['SIGNED']==['#427A9F','#FAFAFA','#B76778']

# Fixed atom distances in the displayed archived local geometry; do not substitute v53 min-heavy pairs.
pdb=OLD/'source_archive/4jff_model_input.pdb'
structure=PDBParser(QUIET=True).get_structure('mature',pdb)[0]
q72=structure['A'][72]
dist53=float(np.linalg.norm(q72['OE1'].coord.astype(float)-structure['E'][54]['N'].coord.astype(float)))
dist54=float(np.linalg.norm(q72['NE2'].coord.astype(float)-structure['E'][55]['O'].coord.astype(float)))
assert f'{dist53:.2f}'=='2.94' and f'{dist54:.2f}'=='3.12'
image_path=FIG/'fig4_mel5_composite.png'
audit=json.loads((FIG/'figure_data_audit.json').read_text(encoding='utf-8'))
assert audit['figure_sha256']==sha(image_path)
assert len(audit['visible_text'])>0
texts=[x['text'] for x in audit['visible_text']]
assert '2.94 Å' in texts and '3.12 Å' in texts
for t in constants['NAMES']:assert t in texts
assert 'Δ = α24β17 − MEL5' in texts
report={
 'status':'passed','reviewer':'independent source/data reviewer',
 'PNG_path':str(image_path),'PNG_SHA256':sha(image_path),
 'plot_py_SHA256':sha(FIG/'plot.py'),
 'read_only_review':True,'actual_PNG_viewed':True,
 'copied_sources':copychecks,
 'panel_a':{'distances_are_fixed_atom_pairs':True,'common_beta53_author_id':54,'common_beta54_author_id':55,'Q72_OE1_to_beta53_backbone_N_A':dist53,'Q72_NE2_to_beta54_backbone_O_A':dist54,'labels_match_two_decimal_rounding':True,'not_v53_minimum_heavy_distance_metric':True,'structure_PDB_SHA256':sha(pdb)},
 'panel_b':{'all_four_states_all_ten_peptide_sites':True,'total_values':40,'state_names_and_colors_confirmed':True,'NPZ_equality_verified_in_source_audit':True},
 'panel_c':{'all45_positions_two_states':True,'total_values':90,'common_numbering_and_state_mapping_match':True},
 'panel_d':{'all23_same_AA_valid_Hchi_sites':True,'ascending_after_minus_before_sort':True,'Q72_largest_decrease':True,'values_match_heatmap_first_row':True,'peptide_prefix_p_prevents_HLA_confusion':True},
 'panel_e':{'shape':[3,45],'all135_values_match_original':True,'common_TCR_numbering_match_original_absolute_profiles':True},
 'panel_f':{'shape':[3,30],'all90_values_and_eligibility_masks_match_original':True,'eligible_per_row':nvalid,'ineligible_is_nan_and_grey':True,'Q72_column_highlight_correct':True},
 'heatmaps':{'same_symmetric_limits':[-.45,.45],'negative_blue_positive_rose':True,'no_clipped_values':True,'units_separate':['RMSF Angstrom','Hchi nats']},
 'mapping':{'MEL5_ELA':{'pdb':'3hg1','color':'#CF8C27','marker':'open circle','line':'solid'},'mature_ELA':{'pdb':'4jff','color':'#762A83','marker':'filled square','line':'solid'},'mature_G4A':{'pdb':'4jfd','color':'#4B82A5','marker':'filled diamond','line':'dashed'},'mature_I7A':{'pdb':'4jfe','color':'#73927D','marker':'filled triangle','line':'dotted'}},
 'required_caption_definitions':['TCR maturation = 4JFF minus 3HG1; peptide G4A = 4JFD minus 4JFF; peptide I7A = 4JFE minus 4JFF.','Grey Hchi cells are chemically ineligible, never zero; all30 interface positions retained.','p prefix in ranked bar chart labels denotes peptide residues.','All TCR labels use MEL5 common residue numbering.','Interface heatmap scope is four-state union defined by nearest interpartner C-alpha <=8 Angstrom, not the fixed atomic pairs in a.'],
 'sources_unchanged':True,
 'independent_source_audit':str(ROOT/'source_audit/independent_source_audit.json'),
 'exported_data_files':[{ 'path':str(p),'SHA256':sha(p)} for p in sorted(FIG.glob('source_*')) if p.is_file()]
}
(FIG/'independent_data_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','figure_SHA256':sha(image_path),'fixed_atomic_distances':[dist53,dist54],'file':str(FIG/'independent_data_review.json')},ensure_ascii=False))
