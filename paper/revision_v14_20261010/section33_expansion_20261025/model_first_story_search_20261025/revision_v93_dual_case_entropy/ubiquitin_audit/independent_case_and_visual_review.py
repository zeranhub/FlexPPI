"""Independent original-cache checks of final plotted ubiquitin values."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from PIL import Image

O=Path(__file__).resolve().parent;R=O.parent
P=R/'figure3'/'source_data'/'all_plotted_values.csv'
PNG=R/'figure3'/'fig3_dual_case_entropy.png'
T=R/'manuscript'/'final_text_contract.json'
ORIGINAL=Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\experiments\corrected_pipeline\neural_exports\evaluation\full\A003R.npz')
z=np.load(ORIGINAL,allow_pickle=False);plots=pd.read_csv(P,float_precision='round_trip');text=json.loads(T.read_text(encoding='utf-8'))
source_idx=np.flatnonzero((z['chain']=='A')&z['valid'][:,1]);assert len(source_idx)==66
ids=z['author_residue_id'][source_idx].astype(int)
checks={}
for label,column in [('MD','target'),('Full','mean')]:
    plotted=plots[(plots.panel=='c')&(plots.series==label)]
    plotted=plotted.assign(residue_id=plotted.identity.astype(int)).sort_values('residue_id')
    assert len(plotted)==66 and not plotted.residue_id.duplicated().any()
    np.testing.assert_array_equal(plotted.residue_id,ids)
    diff=float(np.max(np.abs(plotted.value.to_numpy()-z[column][source_idx,1])))
    assert diff<1e-14
    checks['c_'+label]={'n':66,'raw_cache_max_abs_difference':diff,'all_original_valid_positions_exact':True}
env=pd.read_csv(O/'all119_input_environment.csv',float_precision='round_trip')
six=np.flatnonzero((z['chain']=='A')&(z['amino_acid']=='Q'))
for label,column in [('MD','target'),('Full','mean')]:
    plotted=plots[(plots.panel=='e')&(plots.series==label)]
    assert len(plotted)==6 and set(plotted.identity)=={'Q2','Q31','Q40','Q41','Q49','Q62'}
    diffs=[];sdiffs=[]
    for row in plotted.itertuples():
        resid=int(row.identity[1:]);idx=np.flatnonzero((z['chain']=='A')&(z['author_residue_id']==str(resid)))
        assert len(idx)==1
        diffs.append(abs(row.value-float(z[column][idx[0],1])))
        sasa=env[(env.chain=='A')&(env.input_residue_id==resid)].sidechain_SASA_A2.iloc[0]
        sdiffs.append(abs(row.SASA_A2-sasa))
    assert max(diffs)<1e-14 and max(sdiffs)<1e-12
    checks['e_'+label]={'n':6,'raw_cache_max_abs_difference':max(diffs),'independent_SASA_max_abs_difference':max(sdiffs),'all_six_no_omission':True}
pro=[]
for resid in [19,37,38]:
    idx=int(np.flatnonzero((z['chain']=='A')&(z['author_residue_id']==str(resid)))[0])
    a=plots[(plots.panel=='c')&plots.identity.eq(str(resid))]
    assert len(a)==2 and set(a.series)=={'MD','Full'} and bool(z['valid'][idx,1])
    pro.append({'residue':'P'+str(resid),'source_valid':True,'source_MD_Hchi':float(z['target'][idx,1]),
       'source_Full_Hchi':float(z['mean'][idx,1]),'main_plot_retained':True})
stats={'all66_MAE':float(abs(z['target'][source_idx,1]-z['mean'][source_idx,1]).astype(float).mean()),
  'all66_rho':float(spearmanr(z['target'][source_idx,1],z['mean'][source_idx,1]).statistic),
  'Gln6_MAE':float(abs(z['target'][six,1]-z['mean'][six,1]).astype(float).mean()),
  'Gln6_rho':float(spearmanr(z['target'][six,1],z['mean'][six,1]).statistic)}
manifest=json.loads((R/'pymol'/'render_manifest.json').read_text(encoding='utf-8'))
a=manifest['assets']['ub_Q41_local'];b=manifest['assets']['ub_Q62_local']
np.testing.assert_allclose(a['view'][:9],b['view'][:9],rtol=0,atol=0)
assert a['shared_optical_halfheight_A']==b['shared_optical_halfheight_A']
checks['d_native_structure']={'shared_camera_rotation':True,'shared_optical_halfheight_A':a['shared_optical_halfheight_A'],
 'Q41_whole_neighbor_residues':a['whole_residues'],'Q62_whole_neighbor_residues':b['whole_residues'],
 'Q41_heavy_sticks':a['stick_atom_count'],'Q62_heavy_sticks':b['stick_atom_count'],
 'all_complete_selected_residues_shown':a['all_whole_residue_heavy_atoms_shown_as_sticks'] and b['all_whole_residue_heavy_atoms_shown_as_sticks'],
 'both_unchanged_source_coordinates':a['unchanged_source_coordinates'] and b['unchanged_source_coordinates']}
phrasing='支持局部堆积与溶剂暴露共同影响侧链取向多样性的解释'
result={'reviewer':'independent ubiquitin case subagent','status':'NUMERIC_AND_VISUAL_PASS_WITH_ONE_WORDING_SUGGESTION',
 'reviewed_files':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [PNG,P,T,ORIGINAL]},
 'actual_PNG_opened_and_reviewed':True,'PNG_dimensions_px':Image.open(PNG).size,
 'checks':checks,'Pro_original_valid_zero_positions':pro,'stats_recomputed_from_original_cache':stats,
 'visual_review':{'target_Q41_Q62_labels_correct':True,'c_retains_all66_including_three_Pro_zeros':True,
   'no_Ub_RMSF_panel':True,'all6_Gln_present_without_regression':True,'neighborhoods_visually_complete_and_target_emphasized':True,
   'no_observed_text_clipping_or_panel_overlap':True,
   'minor_design_observation':'c full profile is visually dense but accurately exposes all positions; e no separate legend is acceptable only with shared c symbol mapping and explicit caption'},
 'text_scope_review':{'no_apo_holo_reinterpretation':True,'Ub_original_fold_not_mislabeled_family_OOF':True,
   'nonmonotonic_SASA_caveat_present':True,'no_direct_Q41_Q62_UBA_effect_claim':True,'model_does_not_generate_MD_distributions':True,
   'single_wording_suggestion':{'present':phrasing in '\n'.join(text['section33']['before_figure']),
      'current':phrasing,'why':'共同影响 carries a causal interpretation stronger than this single static-structure/MD comparison; surrounding supports/explanation language softens it but safer wording is available',
      'suggested':'这一差异与紧密堆积约束可访问侧链状态的解释一致，并表明模型能够恢复同类型侧链在不同局部环境中的状态多样性差异'}},
 'scope_limit':'Independent numerical audit here covers Ub panels c/e and metadata/actual image review of d; CaM and 83-record panels were viewed for layout but their numerical sources were not independently recomputed in this review.',
 'original_main_files_modified':False}
(O/'independent_case_and_visual_review.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
