"""Final traceability checks after manual view_image inspection of both masters."""
from pathlib import Path
import sys,json,hashlib,re,shutil,xml.etree.ElementTree as ET
BASE=Path(r'E:\Project\FlexPPI\experiments\_sources')
sys.path.insert(0,str(BASE/'completed_results_20261012/plot_runtime'))
import numpy as np
import pandas as pd
from PIL import Image
OUT=BASE/'candidate_case_atlas_20261021/evidence/non_tcr_figures'
OLD=BASE/'reader_revision_v39_20261016/english_figures/fig6'
UBOLD=BASE/'figure_merge_v43_20261018/figure5/source/original_fig7'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:pd.read_csv(p,float_precision='round_trip')
protocol=json.loads((OLD/'case_model_protocol.json').read_text(encoding='utf-8'))
assert protocol['KIX']['fold']==1 and protocol['ubiquitin']['fold']==0

kix=read(OLD/'data/kix_fold1_all_residues.csv');kix=kix[kix.metric=='rmsf']
km=read(OLD/'data/kix_fold1_case_metrics.csv');km=km[km.metric=='rmsf'].set_index('state')
kp=read(OUT/'C31/all_plotted_points.csv')
assert len(kp)==1400
kixchecks=[]
for state,src in kix.groupby('state'):
    src=src.sort_values('author_residue_number')
    assert len(src)==87
    for panel,kind,key in [('c','MD','md_percentile'),('c','model','model_percentile'),('e','MD','md_value'),('e','model','model_value')]:
        p=kp[(kp.panel==panel)&(kp.state==state)&(kp.kind==kind)].sort_values('residue')
        np.testing.assert_array_equal(p.residue.to_numpy(),src.author_residue_number.to_numpy())
        np.testing.assert_allclose(p.y,src[key],rtol=0,atol=1e-13)
        kixchecks.append(dict(panel=panel,state=state,kind=kind,n=len(p),max_error=float(np.max(np.abs(p.y.to_numpy()-src[key].to_numpy())))))
    truth=set(src.nlargest(18,'md_value').author_residue_number)
    pred=set(src.nlargest(18,'model_value').author_residue_number)
    assert len(truth&pred)==int(km.loc[state,'top_overlap'])
    l12=src[src.author_residue_number.between(614,621)]
    assert len(l12)==8
    for key in ['md','model']:
        v=float(l12[key+'_percentile'].mean())
        assert abs(v-float(km.loc[state,'l12_'+key+'_mean_percentile']))<1e-12

rows=read(UBOLD/'source/all66_residue_summary.csv');defined=read(UBOLD/'experiment/defined_angles_long.csv');hist=pd.read_csv(UBOLD/'source/histograms_501frames.csv')
up=read(OUT/'C32/all_plotted_points.csv');assert len(up)==715
ubchecks=[]
for kind,key in [('MD','frozen_MD_joint_Hchi_nats'),('model','predicted_Hchi_nats')]:
    p=up[(up.panel=='b')&(up.kind==kind)].sort_values('residue_id');src=rows.sort_values('author_residue_id')
    np.testing.assert_array_equal(p.residue_id,src.author_residue_id)
    np.testing.assert_allclose(p.y,src[key],rtol=0,atol=1e-13)
    ubchecks.append(dict(panel='b',kind=kind,n=len(p),source_error_max=float(np.max(np.abs(p.y.to_numpy()-src[key].to_numpy())))))
for chi,n in [(1,66),(2,52),(3,23),(4,10)]:
    p=up[(up.panel=='c')&(up.chi==chi)].sort_values('residue_id');src=defined[defined.chi=='chi'+str(chi)].sort_values('author_residue_id')
    assert len(p)==len(src)==n
    np.testing.assert_array_equal(p.residue_id,src.author_residue_id)
    np.testing.assert_allclose(p.x,src.predicted_Hchi_nats_AA_centered,rtol=0,atol=1e-13)
    np.testing.assert_allclose(p.y,src.fine_entropy_four_origins_mean_nats_AA_centered,rtol=0,atol=1e-13)
    ubchecks.append(dict(panel='c',chi=chi,n=n,all_positions_retained=True))
for site in [2,41]:
    for chi in [1,2,3]:
        p=up[(up.panel=='d')&(up.residue_id==site)&(up.chi==chi)].sort_values('x')
        src=hist[(hist.residue=='Q'+str(site))&(hist.chi=='chi'+str(chi))].sort_values('angle_center_degrees')
        assert len(p)==len(src)==72
        np.testing.assert_allclose(p.x,src.angle_center_degrees,rtol=0,atol=1e-13)
        np.testing.assert_allclose(p.y,src.F_star,rtol=0,atol=1e-13,equal_nan=True)
        good=np.isfinite(src.F_star)&(src.F_star<=1)
        assert list(p.low_Fstar_le1)==list(good)
        ubchecks.append(dict(panel='d',residue='Q'+str(site),chi=chi,n=72,zero_count_bins=int(src.F_star.isna().sum()),lowF_width_degrees=int(good.sum())*5))
gln=rows[rows.amino_acid=='Q'].sort_values('predicted_Hchi_nats')
assert len(gln)==6 and int(gln.iloc[0].author_residue_id)==41 and int(gln.iloc[-1].author_residue_id)==2
assert round(float(gln.iloc[0].predicted_Hchi_nats),3)==2.000
assert round(float(gln.iloc[-1].predicted_Hchi_nats),3)==2.990

manifest=json.loads((OUT/'non_tcr_figure_manifest.json').read_text(encoding='utf-8'))
for caseid,checks,images,domain in [('C31',kixchecks,5,'KIX only'),('C32',ubchecks,1,'Ubiquitin–UBA only')]:
    dest=OUT/caseid;stem='fig_'+caseid.lower()+'_case'
    stats=json.loads((dest/'statistics_and_data_audit.json').read_text(encoding='utf-8'))
    current={ext:sha(dest/(stem+'.'+ext)) for ext in ['png','pdf','svg']}
    assert current==stats['output_sha256']
    svg=(dest/(stem+'.svg')).read_text(encoding='utf-8');tree=ET.fromstring(svg)
    assert len(tree.findall('.//{http://www.w3.org/2000/svg}image'))==images
    assert len(tree.findall('.//{http://www.w3.org/2000/svg}path'))>20
    fonts=[float(x) for x in re.findall(r'font-size:\s*([0-9.]+)px',svg)]
    if not fonts:fonts=[float(x) for x in re.findall(r'font:\s*(?:[0-9]+\s+)?([0-9.]+)px',svg)]
    assert fonts and min(fonts)>=9,(caseid,fonts)
    with Image.open(dest/(stem+'.png')) as im:
        assert im.size==(3242,3631),im.size
        png_size=list(im.size)
    assert stats['text_outside_canvas']==[]
    visual=dict(status='passed',date='2026-10-21',reviewer='case_evidence_20261021',method='Individually opened final PNG using view_image after the final export',
                reviewed_output_sha256=current,canvas_mm=[183,205],png_pixels=png_size,png_dpi=450,minimum_SVG_font_pt=min(fonts),
                language='English',checks={'VV1_no_text_or_legend_occlusion':'pass','VV2_panel_alignment_and_canvas':'pass','VV3_9pt_text_legible':'pass','VV4_semantic_color_consistency':'pass','VV5_complete_source_data_visible':'pass'},
                fixes=['Replaced crowded KIX d abbreviations plus redundant explanation by directly named rotated state labels.'] if caseid=='C31' else [],
                raster_content='Three real molecular renders, the complete87x8 percentile heatmap and its color scale; curves/text are vector.' if caseid=='C31' else 'Only one real molecular render; all curves/scatter/text are vector. Molecular text removed and relabeled in vector.')
    (dest/'visual_QA.json').write_text(json.dumps(visual,ensure_ascii=False,indent=2),encoding='utf-8')
    content=dict(status='passed',case_id=caseid,system=domain,output_sha256=current,direct_plot_source_checks=checks,
                 all_observations_retained=True,no_new_model_training=True,no_new_inference=True,no_new_MD_production=True,
                 model_protocol=protocol['KIX' if caseid=='C31' else 'ubiquitin'],
                 interpretation='Four complete87-position KIX profiles and localization, not recovered state differences or cooperation.' if caseid=='C31' else 'All66 scalar joint-Hchi values, all151 per-angle associations and all432 actual MD bins, not direct predicted angle distributions or binding free energy.',
                 numerical_source_max_tolerance=1e-13,caption_path=str(dest/'caption_zh.md'))
    (dest/'final_content_review.json').write_text(json.dumps(content,ensure_ascii=False,indent=2),encoding='utf-8')
    stats['status']='passed';stats['visual_QA']='passed';stats['direct_source_check_count']=len(checks);stats['source_SVG_font_min_pt']=min(fonts)
    (dest/'statistics_and_data_audit.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf-8')
    for row in manifest['figures']:
        if row['caseid']==caseid:row.update(status='complete',visual_QA=str(dest/'visual_QA.json'),final_content_review=str(dest/'final_content_review.json'),output_sha256=current)
archive=OUT/'source_archive';archive.mkdir(exist_ok=True)
for src in [OLD/'data/kix_fold1_all_residues.csv',OLD/'data/kix_fold1_case_metrics.csv',OLD/'case_model_protocol.json',OLD/'compose_figure.py',UBOLD/'source/all66_residue_summary.csv',UBOLD/'source/histograms_501frames.csv',UBOLD/'experiment/defined_angles_long.csv',UBOLD/'experiment/per_chi_correlations.csv',OLD/'ubiquitin/inputs/A003R_input_structure.pdb',OLD/'kix/prepared_2AGH_model1_ternary.pdb']:
    shutil.copy2(src,archive/src.name)
for name in ['render_non_tcr_structures.py','finalize_non_tcr_figures.py','non_tcr_case_figures_body.py']:
    shutil.copy2(Path(__file__).parent/name,OUT/name)
manifest.update(status='complete',date='2026-10-21',n_figures=2,final_visual_and_data_QA='passed',source_archive_sha256={p.name:sha(p) for p in archive.iterdir() if p.is_file()})
(OUT/'non_tcr_figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(status='complete',figures=2,source_checks={'C31':len(kixchecks),'C32':len(ubchecks)}),ensure_ascii=False))
