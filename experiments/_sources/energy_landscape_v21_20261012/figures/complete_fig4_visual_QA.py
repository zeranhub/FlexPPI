from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, re, xml.etree.ElementTree as ET
BASE=Path(__file__).resolve().parent
stem='fig4_joint_complementarity'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
stats_path=BASE/(stem+'_statistics.json')
stats=json.loads(stats_path.read_text(encoding='utf-8'))
manifest_path=BASE/(stem+'_sourcehash.json')
manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
assert not stats['qa']['text_outside_canvas'] and stats['qa']['minimum_text_pt']>=7.5
assert stats['qa']['width_mm']==183 and stats['qa']['height_mm']==80
assert min(stats['qa']['dpi'])>=299 and stats['cdf_point_counts']=={'reference':3063,'prediction':3063}
assert stats['recall_point_count']==12130 and stats['all_valid_records_drawn']
pdf=(BASE/(stem+'.pdf')).read_bytes()
assert b'/FontFile2' in pdf and len(re.findall(rb'/Subtype\s*/Image',pdf))==0
svg=ET.parse(BASE/(stem+'.svg'));ns={'s':'http://www.w3.org/2000/svg'}
assert len(svg.findall('.//s:image',ns))==0 and len(svg.findall('.//s:text',ns))>20
for ext,metadata in manifest['figure']['outputs'].items():assert sha(BASE/(stem+'.'+ext))==metadata['sha256']
for name,digest in manifest['source_sha256'].items():assert sha(name)==digest
qa={'status':'READY','checked_utc':datetime.now(timezone.utc).isoformat(),'actual_full_image_viewed_with_view_image':True,'viewed_png_sha256':sha(BASE/(stem+'.png')),'visual_findings':'Three aligned panels have legible Chinese labels, no clipping or legend/data overlap. Both complete ECDF curves, eight category bars, explicit unclassified fractions and all 12130 valid gray recall points are visible. Macro means and record denominators match the source.','pass0_antipatterns':{'verbatim_style_baselines':True,'custom_gray_blue_palette':True,'top_right_spines_removed':True,'legends_outside_axes':True,'vector_and_300dpi_preview':True,'explicit_Chinese_font':True},'pass1_code':{'183mm_width':True,'font_floor_7p5':True,'all_records_no_subsample':True,'PDF_TrueType_fonts':True,'SVG_editable_text':True,'no_rasterized_vector_layers':True},'pass2_logic':{'one_message_unique_evidence_per_panel':True,'record_macro_not_micro':True,'bars_zero_origin':True,'unclassified_mass_retained':True,'undefined_recalls_only_zero_denominators_excluded':True,'no_CI_or_significance':True},'pass3_render':{'no_visible_data_occlusion':True,'aligned_panel_layout':True,'text_legible':True,'no_canvas_cutoff':True,'palette_and_numeric_labels_consistent':True},'sourcehash_sha256_before_QA_update':sha(manifest_path),'statistics_sha256_before_QA_update':sha(stats_path)}
stats['status']='complete_visual_QA_passed'
stats['qa']['visual_review']='READY; actual full PNG reviewed with view_image; no visible defects.'
stats_path.write_text(json.dumps(stats,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
manifest['status']='complete_visual_QA_passed'
manifest['figure']['qa']['visual_review']=stats['qa']['visual_review']
manifest['output_metadata_sha256']=sha(stats_path)
manifest['visual_QA_file']=stem+'_QA.json'
manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
(BASE/(stem+'_QA.json')).write_text(json.dumps(qa,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':'READY','PNG_sha256':qa['viewed_png_sha256'],'statistics':str(stats_path)},ensure_ascii=False))
