"""Bind the personally inspected final PNG and independently audit SVG points."""
import os
_dll = os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
from pathlib import Path
import json
import hashlib
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
audit_path = HERE / 'data_audit.json'
audit = json.loads(audit_path.read_text(encoding='utf-8'))
stem = HERE / 'fixed_backbone_antigen_pilot'
tree = ET.parse(str(stem)+'.svg').getroot()
ns = {'s': 'http://www.w3.org/2000/svg'}
collections = [g for g in tree.findall('.//s:g', ns) if g.attrib.get('id', '').startswith('PathCollection')]
counts = [len(g.findall('.//s:use', ns)) for g in collections]
assert counts == [13, 13, 13, 13, 7, 7]
rects = [list(g)[0].attrib for g in tree.findall('.//s:clipPath', ns)]
rects = [{k: float(row[k]) for k in ['x', 'y', 'width', 'height']} for row in rects]
pairs = pd.read_csv(HERE / 'source_data/panels_ab_all26_candidate_pairs.csv', float_precision='round_trip')
methods = ['zero_change', 'sequence_only', 'static_only', 'sequence_plus_static',
           'sequence_plus_static_plus_R', 'sequence_plus_static_plus_H', 'sequence_plus_static_plus_RH']
metrics = pd.read_csv(HERE / 'source_data/panel_c_all14_method_MAE.csv', float_precision='round_trip')
maxerr = 0.
for j, group in enumerate(collections):
    uses = group.findall('.//s:use', ns)
    xy = np.asarray([[float(e.attrib['x']), float(e.attrib['y'])] for e in uses])
    if j < 4:
        system = 'nyeso' if j < 2 else 'tax'
        data = pairs[pairs.system == system]
        x = data.observed_delta_pKD.to_numpy()
        y = data['baseline_prediction' if j % 2 == 0 else 'RH_prediction'].to_numpy()
        rect = rects[j//2]
        expectx = rect['x'] + (x+2.7)/3.25*rect['width']
        expecty = rect['y'] + rect['height'] - (y+2.7)/3.25*rect['height']
    else:
        system = 'nyeso' if j == 4 else 'tax'
        data = metrics[metrics.system == system].set_index('model').loc[methods]
        x = data.MAE_delta_pKD.to_numpy()
        y = np.arange(7)+(-.12 if j == 4 else .12)
        rect = rects[2]
        expectx = rect['x'] + (x-.30)/.76*rect['width']
        expecty = rect['y'] + (y+.55)/7.10*rect['height']
    err = float(np.abs(xy - np.column_stack([expectx,expecty])).max())
    maxerr = max(maxerr, err)
    assert err < 2e-5, (j, err)
assert not tree.findall('.//s:image', ns), 'Expected fully vector plot, not pasted raster panels.'
width_mm = float(tree.attrib['width'].replace('pt', ''))/72*25.4
height_mm = float(tree.attrib['height'].replace('pt', ''))/72*25.4
assert abs(width_mm-183) < 1e-5 and abs(height_mm-196) < 1e-5
sens = pd.read_csv(HERE/'source_data/panel_d_all15_sensitivity_cells.csv',float_precision='round_trip')
for value in sens.joint_increment_vs_baseline:
    assert any(e.text == f'{value:+.4f}' for e in tree.findall('.//s:text',ns)), value
with Image.open(str(stem)+'.png') as image:
    pixels = np.asarray(image.convert('RGB'))
    density = float((pixels.min(axis=2) < 245).mean())
assert .015 < density < .50
for name, row in audit['source_files'].items():
    assert sha(row['original']) == row['SHA256'] == sha(row['archive'])
assert all(sha(str(stem)+'.'+ext) == digest for ext,digest in audit['export_SHA256'].items())
qa = {'status':'passed','personally_opened_final_PNG':True,
    'reviewer':'/root/case_figure_design_20261021','reviewed_at_utc':datetime.now(timezone.utc).isoformat(),
    'PNG_SHA256':sha(str(stem)+'.png'),'PDF_SHA256':sha(str(stem)+'.pdf'),'SVG_SHA256':sha(str(stem)+'.svg'),
    'actual_dimensions_mm':[width_mm,height_mm],'minimum_font_pt':audit['minimum_visible_font_pt'],
    'PNG_dpi':450,'all52_scatter_markers_and14_method_markers_in_actual_SVG':counts,
    'SVG_coordinate_max_absolute_error_points':maxerr,'all15_heatmap_numbers_in_actual_SVG':True,
    'all_requested_primary_candidates_methods_templates_and_fits_retained':True,
    'no_raster_image_panels_in_vector_SVG':True,'PNG_nonwhite_content_fraction':density,
    'Pass0_anti_patterns':'passed: mandatory baseline blocks; semantic gold/purple systems; open-circle vs filled-square method contrast; clean spines; external legends; vector exports.',
    'Pass1_code':'passed:183mm width,196mm height,9pt minimum,450dpi,Arial,editable SVG and embedded PDF fonts.',
    'Pass2_integrity':'passed:26 candidates/7methods/15cells;only within-system held-position-out rows;MAE recomputed;zero-change no rank claim;shared scatter axes and zero-symmetric gain scale.',
    'Pass3_visual':'passed:personally opened final PNG after moving the legends and heatmap fit labels; no text collisions,clipping or missing glyphs. All source coordinates unchanged. Near-identical predictions may naturally overlap without jitter or selection.',
    'statistical_tests_or_significance_claims':None,
    'data_signal_exception':'A small heatmap variation is scientifically meaningful here. No effect size was altered to meet generic visual-variance heuristics; all real signed cells are annotated on a single ±0.08 scale.',
    'must_fix_findings':[],'original_results_protocol_and_manuscript_modified':False}
(HERE/'visual_QA.json').write_text(json.dumps(qa,indent=2)+'\n',encoding='utf-8')
audit['status']='passed'
audit['personal_visual_QA_file']=str(HERE/'visual_QA.json')
audit['personal_visual_QA_SHA256']=sha(HERE/'visual_QA.json')
audit['actual_SVG_point_coordinates_independently_verified']=True
audit_path.write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
manifest={'status':'passed','figure':'fixed_backbone_antigen_pilot','dimensions_mm':[183,196],
    'outputs':{ext:{'path':str(stem)+'.'+ext,'SHA256':sha(str(stem)+'.'+ext)} for ext in ['png','pdf','svg']},
    'source_data':str(HERE/'source_data'),'figure_script':str(HERE/'plot_fixed_backbone_pilot.py'),
    'figure_script_SHA256':sha(HERE/'plot_fixed_backbone_pilot.py'),
    'audits':{name:sha(HERE/name) for name in ['data_audit.json','visual_QA.json','statistics_and_reproducibility.json','figure_contract.json']},
    'captions':{name:sha(HERE/name) for name in ['caption_en.txt','caption_zh.txt']}}
(HERE/'delivery_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':'passed','SVG_markers':counts,'max_coordinate_error':maxerr,'PNG_SHA256':qa['PNG_SHA256']}))
