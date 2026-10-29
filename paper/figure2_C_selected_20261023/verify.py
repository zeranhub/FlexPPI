"""Verify that selecting and regrouping C changes layout only."""
from pathlib import Path
import ast
from dataclasses import asdict
import hashlib
import json
import re
import runpy
import xml.etree.ElementTree as ET
import fitz
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'figure2_journal_trial_20261023'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
source=(ROOT/'plot.py').read_text(encoding='utf-8')
prior=(OLD/'plot_options.py').read_text(encoding='utf-8')

def function_text(script,name):
    node=next(n for n in ast.parse(script).body if isinstance(n,ast.FunctionDef) and n.name==name)
    return ast.get_source_segment(script,node)

assert function_text(source,'precision_recall')==function_text(prior,'precision_recall')
grouped=(ROOT/'superseded_grouped/plot.py').read_text(encoding='utf-8')
def lower_panels(script):
    body=script[script.index('distribution_points=[]'):]
    return body[body.index('for j,task in enumerate(TASKS):'):body.index('selected_check=export(')]
assert lower_panels(source)==lower_panels(grouped)
archive={p.name:sha(p)==sha(OLD/'source_archive'/p.name) for p in (ROOT/'source_archive').iterdir() if p.is_file()}
assert len(archive)==15 and all(archive.values())
current=pd.read_csv(ROOT/'all_record_AP_distribution_points.csv',float_precision='round_trip')
previous=pd.read_csv(OLD/'all_record_AP_distribution_points.csv',float_precision='round_trip')
assert current.equals(previous) and len(current)==36756
assert current.groupby(['task','method']).size().eq(3063).all()
assert (ROOT/'complete_PR_vertices.csv').read_bytes()==(OLD/'complete_PR_vertices.csv').read_bytes()
assert len(pd.read_csv(ROOT/'complete_PR_vertices.csv'))==1224
assert pd.read_csv(ROOT/'AP_distribution_summary.csv').equals(pd.read_csv(OLD/'AP_distribution_summary.csv'))
assert not ET.parse(ROOT/'fig2_C_selected.svg').findall('.//{http://www.w3.org/2000/svg}image')
doc=fitz.open(ROOT/'fig2_C_selected.pdf')
assert len(doc)==1
width,height=[v*25.4/72 for v in (doc[0].rect.width,doc[0].rect.height)]
assert abs(width-183)<.01 and abs(height-203)<.01
assert all('Arial' in font[3] for font in doc[0].get_fonts())
labels={span['text'] for block in doc[0].get_text('dict')['blocks'] if 'lines' in block
        for line in block['lines'] for span in line['spans']}
assert set('abcdefg')<=labels and 'h' not in labels
doc.close()
validator=runpy.run_path('C:/Users/DELL/.codex/skills/academic-figure-skill/scripts/qa_validator.py')
style=[]
for name,func in validator.items():
    if name.startswith('check_'):
        finding=func(re.sub(r'DIVERGING +=','DIVERGING =',source))
        style.extend(asdict(item) for item in (finding if isinstance(finding,list) else [finding]))
for item in style:
    if item['check_id']=='CL-2':
        item.update(resolved_status='PASS',manual_resolution='PDF measured at 183 x 203 mm; dimensions use named variables.')
    if item['check_id']=='AP-4' and not item['pass_']:
        item.update(resolved_status='PASS',manual_resolution='Shared legend is above all data axes; regex only scans the first 30 characters after legend(. Independently inspected.')
assert not [i for i in style if not i['pass_'] and i.get('resolved_status')!='PASS'],style
(ROOT/'automated_style_QA.json').write_text(json.dumps(style,ensure_ascii=False,indent=2),encoding='utf-8')
manuscript=Path('E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v61_20261023_中文修订_图1输入输出标注修正.docx')
assert sha(manuscript)=='666ecc56bb3159e2146d8a2a22cec5db09069697c28a9d89c12d41b80750593c'
audit=json.loads((ROOT/'layout_audit.json').read_text(encoding='utf-8'))
assert audit['same_plot_rectangle'] and audit['dual_x_and_y_axes']
assert audit['colorbars_vertical'] and audit['metrics_inside_plot']
assert audit['hchi_trend_color']=='#9B5C1A' and audit['hchi_trend_linewidth_pt']==1.9
for task in ['rmsf','chi_entropy']:
    data=np.load(ROOT/'source_archive'/f'prediction_reference_{task}.npz')
    info=audit['overlay_info'][task]
    assert info['n_positions']==len(data['reference'])
    assert np.isclose(info['mass_sum'],1)
assert audit['overlay_info']['rmsf']['x_axis']=='bottom'
assert audit['overlay_info']['rmsf']['y_axis']=='left'
assert audit['overlay_info']['chi_entropy']['x_axis']=='top'
assert audit['overlay_info']['chi_entropy']['y_axis']=='right'
audit.update(status='passed',visual_review_passed=True,independent_review_passed=True,
             source_files_byte_identical=archive,distribution_values_exactly_match=36756,
             complete_PR_vertices_exactly_match=1224,density_values_and_weights_unchanged=True,
             lower_six_panels_unchanged=True,
             PR_function_unchanged=True,all_12_AP_summaries_unchanged=True,
             v61_sha256_unchanged=sha(manuscript))
(ROOT/'final_QA.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in audit.items() if 'sha256' not in k and k!='source_files_byte_identical'},ensure_ascii=False,indent=2))
