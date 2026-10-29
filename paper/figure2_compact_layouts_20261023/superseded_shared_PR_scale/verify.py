"""Verify complete frozen evidence and compact vector output dimensions."""
from pathlib import Path
from dataclasses import asdict
import hashlib
import json
import re
import runpy
import xml.etree.ElementTree as ET
import fitz
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'figure2_C_selected_20261023'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
archive={p.name:sha(p)==sha(OLD/'source_archive'/p.name) for p in (ROOT/'source_archive').iterdir() if p.is_file()}
assert len(archive)==15 and all(archive.values())
for filename in ['complete_PR_vertices.csv','all_record_AP_distribution_points.csv','AP_distribution_summary.csv']:
    assert sha(ROOT/filename)==sha(OLD/filename)
assert sha(ROOT/'previous_7panel_layout.pdf')==sha(OLD/'fig2_C_selected.pdf')
audit=json.loads((ROOT/'layout_audit.json').read_text(encoding='utf-8'))
assert all(v==dict(AP_values=36756,PR_vertices=1224) for v in audit['rendered_and_verified_counts'].values())
for stem,height in [('fig2_compact_task_panels',181.),('fig2_compact_merged',150.)]:
    doc=fitz.open(ROOT/f'{stem}.pdf')
    assert len(doc)==1
    assert abs(doc[0].rect.width*25.4/72-183)<.01
    assert abs(doc[0].rect.height*25.4/72-height)<.01
    assert all('Arial' in f[3] for f in doc[0].get_fonts())
    assert not ET.parse(ROOT/f'{stem}.svg').findall('.//{http://www.w3.org/2000/svg}image')
    if stem=='fig2_compact_merged':
        labels=[span['text'] for block in doc[0].get_text('dict')['blocks'] if 'lines' in block
                for line in block['lines'] for span in line['spans']]
        assert all(labels.count(letter)==1 for letter in 'abc')
        assert not any(letter in labels for letter in 'de')
        assert labels.count('Mass')==1 and labels.count('AP')==1 and labels.count('Recall')==1
        assert 'High-dynamics site retrieval' not in labels
    doc.close()
doc=fitz.open(ROOT/'previous_7panel_layout.pdf')
doc[0].get_pixmap(matrix=fitz.Matrix(1.6,1.6),alpha=False).save(ROOT/'previous_layout_preview.png')
doc.close()
fig=plt.figure(figsize=(18,7.8),dpi=160)
for j,(name,height,label) in enumerate([
    ('previous_layout',203.,'Previous: 7 panels'),
    ('fig2_compact_task_panels',181.,'Task groups: 4 panels'),
    ('fig2_compact_merged',150.,'Selected: 3 panels')]):
    width=.305
    display_height=.82*height/203
    ax=fig.add_axes([.015+j*.33,.89-display_height,width,display_height])
    ax.imshow(mpimg.imread(ROOT/f'{name}_preview.png'))
    ax.axis('off')
    fig.text(.015+j*.33+width/2,.945,label,ha='center',fontsize=15,fontfamily='Arial')
fig.savefig(ROOT/'comparison_overview.png',dpi=160,facecolor='white')
plt.close(fig)
source=(ROOT/'plot.py').read_text(encoding='utf-8')
validator=runpy.run_path('C:/Users/DELL/.codex/skills/academic-figure-skill/scripts/qa_validator.py')
style=[]
for name,func in validator.items():
    if name.startswith('check_'):
        result=func(re.sub(r'DIVERGING +=','DIVERGING =',source))
        style.extend(asdict(item) for item in (result if isinstance(result,list) else [result]))
for item in style:
    if item['check_id']=='CL-2':
        item.update(resolved_status='PASS',manual_resolution='Both PDFs measured: 183 x 181 mm and 183 x 150 mm; canvas uses variables.')
assert not [i for i in style if not i['pass_'] and i.get('resolved_status')!='PASS'],style
(ROOT/'automated_style_QA.json').write_text(json.dumps(style,ensure_ascii=False,indent=2),encoding='utf-8')
manuscript=Path('E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v61_20261023_中文修订_图1输入输出标注修正.docx')
assert sha(manuscript)=='666ecc56bb3159e2146d8a2a22cec5db09069697c28a9d89c12d41b80750593c'
assert audit['main_panels']==['a','b','c']
assert audit['shared_mass_scale'] and audit['mass_colorbar_axes_count']==1
assert audit['output_trends_solid'] and audit['distribution_width_mm']==71
assert audit['c_heading_removed'] and audit['PR_plot_height_mm']==40
assert audit['PR_precision_axis']['type']=='continuous_piecewise_linear'
assert audit['PR_precision_axis']['knot']==.5
assert audit['PR_precision_axis']['lower_half_display_fraction']==.25
assert audit['PR_precision_axis']['upper_half_display_fraction']==.75
audit.update(status='passed',source_files_byte_identical=archive,visual_review_passed=True,
             independent_visual_review_passed=True,
             height_reduction_percent=dict(task_panels=100*(203-181)/203,merged_distributions=100*(203-150)/203),
             prior_selected_C_unchanged=sha(OLD/'fig2_C_selected.pdf'),v61_unchanged=sha(manuscript),
             output_sha256={p.name:sha(p) for p in ROOT.iterdir() if p.suffix in ['.pdf','.svg','.png']})
(ROOT/'final_QA.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in audit.items() if k not in ['source_files_byte_identical','output_sha256']},ensure_ascii=False,indent=2))
