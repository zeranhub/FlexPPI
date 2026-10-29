"""Clarify figure 1 inputs and output branches, retaining the v60 palette."""
from pathlib import Path
from copy import deepcopy
import hashlib
import json
import re
import shutil
import subprocess
from lxml import etree as ET

ROOT=Path(__file__).resolve().parents[1]
BASE=Path('E:/Project/FlexPPI/paper/figures_v60/01_模型框架')
OUT=ROOT/'figures_v61'/'01_模型框架'
APP='C:/Program Files/draw.io/draw.io.exe'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    original=(BASE/'fig1_framework.drawio').read_bytes()
    graph=ET.fromstring(original)
    baseline=deepcopy(graph)
    cells={c.get('id'):c for c in graph.findall('.//mxCell')}
    cells['sequence'].set('value','Residue<br>sequence')
    cells['structure-backbone-annotation'].set('value','N–C<sub style="font-size:65%">α</sub>–C;<br>Protein or complex')
    cells['entropy-title'].set('value','Side-chain')
    cells['entropy-title'].find('mxGeometry').set('height','56')
    cells['entropy-metric'].find('mxGeometry').set('y','466')
    # Keep the branch points joined, but leave a clear horizontal run for each arrow.
    cells['head-trunk'].find('mxGeometry/mxPoint[@as="targetPoint"]').set('x','1605')
    for name in ['head-to-rmsf','head-to-hchi']:
        c=cells[name]
        c.set('style',c.get('style').replace('endSize=17;','endSize=12;'))
        geom=c.find('mxGeometry')
        geom.find('mxPoint[@as="sourcePoint"]').set('x','1605')
        geom.find('mxPoint[@as="targetPoint"]').set('x','1643')
        geom.find('Array/mxPoint').set('x','1605')
    allowed={'sequence','structure-backbone-annotation','entropy-title','entropy-metric',
             'head-trunk','head-to-rmsf','head-to-hchi'}
    old={c.get('id'):c for c in baseline.findall('.//mxCell')}
    assert old.keys()==cells.keys()
    for key,c in cells.items():
        if key not in allowed:
            assert ET.tostring(c)==ET.tostring(old[key]),key
    assert len(cells)==208
    output=OUT/'fig1_framework.drawio'
    output.write_bytes(ET.tostring(graph,encoding='UTF-8',xml_declaration=True))
    for ext in ['png','pdf','svg']:
        args=[APP,'--export','--format',ext,'--border','0','--output',str(OUT/('fig1_framework.'+ext)),str(output)]
        if ext=='png':args[1:1]=['--scale','2']
        r=subprocess.run(args,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=90,creationflags=subprocess.CREATE_NO_WINDOW)
        (OUT/('export_'+ext+'.log')).write_text(r.stdout+r.stderr,encoding='utf-8')
        assert r.returncode==0 and (OUT/('fig1_framework.'+ext)).is_file(),ext
    caption=(BASE/'caption_zh.md').read_text(encoding='utf-8')
    caption=caption.replace('氨基酸类型及','残基序列及')
    caption=caption.replace('（c）共享读出层','（c）Hχ为同一残基全部可定义χ角的联合离散状态熵。共享读出层')
    (OUT/'caption_zh.md').write_text(caption,encoding='utf-8')
    audit={'status':'source_passed_visual_pending','baseline':str(BASE),
        'modified_cells':sorted(allowed),'unchanged_cells':len(cells)-len(allowed),
        'palette_images_all_other_labels_and_curves_unchanged':True,
        'labels':{'sequence':'Residue sequence','structure':'Protein or complex','output':'Side-chain entropy Hχ'},
        'branch_geometry':{'junction_x':1605,'arrow_tip_x':1643,'chart_y_axis_x':1660,'horizontal_arrow_length':38,'axis_clearance':17,'arrow_size':12},
        'outputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('fig1_framework.*')}}
    (OUT/'revision_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(OUT)

if __name__=='__main__':main()
