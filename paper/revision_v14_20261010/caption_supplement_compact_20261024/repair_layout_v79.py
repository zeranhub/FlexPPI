"""Allow complete S1 panels to paginate and restore a short plotting definition."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib
import json
import zipfile
from lxml import etree as E

ROOT = Path(__file__).resolve().parent
TARGET = ROOT.parent/'paper_v79_20261024_中文修订_图2图注与补充表精简.docx'
NS = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'm':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
W = '{'+NS['w']+'}'

def text(el):
    return ''.join(el.xpath('.//w:t/text()',namespaces=NS))

def paragraph(template,value):
    el=deepcopy(template)
    properties=el.find('w:r/w:rPr',NS)
    for child in list(el):
        if child.tag!=W+'pPr': el.remove(child)
    run=E.SubElement(el,W+'r')
    if properties is not None: run.append(deepcopy(properties))
    node=E.SubElement(run,W+'t'); node.text=value
    node.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    return el

assert hashlib.sha256(TARGET.read_bytes()).hexdigest()=='79bff4f9667da3855d9b027a913992cfb8340c1ef5dcf8720f6d7aa3245fa33d'
(ROOT/'v79_before_layout_repair.docx').write_bytes(TARGET.read_bytes())
with zipfile.ZipFile(TARGET) as package:
    entries=package.infolist(); parts={entry.filename:package.read(entry.filename) for entry in entries}
root=E.fromstring(parts['word/document.xml']); body=root.find('w:body',NS)
math_before=Counter(E.tostring(el) for el in root.xpath('//m:oMath',namespaces=NS))
tables=body.findall('w:tbl',NS)
row_values=[[[text(c) for c in r.findall('w:tc',NS)] for r in t.findall('w:tr',NS)] for t in tables]
s1_rows=tables[0].findall('w:tr',NS)
assert len(s1_rows)==24
for index in [7,13]:
    for p in s1_rows[index].xpath('.//w:p',namespaces=NS):
        p.find('w:pPr/w:keepNext',NS).set(W+'val','0')
graph=next(p for p in body if text(p).startswith('图2a密度包含3,063条记录的'))
value=text(graph).replace('图2b各方法均含3,063条记录，','图2b各方法均含3,063条记录，组内从左到右为Full、Intra、Static、AA，小提琴表示分布；',1)
assert value!=text(graph)
body.replace(graph,paragraph(graph,value))
window=next(p for p in body if text(p).startswith('固定预测在原汇总标签与第一条共同窗口标签上的核查'))
old='固定预测在原汇总标签与第一条共同窗口标签上的核查保存在配套数据common_window_label_validation.csv中，完整保留各来源的MAE和相关。'
new='原标签与共同窗口的数值核查见配套数据window_label_check.csv，保留各来源的MAE和相关。'
assert text(window).startswith(old)
body.replace(window,paragraph(window,text(window).replace(old,new,1)))
(ROOT/'source_data'/'window_label_check.csv').write_bytes((ROOT/'source_data'/'common_window_label_validation.csv').read_bytes())
assert Counter(E.tostring(el) for el in root.xpath('//m:oMath',namespaces=NS))==math_before
assert [[[text(c) for c in r.findall('w:tc',NS)] for r in t.findall('w:tr',NS)] for t in tables]==row_values
xml=E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
with zipfile.ZipFile(TARGET,'w',zipfile.ZIP_DEFLATED) as package:
    for entry in entries: package.writestr(entry,xml if entry.filename=='word/document.xml' else parts[entry.filename])
audit_path=ROOT/'v79_patch_audit.json'
audit=json.loads(audit_path.read_text(encoding='utf-8'))
audit['target_sha256']=hashlib.sha256(TARGET.read_bytes()).hexdigest()
audit['layout_repair']={'S1_allowed_breaks':'Between complete panels A/B/C only; within-panel rows remain together',
                         'restored_graph_definition':'Full/Intra/Static/AA order and violin distribution in supplement C',
                         'shortened_data_alias':'window_label_check.csv; original long filename also preserved'}
audit['visual_qa']='pending rerender'
audit_path.write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'sha256':audit['target_sha256'],'repair':'PASS'},ensure_ascii=False))
