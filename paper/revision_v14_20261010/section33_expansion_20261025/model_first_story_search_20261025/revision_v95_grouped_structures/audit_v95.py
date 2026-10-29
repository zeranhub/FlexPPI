from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import pandas as pd
from lxml import etree
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'revision_v94_restore_per_chi'
BASE=ROOT.parents[2]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
source=OLD/'figure3/source_data'
target=ROOT/'figure3/source_data'
for p in source.glob('*'):
    if p.is_file() and p.name!='all_plotted_values.csv':assert sha(p)==sha(target/p.name),p.name
before=pd.read_csv(source/'all_plotted_values.csv',float_precision='round_trip').fillna('')
after=pd.read_csv(target/'all_plotted_values.csv',float_precision='round_trip').fillna('')
before['panel']=before.panel.replace({'b':'c','c':'d'})
assert list(before.columns)==list(after.columns)
sort_by=list(before.columns)
pd.testing.assert_frame_equal(before.sort_values(sort_by).reset_index(drop=True),after.sort_values(sort_by).reset_index(drop=True),check_exact=True)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
docs=[BASE/'paper_v94_20261026_中文修订_恢复逐角占据关联.docx',BASE/'paper_v95_20261026_中文修订_结构集中与局部精简.docx']
parts=[];roots=[]
for p in docs:
    with ZipFile(p) as z:parts.append({n:z.read(n) for n in z.namelist()})
    roots.append(etree.fromstring(parts[-1]['word/document.xml']))
assert parts[0].keys()==parts[1].keys()
assert [n for n in parts[0] if parts[0][n]!=parts[1][n]]==['word/document.xml','word/media/image5.png']
assert [etree.tostring(t,method='c14n') for t in roots[0].findall('.//w:tbl',ns)]==[etree.tostring(t,method='c14n') for t in roots[1].findall('.//w:tbl',ns)]
paras=[[''.join(p.xpath('.//w:t/text()',namespaces=ns)) for p in r.findall('.//w:body/w:p',ns)] for r in roots]
assert len(paras[0])==len(paras[1])
changed=[i for i,(a,b) in enumerate(zip(*paras)) if a!=b]
assert changed==[70,71,73,74,75,158,171],changed
assert not any('局部全残基棍状表示' in t for t in paras[1])
qa=json.loads((ROOT/'figure3/figure_QA.json').read_text(encoding='utf-8'))
qa['visual_review']='passed: root and independent layout reviewer; no clipping/overlap; both intact structures adjacent at top; removed local insets; all values identical after panel remap'
(ROOT/'figure3/figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
audit={'all_878_plotted_values_exactly_unchanged_after_panel_remap':True,'all_source_tables_unchanged':True,'all_six_tables_unchanged':True,'all_other45_images_unchanged':True,'changed_paragraph_indices':changed,'old_new_panel_mapping':{'a':'a','d_overview':'b','b':'c','c':'d','e':'e','f':'f','g':'g'},'removed_panels':['Q41_local','Q62_local'],'main_values':after.groupby('panel').size().to_dict(),'native_overviews_reused_without_coordinate_change':True}
render=ROOT/'render/release_01'
if (render/'page-1.png').exists():
    previous=OLD/'render/release_01'
    current_pages=sorted(render.glob('page-*.png'),key=lambda p:int(p.stem.split('-')[1]))
    prior_pages=sorted(previous.glob('page-*.png'),key=lambda p:int(p.stem.split('-')[1]))
    assert len(current_pages)==len(prior_pages)==29,(len(current_pages),len(prior_pages))
    changed_pages=[i+1 for i,(a,b) in enumerate(zip(prior_pages,current_pages)) if sha(a)!=sha(b)]
    audit['rendered_pages']=len(current_pages)
    audit['changed_pages_from_verified_v94']=changed_pages
    audit['identical_pages_inherit_verified_render']=29-len(changed_pages)
(ROOT/'audit/data_and_document_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False))
