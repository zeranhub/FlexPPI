from pathlib import Path
from zipfile import ZipFile
from lxml import etree
import pandas as pd
import json,hashlib,sys
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent;OLD=ROOT.parent/'revision_v95_grouped_structures';BASE=ROOT.parents[2]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for p in (OLD/'figure3/source_data').glob('*'):
    if p.is_file() and p.name!='all_plotted_values.csv':assert sha(p)==sha(ROOT/'figure3/source_data'/p.name),p.name
before=pd.read_csv(OLD/'figure3/source_data/all_plotted_values.csv',float_precision='round_trip').fillna('')
after=pd.read_csv(ROOT/'figure3/source_data/all_plotted_values.csv',float_precision='round_trip').fillna('')
contact=before[before.kind.eq('CaM_contact_residues_4p5A')];assert len(contact)==25
contact.to_csv(ROOT/'figure3/source_data/contact_counts_retained_outside_main.csv',index=False)
before=before[~before.kind.eq('CaM_contact_residues_4p5A')].copy()
before['panel']=before.panel.replace({'c':'b','d':'c','e':'d','f':'e','g':'f'})
assert len(before)==len(after)==853
sort=list(before.columns)
pd.testing.assert_frame_equal(before.sort_values(sort).reset_index(drop=True),after.sort_values(sort).reset_index(drop=True),check_exact=True)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
paths=[BASE/'paper_v95_20261026_中文修订_结构集中与局部精简.docx',BASE/'paper_v96_20261026_中文修订_图3视觉层级重排.docx']
parts=[];trees=[]
for path in paths:
    with ZipFile(path) as z:parts.append({n:z.read(n) for n in z.namelist()})
    trees.append(etree.fromstring(parts[-1]['word/document.xml']))
assert [n for n in parts[0] if parts[0][n]!=parts[1][n]]==['word/document.xml','word/media/image5.png']
assert [etree.tostring(t,method='c14n') for t in trees[0].findall('.//w:tbl',ns)]==[etree.tostring(t,method='c14n') for t in trees[1].findall('.//w:tbl',ns)]
paras=[[''.join(p.xpath('.//w:t/text()',namespaces=ns)) for p in tree.findall('.//w:body/w:p',ns)] for tree in trees]
assert len(paras[0])==len(paras[1])
changed=[i for i,(a,b) in enumerate(zip(*paras)) if a!=b]
assert changed==[70,71,73,74,75,76,79,158,171,185],changed
assert not any('图3g' in p for p in paras[1])
audit=dict(all853_model_MD_and_cohort_values_exactly_unchanged=True,all25_static_contact_counts_retained_outside_main=True,no_model_reference_points_removed=True,all_six_tables_unchanged=True,all_other45_media_unchanged=True,changed_paragraphs=changed,main_plot_counts=after.groupby('panel').size().to_dict(),source_tables_unchanged=True)
render=ROOT/'render/release_02'
if (render/'page-1.png').exists():
    pages=sorted(render.glob('page-*.png'),key=lambda p:int(p.stem.split('-')[1]))
    old_pages=sorted((OLD/'render/release_01').glob('page-*.png'),key=lambda p:int(p.stem.split('-')[1]))
    assert len(pages)==len(old_pages)==29
    audit['rendered_pages']=29
    audit['changed_pages_from_verified_v95']=[i+1 for i,(a,b) in enumerate(zip(old_pages,pages)) if sha(a)!=sha(b)]
(ROOT/'audit/data_and_document_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False))
