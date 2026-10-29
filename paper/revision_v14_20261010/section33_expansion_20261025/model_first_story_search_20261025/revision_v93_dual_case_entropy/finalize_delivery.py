"""Finalize an already built/reviewed DOCX and a compact reproducible figure bundle."""
from pathlib import Path
import hashlib,json,zipfile,csv
from lxml import etree

ROOT=Path(__file__).resolve().parent;BASE=ROOT.parents[2]
DOC=BASE/'paper_v93_20261026_中文修订_双案例与构象占据整合.docx'
SRC=BASE/'paper_v90_20261026_中文修订_图3数据呈现与机制推断.docx'
EXPECTED_DOC='6c603bfff03d033610d228854c05c0cb65bce4ceee3e9f9dae2bf8b3433aab44'
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
assert ROOT.is_relative_to(BASE) and ROOT.name=='revision_v93_dual_case_entropy'
assert sha(DOC)==EXPECTED_DOC
assert sha(SRC)=='dce6f4cd2a23cb30c94b9e710b5a520bf3b6385afeb6b51748e50a4c53d30dd5'
document_audit=json.loads((ROOT/'manuscript/v93_final_document_audit.json').read_text(encoding='utf-8'))
assert document_audit['docx_sha256']==EXPECTED_DOC and document_audit['page_count']==29
assert document_audit['visual_review']['status']=='passed'
old=ROOT/'render/release_04';new=ROOT/'render/release_05'
changed=[i for i in range(1,30) if sha(old/f'page-{i}.png')!=sha(new/f'page-{i}.png')]
assert changed==[1,2,12,13,29]
assert document_audit['figure3_sha256']=='a3b8c18fb3d52e0f5a1143ccd7af4a418940151af1f131a4233f70aae0cfc0bc'

# Final text integrity: captions, sequential sections, actual group-difference wording.
with zipfile.ZipFile(DOC) as z:
    assert z.testzip() is None
    root=etree.fromstring(z.read('word/document.xml'))
    texts=root.xpath('//w:body/w:p//w:t/text()',namespaces={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'})
text=''.join(texts)
for stale in ['高值组实际角区平均宽7.76','移至图S4c、d，不将','非本轮新家族划分的类别分析','图5：','图5d']:
    assert stale not in text,stale
assert '其实际角区较低值组平均宽7.76°/角' in text
assert '表明模型能够恢复同类型侧链在不同环境中的状态多样性差异' in text
assert document_audit['source_table_xml_identical'] and document_audit['unchanged_media_count']==44

contract=json.loads((ROOT/'manuscript/final_text_contract.json').read_text(encoding='utf-8'))
md=['# 新版图3图注与相关正文','',contract['section33']['title'],'']
md.extend(contract['section33']['before_figure']);md.extend(['',contract['section33']['caption'],'',*contract['section33']['after_figure'],'',contract['section34']['title'],'',*contract['section34']['paragraphs']])
(ROOT/'图3图注与3.3_3.4正文.md').write_text('\n\n'.join(md),encoding='utf-8')

audit=dict(status='complete',docx_sha256=EXPECTED_DOC,source_unchanged=True,Word_pages=29,latest_render='release_05',latest_changed_pages_human_reviewed=changed,unchanged_other24_pages_verified_by_SHA=True,figure3_and_caption_same_page=11,TCR_figure4_page=14,S4_and_caption_same_page=29,no_main_figure5=True,sections='3.3 dual cases -> 3.4 angular occupancy -> 3.5 TCR',all908_main_plot_values_audited=True,all24_bootstrap_summaries_audited=True,all_six_tables_and44_other_media_unchanged=True,no_new_model_training_or_MD_simulation=True)
(ROOT/'audit/final_release_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')

exclude={'render','source_media','__pycache__','.git'}
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and not exclude.intersection(p.relative_to(ROOT).parts) and p.suffix not in {'.tmp','.zip','.pdf-lock'} and p.name!='archive_contents_sha256.csv')
assert all(p.resolve().is_relative_to(ROOT.resolve()) for p in files)
manifest=ROOT/'audit/archive_contents_sha256.csv'
with manifest.open('w',encoding='utf-8-sig',newline='') as out:
    writer=csv.writer(out);writer.writerow(['relative_path','bytes','sha256'])
    for p in files:writer.writerow([p.relative_to(ROOT).as_posix(),p.stat().st_size,sha(p)])
files.append(manifest)
archive=BASE/'paper_v93_20261026_图3与S4源数据_PyMOL与审计.zip'
assert not archive.exists(),'Refuse to replace an existing delivery archive'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,Path(ROOT.name)/p.relative_to(ROOT))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None and len(z.infolist())==len(files)
    assert not any('/render/' in x.filename or '/source_media/' in x.filename for x in z.infolist())
delivery=[DOC,*[BASE/f'figure3_v93_20261026_双案例与构象占据.{e}' for e in ['png','pdf','svg']],*[BASE/f'figureS4_v93_20261026_构象占据补充核查.{e}' for e in ['png','pdf','svg']],archive]
result=dict(status='complete',archive_files=len(files),artifacts=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in delivery])
(BASE/'paper_v93_20261026_交付清单.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(status='complete',Word_pages=29,archive_files=len(files),archive_bytes=archive.stat().st_size,docx_sha256=EXPECTED_DOC,archive_sha256=sha(archive)),ensure_ascii=False))
