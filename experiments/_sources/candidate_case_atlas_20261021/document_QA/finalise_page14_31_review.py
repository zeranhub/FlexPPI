"""Bind the manually re-opened corrected page and byte-identical pages to final DOCX."""
from pathlib import Path
from datetime import datetime, timezone
import json, hashlib

QA=Path(__file__).parent
old_path=QA/'page14_31_visual_QA.json'
report=json.loads(old_path.read_text(encoding='utf8'))
binding_path=QA/'render_binding.json'
binding=json.loads(binding_path.read_text(encoding='utf8'))
expected='e479ea9c92e80161216cf418935883245f86f1bfa780044bca0abacce3cc9744'
assert binding['document_sha256']==expected
assert hashlib.sha256(Path(report['document_path']).read_bytes()).hexdigest()==expected
previous_doc=report['document_sha256']
previous_binding=report['render_binding_SHA256']
assert previous_doc=='3285733b085526f8b2942cec1a4043b8a656ffbfc33967b570dba207eaaae2b8'
render=Path(binding['render_directory'])
changed=[]
for row in report['pages']:
    n=row['page'];path=render/f'page-{n}.png'
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest==binding['PNG_SHA256'][str(n)]
    unchanged=digest==row['PNG_SHA256']
    if not unchanged: changed.append(n)
    row['previous_PNG_SHA256']=row['PNG_SHA256']
    row['PNG_SHA256']=digest
    row['path']=str(path)
    row['status']='passed'
    row['reused_by_identical_PNG_bytes']=unchanged
    row['review_method_final']=('Prior individually inspected page is byte-for-byte identical, SHA-256 independently checked.' if unchanged else 'Final rendered page individually reopened with view_image(detail=original), corrected label and complete layout read and verified.')
    row['required_changes']=[]
    if n==18:
        row['notes']='Corrected Hchi scope table reads 同类型完整肽段, 8位. All characters, native three-line table, body paragraphs and pagination remain intact.'
assert changed==[18],changed
report['status']='passed'
report['document_sha256']=expected
report['PDF_SHA256']=binding['PDF_SHA256']
report['render_binding_SHA256']=hashlib.sha256(binding_path.read_bytes()).hexdigest()
report['render_directory']=str(render)
report['reviewed_at_utc']=datetime.now(timezone.utc).isoformat()
report['required_changes']=[]
report['method']='18 pages were individually inspected on render_1. In final render_2, page 18 was individually reopened and reviewed at original resolution; the other 17 PNGs were independently verified byte-identical by SHA-256 and their previous visual conclusions reused.'
report['revision_history']=[{'previous_document_sha256':previous_doc,'previous_render_binding_SHA256':previous_binding,'resolved_issue':'C08_Hchi_scope_label','resolution':'同类型完整界面 changed to 同类型完整肽段','reopened_pages':[18],'byte_identical_pages_reused':[n for n in range(14,32) if n!=18]}]
old_path.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'status':report['status'],'document_sha256':expected,'reopened_pages':changed,'byte_identical_pages_reused':17},ensure_ascii=False))
