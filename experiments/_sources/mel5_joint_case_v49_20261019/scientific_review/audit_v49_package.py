"""Independent read-only native Word package and approved-edit audit."""
from pathlib import Path
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent.parent
QA = ROOT / 'manuscript_QA'
integration_path = QA / 'integration_audit.json'
integration = json.loads(integration_path.read_text(encoding='utf8'))
OLD, DOC, FIG = [Path(integration[k]) for k in ('source','document','figure4')]
NS = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
      'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
text = lambda element: ''.join(n.text or '' for n in element.findall('.//w:t', NS))
xml = lambda element: ET.tostring(element, encoding='utf8')
checks = []
def check(name, passed, detail):
    checks.append({'check':name, 'passed':bool(passed), 'detail':detail})

with zipfile.ZipFile(OLD) as z:
    oldparts = {n:z.read(n) for n in z.namelist()}
with zipfile.ZipFile(DOC) as z:
    newparts = {n:z.read(n) for n in z.namelist()}
oldroot,newroot = [ET.fromstring(p['word/document.xml']) for p in (oldparts,newparts)]
oldbody,newbody = [r.find('w:body',NS) for r in (oldroot,newroot)]
mapping = {int(k):int(v) for k,v in integration['source_to_final_body_index'].items()}
editable = set(integration['editable_source_body_indices'])
inserted = integration['inserted_body_nodes']
check('source_and_target_hashes_match_integration',sha(OLD)==integration['source_sha256'] and sha(DOC)==integration['document_sha256'],
      {'source':sha(OLD),'document':sha(DOC)})
check('unchanged_zip_inventory',oldparts.keys()==newparts.keys(),
      {'added':sorted(newparts.keys()-oldparts.keys()),'removed':sorted(oldparts.keys()-newparts.keys())})
changed_parts=sorted(k for k in oldparts.keys()&newparts.keys() if oldparts[k]!=newparts[k])
check('only_document_and_figure4_parts_changed',changed_parts==['word/document.xml','word/media/image41.png'],changed_parts)
check('complete_injective_source_body_mapping',set(mapping)==set(range(len(oldbody))) and len(set(mapping.values()))==len(oldbody) and len(newbody)==len(oldbody)+1,
      {'source_nodes':len(oldbody),'new_nodes':len(newbody),'mapping_nodes':len(mapping),'inserted_n':len(inserted)})
unexpected=[i for i in range(len(oldbody)) if i not in editable and xml(oldbody[i])!=xml(newbody[mapping[i]])]
check('all_unrelated_body_nodes_preserved',not unexpected,{'changed_protected_nodes':unexpected})
actual_changes=[i for i in range(len(oldbody)) if xml(oldbody[i])!=xml(newbody[mapping[i]])]
check('all_source_node_changes_accounted_for',set(actual_changes)==editable,
      {'changed_source_indices':actual_changes,'approved_source_indices':sorted(editable)})
prose_checks=[]
for change in integration['prose_changes']:
    oldi,newi=change['source_body_index'],change['body_index']
    prose_checks.append({'source_body_index':oldi,'final_body_index':newi,
        'mapping_exact':mapping[oldi]==newi,'before_exact':text(oldbody[oldi])==change['before'],
        'after_exact':text(newbody[newi])==change['after']})
check('all_prose_replacements_exact',all(all(c[k] for k in ('mapping_exact','before_exact','after_exact')) for c in prose_checks),prose_checks)
new_indices=set(range(len(newbody)))-set(mapping.values())
check('single_additional_mechanism_paragraph_exact',len(inserted)==1 and new_indices=={inserted[0]['body_index']} and text(newbody[inserted[0]['body_index']])==inserted[0]['text'],inserted)
for label,path,count in [('native_three_line_tables','w:tbl',10),('native_equations','m:oMath',36)]:
    a,b=oldroot.findall('.//'+path,NS),newroot.findall('.//'+path,NS)
    check(label+'_preserved',len(a)==len(b)==count and all(xml(x)==xml(y) for x,y in zip(a,b)),
          {'source':len(a),'target':len(b),'exact_xml':all(xml(x)==xml(y) for x,y in zip(a,b))})
refs=lambda r:[text(p) for p in r.findall('.//w:p',NS) if re.match(r'^\[\d+\]',text(p))]
oldrefs,newrefs=refs(oldroot),refs(newroot)
check('reference_1_to_27_preserved_28_intentionally_replaced',len(oldrefs)==len(newrefs)==28 and oldrefs[:27]==newrefs[:27] and oldrefs[27]!=newrefs[27] and 'Madura' in newrefs[27] and 'M113.464560' in newrefs[27],
      {'source_references':len(oldrefs),'target_references':len(newrefs),'new_reference28':newrefs[-1]})
comment_tags=('commentRangeStart','commentRangeEnd','commentReference')
comment_counts={tag:[len(r.findall('.//w:'+tag,NS)) for r in (oldroot,newroot)] for tag in comment_tags}
check('comments_and_anchors_preserved',all(a==b for a,b in comment_counts.values()) and all(oldparts[k]==newparts.get(k) for k in oldparts if 'comment' in k.lower()),comment_counts)
images=lambda r:[x.get('{'+NS['r']+'}embed') for x in r.findall('.//a:blip',NS)]
check('all6_active_images_and_order_preserved',images(oldroot)==images(newroot) and len(images(newroot))==6,images(newroot))
check('figure4_embedding_exact',newparts['word/media/image41.png']==FIG.read_bytes() and sha(FIG)==integration['figure4_sha256'],
      {'source_asset_sha256':sha(FIG),'embedded_sha256':hashlib.sha256(newparts['word/media/image41.png']).hexdigest()})
prose=ROOT/'prose_v49.json'
check('approved_prose_contract_current',sha(prose)==integration['prose_contract_sha256'],{'prose_sha256':sha(prose)})
joined='\n'.join(text(p) for p in newroot.findall('.//w:p',NS))
check('former_KRAS_case_removed',not any(t in joined for t in ('KRAS','7OW5','7OW6','JDIa41b1','Poole')),
      {'checked_tokens':['KRAS','7OW5','7OW6','JDIa41b1','Poole']})
report={'status':'passed' if all(c['passed'] for c in checks) else 'failed',
        'reviewer':'independent residue_pair_evidence agent','document':str(DOC),'document_sha256':sha(DOC),
        'source':str(OLD),'source_sha256':sha(OLD),'figure4':str(FIG),'figure4_sha256':sha(FIG),
        'integration_audit_sha256':sha(integration_path),'checks':checks}
(QA/'package_preservation_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
(ROOT/'scientific_review/package_preservation_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'status':report['status'],'failed_checks':[c for c in checks if not c['passed']]},ensure_ascii=False,indent=2))
