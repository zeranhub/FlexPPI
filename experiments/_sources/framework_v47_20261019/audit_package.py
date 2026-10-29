"""Independently compare every native part and all non-editable body nodes."""
from pathlib import Path
import hashlib,json,zipfile,re
from lxml import etree
ROOT=Path(__file__).resolve().parent
QA=ROOT/'manuscript_QA'
integration=json.loads((QA/'integration_audit.json').read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
source=Path(integration['source']);new=Path(integration['document'])
assert sha(source)==integration['source_sha256'] and sha(new)==integration['document_sha256']
def parts(p):
    with zipfile.ZipFile(p) as z:
        assert z.testzip() is None
        return {n:z.read(n) for n in z.namelist()}
a,b=parts(source),parts(new)
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
ar,br=etree.fromstring(a['word/document.xml']),etree.fromstring(b['word/document.xml'])
ab,bb=ar.find('w:body',NS),br.find('w:body',NS)
c14=lambda e:etree.tostring(e,method='c14n')
text=lambda e:''.join(e.xpath('.//w:t/text()',namespaces=NS))
def caption_text(p):
    result=[]
    for child in p:
        if child.tag=='{'+NS['m']+'}oMath':
            token=''.join(child.xpath('.//m:t/text()',namespaces=NS))
            assert token in ('Ri','ti')
            result.append(token[0]+'ᵢ')
        else:result.append(text(child))
    return ''.join(result)
allowed=set(integration['editable_body_indices'])
checks={
 'same_package_parts':set(a)==set(b),
 'only_document_and_figure1_media_changed':{n for n in a if a[n]!=b[n]}=={'word/document.xml','word/media/image1.png'},
 'same_body_node_count':len(ab)==len(bb),
 'all_noneditable_body_nodes_exact':len(ab)==len(bb) and all(c14(x)==c14(y) for i,(x,y) in enumerate(zip(ab,bb)) if i not in allowed),
 'figure1_media_is_final_PNG':b['word/media/image1.png']==Path(integration['figure1']).read_bytes(),
 'all_other_media_exact':all(a[n]==b[n] for n in a if n.startswith('word/media/') and n!='word/media/image1.png'),
 'all_native_tables_exact':list(map(c14,ar.xpath('.//w:tbl',namespaces=NS)))==list(map(c14,br.xpath('.//w:tbl',namespaces=NS))) and len(br.xpath('.//w:tbl',namespaces=NS))==10,
 'all_native_equations_exact':list(map(c14,ar.xpath('.//m:oMath',namespaces=NS)))==list(map(c14,br.xpath('.//m:oMath',namespaces=NS))) and len(br.xpath('.//m:oMath',namespaces=NS))==36,
 'all_comment_parts_exact':all(a[n]==b[n] for n in a if 'comment' in n.lower()),
 'all_comment_anchors_exact':all(list(map(c14,ar.xpath('.//w:'+k,namespaces=NS)))==list(map(c14,br.xpath('.//w:'+k,namespaces=NS))) for k in ['commentRangeStart','commentRangeEnd','commentReference']),
 'all_relationships_exact':all(a[n]==b[n] for n in a if n.endswith('.rels')),
 'all_headers_footers_styles_sections_exact':all(a[n]==b[n] for n in a if any(s in n for s in ['word/header','word/footer','word/styles','word/numbering','word/settings'])) and list(map(c14,ar.xpath('.//w:sectPr',namespaces=NS)))==list(map(c14,br.xpath('.//w:sectPr',namespaces=NS))),
 'all_reference_paragraphs_exact':[text(p) for p in ab if re.match(r'^\[\d+\]',text(p))]==[text(p) for p in bb if re.match(r'^\[\d+\]',text(p))] and len([p for p in bb if re.match(r'^\[\d+\]',text(p))])==28,
 'same_active_images_order':ar.xpath('.//a:blip/@r:embed',namespaces=NS)==br.xpath('.//a:blip/@r:embed',namespaces=NS) and len(br.xpath('.//a:blip',namespaces=NS))==6,
 'caption_matches_final_asset_file':caption_text(bb[integration['caption_body_index']])==(ROOT/'figure1/caption_zh.md').read_text(encoding='utf8').strip(),
 'all_readout_clarifications_exact':all(text(ab[item['body_index']])==item['before'] and text(bb[item['body_index']])==item['after'] for item in integration['prose_changes']),
 'no_new_training_MD_inference':all(integration[k]==0 for k in ['new_MD','new_training','new_inference'])
}
report={'status':'passed' if all(checks.values()) else 'failed','source':str(source),'source_sha256':sha(source),'document':str(new),'document_sha256':sha(new),'figure1_sha256':sha(integration['figure1']),'checks':checks}
(QA/'package_preservation_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(report,ensure_ascii=False))
assert all(checks.values())
