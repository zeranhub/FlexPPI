"""Targeted v90 -> v93 revision. Execute only after revision_content.json is ready.

The separate document-skill edit marker must succeed exactly once before executing
this builder. It is deliberately not invoked here, so rerenders/rebuilds cannot
accidentally count a second artifact operation.
"""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from copy import deepcopy
from lxml import etree
from PIL import Image
import json, hashlib, sys

HERE=Path(__file__).resolve().parent
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
SOURCE=BASE/'paper_v90_20261026_中文修订_图3数据呈现与机制推断.docx'
EXPECTED_SOURCE_SHA='dce6f4cd2a23cb30c94b9e710b5a520bf3b6385afeb6b51748e50a4c53d30dd5'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships','pic':'http://schemas.openxmlformats.org/drawingml/2006/picture'}
def q(prefix,name):return '{'+NS[prefix]+'}'+name
def signature(el):return etree.tostring(el,method='c14n')
def text_of(el):return ''.join(el.xpath('.//w:t/text()',namespaces=NS))
def set_text(p,text):
    old=p.findall('w:r',NS)
    rpr=deepcopy(old[0].find('w:rPr',NS)) if old and old[0].find('w:rPr',NS) is not None else None
    for child in list(p):
        if child.tag!=q('w','pPr'):p.remove(child)
    r=etree.SubElement(p,q('w','r'))
    if rpr is not None:r.append(rpr)
    t=etree.SubElement(r,q('w','t'));t.text=text
    if text.startswith(' ') or text.endswith(' '):t.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    return p
def clone_text(template,text):return set_text(deepcopy(template),text)
def set_paragraph_flag(p,name,value):
    ppr=p.find('w:pPr',NS)
    if ppr is None:ppr=etree.SubElement(p,q('w','pPr'));p.insert(0,ppr)
    el=ppr.find('w:'+name,NS)
    if el is None:el=etree.SubElement(ppr,q('w',name))
    el.set(q('w','val'),'1' if value else '0')
def set_image(node,path,parts,rels,width_mm,name,description):
    rids=node.xpath('.//a:blip/@r:embed',namespaces=NS)
    assert len(rids)==1,rids
    target='word/'+rels[rids[0]]
    parts[target]=Path(path).read_bytes()
    im=Image.open(path);cx=round(width_mm*36000);cy=round(cx*im.height/im.width)
    for n in node.xpath('.//wp:extent|.//a:xfrm/a:ext',namespaces=NS):n.set('cx',str(cx));n.set('cy',str(cy))
    for n in node.xpath('.//wp:docPr|.//pic:cNvPr',namespaces=NS):n.set('name',name);n.set('descr',description)
    return {'media_part':target,'relationship_id':rids[0],'source_png':str(path),'width_mm':width_mm,'height_mm':cy/36000,'sha256':hashlib.sha256(parts[target]).hexdigest()}

def main():
    config_path=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'revision_content.json'
    cfg=json.loads(config_path.read_text(encoding='utf-8'))
    output=Path(cfg.get('output_docx',str(BASE/'paper_v93_20261026_中文修订_双案例与构象占据合图.docx')))
    assert SOURCE.resolve()!=output.resolve()
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED_SOURCE_SHA
    with ZipFile(SOURCE) as z:parts={i.filename:z.read(i.filename) for i in z.infolist()}
    original_parts=parts.copy()
    root=etree.fromstring(parts['word/document.xml']);body=root.find('w:body',NS)
    old_children=list(body);pars=body.findall('w:p',NS)
    assert len(pars)==186 and len(old_children)==193
    before=[signature(p) for p in pars]
    original_tables=[signature(x) for x in body.findall('w:tbl',NS)]
    relroot=etree.fromstring(parts['word/_rels/document.xml.rels'])
    rels={n.get('Id'):n.get('Target') for n in relroot}
    changes={int(i):t for i,t in cfg.get('paragraph_replacements',{}).items()}
    forbidden=set(range(68,91))
    assert not (forbidden & set(changes)), 'Use section33/34 for replaced section text; TCR stays intact except heading.'
    for i,t in changes.items():set_text(pars[i],t)
    set_paragraph_flag(pars[91],'pageBreakBefore',False)
    c33=cfg['section33'];c34=cfg['section34']
    fig3=deepcopy(pars[72]);figS4=deepcopy(pars[89])
    fig_audit=[set_image(fig3,cfg['figure3_png'],parts,rels,cfg.get('figure3_width_mm',164),'Figure 3 Backbone and sidechain case studies with conformation occupancy validation',cfg.get('figure3_alt','CaMKI backbone RMSF, ubiquitin sidechain Hchi, and cohort angle occupancy validation')),
       set_image(figS4,cfg['supplement4_png'],parts,rels,cfg.get('supplement4_width_mm',164),'Figure S4 Supplementary conformation occupancy diagnostics',cfg.get('supplement4_alt','Supplementary paired flexibility, ubiquitin angular entropy and empirical relative free energy diagnostics'))]
    assert fig_audit[0]['media_part']=='word/media/image5.png'
    assert fig_audit[1]['media_part']=='word/media/image40.png'
    pbody=pars[69];phead=pars[68];pcap=pars[73]
    new33=[clone_text(phead,c33['title'])]+[clone_text(pbody,t) for t in c33['before_figure']]
    new33 += [fig3,clone_text(pcap,c33['caption'])]
    after33=[clone_text(pbody,t) for t in c33.get('after_figure',[])]
    if after33:set_paragraph_flag(after33[0],'keepLines',True)
    new33 += after33
    new34=[clone_text(phead,c34['title'])]+[clone_text(pbody,t) for t in c34['paragraphs']]
    tcr=[old_children[i] for i in [77,78,79,82,83,80,81]]
    assert text_of(tcr[0]).startswith('3.4 ')
    set_text(tcr[0],text_of(tcr[0]).replace('3.4 ','3.5 ',1))
    expected_tcr_before=[signature(p) for p in pars[78:84]]
    tail=old_children[91:]
    assert tail[-1].tag==q('w','sectPr')
    s4captions=cfg['supplement4_caption']
    if isinstance(s4captions,str):s4captions=[s.strip() for s in s4captions.split('\n\n') if s.strip()]
    supplement=[figS4]+[clone_text(pars[90],s) for s in s4captions]
    new_children=old_children[:68]+new33+new34+tcr+tail[:-1]+supplement+[tail[-1]]
    body[:]=new_children
    assert original_tables==[signature(x) for x in body.findall('w:tbl',NS)]
    assert expected_tcr_before==[signature(p) for p in pars[78:84]]
    allowed_changes=set(changes)|set(range(68,91))|{91}
    for i,p in enumerate(pars):
        if i not in allowed_changes:assert signature(p)==before[i],f'Unexpected source paragraph modification {i}'
    finalpars=body.findall('w:p',NS)
    full_text='\n'.join(text_of(p) for p in finalpars)
    assert not any(x in full_text for x in ['图5','图 5','Figure 5','Fig. 5']), 'Unresolved main figure 5 reference'
    assert full_text.index(c33['title'])<full_text.index(c34['title'])<full_text.index('3.5 TCR')
    for name in ['图S1','图S2','图S3','图S4']:
        assert name in full_text,name
    parts['word/document.xml']=etree.tostring(root,encoding='UTF-8',xml_declaration=True,standalone=True)
    changed_parts=[k for k in parts if parts[k]!=original_parts[k]]
    assert set(changed_parts)=={'word/document.xml','word/media/image5.png','word/media/image40.png'}
    with ZipFile(output,'w',ZIP_DEFLATED) as z:
        for name,data in parts.items():z.writestr(name,data)
    with ZipFile(output) as z:assert z.testzip() is None
    (HERE/'v93_final_text.txt').write_text(full_text,encoding='utf-8')
    finalinv=[{'paragraph_index':i,'text':text_of(p),'body_child_index':list(body).index(p)} for i,p in enumerate(finalpars)]
    (HERE/'v93_final_paragraphs.json').write_text(json.dumps(finalinv,ensure_ascii=False,indent=2),encoding='utf-8')
    audit={'source':str(SOURCE),'source_sha256':EXPECTED_SOURCE_SHA,'output':str(output),'docx_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
       'config':str(config_path),'changed_zip_parts':changed_parts,'old_paragraphs_replaced_or_moved':list(range(68,91)),
       'additional_source_paragraphs_replaced':sorted(changes),'all_six_tables_unchanged':True,'TCR_body_and_image_caption_unchanged':True,
       'TCR_source_paragraph_order':[77,78,79,82,83,80,81],
       'paragraph_format_changes':{'source_paragraph_91':'pageBreakBefore false, allow Discussion to follow Results without large blank gap','new_first_ubiquitin_paragraph':'keepLines true, avoid two-line fragment below Figure 3'},
       'all_other_source_paragraphs_unchanged':True,'all_other_media_unchanged':True,'figures':fig_audit,
       'new_paragraph_count':len(finalpars),'new_body_child_count':len(body),'main_figure5_references':0}
    (HERE/'v93_build_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(audit,ensure_ascii=False))

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8',errors='replace')
    main()
