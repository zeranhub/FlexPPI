"""Read-only final DOCX content, native table, figure-byte and reference audit."""
from pathlib import Path
import zipfile,json,hashlib,re,posixpath,sys
ROOT=Path(__file__).resolve().parent.parent
BASE=ROOT.parent.parent
DOC=BASE/'FlexPPI_32个潜在案例图册_背景证据与局限_20261021.docx'
MAP=ROOT/'ATLAS_FIGURE_MAP.json'
MAIN=BASE/'paper_v49_20261019_中文修订_MEL5双动态机制实例_三线表.docx'
EXPECTED_MAIN_SHA='5f79ddf7d657dc9b782c9fb432e373a1e7c170c1c88d019e528fbf73c366c222'
sys.path.insert(0,r'E:\Project\FlexPPI\experiments\_sources\completed_results_20261012\plot_runtime')
from lxml import etree as E
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
    'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
    'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
    'pr':'http://schemas.openxmlformats.org/package/2006/relationships'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(b):return hashlib.sha256(b).hexdigest()
def text(node):return ''.join(node.xpath('.//w:t/text()',namespaces=NS))
def norm(s):return re.sub(r'\s+','',s or '')
def val(node,kind):
    k=node.find('w:'+kind,NS) if node is not None else None
    return k.get('{'+NS['w']+'}val') if k is not None else None

def main():
    assert DOC.exists(),f'Document has not been created: {DOC}'
    bundle=json.loads(MAP.read_text(encoding='utf8'))
    narrative=json.loads((ROOT/'evidence/narratives.json').read_text(encoding='utf8'))
    catalog=json.loads((ROOT/'catalog/candidate_inventory.json').read_text(encoding='utf8'))
    bycat={x['case_id']:x for x in catalog['cases']}
    narr={x['catalog_case_id']:x for x in narrative['cases']}
    expected={f'C{i:02}' for i in range(1,33)}
    assert set(narr)==set(bundle['figures'])==expected
    before=sha(MAIN);assert before==EXPECTED_MAIN_SHA==bundle['main_manuscript_sha256']
    docsha=sha(DOC)
    with zipfile.ZipFile(DOC) as z:
        assert z.testzip() is None
        root=E.fromstring(z.read('word/document.xml'))
        relroot=E.fromstring(z.read('word/_rels/document.xml.rels'))
        rels={x.get('Id'):x.get('Target') for x in relroot}
        body=root.find('w:body',NS);children=list(body)
        alltext=text(root)
        heading_positions={}
        for i,node in enumerate(children):
            style=node.find('w:pPr/w:pStyle',NS)
            if style is not None and style.get('{'+NS['w']+'}val')=='Heading1':
                hit=re.match(r'^(C\d{2})\s',text(node))
                if hit:
                    assert hit[1] not in heading_positions,('Duplicate case heading',hit[1])
                    heading_positions[hit[1]]=i
        assert set(heading_positions)==expected,('Case headings',sorted(heading_positions))
        image_records=[]
        for inline in root.xpath('.//wp:inline',namespaces=NS):
            pr=inline.find('wp:docPr',NS);blip=inline.find('.//a:blip',NS)
            assert pr is not None and blip is not None
            hit=re.match(r'^(C\d{2})\s',pr.get('name',''))
            assert hit,('Unidentified image',pr.get('name'))
            cid=hit[1];rid=blip.get('{'+NS['r']+'}embed')
            target=posixpath.normpath(posixpath.join('word',rels[rid]))
            h=digest(z.read(target));expectedh=bundle['figures'][cid]['sha256']['png']
            assert h==expectedh==sha(bundle['figures'][cid]['png']),(cid,'Embedded image bytes differ from final map')
            image_records.append(dict(case_id=cid,rId=rid,package_part=target,SHA256=h,map_source=bundle['figures'][cid]['png']))
        assert len(image_records)==32 and {r['case_id'] for r in image_records}==expected
        assert len(root.xpath('.//a:blip',namespaces=NS))==32
        assert len({x['SHA256'] for x in image_records})==32,'Every case must have its own figure'
        case_records=[]
        hsorted=sorted(heading_positions.items(),key=lambda x:x[1])
        next_major=[i for i,n in enumerate(children) if n.find('w:pPr/w:pStyle',NS) is not None and n.find('w:pPr/w:pStyle',NS).get('{'+NS['w']+'}val')=='Heading1']
        for cid,start in hsorted:
            end=min((j for j in next_major if j>start),default=len(children))
            nodes=children[start:end];case_text=''.join(text(n) for n in nodes)
            source=narr[cid]
            missing=[j for j,p in enumerate(source['paragraphs']) if norm(p) not in norm(case_text)]
            assert not missing,(cid,'Narrative paragraphs missing',missing)
            roles=source.get('paragraph_roles',[])
            assert any(re.search(r'局限|限制|不能|缺陷|未解决',r) for r in roles),(cid,'No explicit limitation role')
            limitpara=source['paragraphs'][-1]
            assert len(limitpara)>50 and re.search(r'不能|不足|未|局限|限制|不等|尚|只|弱|问题|缺',limitpara),(cid,'Empty/vague limitation paragraph')
            pdbs=bycat.get(cid,{}).get('pdbs',source['pdbs'])
            assert all(p.upper() in case_text.upper() for p in pdbs),(cid,'PDB missing',pdbs)
            expected_refs=bundle['case_references'][cid]
            cite='['+','.join(str(n) for n in expected_refs)+']'
            assert cite in case_text,(cid,'Case source citations missing',cite)
            assert '图'+cid in case_text,(cid,'Figure caption absent')
            tables=[n for n in nodes if n.tag=='{'+NS['w']+'}tbl']
            assert len(tables)==1,(cid,'Case scope table')
            st=text(tables[0])
            assert 'RMSF' in st or 'Hχ' in st or 'Hχ' in case_text
            if cid=='C08':assert '同类型完整肽段' in st and '同类型完整界面' not in st,(cid,'Peptide Hchi scope mislabelled as interface',st)
            if cid=='C14':assert '16/19位' in st,(cid,'DMF5/DMF4 interface counts switched or omitted',st)
            if cid=='C17':assert '每态10条' in st,(cid,'CatchBond replicate protocol mislabelled')
            case_records.append(dict(case_id=cid,title=text(children[start]),narrative_paragraphs=len(source['paragraphs']),
                 all_paragraphs_present=True,has_explicit_limitations=True,all_structure_PDBs_present=True,
                 source_reference_numbers=expected_refs,scope_table_text=st,
                 has_one_case_figure=True,figure_SHA256=next(x['SHA256'] for x in image_records if x['case_id']==cid)))
        tables=root.xpath('.//w:tbl',namespaces=NS);assert len(tables)==34,('Native table count',len(tables))
        table_records=[]
        for i,t in enumerate(tables,1):
            bd=t.find('w:tblPr/w:tblBorders',NS)
            assert val(bd,'top')=='single' and val(bd,'bottom')=='single',(i,'Top/bottom rules absent')
            assert all(val(bd,k) in ['nil','none'] for k in ['left','right','insideH','insideV']),(i,'Non-three-line borders')
            rs=t.findall('w:tr',NS)
            assert rs[0].find('w:trPr/w:tblHeader',NS) is not None,(i,'Header not repeated')
            for c in rs[0].findall('w:tc',NS):assert val(c.find('w:tcPr/w:tcBorders',NS),'bottom')=='single',(i,'Header separator absent')
            for row in rs[1:]:
                for c in row.findall('w:tc',NS):
                    for bb in c.findall('w:tcPr/w:tcBorders/*',NS):assert bb.get('{'+NS['w']+'}val') in ['nil','none'],(i,'Extra body cell rules')
            table_records.append(dict(table=i,rows=len(rs),native_three_line=True,repeated_header=True))
        # Complete references and every case's assigned citations are in the package.
        refpos=next(i for i,n in enumerate(children) if text(n)=='参考文献与结构条目')
        last=next(i for i,n in enumerate(children) if i>refpos and text(n)=='发现性结构配对索引')
        reftext=''.join(text(n) for n in children[refpos:last])
        reference_records=[]
        assigned={n for ns in bundle['case_references'].values() for n in ns}
        for ref in bundle['references']:
            assert ref['number'] in assigned
            assert f"[{ref['number']}]" in reftext and norm(ref['title']) in norm(reftext),(ref['number'],'Incomplete reference')
            for field in ['doi','url']:
                if ref.get(field):assert ref[field] in reftext,(ref['number'],'Missing reference identifier',field)
            reference_records.append(dict(number=ref['number'],title=ref['title'],doi=ref.get('doi'),url=ref.get('url'),assigned_to_case=True,full_entry_present=True))
        # Discovery index table preserves every previous pair; it is separate from 32 biological narratives.
        index=tables[-1];pairs=json.loads((ROOT/'catalog/all_discovery_pair_index.json').read_text(encoding='utf8'))
        cells=[text(c) for c in index.findall('.//w:tc',NS)]
        for p in pairs:assert p['a'].upper()+' / '+p['b'].upper() in cells
        bookmarks=root.xpath('.//w:bookmarkStart/@w:name',namespaces=NS)
        assert all(bookmarks.count(cid)==1 for cid in expected)
        assert not root.xpath('.//w:ins|.//w:del',namespaces=NS),'Unexpected unresolved tracked changes'
    after=sha(MAIN);assert before==after==EXPECTED_MAIN_SHA
    assert sha(DOC)==docsha,'DOCX changed during read-only audit'
    report=dict(status='passed',document=str(DOC),document_sha256=docsha,
        figures=32,cases=32,source_narrative_paragraphs=sum(c['narrative_paragraphs'] for c in case_records),
        native_three_line_tables=34,references=len(reference_records),discovery_pairs=len(pairs),
        full_source_narratives_present=True,all_cases_have_explicit_limitations=True,all_32_embedded_images_byte_exact=True,
        reference_entries_and_case_citations_complete=True,main_manuscript=str(MAIN),main_manuscript_sha256=after,
        main_manuscript_unchanged=True,main_manuscript_baseline_sha256=EXPECTED_MAIN_SHA,
        final_figure_map_sha256=sha(MAP),narrative_source_sha256=sha(ROOT/'evidence/narratives.json'),
        case_reviews=case_records,image_reviews=image_records,table_reviews=table_records,reference_reviews=reference_records,
        visual_layout_review='Separate render/page reviews performed by root and other agents; this report audits OOXML, content and exact source bytes only',
        new_training=0,new_inference=0,new_MD=0,document_modified_by_auditor=False)
    out=ROOT/'document_QA/independent_document_audit.json';out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({k:report[k] for k in ['status','document_sha256','figures','cases','native_three_line_tables','references','discovery_pairs','main_manuscript_unchanged']},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
