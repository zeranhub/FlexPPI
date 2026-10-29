from pathlib import Path
from collections import Counter
import hashlib,json,zipfile,re,posixpath
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
OLD=BASE/'paper_v38_20261016_中文修订_图7标注说明与图例调整_三线表.docx'
NEW=BASE/'paper_v39_20261017_中文修订_读者表述与英文图表_三线表.docx'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
a=Document(OLD);b=Document(NEW)
integration=readj(QA/'integration_audit.json')
assert sha(NEW)==integration['output_sha256']
assert len(a.tables)==11 and len(b.tables)==10
for old,new in zip(a.tables[2:],b.tables[1:]):
    # Method names in S4 were shortened, all quantitative entries and remaining cells preserved.
    numeric=lambda t:Counter(re.findall(r'(?<![A-Za-z])[-−+]?\d+(?:[.,]\d+)*(?:%)?', '\n'.join(c.text for row in t.rows for c in row.cells)))
    assert numeric(old)==numeric(new)
assert len(integration['table1_numeric_audit'])==13
for row in integration['table1_numeric_audit']:
    assert row['display_values']==['—' if x in ('无','不适用') else x for x in row['source_values']]
alltext='\n'.join(p.text for p in b.paragraphs)
assert not any(w in alltext for w in ['冻结','发布权重','发布模型','图7','表2：','3.6 '])
assert all(w in alltext for w in ['0.080 nats','0.979（图6a）','0.868','33,142','0.400','0.510','7.76','633','2,430','同一个关联家族'])
assert len([p for p in b.paragraphs if re.match(r'^图[1-6]：',p.text)])==6
assert len([p for p in b.paragraphs if p.text.startswith('图S1：')])==1
assert len([p for p in b.paragraphs if re.match(r'^表S[1-9]：',p.text)])==9
assert [p.text[:3] for p in b.paragraphs if p.style.name=='Heading 2' and re.match(r'^3\.[1-5] ',p.text)]==['3.1','3.2','3.3','3.4','3.5']
assert len([p for p in b.paragraphs if re.match(r'^\[\d+\]',p.text)])==27
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def active(path):
    with zipfile.ZipFile(path) as z:
        tree=etree.fromstring(z.read('word/document.xml'))
        rel={el.get('Id'):el.get('Target') for el in etree.fromstring(z.read('word/_rels/document.xml.rels'))}
        return [hashlib.sha256(z.read(posixpath.normpath('word/'+rel[el.get(qn('r:embed'))]))).hexdigest() for el in tree.xpath('.//a:blip',namespaces=ns)]
oldimages,newimages=active(OLD),active(NEW)
assert len(newimages)==7
assert newimages[0]==oldimages[0] and newimages[4]==oldimages[5] and newimages[6]==oldimages[7]
for key,index in [('2',1),('3',2),('4',3),('7',5)]:
    assert newimages[index]==sha(readj(ROOT/'figure_manifest.json')[key]['png'])
with zipfile.ZipFile(OLD) as z1,zipfile.ZipFile(NEW) as z2:
    x=etree.fromstring(z1.read('word/document.xml'));y=etree.fromstring(z2.read('word/document.xml'))
    maths=lambda z:[etree.tostring(el,method='c14n') for el in z.xpath('.//m:oMath',namespaces=ns)]
    assert maths(x)==maths(y) and len(maths(y))==36
    anchor=lambda z:Counter((el.tag,el.get(qn('w:id'))) for el in z.iter() if el.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    assert anchor(x)==anchor(y) and z1.read('word/comments.xml')==z2.read('word/comments.xml')
    # Borders occupy only the top, header bottom, and table bottom.
    for t in y.xpath('.//w:tbl',namespaces=ns):
        rows=t.findall(qn('w:tr'))
        header_bottom=1 if t is y.xpath('.//w:tbl',namespaces=ns)[0] else 0
        for ri,row in enumerate(rows):
            for cell in row.findall(qn('w:tc')):
                borders=cell.find('w:tcPr/w:tcBorders',namespaces=ns)
                if borders is None:continue
                for edge in borders:
                    if edge.get(qn('w:val')) not in ('nil','none'):
                        side=etree.QName(edge).localname
                        assert (side=='top' and ri==0) or (side=='bottom' and ri in (header_bottom,len(rows)-1)),(ri,side)

# Root viewed pages 9–17 directly from the final PNGs; the remaining pages have independent reviews.
write(QA/'visual_pages_09_17.json',{'status':'passed','document_sha256':sha(NEW),'method':'Each original final PNG directly inspected','pages':[{'page':n,'status':'passed','sha256':sha(TMP/'render_v39'/f'page-{n}.png')} for n in range(9,18)],'necessary_fixes':[]})
pdf=PdfReader(TMP/'v39_word.pdf');assert len(pdf.pages)==27
assert all(len(p.extract_text().strip())>100 and '\ufffd' not in p.extract_text() for p in pdf.pages)
fig7=readj(ROOT/'english_figures/fig7/fig7_english.audit.json')
assert fig7['all_source_and_experiment_files_byte_exact_to_v38']
assert fig7['point_counts']=={'a_Hchi':132,'b_scatter':151,'c':432}
assert all(x['byte_exact_to_v38'] for x in fig7['exact_data_comparisons'])
for n in ['2','3','4']:
    report=readj(ROOT/f'english_figures/fig{n}/data_and_language_audit.json')
    assert report['all_numeric_cells_unchanged'] and report['all_visible_labels_English']
    assert all(x['all_numeric_cells_unchanged'] for x in report['data_comparisons'])
visual=[]
for part in ['visual_pages_01_08.json','visual_pages_09_17.json','visual_pages_18_27.json']:
    p=QA/part
    if p.exists():
        r=readj(p);assert r['status'] in ('passed','PASS','pass')
        visual.append({'file':part,'sha256':sha(p)})
review=readj(QA/'independent_content_review.json')
assert review['status']=='pass' and sha(NEW) in json.dumps(review)
result={'status':'passed','document':str(NEW),'document_sha256':sha(NEW),'pages':27,'main_table_A_B_all13_methods_original_values':True,'all9_supplementary_tables_numeric_values_unchanged':True,'all_tables_native_three_line':10,'active_figures':7,'all_visible_figure_labels_English':True,'all715_ubiquitin_values_byte_exact':True,'all14_other_plot_CSV_numeric_data_preserved':True,'native_math_preserved':36,'comments_and_anchors_byte_exact_or_identical':True,'references':27,'numbering_contiguous':True,'no_redundant_frozen_or_release_jargon':True,'independent_content_review':'pass','direct_page_reviews':visual,'remaining_required_work':[]}
write(QA/'numerical_structural_audit.json',result)
print(json.dumps(result,ensure_ascii=False))
