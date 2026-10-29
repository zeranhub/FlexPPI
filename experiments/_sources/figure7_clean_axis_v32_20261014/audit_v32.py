from pathlib import Path
import hashlib,json,zipfile,posixpath
import pandas as pd
from docx import Document
from lxml import etree
from PIL import Image,ImageChops
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';OLD=ROOT.parent/'figure7_axis_sites_v31_20261014'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
STEM='fig7_clean_axis_six_Gln_and_occupancy'
DOC=BASE/'paper_v32_20261014_中文修订_图7横轴标注整理_三线表.docx'
SOURCE=BASE/'paper_v31_20261014_中文修订_图7简化与六谷氨酰胺轴标注_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json');assert sha(DOC)==integration['document_sha256']
assert sha(SOURCE)==integration['source_sha256']=='a9011e448b28be019d4adb936c84e2ef45cd8196b0705f615b96817c0efdc2b7'
d=Document(DOC);old=Document(SOURCE)
assert [p.text for p in d.paragraphs]==[p.text for p in old.paragraphs]
assert len(d.tables)==11
ns={'a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def active_images(path):
    with zipfile.ZipFile(path) as z:
        rels=etree.fromstring(z.read('word/_rels/document.xml.rels'));mapping={e.get('Id'):e.get('Target') for e in rels}
        doc=etree.fromstring(z.read('word/document.xml'))
        return [hashlib.sha256(z.read(posixpath.normpath('word/'+mapping[e.get('{'+ns['r']+'}embed')]))).hexdigest() for e in doc.xpath('.//a:blip',namespaces=ns)]
before,after=active_images(SOURCE),active_images(DOC)
assert len(before)==len(after)==8 and [i for i,(a,b) in enumerate(zip(before,after)) if a!=b]==[6]
assert after[6]==sha(FIG/(STEM+'.png'))
fa=readj(FIG/(STEM+'.audit.json'));labels=fa['all6_Gln_axis_labels']
assert [r['residue_id'] for r in labels]==[2,31,40,41,49,62]
assert [r['anchor_x_data'] for r in labels]==[2,31,40,41,49,62]
assert len(set(r['y_axes'] for r in labels))==1
assert [r['label_x_data'] for r in labels]==[2,31,38,43,49,62]
assert fa['numeric_axis_ticks']==[1,10,20,30,40,50,60,72] and fa['axis_labels_no_overlap']
newpoints=pd.read_csv(FIG/(STEM+'.points.csv'));oldpoints=pd.read_csv(OLD/'figures/fig7_Hchi_profile_six_Gln_axis_and_occupancy.points.csv')
assert len(newpoints)==1872 and newpoints.equals(oldpoints)
assert pd.read_csv(FIG/(STEM+'.display_widths.csv')).equals(pd.read_csv(OLD/'figures/fig7_Hchi_profile_six_Gln_axis_and_occupancy.display_widths.csv'))
pdf=PdfReader(TMP/'v32_word.pdf');assert len(pdf.pages)==28
assert all(len((page.extract_text() or '').strip())>100 and '\ufffd' not in (page.extract_text() or '') for page in pdf.pages)
comparisons=[]
for i in range(1,29):
    path=TMP/'render_v32'/f'page-{i}.png';prev=TMP/'render_v31'/f'page-{i}.png'
    exact=ImageChops.difference(Image.open(path),Image.open(prev)).getbbox() is None
    comparisons.append({'page':i,'pixel_identical_to_reviewed_v31':exact,'PNG_sha256':sha(path)})
assert [p['page'] for p in comparisons if not p['pixel_identical_to_reviewed_v31']]==[17]
write(QA/'page_pixel_comparison.json',comparisons)
write(QA/'numerical_content_audit.json',{'status':'passed','document_sha256':sha(DOC),'all1872_figure_points_exact':True,'all_text_unchanged':True,'six_Gln_anchor_positions':[2,31,40,41,49,62],'six_Gln_labels_one_row':True,'Q40_Q41_label_offsets_with_correct_anchors':True,'regular_numeric_ticks_one_row_including40':True,'native_tables_math_comments_preserved':integration['all_tables_math_comments_preserved'],'other7_active_images_identical':True,'pages':28})
print(json.dumps({'status':'passed','pages':28,'changed_pages':[17],'unchanged_pages':27},ensure_ascii=False))
