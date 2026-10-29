from pathlib import Path
import hashlib,json,zipfile,posixpath
import numpy as np
import pandas as pd
from docx import Document
from lxml import etree
from PIL import Image,ImageChops
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';OLD=ROOT.parent/'figure7_clean_axis_v32_20261014'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
STEM='fig7_combined_predictions_and_MD_curves'
DOC=BASE/'paper_v33_20261014_中文修订_图7预测合并与自由能曲线_三线表.docx'
SOURCE=BASE/'paper_v32_20261014_中文修订_图7横轴标注整理_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json');assert sha(DOC)==integration['document_sha256']
assert sha(SOURCE)==integration['source_sha256']=='301daed0997726951d5b6a958532606faffd5a74480b29170ceaabca39d95d35'
d=Document(DOC);old=Document(SOURCE);text='\n'.join(p.text for p in d.paragraphs)
assert len(d.tables)==11 and len(d.paragraphs)==len(old.paragraphs)
assert sum(p.text!=q.text for p,q in zip(d.paragraphs,old.paragraphs))==3
assert '图7a右' in text and '图7b–d' in text
assert not any(term in text for term in ['图7c–e','图7b）','同类型逐角占据','完整逐角占据的对应','展开各自χ1–χ3的完整角度占据'])
result=next(p.text for p in d.paragraphs if p.text.startswith('全部六个谷氨酰胺按冻结Hχ预测从低到高排列'))
assert 'Q49、Q2的χ3占据较分散' not in result
assert all(value in text for value in ['0.080 nats','0.105 nats','Q40高估0.368','0.400','0.510','35.41','43.17','7.76'])
cells=lambda doc:[[[c.text for c in r.cells] for r in t.rows] for t in doc.tables]
assert cells(d)==cells(old)
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
assert [r['anchor_x_data'] for r in labels]==[2,31,40,41,49,62]
assert len(set(r['y_axes'] for r in labels))==1
assert fa['occupancy_heatmaps_removed'] and fa['six_Gln_predictions_moved_to_a_right']
newpoints=pd.read_csv(FIG/(STEM+'.points.csv'));oldpoints=pd.read_csv(OLD/'figures/fig7_clean_axis_six_Gln_and_occupancy.points.csv')
assert len(newpoints)==576 and newpoints.groupby('panel').size().to_dict()=={'a_Gln':12,'a_Hchi':132,'b':144,'c':144,'d':144}
assert 'b_probability' not in set(newpoints.panel)
oldpoints=oldpoints[oldpoints.panel!='b_probability'].copy()
oldpoints['panel']=oldpoints.panel.map({'a_Hchi':'a_Hchi','b_Hchi':'a_Gln','c':'b','d':'c','e':'d'})
paired=newpoints.merge(oldpoints,on=['panel','residue_id','kind','chi','x'],suffixes=('_new','_old'),validate='one_to_one')
assert len(paired)==576 and np.allclose(paired.y_new,paired.y_old,rtol=0,atol=0,equal_nan=True)
assert np.allclose(paired.error_new,paired.error_old,rtol=0,atol=0,equal_nan=True)
assert paired.low_Fstar_le1_new.fillna('absent').equals(paired.low_Fstar_le1_old.fillna('absent'))
assert pd.read_csv(FIG/(STEM+'.display_widths.csv')).equals(pd.read_csv(OLD/'figures/fig7_clean_axis_six_Gln_and_occupancy.display_widths.csv'))
assert newpoints.query('panel=="a_Hchi" and kind=="MD_reference"').y.eq(0).sum()==3
pdf=PdfReader(TMP/'v33_word.pdf');assert len(pdf.pages)==28
assert all(len((page.extract_text() or '').strip())>100 and '\ufffd' not in (page.extract_text() or '') for page in pdf.pages)
comparisons=[]
for i in range(1,29):
    path=TMP/'render_v33'/f'page-{i}.png';prev=TMP/'render_v32'/f'page-{i}.png'
    exact=ImageChops.difference(Image.open(path),Image.open(prev)).getbbox() is None
    comparisons.append({'page':i,'pixel_identical_to_reviewed_v32':exact,'PNG_sha256':sha(path)})
write(QA/'page_pixel_comparison.json',comparisons)
write(QA/'numerical_content_audit.json',{'status':'passed','document_sha256':sha(DOC),'all576_retained_points_exact':True,'all66_predictions_and6Gln_retained':True,'all432_original_MD_curve_points_exact':True,'occupancy_heatmaps_and_crossrefs_removed':True,'six_Gln_anchor_positions':[2,31,40,41,49,62],'native_tables_math_comments_preserved':integration['all_tables_math_comments_preserved'],'other7_active_images_identical':True,'changed_text_paragraphs':3,'pages':28})
print(json.dumps({'status':'passed','pages':28,'changed_pages':[p['page'] for p in comparisons if not p['pixel_identical_to_reviewed_v32']]},ensure_ascii=False))
