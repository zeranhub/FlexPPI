from pathlib import Path
import hashlib,json,re,zipfile,posixpath
import numpy as np
import pandas as pd
from docx import Document
from lxml import etree
from PIL import Image,ImageChops
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';OLD=ROOT.parent/'rich_angle_figure_v30_20261014'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
STEM='fig7_Hchi_profile_six_Gln_axis_and_occupancy'
DOC=BASE/'paper_v31_20261014_中文修订_图7简化与六谷氨酰胺轴标注_三线表.docx'
SOURCE=BASE/'paper_v30_20261014_中文修订_全残基逐角占据与自由能验证_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json');assert sha(DOC)==integration['document_sha256']
assert sha(SOURCE)==integration['source_sha256']=='bc85e75a3d6945e123e98910404de27c08d2d30c14bd114694934119ae696d47'
d=Document(DOC);old=Document(SOURCE);text='\n'.join(p.text for p in d.paragraphs)
assert '全部151个可定义χ角' not in text
assert '并在横轴标出全部六个谷氨酰胺位点' in text
caption=next(p.text for p in d.paragraphs if p.text.startswith('图7：'))
assert 'χ1–χ4' not in caption and '预测减参考的误差' not in caption and '灰格' not in caption
assert all(s in text for s in ['图7a','图7b','图7c–e','0.080 nats','0.979','0.868','0.105 nats','Q40高估0.368','35.41','43.17','7.76','0.400','0.510'])
refs=lambda doc:[p.text for p in doc.paragraphs if re.match(r'^\[\d+\]',p.text)]
assert refs(d)==refs(old) and len(refs(d))==27
cells=lambda doc:[[[c.text for c in r.cells] for r in t.rows] for t in doc.tables]
assert cells(d)==cells(old) and len(d.tables)==11
assert len(d.paragraphs)==len(old.paragraphs)
assert sum(p.text!=q.text for p,q in zip(d.paragraphs,old.paragraphs))==2
ns={'a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def active_images(path):
    with zipfile.ZipFile(path) as z:
        rels=etree.fromstring(z.read('word/_rels/document.xml.rels'));mapping={e.get('Id'):e.get('Target') for e in rels}
        doc=etree.fromstring(z.read('word/document.xml'))
        return [hashlib.sha256(z.read(posixpath.normpath('word/'+mapping[e.get('{'+ns['r']+'}embed')]))).hexdigest() for e in doc.xpath('.//a:blip',namespaces=ns)]
before,after=active_images(SOURCE),active_images(DOC)
assert len(before)==len(after)==8 and [i for i,(a,b) in enumerate(zip(before,after)) if a!=b]==[6]
assert after[6]==sha(FIG/(STEM+'.png'))
fa=readj(FIG/(STEM+'.audit.json'))
expected=[2,31,40,41,49,62]
assert [r['residue_id'] for r in fa['all6_Gln_axis_labels']]==expected
assert [r['x_data'] for r in fa['all6_Gln_axis_labels']]==expected
assert all(r['label']=='Q'+str(r['residue_id']) for r in fa['all6_Gln_axis_labels'])
assert fa['all6_Gln_axis_labels'][3]['y_axes']<fa['all6_Gln_axis_labels'][2]['y_axes']
z=np.load(OLD/'experiment/ubiquitin_all66_chi_501frames.npz',allow_pickle=False)
assert set(z['author_residue_id'][z['amino_acid']=='Q'].astype(int))==set(expected)
points=pd.read_csv(FIG/(STEM+'.points.csv'))
previous=pd.read_csv(OLD/'figures/fig7_full_residue_and_chi_occupancy.points.csv')
assert len(points)==1872 and set(points.panel)=={'a_Hchi','b_Hchi','b_probability','c','d','e'}
assert points.groupby('panel').size().to_dict()=={'a_Hchi':132,'b_Hchi':12,'b_probability':1296,'c':144,'d':144,'e':144}
previous=previous[~previous.panel.isin(['a_error','a_width'])]
paired=points.merge(previous,on=['panel','residue_id','kind','chi','x'],suffixes=('_new','_old'),validate='one_to_one')
assert len(paired)==1872 and np.allclose(paired.y_new,paired.y_old,rtol=0,atol=0,equal_nan=True)
for col in ['count','error']:
    assert np.allclose(paired[col+'_new'],paired[col+'_old'],rtol=0,atol=0,equal_nan=True)
assert paired.low_Fstar_le1_new.fillna('absent').equals(paired.low_Fstar_le1_old.fillna('absent'))
assert points.query('panel=="a_Hchi" and kind=="MD_reference"').y.eq(0).sum()==3
assert pd.read_csv(FIG/(STEM+'.display_widths.csv')).equals(pd.read_csv(OLD/'figures/fig7_full_residue_and_chi_occupancy.display_widths.csv'))
pdf=PdfReader(TMP/'v31_word.pdf');assert len(pdf.pages)==28
assert all(len((page.extract_text() or '').strip())>100 and '\ufffd' not in (page.extract_text() or '') for page in pdf.pages)
comparisons=[]
for i in range(1,29):
    path=TMP/'render_v31'/f'page-{i}.png';prev=TMP/'render_v30'/f'page-{i}.png'
    exact=ImageChops.difference(Image.open(path),Image.open(prev)).getbbox() is None
    comparisons.append({'page':i,'pixel_identical_to_reviewed_v30':exact,'PNG_sha256':sha(path)})
write(QA/'page_pixel_comparison.json',comparisons)
write(QA/'numerical_content_audit.json',{'status':'passed','document_sha256':sha(DOC),'all1872_retained_figure_points_exact':True,'all132_a_Hchi_values_exact':True,'all1740_b_to_e_values_exact':True,'a_error_and_width_tracks_removed':True,'all6_Gln_axis_labels_exact_positions':expected,'Q40_Q41_vertical_stagger_only':True,'all3_valid_zero_values_retained':True,'native_tables_math_comments_preserved':integration['all_tables_math_comments_preserved'],'other7_active_images_identical':True,'references_27_identical':True,'changed_text_paragraphs':2,'pages':28})
print(json.dumps({'status':'passed','pages':28,'changed_pages':[p['page'] for p in comparisons if not p['pixel_identical_to_reviewed_v30']],'unchanged_pages':[p['page'] for p in comparisons if p['pixel_identical_to_reviewed_v30']]},ensure_ascii=False))
