from pathlib import Path
import hashlib,json,zipfile,posixpath,sys
RT=Path(r'E:\Project\FlexPPI\experiments\_sources\completed_results_20261012\plot_runtime')
sys.path.insert(0,str(RT))
import numpy as np
import pandas as pd
from docx import Document
from lxml import etree
from PIL import Image,ImageChops
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures'
OLD=ROOT.parent/'figure7_combined_prediction_v33_20261014'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
STEM='fig7_sequence_order_predictions_and_MD_curves';OLDSTEM='fig7_combined_predictions_and_MD_curves'
DOC=BASE/'paper_v34_20261015_中文修订_图7六谷氨酰胺顺序对照_三线表.docx'
SOURCE=BASE/'paper_v33_20261014_中文修订_图7预测合并与自由能曲线_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json')
assert sha(DOC)==integration['document_sha256']
assert sha(SOURCE)==integration['source_sha256']=='750a323bee603d05f448c1deaaf863ad3a3970fac8aef780480837627873d4e7'
d=Document(DOC);old=Document(SOURCE)
assert len(d.tables)==11 and len(d.paragraphs)==len(old.paragraphs)
changed_paragraphs=[i for i,(p,q) in enumerate(zip(d.paragraphs,old.paragraphs)) if p.text!=q.text]
assert changed_paragraphs==[51,100,106],changed_paragraphs
text='\n'.join(p.text for p in d.paragraphs)
for index in [51,100,106]:
    assert '按冻结预测从低到高排列' not in d.paragraphs[index].text
    assert '带符号误差' not in d.paragraphs[index].text
    assert '逐位误差' not in d.paragraphs[index].text
assert 'Q40高估0.368' not in d.paragraphs[100].text
assert '按残基编号排列' in d.paragraphs[100].text
assert all(value in text for value in ['0.080 nats','0.105 nats','0.400','0.510','35.41','43.17','7.76'])
assert d.paragraphs[98].text==old.paragraphs[98].text
assert d.paragraphs[99].text==old.paragraphs[99].text
assert d.paragraphs[109].text==old.paragraphs[109].text
cells=lambda doc:[[[c.text for c in r.cells] for r in t.rows] for t in doc.tables]
assert cells(d)==cells(old)
ns={'a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def active_images(path):
    with zipfile.ZipFile(path) as z:
        rels=etree.fromstring(z.read('word/_rels/document.xml.rels'));mapping={e.get('Id'):e.get('Target') for e in rels}
        xml=etree.fromstring(z.read('word/document.xml'))
        return [hashlib.sha256(z.read(posixpath.normpath('word/'+mapping[e.get('{'+ns['r']+'}embed')]))).hexdigest() for e in xml.xpath('.//a:blip',namespaces=ns)]
before,after=active_images(SOURCE),active_images(DOC)
assert len(before)==len(after)==8
assert [i for i,(a,b) in enumerate(zip(before,after)) if a!=b]==[6]
assert after[6]==sha(FIG/(STEM+'.png'))
fa=readj(FIG/(STEM+'.audit.json'))
labels=fa['all6_Gln_axis_labels']
assert [r['anchor_x_data'] for r in labels]==[2,31,40,41,49,62]
assert len(set(r['y_axes'] for r in labels))==1
newpoints=pd.read_csv(FIG/(STEM+'.points.csv'));oldpoints=pd.read_csv(OLD/'figures'/(OLDSTEM+'.points.csv'))
assert len(newpoints)==len(oldpoints)==576
assert newpoints.groupby('panel').size().to_dict()=={'a_Gln':12,'a_Hchi':132,'b':144,'c':144,'d':144}
key=['panel','residue_id','kind','chi','x']
paired=newpoints.merge(oldpoints,on=key,suffixes=('_new','_old'),validate='one_to_one')
assert len(paired)==576
unchanged=paired[paired.panel!='a_Gln']
assert len(unchanged)==564 and np.allclose(unchanged.y_new,unchanged.y_old,rtol=0,atol=0,equal_nan=True)
assert np.allclose(paired.error_new,paired.error_old,rtol=0,atol=0,equal_nan=True)
assert paired.low_Fstar_le1_new.fillna('absent').equals(paired.low_Fstar_le1_old.fillna('absent'))
expected={site:index for index,site in enumerate([2,31,40,41,49,62])}
assert all(row.y==expected[row.residue_id] for row in newpoints[newpoints.panel=='a_Gln'].itertuples())
assert pd.read_csv(FIG/(STEM+'.display_widths.csv')).equals(pd.read_csv(OLD/'figures'/(OLDSTEM+'.display_widths.csv')))
gln=pd.read_csv(FIG/(STEM+'.all6_sequence_order.csv'))
assert gln.author_residue_id.tolist()==[2,31,40,41,49,62]
oldgln=pd.read_csv(OLD/'figures'/(OLDSTEM+'.all6_prediction_order.csv'))
pd.testing.assert_frame_equal(gln[oldgln.columns].sort_values('author_residue_id').reset_index(drop=True),oldgln.sort_values('author_residue_id').reset_index(drop=True),check_exact=True)
assert np.allclose(gln.signed_error_nats,gln.predicted_Hchi_nats-gln.frozen_MD_joint_Hchi_nats,rtol=0,atol=1e-15)
assert newpoints.query('panel=="a_Hchi" and kind=="MD_reference"').y.eq(0).sum()==3
q40=newpoints.query('panel=="a_Gln" and residue_id==40')
assert len(q40)==2 and np.allclose(q40.error,.3676173686981201,rtol=0,atol=1e-15)
assert abs((gln.predicted_Hchi_nats-gln.frozen_MD_joint_Hchi_nats).abs().mean()-.1050252914428711)<1e-6
pdf=PdfReader(TMP/'v34_word.pdf');assert len(pdf.pages)==28
assert all(len((page.extract_text() or '').strip())>100 and '\ufffd' not in (page.extract_text() or '') for page in pdf.pages)
comparisons=[]
for i in range(1,29):
    path=TMP/'render_v34'/f'page-{i}.png';prev=TMP/'render_v33'/f'page-{i}.png'
    exact=ImageChops.difference(Image.open(path),Image.open(prev)).getbbox() is None
    comparisons.append({'page':i,'pixel_identical_to_reviewed_v33':exact,'PNG_sha256':sha(path)})
write(QA/'page_pixel_comparison.json',comparisons)
write(QA/'numerical_content_audit.json',{'status':'passed','document_sha256':sha(DOC),'all576_values_retained':True,'all564_otherpanel_points_exact':True,'all6_paired_Hchi_and_signed_error_data_unchanged':True,'a_Gln_row_order_only_changed':True,'signed_error_numbers_removed_from_display_only':True,'Q40_retained_in_all_datasets_and_summary':True,'all432_original_MD_curve_points_exact':True,'all6_sequence_order':[2,31,40,41,49,62],'native_tables_math_comments_preserved':integration['all_tables_math_comments_preserved'],'other7_active_images_identical':True,'changed_text_paragraphs':changed_paragraphs,'pages':28})
print(json.dumps({'status':'passed','pages':28,'changed_pages':[p['page'] for p in comparisons if not p['pixel_identical_to_reviewed_v33']]},ensure_ascii=False))
