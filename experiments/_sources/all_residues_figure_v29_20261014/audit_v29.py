from pathlib import Path
import hashlib,json,re,zipfile,posixpath
import numpy as np
import pandas as pd
from docx import Document
from lxml import etree
from PIL import Image,ImageChops
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';OLD=ROOT.parent/'curve_only_v28_20261014'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
STEM='fig7_all_residues_Gln_and_MD_free_energy'
DOC=BASE/'paper_v29_20261014_中文修订_全残基与同类型预测及自由能曲线_三线表.docx'
SOURCE=BASE/'paper_v28_20261014_中文修订_图7三条自由能曲线_三线表.docx'
PRED=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\neural_exports\evaluation\full\A003R.npz')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json');assert sha(DOC)==integration['document_sha256']
assert sha(SOURCE)==integration['source_sha256']=='e86276be39fb476680c0f915d2f155c8be0030599789e035d50a1b83d1dd97ac'
d=Document(DOC);old=Document(SOURCE);text='\n'.join(p.text for p in d.paragraphs)
assert all(s in text for s in ['图7a','图7b','图7c–e','0.080 nats','0.979','0.868','0.105 nats','Q40高估0.368','35.41','43.17','7.76','0.400','0.510'])
assert not any(s in text for s in ['图7a–c','群体范围曲线','范围曲线在E=','覆盖曲线','范围曲线沿用'])
refs=lambda doc:[p.text for p in doc.paragraphs if re.match(r'^\[\d+\]',p.text)]
assert refs(d)==refs(old) and len(refs(d))==27
cells=lambda doc:[[[c.text for c in r.cells] for r in t.rows] for t in doc.tables]
assert cells(d)==cells(old) and len(d.tables)==11
assert len(d.paragraphs)==len(old.paragraphs)
assert sum(p.text!=q.text for p,q in zip(d.paragraphs,old.paragraphs))==4
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def active_images(path):
    with zipfile.ZipFile(path) as z:
        rels=etree.fromstring(z.read('word/_rels/document.xml.rels'));mapping={e.get('Id'):e.get('Target') for e in rels}
        doc=etree.fromstring(z.read('word/document.xml'))
        return [hashlib.sha256(z.read(posixpath.normpath('word/'+mapping[e.get('{'+ns['r']+'}embed')]))).hexdigest() for e in doc.xpath('.//a:blip',namespaces=ns)]
before,after=active_images(SOURCE),active_images(DOC)
assert len(before)==len(after)==8 and [i for i,(a,b) in enumerate(zip(before,after)) if a!=b]==[6]
assert after[6]==sha(FIG/(STEM+'.png'))
z=np.load(PRED,allow_pickle=False);idx=np.flatnonzero((z['chain']=='A')&z['valid'][:,1])
assert len(idx)==66
res=pd.read_csv(FIG/(STEM+'.all66_residue_predictions.csv'));gln=pd.read_csv(FIG/(STEM+'.all6_Gln_predictions.csv'))
assert list(res.residue_index)==list(idx) and res.residue_id.nunique()==66
assert np.array_equal(res.residue_id,z['author_residue_id'][idx].astype(int))
np.testing.assert_allclose(res.MD_Hchi_nats,z['target'][idx,1].astype(float),rtol=0,atol=1e-14)
np.testing.assert_allclose(res.model_Hchi_nats,z['mean'][idx,1].astype(float),rtol=0,atol=1e-14)
assert list(gln.residue_id)==[2,31,40,41,49,62] and list(gln.residue_id)==list(res[res.amino_acid=='Q'].residue_id)
assert res.MD_Hchi_nats.eq(0).sum()==3 and res.loc[res.absolute_error_nats.idxmax(),'residue_id']==40
np.testing.assert_allclose(res.absolute_error_nats,abs(res.model_Hchi_nats-res.MD_Hchi_nats),rtol=0,atol=1e-14)
np.testing.assert_allclose(gln.absolute_error_nats,abs(gln.model_Hchi_nats-gln.MD_Hchi_nats),rtol=0,atol=1e-14)
assert gln.loc[gln.model_Hchi_nats.idxmin(),'residue_id']==41 and gln.loc[gln.model_Hchi_nats.idxmax(),'residue_id']==2
stats=readj(FIG/(STEM+'.statistics.json'))
assert abs(stats['all66_Hchi_MAE_nats']-res.absolute_error_nats.mean())<1e-14
assert abs(stats['all66_Spearman_rho']-res.MD_Hchi_nats.rank().corr(res.model_Hchi_nats.rank()))<1e-14
assert abs(stats['all6_Gln_MAE_nats']-gln.absolute_error_nats.mean())<1e-14
points=pd.read_csv(FIG/(STEM+'.points.csv'));assert len(points)==510
assert points.groupby('panel').size().to_dict()=={'a':66,'b':12,'c':144,'d':144,'e':144}
pa=points.query('panel=="a"').merge(res,on='residue_id',validate='one_to_one')
np.testing.assert_allclose(pa.x,pa.MD_Hchi_nats,rtol=0,atol=1e-14);np.testing.assert_allclose(pa.y,pa.model_Hchi_nats,rtol=0,atol=1e-14)
pb=points.query('panel=="b"').merge(gln,on='residue_id',validate='many_to_one');assert len(pb)==12
for kind,col in [('MD_reference','MD_Hchi_nats'),('model_prediction','model_Hchi_nats')]:
    subset=pb[pb.kind==kind];assert len(subset)==6
    np.testing.assert_allclose(subset.y,subset[col],rtol=0,atol=1e-14)
previous=pd.read_csv(OLD/'figures/fig7_model_selected_three_MD_free_energy_curves.points.csv')
previous['panel']=previous.panel.map({'a':'c','b':'d','c':'e'});previous['residue_id']=previous.residue.str[1:].astype(int)
curves=points[points.panel.isin(list('cde'))]
paired=curves.merge(previous,on=['panel','residue_id','chi','x'],suffixes=('_new','_old'),validate='one_to_one')
assert len(paired)==432 and np.allclose(paired.y_new,paired.y_old,rtol=0,atol=0,equal_nan=True)
assert np.array_equal(paired.low_Fstar_le1_new,paired.low_Fstar_le1_old)
new_widths=pd.read_csv(FIG/(STEM+'.display_widths.csv'))
old_widths=pd.read_csv(OLD/'figures/fig7_model_selected_three_MD_free_energy_curves.display_widths.csv')
width_pairs=new_widths.merge(old_widths,on=['residue','chi'],validate='one_to_one')
assert len(width_pairs)==6 and np.array_equal(width_pairs.origin0_lowFstar_width_degrees,width_pairs.display_origin0_lowFstar_width_degrees)
pdf=PdfReader(TMP/'v29_word.pdf');assert len(pdf.pages)==28
assert all(len((p.extract_text() or '').strip())>100 and '\ufffd' not in (p.extract_text() or '') for p in pdf.pages)
comparisons=[]
for i in range(1,29):
    path=TMP/'render_v29'/f'page-{i}.png';prev=TMP/'render_v28'/f'page-{i}.png'
    exact=ImageChops.difference(Image.open(path),Image.open(prev)).getbbox() is None
    comparisons.append({'page':i,'pixel_identical_to_reviewed_v28':exact,'PNG_sha256':sha(path)})
write(QA/'page_pixel_comparison.json',comparisons)
write(QA/'numerical_content_audit.json',{'status':'passed','document_sha256':sha(DOC),'all66_valid_model_and_reference_values_preserved':True,'all3_valid_zero_values_retained':True,'all6_Gln_including_Q40_largest_error_retained':True,'all510_plotted_data_points_checked':True,'all432_original_MD_curve_points_exact':True,'native_tables_math_comments_preserved':integration['all_tables_math_comments_preserved'],'other7_active_images_identical':True,'references_27_identical':True,'changed_text_paragraphs':4,'pages':28,'stats':stats})
print(json.dumps({'status':'passed','pages':28,'changed_pages':[p['page'] for p in comparisons if not p['pixel_identical_to_reviewed_v28']],'unchanged_pages':[p['page'] for p in comparisons if p['pixel_identical_to_reviewed_v28']]},ensure_ascii=False))
