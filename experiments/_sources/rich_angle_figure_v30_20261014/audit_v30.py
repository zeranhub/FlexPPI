from pathlib import Path
import hashlib,json,re,zipfile,posixpath
import numpy as np
import pandas as pd
from docx import Document
from lxml import etree
from PIL import Image,ImageChops
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';DATA=ROOT/'experiment'
OLD=ROOT.parent/'all_residues_figure_v29_20261014'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
STEM='fig7_full_residue_and_chi_occupancy'
DOC=BASE/'paper_v30_20261014_中文修订_全残基逐角占据与自由能验证_三线表.docx'
SOURCE=BASE/'paper_v29_20261014_中文修订_全残基与同类型预测及自由能曲线_三线表.docx'
PRED=OLD/'source/A003R.npz'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json');assert sha(DOC)==integration['document_sha256']
assert sha(SOURCE)==integration['source_sha256']=='dc39f5dacad17bbb284621083c39834186e87d53a3a73a49d2cf5992006fa904'
d=Document(DOC);old=Document(SOURCE);text='\n'.join(p.text for p in d.paragraphs)
assert all(s in text for s in ['图7a','图7b','图7c–e','0.080 nats','0.979','0.868','0.105 nats','Q40高估0.368','全部151个可定义χ角','χ1–χ3','35.41','43.17','7.76','0.400','0.510'])
assert not any(s in text for s in ['图7a–c','群体范围曲线','范围曲线在E=','覆盖曲线','范围曲线沿用'])
refs=lambda doc:[p.text for p in doc.paragraphs if re.match(r'^\[\d+\]',p.text)]
assert refs(d)==refs(old) and len(refs(d))==27
cells=lambda doc:[[[c.text for c in r.cells] for r in t.rows] for t in doc.tables]
assert cells(d)==cells(old) and len(d.tables)==11
assert len(d.paragraphs)==len(old.paragraphs)
assert sum(p.text!=q.text for p,q in zip(d.paragraphs,old.paragraphs))==3
ns={'a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def active_images(path):
    with zipfile.ZipFile(path) as z:
        rels=etree.fromstring(z.read('word/_rels/document.xml.rels'));mapping={e.get('Id'):e.get('Target') for e in rels}
        doc=etree.fromstring(z.read('word/document.xml'))
        return [hashlib.sha256(z.read(posixpath.normpath('word/'+mapping[e.get('{'+ns['r']+'}embed')]))).hexdigest() for e in doc.xpath('.//a:blip',namespaces=ns)]
before,after=active_images(SOURCE),active_images(DOC)
assert len(before)==len(after)==8 and [i for i,(a,b) in enumerate(zip(before,after)) if a!=b]==[6]
assert after[6]==sha(FIG/(STEM+'.png'))
z=np.load(DATA/'ubiquitin_all66_chi_501frames.npz',allow_pickle=False);p=np.load(PRED,allow_pickle=False)
idx=np.flatnonzero((p['chain']=='A')&p['valid'][:,1]);assert len(idx)==66 and np.array_equal(idx,z['residue_index'])
assert np.array_equal(p['target'][idx,1].astype(float),z['frozen_MD_joint_Hchi_nats'])
assert np.array_equal(p['mean'][idx,1].astype(float),z['predicted_Hchi_nats'])
assert np.array_equal(p['author_residue_id'][idx].astype(int),z['author_residue_id'].astype(int))
assert np.sum(z['chi_defined'])==151 and np.sum(~z['chi_defined'])==113
assert z['chi_defined'].sum(axis=0).tolist()==[66,52,23,10]
assert np.isnan(z['chi_degrees'][:,~z['chi_defined']]).all()
assert np.isnan(z['probability_origin0'][~z['chi_defined']]).all()
assert np.isnan(z['width_per_origin_degrees'][~z['chi_defined']]).all()
assert np.allclose(z['probability_origin0'][z['chi_defined']].sum(axis=1),1,rtol=0,atol=1e-14)
assert (z['frozen_MD_joint_Hchi_nats']==0).sum()==3
points=pd.read_csv(FIG/(STEM+'.points.csv'));assert len(points)==2202
assert points.groupby('panel').size().to_dict()=={'a_Hchi':132,'a_error':66,'a_width':264,'b_Hchi':12,'b_probability':1296,'c':144,'d':144,'e':144}
lookup={int(site):i for i,site in enumerate(z['author_residue_id'])}
for r in points.itertuples(index=False):
    i=lookup[int(r.residue_id)]
    if r.panel in ['a_Hchi','b_Hchi']:
        expected=z['predicted_Hchi_nats'][i] if r.kind=='model_prediction' else z['frozen_MD_joint_Hchi_nats'][i]
        actual=r.y if r.panel=='a_Hchi' else r.x
        assert abs(actual-expected)<1e-14
    elif r.panel=='a_error':assert abs(r.y-(z['predicted_Hchi_nats'][i]-z['frozen_MD_joint_Hchi_nats'][i]))<1e-14
    elif r.panel=='a_width':
        q=int(r.chi)-1
        assert bool(r.chi_defined)==bool(z['chi_defined'][i,q])
        if z['chi_defined'][i,q]:assert abs(r.y-np.mean(z['width_per_origin_degrees'][i,q]))<1e-14
        else:assert np.isnan(r.y)
    elif r.panel=='b_probability':
        q=int(r.chi)-1;k=int(round((r.x+177.5)/5))
        assert abs(r.y-z['probability_origin0'][i,q,k])<1e-14 and r.count==z['counts_origin0'][i,q,k]
        assert bool(r.low_Fstar_le1)==bool(np.isfinite(z['Fstar_origin0'][i,q,k]) and z['Fstar_origin0'][i,q,k]<=1)
order=pd.read_csv(FIG/(STEM+'.all6_prediction_order.csv'))
assert list(order.author_residue_id)==[41,31,40,62,49,2]
assert set(order.author_residue_id)==set(z['author_residue_id'][z['amino_acid']=='Q'].astype(int))
assert points.query('panel=="a_width"').y.notna().sum()==151
assert points.query('panel=="a_Hchi" and kind=="MD_reference"').y.eq(0).sum()==3
previous=pd.read_csv(OLD/'figures/fig7_all_residues_Gln_and_MD_free_energy.points.csv')
curves=points[points.panel.isin(list('cde'))]
old_curves=previous[previous.panel.isin(list('cde'))]
paired=curves.merge(old_curves,on=['panel','residue_id','chi','x'],suffixes=('_new','_old'),validate='one_to_one')
assert len(paired)==432 and np.allclose(paired.y_new,paired.y_old,rtol=0,atol=0,equal_nan=True)
assert np.array_equal(paired.low_Fstar_le1_new,paired.low_Fstar_le1_old)
assert pd.read_csv(FIG/(STEM+'.display_widths.csv')).equals(pd.read_csv(OLD/'figures/fig7_all_residues_Gln_and_MD_free_energy.display_widths.csv'))
datareview=readj(DATA/'independent_verification.json');assert datareview['status']=='passed'
assert datareview['old_six_GLN_maximum_circular_error_degrees']==0.0 and datareview['old_two_site_histogram_rows_checked']==432
for filename,info in datareview['files'].items():assert sha(DATA/filename)==info['sha256']
pdf=PdfReader(TMP/'v30_word.pdf');assert len(pdf.pages)==28
assert all(len((page.extract_text() or '').strip())>100 and '\ufffd' not in (page.extract_text() or '') for page in pdf.pages)
comparisons=[]
for i in range(1,29):
    path=TMP/'render_v30'/f'page-{i}.png';prev=TMP/'render_v29'/f'page-{i}.png'
    exact=ImageChops.difference(Image.open(path),Image.open(prev)).getbbox() is None
    comparisons.append({'page':i,'pixel_identical_to_reviewed_v29':exact,'PNG_sha256':sha(path)})
write(QA/'page_pixel_comparison.json',comparisons)
write(QA/'numerical_content_audit.json',{'status':'passed','document_sha256':sha(DOC),'all66_valid_model_and_reference_values_preserved':True,'all3_valid_zero_values_retained':True,'all151_angles_and113missing_slots_preserved':True,'all6_Gln_all3angles_and72bins_retained':True,'all2202_source_points_checked':True,'all432_original_MD_curve_points_exact':True,'independent_angle_verification_passed':True,'native_tables_math_comments_preserved':integration['all_tables_math_comments_preserved'],'other7_active_images_identical':True,'references_27_identical':True,'changed_text_paragraphs':3,'pages':28})
print(json.dumps({'status':'passed','pages':28,'changed_pages':[p['page'] for p in comparisons if not p['pixel_identical_to_reviewed_v29']],'unchanged_pages':[p['page'] for p in comparisons if p['pixel_identical_to_reviewed_v29']]},ensure_ascii=False))
