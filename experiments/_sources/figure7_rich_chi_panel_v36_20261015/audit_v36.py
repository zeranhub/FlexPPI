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
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';EXP=ROOT/'experiment'
OLD=ROOT.parent/'figure7_all_chi_v35_20261015'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
STEM='fig7_rich_chi_scatter_and_MD_curves';OLDSTEM='fig7_all_chi_associations_and_MD_curves'
DOC=BASE/'paper_v36_20261015_中文修订_图7逐角散点与三子图布局_三线表.docx'
SOURCE=BASE/'paper_v35_20261015_中文修订_图7全残基逐χ占据关联_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
readcsv=lambda p:pd.read_csv(p,float_precision='round_trip')
integration=readj(QA/'integration_audit.json')
assert sha(DOC)==integration['document_sha256']
assert sha(SOURCE)==integration['source_sha256']=='35c3d0a6ae967c1bd76e45fb826245233947f07d0eff56660d7291bd478c9ecf'
d=Document(DOC);old=Document(SOURCE)
assert len(d.tables)==11 and len(d.paragraphs)==len(old.paragraphs)
changed_paragraphs=[i for i,(p,q) in enumerate(zip(d.paragraphs,old.paragraphs)) if p.text!=q.text]
assert changed_paragraphs==[51,100,106],changed_paragraphs
text='\n'.join(p.text for p in d.paragraphs)
assert all('六个谷氨酰胺中' in d.paragraphs[i].text for i in [51,100,106])
assert '仅标出Q2与Q41' in d.paragraphs[51].text
assert '图7a右' not in text and '图7b–d' not in text
assert '（图7b）' in text and '图7c进一步展开' in text
assert '仅标出Q2与Q41' in d.paragraphs[106].text
assert '0.105 nats' not in d.paragraphs[100].text
assert '0.447–0.525' in text and '0.640' in text and '0.127和−0.183' in text
assert all(value in text for value in ['0.080 nats','0.400','0.510','35.41','43.17','7.76'])
assert d.paragraphs[98].text==old.paragraphs[98].text
assert d.paragraphs[99].text==old.paragraphs[99].text
assert all(d.paragraphs[i].text==old.paragraphs[i].text for i in [102,103,104])
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
for ext in ['png','pdf','svg']:assert sha(FIG/(STEM+'.'+ext))==fa['output_sha256'][ext]
# Axis labels are intentionally reduced to the two frozen-prediction case endpoints.
labels=fa['a_labels']
assert fa['panel_letters']==['a','b','c']
if not labels:
    raise AssertionError(('Axis label audit key missing',list(fa)))
assert [r['anchor_x_data'] for r in labels]==[2,41]
assert [r['label_x_data'] for r in labels]==[2,41]
newpoints=readcsv(FIG/(STEM+'.points.csv'));oldpoints=readcsv(OLD/'figures'/(OLDSTEM+'.points.csv'))
assert len(newpoints)==723
assert newpoints.groupby('panel').size().to_dict()=={'a_Hchi':132,'b_scatter':151,'b_summary':8,'c':432}
original=oldpoints[oldpoints.panel!='a_chi_association'].copy()
original['panel']=original.panel.replace({'b':'c','c':'c','d':'c'})
retained=newpoints[newpoints.panel.isin(['a_Hchi','c'])].copy()
key=['panel','residue_id','kind','chi','x']
paired=retained.merge(original,on=key,suffixes=('_new','_old'),validate='one_to_one')
assert len(paired)==564 and np.allclose(paired.y_new,paired.y_old,rtol=0,atol=0,equal_nan=True)
assert paired.low_Fstar_le1_new.fillna('absent').equals(paired.low_Fstar_le1_old.fillna('absent'))
assert pd.read_csv(FIG/(STEM+'.display_widths.csv')).equals(pd.read_csv(OLD/'figures'/(OLDSTEM+'.display_widths.csv')))
assert newpoints.query('panel=="a_Hchi" and kind=="MD_reference"').y.eq(0).sum()==3
# All scatter values are exact copies of the already audited full angle table.
angles=readcsv(EXP/'defined_angles_long.csv');slots=readcsv(EXP/'all_angle_slots_long.csv');stat=readcsv(EXP/'figure_statistics.csv')
assert len(angles)==151 and len(slots)==264
assert angles.groupby('chi').size().to_dict()=={'chi1':66,'chi2':52,'chi3':23,'chi4':10}
assert len(angles.author_residue_id.unique())==66
assert len(slots[~slots.chi_defined])==113
assert (slots.loc[~slots.chi_defined,'fine_entropy_four_origins_mean_nats'].isna()).all()
assert angles[angles.author_residue_id==40].chi.tolist()==['chi1','chi2','chi3']
scatter=newpoints[newpoints.panel=='b_scatter'].copy()
assert scatter.groupby('chi').size().to_dict()=={1:66,2:52,3:23,4:10}
values=[]
for _,row in stat.iterrows():
    group=angles[angles.chi==row.chi]
    assert len(group)==row.n
    actual=scatter[scatter.chi==row.chi_index].sort_values('residue_id')
    source_group=group.sort_values('author_residue_id')
    assert np.array_equal(actual.residue_id,source_group.author_residue_id)
    np.testing.assert_allclose(actual.x,source_group.predicted_Hchi_nats_AA_centered,rtol=0,atol=0)
    np.testing.assert_allclose(actual.y,source_group.fine_entropy_four_origins_mean_nats_AA_centered,rtol=0,atol=0)
    for metric,target,col in [('fine_entropy','fine_entropy_four_origins_mean_nats','rho_entropy_typecentered'),('low_Fstar_width','width_four_origins_mean_degrees','rho_width_typecentered')]:
        x=group.predicted_Hchi_nats_AA_centered.to_numpy()
        y=group[target+'_AA_centered'].to_numpy()
        rho=float(np.corrcoef(pd.Series(x).rank(method='average'),pd.Series(y).rank(method='average'))[0,1])
        assert abs(rho-row[col])<1e-14,(row.chi,metric,rho,row[col])
        summary=newpoints[(newpoints.panel=='b_summary')&(newpoints.kind==metric)&(newpoints.chi==row.chi_index)]
        assert len(summary)==1 and abs(summary.rho.iloc[0]-rho)<1e-14
        assert int(summary.n_defined.iloc[0])==int(row.n)
        values.append({'chi':row.chi,'metric':metric,'n':int(row.n),'rho':rho})
assert stat.loc[stat.chi=='chi4','rho_width_typecentered'].iloc[0]<0
assert abs(stat.rho_entropy_typecentered.iloc[:3].min()-.44718640594313047)<1e-14
assert abs(stat.rho_entropy_typecentered.iloc[:3].max()-.52470355731225293)<1e-14
assert stat.loc[stat.rho_width_typecentered.idxmax(),'chi']=='chi3'
experiment_audit=readj(EXP/'audit.json');assert experiment_audit['status']=='passed' and experiment_audit['checks_passed']==1149
for name,digest in experiment_audit['output_sha256'].items():assert sha(EXP/name)==digest
for name,item in experiment_audit['original_source_provenance'].items():assert sha(item['read_path'])==item['original_sha256']
pdf=PdfReader(TMP/'v36_word.pdf');n_pages=len(pdf.pages)
assert 28<=n_pages<=29,n_pages
assert all(len((page.extract_text() or '').strip())>100 and '\ufffd' not in (page.extract_text() or '') for page in pdf.pages)
comparisons=[]
for i in range(1,n_pages+1):
    path=TMP/'render_v36'/f'page-{i}.png';prev=TMP/'render_v35'/f'page-{i}.png'
    exact=prev.exists() and ImageChops.difference(Image.open(path),Image.open(prev)).getbbox() is None
    comparisons.append({'page':i,'pixel_identical_to_reviewed_v35':exact,'PNG_sha256':sha(path)})
write(QA/'page_pixel_comparison.json',comparisons)
write(QA/'numerical_content_audit.json',{'status':'passed','document_sha256':sha(DOC),'all564_original_curve_profile_values_exact':True,'all8_displayed_associations_independently_recomputed':values,'all151_displayed_scatter_values_exact':True,'panels_a_b_c':True,'all151_defined_angles_retained':True,'missing_slots_NaN':113,'all66_residues_including_3Pro_and_Q40_retained':True,'only_Q2_Q41_axis_labels':True,'all432_original_MD_curve_points_exact':True,'new_MD_or_training':False,'native_tables_math_comments_preserved':integration['all_tables_math_comments_preserved'],'other7_active_images_identical':True,'changed_text_paragraphs':changed_paragraphs,'pages':n_pages,'experimental_checks_passed':1149})
print(json.dumps({'status':'passed','pages':n_pages,'changed_pages':[p['page'] for p in comparisons if not p['pixel_identical_to_reviewed_v35']]},ensure_ascii=False))
