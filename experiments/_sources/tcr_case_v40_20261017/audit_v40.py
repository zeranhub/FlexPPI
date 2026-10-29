from pathlib import Path
from collections import Counter
import hashlib,json,zipfile,re,posixpath,sys
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
sys.path.insert(0,str(ROOT.parent/'completed_results_20261012/plot_runtime'))
import pandas as pd
import numpy as np
OLD=BASE/'paper_v39_20261017_中文修订_读者表述与英文图表_三线表.docx'
NEW=BASE/'paper_v40_20261017_中文修订_TCR识别环境变化实例_三线表.docx'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
read=lambda p:pd.read_csv(p,float_precision='round_trip')
a=Document(OLD);b=Document(NEW)
assert len(a.tables)==len(b.tables)==10
for x,y in zip(a.tables,b.tables):assert etree.tostring(x._tbl,method='c14n')==etree.tostring(y._tbl,method='c14n')
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def active(path):
    with zipfile.ZipFile(path) as z:
        tree=etree.fromstring(z.read('word/document.xml'));rel={el.get('Id'):el.get('Target') for el in etree.fromstring(z.read('word/_rels/document.xml.rels'))}
        return [hashlib.sha256(z.read(posixpath.normpath('word/'+rel[el.get(qn('r:embed'))]))).hexdigest() for el in tree.xpath('.//a:blip',namespaces=ns)]
before,after=active(OLD),active(NEW)
assert len(before)==7 and len(after)==8
assert after[:5]==before[:5] and after[6:]==before[5:] and after[5]==sha(ROOT/'figure6/fig6_tcr_case.png')
with zipfile.ZipFile(OLD) as z1,zipfile.ZipFile(NEW) as z2:
    x=etree.fromstring(z1.read('word/document.xml'));y=etree.fromstring(z2.read('word/document.xml'))
    maths=lambda t:[etree.tostring(el,method='c14n') for el in t.xpath('.//m:oMath',namespaces=ns)]
    assert maths(x)==maths(y) and len(maths(y))==36
    anchors=lambda t:Counter((el.tag,el.get(qn('w:id'))) for el in t.iter() if el.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    assert anchors(x)==anchors(y)
    assert etree.tostring(etree.fromstring(z1.read('word/comments.xml')),method='c14n')==etree.tostring(etree.fromstring(z2.read('word/comments.xml')),method='c14n')
text='\n'.join(p.text for p in b.paragraphs)
assert not any(s in text for s in ['冻结','发布权重','发布模型'])
assert [int(re.match(r'^图(\d+)：',p.text)[1]) for p in b.paragraphs if re.match(r'^图\d+：',p.text)]==list(range(1,8))
assert len([p for p in b.paragraphs if p.text.startswith('图S1：')])==1
assert [int(re.match(r'^\[(\d+)\]',p.text)[1]) for p in b.paragraphs if re.match(r'^\[\d+\]',p.text)]==list(range(1,29))
assert [p.text[:3] for p in b.paragraphs if re.match(r'^3\.\d ',p.text)]==['3.1','3.2','3.3','3.4','3.5','3.6']
assert all(s in text for s in ['0.979（图7a）','33,142','0.510','0.400','7.76','0.945','0.093 Å','1.099','−0.018','Q9的熵增幅未被模型重现'])
assert '[28]' in text and 'doi:10.1074/jbc.M112.447185' in text
p=read(ROOT/'figure6/all13_model_and_MD_changes.csv');rep=read(ROOT/'sb27_sb47/all13_peptide_per_replica_labels.csv')
assert len(p)==13 and p.hchi_valid.sum()==11
for pdb,key in [('2ak4','a'),('4jry','b')]:
    source=np.load(ROOT/f'family_predictions/full/itrace_{pdb}.npz')
    mask=source['chain']=='C';assert mask.sum()==13
    for metric,j,col in [('rmsf',0,'rmsf_A'),('hchi',1,'Hchi_nats')]:
        vals=rep[rep.pdb.eq(pdb)].groupby('position')[col]
        target=source['target'][mask,j]
        assert np.allclose(vals.mean(),target,equal_nan=True,atol=3e-6)
        valid=p[f'{metric}_valid'].to_numpy(bool)
        assert np.array_equal(valid,source['valid'][mask,j])
        assert np.allclose(vals.var(ddof=1).to_numpy()[valid],source['label_variance'][mask,j][valid],atol=3e-6)
        assert np.allclose(p.loc[valid,f'{metric}_full_{key}'],source['mean'][mask,j][valid],atol=1e-7)
        assert np.allclose(p.loc[valid,f'{metric}_md_{key}'],target[valid],atol=1e-7)
    for file in [ROOT/f'sb27_sb47/{pdb}_model_input.pdb',ROOT/f'sb27_sb47/{pdb}_pMHC_only.pdb']:assert file.exists()
    full=(ROOT/f'sb27_sb47/{pdb}_model_input.pdb').read_text()
    retained=''.join(l for l in full.splitlines(True) if l[:6].strip() in ['ATOM','HETATM','TER'] and l[21:22] in 'ABC')+'END\n'
    assert retained==(ROOT/f'sb27_sb47/{pdb}_pMHC_only.pdb').read_text()
    repro=readj(ROOT/'sb27_sb47/context_inference_audit.json')['full_reproduction'][pdb]
    assert repro['max_absolute_error']<2e-5 and repro['input_sha256']==sha(ROOT/f'sb27_sb47/{pdb}_model_input.pdb')
    assert ''.join(source['amino_acid'][mask])=='LPEPLPQGQLTAY'
    if pdb=='2ak4':mhc=source['amino_acid'][source['chain']=='A']
    else:assert np.array_equal(mhc,source['amino_acid'][source['chain']=='A'])
    for metric in ['rmsf','hchi']:
        assert np.allclose(p[f'{metric}_full_delta'],p[f'{metric}_full_b']-p[f'{metric}_full_a'],equal_nan=True)
        assert np.allclose(p[f'{metric}_md_delta'],p[f'{metric}_md_b']-p[f'{metric}_md_a'],equal_nan=True)
        expected=np.sqrt((p[f'{metric}_sd_a']**2+p[f'{metric}_sd_b']**2)/3)
        assert np.allclose(p[f'{metric}_delta_sem'],expected,equal_nan=True)
assert np.isclose(np.corrcoef(p.rmsf_full_delta.rank(),p.rmsf_md_delta.rank())[0,1],.945054945054945)
assert np.isclose(abs(p.rmsf_full_delta-p.rmsf_md_delta).mean(),.09310372059161846)
pairs=read(ROOT/'local_pair_screen/biological_candidate_pairs.csv')
assert pairs.category.eq('same_pMHC_different_TCR').sum()==124
ctx=read(ROOT/'sb27_sb47/context_inference.csv');q=ctx[ctx.position.eq(7)]
qdiff={name:float(g.set_index('pdb').loc['4jry'].Hchi-g.set_index('pdb').loc['2ak4'].Hchi) for name,g in q.groupby('input')}
assert np.isclose(qdiff['full'],1.098828911781311) and np.isclose(qdiff['pMHC_only'],-.017818450927734375)
raw=readj(ROOT/'sb27_sb47/raw_trajectory_audit.json');assert raw['status']=='passed' and raw['raw_topology_and_trajectory_files_unchanged']
figaudit=readj(ROOT/'figure6/figure6_data_audit.json');assert figaudit['status']=='passed'
pdf=PdfReader(TMP/'v40_word.pdf');assert len(pdf.pages)==30
assert all(len(p.extract_text().strip())>100 and '\ufffd' not in p.extract_text() for p in pdf.pages)
report={'status':'passed','document':str(NEW),'document_sha256':sha(NEW),'pages':30,'main_figures':7,'supplementary_figures':1,'all7_existing_images_byte_exact':True,'all10_native_three_line_tables_byte_exact':True,'native_math_byte_exact':36,'original_comments_and_anchors_preserved':True,'references':28,'all_chains_and_peptide_sequences_verified':True,'raw_MD_labels_mean_and_variance_recomputed_all13positions':True,'all13_positions_and11valid_sidechains_retained':True,'complete_pair_screen_available':True,'same_model_reproduction_tolerance':2e-5,'same_model_TCR_removal_preserves_all_pMHC_coordinate_lines':True,'Q7_predicted_delta_full':qdiff['full'],'Q7_predicted_delta_TCR_removed':qdiff['pMHC_only'],'new_training':0,'new_MD_production':0,'new_same_model_inferences':4,'numbering_and_cross_references':'passed','visual_review':'pending'}
write(QA/'numerical_structural_audit.json',report)
print(json.dumps(report,ensure_ascii=False,indent=2))
