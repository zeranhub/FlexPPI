"""Read-only checks of v21 manuscript against independently certified results."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone, timedelta
import hashlib, json, re, zipfile
import pandas as pd
from docx import Document
from docx.oxml.ns import qn
from lxml import etree
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent.parent
QA = ROOT / 'manuscript_QA'
DOCX = BASE / 'paper_v21_20261012_中文修订_构象分布与自由能分析_三线表.docx'
SOURCE = BASE / 'paper_v20_20261012_中文修订_联合模型_三线表.docx'
TEMP = Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
doc, old = Document(DOCX), Document(SOURCE)
checks = []
def check(name, condition, detail=None):
    checks.append({'name': name, 'passed': bool(condition), 'detail': detail})
    assert condition, (name, detail)
def read_csv(relative):
    return pd.read_csv(ROOT / relative)
text = '\n'.join(p.text for p in doc.paragraphs)
audit = json.loads((QA / 'integration_audit.json').read_text(encoding='utf8'))
check('Frozen source and final document hashes', sha(DOCX)==audit['sha256'] and sha(SOURCE)==audit['source_document_sha256'])
for path, digest in audit['scientific_sources_sha256'].items():
    check('Integrated scientific source hash exact', sha(path)==digest, path)
check('Transferred common-residue gzip equals remote certified file', sha(ROOT/'oof_analysis/paired_common_residues.csv.gz')=='1b4061f68901f35a5edd72023b3fa0d5cc907418f2ae6595adefe32042b5b868')
check('Exactly 15 native Word tables', len(doc.tables)==15)
cells = lambda t: [[c.text for c in r.cells] for r in t.rows]
check('All 13 original tables retain every cell value', all(cells(doc.tables[i])==cells(old.tables[i]) for i in range(13)))

classes=read_csv('oof_analysis/reference_class_record_macro.csv')
classes=classes[classes.dataset=='all'].set_index('reference_class')
keys=['low_RMSF_low_Hchi','low_RMSF_high_Hchi','high_RMSF_low_Hchi','high_RMSF_high_Hchi','unclassified']
names=['低RMSF／低Hχ','低RMSF／高Hχ','高RMSF／低Hχ','高RMSF／高Hχ','未分类']
expected=[]
for key,name in zip(keys,names):
    row=classes.loc[key]
    expected.append([name,f'{int(row.common_residues):,}',f'{100*row.reference_fraction_record_macro:.3f}',f'{row.rmsf_mae_record_macro:.3f}',f'{row.hchi_mae_record_macro:.3f}',f'{100*row.rank_class_recall_record_macro:.1f}' if key!='unclassified' else '无'])
check('Every new table S12 value from source CSV', cells(doc.tables[13])[1:]==expected)
check('Reference categories partition full common mask', int(classes.common_residues.sum())==824988)
discordant=classes.loc[keys[1:3]]
check('Discordant category counts and macro/micro percentages', int(discordant.common_residues.sum())==47354 and f'{100*discordant.reference_fraction_record_macro.sum():.3f}%' in text and f'{100*discordant.common_residues.sum()/824988:.3f}%' in text)

joint=read_csv('trajectory_analysis/manuscript_summary_table.csv')
joint=joint[joint.dataset=='all'].set_index('shift_degrees')
nonpro=read_csv('trajectory_analysis/nonpro_chi1_source_summary.csv')
nonpro=nonpro[(nonpro.dataset=='all')&(nonpro.replica_mode=='replica1')].set_index('shift_degrees')
expected=[]
for phase in [0,30,60,90]:
    r,q=joint.loc[phase],nonpro.loc[phase]
    expected.append([f'{phase}°',f'{r.joint_H_mean_nats:.3f} ± {r.joint_H_record_sd_nats:.3f}',f'{q.h_chi1_mean:.3f}',f'{100*q.chi1_boundary_10deg_fraction_mean:.2f}',f'{q.apparent_chi1_changes_mean:.1f}',f'{r.half_joint_JSD_mean_nats:.3f}'])
check('Every new table S13 value from correct non-Pro and joint CSVs',cells(doc.tables[14])[1:]==expected)
for phase in [0,30,60,90]:
    check('Narrative joint partition means exact',f'{joint.loc[phase].joint_H_mean_nats:.3f}' in text,phase)
for phase in [0,60]:
    q=nonpro.loc[phase]
    check('Narrative conventional chi1 means and changes exact',f'{q.h_chi1_mean:.3f}' in text and f'{q.apparent_chi1_changes_mean:.1f}' in text,phase)
    check('Narrative state-conditional variance and window entropy gap',f'{100*joint.loc[phase].joint_CA_between_fraction_mean:.1f}%' in text and f'{joint.loc[phase].full_H_minus_weighted_window_H_mean_nats:.3f}' in text,phase)
ub=read_csv('case_analysis/ubiquitin_bin_sensitivity.csv').set_index(['residue','phase'])
for residue in ['Q41','Q62']:
    r=ub.loc[(residue,60)]
    check('Ubiquitin narrative uses five-window joint and chi1 means',f'{r.joint_window_mean_h:.3f}' in text and f'{r.chi1_window_mean_h:.3f}' in text,residue)

for index,t in enumerate(doc.tables):
    check('No native table grid or shading',not t._tbl.xpath('.//w:shd | ./w:tblPr/w:tblBorders'),index)
    for i,row in enumerate(t.rows):
        for cell in row.cells:
            borders=cell._tc.tcPr.find(qn('w:tcBorders'))
            check('Three-line table borders defined',borders is not None,index)
            values={etree.QName(e).localname:e.get(qn('w:val')) for e in borders}
            visible={k for k,v in values.items() if v not in {'nil','none'}}
            expected=({'top','bottom'} if i==0 else {'bottom'} if i==len(t.rows)-1 else set())
            check('Only top/header/bottom horizontal table lines',visible==expected,(index,i))

ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(DOCX) as b:
    ax,bx=etree.fromstring(a.read('word/document.xml')),etree.fromstring(b.read('word/document.xml'))
    check('User comments preserved byte-for-byte',a.read('word/comments.xml')==b.read('word/comments.xml'))
    anchors=lambda x:Counter((e.tag,e.get(qn('w:id'))) for e in x.iter() if e.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    check('All comment anchors preserved',anchors(ax)==anchors(bx))
    eq=lambda x:[etree.tostring(e,method='c14n') for e in x.xpath('.//m:oMath',namespaces=ns)]
    check('Native math equations preserved exactly',eq(ax)==eq(bx))
    check('Old scientific image parts preserved',all(a.read(n)==b.read(n) for n in a.namelist() if n.startswith('word/media/')))
    check('Exactly 20 inline scientific figures',len(bx.xpath('.//a:blip',namespaces=ns))==20)

figure_captions=[re.match(r'^图\s*(S?\d+)\s*[:：]',p.text).group(1) for p in doc.paragraphs if re.match(r'^图\s*(S?\d+)\s*[:：]',p.text)]
table_captions=[re.match(r'^表\s*(S?\d+)\s*[:：]',p.text).group(1) for p in doc.paragraphs if re.match(r'^表\s*(S?\d+)\s*[:：]',p.text)]
check('Main and supplementary figure numbers consecutive',figure_captions==[str(i) for i in range(1,8)]+['S'+str(i) for i in range(1,14)],figure_captions)
check('Main and supplementary table numbers consecutive',table_captions==['1','2']+['S'+str(i) for i in range(1,14)],table_captions)
for prefix,captions in [('图',figure_captions),('表',table_captions)]:
    identifiers=set(re.findall(prefix+r'\s*(S?\d+)',text))
    check('All figure/table cross-reference numbers resolve',identifiers<=set(captions),(prefix,sorted(identifiers)))
refs=[int(m.group(1)) for p in doc.paragraphs if (m:=re.match(r'^\[(\d+)\]',p.text))]
check('35 unique numbered references',refs==list(range(1,36)))
check('New primary-reference DOIs and print years retained',all(s in text for s in ['10.1021/jp4105823','10.1002/pro.565','10.1021/ja5024783','10.1073/pnas.0910390107','107(3):1088–1093, 2010']))
check('Default joint model and editorial direction retained','FlexPPI采用联合训练作为默认方案' in text and '家族隔离后成绩下降' not in text and '联合训练收益较小' not in text)
check('Methods separate trajectory F* from model outputs','F*由轨迹的所选坐标计算，RMSF与Hχ由模型预测' in text)
check('No claims of learned energy-barrier prediction','预测能垒' not in text and '预测能量势垒' not in text)

pdf=PdfReader(TEMP/'v21_word.pdf')
check('Word render contains 49 pages',len(pdf.pages)==49)
page_text=[p.extract_text() or '' for p in pdf.pages]
check('No blank pages or missing-glyph replacement characters',all(len(s.strip())>20 and '\ufffd' not in s for s in page_text))
reuse=json.loads((QA/'page_pixel_reuse.json').read_text())
viewed=set(range(1,27))|{40,45,46,47,48,49}
coverage=[{'page':p['page'],'review':'visually_inspected_final_png' if p['page'] in viewed else 'pixel_identical_body_to_previously_inspected_v20','v20_page':p['reuse_v20']} for p in reuse]
check('All 49 page bodies visually verified',all(p['page'] in viewed or p['reuse_v20'] is not None for p in reuse))

report={'status':'passed','document':str(DOCX),'document_sha256':sha(DOCX),'verified_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'checks_passed':len(checks),'checks':checks,'tables':15,'figures':20,'pages':49,'page_review':coverage,'review_limits':'Numerical and manuscript audit is distinct from independent experiment verification; reused pages are byte-identical rendered bodies with only changed footer page numbers.'}
(QA/'final_content_and_layout_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'status':report['status'],'checks':len(checks),'pages':49,'tables':15,'figures':20},ensure_ascii=False))
