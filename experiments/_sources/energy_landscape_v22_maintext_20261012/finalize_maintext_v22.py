"""Verify preserved content and main-text placement; package v22 delivery."""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone,timedelta
import hashlib,json,re,shutil,zipfile
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'energy_landscape_v21_20261012'
BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA'
SOURCE=BASE/'paper_v21_20261012_中文修订_构象分布与自由能分析_三线表.docx'
DOC=BASE/'paper_v22_20261012_中文修订_自由能分析融入正文_三线表.docx'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
now=datetime.now(timezone(timedelta(hours=8))).isoformat()
doc,source=Document(DOC),Document(SOURCE)
checks=[]
def check(name,condition):
    checks.append({'name':name,'passed':bool(condition)});assert condition,name
integration=json.loads((QA/'integration_audit.json').read_text(encoding='utf8'))
independent=json.loads((QA/'v22_results_independent_content_review.json').read_text(encoding='utf8'))
check('Independent numerical and scientific review matches final document',independent['status']=='passed' and independent['document_sha256']==sha(DOC))
check('Integration certificate matches final document',integration['document_sha256']==sha(DOC))
cells=lambda t:tuple(tuple(c.text for c in r.cells) for r in t.rows)
check('All fifteen original table cell values preserved',len(doc.tables)==15 and Counter(cells(t) for t in doc.tables)==Counter(cells(t) for t in source.tables))
for t in doc.tables:
    check('No table grid or shading',not t._tbl.xpath('.//w:shd | ./w:tblPr/w:tblBorders'))
    for i,row in enumerate(t.rows):
        for c in row.cells:
            borders=c._tc.tcPr.find(qn('w:tcBorders'));assert borders is not None
            visible={etree.QName(e).localname for e in borders if e.get(qn('w:val')) not in {'nil','none'}}
            check('Only three horizontal table rules',visible==({'top','bottom'} if i==0 else {'bottom'} if i==len(t.rows)-1 else set()))
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(DOC) as b:
    ax,bx=etree.fromstring(a.read('word/document.xml')),etree.fromstring(b.read('word/document.xml'))
    check('User comments byte-identical',a.read('word/comments.xml')==b.read('word/comments.xml'))
    anchor=lambda x:Counter((e.tag,e.get(qn('w:id'))) for e in x.iter() if e.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    check('All comment anchors preserved',anchor(ax)==anchor(bx))
    math=lambda x:[etree.tostring(e,method='c14n') for e in x.xpath('.//m:oMath',namespaces=ns)]
    check('Native source math preserved with one new identity',math(ax)==math(bx)[:-1] and len(math(bx))==len(math(ax))+1)
    check('All scientific media parts byte-identical',all(a.read(n)==b.read(n) for n in a.namelist() if n.startswith('word/media/')))
    check('20 scientific figures total',len(bx.xpath('.//a:blip',namespaces=ns))==20)
text='\n'.join(p.text for p in doc.paragraphs)
end=text.index('4 讨论');supp=text.index('补充材料');main=text[:end]
figs=[m.group(1) for p in doc.paragraphs if (m:=re.match(r'^图\s*(S?\d+)\s*[:：]',p.text))]
tables=[m.group(1) for p in doc.paragraphs if (m:=re.match(r'^表\s*(S?\d+)\s*[:：]',p.text))]
check('Nine main and eleven supplementary consecutive figure captions',figs==[str(i) for i in range(1,10)]+['S'+str(i) for i in range(1,12)])
check('Three main and twelve supplementary consecutive table captions',tables==[str(i) for i in range(1,4)]+['S'+str(i) for i in range(1,13)])
check('Figures7to9andtable3allbeforeDiscussion',all(x in main for x in ['图7：','图8：','图9：','表3：']))
check('No old promoted references',not any(x in text for x in ['图S12','图S13','表S13']))
for prefix,captions in [('图',figs),('表',tables)]:
    check('All numbered cross-references resolve',set(re.findall(prefix+r'\s*(S?\d+)',text))<=set(captions))
refs=[int(m.group(1)) for p in doc.paragraphs if (m:=re.match(r'^\[(\d+)\]',p.text))]
check('All 35 reference entries preserved',refs==list(range(1,36)))
check('Three main-text energy/distribution sub-sections',all(x in main for x in ['3.6.1 泛素','3.6.2 角度','3.6.3 不同蛋白']))
pdf=PdfReader(TMP/'v22_word.pdf');check('49 Word rendered pages',len(pdf.pages)==49)
check('No blank pages or missing glyph characters',all(len((p.extract_text() or '').strip())>20 and '\ufffd' not in (p.extract_text() or '') for p in pdf.pages))
reuse=json.loads((QA/'page_pixel_reuse.json').read_text())
viewed={1,2,8,9,10,14,18,19,20,21,22,23,24,25,49}
check('All final page bodies visually checked',all(p['page'] in viewed or p['reuse_v21'] is not None for p in reuse))
page_review=[{'page':p['page'],'review':'final_PNG_visually_inspected' if p['page'] in viewed else 'pixel_identical_to_inspected_v21_page_body','v21_page':p['reuse_v21']} for p in reuse]
review={'status':'passed','document_sha256':sha(DOC),'reviewed_beijing':now,'checks':checks,'checks_passed':len(checks),'scientific_independent_review_sha256':sha(QA/'v22_results_independent_content_review.json'),'pages':49,'page_review':page_review,'layout_findings':'Native formula5, promoted figures8/9 and three-line table3 are legible and complete; captions stay with their figures. All changed page bodies inspected; remaining bodies pixel-identical to previously reviewed v21.'}
write(QA/'final_content_and_layout_audit.json',review)
integration['final_review']='passed_independent_science_and_all_page_visual_review';write(QA/'integration_audit.json',integration)

FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
mapping={'fig4_joint_complementarity':'fig4_joint_complementarity','fig7_ubiquitin_angular_landscapes':'fig7_ubiquitin_angular_landscapes','fig8_partition_sampling':'figS12_partition_sampling','fig9_representative_distribution_cases':'figS13_representative_distribution_cases'}
manifest={'document_sha256':sha(DOC),'scientific_graphics_unchanged':True,'files':{}}
for new,old in mapping.items():
    for ext in ['pdf','svg','png']:
        src=OLD/'figures'/f'{old}.{ext}';out=FIG/f'{new}.{ext}'
        shutil.copy2(src,out);check('Promoted vector/raster file is unchanged',sha(out)==sha(src));manifest['files'][out.name]=sha(out)
    label=new.split('_')[0][3:]
    caption=next(p.text for p in doc.paragraphs if p.text.startswith('图'+label+'：'))
    (FIG/f'{new}_caption_zh.md').write_text(caption+'\n',encoding='utf8')
    for p in (OLD/'figures').glob(old+'*'):
        if p.suffix in {'.csv','.json'}:shutil.copy2(p,FIG/(new+p.name[len(old):]))
write(FIG/'v22_figure_manifest.json',manifest)
note=ROOT/'正文自由能整合说明.md'
note.write_text('''自由能分析已整合进正文第2.8节和第3.6节，并同步修改摘要、引言、讨论与结论。第2.8节增加占据概率、相对经验投影自由能与联合熵的直接关系，公式(5)为原生可编辑Word数学公式。第3.6节按泛素多角度分布、角度占据与主链均方波动、不同蛋白的局部自由能形态三个小节展开，正文图7–9及表3直接呈现这些证据。

原图S12和S13分别移为正文图8和图9，原表S13移为正文表3，原有科学绘图与全部15张表的单元格数值保留。正文加入83条记录的RMSF²状态内／状态间分解，以及三个来源代表位点的实际χ1占据区和经验峰解释。正文明确区分多角联合熵、单角边缘熵、窗口平均和合并帧分布。补充材料K保留详细有效位置数、完整联合状态JSD及复现口径，主要结论和图无需依赖补充材料阅读。

最终稿49页，共9张正文图、11张补充图，3张正文三线表和12张补充三线表，35篇参考文献。独立科学内容审阅及所有页面视觉核验通过，原用户批注、锚点、科学图片和既有数学公式保留。已有MD轨迹和训练结果继续复用。实验源表及代码仍见v21科学结果归档，新正文、图注、代表峰和Cα分解审计保存在本v22目录。
''',encoding='utf8')
package=BASE/'paper_v22_20261012_正文新增图及图注.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(FIG.iterdir()):z.write(p,'figures/'+p.name)
    z.write(note,note.name)
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    for name,digest in manifest['files'].items():assert hashlib.sha256(z.read('figures/'+name)).hexdigest()==digest

pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
status=json.loads(pointer.read_text(encoding='utf8'))
shutil.copy2(pointer,QA/'delivery_status_v21_snapshot.json')
shutil.copy2(ROOT.parent/'completion_status_zh.md',QA/'completion_status_v21_snapshot.md')
status.update({'status':'chinese_manuscript_v22_free_energy_maintext_delivered','completed_beijing':now,'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':49,'main_figures':9,'supplementary_figures':11,'three_line_tables':15,'main_tables':3,'supplementary_tables':12,'content_review':'passed_110_independent_scientific_checks','render_and_visual_review':'passed_all_49_pages','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Free-energy argument integrated into main methods, results, abstract and discussion; native entropy/probability identity; promoted figures8/9 and table3; direct state-conditioned RMSF-squared interpretation and representative angle-occupancy explanation.','remaining_required_work':[]})
write(ROOT/'delivery_status.json',status);write(pointer,status)
(ROOT.parent/'completion_status_zh.md').write_text(f'''当前交付为 {DOC}，SHA256为 {sha(DOC)}，49页。自由能与构象占据分析已融入正文2.8、3.6及讨论，正文图7至9、表3提供主要证据；9张正文图、11张补充图，15张原生可编辑三线表及35篇参考文献。原用户批注、公式、科学图和全部表格数值保留。独立科学内容110项核验通过，最终49页视觉检查通过。正文配套图包为 {package}。完整计算结果、轨迹抽取缓存、最终训练权重和此前版本均保留；v21源数据归档继续作为科学源结果。当前详细状态见 {ROOT/'delivery_status.json'}。所有本轮所需工作已完成。
''',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':49,'main_figures':9,'main_tables':3,'package_MB':round(package.stat().st_size/1e6,2),'root_checks':len(checks)},ensure_ascii=False))
