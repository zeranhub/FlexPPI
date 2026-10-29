"""Final content, numerical, figure and rendered-page gate; package delivery."""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone,timedelta
import hashlib,json,re,shutil,zipfile
import numpy as np
import pandas as pd
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from PIL import Image,ImageChops
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
EXP=ROOT/'experiment_predictive_occupancy';FIG=ROOT/'figures'
SOURCE=BASE/'paper_v23_20261012_中文修订_论证精简与补充整合_三线表.docx'
DOC=BASE/'paper_v24_20261013_中文修订_模型论证与构象占据验证_三线表.docx'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
now=datetime.now(timezone(timedelta(hours=8))).isoformat()
checks=[]
def check(name,cond):checks.append({'name':name,'passed':bool(cond)});assert cond,name
integration=json.loads((QA/'integration_audit.json').read_text(encoding='utf8'))
numerical=json.loads((QA/'predictive_occupancy_reconstruction_audit.json').read_text(encoding='utf8'))
check('Final document matches builder certificate',sha(DOC)==integration['document_sha256']=='de8e19c4b941dbf0217c2b539e731e5932858c64f1fb81f7c07153e1b45a242e')
check('Frozen v23 source unchanged',sha(SOURCE)==integration['source_sha256']=='04e952548a9e2de9091f138d510238f93e816f5731ad801bee9af5f07fcbf4f8')
check('Independent algorithm reconstructs all203 trajectories and both analyses',numerical['status']=='passed' and numerical['checks_passed']==2108 and numerical['plan_sha256']==sha(EXP/'analysis_plan.json'))
d,old=Document(DOC),Document(SOURCE)
cells=lambda t:tuple(tuple(c.text for c in row.cells) for row in t.rows)
check('Ten existing model-evidence tables unchanged',len(d.tables)==11 and all(cells(d.tables[i])==cells(old.tables[i]) for i in range(11) if i!=8))
check('363 existing native cells plus42 newly computed cells',sum(len(r.cells) for i,t in enumerate(d.tables) if i!=8 for r in t.rows)==363 and sum(len(r.cells) for r in d.tables[8].rows)==42)
tab=pd.read_csv(ROOT/'tableS7_predictive_occupancy.csv')
expect=[]
for _,r in tab.iterrows():expect.append(['全部' if r.scope=='all' else r.scope[1:]+'角',str(int(r.n_records)),f'{int(r.positions):,}',*[f'{r[k]:.3f}' for k in ['Full_Hchi_rho','Intra_Hchi_rho','Full_RMSF_rho','Full_Hchi_high_low']]])
check('Every new tableS7cell agrees with fullprecision certified CSV',cells(d.tables[8])[1:]==tuple(tuple(r) for r in expect))
for i,t in enumerate(d.tables):
    check(f'Table{i+1} three-line native formatting',not t._tbl.xpath('.//w:shd | ./w:tblPr/w:tblBorders') and all({etree.QName(e).localname for e in c._tc.tcPr.find(qn('w:tcBorders')) if e.get(qn('w:val')) not in {'nil','none'}}==({'top','bottom'} if j==0 else {'bottom'} if j==len(t.rows)-1 else set()) for j,r in enumerate(t.rows) for c in r.cells))
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(DOC) as b:
    ax,bx=etree.fromstring(a.read('word/document.xml')),etree.fromstring(b.read('word/document.xml'))
    check('Original user comments byte-identical',a.read('word/comments.xml')==b.read('word/comments.xml'))
    anchors=lambda x:Counter((e.tag,e.get(qn('w:id'))) for e in x.iter() if e.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    math=lambda x:[etree.tostring(e,method='c14n') for e in x.xpath('.//m:oMath',namespaces=ns)]
    check('All user comment anchors and36native math elements preserved',anchors(ax)==anchors(bx) and math(ax)==math(bx) and len(math(bx))==36)
    refs=etree.fromstring(b.read('word/_rels/document.xml.rels'));targets={x.get('Id'):x.get('Target') for x in refs}
    oldmedia={hashlib.sha256(a.read(n)).hexdigest() for n in a.namelist() if n.startswith('word/media/')}
    images=[targets[e.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed')] for e in bx.xpath('.//a:blip',namespaces=ns)]
    digests=[hashlib.sha256(b.read('word/'+n)).hexdigest() for n in images]
    figsha=sha(FIG/'fig8_model_angular_occupancy.png')
    check('8retained images exact and1newfigure8exact',len(digests)==9 and digests.count(figsha)==1 and all(s in oldmedia or s==figsha for s in digests))
text='\n'.join(p.text for p in d.paragraphs)
main=text[:text.index('4 讨论')]
figs=[m.group(1) for p in d.paragraphs if (m:=re.match(r'^图\s*(S?\d+)\s*[:：]',p.text))]
tables=[m.group(1) for p in d.paragraphs if (m:=re.match(r'^表\s*(S?\d+)\s*[:：]',p.text))]
check('Figure captions1to8andS1consecutive',figs==[str(i) for i in range(1,9)]+['S1'])
check('Table captions1to2andS1toS9consecutive',tables==['1','2']+['S'+str(i) for i in range(1,10)])
for prefix,cap in [('图',figs),('表',tables)]:check('Crossreferences'+prefix+'allresolve',set(re.findall(prefix+r'\s*(S?\d+)',text))<=set(cap))
references=[int(m.group(1)) for p in d.paragraphs if (m:=re.match(r'^\[(\d+)\]',p.text))]
check('27referencesconsecutive',references==list(range(1,28)))
cited=set()
for p in d.paragraphs:
    if re.match(r'^\[\d+\]',p.text):continue
    for match in re.finditer(r'\[([\d,]+)\]',p.text):cited.update(int(n) for n in match.group(1).split(','))
check('All27referencescitedwithnoorphans',cited==set(references))
check('Seven supplement modules', [m.group(1) for p in d.paragraphs if (m:=re.match(r'^([A-G])\s',p.text))]==list('ABCDEFG'))
check('Models evaluated in eachresultssection',all(s in main for s in ['3.1 两类动态预测','3.2 同类型残基','3.3 主链与侧链','3.4 ThermoPCD','3.5 KIX','3.6.1 泛素同类型','3.6.2 折外侧链']))
check('No removedpureMDsectionsorcrossrefs',not any(s in text for s in ['图9','3.6.3','半程JSD','189.2','分区敏感性','31,621个非Pro','事后分组校准']))
check('Predictive values alongsideMD interpretation',all(s in main for s in ['0.868','1.890','2.631','2.000','2.684','图7：','自由能']))
check('Actual material weaker results retained',all(s in text for s in ['15.8%','16.5%','六对状态差值','0.474','0.471','0.537','0.568','0.577']))
vals=json.loads((ROOT/'manuscript_values.json').read_text(encoding='utf8'))
for category,predictor,rounded in [('correlation','full_hchi','0.510'),('correlation','intrachain_hchi','0.493'),('correlation','full_rmsf','0.321'),('high_low','full_hchi','0.293'),('high_low','intrachain_hchi','0.281'),('high_low','full_rmsf','0.198')]:
    source=next(r for r in vals['replica1']['all'][category] if r['predictor']==predictor)
    check(f'Certifiedmainvalue{category}{predictor}',f'{source["mean"]:.3f}'==rounded and rounded in main)
high=next(r for r in vals['replica1']['all']['high_low'] if r['predictor']=='full_hchi')
check('83positive recordhighlowcontrasts; exactroundedinterval',high['positive_records']==83 and f'{high["bootstrap_95_low"]:.3f}'=='0.277' and f'{high["bootstrap_95_high"]:.3f}'=='0.307')
points=pd.read_csv(FIG/'fig8_model_angular_occupancy.points.csv')
check('Figure8all498pointsand83recordspermethod',len(points)==498 and points.groupby(['panel','predictor']).sample_id.nunique().eq(83).all())
for panel,filename,column in [('a','replica1_per_record_correlations.csv','rho_type_centered'),('b','replica1_high_low_per_record.csv','high_minus_low_MD')]:
    raw=pd.read_csv(EXP/filename);raw=raw[(raw.scope=='all')&(raw.metric=='fine_entropy_mean_per_chi')]
    merged=points[points.panel==panel].merge(raw,on=['sample_id','predictor'],validate='one_to_one')
    check(f'Figure8panel{panel}pointsnumericallyexact',len(merged)==249 and np.allclose(merged.y,merged[column],rtol=0,atol=1e-14))
figaudit=json.loads((FIG/'fig8_model_angular_occupancy.audit.json').read_text(encoding='utf8'))
figaudit['visual_review']='passed_actual_original_PNG_and_final_Word_page18; allpoints/pairedlines/meanintervals/titles/units/legend readable'
figaudit['statistics_and_reproducibility']='83records;31,?'.replace('31,?','33142sameAApositions; fixed OOF models; existing203caches; rep1primary; type-centeredSpearman; prediction-onlyAAhalfgroupselection; AAequal thenrecordequalcontrast; source-stratified10000bootstrap; nats/angle; no statisticalpvalues; no multiplecomparison significance claims')
write(FIG/'fig8_model_angular_occupancy.audit.json',figaudit)
pdf=PdfReader(TMP/'v24_word.pdf')
check('29renderedWordpages',len(pdf.pages)==29)
check('No blankpages orreplacementglyph',all(len((p.extract_text() or '').strip())>20 and '\ufffd' not in (p.extract_text() or '') for p in pdf.pages))
pages=[]
for i in range(1,30):
    current=TMP/'render_v24'/f'page-{i}.png';previous=TMP/'render_v24_first_visual'/f'page-{i}.png'
    if i in {10,11,16,27}:method='final_PNG_directly_visually_inspected'
    else:
        check(f'Finalpage{i}identicaltopreviouslyvisuallyinspectedv24',ImageChops.difference(Image.open(current),Image.open(previous)).getbbox() is None)
        method='pixel_identical_to_directly_inspected_first_v24_render'
    pages.append({'page':i,'review':method,'PNG_sha256':sha(current)})
write(QA/'final_content_and_layout_audit.json',{'status':'passed','document_sha256':sha(DOC),'reviewed_beijing':now,'checks_passed':len(checks),'checks':checks,'pages':29,'page_review':pages,'numerical_verification_sha256':sha(QA/'predictive_occupancy_reconstruction_audit.json'),'actual_review':'All29pagePNG bodies manually viewed; final4changedpages viewedagain; other25pixel-identical. Maintext read in full; original37?'.replace('37?','36math elements preserved; newexperimentalclaims checkedagainstcertifiedCSV'),'material_limit':'The model estimates scalarRMSF/Hchi; MDcurves explain predicted site differences; fineentropy cohort tests ranking/groupselection,not full density orkinetic prediction.'})
integration['review']='passed_reconstructed_science_model_argument_content_and_all29page_visualQA';integration['pages']=29
write(QA/'integration_audit.json',integration)
note=ROOT/'模型论证与新实验说明.md'
note.write_text('''v24 将全文组织为模型能力的证据链。整体对照验证静态几何对两项估计的作用，完整／链内及类型内比较检验邻域和位置差异，外部队列考察迁移，KIX验证活动区域定位。3.6.1将泛素的全链留出成绩、Q41/Q62预测与参考、真实多角占据及经验自由能投影连续呈现；3.6.2新增固定折外预测驱动的83记录群体验证。自由能保留正文，并用于解释已验证预测差异的实际取向形态。

新增实验沿用83条冻结记录和203个已有角度缓存，原预测、训练标签与权重保持不变。主指标预先固定为每角5°細箱占据熵，对四个1.25°错位原点和可定义角取均值，含Pro；全部33,493个共同位置中，同类型至少四位点保留33,142个。第一条轨迹主分析中，完整Hχ／链内Hχ／完整RMSF的类型去均值相关为0.510／0.493／0.321。完整Hχ预测高低半组真实指标差为0.293 nats/角，83记录均为正值；其95%来源分层重采样区间为0.277–0.307。角数分层、全部可用重复均值和完整配对统计随源数据提供。另一套向量化直方图实现复算全部203轨迹、每种读数相关及所有高低组，2,108项验证通过。新图8展示全部83记录，补充表S7提供完整角数分层。

原分区敏感性图8、条件Cα分解、JSD与参考驱动的三蛋白病例图9撤出稿件，继续在v21–v23源归档保存。删除重复的否定式防御句和没有论证用途的文献，保留固定标签、留出关系、量纲、采样和聚合等必要定义；弱组合召回、发布基线的实际成绩和KIX状态差值仍如实报告。保留的十张既有表共363单元格逐字不变，以新模型验证42单元格替换原纯分区表。

最终Word29页，八张正文图、一张补充图、十一张原生可编辑三线表、七个补充模块和27条重新编号的参考文献。全部用户批注及原生数学元素保持不变，所有页面完成视觉检查。新数据包包含实验计划、每位点和每记录统计、重采样结果、脚本、图8矢量文件和审计；大型原始角度缓存及既有模型权重沿用已保存源归档。本轮未追加MD生产或模型训练。
'''.replace('細箱','细箱'),encoding='utf8')
package=BASE/'paper_v24_20261013_新增验证图与实验数据.zip'
files=[]
for p in EXP.iterdir():
    if p.is_file() and p.suffix in {'.csv','.gz','.json','.py'}:files.append((p,'experiment/'+p.name))
for directory,prefix in [(FIG,'figures'),(QA,'QA')]:
    for p in directory.iterdir():
        if p.is_file():files.append((p,prefix+'/'+p.name))
for p in [note,ROOT/'manuscript_values.json',ROOT/'tableS7_predictive_occupancy.csv',ROOT/'verify_predictive_occupancy.py',ROOT/'prepare_model_figure.py',ROOT/'build_model_argument_v24.py',Path(__file__)]:files.append((p,p.name))
manifest={'document_sha256':sha(DOC),'frozen_source_document_sha256':sha(SOURCE),'note':'New derived experiment source tables and all resampling draws; original large angle caches referenced by frozen hashes and retained in previous scientific archive. Existing training andMD unchanged.','files':{name:sha(p) for p,name in files}}
write(ROOT/'v24_package_manifest.json',manifest)
files.append((ROOT/'v24_package_manifest.json','v24_package_manifest.json'))
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p,name in files:z.write(p,name)
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
status=json.loads(pointer.read_text(encoding='utf8'))
if not (QA/'delivery_status_v23_snapshot.json').exists():
    shutil.copy2(pointer,QA/'delivery_status_v23_snapshot.json');shutil.copy2(ROOT.parent/'completion_status_zh.md',QA/'completion_status_v23_snapshot.md')
status.update({'status':'chinese_manuscript_v24_model_argument_and_occupancy_validation_delivered','completed_beijing':now,'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':29,'figures':9,'main_figures':8,'supplementary_figures':1,'three_line_tables':11,'main_tables':2,'supplementary_tables':9,'supplementary_modules':7,'references':27,'content_review':'passed_model_argument_and_source_verified_claims','new_numerical_review':'passed_2108checks_separate_vectorized_reconstruction_all203cachedtrajectories','render_and_visual_review':'passed_all_29pages','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'new_predictive_occupancy_records':83,'new_predictive_occupancy_eligible_positions':33142,'changes_this_round':'Entirepaperorganizedaroundmodelclaims; ubiquitinpredictionandMDfreeenergyconnected; newfrozenOOFpredictiveangularoccupancycohortexperimentandfigure8; pureMDpartitionandrepresentativecasesremovedfrommanuscript; truthfulweakerresultsretained.','new_MD_production_this_round':0,'new_model_training_this_round':0,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',status);write(pointer,status)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付为 {DOC}，SHA256 {sha(DOC)}，29页。全文围绕模型能力重组，自由能保留正文并连接泛素预测与新83记录的预测驱动构象占据验证。八张正文图、一张补充图、十一张原生可编辑三线表、七个补充模块和27条重新编号引用；十张既有表数值、原用户批注及数学元素保持不变。新实验已以另一实现复算全部203缓存，2108项验证通过；全部29页视觉验证通过。完整新增图与实验源数据为 {package}，原MD与训练源归档继续保留。详细状态见 {ROOT / "delivery_status.json"}。当前要求全部完成。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':29,'main_figures':8,'tables':11,'root_checks':len(checks),'new_experiment_verifications':2108,'package_MB':round(package.stat().st_size/1e6,2),'package_sha256':sha(package)},ensure_ascii=False))
