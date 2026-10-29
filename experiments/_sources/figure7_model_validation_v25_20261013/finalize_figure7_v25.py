from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,re,shutil,zipfile
import numpy as np
import pandas as pd
from docx import Document
from docx.oxml.ns import qn
from lxml import etree
from PIL import Image,ImageChops
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
DOC=BASE/'paper_v25_20261013_中文修订_图7预测对照明确化_三线表.docx'
SOURCE=BASE/'paper_v24_20261013_中文修订_模型论证与构象占据验证_三线表.docx'
CASE=ROOT.parent/'energy_landscape_v21_20261012/case_analysis'
STEM='fig7_ubiquitin_prediction_and_MD_Fstar'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
now=datetime.now(timezone(timedelta(hours=8))).isoformat()
checks=[]
def check(name,value):
    checks.append({'name':name,'passed':bool(value)})
    assert value,name
integration=json.loads((QA/'integration_audit.json').read_text(encoding='utf8'))
check('Final Word matches frozen builder certificate',sha(DOC)==integration['document_sha256'])
check('Previous v24 unchanged',sha(SOURCE)==integration['source_sha256']=='de8e19c4b941dbf0217c2b539e731e5932858c64f1fb81f7c07153e1b45a242e')
d=Document(DOC);old=Document(SOURCE)
check('Only three text paragraphs modified',sum(a.text!=b.text for a,b in zip(d.paragraphs,old.paragraphs))==3 and len(d.paragraphs)==len(old.paragraphs))
check('All405native table cells preserved',sum(len(r.cells) for t in d.tables for r in t.rows)==405 and all([[c.text for c in r.cells] for r in a.rows]==[[c.text for c in r.cells] for r in b.rows] for a,b in zip(d.tables,old.tables)))
for i,t in enumerate(d.tables):
    check(f'Native three-line table {i+1}',not t._tbl.xpath('.//w:shd | ./w:tblPr/w:tblBorders') and all({etree.QName(e).localname for e in c._tc.tcPr.find(qn('w:tcBorders')) if e.get(qn('w:val')) not in {'nil','none'}}==({'top','bottom'} if j==0 else {'bottom'} if j==len(t.rows)-1 else set()) for j,r in enumerate(t.rows) for c in r.cells))
text='\n'.join(p.text for p in d.paragraphs)
figs=[m.group(1) for p in d.paragraphs if (m:=re.match(r'^图\s*(S?\d+)\s*[:：]',p.text))]
tables=[m.group(1) for p in d.paragraphs if (m:=re.match(r'^表\s*(S?\d+)\s*[:：]',p.text))]
check('Figure captions 1 to8 andS1',figs==[str(i) for i in range(1,9)]+['S1'])
check('Table captions 1 to2 andS1 toS9',tables==['1','2']+['S'+str(i) for i in range(1,10)])
for prefix,cap in [('图',figs),('表',tables)]:check('Crossreferences '+prefix+' resolve',set(re.findall(prefix+r'\s*(S?\d+)',text))<=set(cap))
refs=[int(m.group(1)) for p in d.paragraphs if (m:=re.match(r'^\[(\d+)\]',p.text))]
cited=set()
for p in d.paragraphs:
    if re.match(r'^\[\d+\]',p.text):continue
    for match in re.finditer(r'\[([\d,]+)\]',p.text):cited.update(int(n) for n in match.group(1).split(','))
check('All27numbered references cited',refs==list(range(1,28)) and cited==set(refs))
cap=next(p.text for p in d.paragraphs if p.text.startswith('图7：'))
check('New figure7roles are explicit',all(s in cap for s in ['灰色圆点为MD','蓝色方点为模型','（b–d）','黄色实线为MD Q41','紫色虚线为MD Q62','a检验联合熵的预测幅值与位点差异']))
check('Old nine-panel figure7references removed',not any(s in cap for s in ['（a–c）','（d–f）','（g–i）','图顶端']))
check('Prediction values/errors and differences inmaintext',all(s in text for s in ['2.000和2.684','1.890和2.631','0.110和0.053','0.741 nats','0.684 nats','图7a','图7b–d']))
case=pd.read_csv(CASE/'ubiquitin_bin_sensitivity.csv').query('phase==0').set_index('residue').loc[['Q41','Q62']]
points=pd.read_csv(FIG/(STEM+'.prediction_points.csv'))
stats=pd.read_csv(FIG/(STEM+'.statistics.csv')).set_index('residue')
check('Exactly four actual prediction/reference points',len(points)==4 and set(points.method)=={'MD 参考','模型预测'})
for site,row in case.iterrows():
    ref=float(row.frozen_joint_h);pred=float(row.pred_joint_h)
    check(site+' MD label is five-window mean',abs(ref-row.joint_window_mean_h)<1e-12)
    check(site+' points exact',abs(points.query('residue==@site and method=="MD 参考"').y.iloc[0]-ref)<1e-14 and abs(points.query('residue==@site and method=="模型预测"').y.iloc[0]-pred)<1e-14)
    check(site+' absolute error exact',abs(stats.loc[site,'absolute_error_nats']-abs(pred-ref))<1e-14)
    check(site+' five windows1000frames',stats.loc[site,'MD_label_windows']==5 and stats.loc[site,'frames_per_window']==1000)
check('Exact roundedMDdelta',f'{case.frozen_joint_h.iloc[1]-case.frozen_joint_h.iloc[0]:.3f}'=='0.741')
check('Exact roundedmodeldelta',f'{case.pred_joint_h.iloc[1]-case.pred_joint_h.iloc[0]:.3f}'=='0.684')
curves=pd.read_csv(FIG/(STEM+'.all_MD_curve_points.csv'))
original=pd.read_csv(CASE/'ubiquitin_chi1_histograms_and_Fstar.csv')
check('All2592histogram rows saved unchanged',list(curves.columns)==list(original.columns) and len(curves)==len(original)==2592)
for column in curves:
    if original[column].dtype.kind in 'fiu':check('Curve column exact '+column,np.allclose(curves[column],original[column],rtol=0,atol=1e-14,equal_nan=True))
    else:check('Curve column exact '+column,curves[column].equals(original[column]))
for key,part in curves.groupby(['residue','window','chi']):
    p=part.probability.to_numpy();nonzero=p>0
    check('Fstar reconstruction '+str(key),np.allclose(part.F_star[nonzero],-np.log(p[nonzero]/p.max()),rtol=0,atol=1e-12) and part.F_star[~nonzero].isna().all())
pdf=PdfReader(TMP/'v25_word.pdf')
check('29rendered pages with no blank/replacement text',len(pdf.pages)==29 and all(len((p.extract_text() or '').strip())>20 and '\ufffd' not in (p.extract_text() or '') for p in pdf.pages))
page_review=[]
for i in range(1,30):
    current=TMP/'render_v25'/f'page-{i}.png'
    if i in {16,17,18,19}:method='final_PNG_directly_visually_inspected'
    else:
        previous=TMP/'render_v24'/f'page-{i}.png'
        check(f'Page{i}identicaltoalreadyreviewedv24',ImageChops.difference(Image.open(current),Image.open(previous)).getbbox() is None)
        method='pixel_identical_to_reviewed_v24_final_page'
    page_review.append({'page':i,'review':method,'PNG_sha256':sha(current)})
figaudit=json.loads((FIG/(STEM+'.audit.json')).read_text(encoding='utf8'))
for ext in ['png','pdf','svg']:check('Finalfigure output hash '+ext,sha(FIG/(STEM+'.'+ext))==figaudit['output_sha256'][ext])
figaudit['visual_review']='passed_original_PNG_and_final_Word_page17; allfourpoints/values/errors/differences visible; MDlegends placed in emptytopstrip above maximumFstar7.0003; nofooterxlabeloverlap; all2592curvepoints retained'
figaudit['independent_numeric_review']='passed_separate_agent_read_only_check_originalcacheandfivewindowmean; errors0.109592860/0.052805349; MDdifference0.740577956/modeldifference0.683790445'
write(FIG/(STEM+'.audit.json'),figaudit)
integration.update({'layout_review':'passed_all29pages; fourchangedpagesdirectlyviewed,other25pixelidenticaltoverifiedv24','pages':29})
write(QA/'integration_audit.json',integration)
write(QA/'final_content_and_layout_audit.json',{'status':'passed','document_sha256':sha(DOC),'reviewed_beijing':now,'checks_passed':len(checks),'checks':checks,'pages':29,'page_review':page_review,'scope':'Two predefinedsameAAcasepositions verify scalarprediction/reference; MDprojections explainactualangleoccupancy. Generalcohortargument staysinfigure8. No additionaltraining or MD production.'})

note=ROOT/'图7修改说明.md'
note.write_text('''图7原版本展示MD角度时间序列、密度和投影自由能，模型及参考联合熵只写在图顶端。新版将直接预测检验画成a面板：灰圆点为MD参考，蓝方点为冻结第0折模型联合Hχ，逐点标注读数，标注两位点误差及Q62−Q41差值。Q41参考／预测为1.890488441／2.000081301 nats，Q62为2.631066397／2.683871746 nats；绝对误差为0.110和0.053，参考／预测差值为0.741／0.684 nats。

b–d保留同一MD轨迹χ1、χ2、χ3的F*=−ln[p̂/p̂max]，曲线、坐标及图例均标明MD。黄色实线代表Q41，紫色虚线代表Q62；粗线使用5000帧，淡线完整保留五个连续1000帧窗口。全部2592个原始直方图记录保持不变，零计数断线。参考联合熵仍来自原120°联合分区的五窗口均值，经验F*采用5°逐角箱。a说明模型保留两位点数值和高低差异，b–d解释这些位点实际取向范围；图8沿用83记录的群体验证。

仅修改3.6.1的两个段落、图7图片与图注。其余八张图、十一张原生可编辑三线表的405单元格、数学元素、用户批注及27条参考文献保持不变。最终Word仍为29页，修改影响第16–19页，已逐页检查，其余25页与已检查v24逐像素一致。本轮没有追加模型训练或MD生产。

本小包保存新版图7的PNG/PDF/SVG、完整实际作图点、输入源CSV、脚本和审计。83记录实验完整源数据继续保存在paper_v24_20261013_新增验证图与实验数据.zip；已有训练与MD源归档保持原位置。plot_prediction_and_MD_Fstar.py在原项目位置运行，读取上级energy_landscape_v21_20261012/case_analysis的两个输入CSV（也在source目录保存）。生成图片使用numpy/pandas/matplotlib与Microsoft YaHei字体。
''',encoding='utf8')
srcdir=ROOT/'source';srcdir.mkdir(exist_ok=True)
for name in ['ubiquitin_bin_sensitivity.csv','ubiquitin_chi1_histograms_and_Fstar.csv']:shutil.copy2(CASE/name,srcdir/name)
package=BASE/'paper_v25_20261013_图7预测对照图与源数据.zip'
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(srcdir,'source')]:
    for p in directory.iterdir():
        if p.is_file():files.append((p,prefix+'/'+p.name))
for p in [note,ROOT/'prepare_figure7.py',ROOT/'build_figure7_v25.py',Path(__file__)]:files.append((p,p.name))
previous_package=BASE/'paper_v24_20261013_新增验证图与实验数据.zip'
check('Frozen83record experiment source archive unchanged',sha(previous_package)=='0e59866b710a1c902d5b9936435bd5afe3bcc6f8ee4eac2758cd3bb4ffba110f')
manifest={'document':str(DOC),'document_sha256':sha(DOC),'source_document_sha256':sha(SOURCE),'unchanged_model_occupancy_archive':str(previous_package),'unchanged_model_occupancy_archive_sha256':sha(previous_package),'files':{name:sha(p) for p,name in files}}
write(ROOT/'v25_package_manifest.json',manifest)
files.append((ROOT/'v25_package_manifest.json','v25_package_manifest.json'))
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p,name in files:z.write(p,name)
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
status=json.loads(pointer.read_text(encoding='utf8'))
if not (QA/'delivery_status_v24_snapshot.json').exists():
    shutil.copy2(pointer,QA/'delivery_status_v24_snapshot.json')
    shutil.copy2(ROOT.parent/'completion_status_zh.md',QA/'completion_status_v24_snapshot.md')
status.update({'status':'chinese_manuscript_v25_figure7_prediction_reference_explicit_delivered','completed_beijing':now,'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':29,'content_review':'passed_explicit_prediction_reference_and_MDinterpretation','render_and_visual_review':'passed_all29pages_fourchangedpagesdirectlyviewed_other25exactpixels','model_occupancy_source_archive':str(previous_package),'model_occupancy_source_archive_sha256':sha(previous_package),'new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Figure7 panela plots actualmodel/referenceHchi values, errorsandQ62minusQ41difference; panelsb-dretainonlyMD Fstar; 3.6.1andcaptionconnectthetwo; allothermodelresultsunchanged.','new_MD_production_this_round':0,'new_model_training_this_round':0,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',status);write(pointer,status)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付为 {DOC}，SHA256 {sha(DOC)}，29页。图7a现在直接展示联合Hχ预测与MD参考、误差和位点差值；b–d为同一位点的MD自由能投影，3.6.1及图注明确两类面板的证据作用。其余八张图、十一张原生三线表405单元格、数学元素、批注与27条参考文献保持不变。第16–19页视觉验证通过，其余25页与已验证v24逐像素一致。本轮新增图7小包为 {package}；83记录完整实验源数据继续使用 {previous_package}。没有追加训练或MD。当前要求已完成。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':29,'checks':len(checks),'package_MB':round(package.stat().st_size/1e6,2),'package_sha256':sha(package)},ensure_ascii=False))
