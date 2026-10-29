from pathlib import Path
from datetime import datetime,timezone,timedelta
from collections import Counter
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
DOC=BASE/'paper_v27_20261014_中文修订_模型选择与自由能曲线验证_三线表.docx'
LAST=BASE/'paper_v25_20261013_中文修订_图7预测对照明确化_三线表.docx'
V26=ROOT.parent/'free_energy_model_v26_20261013'
CASE=ROOT/'case_curves';POP=ROOT/'population_curves'
STEM='fig7_prediction_selected_MD_free_energy_curves'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
now=datetime.now(timezone(timedelta(hours=8))).isoformat()
checks=[]
def check(name,value):
    checks.append({'name':name,'passed':bool(value)})
    assert value,name
integration=readj(QA/'integration_audit.json')
check('Final Word matches builder certificate',sha(DOC)==integration['document_sha256'])
check('Last delivered v25 unchanged',sha(LAST)=='523290c8297ba2034cf249bd07c57256cdc1920cb0f9530b4acd21d0e6a9335b')
d=Document(DOC);old=Document(LAST)
check('All405native table cells preserved',len(d.tables)==11 and sum(len(r.cells) for t in d.tables for r in t.rows)==405 and all([[c.text for c in r.cells] for r in a.rows]==[[c.text for c in r.cells] for r in b.rows] for a,b in zip(d.tables,old.tables)))
for i,t in enumerate(d.tables):
    check(f'Native three-line table{i+1}',not t._tbl.xpath('.//w:shd | ./w:tblPr/w:tblBorders') and all({etree.QName(e).localname for e in c._tc.tcPr.find(qn('w:tcBorders')) if e.get(qn('w:val')) not in {'nil','none'}}==({'top','bottom'} if j==0 else {'bottom'} if j==len(t.rows)-1 else set()) for j,r in enumerate(t.rows) for c in r.cells))
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(LAST) as a,zipfile.ZipFile(DOC) as b:
    ax,bx=etree.fromstring(a.read('word/document.xml')),etree.fromstring(b.read('word/document.xml'))
    anchors=lambda x:Counter((e.tag,e.get(qn('w:id'))) for e in x.iter() if e.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    math=lambda x:[etree.tostring(e,method='c14n') for e in x.xpath('.//m:oMath',namespaces=ns)]
    check('Original user comments byte-exact',a.read('word/comments.xml')==b.read('word/comments.xml'))
    check('Allcomment anchors and36 native math preserved',anchors(ax)==anchors(bx) and math(ax)==math(bx) and len(math(bx))==36)
    check('8active images',len(bx.xpath('.//a:blip',namespaces=ns))==8)
text='\n'.join(p.text for p in d.paragraphs)
figs=[m.group(1) for p in d.paragraphs if (m:=re.match(r'^图\s*(S?\d+)\s*[:：]',p.text))]
tables=[m.group(1) for p in d.paragraphs if (m:=re.match(r'^表\s*(S?\d+)\s*[:：]',p.text))]
check('Figure captions1to7andS1',figs==[str(i) for i in range(1,8)]+['S1'])
check('Table captions1to2andS1toS9',tables==['1','2']+['S'+str(i) for i in range(1,10)])
for prefix,cap in [('图',figs),('表',tables)]:check('Crossreferences '+prefix+' resolve',set(re.findall(prefix+r'\s*(S?\d+)',text))<=set(cap))
refs=[int(m.group(1)) for p in d.paragraphs if (m:=re.match(r'^\[(\d+)\]',p.text))]
cited=set()
for p in d.paragraphs:
    if re.match(r'^\[\d+\]',p.text):continue
    for m in re.finditer(r'\[([\d,]+)\]',p.text):cited.update(int(n) for n in m.group(1).split(','))
check('All27references cited',refs==list(range(1,28)) and cited==set(refs))
check('Seven supplementary modules', [m.group(1) for p in d.paragraphs if (m:=re.match(r'^([A-G])\s',p.text))]==list('ABCDEFG'))
check('Nooldfigure8crossrefs', '图8' not in text)
check('Originalmaterialweakresults retained',all(s in text for s in ['15.8%','16.5%','六对状态差值','0.474','0.471','0.537','0.568','0.577']))
check('Case andcohort evidenceroles clear',all(s in text for s in ['模型预测确定曲线验证位点','图7a–c','图7d','仅按冻结','原点0°的5°箱','逐点95%','2.990','64.17','27.50','35.41','43.17','7.76','7.26–8.26']))
check('Nooverclaim ofRMFS independent added information', '主链幅度之外的实际占据差异' not in text)

widthaudit=readj(V26/'experiment/independent_verification.json')
curveaudit=readj(POP/'second_implementation_verification.json')
caseaudit=readj(CASE/'independent_csv_verification.json')
check('Independent all203width reconstruction passed',widthaudit['status']=='passed' and widthaudit['trajectories']==203)
check('Independent83rawcachecurve reconstruction passed',curveaudit['status']=='passed' and curveaudit['checks']==10452 and curveaudit['raw_caches_reconstructed']==83)
check('Modelselectedcase independentCSV passed',caseaudit['status']=='passed')
check('Frozen originalprimaryindices curves verified',readj(POP/'frozen_primary_indices_verification.json')['status']=='passed')
inventory=pd.read_csv(CASE/'selection_inventory_with_MD_outcomes.csv').sort_values(['predicted_joint_Hchi_nats','residue_index'])
check('All6candidate predictions retained; modelargminargmax selected',len(inventory)==6 and inventory.iloc[0].residue=='Q41' and inventory.iloc[-1].residue=='Q2' and inventory.iloc[0].selection=='lowest_predicted_Hchi' and inventory.iloc[-1].selection=='highest_predicted_Hchi')
check('MDwiderQ49 not postselected',inventory.set_index('residue').loc['Q49','MD_width_501frames_deg']>inventory.iloc[-1].MD_width_501frames_deg and inventory.set_index('residue').loc['Q49','selection']=='not_selected')
for site,expected in [('Q41','27.50'),('Q2','64.17')]:check('Casefouroriginmean '+site,f'{inventory.set_index("residue").loc[site,"MD_width_501frames_deg"]:.2f}'==expected)
points=pd.read_csv(FIG/(STEM+'.points.csv'))
hist=pd.read_csv(CASE/'histograms_501frames.csv')
curves=pd.read_csv(POP/'population_curve_summary_frozen_primary_indices.csv')
check('All494curvepoints retained',len(points)==494 and points.groupby('panel').size().to_dict()=={'a':144,'b':144,'c':144,'d':62})
for qi,panel in enumerate('abc',1):
    plotted=points[points.panel==panel]
    source=hist[hist.chi=='chi'+str(qi)]
    merged=plotted.merge(source,left_on=['residue','x'],right_on=['residue','angle_center_degrees'],validate='one_to_one')
    check('Casepanel '+panel+' all144actualpoints exact',len(merged)==144 and np.allclose(merged.y,merged.F_star,rtol=0,atol=1e-14,equal_nan=True))
    check('Casepanel '+panel+' lowFstarbars exact',np.array_equal(merged.low_Fstar_le1_x,np.isfinite(merged.F_star)&(merged.F_star<=1)))
for group in ['low','high']:
    p=points[(points.panel=='d')&(points.group==group)].sort_values('x')
    check('Populationcurve '+group+' all31points exact',len(p)==31 and np.allclose(p.y,curves[group+'_mean'],rtol=0,atol=1e-13) and np.allclose(p.CI_low,curves[group+'_bootstrap_95_low'],rtol=0,atol=1e-13) and np.allclose(p.CI_high,curves[group+'_bootstrap_95_high'],rtol=0,atol=1e-13))
    check('Populationcurve '+group+' monotonic',np.diff(p.y).min()>=-1e-12)
e1=curves[np.isclose(curves.E,1)].iloc[0]
check('Primarycutoff fullmean andold95interval exact',abs(e1.difference_mean-7.757741521669362)<1e-14 and abs(e1.difference_bootstrap_95_low-7.26112475437627)<1e-14 and abs(e1.difference_bootstrap_95_high-8.26265812752495)<1e-14)
check('E1all83recordsdifferencepositive',int(e1.positive_difference_records)==83)
print('Numerical, content and source checks passed',len(checks),flush=True)

pdf=PdfReader(TMP/'v27_word.pdf')
check('28renderedpagesnoblankreplacementglyph',len(pdf.pages)==28 and all(len((p.extract_text() or '').strip())>20 and '\ufffd' not in (p.extract_text() or '') for p in pdf.pages))
first=readj(QA/'first_reviewed_v27_sha256.json');pixelmap=readj(QA/'page_pixel_comparison.json')
pages=[]
for i in range(1,29):
    png=TMP/'render_v27'/f'page-{i}.png'
    if i==17:method='FinalPNGdirectlyreviewed_aftercaptionoriginclarification'
    else:
        check(f'Finalpage{i}sameasreviewedfirstv27',sha(png)==first[str(i)])
        match=next(r for r in pixelmap if r['page']==i)['same_reviewed_v26_body_page']
        method='Firstv27directlyreviewed_andfinalPNGbyteidentical' if match is None else f'Bodypixelidenticaltoalreadyreviewedv26page{match};finalPNGbyteidentical'
    pages.append({'page':i,'review':method,'PNG_sha256':sha(png)})
figaudit=readj(FIG/(STEM+'.audit.json'))
for ext in ['png','pdf','svg']:check('Finalfigure hash '+ext,sha(FIG/(STEM+'.'+ext))==figaudit['output_sha256'][ext])
figaudit['visual_review']='passed_originalPNG_andfinalWordpage17; actualMDcurves/lowFstarbars/modelselection/83recordW(E)curves readable; no label orlegend overlap; all494points present'
figaudit['scientific_review']='passed_independentagent; mandatoryorigin0labeladded; actualcurvesareMD; modelpreselectsQ41/Q2fromall6GLN; Q49 widerMD not substituted; pointwise95andfouroriginaverage correctly distinguished'
write(FIG/(STEM+'.audit.json'),figaudit)
integration.update({'review':'passed_modelselection_to_MDcurves_to83recordW(E);all28pages','pages':28})
write(QA/'integration_audit.json',integration)
write(QA/'final_content_and_layout_audit.json',{'status':'passed','document_sha256':sha(DOC),'reviewed_beijing':now,'checks_passed':len(checks),'checks':checks,'pages':28,'page_review':pages,'independent_curve_checks':curveaudit['checks'],'independent_width_check_counts':widthaudit['check_counts'],'review_notes':'Allchangedv26/v27PNGsactuallyviewed; unmodifiedpagebodiespixel-identicaltoearlierreviewedversions. Finalpage17viewedagain; other27 finalPNGsbyteidenticaltofirstreviewedv27. MDcurvesnotmodelcurvepredictions; modeldrivesselectionandcohortgroups.'})

note=ROOT/'模型选择与自由能曲线验证说明.md'
note.write_text('''新版保留角度自由能曲线，并将模型的角色改为可核验的预测选择：冻结第0折模型先在泛素六个有效谷氨酰胺位点中按Hχ排序，预先选出最低Q41与最高Q2，然后才提取MD角度并画曲线。全部候选预测与MD结果保存，Q49实际宽度更大仍未替换为展示位点。a–c展示所选位置的三个χ角MD F*，底部色条直接表示原点0°箱格F*≤1的实际合格角区。正文四原点平均的跨角宽度为Q41 27.50°/角、Q2 64.17°/角。模型预测为2.000和2.990 nats，原五窗口MD联合熵为1.890和3.009 nats。

d展示原83记录、33,142个同类型有效位置中，完全由折外预测形成高低半组后的实际MD角区总宽度W(E)。E=0至3每0.1一点，全部31阈值保留。E=1时低组35.412503、高组43.170245°/角，差7.757742°/角、95%区间7.261125–8.262658，83记录均为正。曲线区间复用原主统计已冻结的seed20261013配对记录抽样索引，10,000次；阴影是逐点区间。原新增seed20261014范围曲线方案与结果另存，展示索引统一有冻结addendum与审计。

主低自由能宽度实验对全部203缓存独立复算39,993项检查、7,803,539个数组元素；范围曲线以另一实现重建83原缓存与33,493×31矩阵，10,452项检查通过，展示用原索引930,000个抽样值再次核对。病例3024直方图箱、144个逐角原点宽度和30个原联合熵窗口独立CSV复算通过。Q41旧角度档与重新提取角度数值容差内一致；原结构、轨迹、模型权重、预测和标签保持不变。全部sourceSHA及采样定义、所有候选而非只两个展示位点随本包保留。

正文3.6.1连接全链成绩、原Q41/Q62模型比较和新模型预选曲线；3.6.2保留自由能范围及细箱熵的群体验证。完整模型与链内模型在新宽度指标的增量区间跨零，正文仅报告各自读数，没有宣称完整模型显著更好；侧链与主链的排序及分组比较据实际结果报告。曲线说明模型筛选到的位置，未声称模型预测完整F*形状。

最终Word28页，七张正文图、一张补充图、十一张原生可编辑三线表、27条参考文献、七个补充模块；405表单元格、用户批注及36原生数学元素保持不变。所有修改页面完成视觉检查，其余页面与已检查版本逐像素一致。原四统计散点图v26是本轮中间稿，没有替代曲线版作为最终交付。本轮未进行新MD生产、模型训练或新推理。

本包包含v26新增宽度实验、v27模型选择的病例角度/曲线、83记录范围曲线/逐点重采样、冻结计划、脚本、图7矢量文件和审计。原83记录细箱熵实验源数据继续使用paper_v24_20261013_新增验证图与实验数据.zip；原始MD轨迹和模型源归档保持原位置。图7脚本在原项目结构中直接运行，使用既有绘图runtime和Microsoft YaHei字体。
''',encoding='utf8')
package=BASE/'paper_v27_20261014_自由能曲线与模型验证源数据.zip'
files=[]
for folder,prefix,recursive in [(V26/'experiment','lowFstar_width_experiment',True),(CASE,'model_selected_case',False),(POP,'population_curves',False),(FIG,'figures',False),(QA,'QA',False)]:
    for p in folder.rglob('*') if recursive else folder.iterdir():
        if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc':files.append((p,prefix+'/'+str(p.relative_to(folder)).replace('\\','/')))
for p in [note,ROOT/'prepare_curve_figure.py',ROOT/'build_v27.py',Path(__file__)]:files.append((p,p.name))
oldpackage=BASE/'paper_v24_20261013_新增验证图与实验数据.zip'
check('Existing83recordentropyarchiveunchanged',sha(oldpackage)=='0e59866b710a1c902d5b9936435bd5afe3bcc6f8ee4eac2758cd3bb4ffba110f')
manifest={'document_sha256':sha(DOC),'last_delivered_document_sha256':sha(LAST),'existing_entropy_source_archive':str(oldpackage),'existing_entropy_source_archive_sha256':sha(oldpackage),'files':{name:sha(p) for p,name in files}}
write(ROOT/'v27_package_manifest.json',manifest)
files.append((ROOT/'v27_package_manifest.json','v27_package_manifest.json'))
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p,name in files:z.write(p,name)
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
status=readj(pointer)
if not (QA/'delivery_status_v25_snapshot.json').exists():
    shutil.copy2(pointer,QA/'delivery_status_v25_snapshot.json')
    shutil.copy2(ROOT.parent/'completion_status_zh.md',QA/'completion_status_v25_snapshot.md')
status.update({'status':'chinese_manuscript_v27_model_selected_free_energy_curves_delivered','completed_beijing':now,'previous_document':str(LAST),'previous_document_sha256':sha(LAST),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'figures':8,'main_figures':7,'supplementary_figures':1,'three_line_tables':11,'main_tables':2,'supplementary_tables':9,'references':27,'content_review':'passed_modelselection_MDcurves_83recordrangevalidation','render_and_visual_review':'passed_all28pages','model_occupancy_source_archive':str(oldpackage),'model_occupancy_source_archive_sha256':sha(oldpackage),'new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'new_free_energy_width_records':83,'new_free_energy_width_eligible_positions':33142,'new_curve_thresholds':31,'new_model_selected_case_sites':['Q41','Q2'],'new_numerical_review':'passed_independent203cachewidth39993checks_83cachescuv10452checks_caseCSV3024bins','changes_this_round':'Preservedangularfreeenergycurveswithmodel-onlypreselectedsameAAextrema; modeldrivesselectionthenMDevaluates; 83recordprediction-definedW(E)curves; maintextfreeenergywidthandfineentropyvalidation; originalscatter-onlyfigurewithdrawnfromfinalmanuscript.','new_MD_production_this_round':0,'new_model_training_this_round':0,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',status);write(pointer,status)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7保留自由能曲线，以原模型预测先从六个谷氨酰胺中选Q41/Q2，MD曲线验证所选取向；d为83记录由折外预测确定的高低组W(E)群体曲线。原宽度、范围曲线和病例已独立复算，405三线表单元格、数学与批注保持不变；全部28页完成视觉验证。完整新增图与源数据包 {package}；83记录细箱熵源数据继续使用 {oldpackage}。本轮没有新增MD或模型训练，当前要求已完成。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'root_checks':len(checks),'package_MB':round(package.stat().st_size/1e6,2),'package_sha256':sha(package)},ensure_ascii=False))
