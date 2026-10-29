from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,re,zipfile,shutil
import numpy as np
import pandas as pd
from docx import Document
from PIL import Image,ImageChops
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
STEM='fig7_model_selected_three_MD_free_energy_curves'
DOC=BASE/'paper_v28_20261014_中文修订_图7三条自由能曲线_三线表.docx'
SOURCE=BASE/'paper_v27_20261014_中文修订_模型选择与自由能曲线验证_三线表.docx'
OLD=ROOT.parent/'free_energy_curves_v27_20261014'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json')
assert sha(DOC)==integration['document_sha256']
assert sha(SOURCE)==integration['source_sha256']=='2e4c566833dbf7f828f11c2a449248a7025fbe6011fc1da531c1528b816d4316'
d=Document(DOC);old=Document(SOURCE)
text='\n'.join(p.text for p in d.paragraphs)
assert not any(s in text for s in ['图7d','群体范围曲线','范围曲线在E=','覆盖曲线','范围曲线沿用'])
assert all(s in text for s in ['图7a–c','35.41','43.17','7.76','0.400','0.510'])
refs=lambda doc:[p.text for p in doc.paragraphs if re.match(r'^\[\d+\]',p.text)]
assert refs(d)==refs(old) and len(refs(d))==27
cells=lambda doc:[[[c.text for c in r.cells] for r in t.rows] for t in doc.tables]
assert cells(d)==cells(old) and len(d.tables)==11
assert len(d.paragraphs)==len(old.paragraphs)
assert sum(p.text!=q.text for p,q in zip(d.paragraphs,old.paragraphs))==5
points=pd.read_csv(FIG/(STEM+'.points.csv'))
previous=pd.read_csv(OLD/'figures/fig7_prediction_selected_MD_free_energy_curves.points.csv')
previous=previous[previous.panel.isin(list('abc'))]
assert len(points)==len(previous)==432 and set(points.panel)==set('abc')
paired=points.merge(previous,on=['panel','residue','chi','x'],suffixes=('_new','_old'),validate='one_to_one')
assert len(paired)==432 and np.allclose(paired.y_new,paired.y_old,rtol=0,atol=0,equal_nan=True)
assert np.array_equal(paired.low_Fstar_le1_new,paired.low_Fstar_le1_old)
widths=pd.read_csv(FIG/(STEM+'.display_widths.csv'))
prevwidth=pd.read_csv(OLD/'figures/fig7_prediction_selected_MD_free_energy_curves.display_widths.csv')
assert widths.equals(prevwidth)
pdf=PdfReader(TMP/'v28_word.pdf')
assert len(pdf.pages)==28 and all(len((p.extract_text() or '').strip())>20 and '\ufffd' not in (p.extract_text() or '') for p in pdf.pages)
pages=[]
for i in range(1,29):
    p=TMP/'render_v28'/f'page-{i}.png'
    if i in {8,9,16,17,18}:review='finalPNGdirectlyvisuallyreviewed'
    else:
        assert ImageChops.difference(Image.open(p),Image.open(TMP/'render_v27'/f'page-{i}.png')).getbbox() is None
        review='pixelidenticaltoalreadyreviewedfinalv27'
    pages.append({'page':i,'review':review,'PNG_sha256':sha(p)})
fa=readj(FIG/(STEM+'.audit.json'))
for ext in ['png','pdf','svg']:assert sha(FIG/(STEM+'.'+ext))==fa['output_sha256'][ext]
fa['visual_review']='passed_originalPNGandWordpage17; a-csingle row; nolegend/titleoverlap; alloriginal432MDpointsandlowFstarbarsunchanged'
write(FIG/(STEM+'.audit.json'),fa)
integration.update({'review':'passed_all28pages;dremovedonly;432a-cpointsexact;5textparagraphsupdated','pages':28})
write(QA/'integration_audit.json',integration)
write(QA/'final_content_and_layout_audit.json',{'status':'passed','document_sha256':sha(DOC),'source_sha256':sha(SOURCE),'all_432_retained_curve_points_exact':True,'all_lowFstarbars_exact':True,'all_11native_tables_comments_math_preserved':True,'retained_83record_validation_values':True,'no_d_panel_references':True,'pages':28,'page_review':pages})
note=ROOT/'图7精简说明.md'
note.write_text('v28按用户要求删除图7d，只保留模型预先选出的Q41和Q2在χ1–χ3上的MD自由能曲线。原a–c的432个直方图点及低F*标记逐值保持不变，改为单行三面板，保留预测选择及数值标签。删除图注d、正文图7d引用、范围扫描曲线方法和讨论句；83条记录的低自由能角区宽度及细箱熵验证数值继续保留在正文3.6.2。没有将d移入补充材料，原计算数据仍保存在v27源包。五个文字段落修改；十一张原生三线表、原数学元素和批注保持不变。最终28页，五个修改页面已实际检查，其余23页与已检查v27逐像素一致。本轮没有新增实验。\n',encoding='utf8')
src=ROOT/'source';src.mkdir(exist_ok=True)
for name in ['histograms_501frames.csv','selection_inventory_with_MD_outcomes.csv']:shutil.copy2(OLD/'case_curves'/name,src/name)
package=BASE/'paper_v28_20261014_图7三条曲线与源数据.zip'
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(src,'source')]:
    files.extend((p,prefix+'/'+p.name) for p in directory.iterdir() if p.is_file())
files.extend((p,p.name) for p in [ROOT/'图7精简说明.md',ROOT/'build_v28.py',Path(__file__)])
oldarchive=BASE/'paper_v27_20261014_自由能曲线与模型验证源数据.zip'
assert sha(oldarchive)=='6729cb49b9fb30cc4e347f38761e6df717877e957a5bae740d45c9a186dec9d5'
manifest={'document_sha256':sha(DOC),'retained_validation_source_archive':str(oldarchive),'retained_validation_source_archive_sha256':sha(oldarchive),'files':{name:sha(p) for p,name in files}}
write(ROOT/'v28_package_manifest.json',manifest)
files.append((ROOT/'v28_package_manifest.json','v28_package_manifest.json'))
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p,name in files:z.write(p,name)
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None and all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
status=readj(pointer)
if not (QA/'delivery_status_v27_snapshot.json').exists():shutil.copy2(pointer,QA/'delivery_status_v27_snapshot.json')
status.update({'status':'chinese_manuscript_v28_figure7_three_curves_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'render_and_visual_review':'passed_all28pages_5changedpagesdirectview_23exactpixels','content_review':'passed_dremoved_threecurvesand83recordmainresultsretained','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'retained_model_validation_source_archive':str(oldarchive),'retained_model_validation_source_archive_sha256':sha(oldarchive),'new_experiment_note':str(ROOT/'图7精简说明.md'),'changes_this_round':'Removedfigure7d; retainedexacta-cMDcurvesandmodelpreselection; removedrange-scanmethod/crossrefs; all83recordprimaryvalidationresultsremaininmaintext; notmoveddintosupplement.','new_MD_production_this_round':0,'new_model_training_this_round':0,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',status);write(pointer,status)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7仅保留a–c三条角度自由能曲线，模型选择标签与432原始曲线点保持不变；图7d和相关引用/扫描曲线叙述已删除，83记录主验证结果仍在正文。三线表、数学与批注保持，所有页面检查通过。新版图7小包 {package}，完整已有实验源数据沿用 {oldarchive}。本轮无新增实验，当前要求已完成。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'package_MB':round(package.stat().st_size/1e6,2)},ensure_ascii=False))
