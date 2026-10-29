from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,shutil,zipfile

ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
DOC=BASE/'paper_v42_20261017_中文修订_蛋白背景与实例动机_三线表.docx'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\render_v42')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
digest=sha(DOC);assert digest=='0d56dd11e3abab0b7e12bbd75ef33c352d949bc305212c01f9eb5509541df887'
audit=read(QA/'v42_package_integrity_audit.json')
assert audit['document_sha256']==digest and audit['all_package_checks_pass']
assert len(read(QA/'v42_page_map.json'))==29
records=[]
for filename,first,last,reviewer in [
 ('visual_pages_1_10.json',1,10,'kras_figure6_build'),
 ('visual_pages_11_20.json',11,20,'kras_manuscript_text_review'),
 ('visual_pages_21_29.json',21,29,'kras_integration_audit')]:
    r=read(QA/filename)
    assert r['status'] in ('pass','passed') and r['document_sha256']==digest
    for n in range(first,last+1):
        records.append({'page':n,'status':'passed','reviewer':reviewer,'sha256':sha(TMP/f'page-{n}.png')})
write(QA/'final_page_visual_review.json',{'status':'passed','document_sha256':digest,'pages':records,
 'method':'All29 final150dpi original page PNGs individually opened by three assigned reviewers; root directly inspected14,16,17,18 and the previous layout iteration.',
 'required_fixes':[],'layout_adjustment':'Figure6 and its unchanged caption follow the KRAS biological setup before the quantitative paragraphs; no large gap and no split caption.'})
integration=read(QA/'integration_audit.json');assert integration['document_sha256']==digest
integration.update(status='complete',pages=29,render_status='passed',
 visual_review='passed_all29pages_individually_reviewed',remaining_required_work=[])
write(QA/'integration_audit.json',integration)
note=ROOT/'修改说明_zh.md'
note.write_text('''v42 蛋白功能背景与实例动机，2026-10-17

在结果3.4、3.5和3.6.1的开头补入功能背景、结合对象关系与模型检验动机。KIX介绍CBP的转录共激活作用及c-Myb/MLL的两个结合表面，接回四种伙伴组成中的L12活动区域定位；KRAS介绍GTP/GDP开关、G12D替换及细胞内蛋白来源肽经HLA呈递后被TCR识别的关系，接回同一工程TCR/HLA识别WT/G12D十肽的动态预测；泛素介绍76残基小蛋白的共价标记作用、EDD的UBA识别模块，接回同种侧链在不同位置的取向多样性比较。

背景事实复用已核对的原始研究[22]、[23]、[28]，对应DOI分别为10.1016/j.jmb.2005.09.059、10.1074/jbc.M705655200和10.1038/s41467-022-32811-1。各蛋白的来源与支持主张见kix_background.json、kras_background.json和ubiquitin_background.json。泛素核对还读取了作者机构保存的[23]正式全文。未把背景中的生物学机制归为本模型发现。

方法2.7仅补清KIX属于CBP及UBA属于EDD的关系，不重复扩写功能背景。正文保留全部定量结果，未增加实验、训练或推理。图6内容、尺寸和图注完全相同，仅移动到KRAS背景段后，使图文顺序连贯并减少空白。

与v41相比，DOCX包内仅word/document.xml改变；五个原段落修改、图6及其图注原样重排。全部8幅可见图及所有media文件、10张原生三线表、36项原生数学、28条参考文献、批注部件和锚点全部保持。29页已由独立隐藏Word实例只读导出，再用文档技能render_docx.py生成150dpi PNG，逐张视觉检查通过。
''',encoding='utf8')
(ROOT/'README_zh.md').write_text('''v42 背景修订材料

当前中文稿位于上两级，文件名paper_v42_20261017_中文修订_蛋白背景与实例动机_三线表.docx。本轮只补实例生物学背景和阅读衔接。

三个background.json保存各自原始来源、逐项支持与提案。build_v42.py为实际文稿编辑脚本；manuscript_QA包含最终包级检查、逐页视觉QA、正文提取和页码记录。export_v42_word.ps1和render_v42.py记录实际渲染方式。修改说明_zh.md记载最终采用范围；原提案中未采用的更详细机制/同义名不构成正文内容。

CURRENT_FIGURE_MAP.json沿用v41全部有效图源；没有新增或重绘图。KRAS图6及其实验源数据仍在kras_case_v41_20261017及原v41源数据包，泛素/自由能等原图源亦保持可复用。当前稿继承这些科学结果及配置，本次不新增计算实验。
''',encoding='utf8')
shutil.copy2(ROOT.parent/'kras_case_v41_20261017/CURRENT_FIGURE_MAP.json',ROOT/'CURRENT_FIGURE_MAP.json')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
snapshot=QA/'delivery_status_v41_snapshot.json'
if not snapshot.exists():shutil.copy2(pointer,snapshot)
prior=read(snapshot);assert prior['document_sha256']==integration['source_sha256']
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ('v42_edit_manifest.json','delivery_status.json'))
manifest={'document':str(DOC),'document_sha256':digest,'previous_document':integration['source'],
 'files':{p.relative_to(ROOT).as_posix():sha(p) for p in files}}
write(ROOT/'v42_edit_manifest.json',manifest)
package=BASE/'paper_v42_20261017_蛋白背景修订与文献核验.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    z.write(ROOT/'v42_edit_manifest.json','v42_edit_manifest.json')
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    for n,d in manifest['files'].items():assert hashlib.sha256(z.read(n)).hexdigest()==d
prior.update(status='chinese_manuscript_v42_protein_background_delivered',
 completed_beijing=datetime.now(timezone(timedelta(hours=8))).isoformat(),
 document=str(DOC),document_sha256=digest,pages=29,
 previous_document=integration['source'],previous_document_sha256=integration['source_sha256'],
 changes_this_round='Added biological background and case motivation for CBP KIX, KRAS G12D neoantigen and ubiquitin–EDD UBA; clarified two Methods object names; unchanged Figure6/caption repositioned after biological setup.',
 content_review='passed_primary_source_background_and_complete_package_preservation_audit',
 render_and_visual_review='passed_all29pages_individually_reviewed',
 background_revision_note=str(note),document_edit_package=str(package),document_edit_package_sha256=sha(package),
 background_references_reused=[22,23,28],
 new_MD_production_this_round=0,new_model_training_this_round=0,new_model_inference_this_round=0,new_MD_analysis_this_round=0,
 remaining_required_work=[])
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {digest}，29页。KIX、KRAS G12D及泛素–EDD UBA三处实例补充功能背景/伙伴关系/模型验证动机；复用原始参考22/23/28。所有图像、10个原生三线表、36项数学、28参考及原批注和锚点保持。图6原样调整位置。29页逐张检查通过。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':digest,'pages':29},ensure_ascii=False))
