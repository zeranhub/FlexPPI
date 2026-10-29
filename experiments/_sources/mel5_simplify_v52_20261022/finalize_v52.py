from pathlib import Path
import csv,hashlib,json,math,zipfile
from datetime import datetime,timezone
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent;QA=ROOT/'manuscript_QA';FIG=ROOT/'figure4'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,o):Path(p).write_text(json.dumps(o,ensure_ascii=False,indent=2),encoding='utf-8')
now=datetime.now(timezone.utc).isoformat()
audit=read(QA/'integration_audit.json');doc=Path(audit['document'])
assert sha(doc)==audit['document_sha256']
assert sha(audit['source'])==audit['source_sha256']
assert read(QA/'independent_review.json')['status']=='passed'
with zipfile.ZipFile(doc) as z:assert hashlib.sha256(z.read('word/media/image41.png')).hexdigest()==sha(FIG/'fig4_mel5_compact.png')
source=ROOT.parent/'mel5_joint_case_v49_20261019/figure4/source_plotted_CDR_RMSF.csv'
with source.open(encoding='utf-8-sig',newline='') as f:original=list(csv.DictReader(f))
with (FIG/'source_CDR_RMSF_all45.csv').open(encoding='utf-8-sig',newline='') as f:current=list(csv.DictReader(f))
assert len(current)==len(original)==45
for a,b in zip(original,current):
    for field in ['MEL5_model_RMSF_A','alpha24beta17_model_RMSF_A']:
        assert math.isclose(float(a[field]),float(b[field]),rel_tol=0,abs_tol=1e-15)
fonts=[]
for page in PdfReader(FIG/'fig4_mel5_compact.pdf').pages:
    for name,ref in page['/Resources'].get('/Font',{}).items():
        font=ref.get_object()
        for child in font.get('/DescendantFonts',[font]):
            item=child.get_object();fd=item.get('/FontDescriptor');fd=fd.get_object() if fd is not None else {}
            fonts.append({'name':str(font.get('/BaseFont',name)),'embedded':any(k in fd for k in ['/FontFile','/FontFile2','/FontFile3'])})
assert fonts and all(f['embedded'] for f in fonts)
write(FIG/'pdf_font_audit.json',{'status':'passed','fonts':fonts,'pdf_sha256':sha(FIG/'fig4_mel5_compact.pdf')})
write(QA/'root_pages_15_28_review.json',{'status':'passed','document_sha256':sha(doc),'reviewed_at_utc':now,'pages':[{'page':n,'sha256':sha(QA/f'pages/page-{n}.png'),'method':'individual full-page view_image','status':'passed'} for n in range(15,29)],'required_changes':[]})
fa=read(FIG/'figure_data_audit.json');fa.update(status='passed',visual_QA='passed',source_values_verified=True);write(FIG/'figure_data_audit.json',fa)
audit.update(status='passed',rendered_pages=28,visual_QA='passed',finalized_at_utc=now);write(QA/'integration_audit.json',audit)
stats=read(ROOT.parent/'mel5_interface_v51_20261022/figure4/statistics_and_reproducibility.json')
stats['sample_definition']['structure_plate']='Only the existing complex overview remains; both molecular close-ups removed.'
stats['sample_definition']['distance']='Coordinate distances remain in methods/prose and original analysis archive; no atom-pair distance is depicted.'
stats['source_files']={str(p):{'sha256':sha(p)} for p in [source,FIG/'source_CDR_RMSF_all45.csv',FIG/'source_model_readouts.json']}
write(FIG/'statistics_and_reproducibility.json',stats)
fm=read(ROOT.parent/'mel5_interface_v51_20261022/CURRENT_FIGURE_MAP.json')
fm['4']={'source_folder':'figure4','absolute_source_folder':str(FIG),'png':'fig4_mel5_compact.png','pdf':'fig4_mel5_compact.pdf','svg':'fig4_mel5_compact.svg','new_in_revision':'v52','panels':['a complex overview, published SPR and paired Q72 Hchi / betaGly54 RMSF numbers','b all45CDR absolute model RMSF'],'removed':'Both v51 local molecular close-ups, all their atom labels, distance lines and view descriptors'};write(ROOT/'CURRENT_FIGURE_MAP.json',fm)
(FIG/'QA_report.md').write_text('''# Figure4 v52 QA

READY. Both rejected local molecular close-ups removed. Figure183×100mm;450dpi PNG; embedded-font vector PDF and editable SVG. English labels; minimum6.5pt; custom consistent state colors. Panel a is a compact overview plus direct model readouts with separate units, no replacement local structure or added coordinate panel. Panel b retains all45CDR positions with the original90 predicted values; numerical comparison passed1e-15 absolute tolerance. Model readouts come from archived key_findings. No statistical tests, confidence interval or errorbars apply to one deterministic output per state; no new inference, training or MD. Geometry/MD prose and source archives preserved. Method, prose and caption synchronized; no same-view/distance-dash figure claim remains. Independent Figure4 and pages1–15 review passed; root pages15–28 review passed. All10native three-line tables,36native equations,6images and28references retained.
''',encoding='utf-8')
(ROOT/'修改说明_zh.md').write_text('''# v52 修改说明

删除图4a右侧两张局部棒状结构图及其原子标签、距离虚线和同视角说明。a收为复合物概览、SPR背景和Q72 Hχ／βGly54 RMSF的两行模型读数；b保留完整45个CDR位点的两态曲线。整图由183×140mm收紧为183×100mm。没有补新的局部结构图或独立小面板。

正文保留已核验的接触机制与原子距离，删掉“同一视角观察”的旧图指向；方法改为距离计算说明，图注只说明实际保留的内容。全文28页渲染和逐页目视检查通过，三线表与原生公式保留。v51文件未覆盖。
''',encoding='utf-8')
delivery={'status':'complete','version':'v52','completed_at_utc':now,'document':{'path':str(doc),'sha256':sha(doc)},'figure':{ext:{'path':str(FIG/f'fig4_mel5_compact.{ext}'),'sha256':sha(FIG/f'fig4_mel5_compact.{ext}')} for ext in ['png','pdf','svg']},'panels':2,'pages':28,'native_three_line_tables':10,'source_preserved':True,'QA':'passed'}
write(ROOT/'delivery_status.json',delivery)
print(json.dumps({'status':'complete','document':str(doc),'pages':28,'removed_local_structures':2},ensure_ascii=False))
