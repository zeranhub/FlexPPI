from pathlib import Path
import csv, hashlib, json, math, zipfile
from datetime import datetime, timezone
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
QA, FIG = ROOT/'manuscript_QA', ROOT/'figure4'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p, data): Path(p).write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
now = datetime.now(timezone.utc).isoformat()
audit = read(QA/'integration_audit.json')
doc = Path(audit['document'])
assert sha(doc) == audit['document_sha256']
assert sha(audit['source']) == audit['source_sha256']
png = FIG/'fig4_mel5_contact_map.png'
assert sha(png) == audit['figure_sha256']
with zipfile.ZipFile(doc) as z:
    assert hashlib.sha256(z.read('word/media/image41.png')).hexdigest() == sha(png)
for filename, key in [('independent_visual_review.json','png_sha256'), ('independent_data_review.json','PNG_sha256')]:
    report = read(FIG/filename)
    assert report['status'] == 'passed' and report[key] == sha(png)
for filename in ['independent_pages_1_14_review.json', 'independent_pages_15_28_science_review.json']:
    report = read(QA/filename)
    assert report['status'] == 'passed' and report['document_sha256'] == sha(doc)
assert len(list((QA/'pages').glob('page-*.png'))) == 28
assert len(PdfReader(QA/'v53_word.pdf').pages) == 28
scope = read(ROOT/'analysis/data_scope_audit.json')
assert scope['status'] == 'passed'
for item in scope['source_sha256'].values():
    assert sha(item['path']) == item['sha256']
with (FIG/'source_all45_CDR_RMSF.csv').open(encoding='utf-8-sig', newline='') as f:
    current = list(csv.DictReader(f))
with Path(scope['source_sha256']['all45_plotted_CDR_RMSF']['path']).open(encoding='utf-8-sig', newline='') as f:
    original = list(csv.DictReader(f))
assert len(current) == len(original) == 45
for a, b in zip(original, current):
    for field in ['MEL5_model_RMSF_A','alpha24beta17_model_RMSF_A']:
        assert math.isclose(float(a[field]),float(b[field]),rel_tol=0,abs_tol=1e-15)
for name, count in [('source_all54_contacts.csv',54),('source_all6_local_Hchi.csv',6),('source_all9_local_RMSF.csv',9)]:
    with (FIG/name).open(encoding='utf-8-sig', newline='') as f:
        assert len(list(csv.DictReader(f))) == count
fonts = []
for page in PdfReader(FIG/'fig4_mel5_contact_map.pdf').pages:
    for name, ref in page['/Resources'].get('/Font',{}).items():
        font = ref.get_object()
        for child in font.get('/DescendantFonts',[font]):
            item = child.get_object()
            fd = item.get('/FontDescriptor')
            fd = fd.get_object() if fd is not None else {}
            fonts.append({'name':str(font.get('/BaseFont',name)), 'embedded':any(k in fd for k in ['/FontFile','/FontFile2','/FontFile3'])})
assert fonts and all(f['embedded'] for f in fonts)
write(FIG/'pdf_font_audit.json', {'status':'passed','fonts':fonts,'pdf_sha256':sha(FIG/'fig4_mel5_contact_map.pdf')})
fa = read(FIG/'figure_data_audit.json')
fa.update(status='passed',visual_QA='passed',source_values_verified=True,finalized_at_utc=now)
write(FIG/'figure_data_audit.json',fa)
audit.update(status='passed',rendered_pages=28,visual_QA='passed',finalized_at_utc=now)
write(QA/'integration_audit.json',audit)
stats = read(ROOT.parent/'mel5_simplify_v52_20261022/figure4/statistics_and_reproducibility.json')
stats['sample_definition'].update(
    structure_plate='One existing complex overview as a locator. Quantitative contact grid contains all54 pairs between six HLA rows and nine TCRbeta48-56 columns.',
    distance='Minimum distance over all nonhydrogen atoms in each residue pair in archived starting structures; <=4A contact, >4 to6A nearby, >6A far. Not fixed-atom distances or MD occupancy.',
    local_selection='Within HLA65-76, retain each site with valid Hchi and minimum distance to TCRbeta48-56 <=6A in either archived structure: R65 K66 K68 Q72 R75 V76. All54 combinations displayed.',
    predictions='One existing deterministic joint-model output per bound state; six row-aligned Hchi sites, nine column-aligned local RMSF sites, and all45 aligned CDR positions.')
stats['source_files'] = {str(p):{'sha256':sha(p)} for p in sorted(FIG.glob('source_*.csv'))}
stats['source_coordinate_audit'] = str(ROOT/'analysis/data_scope_audit.json')
stats['scope'] = 'Paired local structure contacts aligned with model main-chain and side-chain responses. No Hchi-to-RMSF causal direction or quantitative affinity prediction is inferred.'
write(FIG/'statistics_and_reproducibility.json',stats)
fm = read(ROOT.parent/'mel5_simplify_v52_20261022/CURRENT_FIGURE_MAP.json')
fm['4'] = {'source_folder':'figure4','absolute_source_folder':str(FIG),'png':'fig4_mel5_contact_map.png','pdf':'fig4_mel5_contact_map.pdf','svg':'fig4_mel5_contact_map.svg','new_in_revision':'v53','panels':['a all54 residue contact pairs, six row-aligned model Hchi predictions and nine column-aligned model RMSF predictions','b all45 CDR absolute model RMSF'],'removed':'v52 text and numeric summary cards; no separate c/d panels'}
write(ROOT/'CURRENT_FIGURE_MAP.json',fm)
(FIG/'QA_report.md').write_text('''# Figure4 v53 QA

READY. Actual 6 by9 contact grid with all54 residue pairs and108 state glyphs; six HLA Hchi rows and nine TCR RMSF columns aligned; all18 Q72 distances retained. English labels, Arial, 183 by170mm,450dpi PNG and embedded-font PDF plus editable SVG. Independent full-figure visual and numerical reviews passed. All45 CDR positions and90 model RMSF values unchanged. Both acquired and weakened contacts displayed. Geometry is minimum heavy-atom distance in existing starting structures, with the six-site inclusion rule explicit in methods. No newly run inference, training or MD. All28 manuscript pages independently inspected; figure and full caption on page15. Native10 three-line tables,36 equations,6 images and28 references preserved. Source v52 unchanged.
''',encoding='utf-8')
(ROOT/'修改说明_zh.md').write_text('''# v53 修改说明

图4a改为实际残基接触数据图。六个HLA位点逐行排列，TCR beta48–56逐列排列，每格用两种颜色和形状表示两态的近接触关系。左侧Hchi与行对应，上方RMSF与列对应；Q72行给出全部两态距离，并突出beta53和beta54。b保留六个CDR的全部45个模型预测位点。原文字摘要卡片删除。

方法补充六个HLA位点的几何选择规则；正文及图注统一为最短重原子距离，保留接触增加与减弱的共同信息，并将其与模型两侧动态响应对应。未新增训练、推理或MD。全稿28页渲染及逐页审查通过，三线表与原生公式保留，v52未覆盖。
''',encoding='utf-8')
delivery = {'status':'complete','version':'v53','completed_at_utc':now,'document':{'path':str(doc),'sha256':sha(doc)},'figure':{ext:{'path':str(FIG/f'fig4_mel5_contact_map.{ext}'),'sha256':sha(FIG/f'fig4_mel5_contact_map.{ext}')} for ext in ['png','pdf','svg']},'panels':2,'pages':28,'native_three_line_tables':10,'source_preserved':True,'QA':'passed'}
write(ROOT/'delivery_status.json',delivery)
print(json.dumps({'status':'complete','document':str(doc),'pages':28,'figure':str(png)},ensure_ascii=False))
