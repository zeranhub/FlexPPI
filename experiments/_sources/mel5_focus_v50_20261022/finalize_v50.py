from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
from pypdf import PdfReader

ROOT = Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_focus_v50_20261022')
OLD = ROOT.parent / 'mel5_joint_case_v49_20261019'
DOC = ROOT.parent.parent / 'paper_v50_20261022_中文修订_MEL5残基接触网络与图4精简_三线表.docx'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def write(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding='utf-8')

now = datetime.now(timezone.utc).isoformat()
audit = read(ROOT / 'manuscript_QA/integration_audit.json')
assert sha(DOC) == audit['document_sha256']
assert sha(audit['source']) == audit['source_sha256']
fig_review = read(ROOT / 'figure4/independent_final_review.json')
pages_review = read(ROOT / 'manuscript_QA/independent_pages_1_13_review.json')
assert fig_review['status'] == pages_review['status'] == 'passed'
assert pages_review['document_sha256'] == sha(DOC)

fonts = []
reader = PdfReader(ROOT / 'figure4/fig4_mel5_focus.pdf')
for page in reader.pages:
    for name, indirect in page['/Resources'].get('/Font', {}).items():
        font = indirect.get_object()
        descendants = font.get('/DescendantFonts', [font])
        for descendant in descendants:
            item = descendant.get_object()
            descriptor = item.get('/FontDescriptor')
            descriptor = descriptor.get_object() if descriptor is not None else {}
            embedded = any(key in descriptor for key in ('/FontFile', '/FontFile2', '/FontFile3'))
            fonts.append({'name': str(font.get('/BaseFont', name)), 'embedded': embedded, 'subtype': str(item.get('/Subtype'))})
assert fonts and all(font['embedded'] for font in fonts)
write(ROOT / 'figure4/pdf_font_audit.json', {'status':'passed', 'fonts':fonts, 'pdf_sha256':sha(ROOT / 'figure4/fig4_mel5_focus.pdf')})

root_pages = []
for number in range(14, 29):
    path = ROOT / f'manuscript_QA/pages/page-{number}.png'
    root_pages.append({'page':number, 'path':str(path), 'sha256':sha(path), 'inspection':'complete single-page view_image', 'status':'passed'})
write(ROOT / 'manuscript_QA/root_pages_14_28_review.json', {
    'reviewer':'root', 'reviewed_at_utc':now, 'status':'passed',
    'document_sha256':sha(DOC), 'pages':root_pages,
    'observations':'No clipped text, overlapping elements, broken tables, displaced equations or missing content. Revised Figure 4 and its caption fit page 15; section 3.4 remains intact on page 14.',
    'required_changes':[]
})

figure_audit = read(ROOT / 'figure4/figure_data_audit.json')
figure_audit.update(status='passed', visual_QA='passed', independent_numeric_review='passed', old_panel_d_removed=True, finalized_at_utc=now)
write(ROOT / 'figure4/figure_data_audit.json', figure_audit)
audit.update(status='passed', render_QA='passed', rendered_pages=28, independent_scientific_review='passed', finalized_at_utc=now)
write(ROOT / 'manuscript_QA/integration_audit.json', audit)

figure_map = read(OLD / 'CURRENT_FIGURE_MAP.json')
figure_map['4'] = {
    'source_folder':'figure4', 'absolute_source_folder':str(ROOT / 'figure4'),
    'png':'fig4_mel5_focus.png', 'pdf':'fig4_mel5_focus.pdf', 'svg':'fig4_mel5_focus.svg',
    'new_in_revision':'v50', 'panels':[
        'a Q72-betaPhe53/betaGly54 contact network, published SPR affinity and model response',
        'b absolute model RMSF across all 45 positions in all six CDRs',
        'c all 23 interface Hchi changes ranked, paired with Q72 input-context comparison'
    ], 'previous_panel_d':'removed; existing MD corroboration retained briefly in prose'
}
write(ROOT / 'CURRENT_FIGURE_MAP.json', figure_map)

(ROOT / '修改说明_zh.md').write_text('''# v50 修改记录

案例围绕HLA Q72–TCR βPhe53/βGly54跨界面接触网络展开。正文从MEL5与α24β17的生物学背景、SPR亲和力与晶体结构观察出发，利用全部六个CDR的45个位点和全部23个有效界面位点定位响应，再用Q72侧链状态集中与未替换βGly54主链波动降低解释同一接触片区的稳定化。摘要、方法、图注、讨论和结论同步修改。

图4由四个子图压缩为三个。a将真实结构、文献SPR结果及Q72/βGly54模型变化合并；b保留全部45个CDR位点的两态RMSF曲线；c改为全部23个有效界面位点的Hχ变化排序与Q72的受体输入对照。旧d删除，已有轨迹中的局部接触结果简短保留在正文中。图中使用英文，尺寸由183×185 mm收紧为183×157 mm。

本次未在正文加入四态探索或十个肽位点的失败叙述。新增图使用已有预测与轨迹分析，没有新训练、新推理或新MD。Q72输入对照保留各自pMHC坐标，仅移除输入中的TCR，用于检验模型响应对受体环境的依赖。

全部45×2个RMSF及23×2个Hχ数据、Q72输入对照、残基编号、局部接触距离和正文轨迹占据率均经独立复核；全部28页经渲染与逐页目视检查，10张原生三线表、36处原生公式、6张图片与28条参考文献保留。原v49文件未覆盖。
''', encoding='utf-8')

delivery = {
    'status':'complete', 'version':'v50', 'completed_at_utc':now,
    'document':{'path':str(DOC), 'sha256':sha(DOC)},
    'figure':{ext:{'path':str(ROOT / f'figure4/fig4_mel5_focus.{ext}'), 'sha256':sha(ROOT / f'figure4/fig4_mel5_focus.{ext}')} for ext in ('png','pdf','svg')},
    'source_preserved':True, 'native_three_line_tables':10,
    'QA':{'all_28_pages_visually_reviewed':True, 'independent_source_numeric_review':'passed','figure_visual_review':'passed','figure_PDF_fonts_embedded':True},
    'notes':str(ROOT / '修改说明_zh.md')
}
write(ROOT / 'delivery_status.json', delivery)
print(json.dumps({'status':'complete','document':str(DOC),'pages':28,'figure_panels':3,'PDF_embedded_fonts':len(fonts),'document_sha256':sha(DOC)},ensure_ascii=False))
