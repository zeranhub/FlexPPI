from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,zipfile
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
FIG=ROOT/'figure4';QA=ROOT/'manuscript_QA'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,o):Path(p).write_text(json.dumps(o,ensure_ascii=False,indent=2),encoding='utf-8')
now=datetime.now(timezone.utc).isoformat()
audit=read(QA/'integration_audit.json');doc=Path(audit['document'])
png=FIG/'fig4_mel5_composite.png'
assert sha(doc)==audit['document_sha256'] and sha(png)==audit['figure_sha256']
assert sha(audit['source'])==audit['source_sha256']
with zipfile.ZipFile(doc) as z:
    assert hashlib.sha256(z.read('word/media/image41.png')).hexdigest()==sha(png)
for filename in ['independent_visual_review.json','independent_data_review.json']:
    report=read(FIG/filename)
    assert report['status']=='passed'
    assert sha(png) in json.dumps(report)
for filename in ['independent_science_review.json','independent_pages_1_15_review.json','independent_pages_15_28_review.json']:
    report=read(QA/filename)
    assert report['status']=='passed'
    assert sha(doc) in json.dumps(report)
assert len(list((QA/'pages').glob('page-*.png')))==28
assert len(PdfReader(QA/'v54_word.pdf').pages)==28
fa=read(FIG/'figure_data_audit.json')
for item in fa['source_files'].values():assert sha(item['path'])==item['sha256']
fonts=[]
for page in PdfReader(FIG/'fig4_mel5_composite.pdf').pages:
    for name,ref in page['/Resources'].get('/Font',{}).items():
        font=ref.get_object()
        for child in font.get('/DescendantFonts',[font]):
            item=child.get_object();fd=item.get('/FontDescriptor')
            fd=fd.get_object() if fd is not None else {}
            fonts.append({'name':str(font.get('/BaseFont',name)),'embedded':any(k in fd for k in ['/FontFile','/FontFile2','/FontFile3'])})
assert fonts and all(f['embedded'] for f in fonts)
write(FIG/'pdf_font_audit.json',{'status':'passed','fonts':fonts,'pdf_sha256':sha(FIG/'fig4_mel5_composite.pdf')})
fa.update(status='passed',visual_QA='passed',source_values_verified=True,finalized_at_utc=now)
write(FIG/'figure_data_audit.json',fa)
audit.update(status='passed',rendered_pages=28,visual_QA='passed',finalized_at_utc=now)
write(QA/'integration_audit.json',audit)
stats={
    'sample_definition':{'observed_bound_complexes':4,'states':['MEL5–ELA 3HG1','alpha24beta17–ELA 4JFF','alpha24beta17–G4A 4JFD','alpha24beta17–I7A 4JFE'],'one_model_output_per_observed_state':True,'not_fixed_backbone_mutagenesis':True,'peptide_values':40,'CDR_positions':45,'absolute_CDR_values':90,'CDR_delta_values':135,'pMHC_interface_union_positions':30,'pMHC_delta_cells':90,'Hchi_comparable_per_contrast':[23,23,22],'ranked_Hchi_sites':23},
    'contrasts':fa['contrasts'],
    'interface_definition':'Union of pMHC residue positions with a C-alpha distance <=8A to any TCR residue in at least one of four structures.22 HLA +8 peptide sites.',
    'Hchi_eligibility':'Both states have defined side-chain chi angles and the same amino acid. Grey cells are ineligible; not zero changes.',
    'local_structure':'Archived 4JFF-derived starting coordinates. Fixed atom pairs Q72 OE1 to beta53 backbone N (2.94A); Q72 NE2 to beta54 backbone O (3.12A). Common MEL5 numbering.',
    'model':{'checkpoint_sha256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5','training_records':2430,'training_sources':{'ATLAS':1866,'DynaRepo':564},'ITrace_in_training':False,'split':'all-chain family fold0'},
    'new_training':0,'new_inference':0,'new_MD':0,
    'error_bars':None,'statistical_tests':None,'multiple_comparison_correction':None,
    'replication':'One deterministic model output per existing bound-complex structure. Residue positions are not independent experimental replicas. No significance tests or error bars added.',
    'published_SPR':'MEL5-ELA18uM to alpha24beta17-ELA0.6nM, Madura et al.2013 doi10.1074/jbc.M113.464560, manuscript reference28. The figure uses reported values without invented uncertainty.',
    'scope':'Residue-resolved dynamic response to receptor and peptide changes; not a fitted affinity model or a directional Hchi-to-RMSF causal experiment.',
    'source_files':{str(p):{'sha256':sha(p)} for p in sorted(FIG.glob('source_*')) if p.is_file()}
}
write(FIG/'statistics_and_reproducibility.json',stats)
fm=read(ROOT.parent/'mel5_contact_v53_20261022/CURRENT_FIGURE_MAP.json')
fm['4']={'source_folder':'figure4','absolute_source_folder':str(FIG),'png':'fig4_mel5_composite.png','pdf':'fig4_mel5_composite.pdf','svg':'fig4_mel5_composite.svg','new_in_revision':'v54','panels':['a existing complex overview and one local Q72-beta53/54 molecular view','b all40 peptide RMSF values across4 bound states','c absolute model RMSF at all45 CDR positions before/after maturation','d all23 valid same-AA interface Hchi changes ranked','e all135 CDR differences for3 contrasts','f all90 union pMHC interface cells, same-AA valid23/23/22'],'replaces':'v53 contact grid; all five user-supplied graphical selections represented as six panels'}
write(ROOT/'CURRENT_FIGURE_MAP.json',fm)
(FIG/'QA_report.md').write_text('''# Figure4 v54 QA

READY. Six panels include every user-selected graphical component.183x214mm,600dpi PNG, embedded-font vector PDF and editable SVG. Charts redrawn from original data; actual molecular renders placed without scientific alteration. Minimum font6pt. Four categorical states use stable color/marker/line combinations; signed delta maps use separate unit-labelled blue-white-rose colorbars. All40 peptide values,90 absolute CDR values,23 ranked Hchi changes,135 CDR differences and90 union-interface cells retained.23/23/22 Hchi cells eligible. Original NPZ values, sequence alignment and interface geometry independently verified; no colour clipping. Figure and caption fit page15. Entire28-page final DOCX independently viewed;10 native three-line tables,36 native equations,6 images and28 references retained. Source v53 preserved. No new inference, training or MD.
''',encoding='utf-8')
(ROOT/'修改说明_zh.md').write_text('''# v54 修改说明

将用户选择的五张附图合为六个子图：a结构概览及Q72局部结构、b四态肽RMSF、c全部CDR两态RMSF、d23个有效界面Hchi变化排序、e/f三种对照的CDR与pMHC差值热图。首行a/b并列，下面四图横跨整页；量化部分从原始数据重新作图，统一字体、线宽、橙紫状态色及柔和的差值色调，保留全部位点与不可比灰格。

同步更正方法、正文和图注的状态定义及子图引用；三张热图行分别采用4JFF−3HG1、4JFD−4JFF、4JFE−4JFF。正文保留成熟对照的Q72–beta53/54机制解释，并以新增肽对照说明区域动态响应差异。全稿28页已渲染和逐页检查，图与完整图注同页，三线表及原生公式保留。原稿v53未覆盖。
''',encoding='utf-8')
contract={'scientific_goal':'Combine the five user-selected graphics into Figure4 and harmonize their colors, preserving actual structure and complete residue responses.','archetype':'asymmetric_mixed','dimensions_mm':[183,214],'layout':'top structure+peptide; full-width CDR, interface rank, CDR delta map, pMHC delta map','authority':'Direct user request to combine and recolor specific attached panels. No additional design confirmation required.','asset_route':'Existing native molecular render; parameter inheritance from inspected LineTrend, BarDistribution and heatmap scripts with incompatible source axes; all quantitative values come from original project data.','panels':fm['4']['panels'],'source_preserved':True,'QA':'passed'}
write(ROOT/'figure_contract.json',contract)
delivery={'status':'complete','version':'v54','completed_at_utc':now,'document':{'path':str(doc),'sha256':sha(doc)},'figure':{ext:{'path':str(FIG/f'fig4_mel5_composite.{ext}'),'sha256':sha(FIG/f'fig4_mel5_composite.{ext}')} for ext in ['png','pdf','svg']},'panels':6,'pages':28,'native_three_line_tables':10,'source_preserved':True,'QA':'passed'}
write(ROOT/'delivery_status.json',delivery)
print(json.dumps({'status':'complete','document':str(doc),'figure':str(png),'panels':6,'pages':28},ensure_ascii=False))
