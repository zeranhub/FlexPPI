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
png=FIG/'fig4_mel5_compact_joint.png'
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
assert len(PdfReader(QA/'v55_word.pdf').pages)==28
fa=read(FIG/'figure_data_audit.json')
for item in fa['source_files'].values():assert sha(item['path'])==item['sha256']
fonts=[]
for page in PdfReader(FIG/'fig4_mel5_compact_joint.pdf').pages:
    for name,ref in page['/Resources'].get('/Font',{}).items():
        font=ref.get_object()
        for child in font.get('/DescendantFonts',[font]):
            item=child.get_object();fd=item.get('/FontDescriptor')
            fd=fd.get_object() if fd is not None else {}
            fonts.append({'name':str(font.get('/BaseFont',name)),'embedded':any(k in fd for k in ['/FontFile','/FontFile2','/FontFile3'])})
assert fonts and all(f['embedded'] for f in fonts)
write(FIG/'pdf_font_audit.json',{'status':'passed','fonts':fonts,'pdf_sha256':sha(FIG/'fig4_mel5_compact_joint.pdf')})
fa.update(status='passed',visual_QA='passed',source_values_verified=True,finalized_at_utc=now)
write(FIG/'figure_data_audit.json',fa)
audit.update(status='passed',rendered_pages=28,visual_QA='passed',finalized_at_utc=now)
write(QA/'integration_audit.json',audit)
stats={
    'sample_definition':{'observed_bound_complexes':4,'states':['MEL5–ELA 3HG1','alpha24beta17–ELA 4JFF','alpha24beta17–G4A 4JFD','alpha24beta17–I7A 4JFE'],'one_model_output_per_observed_state':True,'not_fixed_backbone_mutagenesis':True,'peptide_values':40,'CDR_positions':45,'absolute_CDR_values':90,'CDR_delta_values':135,'pMHC_interface_union_positions':30,'pMHC_delta_cells':90,'Hchi_comparable_per_contrast':[23,23,22],'ranked_Hchi_plot_removed':True},
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
fm=read(ROOT.parent/'mel5_composite_v54_20261022/CURRENT_FIGURE_MAP.json')
fm['4']={'source_folder':'figure4','absolute_source_folder':str(FIG),'png':'fig4_mel5_compact_joint.png','pdf':'fig4_mel5_compact_joint.pdf','svg':'fig4_mel5_compact_joint.svg','new_in_revision':'v55','panels':['a existing complex overview and local Q72-beta53/54 structure','b all40 four-state peptide RMSF values','c all90 absolute RMSF values at45 CDR positions','d joint response, upper all135 CDR differences and lower all90 pMHC interface cells for three contrasts'],'removed':'v54 d ranked Hchi bar plot','merged':'v54 e and f become the upper/lower blocks of v55 d'}
write(ROOT/'CURRENT_FIGURE_MAP.json',fm)
(FIG/'QA_report.md').write_text("""# Figure4 v55 QA

READY. Four panels a-d. Prior ranked Hchi plot deleted by direct user request; two complete signed maps combined under one d heading with independent metric units and colorbars.183x172mm,600dpi PNG plus embedded-font PDF and editable SVG. Font sizes and all retained data unchanged from v54:40 peptide,90 absolute CDR,135 CDR differences,90 interface cells with valid23/23/22. Source CSVs and array values/masks independently checked against v54. No new inference/training/MD. Entire28-page DOCX independently reviewed, complete figure/caption same page.10 native three-line tables,36 equations,6 images and28 references preserved. Source v54 unchanged.
""",encoding='utf-8')
(ROOT/'修改说明_zh.md').write_text("""# v55 修改说明

按用户要求删除原d子图的23位点Hchi排序柱图，将原e/f两张热图合为新d，CDR RMSF变化位于上部、pMHC Hchi变化位于下部，保留独立单位和色标。a/b/c的内容、物理尺寸、字号及所有保留数据不变，整图由183×214mm收紧至183×172mm。

正文中原e/f引用改为图4d上部/下部；Q72原图d引用指向新d下部，图注删去柱图与p前缀说明。原始数据和v54稿均保留。全稿28页渲染及逐页检查通过，图与完整图注同页，三线表和公式完整。
""",encoding='utf-8')
contract={'scientific_goal':'Remove prior panel d and combine e/f into one d as directly requested.','archetype':'asymmetric_mixed','dimensions_mm':[183,172],'layout':'top structure+peptide, full-width CDR, single joint-response d containing upper CDR and lower pMHC heatmaps','authority':'Direct user instruction; no additional design confirmation needed.','panels':fm['4']['panels'],'source_preserved':True,'removed_panel':'v54 d rank bars','merged_panels':'v54 e/f to v55 d','all_retained_values_unchanged':True,'QA':'passed'}
write(ROOT/'figure_contract.json',contract)
delivery={'status':'complete','version':'v55','completed_at_utc':now,'document':{'path':str(doc),'sha256':sha(doc)},'figure':{ext:{'path':str(FIG/f'fig4_mel5_compact_joint.{ext}'),'sha256':sha(FIG/f'fig4_mel5_compact_joint.{ext}')} for ext in ['png','pdf','svg']},'panels':4,'pages':28,'native_three_line_tables':10,'source_preserved':True,'QA':'passed'}
write(ROOT/'delivery_status.json',delivery)
print(json.dumps({'status':'complete','document':str(doc),'figure':str(png),'panels':4,'pages':28},ensure_ascii=False))
