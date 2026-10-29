from pathlib import Path
import hashlib
import json
import zipfile
from datetime import datetime, timezone
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
QA = ROOT / 'manuscript_QA'
FIG = ROOT / 'figure4'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
read = lambda p: json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,obj): Path(p).write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
now=datetime.now(timezone.utc).isoformat()
audit=read(QA/'integration_audit.json')
doc=Path(audit['document'])
assert sha(doc)==audit['document_sha256']
assert sha(audit['source'])==audit['source_sha256']
figure_review=read(FIG/'independent_geometry_visual_review.json')
pages_review=read(QA/'independent_pages_1_13_science_review.json')
assert figure_review['status']==pages_review['status']=='passed'
assert figure_review['reviewed_figure_sha256']==sha(FIG/'fig4_mel5_interface.png')
with zipfile.ZipFile(doc) as z:
    assert hashlib.sha256(z.read('word/media/image41.png')).hexdigest()==sha(FIG/'fig4_mel5_interface.png')

fonts=[]
for page in PdfReader(FIG/'fig4_mel5_interface.pdf').pages:
    for name, indirect in page['/Resources'].get('/Font',{}).items():
        font=indirect.get_object()
        for child in font.get('/DescendantFonts',[font]):
            item=child.get_object();fd=item.get('/FontDescriptor')
            fd=fd.get_object() if fd is not None else {}
            fonts.append({'name':str(font.get('/BaseFont',name)), 'embedded':any(k in fd for k in ['/FontFile','/FontFile2','/FontFile3'])})
assert fonts and all(f['embedded'] for f in fonts)
write(FIG/'pdf_font_audit.json',{'status':'passed','fonts':fonts,'pdf_sha256':sha(FIG/'fig4_mel5_interface.pdf')})
write(QA/'root_pages_14_28_review.json',{'status':'passed','reviewer':'root','document_sha256':sha(doc),'reviewed_at_utc':now,
    'pages':[{'page':n,'sha256':sha(QA/f'pages/page-{n}.png'),'status':'passed','method':'individual complete page view_image'} for n in range(14,29)],
    'required_changes':[],'observations':'Revised case occupies page14, Figure4/caption page15; all references and native tables intact; no clipping, overlaps, missing glyphs or table damage.'})

stats={'sample_definition':{
    'structure_plate':'One existing simulation-prepared starting structure per state (3HG1-derived MEL5;4JFF-derived alpha24beta17). One HLA1–180 CA alignment and one shared camera/crop.',
    'distance':'Two prespecified named atom pairs measured once in each structure. No mean, spread, confidence interval or hypothesis test applies.',
    'predictions':'One deterministic joint-model output per state; all45 aligned CDR positions plotted. Q72 Hchi and betaG54 RMSF values attached to corresponding residues.',
    'interface_selection':'All23 valid same-amino-acid pMHC interface sites retained in source and examined to locate largest Hchi decrease; no rank panel or input-removal contrast plotted.',
    'prose_MD':'Three existing independent ITrace trajectories per state; equal-weight means for local contact and same-frame two-hydrogen-bond occupancy; no MD profile shown.'},
    'model':{'checkpoint_sha256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5','training_records':2430,'training_sources':{'ATLAS':1866,'DynaRepo':564},'ITrace_in_training':False,'split':'all-chain family fold0','outputs':'RMSF in Angstrom; fixed-chi joint-state entropy Hchi in nats','new_training':0,'new_inference':0,'new_MD':0},
    'source_files':{str(p.relative_to(ROOT)):{'sha256':sha(p)} for p in [FIG/'source_CDR_RMSF_all45.csv',FIG/'source_interface_Hchi_all23.csv',FIG/'source_atom_pair_distances.csv',ROOT/'analysis/geometry_audit.json',ROOT/'structure/dual_state_molecular_view.json']},
    'statistical_tests':None,'multiple_comparison_correction':None,'error_bars':None,
    'scope':'Paired local geometry and model dynamic response. Distances are not direct affinity predictions and no direction of Hchi-to-RMSF causation is inferred.'}
write(FIG/'statistics_and_reproducibility.json',stats)

fa=read(FIG/'figure_data_audit.json');fa.update(status='passed',visual_QA='passed',independent_geometry_QA='passed',finalized_at_utc=now);write(FIG/'figure_data_audit.json',fa)
audit.update(status='passed',render_QA='passed',rendered_pages=28,independent_science_review='passed',finalized_at_utc=now);write(QA/'integration_audit.json',audit)
fm=read(ROOT.parent/'mel5_focus_v50_20261022/CURRENT_FIGURE_MAP.json')
fm['4']={'source_folder':'figure4','absolute_source_folder':str(FIG),'png':'fig4_mel5_interface.png','pdf':'fig4_mel5_interface.pdf','svg':'fig4_mel5_interface.svg','new_in_revision':'v51','panels':['a paired same-view interface structures with named-atom distances and residue-specific model readouts','b all45CDR paired absolute model RMSF'],'removed_panels':['v50 c ranked interface Hchi and TCR-removal input comparison','v49 d independent MD occupancy']};write(ROOT/'CURRENT_FIGURE_MAP.json',fm)

(ROOT/'Genome_Biology参考与构图说明.md').write_text('''# 图4的参考与构图

实际目视参考：Velando等，Genome Biology (2025) 26:260，doi:10.1186/s13059-025-03703-6，Figure 4A和Figure 5。正式出版PDF由作者机构CSHL公开存档读取，并逐页渲染、目视检查第9/10页。Figure4A以轻背景骨架衬托关键残基与接触；Figure5采用相同坐标结构和直接曲线命名。新版借鉴局部主体与标签层次，不复用原文图像。

官方链接：https://link.springer.com/article/10.1186/s13059-025-03703-6/figures/4 与 https://link.springer.com/article/10.1186/s13059-025-03703-6/figures/5 。作者机构正式版PDF：https://repository.cshl.edu/id/eprint/41951/1/10.1186.s13059-025-03703-6.pdf 。

另查阅Song等Genome Biology (2025) 26:283的官方文章与Figure7图注（doi:10.1186/s13059-025-03744-x），其中采用前后局部结构显示特定位点环境变化；以及2019年剪接体论文Figure3b的跨蛋白局部图注（doi:10.1186/s13059-019-1680-9）。这两篇实际图像未成功获取，不将图注阅读记作目视审查。

新版a用双态相同视角和尺度，直接呈现Q72侧链OE1–β53主链N及NE2–βGly54主链O两组几何关系，把两项模型读数标在同一结构中的对应残基上。b保留全部45个CDR位点，提供局部Gly54响应在完整识别环中的位置。旧c全部删除，不另外补一个弱面板。图中只有结构配对线，没有Hχ→RMSF的单向因果箭头。
''',encoding='utf-8')
(ROOT/'修改说明_zh.md').write_text('''# v51 修改记录

图4改为a、b两个子图。a以HLA刚性对齐后的MEL5与α24β17同视角、同尺度局部结构为主体，直接标出Q72、β52、βIle53/βPhe53、βGly54与两对原子距离，把Q72 Hχ和Gly54 RMSF贴到对应残基。小复合物全景仅作定位；b保留六个CDR的全部45个对齐位点曲线。旧c的23位点排序与TCR移除输入图已删除，正文对应叙述删除；完整23位点源表仍保留供追溯。

方法与3.4及图注已同步重写。案例主线是：受体改造伴随βCDR2跨界面接触重排，Q72与β53/β54骨架的两端配对靠近，模型把Q72侧链状态集中和未替换Gly54主链活动收紧对应到同一局部网络。原研究SPR与晶体结构、已有ITrace局部轨迹作为生物学背景和动态佐证；不将两项输出之间的因果方向写成已经证明的结论。

坐标距离逐对复算：Q72 OE1–β53N为5.85→2.94 Å，Q72 NE2–βGly54O为10.58→3.12 Å。原子与残基标签经过独立放大复核。论文全部28页逐页渲染审查通过，10张原生三线表、36处原生公式、6张图片及28条文献保留；v50未覆盖。没有新增训练、模型推理或MD。
''',encoding='utf-8')
(FIG/'QA_report.md').write_text('''# Figure4 v51 QA

READY. AP0–AP7: baseline typography/color/export retained; no default palette, rainbow, four-sided data frames, legend occlusion or raster-only export. Molecular images are real-coordinate PyMOL renders. CL1–CL7:183×140 mm, minimum6pt,450dpi preview, PDF embedded TrueType fonts, SVG editable text,0.6pt axes, outward ticks. VI: two nonredundant panels; same coordinates, camera and crop for both states; all45CDR retained; source23interface positions retained; no synthetic atoms, coordinate editing or causal output arrows. No inferential statistical test or errorbar applies to single deterministic paired model outputs. VV: full preview and local enlarged views independently inspected; critical endpoint labels moved away from N/NE2 atoms; no remaining required fixes. Figures are English. All28 manuscript pages individually inspected, including final Figure4 and caption on page15.

See statistics_and_reproducibility.json for data units, sample definition, split/model and traceability. Independent review is bound to the final PNG SHA256 in independent_geometry_visual_review.json.
''',encoding='utf-8')
delivery={'status':'complete','version':'v51','completed_at_utc':now,'document':{'path':str(doc),'sha256':sha(doc)},
          'figure':{ext:{'path':str(FIG/f'fig4_mel5_interface.{ext}'),'sha256':sha(FIG/f'fig4_mel5_interface.{ext}')} for ext in ['png','pdf','svg']},
          'panels':2,'native_three_line_tables':10,'pages':28,'source_preserved':True,'QA':'passed'}
write(ROOT/'delivery_status.json',delivery)
print(json.dumps({'status':'complete','document':str(doc),'figure_panels':2,'pages':28,'fonts_embedded':len(fonts)},ensure_ascii=False))
