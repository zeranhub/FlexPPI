from pathlib import Path
import hashlib, json, os, zipfile
from datetime import datetime, timezone
from PIL import Image
from pypdf import PdfReader

root = Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_absolute_affinity_20261028')
paper_dir = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
target = paper_dir / 'paper_v126_20261028_中文修订_绝对亲和力对照与肽突变验证.docx'
source = paper_dir / 'paper_v125_20261028_中文修订_MEL5核心优势叙述.docx'
render = root / 'manuscript/render/v126_final'
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def save(path,data):
    Path(path).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
now = datetime.now(timezone.utc).isoformat()
change_path = root/'manuscript/revision_changes.json'
changes = json.loads(change_path.read_text(encoding='utf-8-sig'))
assert sha(source)==changes['source_sha256']
assert sha(target)==changes['target_sha256']
with zipfile.ZipFile(source) as old, zipfile.ZipFile(target) as new:
    media=[n for n in old.namelist() if n.startswith('word/media/')]
    assert len(media)==46 and all(old.read(n)==new.read(n) for n in media)
    changed=[n for n in old.namelist() if old.read(n)!=new.read(n)]
    assert changed==changes['changed_original_parts']
pdf=render/(target.stem+'.pdf')
assert len(PdfReader(str(pdf)).pages)==33
page_records=[]
for n in range(1,34):
    p=render/f'page-{n}.png'
    with Image.open(p) as im: dims=list(im.size)
    page_records.append({'page':n,'sha256':sha(p),'pixels':dims,'visual_review':'pass'})
qa={
 'reviewed_utc':now,'status':'passed','document':str(target),'document_sha256':sha(target),
 'render_method':'Packaged render_docx.py with read-only Microsoft Word COM PDF export; PDF rasterized for page review.',
 'page_count':33,'all_33_pages_visually_reviewed':True,
 'checks':{'Chinese_and_Greek_glyphs':'pass','table_S7_ID_and_numeric_columns':'pass','figure_S5_S7_legibility_and_captions':'pass','clipping_and_text_overlap':'none observed','orphaned_overflow_pages':'none in final','original_46_media_binary_preserved':True,'original_document_parts_changed':changed},
 'repairs_before_final_render':['Table S7 fixed layout and zero first-line indentation to prevent PDB label wrapping.','Appendix F wording shortened with full implementation detail retained in result report; eliminated two overflow pages.','Figure S6 title spacing increased to prevent overlap with descriptive statistics.'],
 'independent_content_check':'Dynamics agent compared final main text and result report with source tables; numeric values and limitations passed.',
 'pages':page_records,'review_pdf_sha256':sha(pdf),
 'distribution':'Review PDF and page images are QA artifacts; final manuscript is delivered as DOCX.'
}
save(root/'manuscript/render_QA.json',qa)
changes['render_validation']='passed: 33 pages rendered and visually reviewed; see render_QA.json'
save(change_path,changes)
figures=[]
for stem in ['figS5_absolute_affinity','figS6_peptide_dynamic_validation','figS7_peptide_contact_redistribution']:
    record={'figure':stem,'visual_review':'pass','files':{}}
    for ext in ['pdf','svg','png']:
        p=root/'figures'/f'{stem}.{ext}'
        assert p.is_file() and p.stat().st_size>0
        record['files'][ext]={'bytes':p.stat().st_size,'sha256':sha(p)}
        if ext=='pdf': assert len(PdfReader(str(p)).pages)==1
        if ext=='png':
            with Image.open(p) as im:
                record['files'][ext]['pixels']=list(im.size)
                record['files'][ext]['dpi']=list(im.info.get('dpi',()))
                assert all(v>299 for v in im.info.get('dpi',()))
    figures.append(record)
figqa={
 'reviewed_utc':now,'status':'passed','width_mm':183,'typography':'Arial; base 8 pt, axes 7 pt, ticks 6.5 pt',
 'exports':'Single-page vector PDF, SVG and 300 dpi PNG for each figure.',
 'visual_checks':{'title_statistics_overlap':'fixed in S6 and final image inspected','axes_units_labels_legend':'pass','manuscript_render_legibility':'pass','cropped_text':'none observed'},
 'data_and_statistics_checks':[
 'S5 retains all four states and both models, independently predicted before comparison; no four-point population performance claim.',
 'S6 retains all 10 peptide and 45 CDR positions per contrast; chemically valid same-type pMHC Hchi positions are 23 and 22, all 30 retained in CSV.',
 'Independent MD trajectory is the repeat unit; frames, residues and cross-state 3x3 differences are not biological repeats.',
 'S6 independent-difference SE recomputation max discrepancy 8.3e-17; all three CDR3beta trajectory means per state shown without artificial pairing.',
 'S7 retains all 880 peptide-neighbor pairs per contrast including zeros; selected local panels disclosed as exploratory.',
 'S7 five selected pairs across nine trajectories independently recomputed, 45 checks passed; both 4 and 4.5 Angstrom thresholds reported.',
 'G4A negative agreement and model underestimation retained; no p-value or free-energy causality claim.',
 'PRODIGY version and ProAffinity-GNN simplified inference correction/input-preparation differences disclosed.'
 ],'figures':figures
}
save(root/'figures/figure_QA.json',figqa)
(root/'figures/figure_QA.md').write_text('''# 图S5–S7最终核验\n\n三图均通过源表、数值、统计定义和实际图像核验，已同时检查论文渲染中的嵌入效果。输出为183 mm双栏宽度的PDF、SVG及300 dpi PNG。S6标题与统计文字重叠已修正，最终无文字裁切、图例遮挡或面板重叠。逐文件尺寸和SHA见figure_QA.json。\n\nS5完整显示四态实验与两种模型，未将相关四态当作独立家族评估。S6保留完整十肽、45个CDR位置和全部有效同型pMHC侧链23/22位；30个预定界面位置均在源表。均值差SE由两态独立三轨迹计算，独立复核最大差8.3×10⁻¹⁷。轨迹、残基、保存帧及3×3差值分别说明，未制造生物重复。S7保留每种对照全部880对含零接触组合，局部案例标明完整筛查后的探索性核查，独立45次距离重算通过，两个距离阈值均展示。\n\nI7A的区域支持、G4A的弱恢复及模型对部分幅度的低估均保留。结合态动态与接触占据不作总结合熵、结合自由能或亲和力因果贡献解释。所有图中数值可从同目录statistics.json及对应CSV追溯。最终正文、说明文档也由动态分析代理作独立只读交叉核查，未发现实质数值或论证问题。\n''',encoding='utf-8')

# Only small source/result artifacts; large runtimes, weights and trajectory caches remain local.
excluded_dirs={'venv','venv_prodigy','package_downloads','esm_feature_cache','esm2_t33_650M_UR50D','__pycache__','render','recomputed_labels'}
excluded_names={'prodigy_form.html','prodigy_frontend.js','ProAffinity-GNN_inference_compatible.py'}
allowed={'.py','.ps1','.md','.txt','.csv','.tsv','.json','.log','.pdb','.pdbqt','.fasta','.fa','.png','.svg','.pdf'}
files=[]
for folder,dirs,names in os.walk(root,followlinks=False):
    dirs[:]=[d for d in dirs if d not in excluded_dirs]
    for name in names:
        p=Path(folder)/name
        if name in excluded_names or name=='result_package_manifest.json' or p.suffix.lower() not in allowed:continue
        files.append(p)
files.sort(key=lambda p:p.relative_to(root).as_posix())
manifest={
 'created_utc':now,'source_root':str(root),'files':[{'path':p.relative_to(root).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files],
 'excluded':'Python environments, package downloads, ESM/GNN weights, feature/recomputed-label/distance caches, raw trajectories, manuscript QA page renders and PDF; no original manuscript or final DOCX included.',
 'manuscript_delivered_separately':str(target),
 'reproduction_scope':'Figures can be regenerated from included CSVs; model reruns require frozen dependencies and weights; trajectory reanalysis requires original ITrace and project label/prediction sources listed in audits.'
}
save(root/'result_package_manifest.json',manifest)
archive=paper_dir/'MEL5_绝对亲和力与肽突变验证_20261028.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as out:
    for p in files+[root/'result_package_manifest.json']:
        out.write(p,arcname=p.relative_to(root).as_posix())
with zipfile.ZipFile(archive) as check:
    assert check.testzip() is None
    assert len(check.namelist())==len(files)+1
    for item in manifest['files']:
        assert hashlib.sha256(check.read(item['path'])).hexdigest()==item['sha256']
print(json.dumps({'docx':str(target),'docx_sha256':sha(target),'pages':33,'original_media_preserved':46,'package':str(archive),'package_bytes':archive.stat().st_size,'package_files':len(files)+1,'package_integrity':'passed','qa':'passed'},ensure_ascii=False,indent=2))
