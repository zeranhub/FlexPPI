from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,zipfile,csv
from PIL import Image
from pypdf import PdfReader
R=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028')
P=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
D=P/'paper_v127_20261028_中文修订_TCR结合模型对照.docx'
render=R/'manuscript/render/v127'
oldrender=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_absolute_affinity_20261028\manuscript\render\v126_final')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
changes=json.loads((R/'manuscript/revision_changes.json').read_text(encoding='utf-8'))
assert sha(D)==changes['target_sha256']
texts=json.loads((render/'page_text.json').read_text(encoding='utf-8'))
newpages=[x['page'] for x in texts if not x['pixel_identical_to_previous_reviewed_page']]
assert newpages==[13,15,17,34,35]
qa=dict(status='passed',reviewed_utc=datetime.now(timezone.utc).isoformat(),document=str(D),document_sha256=sha(D),page_count=35,all_pages_reviewed=True,review_method='Read-only Word COM PDF export through packaged render_docx PNG rasterizer. All five changed/new pages inspected. Other 30 pages are pixel-identical to individually inspected v126 pages.',new_or_changed_pages_visually_inspected=newpages,identical_to_previous_reviewed_pages=[x['page'] for x in texts if x['pixel_identical_to_previous_reviewed_page']],source_previous_QA=r'E:\Project\FlexPPI\experiments\_sources\mel5_absolute_affinity_20261028\manuscript\render_QA.json',layout_checks={'table_S8_headers_numeric_alignment':'pass','figure_S8_and_caption_legibility':'pass','Chinese_Greek_and_math_glyphs':'pass','overlap_clipping_or_overflow':'none observed','original_media_preserved_from_v126':49,'original_tables_fields_comments_sections_preserved':True},scientific_crosscheck='Independent ERGO agent rechecked all 72 merged scores against each model source and main/supplement wording; passed.',pages=[dict(page=x['page'],sha256=sha(render/f"page-{x['page']}.png"),review='new visual inspection' if x['page'] in newpages else 'pixel-identical to prior visual inspection') for x in texts])
save(R/'manuscript/render_QA.json',qa)
changes['render_validation']='passed: 35 pages, five changed/new pages visually inspected, 30 others pixel-identical to prior individually inspected pages'
save(R/'manuscript/revision_changes.json',changes)
figfiles={}
for ext in ['pdf','svg','png']:
 p=R/f'figures/figS8_TCR_binding_scores.{ext}'
 figfiles[ext]={'sha256':sha(p),'bytes':p.stat().st_size}
 if ext=='pdf':assert len(PdfReader(p).pages)==1
 if ext=='png':
  with Image.open(p) as im:
   figfiles[ext]['pixels']=list(im.size);figfiles[ext]['dpi']=list(im.info['dpi']);assert min(im.info['dpi'])>299
save(R/'figures/figure_QA.json',dict(status='passed',files=figfiles,width_mm=183,height_mm=132,visual_review='All 12 labels, six configurations, 72 score labels, SPR binding column, colorbar and footnote inspected in PNG and manuscript page35; no clipping or overlap.',data_checks='All 72 values match merged source and each model raw outputs; all correct and incorrect mutation directions retained; no score cutoff in figure.',statistical_limits='Single TCR family, known training overlaps disclosed; CV models not biological repeats; not an overall performance benchmark; cross-method scales are not calibrated probabilities.',peer_review='Independent model agent verified merged scores and main/supplement claims.'))
(R/'figures/figure_QA.md').write_text('''# 图S8最终核验\n\n12个文献条件、三种方法六种配置和72个评分全部保留。图中SPR阳性与NB对应原文条件，未将弱结合重标阴性，也未以0.5切分图像。已核对原始输出、合并表、每格数字及行列顺序；另一模型代理对全部72值作了独立交叉核查，结果一致。\n\n实际查看300 dpi PNG及文稿第35页后确认标题、行标签、数值、色条和注释均清楚，无重叠与裁切。PDF单页，SVG文字可编辑，183×132 mm。固定0–1色标只表示各配置原始分数，不表示跨方法统一校准。已披露训练覆盖、NetTCR四模型的选择依据和全部20模型敏感性结果。没有按结果选择展示配置，没有将模型成员当作生物重复。文件SHA与尺寸见figure_QA.json。\n''',encoding='utf-8')
# Preserve small official sources in their expected relative paths, excluding weights/raw training.
prune={'venv','__pycache__','archive_source','deps','training','download','render','.git','models','weights'}
allowed={'.py','.ps1','.md','.json','.csv','.txt','.log','.png','.svg','.pdf','.yml'}
excluded={'upstream.zip','upstream_master.zip','Madura2013_original.pdf','tcr3d_class1.html'}
files=[]
for folder,dirs,names in os.walk(R,followlinks=False):
 dirs[:]=[d for d in dirs if d not in prune and not (d=='data' and 'upstream' in Path(folder).parts)]
 for name in names:
  p=Path(folder)/name
  if p.stat().st_size==0 or name in excluded or name=='result_package_manifest.json':continue
  if p.suffix.lower() in allowed or name in ['LICENSE','LICENCE']:files.append(p)
ergo=R/'ergo/download/ERGO-II-master'
for name in ['Predict.py','Trainer.py','Models.py','Loader.py','Sampler.py','README.md','LICENSE','requirements.txt']:
 p=ergo/name
 if p.is_file():files.append(p)
files=sorted(set(files),key=lambda p:p.relative_to(R).as_posix())
manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),files=[{'path':p.relative_to(R).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files],excluded='Runtime/dependency environments, model weights, raw training libraries, upstream ZIP downloads, rendering page images/PDF, raw MD. Final DOCX delivered separately.',reproduction='Plotting and merged-result calculations can run from included CSV. Model re-inference requires official weights/training vocabularies at frozen versions and recorded SHA; download scripts and public source URLs included.',manuscript=str(D),manuscript_sha256=sha(D))
save(R/'result_package_manifest.json',manifest)
Z=P/'MEL5_TCR结合模型对照_20261028.zip'
with zipfile.ZipFile(Z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as out:
 for p in files+[R/'result_package_manifest.json']:out.write(p,arcname=p.relative_to(R).as_posix())
with zipfile.ZipFile(Z) as z:
 assert z.testzip() is None
 for item in manifest['files']:assert hashlib.sha256(z.read(item['path'])).hexdigest()==item['sha256']
print(json.dumps({'document':str(D),'pages':35,'QA':'passed','archive':str(Z),'archive_bytes':Z.stat().st_size,'archive_files':len(files)+1,'integrity':'passed'},ensure_ascii=False,indent=2))
