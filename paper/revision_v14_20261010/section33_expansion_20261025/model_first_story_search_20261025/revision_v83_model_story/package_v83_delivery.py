from pathlib import Path
import json,hashlib,zipfile
base=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010');v=base/'section33_expansion_20261025/model_first_story_search_20261025/revision_v83_model_story'
report=v/'render/publish_comparison.json';d=json.loads(report.read_text(encoding='utf-8'));d['visual_review']='passed: pages10/24 inspected after final wording; other27 pages exactly match reviewed final render';d['render_directory']=str(v/'render/release_publish');report.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
a=v/'manuscript/v83_build_audit.json';p=json.loads(a.read_text(encoding='utf-8'));p['render_visual_QA_pending']=False;p['render_visual_QA']='29 pages; main and supplemental figures, captions, all changed pages reviewed; remaining bodies verified against reviewed v82';p['final_render']=str(v/'render/release_publish');a.write_text(json.dumps(p,ensure_ascii=False,indent=2),encoding='utf-8')
archive=base/'paper_v83_20261025_图3与补充图S3_源文件与PyMOL.zip'
files=[]
for name in ['figure3','supplement_cam','pymol']:
    files.extend(x for x in (v/name).rglob('*') if x.is_file() and '__pycache__' not in x.parts)
files.append(v/'palette_audit/palette_decision.json')
for name in ['build_v83_manuscript.py','v83_text_plan.json','v79_paragraph_inventory.json','v83_build_audit.json']:files.append(v/'manuscript'/name)
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for x in sorted(set(files)):z.write(x,x.relative_to(v))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
summary={'paper':p['output'],'paper_sha256':p['output_sha256'],'archive':str(archive),'archive_bytes':archive.stat().st_size,'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'source_files':len(set(files)),'pages':29,'review':'passed','main_plot_values':568,'supplement_plot_values':900,'final_palette':'NHE1 #B87936; IQCG #2166AC; key Met #3A8D8C'}
(v/'v83_delivery_audit.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False))