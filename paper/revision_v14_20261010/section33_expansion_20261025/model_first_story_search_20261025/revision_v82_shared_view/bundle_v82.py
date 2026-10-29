from pathlib import Path
import json,shutil,zipfile,hashlib
base=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010');root=base/'section33_expansion_20261025/model_first_story_search_20261025';v=root/'revision_v82_shared_view';f=v/'figure3';p=v/'pymol'
shutil.copy2(root/'revision_v80_redesign/pymol/CaM_continuous_context_v80.pse',p/'CaM_local_context_from_v80.pse')
qa=json.loads((f/'figure_QA.json').read_text(encoding='utf-8'));qa['pdf_fonts']='Arial and Arial Bold: embedded CID TrueType, verified with pdffonts';(f/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
archive=base/'paper_v82_20261025_图3统一视角_源文件与PyMOL.zip'
files=[x for x in f.rglob('*') if x.is_file()]
files.extend(x for x in (p/'inputs').glob('*.pdb'))
for name in ['CaM_shared_reference_v82.pse','CaM_local_context_from_v80.pse','render_shared_view_v82.py','render_manifest.json','README_shared_reference_zh.md','whole_NHE1.png','whole_IQCG.png','whole_NHE1_trim.png','whole_IQCG_trim.png']:
    files.append(p/name)
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for x in files:assert x.is_file();z.write(x,x.relative_to(v))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print(json.dumps({'archive':str(archive),'files':len(files),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest()},ensure_ascii=False))