from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, zipfile

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
read = lambda path: json.loads(Path(path).read_text(encoding='utf-8'))

protocol_sha = 'b41817fb7b7536db8644124fc5b5b29450c33bac1fb992f509852a15e5f748f4'
assert sha(ROOT/'protocol.json') == protocol_sha
assert sha(ROOT/'direct_affinity/source_only_protocol.json') == '3f5d9c3cf082900e76cf2462a85ea008c6a2fca316d84ca395ad9bc7775c212d'
functional = read(ROOT/'independent_functional_audit.json')
assert functional['status'] == 'passed'
for name, digest in functional['file_hashes'].items():
    assert sha(ROOT/'functional_screening'/name) == digest
direct = read(ROOT/'direct_affinity/independent_numeric_verification.json')
assert direct['status'] == 'passed'
figure = read(ROOT/'result_figures/figure_manifest.json')
assert figure['visual_qa']['status'] == 'passed'
for item in figure['source']:
    assert sha(item['path']) == item['sha256']
for item in figure['figures']:
    for f in item['files']:
        assert sha(f['path']) == f['sha256']
model = BASE/'local_antigen_screening_20261021/model/epoch_040.pt'
assert sha(model) == 'b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f'
paper = BASE.parent/'paper_v49_20261019_中文修订_MEL5双动态机制实例_三线表.docx'
assert sha(paper) == '5f79ddf7d657dc9b782c9fb432e373a1e7c170c1c88d019e528fbf73c366c222'
files = [p for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name != 'delivery_audit.json' and '__pycache__' not in p.parts]
audit = dict(status='passed', finished_utc=datetime.now(timezone.utc).isoformat(), rule='signed largest peptide RMSF decrease over all 9 sites', functional_independent_audit='passed', direct_SPR_independent_audit='passed', figure_visual_QA='passed', root_final_images_inspected=['max_drop_overview.png','max_drop_all_contexts.png'], root_final_visual_QA='passed: readable labels, complete contexts, correct units, no clipping or collisions', main_manuscript_path=str(paper), main_manuscript_sha256=sha(paper), main_manuscript_modified=False, no_new_training=True, no_new_neural_inference=True, no_new_MD=True, evidence_scope='exploratory retrospective historical model; functional hit retention and separate small selected direct SPR panels', files=[dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size, sha256=sha(p)) for p in files])
(ROOT/'delivery_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2),encoding='utf-8')
archive = BASE/'任一肽位点RMSF降幅_完整复核_20261022.zip'
files.append(ROOT/'delivery_audit.json')
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as package:
    for p in files:
        package.write(p, str(Path(ROOT.name)/p.relative_to(ROOT)))
with zipfile.ZipFile(archive) as package:
    assert package.testzip() is None
    for p in files:
        key = str(Path(ROOT.name)/p.relative_to(ROOT)).replace('\\','/')
        assert hashlib.sha256(package.read(key)).hexdigest() == sha(p)
print(json.dumps(dict(status='passed', archive=str(archive), archive_bytes=archive.stat().st_size, archive_sha256=sha(archive), file_count=len(files), manuscript_unchanged=True), ensure_ascii=False,indent=2))
