"""Package only this analysis directory; never include external raw trajectories."""
from pathlib import Path
import shutil, zipfile, json, hashlib, sys, importlib.metadata

R=Path(__file__).resolve().parent
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
assert R.is_relative_to(BASE) and R.name=='analysis_v1_20261026'
audit=json.loads((R/'audit/final_delivery_audit.json').read_text(encoding='utf-8'))
assert audit['status']=='passed'
versions={}
for package in ['numpy','pandas','scipy','pyarrow','matplotlib','Pillow','scikit-learn','joblib','psutil','threadpoolctl','MDAnalysis','mdtraj']:
    try:versions[package]=importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:versions[package]='not installed in this interpreter'
(R/'audit/runtime_versions.json').write_text(json.dumps(dict(python=sys.version,packages=versions),ensure_ascii=False,indent=2),encoding='utf-8')

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

copies=[]
for ext in ['png','pdf','svg']:
    source=R/f'figures/class_core40_absolute_MAE.{ext}'
    target=BASE/f'section33_20261026_功能肽类别分析.{ext}'
    if target.exists():assert sha(target)==sha(source), f'Refusing to replace different file: {target}'
    else:shutil.copy2(source,target)
    copies.append(dict(file=str(target),bytes=target.stat().st_size,sha256=sha(target)))

files=sorted(p for p in R.rglob('*') if p.is_file() and not any(x in {'_staging','__pycache__'} for x in p.relative_to(R).parts))
files=[p for p in files if p.suffix not in {'.tmp','.zip'} and p.name!='archive_contents_sha256.csv']
assert all(p.resolve().is_relative_to(R.resolve()) for p in files)
manifest=R/'audit/archive_contents_sha256.csv'
import csv
with manifest.open('w',encoding='utf-8-sig',newline='') as out:
    writer=csv.writer(out);writer.writerow(['archive_relative_path','bytes','sha256'])
    for p in files:writer.writerow([str(p.relative_to(R)).replace('\\','/'),p.stat().st_size,sha(p)])
files.append(manifest)
archive=BASE/'section33_20261026_功能肽类别分析_源数据与审计.zip'
assert not archive.exists(), f'Archive already exists: {archive}'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
    for p in files:z.write(p,Path('analysis_v1_20261026')/p.relative_to(R))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    assert len(z.infolist())==len(files)
    assert not any('_staging/' in i.filename or '__pycache__/' in i.filename for i in z.infolist())
result=dict(status='complete',analysis_root=str(R),archive=dict(file=str(archive),files=len(files),bytes=archive.stat().st_size,sha256=sha(archive)),figure_copies=copies,raw_trajectories_included=False,manuscript_Word_edited=False)
(BASE/'section33_20261026_功能肽类别分析_交付清单.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False))
