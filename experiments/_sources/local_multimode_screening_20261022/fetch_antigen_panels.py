"""Fetch the complete independent antigen panels from a pinned primary archive."""
from pathlib import Path
import json,hashlib,urllib.request,time,concurrent.futures
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'antigen_labels';RAW=OUT/'raw';RAW.mkdir(parents=True,exist_ok=True)
OLD=ROOT.parent/'local_antigen_screening_20261021/experimental_labels/raw'
COMMIT='cfd85f3fefa652f0f00d3fd1e58cb2b4256f5ce2'
tree=json.loads((OLD/'github_tree.json').read_text(encoding='utf8'))['tree']
paths=[]
for x in tree:
    p=x['path']
    if x['type']!='blob' or not p.startswith('mutation_effects/mskcc/'):continue
    rel=p.removeprefix('mutation_effects/mskcc/')
    if rel in ['constants.py','data_preprocessing/avidity_curves.xlsx','data_preprocessing/infer_ec50s.py']:
        paths.append(p)
    elif '/' not in rel and p.endswith('.csv') and 'ec50_sat_mut' in p:paths.append(p)
    elif rel.startswith('pdbs/fold') and 'no_template' in p and p.endswith('.pdb'):paths.append(p)
    elif rel.startswith('results/') and p.endswith('.csv') and any(q in p for q in ['proteinmpnn','hermes_py_000']) and 'with_relaxation' not in p:
        paths.append(p)

def fetch(p):
    url=f'https://raw.githubusercontent.com/StatPhysBio/tcr_antigen_design/{COMMIT}/{p}'
    relative=p.removeprefix('mutation_effects/mskcc/')
    if relative.startswith('results/'):
        out=RAW/'static_scores'/(hashlib.sha256(p.encode()).hexdigest()[:20]+'.csv')
    else:out=RAW/relative
    out.parent.mkdir(parents=True,exist_ok=True)
    if not out.exists():
        for retry in range(3):
            try:
                req=urllib.request.Request(url,headers={'User-Agent':'FlexPPI-local-experiment/1.0'})
                out.write_bytes(urllib.request.urlopen(req,timeout=40).read());break
            except Exception:
                if retry==2:raise
                time.sleep(1)
    data=out.read_bytes()
    return {'repository_path':p,'path':str(out),'url':url,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}

if __name__=='__main__':
    records=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for r in pool.map(fetch,paths):records.append(r)
    (OUT/'downloads.json').write_text(json.dumps({'commit':COMMIT,'files':records},ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({'files':len(records),'MB':sum(x['bytes'] for x in records)/1e6,'root':str(OUT)},ensure_ascii=False))
