from retrieve_sources import *
from concurrent.futures import ThreadPoolExecutor
import xml.etree.ElementTree as ET
import zipfile
j=json.loads((RAW/'pettmann_elife_api.json').read_text())
urls={}
def walk(x):
    if isinstance(x,dict):
        if isinstance(x.get('uri'),str) and x['uri'].split('/')[-1] in {'elife-67092-fig1-data2-v3.csv','elife-67092-fig1-data3-v3.csv','elife-67092-supp1-v3.csv','elife-67092-supp2-v3.csv'}:urls[x['uri']]=x['uri'].split('/')[-1]
        for v in x.values():walk(v)
    elif isinstance(x,list):
        for v in x:walk(v)
walk(j)
urls['https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12557503/supplementaryFiles']='pnas2025_supplementaryFiles.zip'
urls['https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2862301/supplementaryFiles']='aleksic2010_supplementaryFiles.zip'
records=[]
def fetch(item):
    url,name=item
    try:
        p=download(url,name)
        valid=not (b'<html' in p.read_bytes()[:1000].lower())
        return {'url':url,'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'valid_non_html':valid}
    except Exception as e:return {'url':url,'error':repr(e)}
with ThreadPoolExecutor(max_workers=6) as pool:records=list(pool.map(fetch,urls.items()))
(ROOT/'primary_downloads.json').write_text(json.dumps(records,indent=2),encoding='utf8')
for rec in records:
    print(rec)
    if rec.get('path') and rec['path'].endswith('.csv'):
        print(Path(rec['path']).read_text(encoding='utf-8-sig')[:7000])
    if rec.get('path') and zipfile.is_zipfile(rec['path']):
        with zipfile.ZipFile(rec['path']) as z:
            print('ZIP FILES:',z.namelist())
            for info in z.infolist():
                target=(RAW/'primary_supplements'/Path(rec['path']).stem/info.filename).resolve()
                assert target.is_relative_to((RAW/'primary_supplements').resolve())
                if not info.is_dir():target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(z.read(info))
