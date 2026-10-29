from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import urllib.request,json,hashlib,zipfile,datetime,xml.etree.ElementTree as ET
ROOT=Path(__file__).parent;RAW=ROOT/'raw';RAW.mkdir(exist_ok=True)
def fetch(item):
    pmcid,kind=item
    suffix='fullText.xml' if kind=='fullTextXML' else 'supplementaryFiles.zip'
    path=RAW/(pmcid+'_'+suffix)
    url=f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/{kind}'
    req=urllib.request.Request(url,headers={'User-Agent':'FlexPPI-DMS-label-audit/1.0'})
    data=urllib.request.urlopen(req,timeout=90).read();path.write_bytes(data)
    record={'url':url,'path':str(path),'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)}
    if kind=='fullTextXML':
        article=ET.fromstring(data)
        record['title']=''.join(article.find('.//article-title').itertext())
        record['supplement_links']=[ET.tostring(x,encoding='unicode') for x in article.findall('.//supplementary-material')]
        record['doi']=[x.text for x in article.findall('.//article-id') if x.attrib.get('pub-id-type')=='doi']
    else:
        assert zipfile.is_zipfile(path)
        with zipfile.ZipFile(path) as z:
            names=z.namelist();record['files']=names
            for x in z.infolist():
                if x.is_dir():continue
                p=(RAW/pmcid/x.filename).resolve();assert p.is_relative_to((RAW/pmcid).resolve())
                p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(z.read(x))
    return record
if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=4) as pool:
        out=list(pool.map(fetch,[(p,k) for p in ['PMC4938795','PMC5114409'] for k in ['fullTextXML','supplementaryFiles']]))
    (ROOT/'initial_source_manifest.json').write_text(json.dumps({'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sources':out},ensure_ascii=False,indent=2),encoding='utf8')
    for r in out:
        print({k:v for k,v in r.items() if k!='supplement_links' and k!='files'})
        if 'files'in r:print('files',r['files'])
        if 'supplement_links'in r:print('supplements',r['supplement_links'])
