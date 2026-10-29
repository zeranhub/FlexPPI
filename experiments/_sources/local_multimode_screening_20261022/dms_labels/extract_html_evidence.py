from pathlib import Path
import sys,re,json
sys.path.insert(0,r'E:\Project\FlexPPI\experiments\_sources\completed_results_20261012\plot_runtime')
from lxml import html
P=Path(__file__).parent;R=P/'raw';allrows=[]
for name in ['PMC4938795_pmc.html','PMC5114409_pmc.html']:
    root=html.fromstring((R/name).read_bytes());art=root.xpath('//article')
    if not art:art=[root]
    node=art[0]
    (R/(name+'.txt')).write_text('\n'.join(' '.join(n.itertext()).strip() for n in node.xpath('.//p|.//h1|.//h2|.//h3|.//h4|.//figcaption')),encoding='utf8')
    links=[]
    for n in node.xpath('.//a[@href]'):
        u=n.get('href')
        if any(q in u.lower() for q in ['/bin/','supp','mmc','geo','sra','rcsb','pdb']):links.append({'text':' '.join(n.itertext()),'href':u})
    tables=[]
    for t in node.xpath('.//table'):
        data=[[' '.join(''.join(c.itertext()).split()) for c in tr.xpath('./th|./td')] for tr in t.xpath('.//tr')]
        tables.append({'rows':data})
    (P/(name.replace('_pmc.html','')+'_tables_links.json')).write_text(json.dumps({'tables':tables,'links':links},ensure_ascii=False,indent=2),encoding='utf8')
    print(name,'tables',len(tables),'links',links)
