from pathlib import Path
import re,html
P=Path(__file__).parent/'raw'
for p in P.glob('*oa.xml'):print(p.name,p.read_text(encoding='utf8')[:400])
for p in P.glob('*pmc.html'):
    text=p.read_text(encoding='utf8')
    links=list(dict.fromkeys(html.unescape(u) for u in re.findall(r'href="([^"]+)"',text) if any(s in u.lower() for s in ['supp','mmc','nihms','xlsx','xls','csv','5e9','4ftv'])))
    print(p.name,'\n'.join(links))
    try:
        from bs4 import BeautifulSoup
        soup=BeautifulSoup(text,'html.parser')
        (p.with_suffix('.txt')).write_text(soup.get_text(' ',strip=True),encoding='utf8')
        for s in soup.find_all(['div','section'],id=re.compile(r'.*(?:sup|supp|app).*')):
            st=s.get_text(' ',strip=True)
            if len(st)<2000:print('SUPP_SECTION',s.get('id'),st)
    except ImportError:pass
