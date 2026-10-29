from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import urllib.request,json,hashlib
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028')
def one(pdb):
 path=R/'structures'/(pdb.lower()+'_raw.pdb')
 if not path.exists():
  data=urllib.request.urlopen('https://files.rcsb.org/download/'+pdb.upper()+'.pdb',timeout=40).read()
  assert b'ATOM  ' in data
  path.write_bytes(data)
 lines=path.read_text().splitlines();chains={}
 for l in lines:
  if l.startswith('ENDMDL'):break
  if l.startswith('ATOM  ') and l[16] in [' ','A']:
   chains.setdefault(l[21],{})[l[22:27]]=l[17:20]
 return {'pdb':pdb,'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'chains':{ch:{'n':len(rs),'first':list(rs.items())[:2],'last':list(rs.items())[-2:]} for ch,rs in chains.items()},'compound':[l.strip() for l in lines if l.startswith('COMPND')]}
with ThreadPoolExecutor(max_workers=4) as pool:out=list(pool.map(one,['1C08','5AYU','6P4A','6P4B','6P4C','6P4D']))
(R/'structure_download_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
