from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import urllib.request,json
import pandas as pd
O=Path(__file__).resolve().parent
f=pd.read_csv(O/'shortlist_region_screen.csv')
q=f[~f.sample_id.str.startswith('itrace_')].drop_duplicates('pdb_id').head(100).copy()
def one(pid):
 try:
  u='https://data.rcsb.org/rest/v1/core/entry/'+pid
  d=json.load(urllib.request.urlopen(u,timeout=20))
  return dict(pdb_id=pid,title=d.get('struct',{}).get('title',''),paper_title=d.get('rcsb_primary_citation',{}).get('title',''),doi=d.get('rcsb_primary_citation',{}).get('pdbx_database_id_PubMed',''),classification=d.get('struct_keywords',{}).get('pdbx_keywords',''),entity_ids=d.get('rcsb_entry_container_identifiers',{}).get('polymer_entity_ids',[]))
 except Exception as e:return dict(pdb_id=pid,error=str(e))
with ThreadPoolExecutor(max_workers=8) as ex:rows=list(ex.map(one,q.pdb_id))
pd.DataFrame(rows).to_csv(O/'top100_rcsb_identity.csv',index=False)
q=q.merge(pd.DataFrame(rows),on='pdb_id',validate='one_to_one');q.to_csv(O/'top100_regions_with_identity.csv',index=False)
print(q[['sample_id','pdb_id','region','n','full_mae','static47_sasa_mae','title']].to_string(index=False),flush=True)
