from retrieve_sources import *
from concurrent.futures import ThreadPoolExecutor
from html import unescape

meta=json.loads((ROOT/'source_candidates.json').read_text())
commit=meta['commit']
paths=[x['path'] for x in meta['files'] if x['path'].endswith(('.csv','constants.py','split_csv.py')) and '/pdbs/' not in x['path'] and 'output.csv' not in x['path']]
paths += [x['path'] for x in meta['files'] if '/pdbs/' in x['path'] and x['path'].endswith('.pdb') and '/fold_' not in x['path']]
def fetch(p):
    url=f'https://raw.githubusercontent.com/{REPO}/{commit}/{p}'
    out=download(url,'github/'+p)
    return {'url':url,'path':str(out),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'bytes':out.stat().st_size}
with ThreadPoolExecutor(max_workers=6) as pool:results=list(pool.map(fetch,paths))
(ROOT/'downloaded_panel_files.json').write_text(json.dumps(results,indent=2),encoding='utf8')
html=(RAW/'pnas_PMC12557503.html').read_text(encoding='utf8')
links=[unescape(x) for x in re.findall(r'href="([^"]+)"',html)]
print('Download count:',len(results))
print('Supplement/PDF links:')
print('\n'.join(x for x in links if any(s in x for s in ['supp','2504783','pnas.','ref47','ref48'])))
for p in ['mutation_effects/nyeso/nyeso_peptide_kd.csv','mutation_effects/nyeso/nyeso_peptide_kd_closest.csv','mutation_effects/tax/tax_peptide_kd.csv','mutation_effects/tax/tax_peptide_kd_closest.csv','mutation_effects/nyeso/split_csv.py','mutation_effects/tax/split_csv.py']:
    print('\nFILE',p);print((RAW/'github'/p).read_text())
