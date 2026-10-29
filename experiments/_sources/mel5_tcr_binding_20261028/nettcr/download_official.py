import subprocess,json,pathlib,concurrent.futures,hashlib,urllib.request,time
root=pathlib.Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\nettcr')
repo=root/'upstream';repo.mkdir(parents=True,exist_ok=True)
meta=json.loads(subprocess.check_output(['gh','api','repos/mnielLab/NetTCR-2.2/commits/main']))
commit=meta['sha'];(root/'upstream_commit.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
tree=json.loads(subprocess.check_output(['gh','api',f'repos/mnielLab/NetTCR-2.2/git/trees/{commit}?recursive=1']))
(root/'upstream_tree.json').write_text(json.dumps(tree,indent=2),encoding='utf-8')
files=['README.md','environment.yml','src/predict.py','src/keras_utils.py','src/make_webserver_prediction.py','src/nettcr_archs.py','data/small_example.csv','data/nettcr_2_2_full_dataset.csv','data/nettcr_2_2_limited_dataset.csv','academic_software_license_agreement.pdf']
files += [x['path'] for x in tree['tree'] if x['type']=='blob' and x['path'].startswith('models/nettcr_2_2_pan/checkpoint/')]
records=[]
def fetch(name):
 p=repo/name;p.parent.mkdir(parents=True,exist_ok=True)
 url=f'https://raw.githubusercontent.com/mnielLab/NetTCR-2.2/{commit}/{name}'
 for attempt in range(3):
  try:
   with urllib.request.urlopen(url,timeout=60) as response:b=response.read()
   p.write_bytes(b);return {'path':name,'source':url,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
  except Exception:
   if attempt==2:raise
   time.sleep(2)
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
 for record in pool.map(fetch,files):records.append(record)
(root/'download_audit.json').write_text(json.dumps({'repository':'https://github.com/mnielLab/NetTCR-2.2','commit':commit,'files':records},indent=2),encoding='utf-8')
print(json.dumps({'commit':commit,'files':len(records),'bytes':sum(x['bytes'] for x in records)},indent=2))
