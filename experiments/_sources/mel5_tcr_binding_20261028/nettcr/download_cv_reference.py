from pathlib import Path
import json,urllib.request,hashlib
root=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\nettcr')
a=json.loads((root/'download_audit.json').read_text())
for name in ['models/nettcr_2_2_pan/cv_pred_df.csv','src/train_nettcr_2_2_pan.py']:
 u=f"https://raw.githubusercontent.com/mnielLab/NetTCR-2.2/{a['commit']}/{name}"
 with urllib.request.urlopen(u,timeout=60) as r:b=r.read()
 p=root/'upstream'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
 a['files'].append({'path':name,'source':u,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
(root/'download_audit.json').write_text(json.dumps(a,indent=2))
print('downloaded CV reference and training script')
