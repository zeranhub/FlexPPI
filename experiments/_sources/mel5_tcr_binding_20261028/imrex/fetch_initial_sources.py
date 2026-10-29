import subprocess,json,base64,pathlib,urllib.request,zipfile
r=pathlib.Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\imrex');r.mkdir(parents=True,exist_ok=True)
for name in ['README.md','models/pretrained/README.md','src/scripts/predict/predict.py','requirements.txt']:
    data=json.loads(subprocess.check_output(['gh','api','repos/pmoris/ImRex/contents/'+name]))
    p=r/'upstream'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(base64.b64decode(data['content']))
    print('\nFILE',name,'\n',p.read_text()[:12000])
