import requests,json,pathlib
r=pathlib.Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\imrex')
for url,name in [('https://zenodo.org/api/records/3973547','zenodo_record.json'),('https://codeload.github.com/pmoris/ImRex/zip/d3975b744c030c4bc67c942c82a32b1c83857c4e','upstream.zip')]:
    try:
        x=requests.get(url,timeout=45);x.raise_for_status();(r/name).write_bytes(x.content)
        print(name,x.status_code,len(x.content))
        if name.endswith('.json'): print([(y['key'],y.get('size'),y.get('links')) for y in x.json().get('files',[])])
    except Exception as e:print(type(e).__name__,str(e))
