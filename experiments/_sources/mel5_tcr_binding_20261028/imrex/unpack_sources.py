from pathlib import Path
import zipfile
r=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\imrex')
with zipfile.ZipFile(r/'upstream.zip') as z:
 print('entries',len(z.namelist()))
 for x in z.infolist():
  if any(s in x.filename for s in ['models/pretrained/','src/bio/','src/neural/evaluation.py']):print(x.filename,x.file_size)
 z.extractall(r/'archive_source')
