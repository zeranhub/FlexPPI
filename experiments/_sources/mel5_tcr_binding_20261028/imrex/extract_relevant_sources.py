from pathlib import Path
import zipfile,hashlib,json
r=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\imrex')
records=[]
with zipfile.ZipFile(r/'upstream.zip') as z:
 for item in z.infolist():
  name=item.filename.split('/',1)[-1]
  if not name or item.is_dir():continue
  if not (name.startswith('src/') or name.startswith('models/pretrained/') or name in ['README.md','LICENSE','environment.yml','requirements.txt','setup.py']): continue
  p=r/'upstream'/name;p_long=Path('\\\\?\\'+str(p));p_long.parent.mkdir(parents=True,exist_ok=True);blob=z.read(item);p_long.write_bytes(blob)
  records.append({'path':name,'sha256':hashlib.sha256(blob).hexdigest(),'bytes':len(blob)})
  if name.endswith('.h5'):
   tag='downsample' if 'trbmhcidown-' in name else 'full'
   p2=r/'weights'/(tag+'.h5');p2.parent.mkdir(exist_ok=True);p2.write_bytes(blob)
  if name.endswith('train.csv'):
   tag='downsample' if 'trbmhcidown-' in name else 'full'
   p2=r/'training'/('train_'+tag+'.csv');p2.parent.mkdir(exist_ok=True);p2.write_bytes(blob)
(r/'source_file_hashes.json').write_text(json.dumps({'commit':'d3975b744c030c4bc67c942c82a32b1c83857c4e','archive_sha256':hashlib.sha256((r/'upstream.zip').read_bytes()).hexdigest(),'files':records},indent=2),encoding='utf-8')
print(len(records),'source/model files extracted')
for name in ['src/bio/feature_builder.py','src/processing/image_padding.py','src/neural/evaluation.py','src/bio/peptide.py']:
 print('\nFILE',name,'\n',(r/'upstream'/name).read_text())
