from pathlib import Path
from PIL import Image
from pypdf import PdfReader
import hashlib,json
tmp=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
def key(path):
    im=Image.open(path).convert('RGB')
    return hashlib.sha256(im.crop((0,0,im.width,im.height-120)).tobytes()).hexdigest()
old={key(p):p.name for p in (tmp/'render').glob('page-*.png')}
items=[]
page_count=len(PdfReader(tmp/'v16_word.pdf').pages)
for path in sorted((tmp/'render_v16').glob('page-*.png'),key=lambda p:int(p.stem.split('-')[1])):
    if int(path.stem.split('-')[1])>page_count:continue
    match=old.get(key(path))
    items.append(dict(page=int(path.stem.split('-')[1]),prior_matching_content=match))
(tmp/'v16_pixel_comparison.json').write_text(json.dumps(items,indent=2),encoding='utf-8')
print(json.dumps(items))
