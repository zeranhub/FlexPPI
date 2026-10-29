from pathlib import Path
import json
from PIL import Image
OUT=Path(__file__).resolve().parent
m=json.loads((OUT/'render_manifest.json').read_text())
w=max(a['trim']['width'] for a in m['assets'].values());h=max(a['trim']['height'] for a in m['assets'].values())
for name,a in m['assets'].items():
 im=Image.open(OUT/a['file']).convert('RGBA');b=a['alpha_bbox'];cx=(b[0]+b[2])/2;cy=(b[1]+b[3])/2
 x=round(cx-w/2);y=round(cy-h/2);crop=[x,y,x+w,y+h]
 fn=name+'_matched_scale_trim.png';im.crop(crop).save(OUT/fn)
 a['matched_scale_trim']={'file':fn,'width':w,'height':h,'crop_box':crop,'landmarks':{k:[p[0]-x,p[1]-y] for k,p in a['landmarks'].items()}}
(OUT/'render_manifest.json').write_text(json.dumps(m,indent=2),encoding='utf-8')
print('Equal pixel and physical scale crops:',w,h)