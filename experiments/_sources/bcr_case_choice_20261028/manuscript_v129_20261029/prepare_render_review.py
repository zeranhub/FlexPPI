from pathlib import Path
import json, hashlib, re
from PIL import Image, ImageChops, ImageDraw, ImageFont
from pypdf import PdfReader
R=Path(__file__).resolve().parent
OLD=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\manuscript\v128\render')
NEW=R/'render_final'
reader=PdfReader(next(NEW.glob('*.pdf')))
old_pages={}
for p in OLD.glob('page-*.png'):
    image=Image.open(p).convert('RGB')
    crop=image.crop((0,0,image.width,image.height-110))
    digest=hashlib.sha256(crop.tobytes()).hexdigest()
    old_pages.setdefault(digest,[]).append(int(p.stem.split('-')[1]))
results=[]
for i,page in enumerate(reader.pages,1):
    p=NEW/f'page-{i}.png'
    image=Image.open(p).convert('RGB')
    crop=image.crop((0,0,image.width,image.height-110))
    matches=old_pages.get(hashlib.sha256(crop.tobytes()).hexdigest(),[])
    results.append({'page':i,'body_pixel_identical_to_previous_pages':matches,'text':page.extract_text(),
                    'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(NEW/'page_text.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
for row in results:
    tags=[x for x in ['B5209B','3.4','R101','[37]','H B5209B','输入与参考','MD参考与采样范围','完整区域结果','解释范围与复现'] if x in row['text']]
    print(row['page'],'prior body identical',row['body_pixel_identical_to_previous_pages'],'tags',tags)
font=ImageFont.truetype(r'C:\Windows\Fonts\arial.ttf',22)
pages=list(range(1,len(results)+1))
for start in range(0,len(pages),6):
    group=pages[start:start+6]
    sheet=Image.new('RGB',(1800,2650),'#e0e0e0')
    d=ImageDraw.Draw(sheet)
    for j,i in enumerate(group):
        im=Image.open(NEW/f'page-{i}.png').convert('RGB')
        im.thumbnail((880,835))
        x=(j%2)*900+10;y=(j//2)*880+35
        sheet.paste(im,(x,y))
        d.text((x,y-30),f'Page {i}',font=font,fill='black')
    sheet.save(NEW/f'contact-{start+1:02d}-{start+len(group):02d}.png')
print('pages',len(results))
