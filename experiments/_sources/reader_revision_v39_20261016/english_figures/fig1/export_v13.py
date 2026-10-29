from pathlib import Path
import os
import subprocess
import json
import fitz
from PIL import Image

HERE = Path(__file__).resolve().parent
EXE = os.environ.get('FLEXPPI_DRAWIO_EXE','C:/Program Files/draw.io/draw.io.exe')
SOURCE = HERE / 'FlexPPI_main_v13.drawio'
commands = []
for ext in ['png','svg','pdf']:
    output = HERE / ('FlexPPI_main_v13_native.pdf' if ext=='pdf' else f'FlexPPI_main_v13.{ext}')
    args = [EXE,'-x','-f',ext,'--crop','-o',str(output),str(SOURCE)]
    if ext == 'svg':
        args.insert(4,'--embed-diagram')
    r = subprocess.run(args,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=120)
    assert r.returncode==0, r.stdout+r.stderr
    commands.append(dict(command=args,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr))
args = [EXE,'-x','-f','png','--crop','-s','1.35','-o',str(HERE/'FlexPPI_main_v13_300dpi.png'),str(SOURCE)]
r = subprocess.run(args,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=120)
assert r.returncode == 0, r.stdout+r.stderr
commands.append(dict(command=args,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr))
im = Image.open(HERE/'FlexPPI_main_v13_300dpi.png')
im.save(HERE/'FlexPPI_main_v13_300dpi.png',dpi=(300,300))
src = fitz.open(HERE/'FlexPPI_main_v13_native.pdf')
assert len(src)==1
w = 183/25.4*72
h = w*src[0].rect.height/src[0].rect.width
dest = fitz.open();p=dest.new_page(width=w,height=h)
p.show_pdf_page(p.rect,src,0)
dest.save(HERE/'FlexPPI_main_v13_183mm.pdf',no_new_id=True,garbage=4,deflate=True)
pdf = fitz.open(HERE/'FlexPPI_main_v13_183mm.pdf')
pdf[0].get_pixmap(matrix=fitz.Matrix(2.1,2.1)).save(HERE/'pdf_page_v13.png')
preview = Image.open(HERE/'FlexPPI_main_v13.png')
# Native output includes a one-pixel margin; crops are previews of fixed panel regions.
for name, rect in [('panel_a_preview',(20,85,1580,667)),('panel_b_preview',(26,694,586,1250)),('panel_c_preview',(593,694,1100,1250))]:
    preview.crop(rect).save(HERE/(name+'.png'))
# Export the same editable panel at a useful inspection size.
for ext in ['png','svg','pdf']:
    args = [EXE,'-x','-f',ext,'--crop','-o',str(HERE/f'FlexPPI_local_frames_v13.{ext}'),str(HERE/'FlexPPI_local_frames_v13.drawio')]
    if ext == 'png':
        args[4:4] = ['-s','2']
    elif ext == 'svg':
        args.insert(4,'--embed-diagram')
    r = subprocess.run(args,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=120)
    assert r.returncode == 0, r.stdout+r.stderr
    commands.append(dict(command=args,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr))
spans=[s for b in pdf[0].get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if s['text'].strip()]
report=dict(commands=commands,pdf_width_mm=183,pdf_height_mm=h/72*25.4,
            min_pdf_font_pt=min(s['size'] for s in spans),png_size=list(preview.size),
            vector_text=True)
(HERE/'export_receipt.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='commands'}))
