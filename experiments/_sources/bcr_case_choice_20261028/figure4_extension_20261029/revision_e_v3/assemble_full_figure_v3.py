"""Append the antibody preview to unchanged TCR vector panels a-d.

No figure is embedded in Word. This creates a separate review artifact only.
"""
from pathlib import Path
from copy import deepcopy
import hashlib, json, re, shutil
from pypdf import PdfReader, PdfWriter, PageObject, Transformation
from PIL import Image
from lxml import etree

ROOT=Path(__file__).resolve().parent
PAPER=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
BASE=PAPER/'figure4_v124_20261027_紧凑蓝橙灰'
EXT=ROOT/'figure4e_g_6A77_preview_v3'
OUT=ROOT/'figure4_with_6A77_preview_v3'
GAP=7.
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

old=PdfReader(str(BASE)+'.pdf').pages[0]
new=PdfReader(str(EXT)+'.pdf').pages[0]
width=float(old.mediabox.width);height=float(old.mediabox.height)
scale=width/float(new.mediabox.width)
new_height=float(new.mediabox.height)*scale
page=PageObject.create_blank_page(width=width,height=height+GAP+new_height)
page.merge_transformed_page(new,Transformation().scale(scale))
page.merge_transformed_page(old,Transformation().translate(tx=0,ty=new_height+GAP))
writer=PdfWriter();writer.add_page(page)
writer.add_metadata({'/Title':'Figure 4 TCR and antibody dynamics preview','/Subject':'Separate review artifact; not inserted into the manuscript'})
with Path(str(OUT)+'.pdf').open('wb') as f:writer.write(f)

svg_ns='http://www.w3.org/2000/svg'
old_svg=etree.parse(str(BASE)+'.svg').getroot()
new_svg=etree.parse(str(EXT)+'.svg').getroot()
old_view=[float(v) for v in old_svg.attrib['viewBox'].split()]
new_view=[float(v) for v in new_svg.attrib['viewBox'].split()]
svg_scale=old_view[2]/new_view[2]
total_height=old_view[3]+GAP+new_view[3]*svg_scale
canvas=etree.Element('{'+svg_ns+'}svg',nsmap={None:svg_ns,'xlink':'http://www.w3.org/1999/xlink'})
canvas.set('width',f'{old_view[2]:.6f}pt');canvas.set('height',f'{total_height:.6f}pt')
canvas.set('viewBox',f'0 0 {old_view[2]:.6f} {total_height:.6f}')
back=etree.SubElement(canvas,'{'+svg_ns+'}rect',x='0',y='0',width='100%',height='100%',fill='white')
def namespace_ids(node,prefix):
    names={el.attrib['id']:prefix+el.attrib['id'] for el in node.iter() if 'id' in el.attrib}
    for el in node.iter():
        for key,value in list(el.attrib.items()):
            if key=='id':el.set(key,names[value]);continue
            if value.startswith('#') and value[1:] in names:value='#'+names[value[1:]]
            value=re.sub(r'url\(#([^)]*)\)',lambda m:'url(#'+names.get(m.group(1),m.group(1))+')',value)
            el.set(key,value)
    return node
for node,prefix,y,sc in [(old_svg,'tcr_',0,1.),(new_svg,'ab_',old_view[3]+GAP,svg_scale)]:
    node=namespace_ids(deepcopy(node),prefix)
    group=etree.SubElement(canvas,'{'+svg_ns+'}g',transform=f'translate(0 {y:.6f}) scale({sc:.9f})')
    for child in node:group.append(child)
etree.ElementTree(canvas).write(str(OUT)+'.svg',encoding='utf-8',xml_declaration=True)

upper_image=Image.open(str(BASE)+'.png')
source_dpi=round(upper_image.info.get('dpi',(600,600))[0])
upper=upper_image.convert('RGB')
lower=Image.open(str(EXT)+'.png').convert('RGB')
if lower.width!=upper.width:
    lower=lower.resize((upper.width,round(lower.height*upper.width/lower.width)),Image.Resampling.LANCZOS)
gap_pixels=round(GAP/72*source_dpi)
whole=Image.new('RGB',(upper.width,upper.height+gap_pixels+lower.height),'white')
whole.paste(upper,(0,0));whole.paste(lower,(0,upper.height+gap_pixels))
whole.save(str(OUT)+'.png',dpi=(source_dpi,source_dpi))
report={'status':'assembled_pending_visual_QA','PDF_size_mm':[width/72*25.4,(height+GAP+new_height)/72*25.4],
        'PNG_size':list(whole.size),'PNG_dpi':source_dpi,'original_TCR_panels':'unchanged vector pages/SVG children, unchanged PNG pixels',
        'pixel_check_old_TCR_equal':whole.crop((0,0,upper.width,upper.height)).tobytes()==upper.tobytes(),
        'new_PDF_scale':scale,'gap_pt':GAP,
        'input_sha256':{str(p)+'.'+s:sha(str(p)+'.'+s) for p in [BASE,EXT] for s in ['pdf','svg','png']},
        'outputs':{s:{'path':str(OUT)+'.'+s,'sha256':sha(str(OUT)+'.'+s)} for s in ['pdf','svg','png']},
        'manuscript_image_inserted':False}
(ROOT/'assembled_figure_provenance.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
for s in ['pdf','svg','png']:
    shutil.copy2(str(OUT)+'.'+s,PAPER/f'figure4_20261029_TCR与抗体合并预览_v3.{s}')
html='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>图4 TCR与抗体动态预览</title><style>body{font:16px Arial,sans-serif;color:#222;max-width:1350px;margin:28px auto;padding:0 20px}img{max-width:100%;height:auto}p{line-height:1.6}a{color:#236ab6}</style><h2>图4 TCR与抗体动态 图e重绘版预览</h2><p>图e沿用图a的渲染风格：左侧为6A77复合物，右侧为R96和R101局部接触；a–d、f和g保留。图稿尚未嵌入论文。</p><p><a href="figure4e_g_20261029_抗体动态预览_v3.svg">单独查看新增e–g</a> · <a href="figure4_20261029_TCR与抗体合并预览_v3.pdf">PDF</a> · <a href="figure4_20261029_TCR与抗体合并预览_v3.svg">可编辑SVG</a></p><img src="figure4_20261029_TCR与抗体合并预览_v3.svg"></html>'''
(PAPER/'figure4_20261029_TCR与抗体合并预览_v3.html').write_text(html,encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
