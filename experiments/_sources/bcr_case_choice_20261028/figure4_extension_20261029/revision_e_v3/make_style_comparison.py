from pathlib import Path
from PIL import Image
root=Path(__file__).resolve().parent
paper=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
mm=lambda n:round(n/25.4*600)
old=Image.open(paper/'figure4_v124_20261027_紧凑蓝橙灰.png').convert('RGB')
a=old.crop((0,0,mm(116),mm(51.5)))
e=Image.open(root/'figure4e_redraw_v3.png').convert('RGB')
canvas=Image.new('RGB',(a.width+mm(4)+e.width,max(a.height,e.height)),'white')
canvas.paste(a,(0,0));canvas.paste(e,(a.width+mm(4),0))
canvas.save(paper/'figure4a_e_20261029_风格对照_v3.png',dpi=(600,600))
canvas.save(root/'figure4a_e_style_comparison_v3.png',dpi=(600,600))
