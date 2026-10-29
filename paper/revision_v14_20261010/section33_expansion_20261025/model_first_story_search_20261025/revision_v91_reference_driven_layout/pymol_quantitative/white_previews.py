from pathlib import Path
from PIL import Image,ImageDraw,ImageFont
r=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v91_reference_driven_layout\pymol_quantitative')
for f in ['1MXE_RMSF_axis000','1MXE_Hchi_axis000','1MXE_Hchi_axis045','1MXE_Hchi_axis090','1MXE_Hchi_axis135']:
 im=Image.open(r/(f+'.png')).convert('RGBA');bg=Image.new('RGBA',im.size,'white');bg.alpha_composite(im);bg.convert('RGB').save(r/(f+'_white_preview.png'))