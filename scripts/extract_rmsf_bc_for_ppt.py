"""Extract panels B/C from the original SVG/PDF, preserving geometry and data."""
import ast
import base64
import json
import math
from pathlib import Path
import re
import shutil

import fitz
from lxml import etree as ET
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'output/flexppi_v3_overall_capability'
BUILD=ROOT/'output/rmsf_editable_BC_ppt_build'
BUILD.mkdir(exist_ok=True)
(BUILD/'media').mkdir(exist_ok=True)
svg=SOURCE/'fig_v3_overall_rmsf_capability.svg'
root=ET.parse(str(svg)).getroot()
ids={e.get('id'):e for e in root.iter() if e.get('id')}
ns={'s':root.nsmap[None]}
pdf=fitz.open(SOURCE/'fig_v3_overall_rmsf_capability.pdf')

# Reuse the verified SVG parser functions without executing the A-only exporter.
parser_source=ast.parse((ROOT/'scripts/extract_rmsf_svg_for_ppt.py').read_text(encoding='utf-8'))
functions=[n for n in parser_source.body if isinstance(n,ast.FunctionDef)]
exec(compile(ast.Module(body=functions,type_ignores=[]),'shared_svg_parser','exec'))

all_items=[]
panel_summaries=[]
for letter,group_id,crop,offset in [
    ('B','axes_2',[280.,40.,512.,229.],0.),
    ('C','axes_3',[12.,252.,244.,441.],252.),
]:
    x_crop,y_crop,x_right,y_bottom=crop
    items=[]
    walk(ids[group_id])
    for block in pdf[0].get_text('dict')['blocks']:
        if block['type']!=0: continue
        for line in block['lines']:
            for span in line['spans']:
                x,y=span['origin']; x0,y0,x1,y1=span['bbox']; size=span['size']
                if not (x_crop<=(x0+x1)/2<=x_right and y_crop<=(y0+y1)/2<=y_bottom): continue
                angle=math.degrees(math.atan2(line['dir'][1],line['dir'][0]))
                width=(x1-x0 if abs(angle)<1 else y1-y0)+size*.18
                height=size*1.2
                local_center=np.array([width/2,-size*.906+height/2])
                a=math.radians(angle); rotation=np.array([[math.cos(a),-math.sin(a)],[math.sin(a),math.cos(a)]])
                center=np.array([x,y])+rotation@local_center
                items.append(dict(kind='text',name=f'text_{len(items)}',text=span['text'],
                                  bounds=[float(center[0]-width/2),float(center[1]-height/2),width,height],
                                  rotation=angle,size=size,bold='Bold' in span['font'],font='Arial'))
    for item in items:
        item['name']=letter+'_'+item['name']; item['panel']=letter
        if item['kind']=='path':
            for command in item['commands']:
                for key in ['moveTo','lineTo']:
                    if key in command:
                        command[key]['x']+=offset-x_crop
                        command[key]['y']-=y_crop
        else:
            item['bounds'][0]+=offset-x_crop
            item['bounds'][1]-=y_crop
    texts=[i['text'] for i in items if i['kind']=='text']
    for value in [letter,'ATLAS','DynaRepo','ITrace','n = 1,929','n = 895','n = 244']:
        if value not in texts: raise ValueError(f'{letter}: missing {value}')
    if sum(i['kind']=='image' for i in items)!=1: raise ValueError('Unexpected point layer count')
    panel_summaries.append(dict(panel=letter,n_objects=len(items),n_texts=len(texts),n_images=1,source_crop=crop))
    all_items.extend(items)
model=dict(source_svg=str(svg),source_size=[484.,189.],slide_size_points=[960,540],
           figure_height_points=900/484*189,items=all_items,panels=panel_summaries)
(BUILD/'figure_primitives.json').write_text(json.dumps(model,indent=2),encoding='utf-8')

# Retain the working native-slide construction and typography used for panel A.
builder=(ROOT/'output/rmsf_editable_ppt_build/build_editable_a.mjs').read_text(encoding='utf-8')
builder=builder.replace("color(s.fill||'#000000',s.opacity)","color(s.fill||'#000000',s['fill-opacity']||s.opacity)")
builder=builder.replace("color(s.stroke,s.opacity)","color(s.stroke,s['stroke-opacity']||s.opacity)")
builder=builder.replace("alt:'Exact source RMSF density layer'","alt:item.name")
start=builder.index('slide.speakerNotes.textFrame.setText(')
end=builder.index('\nconst candidatePath',start)
notes='FlexPPI-v3 panels B and C. All 3,068 held-out protein samples: ATLAS 1,929, DynaRepo 895, ITrace 244. B: per-protein Pearson r. C: per-protein MAE in angstroms on a logarithmic axis. Each point represents one sample. White marker: median. Black bar: interquartile range. Source: validated frozen epoch-40 out-of-fold predictions. Exact source point layers and vector geometry reused.'
builder=builder[:start]+'slide.speakerNotes.textFrame.setText('+json.dumps(notes)+');'+builder[end:]
(BUILD/'build_editable_bc.mjs').write_text(builder,encoding='utf-8')

prepare=(ROOT/'scripts/prepare_rmsf_ppt_candidate.py').read_text(encoding='utf-8')
prepare=prepare.replace("ROOT=Path(__file__).resolve().parents[1]",f'ROOT=Path({str(ROOT)!r})')
prepare=prepare.replace('rmsf_editable_ppt_build','rmsf_editable_BC_ppt_build').replace('rmsf_editable_A_source','rmsf_editable_BC_source')
prepare=prepare.replace('editable_A','editable_BC').replace('pages/A.page','pages/BC.page').replace("pages/A.page","pages/BC.page")
prepare=prepare.replace('scale=520/model[\'source_size\'][1]; dx=(960-model[\'source_size\'][0]*scale)/2; dy=10',
                        "scale=model['figure_height_points']/model['source_size'][1]; dx=(960-model['source_size'][0]*scale)/2; dy=(540-model['figure_height_points'])/2")
prepare=prepare.replace("if style.get('opacity'): e['opacity']=float(style['opacity'])", "if style.get('opacity') or style.get('fill-opacity'): e['opacity']=float(style.get('opacity',style.get('fill-opacity')))" )
prepare=prepare.replace('FlexPPI RMSF dataset contours','FlexPPI RMSF panels B and C')
(BUILD/'prepare_candidate.py').write_text(prepare,encoding='utf-8')
finalizer=(ROOT/'output/rmsf_editable_ppt_build/finalize_editable_a.mjs').read_text(encoding='utf-8')
finalizer=finalizer.replace('rmsf_editable_A_source','rmsf_editable_BC_source').replace('editable_A','editable_BC')
(BUILD/'finalize_editable_bc.mjs').write_text(finalizer,encoding='utf-8')
print(json.dumps({'build':str(BUILD),'panels':panel_summaries,'total_objects':len(all_items)}))
