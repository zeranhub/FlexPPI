"""Extract existing scientific artwork without recomputing any data."""
import base64
import json
import math
import re
from pathlib import Path

import fitz
import numpy as np
from lxml import etree as ET

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'output/flexppi_v3_overall_capability_dataset_contours'
BUILD=ROOT/'output/rmsf_editable_ppt_build'
BUILD.mkdir(exist_ok=True)
(BUILD/'media').mkdir(exist_ok=True)
SVG=SOURCE/'fig_A_rmsf_dataset_contours_enlarged.svg'
tree=ET.parse(str(SVG)); root=tree.getroot()
ids={e.get('id'):e for e in root.iter() if e.get('id')}


def transform(text):
    m=np.eye(3)
    for name,values in re.findall(r'(\w+)\(([^)]*)\)',text or ''):
        v=[float(x) for x in re.findall(r'[-+]?(?:\d*\.\d+|\d+)(?:e[-+]?\d+)?',values)]
        t=np.eye(3)
        if name=='translate': t[:2,2]=[v[0],v[1] if len(v)>1 else 0]
        elif name=='scale': t[0,0]=v[0]; t[1,1]=v[1] if len(v)>1 else v[0]
        elif name=='rotate':
            a=math.radians(v[0]); t[:2,:2]=[[math.cos(a),-math.sin(a)],[math.sin(a),math.cos(a)]]
            if len(v)>1:
                p=np.array(v[1:3]); t[:2,2]=p-t[:2,:2]@p
        else: raise ValueError(name)
        m=m@t
    return m


def path_commands(d,m):
    tokens=re.findall(r'[MLQCZmlqcz]|[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?',d)
    i=0; commands=[]; current=np.zeros(2); start=current.copy()
    def point(p): return (m@np.r_[p,1])[:2].tolist()
    while i<len(tokens):
        cmd=tokens[i].upper(); i+=1
        if cmd=='Z': commands.append({'close':{}}); current=start.copy(); continue
        count={'M':2,'L':2,'Q':4,'C':6}[cmd]
        values=np.array([float(t) for t in tokens[i:i+count]]).reshape(-1,2); i+=count
        if cmd in ['M','L']:
            current=values[0]; commands.append({('moveTo' if cmd=='M' else 'lineTo'):dict(zip(['x','y'],point(current)))})
            if cmd=='M': start=current.copy()
        else:
            p0=current.copy()
            for t in np.linspace(0,1,13)[1:]:
                if cmd=='Q': p=(1-t)**2*p0+2*(1-t)*t*values[0]+t*t*values[1]
                else: p=(1-t)**3*p0+3*(1-t)**2*t*values[0]+3*(1-t)*t*t*values[1]+t**3*values[2]
                commands.append({'lineTo':dict(zip(['x','y'],point(p)))})
            current=values[-1]
    return commands


items=[]
def walk(e,m=None,style=None,group=''):
    if not isinstance(e.tag,str): return
    tag=ET.QName(e).localname
    if tag in ['defs','metadata','clipPath','text','tspan']: return
    m=(np.eye(3) if m is None else m)@transform(e.get('transform'))
    style=dict(style or {})
    style.update(dict(x.strip().split(':',1) for x in e.get('style','').split(';') if ':' in x))
    style={k:v.strip() for k,v in style.items()}
    group=e.get('id',group)
    if tag=='use':
        offset=np.eye(3); offset[:2,2]=[float(e.get('x',0)),float(e.get('y',0))]
        walk(ids[e.get('{http://www.w3.org/1999/xlink}href')[1:]],m@offset,style,group)
    elif tag=='path':
        commands=path_commands(e.get('d'),m)
        items.append(dict(kind='path',name=f'{group}_{len(items)}',commands=commands,style=style))
    elif tag=='image':
        data=base64.b64decode(e.get('{http://www.w3.org/1999/xlink}href').split(',',1)[1])
        filename=f'media/{group}.png'; (BUILD/filename).write_bytes(data)
        x,y,w,h=[float(e.get(k)) for k in ['x','y','width','height']]
        p=(m@np.array([[x,y,1],[x+w,y+h,1]]).T).T[:,:2]
        items.append(dict(kind='image',name=group,file=filename,bounds=[*p.min(axis=0),*np.abs(p[1]-p[0])],flip_vertical=bool(m[1,1]<0)))
    else:
        for child in e: walk(child,m,style,group)

walk(root)
pdf=fitz.open(SOURCE/'fig_A_rmsf_dataset_contours_enlarged.pdf')
for block in pdf[0].get_text('dict')['blocks']:
    if block['type']!=0: continue
    for line in block['lines']:
        for span in line['spans']:
            x,y=span['origin']; x0,y0,x1,y1=span['bbox']; size=span['size']
            angle=math.degrees(math.atan2(line['dir'][1],line['dir'][0]))
            # Text baseline and measured span width come from the original PDF.
            width=(x1-x0 if abs(angle)<1 else y1-y0)+size*.18
            height=size*1.2
            local_center=np.array([width/2,-size*.906+height/2])
            a=math.radians(angle); rotation=np.array([[math.cos(a),-math.sin(a)],[math.sin(a),math.cos(a)]])
            center=np.array([x,y])+rotation@local_center
            items.append(dict(kind='text',name=f'text_{len(items)}',text=span['text'],
                              bounds=[float(center[0]-width/2),float(center[1]-height/2),width,height],
                              rotation=angle,size=size,bold='Bold' in span['font'],font='Arial'))

model=dict(source_svg=str(SVG),source_size=[311.811024,317.480315],slide_size_points=[960,540],
           figure_height_points=520,items=items)
(BUILD/'figure_primitives.json').write_text(json.dumps(model,indent=2),encoding='utf-8')
print(json.dumps({'primitives':len(items),'images':sum(i['kind']=='image' for i in items),'texts':sum(i['kind']=='text' for i in items),'build':str(BUILD)}))
