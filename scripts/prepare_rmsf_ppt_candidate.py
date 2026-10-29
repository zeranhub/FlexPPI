"""Preserve SVG dash lengths in native OOXML lines and package the source project."""
import json
import re
import shutil
import zipfile
from pathlib import Path
from lxml import etree as ET

ROOT=Path(__file__).resolve().parents[1]
BUILD=ROOT/'output/rmsf_editable_ppt_build'
model=json.loads((BUILD/'figure_primitives.json').read_text())
ns={'a':'http://schemas.openxmlformats.org/drawingml/2006/main','p':'http://schemas.openxmlformats.org/presentationml/2006/main'}
styles={item['name']:item['style'] for item in model['items'] if item['kind']=='path'}
with zipfile.ZipFile(BUILD/'candidate.pptx') as zin, zipfile.ZipFile(BUILD/'candidate_exact.pptx','w',zipfile.ZIP_DEFLATED) as zout:
    for info in zin.infolist():
        data=zin.read(info.filename)
        if info.filename=='ppt/slides/slide1.xml':
            root=ET.fromstring(data)
            for shape in root.findall('.//p:sp',ns):
                identity=shape.find('p:nvSpPr/p:cNvPr',ns)
                style=styles.get(identity.get('name'),{})
                if 'stroke-dasharray' not in style: continue
                line=shape.find('p:spPr/a:ln',ns)
                for old in list(line):
                    if ET.QName(old).localname in ['prstDash','custDash']: line.remove(old)
                dash=ET.SubElement(line,'{'+ns['a']+'}custDash')
                lengths=[float(v) for v in style['stroke-dasharray'].split(',')]
                width=float(style.get('stroke-width',1))
                for d,space in zip(lengths[::2],lengths[1::2]):
                    ET.SubElement(dash,'{'+ns['a']+'}ds',d=str(round(d/width*100000)),sp=str(round(space/width*100000)))
            data=ET.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
        zout.writestr(info,data)

# JSON is a YAML 1.2 subset. The standalone PPTD source keeps all dependencies
# together; the existing multi-slide deck itself is edited natively in PowerPoint.
project=ROOT/'output/rmsf_editable_A_source'
(project/'pages').mkdir(parents=True,exist_ok=True)
(project/'media').mkdir(exist_ok=True)
scale=520/model['source_size'][1]; dx=(960-model['source_size'][0]*scale)/2; dy=10
elements=[]
for item in model['items']:
    e={'elementId':item['name']}
    if item['kind']=='path':
        pts=[c[k] for c in item['commands'] for k in ['moveTo','lineTo'] if k in c]
        x=min(p['x'] for p in pts); y=min(p['y'] for p in pts)
        w=max(.01,max(p['x'] for p in pts)-x); h=max(.01,max(p['y'] for p in pts)-y)
        d=[]
        for c in item['commands']:
            if 'close' in c: d.append('Z')
            else:
                k=next(iter(c)); p=c[k]; d.append(('M' if k=='moveTo' else 'L')+f"{p['x']-x:.6f},{p['y']-y:.6f}")
        style=item['style']; e.update(elementType='shape',shapeName='custom',viewBox=[w,h],path=' '.join(d),bounds=[dx+x*scale,dy+y*scale,w*scale,h*scale])
        if style.get('fill','black')!='none': e['fill']={'type':'solid','color':style.get('fill','#000000')}
        if style.get('stroke'): e['border']={'style':'dash' if style.get('stroke-dasharray') else 'solid','width':float(style.get('stroke-width',1))*scale,'color':style['stroke']}
        if style.get('opacity'): e['opacity']=float(style['opacity'])
    else:
        x,y,w,h=item['bounds']; e['bounds']=[dx+x*scale,dy+y*scale,w*scale,h*scale]
        if item['kind']=='text':
            e.update(elementType='text',rotation=item['rotation'],content={'text':item['text'],'fontFamily':'Arial','fontSize':item['size']*scale,'bold':item['bold'],'align':['left','top'],'color':'#000000'})
        else:
            e.update(elementType='image',src=item['file'],flip=[False,item['flip_vertical']],fit={'mode':'contain'})
            shutil.copy2(BUILD/item['file'],project/item['file'])
    elements.append(e)
(project/'editable_A.pptd').write_text(json.dumps({'version':'v2','title':'FlexPPI RMSF dataset contours','size':[960,540],'pages':['pages/A.page']},indent=2),encoding='utf-8')
(project/'pages/A.page').write_text(json.dumps({'background':{'type':'solid','color':'#FFFFFF'},'elements':elements},indent=2),encoding='utf-8')
shutil.copy2(BUILD/'figure_primitives.json',project/'exact_source_primitives.json')
print(json.dumps({'candidate':str(BUILD/'candidate_exact.pptx'),'pptd_source':str(project)}))
