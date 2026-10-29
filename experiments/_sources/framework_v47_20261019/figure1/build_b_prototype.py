"""Native Draw.io prototype: real ribbon and true 3-D backbone inset.

The full Fig. 1 is to retain its original five-stage top row + b/c/d bottom row.
"""
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import base64,json,shutil,subprocess,hashlib

HERE=Path(__file__).resolve().parent
S=HERE/'source_archive';S.mkdir(exist_ok=True)
OLD=Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\figures\main_revision_v12_20261008')
for f in ['local_inset.png','inset_receipt.json']:
    if not (S/f).exists():shutil.copy2(OLD/'assets'/f,S/f)
whole=json.loads((HERE/'assets3d/consistent_a_b_receipt.json').read_text(encoding='utf8'))['b']
inset=json.loads((S/'inset_receipt.json').read_text(encoding='utf8'))
frames=json.loads((S/'local_frame_candidates.json').read_text(encoding='utf8'))['all_204_local_frames']
frames={r['residue_number']:r for r in frames}
root=ET.Element('root');ET.SubElement(root,'mxCell',id='0');ET.SubElement(root,'mxCell',id='1',parent='0')
BLUE='#54788F';PURPLE='#81709D';INK='#28363E';GREY='#A6B4BC';LIGHT='#D5DEE3'

def cell(ident,value,box,style):
    c=ET.SubElement(root,'mxCell',id='b47_'+ident,value=value,vertex='1',parent='1',style=style)
    g=ET.SubElement(c,'mxGeometry',attrib={'as':'geometry'})
    for k,v in zip(['x','y','width','height'],box):g.set(k,str(float(v)))
    return c

def label(ident,value,box,color=INK,size=28,bold=False,align='center',background='none'):
    return cell(ident,value,box,f'text;html=0;whiteSpace=wrap;overflow=visible;fillColor={background};strokeColor=none;fontFamily=Arial;fontSize={size};fontColor={color};fontStyle={int(bold)};align={align};verticalAlign=middle;spacing=0;')

def dot(ident,p,r=4,c=PURPLE):
    return cell(ident,'',[p[0]-r,p[1]-r,2*r,2*r],f'ellipse;fillColor={c};strokeColor=#FFFFFF;strokeWidth=1.8;')

def wire(ident,points,color=PURPLE,width=2,arrow=True,dashed=False,halo=False):
    c=ET.SubElement(root,'mxCell',id='b47_'+ident,parent='1',edge='1',
        style=f'edgeStyle=none;html=0;rounded=0;strokeColor={color};strokeWidth={width};endArrow={"block" if arrow else "none"};endFill=1;endSize=5;startArrow=none;dashed={int(dashed)};dashPattern=4 3;')
    g=ET.SubElement(c,'mxGeometry',relative='1',attrib={'as':'geometry'})
    for key,p in [('sourcePoint',points[0]),('targetPoint',points[-1])]:
        ET.SubElement(g,'mxPoint',x=str(float(p[0])),y=str(float(p[1])),attrib={'as':key})
    if len(points)>2:
        a=ET.SubElement(g,'Array',attrib={'as':'points'})
        for p in points[1:-1]:ET.SubElement(a,'mxPoint',x=str(float(p[0])),y=str(float(p[1])))
    return c

def image(ident,path,box):
    data=base64.b64encode(Path(path).read_bytes()).decode()
    return cell(ident,'',box,'shape=image;imageAspect=0;aspect=fixed;image=data:image/png,'+data+';')

label('title','b   Global and local frames',[15,9,535,48],size=36,bold=True,align='left')
DY=61
def by(p):return np.array(p)+[0,DY]

# Complete real secondary-structure ribbon is the primary visual object.
protein=[15,50+DY,390,390*whole['image_size'][1]/whole['image_size'][0]]
image('whole_ribbon',HERE/'assets3d/chain_d_consistent.png',protein)
label('chain','Chain D · 204 residues',[15,3+DY,310,36],color='#667C89',size=28,align='left')
camera=np.array(whole['camera_matrix']);projection=camera[:2].copy();projection[1]*=-1
def pworld(number):return np.array(protein[:2])+np.array(whole['pixels'][f'D{number}_origin'])*protein[2]/whole['image_size'][0]
positions={}
for n,name in [(67,'j'),(166,'i')]:
    o=pworld(n);positions[n]=o;R=np.array(frames[n]['R_columns_xyz'])
    # A white hairline under each short true projected axis improves separation
    # from the real ribbon without introducing RGB axes over the whole protein.
    for k,axis in enumerate('xyz'):
        end=o+35*(projection@R[:,k])
        wire(f'{name}_axis_halo_{axis}',[o,end],color='#FFFFFF',width=5.5)
        wire(f'{name}_axis_{axis}',[o,end],color=PURPLE,width=2.2)
    dot(name+'_origin',o,r=4.3)
    shift=np.array([-53,-9]) if name=='j' else np.array([-5,29])
    label(name+'_name',name,[*(o+shift),40,34],color=PURPLE,size=31,bold=True,background='#FFFFFF')

# The inset is a true 12-atom peptide backbone (D165–D167), not a flat diagram.
ib=[335,211+DY,185,185*inset['image_size'][1]/inset['image_size'][0]]
card=cell('inset_card','',[ib[0]-7,ib[1]-50,ib[2]+14,ib[3]+55],
          'rounded=1;arcSize=7;fillColor=#FFFFFF;strokeColor=#DCE3E7;strokeWidth=1.2;')
image('atomic_inset',S/'local_inset.png',ib)
def atom(key):return np.array(ib[:2])+np.array(inset['pixels'][key])*ib[2]/inset['image_size'][0]
origin=atom('origin')
# Connect the selected real Cα to its magnification across the white gap.
o=positions[166]
cell('target_ring','',[o[0]-12,o[1]-12,24,24],'ellipse;fillColor=none;strokeColor=#A594B8;strokeWidth=1.8;')
wire('zoom_leader',[o+[13,0],[323,o[1]],[323,origin[1]],origin+[-17,0]],color='#A594B8',width=1.5,arrow=False,dashed=True)
label('inset_name','Residue i',[331,166+DY,192,34],color=PURPLE,size=28)
for axis in 'xyz':
    end=atom(axis)
    wire('inset_axis_halo_'+axis,[origin,end],color='#FFFFFF',width=5.5)
    wire('inset_axis_'+axis,[origin,end],color=PURPLE,width=2.0)
    outward=(end-origin)/np.linalg.norm(end-origin)
    pos=end+outward*13
    label('inset_axis_name_'+axis,axis,[pos[0]-12,pos[1]-18,26,34],color=PURPLE,size=28)
for key,name,offset in [('origin','Cα',[-61,7]),('N','N',[-76,-30]),('C','C',[49,-23])]:
    p=atom(key)+offset
    label('atom_'+key,name,[p[0],p[1],45,34],size=28,color='#4D6070')
    if key=='N':wire('atom_leader_N',[[p[0]+44,p[1]+23],atom(key)+[-7,-4]],color='#B2BEC6',width=1.2,arrow=False)
    if key=='C':wire('atom_leader_C',[[p[0],p[1]+18],atom(key)+[9,-2]],color='#B2BEC6',width=1.2,arrow=False)

# A small global triad stays in white space and uses the same physical camera.
g=by([123,360])
for k,axis in enumerate('XYZ'):
    v=projection[:,k]
    end=g+60*v
    wire('global_'+axis,[g,end],color=GREY,width=1.7)
    v=v/np.linalg.norm(v);p=end+v*12
    label('global_name_'+axis,axis,[p[0]-13,p[1]-16,28,32],color='#8D9DA6',size=28)
dot('global_origin',g,r=3,c='#AAB7BF')
label('global_name','Global axes',[24,454+DY,210,34],color='#7B909D',size=28)
label('equation','pᵢ = Rᵢqᵢ + tᵢ',[260,450+DY,290,36],size=30)

mx=ET.Element('mxfile',host='app.diagrams.net',version='24.7.17',type='device')
d=ET.SubElement(mx,'diagram',id='b47',name='Real 3D local frames')
m=ET.SubElement(d,'mxGraphModel',page='1',pageWidth='565',pageHeight='575',grid='0',background='#FFFFFF');m.append(root)
tree=ET.ElementTree(mx);ET.indent(tree,space='  ')
source=HERE/'fig1b_prototype.drawio';tree.write(source,encoding='utf-8',xml_declaration=True)
exe='C:/Program Files/draw.io/draw.io.exe'
for ext in ['png','svg','pdf']:
    command=[exe,'-x','-f',ext,'--crop','-o',str(HERE/('fig1b_prototype.'+ext)),str(source)]
    if ext=='png':command[4:4]=['-s','2']
    proc=subprocess.run(command,capture_output=True,text=True,encoding='utf8',errors='replace',timeout=120)
    assert proc.returncode==0,proc.stdout+proc.stderr
print(source)
