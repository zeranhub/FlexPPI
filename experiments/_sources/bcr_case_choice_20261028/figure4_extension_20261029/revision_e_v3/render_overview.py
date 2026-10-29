"""Actual 6A77 coordinates, using Figure 4a's exact palette and lighting setup."""
from pathlib import Path
import json, os, math, hashlib
import numpy as np
from PIL import Image, ImageDraw
_dll = os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import pymol
pymol.finish_launching(['pymol', '-cq'])
from pymol import cmd
ROOT = Path(__file__).resolve().parent
SOURCE = Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\dataset\structures\A00LZ.pdb')
cmd.load(str(SOURCE), 'fab')
original = np.array(cmd.get_coords('fab'), dtype=float)
xyz = lambda s: np.array(cmd.get_coords('fab and '+s), dtype=float)
heavy = xyz('chain C and resi 1-110 and name CA').mean(0)
light = xyz('chain B and resi 1-110 and name CA').mean(0)
base = xyz('chain B+C and resi 140-190 and name CA').mean(0)
tip = xyz('chain C and resi 95-109 and name CA').mean(0)
y = tip-base; y /= np.linalg.norm(y)
x = heavy-light; x -= y*np.dot(x,y); x /= np.linalg.norm(x)
z = np.cross(x,y); z /= np.linalg.norm(z)
x = np.cross(y,z)
basis = np.column_stack([x,y,z])
center = xyz('name CA').mean(0)
for k,v in [('orthoscopic',1),('ray_shadows',0),('specular',.15),('antialias',2),('cartoon_fancy_helices',1),('cartoon_smooth_loops',1),('ray_opaque_background',0)]:
    cmd.set(k,v)
for name,color in [('heavy','#2166AC'),('light','#92999F'),('antigen','#BAC1C7'),('focus','#B87936')]:
    cmd.set_color(name,[int(color[i:i+2],16)/255 for i in (1,3,5)])
cmd.bg_color('white')
for k,v in json.loads((ROOT/'panel_a_render_style.json').read_text())['settings'].items():cmd.set(k,v)
cmd.hide('everything','all');cmd.show('cartoon','fab')
cmd.color('heavy','fab and chain C');cmd.color('light','fab and chain B');cmd.color('antigen','fab and chain A')
cmd.color('focus','fab and chain C and resi 95-109')
W,H = 1000,1600
records=[]
for angle in [-25,0,25]:
    a=math.radians(angle)
    yaw=np.array([[math.cos(a),0,math.sin(a)],[0,1,0],[-math.sin(a),0,math.cos(a)]])
    transform=basis@yaw
    cmd.load_coords((original-center)@transform,'fab')
    cmd.viewport(W,H);cmd.set_view((1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,-140.,0.,0.,0.,60.,220.,0.))
    cmd.zoom('fab and name CA',buffer=2,complete=1)
    cmd.set('orthoscopic',1)
    filename=f'overview_{angle:+03d}.png'
    cmd.png(str(ROOT/filename),W,H,dpi=600,ray=1)
    cmd.save(str(ROOT/f'overview_{angle:+03d}.pse'))
    records.append({'angle':angle,'png':filename,'view':list(cmd.get_view()),'basis':transform.tolist(),'center':center.tolist()})
(ROOT/'overview_views.json').write_text(json.dumps({'source':str(SOURCE),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'size':[W,H],'views':records},indent=2),encoding='utf-8')
sheet=Image.new('RGB',(900,520),'white');draw=ImageDraw.Draw(sheet)
for i,rec in enumerate(records):
    im=Image.open(ROOT/rec['png']).convert('RGBA');bg=Image.new('RGBA',im.size,'white');bg.alpha_composite(im);im=bg.convert('RGB');im.thumbnail((300,480))
    sheet.paste(im,(i*300,30));draw.text((i*300+10,8),rec['png'],fill='black')
sheet.save(ROOT/'overview_contact_sheet.png')
cmd.quit()
