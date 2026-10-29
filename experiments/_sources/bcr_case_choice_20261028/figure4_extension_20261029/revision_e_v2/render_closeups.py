# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (e) Molecular structure -> cross-type inherit from existing real-coordinate
#     mel5_interface_v51_20261022/structure/render_pair.py -> param inherit
# No molecular structure production asset exists in the skill directory map.
# Quantitative data are not encoded in this molecular view. Existing ray,
# orthographic camera, restrained palette and projection QA are inherited.
# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})
# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"
# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

from pathlib import Path
import os,sys,json,hashlib,math
import numpy as np
from PIL import Image,ImageDraw
_dll=os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd
ROOT=Path(__file__).resolve().parent
SOURCE=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\dataset\structures\A00LZ.pdb')
cmd.load(str(SOURCE),'fab')
original=np.asarray(cmd.get_coords('fab'),float)
W,H=1800,900
def xyz(sel):return np.asarray(cmd.get_coords('fab and '+sel),float)[0]
cmd.set_color('blue',[.137,.416,.714])
cmd.set_color('orange',[.745,.490,.208])
cmd.set_color('lightgray',[.62,.65,.68])
cmd.set_color('contextblue',[.74,.82,.90])
cmd.set_color('softred',[.82,.24,.23])
cmd.set_color('softnitrogen',[.21,.32,.70])
for k,v in [('orthoscopic',1),('ray_shadows',0),('specular',.15),('antialias',2),('ray_opaque_background',0),('depth_cue',0),('ray_trace_fog',0),('ambient',.53),('direct',.65),('reflect',.08),('ray_trace_mode',0),('cartoon_loop_radius',.22),('cartoon_tube_radius',.24),('cartoon_smooth_loops',1),('cartoon_flat_sheets',1),('cartoon_sampling',16),('stick_quality',24)]:cmd.set(k,v)
cmd.bg_color('white')
specs={
 'R96':{'primary':96,'partner':'chain C and resi 108','a':'chain C and resi 96 and name NH1','b':'chain C and resi 108 and name OD2',
         'context':'chain C and resi 25-33+94-110','loop':'chain C and resi 94-110','stick_extra':'chain C and resi 30+107+109'},
 'R101':{'primary':101,'partner':'chain B and resi 49','a':'chain C and resi 101 and name O','b':'chain B and resi 49 and name OH',
         'context':'(chain C and resi 97-106) or (chain B and resi 46-53)','loop':'chain C and resi 97-106','stick_extra':'none'}
}
records=[]
for name,s in specs.items():
    cmd.load_coords(original,'fab')
    ca=xyz(f'chain C and resi {s["primary"]} and name CA');cz=xyz(f'chain C and resi {s["primary"]} and name CZ')
    a=xyz(s['a']);b=xyz(s['b'])
    x=cz-ca;x/=np.linalg.norm(x)
    y=b-a;y-=x*np.dot(y,x);y/=np.linalg.norm(y)
    z=np.cross(x,y);z/=np.linalg.norm(z);y=np.cross(z,x)
    basis=np.column_stack([x,y,z]);center=(ca+cz+a+b)/4
    for angle in [-30,0,30]:
        theta=math.radians(angle)
        tilt=np.array([[1,0,0],[0,math.cos(theta),-math.sin(theta)],[0,math.sin(theta),math.cos(theta)]])
        transform=basis@tilt
        cmd.load_coords((original-center)@transform,'fab')
        cmd.hide('everything','all');cmd.delete('contact')
        cmd.color('contextblue','fab and chain C');cmd.color('lightgray','fab and chain B')
        cmd.show('cartoon','fab and ('+s['context']+')')
        cmd.cartoon('automatic','fab');cmd.set('cartoon_transparency',.27,'fab')
        cmd.color('blue','fab and ('+s['loop']+')')
        cmd.set('cartoon_transparency',.08,'fab and ('+s['loop']+')')
        primary=f'fab and chain C and resi {s["primary"]} and not hydro'
        partner='fab and ('+s['partner']+') and not hydro'
        cmd.show('sticks',primary+' or ('+partner+')')
        cmd.set('stick_radius',.18,primary+' or ('+partner+')')
        cmd.color('blue',primary);cmd.color('orange',primary+' and not name N+CA+C+O+OXT')
        cmd.color('lightgray',partner)
        if name=='R96':
            extra='fab and ('+s['stick_extra']+') and not hydro'
            cmd.show('sticks',extra);cmd.set('stick_radius',.10,extra);cmd.set('stick_transparency',.62,extra)
            cmd.color('contextblue',extra)
        key=primary+' or ('+partner+')'
        cmd.color('softred','('+key+') and elem O');cmd.color('softnitrogen','('+key+') and elem N')
        cmd.distance('contact','fab and '+s['a'],'fab and '+s['b'])
        cmd.hide('labels','all');cmd.set('dash_color','gray35');cmd.set('dash_radius',.045);cmd.set('dash_gap',.16);cmd.set('dash_length',.13)
        cmd.viewport(W,H);cmd.set_view((1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,-70.,0.,0.,0.,30.,110.,1.))
        cmd.zoom(primary+' or ('+partner+')',buffer=1.6,complete=1)
        cmd.set('orthoscopic',1)
        view=list(cmd.get_view())
        filename=f'{name}_{angle:+03d}.png'
        cmd.png(str(ROOT/filename),W,H,ray=1,dpi=600)
        anchors={k:xyz(sel).tolist() for k,sel in {'CA':f'chain C and resi {s["primary"]} and name CA','CZ':f'chain C and resi {s["primary"]} and name CZ','contact_a':s['a'],'contact_b':s['b']}.items()}
        records.append(dict(name=name,angle=angle,png=filename,view=view,anchors_camera=anchors,basis=transform.tolist(),center=center.tolist(),distance=float(cmd.get_distance('fab and '+s['a'],'fab and '+s['b']))))
        print('rendered',filename,flush=True)
(ROOT/'candidate_views.json').write_text(json.dumps({'source':str(SOURCE),'sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'size':[W,H],'views':records},indent=2),encoding='utf-8')
sheet=Image.new('RGB',(1500,560),'white');draw=ImageDraw.Draw(sheet)
for i,rec in enumerate(records):
    im=Image.open(ROOT/rec['png']).convert('RGBA');bg=Image.new('RGBA',im.size,'white');bg.alpha_composite(im);im=bg.convert('RGB');im.thumbnail((500,250))
    col=i%3;row=i//3;sheet.paste(im,(col*500,row*280+25));draw.text((col*500+12,row*280+5),rec['png'],fill='black')
sheet.save(ROOT/'candidate_contact_sheet.png')
cmd.quit()

