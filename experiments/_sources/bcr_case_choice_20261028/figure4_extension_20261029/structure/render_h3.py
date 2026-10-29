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
import hashlib, json, math, os
import numpy as np
from PIL import Image
_dll=os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd

ROOT=Path(__file__).resolve().parent
SOURCE=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\dataset\structures\A00LZ.pdb')
W,H=2200,2000
cmd.load(str(SOURCE),'fab')
original=np.array(cmd.get_coords('fab'),float)

def xyz(sel):
    return np.asarray(cmd.get_coords('fab and '+sel),float)

r96=xyz('chain C and resi 96 and name CZ')[0]
r101=xyz('chain C and resi 101 and name CZ')[0]
ca95=xyz('chain C and resi 95 and name CA')[0]
ca109=xyz('chain C and resi 109 and name CA')[0]
mid=xyz('chain C and resi 100-103 and name CA').mean(0)
x=r101-r96;x/=np.linalg.norm(x)
y=mid-(ca95+ca109)/2;y-=x*np.dot(y,x);y/=np.linalg.norm(y)
z=np.cross(x,y);z/=np.linalg.norm(z)
y=np.cross(z,x)
center=(r96+r101+ca95+ca109+mid)/5
rotation=np.column_stack([x,y,z])

for setting,value in [('orthoscopic',1),('ray_shadows',0),('specular',.12),('antialias',2),('cartoon_fancy_helices',1),('cartoon_smooth_loops',1),('ray_opaque_background',0),('depth_cue',0),('ray_trace_fog',0),('ambient',.42),('direct',.70),('reflect',.10),('ray_trace_mode',0)]:
    cmd.set(setting,value)
cmd.bg_color('white')
cmd.set_color('heavy_blue',[35/255,106/255,182/255])
cmd.set_color('highlight',[190/255,125/255,53/255])
cmd.set_color('partner_gray',[.60,.63,.65])
cmd.set_color('antigen_gray',[.77,.79,.80])
cmd.set_color('blue_context',[.53,.68,.83])

def show_scene(include_antigen=True):
    cmd.hide('everything','all')
    cmd.color('heavy_blue','fab and chain C')
    cmd.color('partner_gray','fab and chain B')
    cmd.color('antigen_gray','fab and chain A')
    cmd.select('h3','fab and chain C and resi 95-109')
    cmd.select('nearheavy','fab and chain C and resi 25-33+94-112')
    cmd.select('nearlight','fab and chain B and resi 30-34+47-53')
    cmd.select('nearantigen','fab and chain A and resi 57-60+65-68')
    cmd.show('cartoon','nearheavy or nearlight')
    cmd.color('blue_context','nearheavy')
    cmd.color('heavy_blue','h3')
    cmd.cartoon('tube','nearheavy or nearlight')
    cmd.set('cartoon_tube_radius',.18)
    cmd.set('cartoon_transparency',.42,'nearheavy')
    cmd.set('cartoon_transparency',.25,'nearlight')
    cmd.set('cartoon_transparency',.02,'h3')
    if include_antigen:
        cmd.show('cartoon','nearantigen')
        cmd.cartoon('tube','nearantigen')
        cmd.set('cartoon_transparency',.50,'nearantigen')
    key='fab and ((chain C and resi 96+101+108) or (chain B and resi 49)) and not hydro'
    cmd.show('sticks',key)
    cmd.set('stick_radius',.17,key)
    cmd.color('highlight','fab and chain C and resi 96+101')
    cmd.color('oxygen',key+' and elem O')
    cmd.color('nitrogen',key+' and elem N')
    # Use selected N/O colors only on atom tips, preserving a small semantic palette.
    # Delete old measurements before each rigid camera transform; otherwise
    # PyMOL retains atom coordinates from earlier candidate views as extra dashes.
    cmd.delete('r101_backbone_contact')
    cmd.delete('r96_internal_contact')
    cmd.distance('r101_backbone_contact','fab and chain C and resi 101 and name O','fab and chain B and resi 49 and name OH')
    cmd.distance('r96_internal_contact','fab and chain C and resi 96 and name NH1','fab and chain C and resi 108 and name OD2')
    cmd.show('dashes','r101_backbone_contact')
    cmd.show('dashes','r96_internal_contact')
    cmd.hide('labels','all')
    cmd.set('dash_color','gray35')
    cmd.set('dash_radius',.055)
    cmd.set('dash_gap',.20)
    cmd.set('dash_length',.15)

views=[]
for angle in [0,30,-30]:
    a=math.radians(angle)
    tilt=np.array([[1,0,0],[0,math.cos(a),-math.sin(a)],[0,math.sin(a),math.cos(a)]])
    world_to_camera=rotation@tilt
    cmd.load_coords((original-center)@world_to_camera,'fab')
    show_scene(True)
    cmd.viewport(W,H)
    cmd.set_view((1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,-65.,0.,0.,0.,25.,105.,0.))
    cmd.zoom('(h3 or (fab and chain B and resi 49) or (fab and chain C and resi 108)) and not hydro',buffer=2.6,complete=1)
    cmd.set('orthoscopic',1)
    view=list(cmd.get_view())
    fname=f'H3_local_angle_{angle:+03d}.png'
    cmd.png(str(ROOT/fname),W,H,dpi=600,ray=1)
    views.append({'angle':angle,'PNG':fname,'view':view,'world_to_camera':world_to_camera.tolist()})
    print('rendered',fname,flush=True)

# Initial camera is selected without changing molecular coordinates, cropping,
# or any atomic-distance calculation. All candidate views are retained.
chosen=views[0]
world_to_camera=np.array(chosen['world_to_camera'])
cmd.load_coords((original-center)@world_to_camera,'fab')
show_scene(False);cmd.set_view(chosen['view'])
cmd.png(str(ROOT/'H3_local.png'),W,H,dpi=600,ray=1)
cmd.save(str(ROOT/'H3_local_scene.pse'))
probes=[np.array([0.,0.,0.]),np.array([3.,0.,0.]),np.array([0.,3.,0.]),np.array([-3.,-3.,3.])]
colors=['red','green','blue','yellow']
cmd.disable('all')
for i,(p,c) in enumerate(zip(probes,colors)):
    cmd.pseudoatom(f'probe{i}',pos=p.tolist())
    cmd.show('spheres',f'probe{i}')
    cmd.set('sphere_scale',.15,f'probe{i}')
    cmd.color(c,f'probe{i}')
cmd.set('ambient',1);cmd.set('direct',0);cmd.set('specular',0)
cmd.set_view(chosen['view'])
cmd.png(str(ROOT/'projection_calibration.png'),W,H,dpi=600,ray=1)
arr=np.array(Image.open(ROOT/'projection_calibration.png').convert('RGBA'))
r,g,b,a=[arr[:,:,k].astype(int) for k in range(4)]
masks=[(r>180)&(g<70)&(b<70)&(a>150),(g>150)&(r<70)&(b<70)&(a>150),(b>180)&(r<70)&(g<70)&(a>150),(r>180)&(g>180)&(b<70)&(a>150)]
screen=[]
for mask in masks:
    yy,xx=np.nonzero(mask)
    if not len(xx):raise RuntimeError('Projection probe not visible')
    screen.append([float(xx.mean()),float(yy.mean())])
mapping=np.linalg.solve(np.column_stack([probes,np.ones(4)]),np.array(screen))
cmd.disable('all');holdout=np.array([1.,-1.2,.7])
cmd.pseudoatom('holdout',pos=holdout.tolist());cmd.show('spheres','holdout');cmd.set('sphere_scale',.15,'holdout');cmd.color('magenta','holdout');cmd.set_view(chosen['view'])
cmd.png(str(ROOT/'projection_holdout.png'),W,H,dpi=600,ray=1)
arr=np.array(Image.open(ROOT/'projection_holdout.png').convert('RGBA'))
r,g,b,a=[arr[:,:,k].astype(int) for k in range(4)]
yy,xx=np.nonzero((r>180)&(b>180)&(g<70)&(a>150))
observed=np.array([xx.mean(),yy.mean()]);predicted=np.r_[holdout,1.]@mapping
error=float(np.linalg.norm(observed-predicted))
assert error<2.0,error
cmd.enable('fab')

anchors={
 'R96_CA':'chain C and resi 96 and name CA',
 'R96_CZ':'chain C and resi 96 and name CZ',
 'R96_NH1':'chain C and resi 96 and name NH1',
 'R101_CA':'chain C and resi 101 and name CA',
 'R101_CZ':'chain C and resi 101 and name CZ',
 'R101_O':'chain C and resi 101 and name O',
 'D108_OD2':'chain C and resi 108 and name OD2',
 'Light_Y49_OH':'chain B and resi 49 and name OH',
 'H3_tip_CA':'chain C and resi 102 and name CA',
 'Heavy_near_R96':'chain C and resi 30 and name CA',
 'Light_context':'chain B and resi 50 and name CA',
 'Antigen_context':'chain A and resi 67 and name CA',
}
projected={k:(np.r_[xyz(s)[0],1.]@mapping).tolist() for k,s in anchors.items()}
metadata={'status':'rendered_pending_visual_selection','source':str(SOURCE),
 'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
 'all_atoms_rendered_from_source_coordinates_with_rigid_camera_transform_only':True,
 'size':[W,H],'PNG':str(ROOT/'H3_local.png'),'coordinates_origin':'upper-left PNG pixel','candidate_views':views,
 'selected_angle':chosen['angle'],'screen_projection':{'xyz1_camera_to_screen':mapping.tolist(),'holdout_error_pixels':error},
 'anchors_screen_pixels':projected,'center_original':center.tolist(),
 'distances_A':{
  'R101_O_Light_Y49_OH':float(cmd.get_distance('fab and chain C and resi 101 and name O','fab and chain B and resi 49 and name OH')),
  'R96_NH1_D108_OD2':float(cmd.get_distance('fab and chain C and resi 96 and name NH1','fab and chain C and resi 108 and name OD2'))},
 'scope':'One static WT Fab-antigen coordinate view, focused on H-CDR3, selected heavy-chain context and light-chain Y49. Antigen background omitted from final local view to avoid projected overlap with the solvent-facing R101 sidechain. Dashes show selected short distances, not independently established hydrogen bonds or affinity contributions.',
 'colors':{'heavy':'#236AB6','R96_R101':'#BE7D35','light':'#999FA6','antigen':'#C4C9CC','oxygen':'PyMOL oxygen','nitrogen':'PyMOL nitrogen'},
 'render_settings':{'orthographic':True,'ray_shadows':False,'depth_cue':False,'transparent_background':True,'ray_dpi':600}}
(ROOT/'H3_molecular_view.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
print(json.dumps({'projection_error':error,'anchors':projected,'distances':metadata['distances_A']},indent=2),flush=True)
cmd.quit()
