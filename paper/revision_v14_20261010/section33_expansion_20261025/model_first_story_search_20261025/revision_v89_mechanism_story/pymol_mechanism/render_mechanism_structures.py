# Academic Figure Skill Asset Confirmation
# Molecular structure -> v88/pymol/render_1MXE_v85.py -> copy-first PyMOL style/camera.
# Exact observed PDB atoms only. Figure annotations are composed separately as vector text.
from pathlib import Path
import os,json,hashlib,shutil
import numpy as np
from PIL import Image
ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
S=Path(r'E:/Project/FlexPPI/paper/revision_v14_20261010/section33_expansion_20261025/model_first_story_search_20261025')
OUT=Path(__file__).resolve().parent;INP=OUT/'inputs';INP.mkdir(exist_ok=True)
OLD=S/'revision_v88_two_panel_story/source_analysis/M316_static_context/inputs'
for fn in ['1A06_rigidly_aligned_to_1MXE.pdb','1MXE_CaM_native4_147_CaMKI294_318.pdb','1MXE_camera_manifest.json']:shutil.copyfile(OLD/fn,INP/fn)
VIEW=json.loads((INP/'1MXE_camera_manifest.json').read_text())['assets']['overview']['view']
COLORS={'kinase':'#BAC1C7','N_lobe':'#8A9299','C_lobe':'#BAC1C7','linker':'#D9DDE0','peptide':'#8B72A5','key':'#3A8D8C','contact':'#767F87'}
def rgb(h):return [int(h[i:i+2],16)/255 for i in (1,3,5)]
def project(c,sel,w,h):
 xyz=np.asarray(c.get_coords(sel)).mean(0);v=np.array(c.get_view());q=v[:9].reshape(3,3).T@(xyz-v[12:15])+v[9:12]
 hh=-v[11]*np.tan(np.deg2rad(float(c.get('field_of_view'))/2));hw=hh*w/h
 return [float((.5+q[0]/(2*hw))*w),float((.5-q[1]/(2*hh))*h)]
def style(c,code):
 obj='kinase' if code=='1A06' else 'bound';ch='A' if code=='1A06' else 'B'
 c.hide('everything',obj);c.show('cartoon',obj+' and polymer.protein')
 if code=='1A06':
  c.color('kinase',obj);c.set('cartoon_color','kinase',obj)
  c.color('peptide',f'{obj} and resi 294-316');c.set('cartoon_color','peptide',f'{obj} and resi 294-316')
  # Show observed F298 neighbours as thin neutral sticks while retaining full kinase context.
  sel=f'{obj} and resi 104+184+210+214+216 and not name N+C+O'
  c.show('sticks',sel);c.color('contact',sel);c.set('stick_radius',.12,sel)
 else:
  for col,sel in [('N_lobe','chain A and resi 4-75'),('linker','chain A and resi 76-81'),('C_lobe','chain A and resi 82-147'),('peptide','chain B')]:
   c.color(col,f'{obj} and ({sel})');c.set('cartoon_color',col,f'{obj} and ({sel})')
 for site in [296,298,303,316]:
  sel=f'{obj} and chain {ch} and resi {site} and not name N+C+O'
  c.show('sticks',sel);c.color('key',sel);c.set('stick_radius',.27,sel)
 # Surface purposely absent: full observed cartoons and pockets remain spatially readable.
 return obj,ch

def render(c,name,obj,ch,w=1800,h=1500):
 c.viewport(w,h);c.png(str(OUT/(name+'.png')),width=w,height=h,dpi=800,ray=1,quiet=1)
 im=Image.open(OUT/(name+'.png')).convert('RGBA');bb=im.getchannel('A').getbbox();pad=45
 crop=[max(0,bb[0]-pad),max(0,bb[1]-pad),min(w,bb[2]+pad),min(h,bb[3]+pad)]
 im.crop(crop).save(OUT/(name+'_trim.png'))
 marks={}
 for site,aa,atom in [(296,'K','NZ'),(298,'F','CZ'),(303,'W','CZ2'),(316,'M','CE')]:marks[f'{aa}{site}']=f'{obj} and chain {ch} and resi {site} and name {atom}'
 if obj=='kinase':
  marks.update(kinase_center=f'{obj} and resi 10-279 and name CA',F298_pocket=f'{obj} and resi 104+184+210+214+216 and name CA',G214=f'{obj} and resi 214 and name O',F104=f'{obj} and resi 104 and name CE2',I210=f'{obj} and resi 210 and name CG2',regulatory_center=f'{obj} and resi 294-316 and name CA',E292=f'{obj} and resi 292 and name OE2')
 else:
  marks.update(CaM_N=f'{obj} and chain A and resi 4-75 and name CA',CaM_C=f'{obj} and chain A and resi 82-147 and name CA',CaM_E127=f'{obj} and chain A and resi 127 and name OE2',regulatory_center=f'{obj} and chain B and name CA')
 lm={k:project(c,s,w,h) for k,s in marks.items()}
 entry={'name':name,'file':name+'.png','trim_file':name+'_trim.png','width':w,'height':h,'view':list(c.get_view()),'field_of_view':float(c.get('field_of_view')),'alpha_bbox':list(bb),'crop_box':crop,'trim_width':crop[2]-crop[0],'trim_height':crop[3]-crop[1],'landmarks':lm,'trim_landmarks':{k:[v[0]-crop[0],v[1]-crop[1]] for k,v in lm.items()}}
 c.scene(name,'store');c.save(str(OUT/(name+'.pse')))
 print('RENDERED',name,flush=True);return entry

def main():
 manifest={'palette':COLORS,'source_files':{},'assets':{},'scope':'Exact observed structures, not conformational interpolation. 1A06 had already been rigidly aligned to 1MXE by prior work; this renderer preserves all loaded atomic coordinates. Only camera zoom/crop and representation differ. 1MXE is isolated regulatory peptide with CaM, not activated full kinase. No model/MD prediction for 1A06.'}
 with pymol2.PyMOL() as p:
  c=p.cmd
  for k,v in dict(ray_opaque_background=0,opaque_background=0,ray_trace_mode=0,ray_shadows=0,antialias=3,orthoscopic=1,depth_cue=0,specular=.12,shininess=35,ambient=.65,direct=.35,reflect=.12,light_count=2,cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,cartoon_smooth_loops=1,cartoon_sampling=20,cartoon_loop_radius=.24,cartoon_tube_radius=.24,cartoon_oval_length=1.25,cartoon_oval_width=.28,max_threads=4,valence=0).items():c.set(k,v)
  for k,v in COLORS.items():c.set_color(k,rgb(v))
  c.bg_color('white')
  coords={};zview=None
  for code,fn in [('1A06','1A06_rigidly_aligned_to_1MXE.pdb'),('1MXE','1MXE_CaM_native4_147_CaMKI294_318.pdb')]:
   obj='kinase' if code=='1A06' else 'bound';c.load(str(INP/fn),obj);c.dss(obj);coords[obj]=c.get_coords(obj).copy();manifest['source_files'][code]={'file':'inputs/'+fn,'sha256':hashlib.sha256((INP/fn).read_bytes()).hexdigest(),'loaded_atom_count':len(coords[obj])}
   c.disable('all');c.enable(obj);obj,ch=style(c,code)
   c.viewport(1800,1500);c.set_view(VIEW);c.zoom(obj,buffer=2.0,complete=1);v=list(c.get_view());v[:9]=VIEW[:9];v[15]=1.;v[16]=-v[11]+200.;c.set_view(v)
   if code=='1A06':zview=list(v)
   manifest['assets'][code+'_overview']=render(c,code+'_overview',obj,ch)
   if code=='1MXE':
    vv=list(v);vv[9:12]=zview[9:12];vv[15]=1.;vv[16]=-vv[11]+200.;c.set_view(vv)
    manifest['assets'][code+'_matched_scale']=render(c,code+'_matched_scale',obj,ch)
   # Separate unlabelled regulatory segment uses the same observed atom conformation and rotation.
   only=obj+'_segment';c.create(only,f'{obj} and chain {ch} and resi 294-318');c.disable('all');c.enable(only)
   c.hide('everything',only);c.show('cartoon',only);c.color('peptide',only);c.set('cartoon_color','peptide',only)
   for site in [296,298,303,316]:
    ss=f'{only} and resi {site} and not name N+C+O';c.show('sticks',ss);c.color('key',ss);c.set('stick_radius',.27,ss)
   c.set_view(VIEW);c.zoom(only,buffer=2,complete=1);v=list(c.get_view());v[:9]=VIEW[:9];v[15]=1.;v[16]=-v[11]+200.;c.set_view(v)
   # Rendering helper expects original chain and derives context selections; temporary object alias avoids missing marks.
   c.viewport(1800,1000);c.png(str(OUT/(code+'_regulatory_segment.png')),width=1800,height=1000,dpi=800,ray=1,quiet=1)
   im=Image.open(OUT/(code+'_regulatory_segment.png')).convert('RGBA');bb=im.getchannel('A').getbbox();im.crop((max(0,bb[0]-30),max(0,bb[1]-30),min(1800,bb[2]+30),min(1000,bb[3]+30))).save(OUT/(code+'_regulatory_segment_trim.png'));c.save(str(OUT/(code+'_regulatory_segment.pse')));c.delete(only)
  manifest['geometry_audit']={obj:{'count':len(before),'exactly_unchanged':bool(np.array_equal(before,c.get_coords(obj))),'max_abs_coordinate_change_A':float(abs(before-c.get_coords(obj)).max())} for obj,before in coords.items()}
  assert all(a['exactly_unchanged'] for a in manifest['geometry_audit'].values())
  manifest['same_rotation_all_overviews']=all(np.allclose(a['view'][:9],VIEW[:9],atol=1e-7) for a in manifest['assets'].values())
  c.disable('all');c.enable('kinase');c.set_view(zview);c.save(str(OUT/'both_observed_structures.pse'))
 (OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 print('COMPLETE',flush=True)
if __name__=='__main__':main()