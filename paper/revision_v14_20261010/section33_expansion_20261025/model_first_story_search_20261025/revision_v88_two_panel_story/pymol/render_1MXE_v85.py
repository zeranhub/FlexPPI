# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) Molecular structure → revision_v83_model_story/pymol/render_palette_v83.py → param inherit
# RULE: Scientific atomic coordinates are rendered in PyMOL, with no generated atoms.
from pathlib import Path
import os,json,hashlib
import numpy as np
ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
if hasattr(os,'add_dll_directory'):
 DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
from PIL import Image
OUT=Path(__file__).resolve().parent
INPUT=OUT/'inputs/1MXE_CaM_native4_147_CaMKI294_318.pdb'
COLORS={'N_lobe':'#8A9299','C_lobe':'#BAC1C7','linker':'#D9DDE0','partner':'#8B72A5','key':'#3A8D8C','sulfur':'#B18F31','oxygen':'#A83D3D','nitrogen':'#447BA7','context':'#A6AEB4'}
def rgb(h):return [int(h[i:i+2],16)/255 for i in (1,3,5)]
def project(c,sel,w,h):
 xyz=c.get_coords(sel).mean(0);v=np.array(c.get_view());q=v[:9].reshape(3,3).T@(xyz-v[12:15])+v[9:12]
 hh=-v[11]*np.tan(np.deg2rad(float(c.get('field_of_view'))/2));hw=hh*w/h
 return [float((.5+q[0]/(2*hw))*w),float((.5-q[1]/(2*hh))*h)]
def style(c,obj,site=None):
 c.hide('everything',obj);c.show('cartoon',obj+' and polymer.protein')
 for color,sel in [('N_lobe','chain A and resi 4-75'),('linker','chain A and resi 76-81'),('C_lobe','chain A and resi 82-147'),('partner','chain B')]:
  c.color(color,obj+' and ('+sel+')');c.set('cartoon_color',color,obj+' and ('+sel+')')
 sites=[297,303,308,316] if site is None else [site]
 for s in sites:
  key=f'{obj} and chain B and resi {s}';c.show('sticks',key+' and not name N+C+O');c.set('stick_radius',.24 if site is None else .23,key)
  c.color('key',key+' and elem C');c.color('nitrogen',key+' and elem N');c.color('sulfur',key+' and elem S');c.color('oxygen',key+' and elem O')
 if site is not None:
  focus=f'{obj} and chain B and resi {site}'
  fxyz=np.asarray(c.get_coords(f'({focus}) and not name N+CA+C+O'))
  contact=sorted(set((a.chain,a.resi,a.resn) for a in c.get_model(f'{obj} and chain A').atom if np.min(np.linalg.norm(fxyz-np.asarray(a.coord),axis=1))<=4.5),key=lambda t:int(t[1]))
  neighbor=f'({obj} and chain A and resi '+ '+'.join(t[1] for t in contact)+')'
  c.show('sticks',neighbor+' and not name N+C+O');c.set('stick_radius',.12,neighbor);c.set('stick_transparency',.12,neighbor)
  for color,element in [('oxygen','O'),('nitrogen','N'),('sulfur','S')]:c.color(color,neighbor+' and elem '+element)
  c.set('cartoon_transparency',.50,obj+' and chain A');c.set('cartoon_transparency',.06,obj+' and chain B')
  return contact
def render(c,name,w,h,marks):
 c.viewport(w,h);c.png(str(OUT/(name+'.png')),width=w,height=h,dpi=600,ray=1,quiet=1)
 im=Image.open(OUT/(name+'.png')).convert('RGBA');bbox=list(im.getchannel('A').getbbox())
 entry={'file':name+'.png','width':w,'height':h,'dpi':600,'view':list(c.get_view()),'field_of_view':float(c.get('field_of_view')),'alpha_bbox':bbox,'landmarks':{k:project(c,s,w,h) for k,s in marks.items()}}
 
 if name=='overview':
  bb=[max(0,bbox[0]-40),max(0,bbox[1]-40),min(w,bbox[2]+40),min(h,bbox[3]+40)]
  im.crop(bb).save(OUT/(name+'_trim.png'))
  entry['trim']={'file':name+'_trim.png','crop_box':bb,'width':bb[2]-bb[0],'height':bb[3]-bb[1],'landmarks':{k:[p[0]-bb[0],p[1]-bb[1]] for k,p in entry['landmarks'].items()}}
 c.scene(name,'store');c.save(str(OUT/(name+'.pse')));print('RENDERED',name,flush=True)
 return entry

def main():
 manifest={'palette':COLORS,'assets':{},'scope':'Exact family-fold1 input heavy atoms. Native numbering only; no atom added, rebuilt or moved. All local views preserve overview camera rotation; zoom/crop changes only. Full-chain cartoons remain present in local PSEs.','input_sha256':hashlib.sha256(INPUT.read_bytes()).hexdigest()}
 with pymol2.PyMOL() as pm:
  c=pm.cmd
  for k,v in dict(ray_opaque_background=0,opaque_background=0,ray_trace_mode=0,ray_shadows=0,antialias=3,orthoscopic=1,depth_cue=0,specular=.12,shininess=35,ambient=.65,direct=.35,reflect=.12,light_count=2,cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,cartoon_smooth_loops=1,cartoon_sampling=20,cartoon_loop_radius=.24,cartoon_tube_radius=.24,cartoon_oval_length=1.25,cartoon_oval_width=.28,max_threads=4,stick_radius=.17,valence=0).items():c.set(k,v)
  for k,v in COLORS.items():c.set_color(k,rgb(v))
  c.bg_color('white');c.load(str(INPUT),'complex');c.dss('complex');style(c,'complex')
  before=np.array(c.get_coords('complex')).copy()
  c.viewport(2100,1500);c.orient('complex');c.turn('z',-20);c.turn('x',10);c.zoom('complex',buffer=1.5,complete=1)
  view=c.get_view()
  marks={f'{aa}{n}':f'complex and chain B and resi {n} and name '+('CZ2' if n==303 else 'CE' if n==316 else 'CG') for aa,n in [('W',303),('M',316),('N',297),('N',308)]}
  marks.update(CaM_N='complex and chain A and resi 4-75 and name CA',CaM_C='complex and chain A and resi 82-147 and name CA',peptide='complex and chain B and name CA')
  manifest['assets']['overview']=render(c,'overview',2100,1500,marks)
  rotation=np.array(view[:9]);local_contact={}
  for site in [303,316]:
   c.disable('all');obj='local_'+str(site);c.create(obj,'complex');c.enable(obj);local_contact[str(site)]=style(c,obj,site)
   c.set_view(view);c.viewport(1200,1500);c.zoom(f'{obj} and chain B and resi {site}',buffer=7,complete=1)
   # Fix an identical optical scale for both 1–14 anchor pockets.
   lv=list(c.get_view());lv[:9]=view[:9]
   if site==303:local_scale=lv[9:12];local_clip=[lv[15],lv[16]]
   else:lv[9:12]=local_scale
   # Keep the depth slab broad: local zoom must not cut through helices.
   lv[15]=1.;lv[16]=-lv[11]+150.
   c.set_view(lv)
   lm={f'{aa}{n}':f'{obj} and chain B and resi {n} and name '+('CZ2' if n==303 else 'CE' if n==316 else 'CG') for aa,n in [('W',303),('M',316),('N',297),('N',308)]}
   for chain,resn,residue in local_contact[str(site)]:lm['CaM_'+residue+resn]=f'{obj} and chain {chain} and resi {resn} and name CA'
   manifest['assets']['local_'+str(site)]=render(c,'local_'+str(site),1200,1500,lm)
   c.delete(obj)
  after=np.array(c.get_coords('complex'));assert np.array_equal(before,after)
  manifest['coordinate_render_audit']={'atom_count':len(before),'exactly_unchanged':True,'max_coordinate_change_A':float(abs(before-after).max())}
  manifest['local_contact_residues_within_4p5A_of_anchor_sidechain']=local_contact
  manifest['pymol_version']=str(c.get_version()[0])
  manifest['camera_rotation_consistency']=bool(all(np.allclose(entry['view'][:9],view[:9],atol=1e-7) for entry in manifest['assets'].values()))
  manifest['local_scale_consistency']=bool(np.allclose(manifest['assets']['local_303']['view'][9:12],manifest['assets']['local_316']['view'][9:12]))
 (OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 print('COMPLETE',str(OUT),flush=True)
if __name__=='__main__':main()