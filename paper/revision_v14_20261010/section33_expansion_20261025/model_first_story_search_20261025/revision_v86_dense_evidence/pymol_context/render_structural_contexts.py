# Academic Figure Skill Asset Confirmation
# (f) Molecular structural context → v85 PyMOL component palette and transparent native rendering → param inherit
# Scientific coordinates only; no atoms or unresolved residues are generated.
from pathlib import Path
import os,json,hashlib,csv
import numpy as np
ENV=Path('D:/Anaconda/envs/pvb_cpu');os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
if hasattr(os,'add_dll_directory'):DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
from PIL import Image,ImageDraw
OUT=Path(__file__).resolve().parent
COLORS={'background':'#BAC1C7','regulatory':'#8B72A5','key':'#3A8D8C','oxygen':'#A83D3D','nitrogen':'#447BA7','sulfur':'#B18F31'}
def rgb(h):return[int(h[k:k+2],16)/255 for k in [1,3,5]]
def project(c,sel,w,h):
 xyz=c.get_coords(sel).mean(0);v=np.array(c.get_view());q=v[:9].reshape(3,3).T@(xyz-v[12:15])+v[9:12];hh=-v[11]*np.tan(np.deg2rad(float(c.get('field_of_view'))/2));hw=hh*w/h
 return[float((.5+q[0]/(2*hw))*w),float((.5-q[1]/(2*hh))*h)]
def style(c,obj,chain,reg):
 c.hide('everything',obj);c.show('cartoon',obj);c.color('background',obj);c.set('cartoon_color','background',obj)
 c.color('regulatory',reg);c.set('cartoon_color','regulatory',reg)
 for site in [303,316]:
  focus=f'{obj} and chain {chain} and resi {site}';c.show('sticks',focus+' and not name N+C+O');c.set('stick_radius',.30,focus)
  c.color('key',focus+' and elem C');c.color('oxygen',focus+' and elem O');c.color('nitrogen',focus+' and elem N');c.color('sulfur',focus+' and elem S')

def contact_metrics(c,obj,chain,background):
 rows=[]
 for site in [303,316]:
  focus=f'{obj} and chain {chain} and resi {site} and not name N+CA+C+O';fa=c.get_model(focus).atom;fxyz=np.array([a.coord for a in fa]);pa=c.get_model(background).atom
  dist=np.linalg.norm(np.array([a.coord for a in pa])[:,None,:]-fxyz[None,:,:],axis=2)
  contacts=sorted(set((a.chain,int(a.resi),a.resn) for a,d in zip(pa,dist.min(axis=1)) if d<=4.5),key=lambda x:x[1]);i,j=np.unravel_index(np.argmin(dist),dist.shape)
  rows.append({'state':obj,'site':site,'min_anchor_sidechain_to_background_heavy_distance_A':float(dist[i,j]),'closest_anchor_atom':fa[j].name,'closest_background_atom':f'{pa[i].chain}:{pa[i].resn}{pa[i].resi}:{pa[i].name}','background_residues_within4p5A':contacts,'contact_residue_count_4p5A':len(contacts)})
 return rows

def sasa(c,obj,chain,reg):
 original_b=[a.b for a in c.get_model(obj).atom]
 c.set('dot_solvent',1);c.set('solvent_radius',1.4);c.set('dot_density',4)
 c.get_area(obj,load_b=1);full={site:sum(a.b for a in c.get_model(f'{obj} and chain {chain} and resi {site} and not name N+CA+C+O').atom) for site in [303,316]}
 tmp=obj+'_isolated_regulatory';c.create(tmp,reg);c.get_area(tmp,load_b=1);isol={site:sum(a.b for a in c.get_model(f'{tmp} and chain {chain} and resi {site} and not name N+CA+C+O').atom) for site in [303,316]};c.delete(tmp)
 c.alter(obj,'b=next(original_b)',space={'original_b':iter(original_b)})
 return {site:{'sidechain_SASA_whole_context_A2':full[site],'sidechain_SASA_isolated_regulatory_fragment_A2':isol[site],'sidechain_deltaSASA_buried_by_background_A2':isol[site]-full[site],'sidechain_fraction_remaining_exposed':full[site]/isol[site] if isol[site] else None} for site in [303,316]}

def main():
 manifest={'source_urls':{'1A06':'https://files.rcsb.org/download/1A06.pdb'},'limitations':['1A06 is a literature crystal structure reference. No new model or MD was run on this structure.','1MXE is a CaM-bound regulatory peptide, not an intact activated kinase.','1A06 has no experimental coordinates for residues317–318; they are not generated.','Static SASA/contact distance are structural descriptors, not free-energy or temporal sequence evidence.'],'palette':COLORS,'assets':{}}
 with pymol2.PyMOL() as pm:
  c=pm.cmd
  for k,v in dict(ray_opaque_background=0,opaque_background=0,ray_trace_mode=0,ray_shadows=0,antialias=3,orthoscopic=1,depth_cue=0,specular=.12,shininess=35,ambient=.65,direct=.35,reflect=.12,light_count=2,cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,cartoon_smooth_loops=1,cartoon_sampling=20,cartoon_loop_radius=.24,cartoon_tube_radius=.24,cartoon_oval_length=1.25,cartoon_oval_width=.28,max_threads=4,stick_radius=.17,valence=0).items():c.set(k,v)
  for k,v in COLORS.items():c.set_color(k,rgb(v))
  c.bg_color('white');c.load(str(OUT/'inputs/1A06_observed_native.pdb'),'autoinhibited');c.load(str(OUT/'inputs/1MXE_model_input_native.pdb'),'CaM_bound');c.dss('all')
  cfg={'autoinhibited':('A','autoinhibited and chain A and resi 294-316','autoinhibited and chain A and not resi 294-316'),'CaM_bound':('B','CaM_bound and chain B','CaM_bound and chain A')}
  manifest['structural_descriptors']={}
  for obj,(chain,reg,background) in cfg.items():
   style(c,obj,chain,reg)
   if obj=='autoinhibited':
    c.set('cartoon_transparency',.43,background);c.set('cartoon_transparency',0,reg)
    contact_sel='autoinhibited and chain A and resi 45+99'
    c.show('sticks',contact_sel+' and not name N+C+O');c.color('background',contact_sel);c.set('stick_radius',.24,contact_sel)
   manifest['structural_descriptors'][obj]={'contacts':contact_metrics(c,obj,chain,background),'SASA':sasa(c,obj,chain,reg)}
  bef=np.array(c.get_coords('autoinhibited'));coords0=bef.copy();bound0=np.array(c.get_coords('CaM_bound'))
  fit=c.pair_fit('autoinhibited and chain A and resi 294-316 and name CA','CaM_bound and chain B and resi 294-316 and name CA');aft=np.array(c.get_coords('autoinhibited'))
  manifest['rigid_fit_23_common_CA_RMSD_A']=float(fit)
  invariant=float(np.max(np.abs(np.linalg.norm(coords0[:,None,:]-coords0[None,:,:],axis=-1)-np.linalg.norm(aft[:,None,:]-aft[None,:,:],axis=-1))))
  manifest['rigid_body_audit']={'autoinhibited_atom_count':len(coords0),'max_pairwise_distance_change_A':invariant,'CaM_bound_coordinates_unchanged':bool(np.array_equal(bound0,c.get_coords('CaM_bound'))),'fit_did_not_rebuild_atoms':True}
  c.save(str(OUT/'inputs/1A06_rigidly_aligned_to_1MXE.pdb'),'autoinhibited')
  c.viewport(1500,700);c.orient('CaM_bound and chain B');c.turn('x',15);c.turn('y',25)
  rot=list(c.get_view()[:9]);manifest['camera_common_rotation']=rot
  for obj,(chain,reg,background) in cfg.items():
   c.disable('all');c.enable(obj);c.zoom(obj,buffer=.8,complete=1);v=list(c.get_view());v[:9]=rot;v[15]=1;v[16]=-v[11]+180;c.set_view(v)
   marks={f'{aa}{site}':f'{obj} and chain {chain} and resi {site} and name '+('CZ2' if site==303 else 'CE') for aa,site in [('W',303),('M',316)]}
   marks.update(regulatory=reg+' and name CA',background=background+' and name CA')
   if obj=='autoinhibited':marks.update(L45='autoinhibited and chain A and resi 45 and not name N+CA+C+O',S99='autoinhibited and chain A and resi 99 and not name N+CA+C+O')
   c.png(str(OUT/(obj+'.png')),width=1500,height=700,dpi=600,ray=1,quiet=1);im=Image.open(OUT/(obj+'.png')).convert('RGBA');bbox=list(im.getchannel('A').getbbox());pad=15;bb=[max(0,bbox[0]-pad),max(0,bbox[1]-pad),min(1500,bbox[2]+pad),min(700,bbox[3]+pad)]
   # Retain full canvas because paired states need a consistent landscape format.
   c.scene(obj,'store');c.save(str(OUT/(obj+'.pse')))
   landmarks={key:project(c,sel,1500,700) for key,sel in marks.items()}
   im.crop(bb).save(OUT/(obj+'_trim.png'))
   manifest['assets'][obj]={'file':obj+'.png','width':1500,'height':700,'view':list(c.get_view()),'alpha_bbox':bbox,'landmarks':landmarks,'trim':{'file':obj+'_trim.png','crop_box':bb,'width':bb[2]-bb[0],'height':bb[3]-bb[1],'landmarks':{key:[xy[0]-bb[0],xy[1]-bb[1]] for key,xy in landmarks.items()}}} 
   print('RENDERED',obj,flush=True)
  manifest['pymol_version']=str(c.get_version()[0]);manifest['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (OUT/'inputs').glob('*.pdb')}
 (OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 sheet=Image.new('RGB',(1500,1400),'white')
 for i,obj in enumerate(cfg):
  im=Image.open(OUT/(obj+'.png')).convert('RGBA');sheet.paste(im,(0,700*i),im)
 sheet.save(OUT/'QA_two_structural_contexts.png')
 print(json.dumps(manifest['structural_descriptors'],indent=2));print('RMSD',fit,'pairwise',invariant)
if __name__=='__main__':main()