# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# Molecular quantitative structure -> prior V90 PyMOL renderer -> native PyMOL style inheritance
# Data values are mapped to exact observed coordinates; no shape interpolation, no model/MD run.
from pathlib import Path
import csv, json, os, hashlib, math
import numpy as np
from PIL import Image
from matplotlib import colormaps
ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
OUT=Path(__file__).resolve().parent
INP=OUT/'inputs'; PDB=INP/'1MXE_CaM_native4_147_CaMKI294_318.pdb'
CSV=INP/'all25_dynamics_and_static_lobe_contacts.csv'
ROWS=list(csv.DictReader(CSV.open(encoding='utf-8-sig')))
KEYS=[294,296,297,298,303,304,308,316,318]
PREVIEW=os.environ.get("V91_PREVIEW","0")=="1"
WIDTH,HEIGHT=(1300,800) if PREVIEW else (2600,1600)
ANGLES=[int(x) for x in os.environ.get("V91_ANGLES","0,45,90,135").split(",")]
FOV=20.0
SURFACE_TRANSPARENCY=0.92
CARTOON_TRANSPARENCY=0.84
STYLES={'CaM':'#D2D7DA','backbone':'#9DA4AA','undefined':'#A7ADB2','anchor':'#A7ADB2'}

def rgb(h):return [int(h[i:i+2],16)/255 for i in (1,3,5)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def projection(c,coords,w=WIDTH,h=HEIGHT):
 xyz=np.asarray(coords).reshape(-1,3).mean(0);v=np.asarray(c.get_view());q=v[:9].reshape(3,3).T@(xyz-v[12:15])+v[9:12]
 hh=-v[11]*math.tan(math.radians(float(c.get('field_of_view'))/2));hw=hh*w/h
 return [float((.5+q[0]/(2*hw))*w),float((.5-q[1]/(2*hh))*h)]
def selcent(c,sel):
 a=c.get_coords(sel)
 return None if a is None else np.asarray(a).mean(0)
def set_style(c):
 for k,v in dict(ray_opaque_background=0,opaque_background=0,ray_trace_mode=0,ray_shadows=0,antialias=1,orthoscopic=1,field_of_view=FOV,depth_cue=0,specular=.12,shininess=30,ambient=.72,direct=.30,reflect=.08,light_count=2,cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,cartoon_smooth_loops=1,cartoon_sampling=24,cartoon_loop_radius=.22,cartoon_tube_radius=.22,cartoon_oval_length=1.25,cartoon_oval_width=.25,max_threads=4,valence=0,transparency_mode=2,two_sided_lighting=1,surface_quality=1,stick_quality=24).items(): c.set(k,v)
 for k,v in STYLES.items():c.set_color(k,rgb(v))
 c.bg_color('white')
 c.hide('everything','all');c.show('cartoon','cam')
 c.color('CaM','cam');c.set('cartoon_color','CaM','cam');c.set('transparency',SURFACE_TRANSPARENCY,'cam');c.set('cartoon_transparency',CARTOON_TRANSPARENCY,'cam')
 c.show('cartoon','peptide');c.set('cartoon_transparency',0,'peptide');c.set('cartoon_color','backbone','peptide');c.color('backbone','peptide')

def metricstyle(c,metric):
 c.hide('sticks','peptide');c.color('backbone','peptide');c.set('cartoon_color','backbone','peptide');c.set('cartoon_transparency',0.35 if metric=='Hchi' else 0.,'peptide')
 col='pred_rmsf' if metric=='RMSF' else 'pred_chi_entropy';vmax=2.3 if metric=='RMSF' else 4.3
 for row in ROWS:
  n=int(row['CaMKI_native_position']);aa=row['amino_acid'];valid=row['valid_chi_entropy'].lower()=='true';c.set('cartoon_color','backbone',f'peptide and resi {n}')
  val=float(row[col]) if metric=='RMSF' or valid else None
  ss=f'peptide and resi {n}'
  if val is not None:
   cc=f'{metric}_{n}';color=list(colormaps['cividis'](val/vmax)[:3]);c.set_color(cc,color)
   if metric=='RMSF': c.set('cartoon_color',cc,ss);c.color(cc,ss)
   else: c.color(cc,ss+' and not name N+C+O')
  elif metric=='Hchi':c.color('undefined',ss+' and not name N+C+O')
 if metric=='Hchi':
  c.show('sticks','peptide and not name N+C+O');c.set('stick_radius',.21,'peptide');c.set('stick_ball',0,'peptide')
 else:
  c.show('sticks','peptide and resi 303+316 and not name N+C+O');c.color('anchor','peptide and resi 303+316 and not name N+C+O');c.set('stick_radius',.12,'peptide')
 return dict(value_column=col,minimum=0,maximum=vmax,colormap='matplotlib cividis',units='A' if metric=='RMSF' else 'nats')

def camera(c,angle):
 ca=c.get_coords('peptide and name CA');center=ca.mean(0);x=np.linalg.svd(ca-center)[2][0]
 if np.dot(x,ca[-3:].mean(0)-ca[:3].mean(0))<0:x=-x
 g=c.get_coords('cam and name CA').mean(0)-center;y=g-x*np.dot(g,x);y/=np.linalg.norm(y);z=np.cross(x,y);z/=np.linalg.norm(z)
 a=math.radians(angle);yr=y*math.cos(a)+z*math.sin(a);zr=np.cross(x,yr)
 R=np.column_stack([x,yr,zr]);allxyz=c.get_coords('cam or peptide');q=(allxyz-center)@R
 # Fixed extent includes the full CaM context but deliberately centers the peptide.
 # Entire untransformed structure is present, and both metric renderings reuse this exact view.
 hspan=max(np.max(abs(q[:,1]))*1.05,np.max(abs(q[:,0]))/(WIDTH/HEIGHT)*1.05)
 dist=hspan/math.tan(math.radians(FOV/2))
 v=list(R.ravel())+[0.,0.,-float(dist)]+list(center)+[1.,float(dist+150),-20.]
 c.set_view(v);return v,{'long_axis_vector':x.tolist(),'y_vector':yr.tolist(),'z_vector':zr.tolist(),'center':center.tolist(),'angle_about_peptide_axis_degrees':angle,'half_height_A':float(hspan)}

def render(c,name,metric,view,info):
 c.set_view(view);c.viewport(WIDTH,HEIGHT);metricinfo=metricstyle(c,metric)
 c.png(str(OUT/(name+'.png')),width=WIDTH,height=HEIGHT,dpi=800,ray=1,quiet=1)
 c.save(str(OUT/(name+'.pse')))
 im=Image.open(OUT/(name+'.png')).convert('RGBA')
 bb=im.getchannel('A').getbbox();pad=35;crop=[max(0,bb[0]-pad),max(0,bb[1]-pad),min(WIDTH,bb[2]+pad),min(HEIGHT,bb[3]+pad)]
 # Shared crop is assigned once after both metrics render; full original PNG remains the registration master.
 lm={}
 for row in ROWS:
  n=int(row['CaMKI_native_position']);ss=f'peptide and resi {n}';ca=selcent(c,ss+' and name CA');sc=selcent(c,ss+' and not name N+CA+C+O');lm[str(n)]={'name':row['amino_acid']+str(n),'CA':projection(c,ca),'sidechain_center':projection(c,sc),'CA_xyz_A':ca.tolist(),'sidechain_center_xyz_A':sc.tolist(),'sidechain_heavy_atom_count':int(c.count_atoms(ss+' and not name N+CA+C+O'))}
 entry=dict(name=name,metric=metric,filename=name+'.png',pse=name+'.pse',width=WIDTH,height=HEIGHT,view=list(c.get_view()),field_of_view=float(c.get('field_of_view')),alpha_bbox=list(bb),suggested_crop=crop,metric_mapping=metricinfo,camera_construction=info,landmarks=lm,alpha_mode='transparent background; context uses actual partial alpha',geometry='original atom coordinates retained')
 print('RENDERED',name,flush=True);return entry

def main():
 manifest={'purpose':'Full model RMSF and sidechain Hchi on actual 1MXE peptide, same observed bound coordinates, no inferred ensemble or free-state structure','source_files':{str(p.name):sha(p) for p in [PDB,CSV]},'settings':dict(CaM_surface_transparency=SURFACE_TRANSPARENCY,CaM_cartoon_transparency=CARTOON_TRANSPARENCY,context_colors=STYLES),'assets':{},'residue_mapping':[]}
 assert len(ROWS)==25 and [int(r['CaMKI_native_position']) for r in ROWS]==list(range(294,319))
 assert sum(r['valid_chi_entropy'].lower()=='true' for r in ROWS)==21
 for row in ROWS:
  n=int(row['CaMKI_native_position']);valid=row['valid_chi_entropy'].lower()=='true'
  manifest['residue_mapping'].append(dict(native_residue=n,amino_acid=row['amino_acid'],RMSF=float(row['pred_rmsf']),Hchi=float(row['pred_chi_entropy']) if valid else None,Hchi_valid=valid,RMSF_RGB=list(colormaps['cividis'](float(row['pred_rmsf'])/2.3)[:3]),Hchi_RGB=list(colormaps['cividis'](float(row['pred_chi_entropy'])/4.3)[:3]) if valid else rgb(STYLES['undefined'])))
 with pymol2.PyMOL() as p:
  c=p.cmd;c.load(str(PDB),'observed');original=c.get_coords('observed').copy();c.create('cam','observed and chain A');c.create('peptide','observed and chain B');c.disable('observed');c.dss('cam or peptide');set_style(c)
  assert c.count_atoms('peptide and name CA')==25
  assert c.count_atoms('cam and name CA')==144
  for angle in ANGLES:
   view,info=camera(c,angle)
   for metric in ['RMSF','Hchi']:
    name=f'1MXE_{metric}_axis{angle:03d}'
    manifest['assets'][name]=render(c,name,metric,view,info)
  current=c.get_coords('observed');manifest['geometry_audit']={'observed_atoms':len(original),'original_atom_coordinates_identical':bool(np.array_equal(original,current)),'max_coordinate_deviation_A':float(np.max(abs(original-current))),'copies_coordinates_identical':bool(np.array_equal(np.vstack([c.get_coords('cam'),c.get_coords('peptide')]),original)),'PDB_CA_count_CaM':144,'PDB_CA_count_CaMKI':25,'native_sequence_294_318':''.join(r['amino_acid'] for r in ROWS),'n_RMSF_values':25,'n_valid_Hchi_values':21,'Hchi_undefined_residues':[int(r['CaMKI_native_position']) for r in ROWS if r['valid_chi_entropy'].lower()!='true'],'same_metric_view_each_angle':all(np.allclose(manifest['assets'][f'1MXE_RMSF_axis{a:03d}']['view'],manifest['assets'][f'1MXE_Hchi_axis{a:03d}']['view']) for a in ANGLES),'no_new_inference_or_MD':True}
  assert manifest['geometry_audit']['original_atom_coordinates_identical'] and manifest['geometry_audit']['copies_coordinates_identical']
 for angle in ANGLES:
  es=[manifest['assets'][f'1MXE_{m}_axis{angle:03d}'] for m in ['RMSF','Hchi']]
  crop=[min(e['alpha_bbox'][0] for e in es)-25,min(e['alpha_bbox'][1] for e in es)-25,max(e['alpha_bbox'][2] for e in es)+25,max(e['alpha_bbox'][3] for e in es)+25];crop=[max(0,crop[0]),max(0,crop[1]),min(WIDTH,crop[2]),min(HEIGHT,crop[3])]
  for e in es:
   e['shared_crop_box']=crop;e['trim_file']=e['name']+'_trim.png';e['trim_width']=crop[2]-crop[0];e['trim_height']=crop[3]-crop[1];e['trim_landmarks']={n:{key:[p[0]-crop[0],p[1]-crop[1]] for key,p in d.items() if key in ['CA','sidechain_center']} for n,d in e['landmarks'].items()};Image.open(OUT/e['filename']).crop(crop).save(OUT/e['trim_file'])
 (OUT/('render_manifest_all_candidates.json' if PREVIEW else 'render_manifest_final.json')).write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8')
 print('COMPLETE',flush=True)
if __name__=='__main__':main()