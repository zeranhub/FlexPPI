# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# Native molecular overview and residue-neighborhood sticks -> V91 PyMOL style inheritance.
# Exact observed atoms only; identical rigid rotation for overview and both close-ups.
from pathlib import Path
import csv, json, os, math, hashlib
import numpy as np
from PIL import Image
ENV=Path('D:/Anaconda/envs/pvb_cpu');os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','');DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
OUT=Path(__file__).resolve().parent;INP=OUT/'inputs';PDB=INP/'1MXE_CaM_native4_147_CaMKI294_318.pdb'
CAM=json.loads((INP/'base_camera_rotation.json').read_text());R=np.asarray(CAM['rotation']).reshape(3,3);FOV=20.
SCALE=int(os.environ.get('V92_RENDER_SCALE','1'));OW,OH=1400*SCALE,1077*SCALE;ZW,ZH=1400*SCALE,1232*SCALE
COLORS={'CaM':'#BAC1C7','peptide':'#8B72A5','context':'#D5DADD','neighbor':'#B9C0C5','target':'#3A8D8C'}
AA={'ILE':'I','LYS':'K','ASN':'N','PHE':'F','ALA':'A','SER':'S','TRP':'W','GLN':'Q','THR':'T','VAL':'V','ARG':'R','HIS':'H','MET':'M'}

def rgb(h):return [int(h[i:i+2],16)/255 for i in [1,3,5]]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def projected(xyz,view,w,h):
 v=np.asarray(view);q=(np.asarray(xyz)-v[12:15])@v[:9].reshape(3,3)+v[9:12];hh=-v[11]*math.tan(math.radians(FOV/2));hw=hh*w/h
 return np.column_stack([(0.5+q[:,0]/(2*hw))*w,(0.5-q[:,1]/(2*hh))*h])
def makeview(center,halfheight):
 d=halfheight/math.tan(math.radians(FOV/2));return list(R.ravel())+[0.,0.,-float(d)]+[float(x) for x in center]+[1.,float(d+200),-20.]
def setup(c):
 for k,v in dict(ray_opaque_background=0,opaque_background=0,ray_trace_mode=0,ray_shadows=0,antialias=1,orthoscopic=1,field_of_view=FOV,depth_cue=0,specular=.12,shininess=30,ambient=.72,direct=.30,reflect=.08,light_count=2,cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,cartoon_smooth_loops=1,cartoon_sampling=24,cartoon_loop_radius=.22,cartoon_tube_radius=.22,cartoon_oval_length=1.25,cartoon_oval_width=.25,max_threads=4,valence=0,transparency_mode=2,two_sided_lighting=1,stick_quality=24).items():c.set(k,v)
 for k,v in COLORS.items():c.set_color(k,rgb(v))
 c.bg_color('white');c.load(str(PDB),'bound');c.dss('bound');c.hide('everything','bound')
def landmarks(c,view,w,h):
 lm={}
 for n in range(294,319):
  ss=f'bound and chain B and resi {n}';ca=c.get_coords(ss+' and name CA');sc=c.get_coords(ss+' and not name N+CA+C+O');name=AA[c.get_model(ss+' and name CA').atom[0].resn]+str(n)
  lm[name]={'CA':projected(ca.mean(0)[None,:],view,w,h)[0].tolist(),'sidechain_center':projected(sc.mean(0)[None,:],view,w,h)[0].tolist(),'CA_xyz_A':ca.mean(0).tolist(),'sidechain_center_xyz_A':sc.mean(0).tolist()}
 return lm
def render(c,key,view,w,h,scope,selection='bound'):
 c.set_view(view);c.viewport(w,h);c.png(str(OUT/(key+'.png')),width=w,height=h,dpi=800,ray=1,quiet=1);c.save(str(OUT/(key+'.pse')))
 im=Image.open(OUT/(key+'.png')).convert('RGBA');bg=Image.new('RGBA',im.size,'white');bg.alpha_composite(im);bg.convert('RGB').save(OUT/(key+'_white_preview.png'))
 px=projected(c.get_coords(selection),view,w,h);inside=(px[:,0]>=0)&(px[:,0]<w)&(px[:,1]>=0)&(px[:,1]<h)
 assert inside.all(),(key,'selected atom projection clipped')
 return {'png':key+'.png','white_preview':key+'_white_preview.png','pse':key+'.pse','width':w,'height':h,'view':list(c.get_view()),'crop':[0,0,w,h],'field_of_view':FOV,'named_landmarks':landmarks(c,view,w,h),'all_displayed_scope_atom_projections_inside_frame':bool(inside.all()),'minimum_projection_margin_px':{'left':float(px[:,0].min()),'right':float(w-px[:,0].max()),'top':float(px[:,1].min()),'bottom':float(h-px[:,1].max())},'scope':scope,'PNG_sha256':sha(OUT/(key+'.png'))}

def main():
 manifest={'source':{'PDB':'inputs/'+PDB.name,'PDB_sha256':sha(PDB),'camera_rotation_source':'inputs/base_camera_rotation.json'},'assets':{},'neighborhood_definition':'For each target N297/N308: all observed protein residues having any heavy atom <=6.0 A from any heavy sidechain atom of that target. All heavy atoms of each selected residue are shown as sticks; selection includes both CaM chain A and CaMKI peptide chain B. Whole residues are never truncated to only contact atoms.','colors':COLORS,'overview_data_mapping':'None: chain/component colors only. Quantitative predictions belong to panel c.','render_scale':SCALE,'geometry_audit':{}}
 with pymol2.PyMOL() as p:
  c=p.cmd;setup(c);src=c.get_coords('bound').copy();model=c.get_model('bound');atoms=model.atom;xyz=np.asarray([a.coord for a in atoms]);neighs={}
  # Determine selected whole residues once from the exact source coordinates.
  for target in [297,308]:
   sc=c.get_coords(f'bound and chain B and resi {target} and not name N+CA+C+O');center=sc.mean(0);dist=np.linalg.norm(xyz[:,None,:]-sc[None,:,:],axis=2).min(1);resids=sorted(set((a.chain,int(a.resi),a.resn) for a,d in zip(atoms,dist) if d<=6.0),key=lambda x:(x[0],x[1]));clauses=[]
   for ch in sorted(set(t[0] for t in resids)):clauses.append(f"(chain {ch} and resi {'+'.join(str(t[1]) for t in resids if t[0]==ch)})")
   sel='bound and ('+' or '.join(clauses)+')';near=c.get_coords(sel);q=(near-center)@R
   neighs[target]={'selection':sel,'center':center,'coords':near,'whole_residues':[{'chain':ch,'native_residue':n,'amino_acid_3letter':rn,'n_heavy_atoms':c.count_atoms(f'bound and chain {ch} and resi {n}')} for ch,n,rn in resids],'n_CaM_residues':sum(ch=='A' for ch,n,rn in resids),'n_peptide_residues':sum(ch=='B' for ch,n,rn in resids),'all_complete_residue_heavy_atoms':len(near),'minimum_CaM_distance_to_target_sidechain_A':float(dist[np.array([a.chain=='A' for a in atoms])].min()),'needed_half_height_A':float(max(abs(q[:,1]).max(),abs(q[:,0]).max()/(ZW/ZH))+1.6)}
  halfheight=max(n['needed_half_height_A'] for n in neighs.values());manifest['local_shared_half_height_A']=halfheight;manifest['local_pixels_per_A_at_height']=ZH/(2*halfheight)
  # Overview: all-cartoon, no isolated sticks, whole complex is within the frame.
  c.show('cartoon','bound');c.color('CaM','bound and chain A');c.color('peptide','bound and chain B');c.set('cartoon_color','CaM','bound and chain A');c.set('cartoon_color','peptide','bound and chain B');c.set('cartoon_transparency',.32,'bound and chain A');c.set('cartoon_transparency',0.,'bound and chain B')
  center=xyz.mean(0);q=(xyz-center)@R;shift=(q.min(0)+q.max(0))/2;center=center+R@shift;q=(xyz-center)@R;ohalf=max(abs(q[:,1]).max(),abs(q[:,0]).max()/(OW/OH))+1.5;view=makeview(center,ohalf)
  manifest['assets']['overview']=render(c,'overview',view,OW,OH,'Complete 1MXE CaM + complete 25-residue CaMKI peptide. All cartoon, zero sticks; conventional component colors only.')
  manifest['assets']['overview'].update(CaM_cartoon_transparency=.32,peptide_cartoon_transparency=0.,stick_count=0)
  assert np.array_equal(src,c.get_coords('bound'))
  for target in [297,308]:
   n=neighs[target];sel=n['selection'];t=f'bound and chain B and resi {target}';c.hide('everything','bound');c.show('cartoon','bound');c.color('context','bound');c.set('cartoon_color','context','bound and chain A');c.set('cartoon_color','context','bound and chain B');c.set('cartoon_transparency',.86,'bound and chain A');c.set('cartoon_transparency',.86,'bound and chain B');c.show('sticks',sel);c.color('neighbor',sel);c.color('target',t);c.set('stick_radius',.12,sel);c.set('stick_radius',.19,t);c.set('stick_transparency',0.,'bound')
   # All atoms of target Asn use the same teal; no element coloring is applied.
   key=f'zoom_N{target}';view=makeview(n['center'],halfheight);entry=render(c,key,view,ZW,ZH,f'6 A whole-residue neighborhood centered on N{target}; full selected residues are sticks on faint cartoon background. Distal cartoon may extend beyond this local frame.',sel)
   entry.update(target=f'N{target}',selection=sel,cutoff_A=6.0,selection_reference_atoms='target sidechain heavy atoms; N/CA/C/O excluded',whole_residues=n['whole_residues'],n_CaM_residues=n['n_CaM_residues'],n_peptide_residues=n['n_peptide_residues'],n_all_complete_residue_heavy_atoms=n['all_complete_residue_heavy_atoms'],minimum_CaM_distance_to_target_sidechain_A=n['minimum_CaM_distance_to_target_sidechain_A'],background_cartoon_transparency=.86,all_neighbor_stick_radius_A=.12,target_stick_radius_A=.19,all_selected_residue_heavy_atoms_in_stick_selection=True,no_surface=True,shared_optical_half_height_A=halfheight)
   manifest['assets'][key]=entry;assert np.array_equal(src,c.get_coords('bound'));print('RENDERED',key,flush=True)
  manifest['geometry_audit']={'all_original_loaded_atom_coordinates_identical_after_all_renders':bool(np.array_equal(src,c.get_coords('bound'))),'maximum_coordinate_change_A':float(abs(src-c.get_coords('bound')).max()),'source_heavy_atom_count':len(src),'CaM_CA_count':c.count_atoms('bound and chain A and name CA'),'CaMKI_CA_count':c.count_atoms('bound and chain B and name CA'),'local_rotation_identical':bool(np.array_equal(manifest['assets']['zoom_N297']['view'][:9],manifest['assets']['zoom_N308']['view'][:9])),'local_optical_distance_identical':manifest['assets']['zoom_N297']['view'][11]==manifest['assets']['zoom_N308']['view'][11],'no_new_model_or_MD':True}
 (OUT/'render_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8');print('COMPLETE',flush=True)
if __name__=='__main__':main()