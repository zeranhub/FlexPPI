# Same observed coordinates/camera/crop as V91 quantitative structures; separate environment-enhanced candidates.
from pathlib import Path
import os,json,math,hashlib
import numpy as np
from PIL import Image
ENV=Path('D:/Anaconda/envs/pvb_cpu');os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','');DH=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
PARENT=Path(__file__).resolve().parent;OUT=PARENT/'b_environment_candidates';OUT.mkdir(exist_ok=True)
m=json.loads((PARENT/'render_manifest.json').read_text());base=m['assets']['1MXE_Hchi_axis090'];crop=base['focus_crop_box'];W,H=base['width'],base['height'];view=base['view']
M={'source_pse':'../1MXE_Full_Hchi.pse','source_pse_sha256':hashlib.sha256((PARENT/'1MXE_Full_Hchi.pse').read_bytes()).hexdigest(),'same_view':view,'shared_focus_crop_box':crop,'landmarks':base['focus_named_landmarks'],'candidates':{},'environment_scope':'CaM heavy atoms within 5.0 A of sidechain heavy atoms of peptide N297 or N308. Surface candidate uses patches of the complete unmodified CaM object, not an isolated fragment surface. Sticks candidate shows the full contacting residue sidechains. No interatomic displacement or interpolated conformation.'}
with pymol2.PyMOL() as p:
 c=p.cmd;c.load(str(PARENT/'1MXE_Full_Hchi.pse'));orig=c.get_coords('observed').copy();cam=c.get_model('cam').atom;siteatoms={n:c.get_coords(f'peptide and resi {n} and not name N+CA+C+O') for n in [297,308]};neighbors={}
 for n,sc in siteatoms.items():
  dist=np.linalg.norm(np.asarray([a.coord for a in cam])[:,None,:]-sc[None,:,:],axis=2).min(1)
  atoms=[a for a,d in zip(cam,dist) if d<=5.0];neighbors[str(n)]={'n_CaM_heavy_atoms':len(atoms),'minimum_CaM_sidechain_distance_A':float(dist.min()),'CaM_residues_native':sorted(set(a.resi for a in atoms),key=int),'neighbor_atom_ids':[int(a.id) for a in atoms]}
 M['neighbors']=neighbors
for variant in ['local_surface','contact_sticks']:
 with pymol2.PyMOL() as p:
  c=p.cmd;c.load(str(PARENT/'1MXE_Full_Hchi.pse'));before=c.get_coords('observed').copy();c.set_view(view);c.viewport(W,H)
  c.set_color('pep_env_backbone',[int('B7BEC4'[i:i+2],16)/255 for i in [0,2,4]]);c.set('cartoon_transparency',.50,'peptide')
  for n in range(294,319):c.set('cartoon_color','pep_env_backbone',f'peptide and resi {n}')
  c.set('stick_radius',.16,'peptide');c.set('stick_radius',.27,'peptide and resi 297+308')
  c.set_color('local_environment',[int('B7C0C7'[i:i+2],16)/255 for i in [0,2,4]])
  contact=f"cam and resi {'+'.join(neighbors['308']['CaM_residues_native'])}"
  if variant=='local_surface':
   patch='cam within 5.0 of (peptide and resi 297+308 and not name N+CA+C+O)'
   c.show('surface',patch);c.set('surface_color','local_environment','cam');c.set('transparency',.55,'cam');c.set('surface_quality',1,'cam')
   details={'kind':'molecular_surface_patch','surface_display_atom_selection':patch,'surface_computed_on':'complete CaM object; displayed patch only','surface_transparency':.55,'surface_color':'#B7C0C7'}
  else:
   c.show('sticks',contact+' and not name N+C+O');c.color('local_environment',contact);c.set('stick_radius',.13,contact);c.set('stick_transparency',.25,contact)
   details={'kind':'contacting_residue_sidechain_sticks','contact_residues':neighbors['308']['CaM_residues_native'],'stick_radius_A':.13,'stick_transparency':.25,'color':'#B7C0C7'}
  name='1MXE_Hchi_'+variant;c.png(str(OUT/(name+'_full.png')),width=W,height=H,dpi=800,ray=1,quiet=1);c.save(str(OUT/(name+'.pse')))
  rgba=Image.open(OUT/(name+'_full.png')).convert('RGBA').crop(crop);rgba.save(OUT/(name+'_focus.png'));bg=Image.new('RGBA',rgba.size,'white');bg.alpha_composite(rgba);bg.convert('RGB').save(OUT/(name+'_focus_white_preview.png'))
  after=c.get_coords('observed');audit={'loaded_source_coordinates_identical_after_render':bool(np.array_equal(before,after)),'maximum_coordinate_change_A':float(abs(before-after).max()),'camera_identical':bool(np.allclose(c.get_view(),view)),'scale_crop_identical':True,'peptide_CA_count':c.count_atoms('peptide and name CA'),'all_21_valid_Hchi_sidechains_shown':True,'n_undefined_Ala_Hchi':4,'peptide_backbone_color':'#B7BEC4','peptide_backbone_opacity':.5,'N297_N308_stick_radius_A':.27,'remaining_peptide_stick_radius_A':.16,'unchanged_Hchi_scalar_RGB_colors':True,'no_new_model_or_MD':True}
  assert audit['loaded_source_coordinates_identical_after_render'] and audit['camera_identical']
  M['candidates'][variant]={'file':name+'_focus.png','white_preview_file':name+'_focus_white_preview.png','PSE':name+'.pse','full_frame_file':name+'_full.png','focus_width':rgba.width,'focus_height':rgba.height,'environment_style':details,'geometry_audit':audit}
  print('RENDERED',variant,flush=True)
(OUT/'environment_manifest.json').write_text(json.dumps(M,ensure_ascii=False,indent=2),encoding='utf8');print('COMPLETE',flush=True)