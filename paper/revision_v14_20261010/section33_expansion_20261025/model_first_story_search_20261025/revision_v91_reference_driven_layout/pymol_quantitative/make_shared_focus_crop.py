from pathlib import Path
import json,math,hashlib
import numpy as np
from PIL import Image
OUT=Path(__file__).resolve().parent
m=json.loads((OUT/'render_manifest.json').read_text());audit=json.loads((OUT/'geometry_data_color_audit.json').read_text())
e0=m['assets']['1MXE_Hchi_axis090'];v=np.array(e0['view']);R=v[:9].reshape(3,3);xyz=[];atoms=[]
for l in (OUT/'inputs/1MXE_CaM_native4_147_CaMKI294_318.pdb').read_text().splitlines():
 if l.startswith('ATOM') and l[21]=='B':
  xyz.append([float(l[30:38]),float(l[38:46]),float(l[46:54])]);atoms.append({'native_residue':int(l[22:26]),'atom_name':l[12:16].strip()})
xyz=np.array(xyz);q=(xyz-v[12:15])@R+v[9:12];w,h=e0['width'],e0['height'];hh=-v[11]*math.tan(math.radians(e0['field_of_view']/2));hw=hh*w/h
px=np.column_stack([(0.5+q[:,0]/(2*hw))*w,(0.5-q[:,1]/(2*hh))*h])
lo=px.min(0)-70;hi=px.max(0)+70
# Widen vertical context to a 1.95 ratio if that only adds context; never clip any peptide atom.
width=hi[0]-lo[0];height=hi[1]-lo[1]
if width/height>1.95:
 add=(width/1.95-height)/2;lo[1]-=add;hi[1]+=add
crop=[max(0,int(np.floor(lo[0]))),max(0,int(np.floor(lo[1]))),min(w,int(np.ceil(hi[0]))),min(h,int(np.ceil(hi[1])))]
inside=(px[:,0]>=crop[0])&(px[:,0]<crop[2])&(px[:,1]>=crop[1])&(px[:,1]<crop[3]);assert inside.all()
for metric in ['RMSF','Hchi']:
 e=m['assets'][f'1MXE_{metric}_axis090'];rgba=Image.open(OUT/e['filename']).convert('RGBA');f=f'1MXE_Full_{metric}_focus.png';rgba.crop(crop).save(OUT/f)
 im=rgba.crop(crop);bg=Image.new('RGBA',im.size,'white');bg.alpha_composite(im);wf=f'1MXE_Full_{metric}_focus_white_preview.png';bg.convert('RGB').save(OUT/wf)
 lm={n:{k:[p[0]-crop[0],p[1]-crop[1]] for k,p in r.items() if k in ['CA','sidechain_center']} for n,r in e['landmarks'].items()}
 e['focus_file']=f;e['focus_white_preview_file']=wf;e['focus_crop_box']=crop;e['focus_width']=crop[2]-crop[0];e['focus_height']=crop[3]-crop[1];e['focus_landmarks']=lm
 e['focus_scope']='Local CaM support around the full 25-residue peptide. Distal pale CaM is intentionally outside crop; not a whole-complex overview.'
 m['delivery_assets'][metric].update(rgba_focus=f,white_preview_focus=wf,focus_landmarks=lm,focus_crop_box=crop,focus_width=e['focus_width'],focus_height=e['focus_height'])
 border={'left':float((px[:,0]-crop[0]).min()),'right':float((crop[2]-px[:,0]).min()),'top':float((px[:,1]-crop[1]).min()),'bottom':float((crop[3]-px[:,1]).min())}
 audit['focus_crop']={'shared_exact_crop_box':crop,'width_px':crop[2]-crop[0],'height_px':crop[3]-crop[1],'aspect':(crop[2]-crop[0])/(crop[3]-crop[1]),'entire_25_residue_peptide_heavy_atom_projection_retained':bool(inside.all()),'n_peptide_heavy_atoms_tested':len(atoms),'n_native_residues_tested':len(set(a['native_residue'] for a in atoms)),'minimum_atom_projection_margins_px':border,'whole_complex_overview':False,'distal_CaM_crop_allowed_and_declared':True,'no_rescaling_or_warping_of_native_render':True}
m['focus_scope_note']='Preferred horizontal crop retains every real peptide heavy atom projection plus >=70 px margins. Only distal CaM context is excluded; same native camera, pixel scale and rectangle for both metrics.'
(OUT/'render_manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf8');(OUT/'geometry_data_color_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps(audit['focus_crop'],indent=2))