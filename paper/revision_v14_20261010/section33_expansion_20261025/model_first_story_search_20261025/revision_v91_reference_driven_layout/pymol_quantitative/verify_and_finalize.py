from pathlib import Path
import os,csv,json,hashlib,shutil,math
import numpy as np
from PIL import Image
from matplotlib import colormaps
ENV=Path('D:/Anaconda/envs/pvb_cpu');os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','');DH=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
OUT=Path(__file__).resolve().parent
m=json.loads((OUT/'render_manifest_final.json').read_text())
rows=list(csv.DictReader((OUT/'inputs/all25_dynamics_and_static_lobe_contacts.csv').open(encoding='utf8')))
PDB=OUT/'inputs/1MXE_CaM_native4_147_CaMKI294_318.pdb'
def sh(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sourceatoms():
 a={}
 for l in PDB.read_text().splitlines():
  if l.startswith('ATOM'):a[(l[21],int(l[22:26]),l[12:16].strip())]=[float(l[30:38]),float(l[38:46]),float(l[46:54])]
 return a
orig=sourceatoms();audit={'purpose':'Independent data/color/geometry check of final 90-degree PSE and PNG assets','n_input_rows':len(rows),'n_RMSF':25,'n_valid_Hchi':sum(r['valid_chi_entropy'].lower()=='true' for r in rows),'source_sha256':{p.name:sh(p) for p in [PDB,OUT/'inputs/all25_dynamics_and_static_lobe_contacts.csv']},'scenes':{},'no_new_model_or_MD':True,'camera_and_scale_identical':np.allclose(m['assets']['1MXE_RMSF_axis090']['view'],m['assets']['1MXE_Hchi_axis090']['view']),'shared_crop_identical':m['assets']['1MXE_RMSF_axis090']['shared_crop_box']==m['assets']['1MXE_Hchi_axis090']['shared_crop_box']}
for metric in ['RMSF','Hchi']:
 e=m['assets'][f'1MXE_{metric}_axis090'];im=Image.open(OUT/e['filename']);bg=Image.new('RGBA',im.size,'white');bg.alpha_composite(im);white=e['name']+'_white_preview.png';bg.convert('RGB').save(OUT/white);e['white_preview_file']=white;e['white_preview_uses']='Standard alpha compositing onto white; original PNG retains native RGBA transparency.'
 alias='1MXE_Full_'+metric
 for ext in ['png','pse']:shutil.copy2(OUT/(e['name']+'.'+ext),OUT/(alias+'.'+ext))
 shutil.copy2(OUT/e['trim_file'],OUT/(alias+'_trim.png'))
 whiteim=Image.open(OUT/white);whiteim.crop(e['shared_crop_box']).save(OUT/(alias+'_trim_white_preview.png'))
 e['delivery_RGBA']=alias+'.png';e['delivery_PSE']=alias+'.pse';e['delivery_trim_RGBA']=alias+'_trim.png';e['delivery_trim_white_preview']=alias+'_trim_white_preview.png';e['png_sha256']=sh(OUT/e['filename'])
 with pymol2.PyMOL() as p:
  c=p.cmd;c.load(str(OUT/e['pse']));atoms=c.get_model('cam or peptide').atom;now={(a.chain,int(a.resi),a.name):np.asarray(a.coord) for a in atoms};maxd=max(float(np.max(np.abs(now[k]-orig[k]))) for k in orig)
  colors=[];c.iterate('peptide and not name N+CA+C+O','colors.append((int(resi),name,color))',space={'colors':colors})
  issues=[];cr=[]
  for row in rows:
   n=int(row['CaMKI_native_position']);valid=row['valid_chi_entropy'].lower()=='true';value=float(row['pred_rmsf' if metric=='RMSF' else 'pred_chi_entropy']);expected=np.asarray(colormaps['cividis'](value/(2.3 if metric=='RMSF' else 4.3))[:3])
   if metric=='Hchi':
    if not valid:expected=np.array([int('A7ADB2'[i:i+2],16)/255 for i in [0,2,4]])
    observed=[np.asarray(c.get_color_tuple(k)) for rn,name,k in colors if rn==n]
    err=max(float(np.max(abs(t-expected))) for t in observed)
    if err>1e-6:issues.append(n)
    cr.append(dict(native_residue=n,valid_Hchi=valid,n_colored_sidechain_heavy_atoms=len(observed),maximum_RGB_error=err))
   else:
    ncols=[];c.iterate(f'peptide and resi {n} and name N','ncols.append(color)',space={'ncols':ncols});raw=ncols[0]
    observed=np.asarray(c.get_color_tuple(int(raw)))
    err=float(np.max(abs(observed-expected)))
    if err>1e-6:issues.append(n)
    cr.append(dict(native_residue=n,base_N_atom_color_id=raw,maximum_RGB_error=err))
  audit['scenes'][metric]={'rendered_atoms_same_keys_as_source':set(now)==set(orig),'maximum_abs_coordinate_difference_from_decimal_PDB_A':maxd,'all_coordinates_within_float32_roundoff_5e-6A':maxd<5e-6,'source_vs_saved_loaded_geometry_identical':m['geometry_audit']['original_atom_coordinates_identical'],'all_25_CA_present':c.count_atoms('peptide and name CA')==25,'all_144_CaM_CA_present':c.count_atoms('cam and name CA')==144,'residue_color_mismatches':issues,'per_residue_color_audit':cr,'peptide_backbone_Hchi_color':'#9DA4AA' if metric=='Hchi' else 'cividis scalar per residue','peptide_cartoon_transparency':float(c.get('cartoon_transparency','peptide')),'CaM_cartoon_transparency':float(c.get('cartoon_transparency','cam')),'all_21_valid_sidechains_selected':metric!='Hchi' or len(set(rn for rn,an,cl in colors if rn not in [299,306,309,311]))==21}
  assert set(now)==set(orig) and maxd<5e-6 and not issues
for row,entry in zip(rows,m['residue_mapping']):
 assert int(row['CaMKI_native_position'])==entry['native_residue']
 assert float(row['pred_rmsf'])==entry['RMSF']
 if entry['Hchi_valid']:assert float(row['pred_chi_entropy'])==entry['Hchi']
assert audit['camera_and_scale_identical'] and audit['shared_crop_identical']
m['selected_angle_degrees']=90;m['rendering_note']='Actual observed 1MXE peptide, N left/C right, rigid camera only. Pale CaM cartoon alpha 0.16 context. Hchi backbone neutral #9DA4AA at opacity 0.65. 4 alanines Hchi undefined and neutral gray. Full 25/21 values used.'
m['geometry_and_data_audit_file']='geometry_data_color_audit.json';m['delivery_assets']={metric:{'rgba':m['assets'][f'1MXE_{metric}_axis090']['delivery_trim_RGBA'],'white_preview':m['assets'][f'1MXE_{metric}_axis090']['delivery_trim_white_preview'],'pse':m['assets'][f'1MXE_{metric}_axis090']['delivery_PSE'],'full_frame_landmarks':m['assets'][f'1MXE_{metric}_axis090']['landmarks'],'cropped_landmarks':m['assets'][f'1MXE_{metric}_axis090']['trim_landmarks']} for metric in ['RMSF','Hchi']}
(OUT/'render_manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf8');(OUT/'geometry_data_color_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k!='scenes'},ensure_ascii=False,indent=2));print('FINAL COMPLETE')