from pathlib import Path
import os,json,hashlib,math
import numpy as np
ENV=Path('D:/Anaconda/envs/pvb_cpu');os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','');DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
OUT=Path(__file__).resolve().parent;PDB=OUT/'inputs/1MXE_CaM_native4_147_CaMKI294_318.pdb'
def rgb(h):return np.asarray([int(h[i:i+2],16)/255 for i in [1,3,5]])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atoms(c,sel):
 return {(a.chain,int(a.resi),a.name):np.asarray(a.coord) for a in c.get_model(sel).atom}
def check():
 m=json.loads((OUT/'render_manifest.json').read_text());A={'source_PDB_sha256':sha(PDB),'renderer_sha256':sha(OUT/'render_v92_structures.py'),'assets':{},'scientific_scope':'Experimental 1MXE observed geometry; structural localization and same-rule residue environments only. No RMSF/Hchi is encoded in these images. No new model, MD, free state, or mechanistic path.'}
 with pymol2.PyMOL() as p:
  c=p.cmd;c.load(str(PDB),'source');src=atoms(c,'source');model=c.get_model('source').atom;xyz=np.asarray([a.coord for a in model]);sc={n:c.get_coords(f'source and chain B and resi {n} and not name N+CA+C+O') for n in [297,308]};expected={}
  for n in [297,308]:
   d=np.linalg.norm(xyz[:,None,:]-sc[n][None,:,:],axis=2).min(1);resids=set((a.chain,int(a.resi)) for a,dd in zip(model,d) if dd<=6.0);expected[n]={k for k in src if k[:2] in resids}
 for key,e in m['assets'].items():
  with pymol2.PyMOL() as p:
   c=p.cmd;c.load(str(OUT/e['pse']));now=atoms(c,'bound');reps=atoms(c,'bound and rep sticks');colors=[];c.iterate('bound','colors.append((chain,int(resi),name,color))',space={'colors':colors});colorMap={(ch,n,name):np.asarray(c.get_color_tuple(col)) for ch,n,name,col in colors}
   geometry=(set(src)==set(now) and all(np.array_equal(src[k],now[k]) for k in src));maxd=max(float(abs(src[k]-now[k]).max()) for k in src);record={'all_source_atoms_retained':set(src)==set(now),'source_heavy_atoms':len(src),'loaded_source_coordinates_exactly_equal':geometry,'maximum_coordinate_change_A':maxd,'surface_atom_count':c.count_atoms('bound and rep surface'),'same_saved_camera':bool(np.array_equal(c.get_view(),e['view'])),'PNG_sha256':sha(OUT/e['png']),'PSE_sha256':sha(OUT/e['pse'])}
   if key=='overview':
    err=max(float(abs(colorMap[k]-rgb('#BAC1C7' if k[0]=='A' else '#8B72A5')).max()) for k in colorMap);record.update(stick_atom_count=len(reps),pure_cartoon_no_sticks=len(reps)==0,maximum_component_RGB_error=err,all_144_CaM_and_25_peptide_CA_present=c.count_atoms('bound and chain A and name CA')==144 and c.count_atoms('bound and chain B and name CA')==25);assert len(reps)==0 and err<1e-6
   else:
    n=int(e['target'][1:]);ex=expected[n];selection_complete=set(reps)==ex;target={k for k in reps if k[0]=='B' and k[1]==n};target_err=max(float(abs(colorMap[k]-rgb('#3A8D8C')).max()) for k in target);neighbor_err=max(float(abs(colorMap[k]-rgb('#B9C0C5')).max()) for k in reps if k not in target);listres={(r['chain'],r['native_residue']) for r in e['whole_residues']};record.update(cutoff_A=6.0,native_target=n,sticks_represented_atom_count=len(reps),expected_whole_residue_heavy_atom_count=len(ex),all_expected_whole_residue_atoms_in_sticks_representation=selection_complete,no_sticks_outside_selected_whole_residues=set(reps)<=ex,manifest_whole_residue_list_matches={k[:2] for k in ex}==listres,n_CaM_residues=len({k[:2] for k in ex if k[0]=='A'}),n_peptide_residues=len({k[:2] for k in ex if k[0]=='B'}),target_heavy_atom_count=len(target),maximum_target_RGB_error=target_err,maximum_neighbor_RGB_error=neighbor_err,no_element_colors=target_err<1e-6 and neighbor_err<1e-6,all_selected_atom_projections_inside_frame=e['all_displayed_scope_atom_projections_inside_frame'],no_target_or_neighbor_atom_projection_clipped=e['all_displayed_scope_atom_projections_inside_frame']);assert selection_complete and target_err<1e-6 and neighbor_err<1e-6 and record['manifest_whole_residue_list_matches']
   assert geometry and record['surface_atom_count']==0 and record['same_saved_camera'];A['assets'][key]=record
 v1=np.array(m['assets']['zoom_N297']['view']);v2=np.array(m['assets']['zoom_N308']['view']);A['both_local_views_same_rotation']=bool(np.array_equal(v1[:9],v2[:9]));A['both_local_views_same_optical_distance']=bool(v1[11]==v2[11]);A['both_local_views_same_dimensions']=(m['assets']['zoom_N297']['width'],m['assets']['zoom_N297']['height'])==(m['assets']['zoom_N308']['width'],m['assets']['zoom_N308']['height']);assert A['both_local_views_same_rotation'] and A['both_local_views_same_optical_distance'] and A['both_local_views_same_dimensions']
 (OUT/'independent_geometry_representation_color_audit.json').write_text(json.dumps(A,ensure_ascii=False,indent=2),encoding='utf8');m['independent_audit']='independent_geometry_representation_color_audit.json';m['audit_status']='passed: exact loaded coordinates, actual sticks representation completeness, colors, common local optics and frame retention';(OUT/'render_manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps(A,ensure_ascii=False,indent=2))
if __name__=='__main__':check()