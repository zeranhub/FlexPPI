"""Verify the adopted DSS ribbon against its PDB, PSE, pixels and original camera."""
from pathlib import Path
import os,json,hashlib
_dll=os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import numpy as np
from PIL import Image
import pymol
from pymol import cmd
ROOT=Path(__file__).resolve().parent
S=ROOT/'source_archive'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
coords_sha=lambda a:hashlib.sha256(np.ascontiguousarray(a,dtype='<f8').tobytes()).hexdigest()
receipt=json.loads((ROOT/'chain_d_ribbon_receipt.json').read_text(encoding='utf-8'))
old=json.loads((S/'whole_chain_receipt.json').read_text(encoding='utf-8'))
projection=json.loads((S/'projection_receipt.json').read_text(encoding='utf-8'))
assert sha(ROOT/'chain_d_ribbon.png')==receipt['image_sha256']
assert sha(S/'itrace_4g9f.pdb')==receipt['source_pdb_sha256']==old['source_sha256']
pymol.finish_launching(['pymol','-cq'])
cmd.load(str(S/'itrace_4g9f.pdb'),'original')
original=cmd.get_coords('original').copy()
original_d=cmd.get_coords('original and chain D').copy()
cmd.delete('all')
cmd.load(str(ROOT/'chain_d_ribbon.pse'))
scene=cmd.get_coords('assembly')
scene_d=cmd.get_coords('assembly and chain D')
assert np.array_equal(original,scene)
assert np.array_equal(original_d,scene_d)
ca=cmd.get_coords('assembly and chain D and name CA')
assert len(ca)==204
assert coords_sha(ca)==old['loaded_CA_xyz_sha256']
model=cmd.get_model('assembly and chain D and name CA')
counts={k:sum(a.ss==k for a in model.atom) for k in ['H','S','L','']}
assert counts==receipt['secondary_structure_counts']=={'H':6,'S':80,'L':117,'':1}
assert np.max(np.abs(np.array(cmd.get_view())-np.array(old['view'])))<1e-6
assert receipt['camera_view']==old['view']
assert cmd.count_atoms('assembly and chain D and rep cartoon')>0
raw=Image.open(ROOT/'chain_d_ribbon_raw.png').convert('RGBA')
cropped=Image.open(ROOT/'chain_d_ribbon.png').convert('RGBA')
assert np.array_equal(np.array(cropped),np.array(raw.crop(receipt['crop'])))
assert cropped.size==tuple(receipt['image_size'])
assert max(cropped.size)>=1000
for key,oldpixel in projection['pixels'].items():
    target=[oldpixel[0]+old['crop'][0]-receipt['crop'][0],oldpixel[1]+old['crop'][1]-receipt['crop'][1]]
    assert receipt['pixels'][key]==target
assert cmd.get_setting_float('cartoon_transparency')==0
report={'status':'passed','reviewer':'fig1_method_audit','adopted_from':'fig1_framework_redraw candidate; independently checked, not rerendered concurrently',
        'runtime':r'D:\Anaconda\envs\pvb_cpu\python.exe','dll_directory':r'D:\Anaconda\envs\pvb_cpu\Library\bin',
        'pymol_version':cmd.get_version()[0],'source_pdb':str(S/'itrace_4g9f.pdb'),'source_pdb_sha256':sha(S/'itrace_4g9f.pdb'),
        'png':str(ROOT/'chain_d_ribbon.png'),'png_sha256':sha(ROOT/'chain_d_ribbon.png'),
        'pse':str(ROOT/'chain_d_ribbon.pse'),'pse_sha256':sha(ROOT/'chain_d_ribbon.pse'),
        'script':str(ROOT/'render_chain_d.py'),'script_sha256':sha(ROOT/'render_chain_d.py'),
        'full_assembly_atom_count':len(original),'chain_D_atom_count':len(original_d),'chain_D_residue_count':204,
        'coordinate_hash_format':'C-contiguous little-endian float64 XYZ bytes in PyMOL atom order',
        'all_source_coordinates_sha256':coords_sha(original),'all_PSE_coordinates_sha256':coords_sha(scene),
        'chain_D_source_coordinates_sha256':coords_sha(original_d),'chain_D_PSE_coordinates_sha256':coords_sha(scene_d),
        'all_input_and_saved_coordinates_identical':True,'input_coordinates_changed':False,
        'secondary_structure_method':'PyMOL dss on actual chain D; automatic regular cartoon; no SS manually assigned',
        'secondary_structure_counts':counts,
        'secondary_structure_per_CA':[{'chain':a.chain,'residue':a.resi,'name':a.resn,'ss':a.ss} for a in model.atom],
        'image_size':list(cropped.size),'image_original_alpha':True,'cartoon_transparency':0.0,
        'camera_matches_original_v13':True,'crop_adjusted_probe_positions_verified':True,
        'original_crop':old['crop'],'new_crop':receipt['crop'],'pixels':receipt['pixels'],
        'visual_review':'PNG individually opened. Real beta arrows, loops, short helix and gold D166 are visible; soft opaque blue material with transparent outer background, no forced tube or faded image alpha.',
        'new_training':0,'new_inference':0,'new_MD':0}
(ROOT/'independent_asset_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:report[k] for k in ['status','png','image_size','secondary_structure_counts','all_input_and_saved_coordinates_identical']}),flush=True)
cmd.quit()
