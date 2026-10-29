"""Render unchanged input assembly and chain D with the same material.

Run using D:/Anaconda/envs/pvb_cpu/python.exe. No coordinates or model values
are generated. Colours are inherited from the original Figure 1 source.
"""
from pathlib import Path
import os,json,hashlib,shutil,csv
_dll=os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import numpy as np
from PIL import Image
import pymol
from pymol import cmd

OUT=Path(__file__).resolve().parent
S=OUT/'source_archive'
OLD=Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\figures')
V1=OLD/'reference_rebuild_20261005/flexppi_main'
V13=OLD/'main_revision_v13_20261008'
V9=OLD/'main_revision_v9_20261006'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
xyz_sha=lambda a:hashlib.sha256(np.ascontiguousarray(a,dtype='<f8').tobytes()).hexdigest()
copies={
    V1/'assets/render_source_manifest.json':'assembly_render_source_manifest.json',
    V1/'source_data/itrace_4g9f_prediction.npz':'itrace_4g9f_prediction.npz',
    V1/'source_data/itrace_4g9f_label.npz':'itrace_4g9f_label.npz',
    V9/'rmsf_structure/all_840_rmsf_colors.csv':'all_840_rmsf_colors_original.csv',
    V9/'rmsf_structure/render_receipt.json':'rmsf_coloring_original_receipt.json',
    V9/'render_colored_rmsf.py':'render_colored_rmsf_original.py',
    V13/'assets/whole_chain_receipt.json':'whole_chain_receipt.json',
    V13/'assets/projection_receipt.json':'projection_receipt.json',
    V13/'source_data/itrace_4g9f.pdb':'itrace_4g9f.pdb',
}
S.mkdir(exist_ok=True)
for src,name in copies.items():shutil.copy2(src,S/name)
manifest=json.loads((S/'assembly_render_source_manifest.json').read_text(encoding='utf-8'))
b_original=json.loads((S/'whole_chain_receipt.json').read_text(encoding='utf-8'))
projection=json.loads((S/'projection_receipt.json').read_text(encoding='utf-8'))
assert sha(S/'itrace_4g9f.pdb')==manifest['pdb']['sha256']==b_original['source_sha256']
palette=manifest['chains']
assert palette=={'A':'#2166AC','B':'#F1A340','C':'#B2182B','D':'#1B7837','E':'#762A83'}

STYLE={
    'cartoon_fancy_helices':1,'cartoon_flat_sheets':1,'cartoon_smooth_loops':1,
    'cartoon_highlight_color':-1,'cartoon_discrete_colors':1,'cartoon_sampling':24,
    'cartoon_loop_radius':.18,'cartoon_tube_radius':.18,
    'cartoon_rect_length':1.2,'cartoon_rect_width':.19,'cartoon_transparency':0,
    'ray_opaque_background':0,'orthoscopic':1,'antialias':2,
    'ray_shadows':0,'depth_cue':0,'ambient':.55,'direct':.40,
    'specular':.10,'shininess':20,'ray_trace_mode':0,'ray_trace_gain':0,
}
def style():
    for k,v in STYLE.items():cmd.set(k,v)
    cmd.bg_color('white')
    cmd.cartoon('automatic','assembly')

def render(name,selection,view,size,pad=45):
    cmd.hide('everything')
    cmd.show('cartoon',selection)
    cmd.viewport(*size)
    cmd.set_view(view)
    assert np.max(np.abs(np.array(cmd.get_view())-np.array(view)))<1e-6
    raw=OUT/(name+'_raw.png')
    cmd.png(str(raw),width=size[0],height=size[1],dpi=600,ray=1)
    im=Image.open(raw).convert('RGBA')
    bb=im.getchannel('A').getbbox()
    assert bb
    crop=[max(0,bb[0]-pad),max(0,bb[1]-pad),min(size[0],bb[2]+pad),min(size[1],bb[3]+pad)]
    image=im.crop(crop)
    image.save(OUT/(name+'.png'),dpi=(600,600))
    cmd.save(str(OUT/(name+'.pse')))
    return {'png':str(OUT/(name+'.png')),'png_sha256':sha(OUT/(name+'.png')),
            'pse':str(OUT/(name+'.pse')),'pse_sha256':sha(OUT/(name+'.pse')),
            'raw_size':list(size),'alpha_bbox':list(bb),'crop':crop,'image_size':list(image.size),'camera_view':list(view)}

pymol.finish_launching(['pymol','-cq'])
cmd.load(str(S/'itrace_4g9f.pdb'),'assembly')
xyz=cmd.get_coords('assembly').copy()
ca=cmd.get_model('assembly and name CA').atom
assert len(ca)==840
chain_counts={c:sum(a.chain==c for a in ca) for c in sorted(palette)}
assert chain_counts=={'A':276,'B':100,'C':10,'D':204,'E':250}
cmd.dss('assembly')
style()
for chain,hexcolor in palette.items():
    rgb=[int(hexcolor[k:k+2],16)/255 for k in (1,3,5)]
    cmd.set_color('original_chain_'+chain,rgb)
    cmd.color('original_chain_'+chain,'assembly and chain '+chain)
models=cmd.get_model('assembly and name CA').atom
ss={c:{k:sum(a.chain==c and a.ss==k for a in models) for k in ['H','S','L','']} for c in palette}
all_coords_hash=xyz_sha(xyz)
a=render('assembly_input_consistent','assembly',manifest['view'],[2400,2100])
assert np.array_equal(xyz,cmd.get_coords('assembly'))
cmd.set_color('original_target_gold',[.88,.64,.28])
cmd.color('original_target_gold','assembly and chain D and resi 166')
b=render('chain_d_consistent','assembly and chain D',b_original['view'],[2600,2200])
assert np.array_equal(xyz,cmd.get_coords('assembly'))
b['pixels']={k:[v[0]+b_original['crop'][0]-b['crop'][0],v[1]+b_original['crop'][1]-b['crop'][1]] for k,v in projection['pixels'].items()}
b['camera_matrix']=projection['PyMOL_rotation_convention']['camera_matrix']
b['camera_nominal_projection_px_per_A']=projection['camera_nominal_projection_px_per_A']
b['target_residue_gold']={'chain':'D','residue':'166','rgb':[.88,.64,.28]}
report={'status':'passed','source_pdb':str(S/'itrace_4g9f.pdb'),'source_pdb_sha256':sha(S/'itrace_4g9f.pdb'),
        'pymol_version':cmd.get_version()[0],'runtime':r'D:\Anaconda\envs\pvb_cpu\python.exe',
        'chain_palette_source':str(S/'assembly_render_source_manifest.json'),'chain_palette':palette,
        'chain_lengths':chain_counts,'secondary_structure_assignment':'PyMOL dss from unchanged full assembly coordinates',
        'secondary_structure_per_chain':ss,'style_identical_a_b':True,'shared_material_settings':STYLE,
        'a':a,'b':b,'coordinate_hash_format':'little-endian float64 C-contiguous PyMOL XYZ bytes',
        'source_XYZ_sha256':all_coords_hash,'after_XYZ_sha256':xyz_sha(cmd.get_coords('assembly')),
        'all_atom_coordinates_unchanged':True,'new_inference':0,'new_MD':0,
        'note':'a assembly is chain-coloured; b uses exactly the same original chain-D green and only highlights D166 in the original target gold. DSS and material are shared; cameras retain their respective original Figure 1 views.'}
(OUT/'consistent_a_b_receipt.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
# A review-only side-by-side preview; this does not replace either native asset.
ims=[Image.open(OUT/n).convert('RGBA') for n in ['assembly_input_consistent.png','chain_d_consistent.png']]
sheet=Image.new('RGBA',(2000,850),'white')
for i,im in enumerate(ims):
    im.thumbnail((950,760))
    sheet.alpha_composite(im,(i*1000+(1000-im.width)//2,(850-im.height)//2))
sheet.convert('RGB').save(OUT/'consistent_a_b_preview.png')
print(json.dumps({'status':'passed','a':a['png'],'b':b['png'],'pair_preview':str(OUT/'consistent_a_b_preview.png'),'chain_palette':palette,'dss_chain_D':ss['D']}),flush=True)
cmd.quit()
