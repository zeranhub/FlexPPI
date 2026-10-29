"""Re-render a's prediction structure in the a/b shared material.

Every per-residue RGB and RMSF is copied from the original audited v9 source.
The legacy asset name flexibility_envelope was retained in the native diagram,
but the active bytes are the v9 constant-width coloured cartoon, not a putty
or uncertainty envelope. No amplitude radius or geometric rescaling is used.
"""
from pathlib import Path
import os,json,hashlib,csv,shutil
_dll=os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import numpy as np
from PIL import Image
import pymol
from pymol import cmd
OUT=Path(__file__).resolve().parent
S=OUT/'source_archive'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
xyz_sha=lambda a:hashlib.sha256(np.ascontiguousarray(a,dtype='<f8').tobytes()).hexdigest()
consistent=json.loads((OUT/'consistent_a_b_receipt.json').read_text(encoding='utf-8'))
original=json.loads((S/'rmsf_coloring_original_receipt.json').read_text(encoding='utf-8'))
style=consistent['shared_material_settings']
with (S/'all_840_rmsf_colors_original.csv').open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
assert len(rows)==840
with np.load(S/'itrace_4g9f_prediction.npz',allow_pickle=False) as z:
    data={k:z[k].copy() for k in z.files}
assert data['mean'].shape==(840,2) and data['valid'][:,0].sum()==840
assert sha(S/'itrace_4g9f_prediction.npz')==original['source_sha256']['itrace_4g9f_prediction.npz']
legend_src=OUT.parents[3]/'framework_v46_20261019/figure1/source_archive/rmsf_original_scale.json'
# The native 64-cell colourbar stays separate from the continuous per-residue map.
if not legend_src.exists():
    legend_src=Path(r'E:\Project\FlexPPI\experiments\_sources\framework_v46_20261019\figure1\source_archive\rmsf_original_scale.json')
assert legend_src.is_file()
shutil.copy2(legend_src,S/'rmsf_original_64_cell_legend.json')
anchors=np.array([[int(h[k:k+2],16)/255 for k in [1,3,5]] for h in original['palette_anchors_hex']])
assert original['palette_anchors_hex']==['#1D5C93','#57B9BB','#F6D77B','#DD8B50','#B34342']
pymol.finish_launching(['pymol','-cq'])
cmd.load(str(S/'itrace_4g9f.pdb'),'assembly')
xyz=cmd.get_coords('assembly').copy()
assert xyz_sha(xyz)==consistent['source_XYZ_sha256']
lookup={(a.chain,a.resi):a for a in cmd.get_model('assembly and name CA').atom}
assert len(lookup)==840
cmd.dss('assembly')
for k,v in style.items():cmd.set(k,v)
cmd.bg_color('white')
cmd.cartoon('automatic','assembly')
cmd.hide('everything')
cmd.show('cartoon','assembly')
max_rgb_error=0.0
mapped=[]
for i,row in enumerate(rows):
    assert int(row['global_index'])==i
    key=(str(data['chain'][i]),str(data['author_residue_id'][i]))
    assert key==(row['chain'],row['author_residue_id'])
    value=float(data['mean'][i,0])
    assert value==float(row['rmsf_prediction_A'])
    assert str(data['amino_acid'][i])==row['amino_acid']
    atom=lookup[key]
    assert np.array_equal(atom.coord,[float(row[c]) for c in ['ca_x_A','ca_y_A','ca_z_A']])
    rgb=np.array([float(row[c]) for c in ['red','green','blue']])
    k=min(int(value),3)
    expected=anchors[k]*(1-(value-k))+anchors[k+1]*(value-k)
    assert np.array_equal(rgb,expected)
    col='preserved_rmsf_'+str(i)
    cmd.set_color(col,rgb.tolist())
    selection=f'assembly and chain {key[0]} and resi {key[1]}'
    assert cmd.count_atoms(selection+' and name CA')==1
    cmd.color(col,selection)
    stored=np.array(cmd.get_color_tuple(col))
    error=float(np.max(np.abs(rgb-stored)))
    max_rgb_error=max(max_rgb_error,error)
    assert error<3.1e-8
    mapped.append({'global_index':i,'chain':key[0],'author_residue_id':key[1],
                   'rmsf_prediction_A':value,'RGB_original':rgb.tolist(),'RGB_input_to_PyMOL':rgb.tolist(),
                   'RGB_stored_float32':stored.tolist(),'CA_source_XYZ_A':list(atom.coord)})
view=consistent['a']['camera_view']
size=consistent['a']['raw_size']
cmd.viewport(*size)
cmd.set_view(view)
cmd.png(str(OUT/'assembly_prediction_consistent_raw.png'),width=size[0],height=size[1],dpi=600,ray=1)
assert np.array_equal(xyz,cmd.get_coords('assembly'))
assert np.max(np.abs(np.array(cmd.get_view())-np.array(view)))<1e-6
cmd.save(str(OUT/'assembly_prediction_consistent.pse'))
im=Image.open(OUT/'assembly_prediction_consistent_raw.png').convert('RGBA')
crop=consistent['a']['crop']
alpha=im.getchannel('A').getbbox()
assert alpha and alpha[0]>=crop[0] and alpha[1]>=crop[1] and alpha[2]<=crop[2] and alpha[3]<=crop[3]
im.crop(crop).save(OUT/'assembly_prediction_consistent.png',dpi=(600,600))
report={'status':'passed','source_pdb_sha256':sha(S/'itrace_4g9f.pdb'),
        'source_prediction_sha256':sha(S/'itrace_4g9f_prediction.npz'),
        'original_active_png_sha256':original['outputs']['predicted_rmsf_colored.png'],
        'original_active_geometry_definition':original['geometry_definition'],
        'original_active_is_amplitude_envelope_or_putty':False,
        'new_geometry_definition':'Identical unchanged static PDB coordinates with automatic DSS secondary structure and a/b shared constant-thickness cartoon material. Colour alone represents raw predicted RMSF.',
        'all_840_original_prediction_values_exact':True,'all_840_original_input_RGB_exact':True,
        'rendered_residues':840,'valid_RMSF':840,'supervised_Hchi_valid':int(data['valid'][:,1].sum()),
        'chain_lengths':consistent['chain_lengths'],'original_colormap_anchors_hex':original['palette_anchors_hex'],
        'original_colormap_anchor_positions_A':original['palette_anchor_positions_A'],'original_range_A':[0,4],
        'original_64_cell_legend':str(S/'rmsf_original_64_cell_legend.json'),
        'original_RGB_CSV_sha256':sha(S/'all_840_rmsf_colors_original.csv'),
        'maximum_PyMOL_float32_storage_RGB_error':max_rgb_error,
        'all_source_coordinates_sha256':xyz_sha(xyz),'all_after_coordinates_sha256':xyz_sha(cmd.get_coords('assembly')),
        'input_coordinates_unchanged':True,'style_identical_to_a_b':True,'shared_material_settings':style,
        'camera_view':view,'crop':crop,'raw_size':size,'image_size':list(Image.open(OUT/'assembly_prediction_consistent.png').size),
        'png':str(OUT/'assembly_prediction_consistent.png'),'png_sha256':sha(OUT/'assembly_prediction_consistent.png'),
        'pse':str(OUT/'assembly_prediction_consistent.pse'),'pse_sha256':sha(OUT/'assembly_prediction_consistent.pse'),
        'new_inference':0,'new_training':0,'new_MD':0}
(OUT/'consistent_prediction_receipt.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'all_840_preserved_RGB_audit.json').write_text(json.dumps(mapped,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','png':report['png'],'rendered_residues':840,'all_input_RGB_exact':True,'max_float32_storage_error':max_rgb_error}),flush=True)
cmd.quit()
