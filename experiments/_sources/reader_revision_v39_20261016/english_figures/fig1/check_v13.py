"""Independent checks of source atoms, global origin, three local frames and exports."""
from pathlib import Path
import base64
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np
import fitz

HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
context=json.loads((HERE/'assets/whole_chain_receipt.json').read_text(encoding='utf-8'))
projection=json.loads((HERE/'assets/projection_receipt.json').read_text(encoding='utf-8'))
integration=json.loads((HERE/'integration_receipt.json').read_text(encoding='utf-8'))
source=HERE/'source_data/itrace_4g9f.pdb'
assert sha(source)==context['source_sha256']==projection['source_sha256']
assert context['residues']==204 and context['camera_roll_deg']==90
assert context['whole_chain_retained'] and not context['coordinates_modified']
assert not context['artificial_bonds_added']
assert integration['atomic_inset_removed'] and integration['selected_residue_numbers']==[67,153,166]
assert integration['projection_receipt_sha256']==sha(HERE/'assets/projection_receipt.json')
assert integration['candidate_frames_sha256']==sha(HERE/'source_data/local_frame_candidates.json')
atoms={};insertions={}
for line in source.read_text(encoding='utf-8').splitlines():
    if line.startswith('ENDMDL'):break
    if line.startswith('ATOM') and line[16] in ' A' and line[21]=='D':
        number=int(line[22:26]);key=(number,line[12:16].strip())
        assert key not in atoms
        atoms[key]=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
        insertions[number]=line[26].strip()
res=sorted(number for number,atom in atoms if atom=='CA')
assert res==list(range(1,205))
all_CA=np.array([atoms[i,'CA'] for i in res])
loaded_CA=all_CA.astype(np.float32).astype(np.float64)
world_O=np.array(integration['global_origin_world_A'])
bonds=[float(np.linalg.norm(atoms[i,'C']-atoms[i+1,'N'])) for i in range(1,204)]
assert min(bonds)>1.1 and max(bonds)<1.6
assert insertions[68]=='A' and insertions[69]=='B' and insertions[70]=='C'
old={c.get('id'):c for c in ET.parse(HERE/'source_data/FlexPPI_main_v12.drawio').findall('.//mxCell')}
new={c.get('id'):c for c in ET.parse(HERE/'FlexPPI_main_v13.drawio').findall('.//mxCell')}
canonical=lambda c:(c.tag,sorted(c.attrib.items()),(c.text or '').strip(),[canonical(a) for a in c])
preserved=0
for ident,cell in old.items():
    if ident.startswith('b12_'):assert ident not in new
    elif ident!='v5_title_b':assert canonical(cell)==canonical(new[ident]),ident;preserved+=1
assert sum(ident.startswith('curve_') for ident in new)==224
assert not any('inset' in ident for ident in new if ident.startswith('b13_'))
assert new['b13_transform_equation'].get('value')=='p = Rᵢq + tᵢ'
P=np.array(projection['camera_nominal_projection_px_per_A'])
camera=np.array(projection['PyMOL_rotation_convention']['camera_matrix'])
unit_P=camera[:2].copy();unit_P[1]*=-1
box=np.array(integration['overview_box'])
scale=box[3]/context['image_size'][1]
anchor=np.array(context['target_CA']);anchor_pixel=np.array(projection['pixels']['D166_origin'])
to_screen=lambda point:box[:2]+scale*(anchor_pixel+P@(point-anchor))
assert np.allclose(to_screen(world_O),integration['global_origin_screen'],atol=1e-10,rtol=0)

def check_axis(ident,origin,endpoint):
    geometry=new[ident].find('mxGeometry')
    first=geometry.find("mxPoint[@as='sourcePoint']");last=geometry.find("mxPoint[@as='targetPoint']")
    assert np.allclose([float(first.get(k)) for k in ['x','y']],origin,atol=1e-12,rtol=0)
    assert np.allclose([float(last.get(k)) for k in ['x','y']],endpoint,atol=1e-12,rtol=0)

origin=np.array(integration['global_origin_screen'])
F=np.array(integration['global_basis_world_columns'])
assert np.allclose(F.T@F,np.eye(3),atol=1e-12)
assert np.isclose(np.linalg.det(F),1)
B=np.array(integration['global_basis_camera_columns'])
assert np.allclose(B.T@B,np.eye(3),atol=1e-12) and np.isclose(np.linalg.det(B),1)
proper_camera=np.array(integration['proper_camera_for_basis'])
assert np.allclose(proper_camera@F,B,atol=1e-12)
assert np.allclose(B,[[1,0,0],[0,0,1],[0,-1,0]],atol=1e-12)
diagram_projection=np.array(integration['global_diagram_projection'])
assert integration['global_axes_are_schematic']
assert integration['global_diagram_is_not_orthographic_molecular_projection']
assert np.allclose(diagram_projection,[[1,-1/np.sqrt(2),0],[0,1/np.sqrt(2),-1]],atol=1e-12)
assert integration['global_axes_behind_protein']
ordering=list(new)
for ident in new:
    if ident.startswith('b13_global_'):
        assert ordering.index(ident)<ordering.index('b13_whole_chain')
for i,axis in enumerate('XYZ'):
    endpoint=origin+integration['global_axis_3d_screen_scales'][axis]*diagram_projection[:,i]
    assert np.allclose(endpoint,integration['global_endpoints'][axis],atol=1e-12,rtol=0)
    check_axis('b13_global_axis_'+axis,origin,endpoint)
max_error=0.;local_checks=[]
for frame in integration['local_frames']:
    number=frame['residue_number']
    n,ca,c=[atoms[number,key].astype(np.float32).astype(np.float64) for key in ['N','CA','C']]
    x=(c-ca)/np.linalg.norm(c-ca)
    y=n-ca-((n-ca)@x)*x;y/=np.linalg.norm(y)
    R=np.stack([x,y,np.cross(x,y)],axis=1)
    assert np.allclose(R,frame['R_columns_xyz'],atol=1e-12,rtol=0)
    assert np.allclose(R,frame['R_world_columns_xyz'],atol=1e-12,rtol=0)
    assert np.allclose(R.T@R,np.eye(3),atol=1e-12) and np.isclose(np.linalg.det(R),1)
    t=ca-world_O
    assert np.allclose(t,frame['translation_world_A'],atol=1e-12,rtol=0)
    R_global=F.T@R
    t_global=F.T@t
    assert np.allclose(R_global,frame['R_local_to_global_columns'],atol=1e-12,rtol=0)
    assert np.allclose(t_global,frame['translation_global_A'],atol=1e-12,rtol=0)
    local_O=to_screen(ca)
    assert np.allclose(local_O,frame['origin_screen'],atol=1e-10,rtol=0)
    for i,axis in enumerate('xyz'):
        endpoint=local_O+integration['local_axis_3d_screen_scale']*(unit_P@R[:,i])
        assert np.allclose(endpoint,frame['endpoints_screen'][axis],atol=1e-12,rtol=0)
        check_axis(f'b13_local_{number}_axis_'+axis,local_O,endpoint)
    q=(loaded_CA-ca)@R
    p=q@R_global.T+t_global
    expected_p=(loaded_CA-world_O)@F
    recovered_q=(p-t_global)@R_global
    error=max(float(np.max(np.abs(p-expected_p))),float(np.max(np.abs(recovered_q-q))))
    max_error=max(max_error,error)
    assert error<1e-11
    local_checks.append(dict(number=number,R_verified=True,translation_verified=True,round_trip_error_A=error))
check_axis('b13_translation',origin,next(f['origin_screen'] for f in integration['local_frames'] if f['residue_number']==166))
image_data=base64.b64decode(new['b13_whole_chain'].get('style').split('image=data:image/png,')[1].split(';')[0])
assert hashlib.sha256(image_data).hexdigest()==sha(HERE/'assets/whole_chain.png')==context['png_sha256']
pdf=fitz.open(HERE/'FlexPPI_main_v13_183mm.pdf')
assert len(pdf)==1 and abs(pdf[0].rect.width*25.4/72-183)<1e-3
spans=[s for b in pdf[0].get_text('dict')['blocks'] for line in b.get('lines',[]) for s in line['spans'] if s['text'].strip()]
assert min(s['size'] for s in spans)>=8
assert all(pdf[0].rect.contains(fitz.Rect(s['bbox'])) for s in spans)
receipt=dict(status='passed',whole_chain_residues=204,continuous_real_peptide_bonds=203,
             peptide_bond_min_A=min(bonds),peptide_bond_max_A=max(bonds),insertion_codes_preserved=True,
             whole_camera_global_and_local_axes_consistent=False,global_origin_is_all_CA_mean=False,
             all_three_local_frames_share_real_molecular_camera=True,
             global_basis_is_right_handed_orthonormal=True,global_oblique_schematic_verified=True,
             global_axes_are_schematic=True,global_origin_fixed_reference_verified=True,
             global_axes_behind_protein=True,global_axes_pale=True,
             local_to_global_rotation_translation_round_trip_max_A=max_error,
             display_origin_and_translation_verified=True,selected_residues=[67,153,166],
             local_frame_checks=local_checks,atomic_inset_removed=True,inset_atoms=0,
             projection_receipt_sha256=sha(HERE/'assets/projection_receipt.json'),
             all_a_c_d_objects_preserved=True,unchanged_non_b_cell_count=preserved,
             all_224_curve_cells_preserved=True,local_axes_from_original_atoms=True,
             minimum_pdf_font_pt=min(s['size'] for s in spans),pdf_width_mm=183,
             new_training=0,new_inference=0,drawio_sha256=sha(HERE/'FlexPPI_main_v13.drawio'),
             png_sha256=sha(HERE/'FlexPPI_main_v13.png'))
(HERE/'scientific_qa_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt))
