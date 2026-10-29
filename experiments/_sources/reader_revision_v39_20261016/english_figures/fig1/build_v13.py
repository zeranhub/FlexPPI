"""One protein, a shared global frame, three true residue frames, and one transform."""
from pathlib import Path
import base64
import copy
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'source_data/FlexPPI_main_v12.drawio'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
context=json.loads((HERE/'assets/whole_chain_receipt.json').read_text(encoding='utf-8'))
projection=json.loads((HERE/'assets/projection_receipt.json').read_text(encoding='utf-8'))
frames=json.loads((HERE/'source_data/local_frame_candidates.json').read_text(encoding='utf-8'))
assert context['camera_roll_deg']==90 and context['residues']==204
assert projection['source_sha256']==frames['source_sha256']==context['source_sha256']
assert projection['whole_chain_png_sha256']==context['png_sha256']==sha(HERE/'assets/whole_chain.png')
doc=ET.parse(SOURCE)
root=doc.find('.//root')
before={c.get('id'):ET.tostring(c) for c in root}
removed=[]
for cell in list(root):
    if cell.get('id','').startswith('b12_'):
        removed.append(cell.get('id'));root.remove(cell)
cells={c.get('id'):c for c in root}
cells['v5_title_b'].set('value','b   Global and local frames')

def box(ident,value,bounds,style):
    cell=ET.SubElement(root,'mxCell',id='b13_'+ident,value=value,vertex='1',parent='1',style=style)
    g=ET.SubElement(cell,'mxGeometry',attrib={'as':'geometry'})
    for key,value in zip(('x','y','width','height'),bounds):g.set(key,str(value))
    return cell

def text(ident,value,bounds,color='#263238',size=28,bold=False):
    return box(ident,value,bounds,f'text;html=0;whiteSpace=wrap;fillColor=none;strokeColor=none;fontFamily=Arial;fontSize={size};fontColor={color};fontStyle={int(bold)};align=center;verticalAlign=middle;spacing=0;')

def line(ident,start,end,color,width=2.6,arrow=True,end_size=5):
    cell=ET.SubElement(root,'mxCell',id='b13_'+ident,parent='1',edge='1',style=f'edgeStyle=none;html=0;strokeColor={color};strokeWidth={width};endArrow={"block" if arrow else "none"};endFill=1;endSize={end_size};rounded=0;')
    geom=ET.SubElement(cell,'mxGeometry',relative='1',attrib={'as':'geometry'})
    for name,point in [('sourcePoint',start),('targetPoint',end)]:
        ET.SubElement(geom,'mxPoint',x=str(point[0]),y=str(point[1]),attrib={'as':name})
    return cell

width,height=context['image_size']
scale=min(505/width,350/height)
overview_box=[50+(505-width*scale)/2,795+(350-height*scale)/2,width*scale,height*scale]
protein_style='shape=image;imageAspect=0;aspect=fixed;image=data:image/png,'+base64.b64encode((HERE/'assets/whole_chain.png').read_bytes()).decode()+';'
P=np.array(projection['camera_nominal_projection_px_per_A'],dtype=float)
camera=np.array(projection['PyMOL_rotation_convention']['camera_matrix'],dtype=float)
unit_projection=camera[:2].copy();unit_projection[1]*=-1
# The large triad is a conventional oblique coordinate diagram, matching the
# user's marked reference: X right, Z up, and the depth axis Y left-down.
# Its pictorial projection is kept separate from the actual orthographic
# projection of the molecule and residue-local axes. No protein is distorted.
U,_,Vt=np.linalg.svd(camera)
proper_camera=U@Vt
assert np.linalg.det(proper_camera)>0
B=np.array([[1.,0.,0.],[0.,0.,1.],[0.,-1.,0.]])
F=proper_camera.T@B
assert np.allclose(F.T@F,np.eye(3),atol=1e-12)
assert np.isclose(np.linalg.det(F),1)
anchor_world=np.array(context['target_CA'])
anchor_pixel=np.array(projection['pixels']['D166_origin'])
all_frames={frame['residue_number']:frame for frame in frames['all_204_local_frames']}
all_CA=np.array([all_frames[i]['CA_world_A'] for i in range(1,205)])
def to_screen(point):
    return np.array(overview_box[:2])+scale*(anchor_pixel+P@(np.array(point)-anchor_world))
global_origin=np.array([260.,955.])
reference_world=all_CA.mean(axis=0)
world_origin=reference_world+np.linalg.pinv(P)@((global_origin-to_screen(reference_world))/scale)
assert np.allclose(to_screen(world_origin),global_origin,atol=1e-10)
selected=[67,153,166]
local={}
for number in selected:
    frame=all_frames[number]
    CA=np.array(frame['CA_world_A'])
    R=np.array(frame['R_columns_xyz'])
    local[number]=dict(residue_number=number,residue_name=frame['residue_name'],
                       residue_identifier=frame['residue_identifier'],CA_world_A=CA.tolist(),
                       R_columns_xyz=R.tolist(),R_world_columns_xyz=R.tolist(),
                       R_local_to_global_columns=(F.T@R).tolist(),
                       translation_world_A=(CA-world_origin).tolist(),
                       translation_global_A=(F.T@(CA-world_origin)).tolist(),
                       origin_screen=to_screen(CA).tolist())

colors=['#BB1E36','#23834D','#266CB5']
global_colors=['#D4A1AA','#A2CDB6','#9CBBDD']
global_scale=135
global_scales={'X':125,'Y':120,'Z':135}
diagram_projection=np.array([[1.,-1/np.sqrt(2),0.],[0.,1/np.sqrt(2),-1.]])
global_endpoints={}
for index,axis in enumerate('XYZ'):
    endpoint=global_origin+global_scales[axis]*diagram_projection[:,index]
    global_endpoints[axis]=endpoint.tolist()
    line('global_axis_'+axis,global_origin,endpoint,global_colors[index],width=4.0,end_size=6)
    offset=16*(endpoint-global_origin)/np.linalg.norm(endpoint-global_origin)
    label_box=[endpoint[0]+offset[0]-16,endpoint[1]+offset[1]-18,32,36]
    if axis=='Y':label_box=[endpoint[0]-20,endpoint[1]+30,32,36]
    text('global_label_'+axis,axis,label_box,global_colors[index],32,True)
box('global_origin_dot','',[global_origin[0]-3.8,global_origin[1]-3.8,7.6,7.6],'ellipse;fillColor=#A2B1BB;strokeColor=none;')
text('global_origin_label','O',[global_origin[0]-40,global_origin[1]-35,29,34],color='#92A5B1',size=29)
# The whole global triad is behind the transparent protein render.
box('whole_chain','',overview_box,protein_style)
line('translation',global_origin,local[166]['origin_screen'],'#967342',width=1.8,end_size=4)
mid=(global_origin+np.array(local[166]['origin_screen']))/2
text('translation_label','tᵢ',[mid[0]-42,mid[1]-12,36,33],color='#866223',size=28)

local_scale=38
labels={67:('Gln67',[-34,-61]),153:('Asp153',[-51,32]),166:('Met166',[-46,35])}
for number,frame in local.items():
    origin=np.array(frame['origin_screen']);R=np.array(frame['R_columns_xyz'])
    endpoints={}
    for index,axis in enumerate('xyz'):
        endpoint=origin+local_scale*(unit_projection@R[:,index])
        endpoints[axis]=endpoint.tolist()
        line(f'local_{number}_halo_'+axis,origin,endpoint,'#FFFFFF',width=5,end_size=4)
        line(f'local_{number}_axis_'+axis,origin,endpoint,colors[index],width=2.7,end_size=4)
        if number==166:
            offset=13*(endpoint-origin)/np.linalg.norm(endpoint-origin)
            text(f'local_{number}_label_'+axis,axis,[endpoint[0]+offset[0]-12,endpoint[1]+offset[1]-16,24,32],colors[index],26)
    frame['endpoints_screen']=endpoints
    name,offset=labels[number]
    bounds=[origin[0]+offset[0],origin[1]+offset[1],100,34]
    frame['residue_label_bounds']=bounds
    text(f'local_{number}_residue',name,bounds,color='#5B7180',size=27)
    box(f'local_{number}_origin_dot','',[origin[0]-3.1,origin[1]-3.1,6.2,6.2],'ellipse;fillColor=#D9A547;strokeColor=#FFFFFF;strokeWidth=1;')
text('chain_name','Chain D · 204',[39,1191,270,35],color='#607D8B',size=27)
text('transform_equation','p = Rᵢq + tᵢ',[312,1188,261,39],size=31)

panel_root=ET.Element('root');ET.SubElement(panel_root,'mxCell',id='0');ET.SubElement(panel_root,'mxCell',id='1',parent='0')
for cell in root:
    if cell.get('id','').startswith('b13_'):
        clone=copy.deepcopy(cell);g=clone.find('mxGeometry')
        if g is not None:
            for key,shift in [('x',25),('y',751)]:
                if key in g.attrib:g.set(key,str(float(g.get(key))-shift))
            for point in g.iter('mxPoint'):
                for key,shift in [('x',25),('y',751)]:
                    if key in point.attrib:point.set(key,str(float(point.get(key))-shift))
        panel_root.append(clone)
mx=ET.Element('mxfile',host='app.diagrams.net',version='24.7.17',type='device')
diagram=ET.SubElement(mx,'diagram',id='b13',name='Global and residue frames')
graph=ET.SubElement(diagram,'mxGraphModel',page='1',pageWidth='570',pageHeight='520',grid='0',background='#FFFFFF');graph.append(panel_root)
panel_doc=ET.ElementTree(mx);ET.indent(panel_doc,space='  ')
panel_doc.write(HERE/'FlexPPI_local_frames_v13.drawio',encoding='utf-8',xml_declaration=True)
preserved=[c.get('id') for c in root if c.get('id') in before and ET.tostring(c)==before[c.get('id')]]
assert all(c.get('id') in preserved for c in root if not c.get('id','').startswith('b13_') and c.get('id')!='v5_title_b')
assert sum(key.startswith('curve_') for key in preserved)==224
ET.indent(doc,space='  ');doc.write(HERE/'FlexPPI_main_v13.drawio',encoding='utf-8',xml_declaration=True)
receipt=dict(source_v12_sha256=sha(SOURCE),removed_b_cells=removed,preserved_cells=preserved,
             overview_box=overview_box,displayed_chain='D',complete_residues=204,atomic_inset_removed=True,
             all_other_panels_exactly_preserved=True,projection_receipt_sha256=sha(HERE/'assets/projection_receipt.json'),
             candidate_frames_sha256=sha(HERE/'source_data/local_frame_candidates.json'),
             global_origin_definition='fixed reference point projected at the user-marked left-of-center origin',
             global_origin_world_A=world_origin.tolist(),global_origin_screen=global_origin.tolist(),
             global_basis_world_columns=F.tolist(),global_basis_camera_columns=B.tolist(),
             proper_camera_for_basis=proper_camera.tolist(),
             global_coordinate_definition='p = F.T @ (world_point - O)',
             global_axis_convention='conventional oblique schematic: Z up, X right, Y left-down',
             global_diagram_projection=diagram_projection.tolist(),
             global_axes_are_schematic=True,
             global_diagram_is_not_orthographic_molecular_projection=True,
             global_axes_behind_protein=True,global_axis_colors=global_colors,
             global_axis_3d_screen_scale=global_scale,global_axis_3d_screen_scales=global_scales,
             global_axis_lengths_are_schematic=True,global_endpoints=global_endpoints,
             local_axis_3d_screen_scale=local_scale,selected_residue_numbers=selected,
             local_frames=list(local.values()),local_to_global_equation='p = R_i q + t_i',
             global_points_relative_to_display_origin=True,world_axes_and_local_axes_camera_consistent=False,
             all_three_local_frames_share_real_molecular_camera=True,
             new_training=0,new_inference=0)
(HERE/'integration_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps({key:value for key,value in receipt.items() if key not in ['removed_b_cells','preserved_cells']}))
