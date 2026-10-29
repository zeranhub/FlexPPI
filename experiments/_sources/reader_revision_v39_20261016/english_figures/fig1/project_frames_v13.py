"""True global and three residue-local projections at the final v13 camera."""
from pathlib import Path
import hashlib,json
import numpy as np
from PIL import Image
import pymol
from pymol import cmd
HERE=Path(__file__).resolve().parent
OUT=HERE/"assets"; SOURCE=HERE/"source_data"
PROBES=SOURCE/"v13_projection_probes"
PROBES.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
PDB=SOURCE/"itrace_4g9f.pdb"
FROZEN=OUT/"whole_chain_receipt.json"
frozen=json.loads(FROZEN.read_text(encoding="utf-8"))
candidate=json.loads((SOURCE/"local_frame_candidates.json").read_text(encoding="utf-8"))
assert sha(PDB)==frozen["source_sha256"]==candidate["source_sha256"]
assert sha(OUT/"whole_chain.png")==frozen["png_sha256"]
assert frozen["residues"]==204 and frozen["camera_roll_deg"]==90
W,H=frozen["raw_image_size"]; crop=frozen["crop"]; view=frozen["view"]
selected=[67,153,166]
pymol.finish_launching(["pymol","-cq"])
cmd.load(str(PDB),"protein");cmd.hide("everything")
cmd.select("Dchain","protein and polymer.protein and chain D")
CA=np.asarray(cmd.get_coords("Dchain and name CA"),dtype="<f8")
original=np.asarray(cmd.get_coords("Dchain"),dtype=np.float64)
assert len(CA)==204
assert hashlib.sha256(CA.tobytes()).hexdigest()==frozen["loaded_CA_xyz_sha256"]
O=CA.mean(axis=0)
frames={}
for number in selected:
 stored=candidate["all_204_local_frames"][number-1]
 assert stored["residue_number"]==number
 resi=stored["residue_identifier"]
 prefix=f"protein and chain D and resi {resi}"
 t=np.array(cmd.get_atom_coords(prefix+" and name CA"),dtype=np.float64)
 n=np.array(cmd.get_atom_coords(prefix+" and name N"),dtype=np.float64)
 c=np.array(cmd.get_atom_coords(prefix+" and name C"),dtype=np.float64)
 x=c-t;x/=np.linalg.norm(x)
 y=n-t;y-=y.dot(x)*x;y/=np.linalg.norm(y)
 z=np.cross(x,y);R=np.column_stack((x,y,z))
 assert np.max(np.abs(R-np.asarray(stored["R_columns_xyz"])))<1e-12
 frames[number]=dict(R=R,t=t,N=n,C=c,identifier=resi,name=stored["residue_name"])
cmd.viewport(W,H);cmd.bg_color("white")
cmd.set("ray_opaque_background",0);cmd.set("orthoscopic",1);cmd.set("antialias",2)
cmd.set("ray_shadows",0);cmd.set("depth_cue",0);cmd.set("ambient",.68)
cmd.set("specular",.05);cmd.set("ray_trace_mode",1)
cmd.set("ray_trace_color",[.53,.65,.72]);cmd.set("ray_trace_gain",.022)
cmd.set_view(view)
probe_records={}
# Saved probe images are valid only for this exact frozen view/crop/source.
signature=dict(source_sha256=sha(PDB),whole_chain_receipt_sha256=sha(FROZEN),raw_image_size=[W,H],crop=crop,view=view)
signature_path=PROBES/"signature.json"
if signature_path.exists():
 assert json.loads(signature_path.read_text(encoding="utf-8"))==signature
else:
 signature_path.write_text(json.dumps(signature,indent=2),encoding="utf-8")
def measure(key,xyz):
 path=PROBES/(key+".png")
 if not path.exists():
  cmd.pseudoatom("single_probe",pos=np.asarray(xyz).tolist(),vdw=.08)
  cmd.show("spheres","single_probe");cmd.color("red","single_probe")
  cmd.set_view(view);cmd.png(str(path),width=W,height=H,dpi=600,ray=1)
  cmd.delete("single_probe")
 bb=Image.open(path).convert("RGBA").getchannel("A").getbbox()
 assert bb is not None,key
 raw=[(bb[0]+bb[2])/2,(bb[1]+bb[3])/2]
 pixel=np.array([raw[0]-crop[0],raw[1]-crop[1]])
 probe_records[key]=dict(world_xyz_A=np.asarray(xyz).tolist(),raw_screen_xy_px=raw,cropped_screen_xy_px=pixel.tolist(),alpha_bbox=list(bb),png_sha256=sha(path))
 return pixel
origin=measure("global_O",O)
world_basis=np.eye(3); length=8.
global_probes={}
for i,axis in enumerate("XYZ"):
 plus=measure("global_"+axis+"_plus8A",O+length*world_basis[:,i])
 minus=measure("global_"+axis+"_minus8A",O-length*world_basis[:,i])
 global_probes[axis]=(plus,minus)
measured=np.column_stack([(global_probes[a][0]-global_probes[a][1])/(2*length) for a in "XYZ"])
flat=np.array(view[:9]).reshape(3,3); choices=[]
for label,M in [("reshape_as_rows",flat),("reshape_then_transpose",flat.T)]:
 nominal=M[:2].copy();nominal[1]*=-1
 scale=float(np.sum(nominal*measured)/np.sum(nominal*nominal))
 residual=float(np.max(np.abs(measured-scale*nominal)))
 choices.append(dict(label=label,camera_matrix=M.tolist(),scale_px_per_A=scale,max_abs_residual_px_per_A=residual))
chosen=min(choices,key=lambda a:a["max_abs_residual_px_per_A"])
assert chosen["label"]=="reshape_then_transpose" and chosen["max_abs_residual_px_per_A"]<.07
M=np.array(chosen["camera_matrix"]);P=M[:2].copy();P[1]*=-1;P*=chosen["scale_px_per_A"]
def axis_record(vector,plus,minus,center):
 delta=(plus-minus)/(2*length)
 return dict(world_unit_vector=np.asarray(vector).tolist(),cropped_origin_xy_px=np.asarray(center).tolist(),plus8A_screen_xy_px=plus.tolist(),minus8A_screen_xy_px=minus.tolist(),screen_delta_xy_px_per_A=delta.tolist(),screen_unit_direction_xy=(delta/np.linalg.norm(delta)).tolist(),projected_length_px_per_A=float(np.linalg.norm(delta)))
global_axes={a:axis_record(world_basis[:,i],*global_probes[a],origin) for i,a in enumerate("XYZ")}
print("Global O and axes probes completed",flush=True)
local_records=[]
pixels={"global_O":origin.tolist()}
for number in selected:
 frame=frames[number];R=frame["R"];t=frame["t"]
 local_origin=measure(f"D{number}_origin",t)
 pixels[f"D{number}_origin"]=local_origin.tolist()
 axes={}
 for i,axis in enumerate("xyz"):
  plus=measure(f"D{number}_{axis}_plus8A",t+length*R[:,i])
  minus=measure(f"D{number}_{axis}_minus8A",t-length*R[:,i])
  row=axis_record(R[:,i],plus,minus,local_origin)
  estimate=P@R[:,i]
  row["nominal_camera_delta_px_per_A"]=estimate.tolist()
  row["max_abs_matrix_residual_px_per_A"]=float(np.max(np.abs(np.array(row["screen_delta_xy_px_per_A"])-estimate)))
  assert row["max_abs_matrix_residual_px_per_A"]<.08
  axes[axis]=row
 estimate=origin+P@(t-O)
 error=float(np.max(np.abs(estimate-local_origin)))
 assert error<1.25,(number,error)
 local_records.append(dict(id=f"D{number}",chain="D",residue_number=number,residue_identifier=frame["identifier"],residue_name=frame["name"],
  N_world_A=frame["N"].tolist(),CA_world_A=t.tolist(),C_world_A=frame["C"].tolist(),
  local_R_columns_xyz=R.tolist(),CA_global_relative_to_O_A=(t-O).tolist(),
  cropped_origin_xy_px=local_origin.tolist(),analytical_origin_xy_px=estimate.tolist(),origin_projection_max_abs_error_px=error,
  axes=axes,local_to_PDB_world="p_world = R_i @ p_local + CA_i",
  local_to_global_centered_at_O="p_global = R_i @ p_local + (CA_i-O)"))
 print("Local probes completed:",number,flush=True)
all_CA=origin+(CA-O)@P.T
all_frames=candidate["all_204_local_frames"]
all_records=[]
for i,xyz in enumerate(CA):
 fr=all_frames[i];R=np.array(fr["R_columns_xyz"])
 assert np.array_equal(xyz,np.array(fr["CA_world_A"]))
 all_records.append(dict(chain="D",residue_number=i+1,residue_identifier=fr["residue_identifier"],residue_name=fr["residue_name"],
  world_CA_xyz_A=xyz.tolist(),global_CA_relative_to_O_A=(xyz-O).tolist(),
  cropped_screen_xy_px=all_CA[i].tolist(),projection_method="analytical orthographic matrix calibrated against actual probes",
  local_R_columns_xyz=R.tolist(),local_basis_screen_delta_px_per_A=(P@R).tolist()))
assert np.array_equal(np.asarray(cmd.get_coords("Dchain"),dtype=np.float64),original)
report=dict(status="real_final_roll90_global_and_three_local_frames_verified",
 source_sha256=sha(PDB),whole_chain_receipt_sha256=sha(FROZEN),whole_chain_png_sha256=sha(OUT/"whole_chain.png"),
 raw_image_size=[W,H],image_size=frozen["image_size"],crop=crop,view=view,camera_roll_deg=90,
 chain="D",residues=204,residue_numbers=list(range(1,205)),loaded_CA_xyz_sha256=frozen["loaded_CA_xyz_sha256"],
 coordinates_modified=False,artificial_bonds_added=False,probe_offset_length_A=length,
 global_origin_world_A=O.tolist(),global_origin_definition="Arithmetic mean of the frozen 204 chain-D CA coordinates",
 global_origin_is_real_3D_coordinate=True,global_O_cropped_screen_xy_px=origin.tolist(),
 global_axes=global_axes,pixels=pixels,
 local_frames=local_records,local_frames_by_id={r["id"]:r for r in local_records},
 local_axes=next(r["axes"] for r in local_records if r["residue_number"]==166),
 local_R_columns_xyz=frames[166]["R"].tolist(),
 selected_local_residue_numbers=selected,
 local_R_convention="R_i=[x_i y_i z_i] has basis vectors as columns, all stored in PDB-world directions",
 global_coordinate_definition="p_global = p_world - O; global XYZ directions equal PDB-world unit vectors",
 local_to_global="p_global = R_i @ p_local + t_i; t_i = CA_i - O",
 global_to_local="p_local = R_i.T @ (p_global - t_i)",
 PyMOL_rotation_convention=chosen,camera_matrix=M.tolist(),
 PyMOL_rotation_other_convention_checks=choices,
 world_to_camera_column_convention="camera_vector = reshape(view[:9],3,3).T @ world_vector",
 camera_nominal_projection_px_per_A=P.tolist(),measured_world_projection_px_per_A=measured.tolist(),
 screen_y_convention="pixels grow downward",orthographic=True,
 all_204_CA_pixels=all_records,all_204_CA_projection_method="Analytical camera projection calibrated by real global +/-8 A probes and checked at selected CA anchors",
 actual_CA_anchor_count=3,maximum_CA_anchor_matrix_error_px=max(r["origin_projection_max_abs_error_px"] for r in local_records),
 probe_records=probe_records,new_training=0,new_inference=0,PyMOL_version=cmd.get_version()[0])
(OUT/"projection_receipt.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
print("V13 projection receipt ready",flush=True)
print("O pixel:",origin.tolist(),"local pixels:",pixels,flush=True)
cmd.quit()

