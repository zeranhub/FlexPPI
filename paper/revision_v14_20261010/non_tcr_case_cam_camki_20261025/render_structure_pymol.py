"""PyMOL render of the 1MXE CaM–CaMKI complex for figure panel a.

Run with the local PyMOL environment:
    E:\\Project\\FlexPPI\\tmp\\mamba_root\\envs\\pymol\\Scripts\\pymol.bat -cq render_structure_pymol.py

The script writes work/cam_pymol_raw.png; run crop_structure.py afterwards to
produce the square cam_structure.png used by make_case_figure_v3.py.
"""
from pathlib import Path

from pymol import cmd

CASE = Path(r"E:\Project\FlexPPI\paper\revision_v14_20261010\non_tcr_case_cam_camki_20261025")
WORK = CASE / "work"
WORK.mkdir(exist_ok=True)
PDB = r"E:\Project\FlexPPI\experiments\_sources\baselines\cohort\structures\record_aba5e1f3db118b0aba8895de.pdb"

cmd.load(PDB, "cam")
cmd.remove("solvent")
cmd.hide("everything")
cmd.bg_color("white")
cmd.show("cartoon", "cam")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops", 1)
cmd.set("cartoon_transparency", 0.0)
cmd.set_color("cam_gray", [0.80, 0.80, 0.80])
cmd.set_color("cam_orange", [0.85, 0.55, 0.25])
cmd.color("cam_gray", "chain A")
cmd.color("cam_orange", "chain B")
cmd.set("ray_shadows", 1)
cmd.set("ray_opaque_background", 1)
cmd.set("ambient", 0.35)
cmd.set("direct", 0.55)
cmd.set("specular", 0.15)
cmd.set("shininess", 10)
cmd.set("antialias", 2)
cmd.set("ray_trace_mode", 1)
cmd.set("ray_trace_gain", 0.12)
cmd.orient("cam")
cmd.turn("x", -20)
cmd.zoom("cam", 1.0)
cmd.ray(2400, 2400)
cmd.png(str(WORK / "cam_pymol_raw.png"), dpi=300)
print("wrote", WORK / "cam_pymol_raw.png")
