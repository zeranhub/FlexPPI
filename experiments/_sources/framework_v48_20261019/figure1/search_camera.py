"""Choose a close-up camera whose local frame reads clearly in projection.

Searches rotations and scores them by (i) how much of each local axis survives
projection, (ii) how far apart the three projected axes are, and (iii) how
compact the backbone trace looks. Writes the chosen angles to a JSON file.
"""
from pathlib import Path
import os, json, itertools
os.add_dll_directory(r"D:\Anaconda\envs\pvb_cpu\Library\bin")
import numpy as np
import pymol
from pymol import cmd

HERE = Path(__file__).resolve().parent
OUT = HERE / "render" / "camera_search.json"

pymol.finish_launching(["pymol", "-cq"])
cmd.load(str(HERE / "assets" / "itrace_4g9f.pdb"), "asm")
cmd.dss("asm")

frames = {}
for resi in (165, 166, 167):
    sel = f"asm and chain D and resi {resi}"
    ca = np.array(cmd.get_atom_coords(f"{sel} and name CA"), float)
    n = np.array(cmd.get_atom_coords(f"{sel} and name N"), float)
    c = np.array(cmd.get_atom_coords(f"{sel} and name C"), float)
    x = (c - ca) / np.linalg.norm(c - ca)
    y = n - ca
    y = y - np.dot(y, x) * x
    y /= np.linalg.norm(y)
    frames[resi] = {"ca": ca, "x": x, "y": y, "z": np.cross(x, y)}

ca_trace = np.array([cmd.get_atom_coords(f"asm and chain D and resi {r} and name CA")
                     for r in range(163, 170)], float)
f = frames[166]


def score(ax_deg, ay_deg, az_deg):
    cmd.reset()
    for axis, ang in (("x", ax_deg), ("y", ay_deg), ("z", az_deg)):
        if ang:
            cmd.turn(axis, ang)
    rot = np.array(cmd.get_view()[:9], float).reshape(3, 3)
    proj = {}
    for key in ("x", "y", "z"):
        v = f[key] / np.linalg.norm(f[key])
        cam = rot.T @ v
        p = np.array([cam[0], -cam[1]])
        proj[key] = p / (np.linalg.norm(p) + 1e-9) * np.linalg.norm(cam[:2])
    lens = {k: float(np.linalg.norm(v)) for k, v in proj.items()}
    angs = []
    for a, b in itertools.combinations(proj, 2):
        u, v = proj[a], proj[b]
        cosang = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-9)
        angs.append(float(np.degrees(np.arccos(np.clip(cosang, -1, 1)))))
    # backbone trace extent in projection
    cams = np.array([rot.T @ (p - f["ca"]) for p in ca_trace])
    span = np.ptp(cams[:, :2], axis=0)
    aspect = span[0] / (span[1] + 1e-9)
    s = min(lens.values()) * 2.0 + (min(angs) - 25) / 45.0
    s -= 0.7 * abs(np.log(aspect))
    return float(s), lens, angs, float(aspect)


best = None
for ax_deg in range(-90, 91, 15):
    for ay_deg in range(-180, 180, 15):
        for az_deg in range(0, 180, 20):
            s, lens, angs, aspect = score(ax_deg, ay_deg, az_deg)
            if best is None or s > best[0]:
                best = (s, (ax_deg, ay_deg, az_deg), lens, angs, aspect)
s, angles, lens, angs, aspect = best
OUT.write_text(json.dumps({"score": s, "turns": angles,
                           "projected_axis_lengths": lens,
                           "min_angle_deg": min(angs), "aspect": aspect}, indent=1))
cmd.quit()
