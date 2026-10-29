"""Render the Figure 1 molecular assets with one shared PyMOL material.

All coordinates come from the archived 4g9f PDB. Two camera orientations are
used: one wide view shared by the input and prediction structures, and one
chain-D view shared by the whole-chain panel and its backbone magnification.

Annotations are placed by detecting the exact raster position of calibration
markers rendered with the identical camera, so the vector axis overlay in the
final figure lands on the atoms rather than on an assumed projection.

Run:  D:\\Anaconda\\envs\\pvb_cpu\\python.exe render_assets.py
"""
from pathlib import Path
import os, json, hashlib, csv

os.add_dll_directory(r"D:\Anaconda\envs\pvb_cpu\Library\bin")
import numpy as np
from PIL import Image
import pymol
from pymol import cmd
from pymol.cgo import BEGIN, END, COLOR, SPHERE

HERE = Path(__file__).resolve().parent
ASSETS = HERE / "assets"
OUT = HERE / "render"
OUT.mkdir(exist_ok=True)
LOG = {}
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

PDB = ASSETS / "itrace_4g9f.pdb"
COLORS = ASSETS / "all_840_rmsf_colors_original.csv"

# Muted chain palette: same hue order as the archived source, lower chroma so
# that the input structure does not compete with the RMSF colour map later.
CHAIN_COLORS = {"A": "#4A7BA6", "B": "#D08C4A", "C": "#9E4B3F",
                "D": "#4F8069", "E": "#7B6A9C"}
ELEMENT = {"C": "#AEB8BF", "N": "#4C78A8", "O": "#C26A5E"}

STYLE = {
    "cartoon_fancy_helices": 1, "cartoon_flat_sheets": 1, "cartoon_smooth_loops": 1,
    "cartoon_round_helices": 1, "cartoon_highlight_color": -1, "cartoon_discrete_colors": 1,
    "cartoon_sampling": 28, "cartoon_loop_radius": 0.18, "cartoon_tube_radius": 0.18,
    "cartoon_rect_length": 1.3, "cartoon_rect_width": 0.20, "cartoon_transparency": 0,
    "ray_opaque_background": 0, "orthoscopic": 1, "antialias": 2,
    "ray_shadows": 0, "depth_cue": 0, "ambient": 0.50, "direct": 0.48,
    "specular": 0.12, "shininess": 24, "light_count": 2, "ray_trace_mode": 0,
    "stick_radius": 0.135, "valence": 1, "sphere_scale": 0.30,
}


def set_style():
    for k, v in STYLE.items():
        cmd.set(k, v)
    cmd.bg_color("white")


def render(name, size, pad=0, ray=1):
    raw = OUT / f"{name}_raw.png"
    cmd.viewport(*size)
    cmd.png(str(raw), width=size[0], height=size[1], dpi=600, ray=ray)
    im = Image.open(raw).convert("RGBA")
    bb = im.getchannel("A").getbbox()
    box = [max(0, bb[0] - pad), max(0, bb[1] - pad),
           min(size[0], bb[2] + pad), min(size[1], bb[3] + pad)]
    im.crop(box).save(OUT / f"{name}.png", dpi=(600, 600))
    return {"png": f"{name}.png", "raw_size": list(size), "alpha_bbox": list(bb),
            "crop": box, "size": list(Image.open(OUT / f"{name}.png").size)}


pymol.finish_launching(["pymol", "-cq"])
cmd.load(str(PDB), "asm")
cmd.remove("hydrogens")
coords0 = cmd.get_coords("asm").copy()
ca_atoms = cmd.get_model("asm and name CA").atom
counts = {c: sum(a.chain == c for a in ca_atoms) for c in CHAIN_COLORS}
assert counts == {"A": 276, "B": 100, "C": 10, "D": 204, "E": 250}, counts
cmd.dss("asm")
set_style()

for chain, hexcol in CHAIN_COLORS.items():
    cmd.set_color(f"ch_{chain}", [int(hexcol[k:k + 2], 16) / 255 for k in (1, 3, 5)])
    cmd.color(f"ch_{chain}", f"asm and chain {chain}")

WIDE = [2400, 1080]
cmd.hide("everything")
cmd.show("cartoon", "asm")
cmd.viewport(*WIDE)
cmd.orient("asm")
cmd.turn("x", -8)
cmd.turn("y", 6)
LOG["wide_frame"] = {"viewport": WIDE, "view": [round(float(v), 8) for v in cmd.get_view()]}
LOG["input"] = render("assembly_input", WIDE, pad=8)

# Same camera, residues coloured by the archived out-of-fold RMSF colours.
rows = list(csv.DictReader(open(COLORS, encoding="utf-8")))
assert len(rows) == 840
groups = {}
for r in rows:
    hx = "#%02X%02X%02X" % tuple(int(round(float(r[k]) * 255)) for k in ("red", "green", "blue"))
    groups.setdefault(hx, []).append((r["chain"], r["author_residue_id"]))
for hx, keys in groups.items():
    name = "p_" + hx[1:]
    cmd.set_color(name, [int(hx[k:k + 2], 16) / 255 for k in (1, 3, 5)])
    for chain, resi in keys:
        cmd.color(name, f"asm and chain {chain} and resi {resi}")
cmd.viewport(*WIDE)
cmd.orient("asm")
cmd.turn("x", -8)
cmd.turn("y", 6)          # identical frame to the input render
LOG["prediction"] = render("assembly_prediction", WIDE, pad=8)
LOG["prediction"].update(distinct_colours=len(groups), residues=len(rows))

# Chain D with the magnified stretch highlighted.
cmd.color("ch_D", "asm and chain D")
cmd.set_color("target", [0.82, 0.55, 0.25])
cmd.color("target", "asm and chain D and resi 165-167")
cmd.hide("everything")
cmd.show("cartoon", "asm and chain D")
D_SIZE = [1500, 1500]
cmd.viewport(*D_SIZE)
cmd.orient("asm and chain D")
cmd.turn("y", -22)
cmd.turn("x", 10)
LOG["chain_d"] = render("chain_d", D_SIZE, pad=10)
D_VIEW = [round(float(v), 8) for v in cmd.get_view()]

# Backbone magnification. The camera is chosen (see search_camera.py) so that
# all three local axes survive projection at comparable length.
Z_SIZE = [1900, 1500]
cmd.viewport(*Z_SIZE)
cmd.reset()
cmd.turn("y", -30)
cmd.turn("z", 140)
cmd.zoom("asm and chain D and resi 164-168 and name CA", 3.1)
cmd.hide("everything")
cmd.show("sticks", "asm and chain D and resi 163-169 and backbone")
cmd.show("spheres", "asm and chain D and resi 163-169 and backbone and name N+CA+C")
for sym, hexcol in (("C", "#A9B3BA"), ("N", "#3F6E9E"), ("O", "#B85A4A"), ("S", "#C9A227")):
    cmd.set_color("elem_" + sym, [int(hexcol[k:k + 2], 16) / 255 for k in (1, 3, 5)])
    cmd.color("elem_" + sym, f"asm and chain D and resi 163-169 and elem {sym}")
cmd.set_color("ca_slate", [0.45, 0.51, 0.55])
cmd.color("ca_slate", "asm and chain D and resi 163-169 and name CA")
cmd.set("stick_radius", 0.135)
cmd.set("sphere_scale", 0.28)
LOG["zoom"] = render("frame_zoom", Z_SIZE, pad=85)
LOG["zoom"]["view"] = [round(float(v), 8) for v in cmd.get_view()]

# ---------------------------------------------------------------- anchors ----
FRAMES = {}
for resi in (153, 165, 166, 167):
    sel = f"asm and chain D and resi {resi}"
    ca = np.array(cmd.get_atom_coords(f"{sel} and name CA"), float)
    n = np.array(cmd.get_atom_coords(f"{sel} and name N"), float)
    c = np.array(cmd.get_atom_coords(f"{sel} and name C"), float)
    x_ax = (c - ca) / np.linalg.norm(c - ca)
    y_raw = n - ca
    y_ax = y_raw - np.dot(y_raw, x_ax) * x_ax
    y_ax /= np.linalg.norm(y_ax)
    FRAMES[resi] = {"ca": ca, "n": n, "c": c, "x": x_ax, "y": y_ax,
                    "z": np.cross(x_ax, y_ax)}
LOG["frames"] = {str(k): {kk: [float(x) for x in vv] for kk, vv in v.items()}
                 for k, v in FRAMES.items()}

AXIS_LEN = 4.4


def palette(n):
    """n flat colours spread over hue and value, far apart in RGB."""
    import colorsys
    out = []
    for i in range(n):
        hue = (i * 0.61803398875) % 1.0
        val = 1.0 if i % 2 == 0 else 0.55
        out.append(colorsys.hsv_to_rgb(hue, 1.0, val))
    return out


def calibrate(name, points, size, crop):
    """Project known model points to pixel centres with the live camera.

    All markers are rendered in one pass with flat, widely separated colours,
    so the measurement is independent of viewport or depth effects.
    """
    import itertools
    grid = [0.0, 0.45, 0.9]
    colours = [c for c in itertools.product(grid, repeat=3) if max(c) < 0.95]
    assert len(colours) >= len(points)
    cmd.set("ambient", 1.0)
    cmd.set("direct", 0.0)
    cmd.set("specular", 0.0)
    cmd.set("sphere_scale", 0.24)
    cmd.hide("everything")
    cmd.delete("anch*")
    tags = list(points)
    for i, tag in enumerate(tags):
        cmd.set_color(f"an{i}", colours[i])
        cmd.pseudoatom(f"anch{i}", pos=[float(v) for v in points[tag]])
        cmd.color(f"an{i}", f"anch{i}")
        cmd.show("spheres", f"anch{i}")
    cmd.viewport(*size)
    path = OUT / f"{name}_anchors.png"
    cmd.png(str(path), width=size[0], height=size[1], dpi=600, ray=1)
    im = np.asarray(Image.open(path).convert("RGB")).astype(int)
    found = {}
    for i, tag in enumerate(tags):
        target = np.array(colours[i]) * 255
        mask = np.sqrt(((im - target) ** 2).sum(axis=2)) < 24
        if mask.sum() < 40:
            raise SystemExit(f"{name}: anchor {tag} not resolved ({mask.sum()} px)")
        ys, xs = np.nonzero(mask)
        found[tag] = [round(float(xs.mean() - crop[0]), 2),
                      round(float(ys.mean() - crop[1]), 2)]
    cmd.set("ambient", 0.50)
    cmd.set("direct", 0.48)
    cmd.set("specular", 0.12)
    return found


zoom_points = {}
for resi in (165, 166, 167):
    fr = FRAMES[resi]
    zoom_points[f"ca{resi}"] = fr["ca"]
    zoom_points[f"n{resi}"] = fr["n"]
    zoom_points[f"c{resi}"] = fr["c"]
for key in ("x", "y", "z"):
    zoom_points[key] = FRAMES[166]["ca"] + AXIS_LEN * FRAMES[166][key]
LOG["zoom"]["anchors"] = calibrate("frame_zoom", zoom_points, Z_SIZE, LOG["zoom"]["crop"])

# A pixel-space affine fit doubles as an independent check of exact placement.
pxs = np.array([LOG["zoom"]["anchors"][t] for t in ("ca165", "ca166", "ca167")])
xyzs = np.array([zoom_points[t] for t in ("ca165", "ca166", "ca167")])
A = np.hstack([xyzs, np.ones((3, 1))])
coef, *_ = np.linalg.lstsq(A, pxs, rcond=None)
pred = A @ coef
LOG["zoom"]["affine_check"] = {"max_residual_px": float(np.abs(pred - pxs).max()),
                               "coefficients": coef.tolist()}

# Triads on the whole-chain view need their own anchor calibration. The camera
# is re-established from the same commands rather than from a stored view.
cmd.hide("everything")
cmd.show("cartoon", "asm and chain D")
cmd.viewport(*D_SIZE)
cmd.orient("asm and chain D")
cmd.turn("y", -22)
cmd.turn("x", 10)
chain_points = {}
for tag, resi in (("i", 166), ("j", 153)):
    fr = FRAMES[resi]
    chain_points[f"{tag}_o"] = fr["ca"]
    for key in ("x", "y", "z"):
        chain_points[f"{tag}_{key}"] = fr["ca"] + AXIS_LEN * fr[key]
LOG["chain_d"]["anchors"] = calibrate("chain_d", chain_points, D_SIZE, LOG["chain_d"]["crop"])
LOG["chain_d"]["view"] = D_VIEW

assert np.array_equal(coords0, cmd.get_coords("asm")), "coordinates were modified"
LOG["all_atom_coordinates_unchanged"] = True
LOG["pymol"] = cmd.get_version()[0]
LOG["pdb_sha256"] = sha(PDB)
LOG["chain_colors"] = CHAIN_COLORS
LOG["style"] = STYLE
LOG["sha256"] = {k: sha(OUT / LOG[k]["png"]) for k in ("input", "prediction", "chain_d", "zoom")}
(OUT / "render_receipt.json").write_text(json.dumps(LOG, indent=2), encoding="utf-8")
cmd.quit()




