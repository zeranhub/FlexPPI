"""Calibrate PyMOL orthoscopic projection against an archived receipt.

The archived local-inset receipt stores, for a known camera view, both the
3D coordinates of selected atoms and their pixel positions in the rendered
image. If a closed-form projection reproduces those pixels, the same formula
can be used to draw exact 2D axis overlays on new renders.
"""
import json
import numpy as np

RECEIPT = (r"E:\Project\FlexPPI\experiments\_sources"
           r"\framework_v47_20261019\figure1\source_archive\inset_receipt.json")
r = json.load(open(RECEIPT, encoding="utf-8"))
view = np.array(r["view"], dtype=float)
rot = view[:9].reshape(3, 3)
pos = view[9:12]
center = view[12:15]
dist, fov, ortho = view[15], view[16], view[17]
print("rot=\n", np.round(rot, 4))
print("pos", pos, "center", center, "dist", dist, "fov", fov, "ortho", ortho)
print("image", r["image_size"], "crop", r["crop"])

pts3 = r["points"]
px = r["pixels"]
names = list(pts3)

# candidate: screen = k * (rot @ (p - center))[0:2], pixel = origin + [sx, -sy]
for name in names:
    p = np.array(pts3[name]) - np.array(pts3["origin"])
    c = rot @ p
    print(f"{name:7s} 3d->cam {np.round(c,3)}  pixel-delta "
          f"{np.round(np.array(px[name])-np.array(px['origin']),2)}")

# least squares for scale on x and y separately
for k, nm in enumerate(["x", "y"]):
    a, b = [], []
    for name in names:
        p = np.array(pts3[name]) - np.array(pts3["origin"])
        c = rot @ p
        a.append(c[k])
        b.append((np.array(px[name]) - np.array(px['origin']))[k] * (1 if k == 0 else -1))
    a, b = np.array(a), np.array(b)
    scale = float(a @ b / (a @ a)) if a @ a else float('nan')
    print(f"axis {k}: least-squares px/A = {scale:.6f}, residual {np.round(b - scale*a,2)}")
