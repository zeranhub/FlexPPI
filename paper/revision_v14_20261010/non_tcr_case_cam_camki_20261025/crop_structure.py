"""Crop the PyMOL render to the square structure panel used in figure v3."""
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
raw = HERE / "work" / "cam_pymol_raw.png"
out = HERE / "cam_structure.png"

im = Image.open(raw).convert("RGB")
arr = np.asarray(im)
mask = np.any(arr < 245, axis=2)
ys, xs = np.where(mask)
margin = 45
box = (max(0, xs.min() - margin), max(0, ys.min() - margin),
       min(arr.shape[1], xs.max() + margin), min(arr.shape[0], ys.max() + margin))
crop = im.crop(box)
w, h = crop.size
side = max(w, h)
canvas = Image.new("RGB", (side, side), "white")
canvas.paste(crop, ((side - w) // 2, (side - h) // 2))
canvas.resize((1600, 1600), Image.LANCZOS).save(out)
print("wrote", out)
