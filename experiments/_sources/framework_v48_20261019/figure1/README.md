# Figure 1 masters (v48)

`fig1_framework.pdf` / `.svg` / `.png` are the current Figure 1 masters. The PNG is
600 dpi and is the copy embedded in the Word file; the PDF and SVG are vector,
with Arial text kept as text.

Everything is drawn at the final printed size (170 mm x 152 mm); the Word image is
inserted at exactly 170 mm, so the 8 pt body text in the canvas is 8 pt on the page.

## Rebuild

```
D:\Anaconda\envs\pvb_cpu\python.exe render_assets.py   # PyMOL renders + pixel anchors
D:\Anaconda\envs\pvb_cpu\python.exe search_camera.py   # (optional) close-up camera search
D:\Anaconda\envs\pvb_cpu\python.exe build_figure1.py   # vector + raster masters
python build_v48.py                                     # Word replacement (bundled python)
```

`render_assets.py` reads `assets/itrace_4g9f.pdb` and the archived per-residue
RMSF colours, renders four views with one material, and writes
`render/render_receipt.json`. That receipt stores the pixel position of every
annotation anchor (N, Cα, C and the three local axes of residues 165-167 and the
two chain-D frames), obtained by rendering calibration markers with the live
camera. The drawing script uses these positions directly, so the axis overlay in
panel b lands on the atoms instead of assuming a projection.

`figlib.py` holds the millimetre canvas and the flat vector primitives
(rounded boxes, arrows, glyphs, colour bar, image placement).
`build_figure1.py` holds the layout and all scientific labels.

## Scientific inputs

| item | source |
| --- | --- |
| 4g9f coordinates | `assets/itrace_4g9f.pdb` (archived, SHA in receipt) |
| per-residue RMSF colours | `assets/all_840_rmsf_colors_original.csv` (original 840 RGB) |
| RMSF / Hχ profiles | `assets/curve_local_means.csv` (raw and 11-residue means) |
| MD snapshots | `assets/d_snapshot_*.png` (original renders, unchanged) |

No new training, inference or MD. Colours of the prediction structure are the
archived per-residue RGB values; the colour bar is rebuilt from the same table so
bar and structure share one ramp.

## Checks

`figure_audit.json` records the canvas size, the minimum font size (8 pt), the
render hashes and a hash of every file in this delivery folder.
