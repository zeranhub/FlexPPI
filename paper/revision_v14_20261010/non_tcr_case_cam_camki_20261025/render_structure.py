"""Render the 1MXE CaM–CaMKI cartoon used in figure panel a.

The renderer is 3Dmol.js driven by a local copy of the library and captured with
headless Microsoft Edge. No network access and no PyMOL dependency.

Run with the project Python:
    python render_structure.py
"""
from pathlib import Path
import subprocess

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
PDB = Path(r"E:\Project\FlexPPI\experiments\_sources\baselines\cohort\structures\record_aba5e1f3db118b0aba8895de.pdb")
EDGE = Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")
VENDOR = HERE / "3Dmol-min.js"
RAW = WORK / "cam_viewer_raw.png"
HTML = WORK / "cam_viewer.html"
OUT = HERE / "cam_structure.png"

pep_iface = [148, 149, 150, 151, 153, 154, 155, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167]
cam_iface = [7, 10, 11, 14, 15, 31, 32, 35, 47, 50, 51, 68, 71, 75, 80, 81, 83, 84, 87, 116, 119, 120, 123, 140, 141]


def rotated_pdb(angle_deg: float = 45.0, tilt_deg: float = 12.0) -> str:
    lines = PDB.read_text().splitlines()
    coords = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                       for l in lines if l.startswith(("ATOM  ", "HETATM"))])
    center = coords.mean(0)
    t = np.deg2rad(angle_deg)
    ry = np.array([[np.cos(t), 0, np.sin(t)], [0, 1, 0], [-np.sin(t), 0, np.cos(t)]])
    tx = np.deg2rad(tilt_deg)
    rx = np.array([[1, 0, 0], [0, np.cos(tx), -np.sin(tx)], [0, np.sin(tx), np.cos(tx)]])
    rot = rx @ ry
    out = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            p = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) - center
            p = rot @ p + center
            line = line[:30] + f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}" + line[54:]
        out.append(line)
    return "\n".join(out)


def main() -> None:
    assert VENDOR.is_file(), VENDOR
    assert EDGE.is_file(), EDGE
    pdb = rotated_pdb()
    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<script src="file:///{VENDOR.as_posix()}"></script>
<style>html,body{{margin:0;padding:0;background:#fff;}} #v{{width:1800px;height:1400px;}}</style>
</head><body><div id="v"></div>
<script>
var viewer = $3Dmol.createViewer(document.getElementById('v'), {{backgroundColor:'white'}});
viewer.addModel(`{pdb}`, 'pdb');
viewer.setStyle({{chain:'A'}}, {{cartoon:{{color:'#DDE1E5', quality:'high'}}}});
viewer.setStyle({{chain:'B'}}, {{cartoon:{{color:'#E8A05A', quality:'high'}}}});
viewer.setStyle({{chain:'A', resi:{cam_iface}}}, {{cartoon:{{color:'#9FC0DA', quality:'high'}}}});
viewer.setStyle({{chain:'B', resi:{pep_iface}}}, {{cartoon:{{color:'#C4702A', quality:'high'}}}});
viewer.zoomTo();
viewer.zoom(0.96);
viewer.render();
</script></body></html>"""
    HTML.write_text(html, encoding="utf-8")
    subprocess.run([
        str(EDGE), "--headless=new", "--disable-gpu", "--use-gl=swiftshader",
        "--hide-scrollbars", "--virtual-time-budget=3000",
        "--window-size=1800,1400", f"--screenshot={RAW}", "file:///" + HTML.as_posix(),
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    im = Image.open(RAW).convert("RGB")
    arr = np.asarray(im)
    mask = np.any(arr < 245, axis=2)
    ys, xs = np.where(mask)
    margin = 40
    box = (max(0, xs.min() - margin), max(0, ys.min() - margin),
           min(arr.shape[1], xs.max() + margin), min(arr.shape[0], ys.max() + margin))
    crop = im.crop(box)
    w, h = crop.size
    side = max(w, h)
    canvas = Image.new("RGB", (side, side), "white")
    canvas.paste(crop, ((side - w) // 2, (side - h) // 2))
    canvas.resize((1400, 1400), Image.LANCZOS).save(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
