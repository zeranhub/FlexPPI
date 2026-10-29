"""Contact sheet of candidate chain-D camera orientations (design aid only)."""
from pathlib import Path
import os, json
os.add_dll_directory(r"D:\Anaconda\envs\pvb_cpu\Library\bin")
import pymol
from pymol import cmd
from PIL import Image

HERE = Path(__file__).resolve().parent
OUT = HERE / "render" / "views"
OUT.mkdir(parents=True, exist_ok=True)

pymol.finish_launching(["pymol", "-cq"])
cmd.load(str(HERE / "assets" / "itrace_4g9f.pdb"), "asm")
cmd.bg_color("white")
cmd.set("ray_opaque_background", 0)
cmd.set("orthoscopic", 1)
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)
cmd.set("ambient", 0.5)
cmd.set("direct", 0.48)
cmd.set("specular", 0.12)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_flat_sheets", 1)
cmd.dss("asm")
cmd.hide("everything")
cmd.show("cartoon", "asm and chain D")
cmd.set_color("dgreen", [0.31, 0.50, 0.41])
cmd.color("dgreen", "asm and chain D")
cmd.set_color("tgt", [0.82, 0.55, 0.25])
cmd.color("tgt", "asm and chain D and resi 165-167")

VIEWS = [(0, 0), (0, 40), (0, -40), (25, 0), (-25, 0), (20, 60), (-20, -60),
         (35, 120), (-35, -120), (60, 0), (0, 90), (0, -90)]
tiles = []
for i, (tx, ty) in enumerate(VIEWS):
    cmd.viewport(700, 700)
    cmd.orient("asm and chain D")
    if ty:
        cmd.turn("y", ty)
    if tx:
        cmd.turn("x", tx)
    cmd.zoom("asm and chain D", 3)
    f = OUT / f"view_{i:02d}.png"
    cmd.png(str(f), width=700, height=700, dpi=150, ray=1)
    tiles.append((i, tx, ty, f))
sheet = Image.new("RGB", (700 * 4, 700 * 3), "white")
for n, (i, tx, ty, f) in enumerate(tiles):
    im = Image.open(f).convert("RGB")
    sheet.paste(im, ((n % 4) * 700, (n // 4) * 700))
sheet.save(OUT / "contact.png")
(OUT / "views.json").write_text(json.dumps(
    [{"index": i, "turn_x": tx, "turn_y": ty} for i, tx, ty, _ in tiles], indent=1), encoding="utf-8")
cmd.quit()
