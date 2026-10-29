"""Compare two canonical PNG render directories, including exact shifted matches."""
from pathlib import Path
import hashlib
import json
import sys
from PIL import Image, ImageChops, ImageDraw

old_folder, new_folder = map(Path, sys.argv[1:3])
def page_paths(folder):
    return sorted(folder.glob("page-*.png"), key=lambda p: int(p.stem.split("-")[-1]))
def identity(image, hide_footer=False):
    image = image.convert("RGB")
    if hide_footer:
        image = image.copy()
        ImageDraw.Draw(image).rectangle((round(image.width*.44), round(image.height*.94), round(image.width*.56), round(image.height*.985)), fill="white")
    return hashlib.sha256(str(image.size).encode() + image.tobytes()).hexdigest()
reference = []
for p in page_paths(old_folder):
    im=Image.open(p)
    reference.append({"page":int(p.stem.split("-")[-1]), "whole":identity(im), "body":identity(im,True)})
rows=[]
for p in page_paths(new_folder):
    page=int(p.stem.split("-")[-1]); im=Image.open(p)
    whole, body=identity(im),identity(im,True)
    match=next((r for r in reference if r["whole"]==whole),None)
    body_match=next((r for r in reference if r["body"]==body),None)
    rows.append({"new_page":page,"whole_pixel_match_old_page":match["page"] if match else None,"body_pixel_match_old_page":body_match["page"] if body_match else None})
result={"reference_directory":str(old_folder),"new_directory":str(new_folder),"old_page_count":len(reference),"new_page_count":len(rows),"pages":rows}
(new_folder / "render_run_comparison.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
print(json.dumps(result))
