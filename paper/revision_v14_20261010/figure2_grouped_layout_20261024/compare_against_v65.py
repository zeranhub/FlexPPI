"""Compare a final v67 render with the fully inspected final v65 render."""
from pathlib import Path
import difflib
import json
import re
import sys
from PIL import Image, ImageChops, ImageDraw

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
base = Path(__file__).resolve().parent
old_folder = base.parent / "figure2_relayout_20261024" / "render_v65_final"
new_folder = Path(sys.argv[1]) if len(sys.argv) > 1 else base / "render_v67_final"
old = json.loads((old_folder / "page_text.json").read_text(encoding="utf-8"))
new = json.loads((new_folder / "page_text.json").read_text(encoding="utf-8"))

def norm(text):
    return re.sub(r"\s+", "", re.sub(r"^\d+\s*", "", text))

results = []
for q in new:
    qnorm = norm(q["text"])
    candidates = []
    for p in old:
        score = difflib.SequenceMatcher(None, norm(p["text"]), qnorm).ratio()
        candidates.append((score, p))
    score, best = max(candidates, key=lambda row: (row[0], -abs(row[1]["page"] - q["page"])))
    a = Image.open(old_folder / f'page-{best["page"]}.png').convert("RGB")
    b = Image.open(new_folder / f'page-{q["page"]}.png').convert("RGB")
    equal = a.size == b.size and ImageChops.difference(a, b).getbbox() is None
    # Word's centered page-number footer changes when content moves by one page.
    # Hide only the verified footer-number box, retaining all document body pixels.
    aa, bb = a.copy(), b.copy()
    if aa.size == bb.size:
        footer_box = (round(aa.width * .44), round(aa.height * .94), round(aa.width * .56), round(aa.height * .985))
        ImageDraw.Draw(aa).rectangle(footer_box, fill="white")
        ImageDraw.Draw(bb).rectangle(footer_box, fill="white")
    body_equal = aa.size == bb.size and ImageChops.difference(aa, bb).getbbox() is None
    row = {"v65_page": best["page"], "v67_page": q["page"], "text_similarity": round(score, 6), "png_exact_equal": equal, "png_body_exact_equal_excluding_footer_number": body_equal}
    results.append(row)
    print(json.dumps(row))
(new_folder / "page_comparison_v65.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
