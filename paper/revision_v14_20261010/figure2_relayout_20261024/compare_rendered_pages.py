from pathlib import Path
import difflib
import json
import re
import sys
from PIL import Image, ImageChops

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
base = Path(__file__).resolve().parent
old_folder = base / "render_reference"
new_folder = Path(sys.argv[1]) if len(sys.argv) > 1 else base / "render_v65"
old = json.loads((old_folder / "page_text.json").read_text(encoding="utf-8"))
new = json.loads((new_folder / "page_text.json").read_text(encoding="utf-8"))

def norm(text):
    return re.sub(r"\s+", "", re.sub(r"^\d+\s*", "", text))

results = []
for p in old:
    pnorm = norm(p["text"])
    candidates = []
    for q in new:
        qnorm = norm(q["text"])
        score = difflib.SequenceMatcher(None, pnorm, qnorm).ratio()
        candidates.append((score, q))
    score, best = max(candidates, key=lambda row: (row[0], -abs(row[1]["page"] - p["page"])))
    a = Image.open(old_folder / f'page-{p["page"]}.png').convert("RGB")
    b = Image.open(new_folder / f'page-{best["page"]}.png').convert("RGB")
    equal = a.size == b.size and ImageChops.difference(a, b).getbbox() is None
    row = {"reference_page": p["page"], "v65_page": best["page"], "text_similarity": round(score, 6), "png_exact_equal": equal}
    results.append(row)
    print(json.dumps(row))
(new_folder / "page_comparison.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
