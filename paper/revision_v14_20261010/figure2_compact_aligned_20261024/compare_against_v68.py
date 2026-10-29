"""Compare final v69 PNG pages with the fully inspected final v68 render."""
from pathlib import Path
import json
import sys
from PIL import Image, ImageChops

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
base = Path(__file__).resolve().parent
old_folder = base.parent / "figure2_title_spacing_20261024" / "render_v68_final"
new_folder = Path(sys.argv[1]) if len(sys.argv) > 1 else base / "render_v69_final"
old_pages = sorted(old_folder.glob("page-*.png"), key=lambda p: int(p.stem.split("-")[-1]))
new_pages = sorted(new_folder.glob("page-*.png"), key=lambda p: int(p.stem.split("-")[-1]))
rows = []
for image_path in new_pages:
    page = int(image_path.stem.split("-")[-1])
    old_path = old_folder / image_path.name
    if old_path.exists():
        old_image = Image.open(old_path).convert("RGB")
        new_image = Image.open(image_path).convert("RGB")
        equal = old_image.size == new_image.size and ImageChops.difference(old_image, new_image).getbbox() is None
    else:
        equal = False
    row = {"v68_page": page if old_path.exists() else None, "v69_page": page, "png_exact_equal": equal}
    rows.append(row)
    print(json.dumps(row))
report = {"v68_page_count": len(old_pages), "v69_page_count": len(new_pages), "pages": rows}
(new_folder / "page_comparison_v68.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
