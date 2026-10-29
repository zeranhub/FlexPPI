from pathlib import Path
import json
import sys
from pypdf import PdfReader
from PIL import Image, ImageOps, ImageDraw

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
folder = Path(sys.argv[1])
pdf = next(folder.glob("*.pdf"))
reader = PdfReader(pdf)
pages = []
for index, page in enumerate(reader.pages, 1):
    text = page.extract_text()
    pages.append({"page": index, "text": text})
    (folder / f"page-{index}.txt").write_text(text, encoding="utf-8")
(folder / "page_text.json").write_text(json.dumps(pages, ensure_ascii=False, indent=2), encoding="utf-8")
for p in pages:
    if any(token in p["text"] for token in ("图2", "图 2", "图S2", "图 S2", "图S3", "图 S3", "Figure 2", "Figure S2", "Figure S3")):
        print(json.dumps({"page": p["page"], "text": p["text"]}, ensure_ascii=False))
print(f"Total PDF pages: {len(pages)}")
for offset in range(0, len(pages), 8):
    batch = pages[offset:offset + 8]
    thumb_w, thumb_h = 510, 720
    sheet = Image.new("RGB", (thumb_w * 4, (thumb_h + 32) * 2), "#eeeeee")
    draw = ImageDraw.Draw(sheet)
    for i, p in enumerate(batch):
        page_image = Image.open(folder / f'page-{p["page"]}.png').convert("RGB")
        page_image.thumbnail((thumb_w, thumb_h))
        x = (i % 4) * thumb_w + (thumb_w - page_image.width) // 2
        y = (i // 4) * (thumb_h + 32) + 30
        sheet.paste(page_image, (x, y))
        draw.text(((i % 4) * thumb_w + 8, (i // 4) * (thumb_h + 32) + 8), f'Page {p["page"]}', fill="black")
    sheet.save(folder / f'contact-{offset + 1}-{offset + len(batch)}.png')
