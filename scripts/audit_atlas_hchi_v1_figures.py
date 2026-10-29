"""Check final page size, embedded fonts, editable SVG text and PNG dimensions."""
import json
from pathlib import Path

import fitz
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/flexppi_atlas_hchi_v1_extended_figures"


def run():
    source = json.loads((OUT / "statistics_and_provenance.json").read_text())
    results = []
    for figure in source["figures"]:
        stem = figure["stem"]
        doc = fitz.open(OUT / (stem + ".pdf"))
        if len(doc) != 1:
            raise ValueError("Unexpected multipage figure")
        page = doc[0]
        width, height = page.rect.width * 25.4 / 72, page.rect.height * 25.4 / 72
        spans = [span for block in page.get_text("dict")["blocks"] if "lines" in block
                 for line in block["lines"] for span in line["spans"] if span["text"].strip()]
        clipped = [s["text"] for s in spans if s["bbox"][0] < -.1 or s["bbox"][1] < -.1
                   or s["bbox"][2] > page.rect.width + .1 or s["bbox"][3] > page.rect.height + .1]
        font_checks = []
        for font in page.get_fonts(full=True):
            name, ext, kind, payload = doc.extract_font(font[0])
            font_checks.append(dict(name=name, extension=ext, embedded_bytes=len(payload)))
        with Image.open(OUT / (stem + ".png")) as png:
            size, dpi = png.size, png.info.get("dpi", (0, 0))
        svg = (OUT / (stem + ".svg")).read_text(encoding="utf-8")
        result = dict(stem=stem, width_mm=width, height_mm=height,
                      minimum_pdf_font_pt=min(s["size"] for s in spans),
                      clipped_text=clipped, fonts=font_checks, png_size=size, png_dpi=dpi,
                      svg_editable_text="<text" in svg, svg_raster_image_count=svg.count("<image"))
        result["passed"] = (abs(width - 183) < .05 and height <= 247 and not clipped
                            and result["minimum_pdf_font_pt"] >= 7.9
                            and all(f["embedded_bytes"] > 0 for f in font_checks)
                            and min(dpi) >= 599 and result["svg_editable_text"]
                            and result["svg_raster_image_count"] == 0)
        results.append(result)
        doc.close()
    report = dict(passed=all(r["passed"] for r in results), figures=results,
                  note="This is mechanical PDF/SVG/PNG QA. Visual layout still requires image inspection.")
    (OUT / "render_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
