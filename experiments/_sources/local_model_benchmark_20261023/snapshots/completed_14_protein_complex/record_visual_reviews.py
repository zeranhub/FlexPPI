"""Record observations made from the two actual exported PNGs on 2026-10-24.

Run only after viewing those current images. This records hashes; it does not
perform visual inspection itself, and must not certify an unseen future render.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
notes = {
    "benchmark_current.png": [
        "Viewed the actual 2161 x 1665 PNG at original resolution after the final render. Two biological rows are Proteins (2203 records, 1500 families) and Protein complexes (860 records, 261 families), with panels a-f and 14 fixed method/configuration rows in each group.",
        "Blue, ochre and teal are visibly retained. The plot is not rendered in grayscale; direct labels, circle/diamond/triangle/square markers and fixed row positions supplement color.",
        "The final MAE axes show 0.0, 0.5 and 1.0 ticks. The previous awkward endpoint label is removed. No clipped labels or overlap among method names, axis text, group headings or the footer/legend was observed.",
        "Points and family-cluster intervals fit within each axis. Protein-group intervals are often smaller than their markers; this is not a claim that intervals are absent. Larger complex-group intervals remain visible.",
        "All six nonphysical-MAE cells for FastProtFlex, ANM and GNM visibly say relative score. No physical error values were fabricated for them.",
        "The footer visibly states 14 completed entries, same 3063 records and mixed training exposure. It is on separate lines above the bottom color legend. This does not claim equal training or completion of deferred heavy models.",
    ],
    "benchmark_distributions_current.png": [
        "Viewed the actual 2161 x 921 PNG at original resolution after the final render. Two panels use the same two biological labels and n=2203/860.",
        "The six displayed configurations are FlexPPI, BackFlip-2.1 assembly, BackFlip-2.1 chains, DynaProt-M, FastProtFlex and ANM. All six names fit in the external two-row legend.",
        "Blue and teal curves and four ochre curves retain named solid/dashed/dotted/dash-dot patterns. ANM has its separate dash-dot-dot pattern; this supplements color even where some protein-group curves overlap.",
        "Titles, axis labels, tick labels, legend and training-exposure footer are visible without observed clipping or overlap. Negative-rank observations remain within the data-driven full x range; no local zoom or nonlinear axis is used.",
        "The distribution shape is descriptive across all records; the visual does not imply that the 244 related TCR-pMHC records are independent sequence families.",
    ],
}
reviews = {name: {"sha256": hashlib.sha256((ROOT / "figures" / name).read_bytes()).hexdigest(),
    "reviewed_utc": datetime.now(timezone.utc).isoformat(), "status": "inspected_current_hash_with_notes",
    "observations": observations} for name, observations in notes.items()}
(ROOT / "visual_reviews.json").write_text(json.dumps(reviews, indent=2), encoding="utf-8")
print(json.dumps({name: item["sha256"] for name, item in reviews.items()}, indent=2))
