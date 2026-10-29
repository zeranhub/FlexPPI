"""Record current-page checks and exact inherited page-body render checks."""
from pathlib import Path
import hashlib
import json
from PIL import Image, ImageDraw
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
CURRENT = ROOT / 'render_v72_final_r2'
old = json.loads((ROOT / 'render_v72_final' / 'render_run_comparison.json').read_text())
new = json.loads((CURRENT / 'render_run_comparison.json').read_text())
current_images = sorted(CURRENT.glob('page-*.png'), key=lambda p: int(p.stem.split('-')[-1]))
assert len(current_images) == 32
assert [p['new_page'] for p in new['pages'] if p['body_pixel_match_old_page'] is None] == [10, 12, 18]

old_manual = {9, 10, 11, 12, 13, 15, 16, 17, 19, 20, 33}
reviewed = []
for p in new['pages']:
    page = p['new_page']
    prior = p['body_pixel_match_old_page']
    if prior is None:
        assert page in {10, 12, 18}
        evidence = 'Current full-size PNG visually inspected by root'
    else:
        prior_source = old['pages'][prior - 1]['body_pixel_match_old_page']
        if prior_source is not None:
            evidence = f'Exact pixel body match to prior v72 p{prior}, and verified v71 p{prior_source}'
        else:
            assert prior in old_manual
            evidence = f'Exact pixel body match to prior v72 p{prior}, inspected at full size'
    reviewed.append({'page': page, 'body_review_evidence': evidence})

pdf = PdfReader(next(CURRENT.glob('*.pdf')))
for i, page in enumerate(pdf.pages, 1):
    lines = [s.strip() for s in page.extract_text().splitlines() if s.strip()]
    assert str(i) in lines, (i, lines[:4])
    reviewed[i-1]['page_number_matches_position'] = True

cell_w, cell_h = 230, 105
strip = Image.new('RGB', (8 * cell_w, 4 * cell_h), '#eeeeee')
draw = ImageDraw.Draw(strip)
for i, path in enumerate(current_images):
    im = Image.open(path).convert('RGB')
    crop = im.crop((round(im.width * .44), round(im.height * .94),
                    round(im.width * .56), round(im.height * .985)))
    x = (i % 8) * cell_w + (cell_w - crop.width) // 2
    y = (i // 8) * cell_h + 20
    strip.paste(crop, (x, y))
    draw.text(((i % 8) * cell_w + 8, (i // 8) * cell_h + 4), f'Page {i + 1}', fill='black')
strip.save(CURRENT / 'footer_numbers_1_to_32.png')

target = ROOT.parent / 'paper_v72_20261024_中文修订_结果五节重构.docx'
report = {
    'target_sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
    'page_count': 32,
    'new_page_bodies_inspected': [10, 12, 18],
    'page_bodies_retaining_exact_previous_visual_review': 29,
    'pages': reviewed,
    'repaired_orphan_pages': ['3.3 final paragraph', 'CPU timing paragraph'],
    'footer_strip': str(CURRENT / 'footer_numbers_1_to_32.png'),
    'numeric_repairs_from_unrounded_sources': {'Full_minus_FastProtFlex_AUROC': 0.015, 'single_interface_Full_minus_BackFlip21_rho': 0.069},
    'final_footers_visual_status': 'pending',
}
(CURRENT / 'visual_qa_audit.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'pages': 32, 'target_sha256': report['target_sha256'], 'body_review': 'PASS', 'footers': 'pending'}, ensure_ascii=False))
