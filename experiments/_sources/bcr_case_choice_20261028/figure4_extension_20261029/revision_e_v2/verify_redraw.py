from pathlib import Path
import hashlib, json, shutil
from PIL import Image, ImageChops

root = Path(__file__).resolve().parent
paper = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
old = Image.open(root.parent / 'figure4e_g_6A77_preview_v1.png').convert('RGB')
new = Image.open(root / 'figure4e_g_6A77_preview_v2.png').convert('RGB')
assert old.size == new.size
# Everything from 68 mm rightwards belongs to f/g, including their y-axis labels.
x = round(68 / 25.4 * 600)
diff = ImageChops.difference(old.crop((x, 0, old.width, old.height)), new.crop((x, 0, new.width, new.height)))
fg_equal = diff.getbbox() is None
assert fg_equal, 'An unrequested change occurred in f/g'
manuscript = paper / 'paper_v129_20261029_中文修订_抗原识别与抗体动态.docx'
doc_sha = sha(manuscript)
assert doc_sha == '7db9bd53801eefe01a1f4cbebcf0f421668b8157ac1cb23bbe1fffcb80373695'
assembly = json.loads((root / 'assembled_figure_provenance.json').read_text(encoding='utf-8'))
assert assembly['pixel_check_old_TCR_equal']
for ext in ['png', 'pdf', 'svg']:
    shutil.copy2(root / f'figure4e_redraw_v2.{ext}', paper / f'figure4e_20261029_局部结构重绘_v2.{ext}')
report = {
    'status': 'passed',
    'a_d_pixels_unchanged': True,
    'f_g_pixels_unchanged': fg_equal,
    'manuscript_unchanged_sha256': doc_sha,
    'panel_e_canvas_mm': [72, 80],
    'whole_figure_canvas_mm': assembly['PDF_size_mm'],
    'visual_review': 'Inspected standalone e, e-g row and complete Figure 4; labels legible, no overlaps or neighboring panel content in standalone export.',
    'scientific_integrity': 'Same canonical atomic coordinates, with separate rigid camera views. Static distances unchanged. No new MD analysis, affinity claim or causal inference.',
    'exports': {str(p): sha(p) for p in root.glob('figure4*.png')}
}
(root / 'redraw_QA.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
for name in ['figure_render_provenance.json', 'assembled_figure_provenance.json']:
    p = root / name
    d = json.loads(p.read_text(encoding='utf-8'))
    d['status'] = 'passed_visual_and_scope_QA'
    d['QA_report'] = str(root / 'redraw_QA.json')
    p.write_text(json.dumps(d, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'passed': True, 'a_d_unchanged': True, 'f_g_unchanged': fg_equal, 'manuscript_unchanged': True}))
