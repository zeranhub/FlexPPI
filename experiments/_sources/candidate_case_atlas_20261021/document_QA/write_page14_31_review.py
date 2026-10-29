"""Record the actual individual-page visual review; never modify the DOCX."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone

QA = Path(__file__).parent
binding_path = QA / 'render_binding.json'
binding = json.loads(binding_path.read_text(encoding='utf-8'))
expected_doc = '3285733b085526f8b2942cec1a4043b8a656ffbfc33967b570dba207eaaae2b8'
assert binding['document_sha256'] == expected_doc
render = Path(binding['render_directory'])
case_names = {
    14: 'C06 MAGE-A10: text, native scope table, old/new scope linkage',
    15: 'C06 MAGE-A10: complete figure and caption',
    16: 'C07 AHIII12.2 MHC variants: text and native scope table',
    17: 'C07 AHIII12.2: all three quantitative contrasts, complete figure and caption',
    18: 'C08 1G4 C9V: text and native scope table',
    19: 'C08 1G4 C9V: full peptide RMSF and eight same-AA peptide Hchi sites',
    20: 'C09 NRAS: text and native scope table',
    21: 'C09 NRAS: complete figure and caption',
    22: 'C10 PRAME cross-HLA: text and native scope table',
    23: 'C10 PRAME cross-HLA: complete figure and caption',
    24: 'C11 KRAS G12V maturation: text and native scope table',
    25: 'C11 KRAS G12V maturation: complete figure and caption',
    26: 'C12 PRAME S2 to S2-198: text, native scope table, old/new scope linkage',
    27: 'C12 PRAME S2 to S2-198: complete figure and caption',
    28: 'C13 MEL5 decamer/nonamer: text and native scope table',
    29: 'C13 MEL5 decamer/nonamer: complete figure and caption',
    30: 'C14 DMF5/DMF4: text and table, 16/19 Hchi sites correctly associated',
    31: 'C14 DMF5/DMF4: both receptor contrasts; DMF5 n=16 and DMF4 n=19',
}
issue = {
    'id': 'C08_Hchi_scope_label',
    'page': 18,
    'case_id': 'C08',
    'severity': 'minor_semantic_scope_label',
    'current': '同类型完整界面',
    'required': '同类型完整肽段',
    'reason': 'Figure c and its caption contain all eight same-amino-acid peptide positions P1-P8, not a pMHC interface subset. Numeric values and figure are correct.',
}
pages = []
for n in range(14, 32):
    path = render / f'page-{n}.png'
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    assert sha == binding['PNG_SHA256'][str(n)], n
    pages.append({
        'page': n,
        'path': str(path),
        'PNG_SHA256': sha,
        'status': 'needs_scope_label_fix' if n == 18 else 'passed',
        'matches_render_binding': True,
        'content': case_names[n],
        'individual_original_image_inspected': True,
        'checks': {
            'missing_or_corrupt_characters': 'none_observed',
            'text_or_figure_clipping': 'none_observed',
            'text_overlap': 'none_observed',
            'native_three_line_table': 'passed' if n % 2 == 0 else 'not_present',
            'figure_caption_contained_on_page': 'passed' if n % 2 else 'not_present',
            'paragraph_or_table_pagination': 'passed',
            'unintended_blank_page': False,
            'page_number_visible': True,
        },
        'notes': ('Only the Hchi scope-table wording needs a semantic correction; layout and data agree.' if n == 18
                  else 'Text/table and intended separate figure page are complete; lower whitespace is the case-by-case page arrangement.' if n % 2 == 0
                  else 'English figure labels and complete Chinese caption remain readable and contained on the same page.'),
        'required_changes': [issue] if n == 18 else [],
    })
report = {
    'status': 'needs_scope_label_fix',
    'layout_status': 'passed',
    'reviewer': '/root/case_catalog_20261021',
    'reviewed_at_utc': datetime.now(timezone.utc).isoformat(),
    'document_path': str(QA.parents[1] / 'FlexPPI_32个潜在案例图册_背景证据与局限_20261021.docx'),
    'document_sha256': expected_doc,
    'PDF_SHA256': binding['PDF_SHA256'],
    'render_binding_SHA256': hashlib.sha256(binding_path.read_bytes()).hexdigest(),
    'render_directory': str(render),
    'assigned_page_range': [14, 31],
    'page_count_reviewed': len(pages),
    'method': 'All 18 page PNGs individually opened using view_image detail=original and read at full readable resolution; no montage substituted for page reading.',
    'scope': ['Chinese glyphs', 'clipping and overlap', 'figure/caption containment', 'native three-line tables', 'paragraph/table pagination', 'unintended blank pages', 'scope-table consistency with figure and caption'],
    'additional_scope_check': 'Read C05-C17 figure audits and build scope_rows: C08 is the only peptide-Hchi special case needing this label correction. C14 DMF5=16 and DMF4=19 explicitly verified in text, figure and caption.',
    'pages': pages,
    'required_changes': [issue],
}
# QA is under ATLAS under experiment_completion under BASE.
report['document_path'] = str(QA.parents[2] / 'FlexPPI_32个潜在案例图册_背景证据与局限_20261021.docx')
out = QA / 'page14_31_visual_QA.json'
out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(out)
print(report['status'], len(pages), 'pages')
