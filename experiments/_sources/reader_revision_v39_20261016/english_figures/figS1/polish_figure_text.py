"""Revise scientific SVG labels; all non-text drawing elements stay identical."""
from pathlib import Path
import datetime
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
import fitz

BASE = Path(__file__).resolve().parents[2]
OUT = BASE / 'figures/writing_revision_20261004'
OUT.mkdir(exist_ok=True)
NS = '{http://www.w3.org/2000/svg}'
ET.register_namespace('', 'http://www.w3.org/2000/svg')
ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
sha = lambda b: hashlib.sha256(b).hexdigest()
original_bundle = BASE / 'deliverables/FlexPPI_manuscript_editable_figures_20261004.zip'
if original_bundle.is_file():
    with zipfile.ZipFile(original_bundle) as z:
        manifest = json.loads(z.read('EDITABLE_BUNDLE_MANIFEST.json'))
    svg_paths = [row['path'] for row in manifest['files'] if row['path'].endswith('.svg')]
else:
    # The revised editable bundle includes these original SVG/PDF files.
    saved = json.loads((BASE / 'qa/writing_revision_20261004/figure_text_changes.json').read_text(encoding='utf-8'))
    svg_paths = [row['original_svg'] for row in saved['figures']]
assert len(svg_paths) == 18

replacements = {
    'No ordinary intrachain-five cost was measured. All others output RMSF + Hchi; equal-output forward speed is not compared.': 'Timed models: full 1/5, intrachain 1, static 1/5, and BackFlip2.',
    'All 32 records share one family component; timing repeats are not independent biological replicates. No CI or p-value.': 'Calls include parsing, features, prediction, conversion and compressed output.',
    'Every measured round; no best-run choice': 'Three measured rounds',
    'e,f: round 0 / 1 / 2 = circle / square / triangle. f lines span observed rounds, not confidence intervals.': 'e,f: round 0 / 1 / 2 = circle / square / triangle. Lines span the three observed rounds.',
    'Same CPU2 host; 32 records, one family component. Call timing includes compressed NPZ; no family CI or p-value.': 'Xeon Gold6342 CPU, two threads. Calls include parsing, prediction and compressed output.',
    'Positive change means worse MAE under the strict protocol; negative / zero / positive values retained.': 'Ordinate: strict-model MAE minus ordinary five-model MAE for each record.',
    'Training membership and ensemble size both change; this contrast is not the causal effect of family purge.': 'Ordinary: five-model average. Strict: one model trained after removal of related records.',
    'BF2 repeats the same published-weight prediction, not strict retraining. One family; no CI / p values.': 'BF2 uses the same published weights in both protocols. Large symbols: record-equal means.',
    'Descriptive protocol shift: strict single − ordinary five-fold mean': 'External validation: two training protocols',
    'Prespecified Hχ cases: distinct least-favorable profile': 'Residue profiles selected by Hχ MAE difference',
    'Prespecified least-favorable and lower-median profiles': 'Residue profiles selected by RMSF MAE difference',
    'Rank 0: relative failure': 'Rank 0 by MAE difference',
    'Rank 0: smallest advantage': 'Rank 0 by MAE difference',
    'Rank 15: advantage (lower median)': 'Rank 15 by MAE difference',
    'Negative = relative failure, zero = tie, positive = advantage. Label / model variation is not a family CI.': 'Columns: RMSF (left), Hχ (right). Each row shows the same structure.',
    'Curves keep every position and author ID; chains drawn separately; invalid labels stay gaps, no interpolation.': 'Ticks give chain and author residue IDs. Curves are separated at chain and label gaps.',
    'Cases fixed before scores: ordinary matched_static − full MAE ascending, tie by ID; rank0 / rank15.': 'Selection: static MAE minus full MAE, ascending with ID tie-break; ranks 0 and 15.',
    'Frozen target': 'MD reference',
    'Source-out: ordinary matched five-fold means': 'External validation: ordinary five-model means',
    'All32 records retained for MAE; correlation points use the five / four method joint finite population.': 'MAE and Spearman are calculated within each record on valid residue positions.',
    'Diamonds: record-equal means; thin lines join the same record, not a time course. One family; no CI / p.': 'Diamonds: record-equal means. Thin lines connect methods for the same record.',
    'BF2 is a published-weight reference, only RMSF; its training membership differs from the matched models.': 'Full, intrachain, static and AA models use matched training records; BF2 uses published weights.',
    'Source-out: fixed-fold seeds and strict single models': 'External validation: seeds and strict single models',
    'Seeds use all32, same fold0; dashed / dotted horizontal lines are single member0 static / AA references.': 'Seeds use the same fold0 records. Dashed / dotted lines: single-model static / AA means.',
    'Strict neural and controls each use one member / 2431 training records; BF2 is not a strictly refitted model.': 'Strict models use 2431 training records. Large circles / diamonds: record-equal means.',
    'Ordinary-to-strict changes membership and ensemble size together; no isolation-only causal claim, CI or p.': 'Ordinary seed models and strict models are scored separately on the same 32 records.',
    'All records and valid zeros retained. A/G: χ undefined (hatched). SD is not SEM.': 'A/G: χ undefined (hatched). SD describes variation across windows or replicas.',
    'Observed SD and evidence for the per-residue observation count n': 'Observed SD by per-residue observation count n',
    'Side-chain diagnostics: every amino-acid / record group': 'Within-type prediction errors',
    'A/G Hχ remains NaN. Intrachain: solid orange; static: dashed green in c–d. Full support census is in CSV.': 'A/G have no Hχ labels. c,d: intrachain in solid orange; static in dashed green.',
    'χ count affects target scale; curves do not establish intrinsic chemical difficulty. No CI or p values.': 'MAE is averaged within each record and amino-acid type, then equally across records.',
    'Corrected OOF: common labels and held-out members': 'Five-fold prediction of RMSF and side-chain state entropy',
    'MAE retains every record. Correlations use the joint finite intersection of all five / four methods.': 'MAE and correlations use the same valid residue labels for each target.',
    'Support census preserves own-defined and excluded records; no CI or p values.': 'Each record is predicted by its held-out fold model. BF2 provides RMSF only.',
    'Corrected OOF: partner ablation and source sensitivity': 'Assembly geometry and source-specific prediction errors',
    'Each record uses its held-out fold; all 3063 negative / zero / positive ablation differences remain.': 'All 3063 records are predicted by their held-out fold models.',
    'Ablation is descriptive and does not establish causal partner effects. No CI or p values.': 'a,b: paired model errors. c,d: separate averages for ATLAS, DynaRepo and ITrace.',
    'Points are seed7/11/23; no best-seed selection. All error records retained.': 'Points: seed7/11/23 models trained on the same fold0 records.',
    'Correlation uses all six seeds plus two probes jointly defined; original supports remain in CSV.': 'a,b: record-equal means by source. c,d: differences for the same record and seed.',
    'Sample SD / span across three training seeds is descriptive; no CI or p value.': 'Static and AA references use the same fold0 training records.',
    'Fixed fold0: correlation on common support': 'Fixed fold0: within-record Spearman correlation',
    'All612 records and adverse values retained; no ratio, threshold, CI, p value or significance claim.': 'Each point represents one of the 612 held-out records; colours identify data sources.',
    'Span is max−min of three separately scored members; ordinate is mean MAE gain, not seed-ensemble MAE.': 'Span: maximum minus minimum MAE of three models. Gain: static MAE minus their mean MAE.',
    'relative failure; rmsf gain=-1.44 Å': 'RMSF MAE difference=-1.44 Å',
    'lower median; rmsf gain=+0.0713 Å': 'RMSF MAE difference=+0.0713 Å',
    'relative failure; chi_entropy gain=-0.022 nats': 'Hχ MAE difference=-0.022 nats',
    'Case rule: three-full-seed mean MAE versus matched static; rank0 / rank305 by RMSF.': 'Selection: static MAE minus mean full-model MAE; RMSF ranks 0 and 305.',
    'Extra Hχ least-favourable case is preserved when distinct. Invalid labels are gaps; no interpolation.': 'Curves show valid residue positions and are separated at chain or label gaps.',
    'All chains, positions and six seeds remain in source CSV; profiles do not imply calibrated uncertainty.': 'Columns compare full and intrachain models. Line styles identify the three training seeds.',
    'Remove all linked training members': 'Record length and sequence-group size',
    'External source shares one component': 'Sequence groups and external records',
}
reports = []
for relative in svg_paths:
    source = BASE / relative
    original_bytes = source.read_bytes()
    root = ET.fromstring(original_bytes)
    original_elements = list(root.iter())
    non_text = [(e.tag, sorted(e.attrib.items())) for e in original_elements if e.tag not in {NS+'text', NS+'tspan'}]
    changed = []
    for element in root.iter(NS+'text'):
        value = ''.join(element.itertext())
        replacement = replacements.get(value)
        if source.name == 'fig_seed_profiles_chi_entropy_2.svg' and value == 'Case rule: three-full-seed mean MAE versus matched static; rank0 / rank305 by RMSF.':
            replacement = 'Selection: static MAE minus mean full-model MAE; Hχ rank 0.'
        if replacement is not None:
            assert re.findall(r'[+-]?\d+(?:\.\d+)?', value) == re.findall(r'[+-]?\d+(?:\.\d+)?', replacement) or not re.fullmatch(r'[+-]?\d+(?:\.\d+)?', value)
            for child in list(element):
                element.remove(child)
            element.text = replacement
            changed.append({'old': value, 'new': replacement})
    for element in root.iter(NS+'text'):
        value = ''.join(element.itertext())
        assert not re.search(r'relative failure|no CI|no best|adverse values|no isolation|not a family CI|not confidence|do not establish|does not establish|not SEM|not the causal|no family CI', value, re.I), (relative, value)
    current_non_text = [(e.tag, sorted(e.attrib.items())) for e in root.iter() if e.tag not in {NS+'text', NS+'tspan'}]
    assert non_text == current_non_text
    svg = OUT / source.name
    ET.ElementTree(root).write(svg, encoding='utf-8', xml_declaration=True)
    pdf = svg.with_suffix('.pdf')
    png = svg.with_suffix('.png')
    # Preserve the original PDF vector drawing as well. Change only native text
    # objects; transparent redaction leaves all paths and images in place.
    original_pdf = source.with_suffix('.pdf')
    doc = fitz.open(original_pdf)
    page = doc[0]
    original_drawings = page.get_drawings()
    normalise = lambda s: re.sub(r'\s+', ' ', s.replace('−', '-')).strip()
    change_map = {normalise(c['old']): c['new'] for c in changed}
    lines = [line for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines']]
    edits = []
    for line in lines:
        value = ''.join(span['text'] for span in line['spans'])
        key = normalise(value)
        if key not in change_map:
            continue
        replacement = change_map[key]
        first = line['spans'][0]
        rect = fitz.Rect(line['bbox'])
        # Use the source SVG anchor to keep centered titles centered.
        source_elements = [e for e in ET.fromstring(original_bytes).iter(NS+'text') if normalise(''.join(e.itertext())) == key]
        anchor = 'start'
        if source_elements:
            style = source_elements[0].attrib.get('style', '')
            if 'text-anchor: middle' in style:
                anchor = 'middle'
            elif 'text-anchor: end' in style:
                anchor = 'end'
        fontfile = Path('C:/Windows/Fonts/arialbd.ttf' if first['flags'] & 16 else 'C:/Windows/Fonts/arial.ttf')
        font = fitz.Font(fontfile=str(fontfile))
        size = first['size']
        new_width = font.text_length(replacement, fontsize=size)
        x, y = first['origin']
        if anchor == 'middle':
            x = (rect.x0 + rect.x1 - new_width) / 2
        elif anchor == 'end':
            x = rect.x1 - new_width
        assert x >= -0.01 and x + new_width <= page.rect.width + 0.01, (source.name, replacement, x, new_width)
        page.add_redact_annot(rect, fill=False, cross_out=False)
        edits.append({'old': value, 'new': replacement, 'origin': (x, y), 'size': size,
                      'fontfile': str(fontfile), 'fontname': 'polish_bold' if first['flags'] & 16 else 'polish_regular',
                      'color': ((first['color'] >> 16 & 255)/255, (first['color'] >> 8 & 255)/255, (first['color'] & 255)/255)})
    assert len(edits) == len(changed), (source.name, len(edits), len(changed), set(change_map)-{normalise(e['old']) for e in edits})
    page.apply_redactions(images=0, graphics=0, text=0)
    for edit in edits:
        page.insert_font(fontname=edit['fontname'], fontfile=edit['fontfile'])
        page.insert_text(edit['origin'], edit['new'], fontname=edit['fontname'], fontsize=edit['size'], color=edit['color'])
    def drawing_signature(drawings):
        return [{k:v for k,v in d.items() if k not in {'seqno', 'layer'}} for d in drawings]
    assert drawing_signature(original_drawings) == drawing_signature(page.get_drawings()), (source.name, 'PDF drawings changed')
    pdf.write_bytes(doc.tobytes(garbage=4, deflate=True))
    doc.close()
    doc = fitz.open(pdf)
    doc[0].get_pixmap(matrix=fitz.Matrix(300/72, 300/72), alpha=False).save(png)
    width_mm = doc[0].rect.width * 25.4 / 72
    assert len(doc) == 1 and abs(width_mm - 183) < 0.01, (source.name, width_mm)
    doc.close()
    reports.append({'original_svg': relative, 'original_svg_sha256': sha(original_bytes),
                    'new_svg': svg.relative_to(BASE).as_posix(), 'svg_sha256': sha(svg.read_bytes()),
                    'pdf_sha256': sha(pdf.read_bytes()), 'png_sha256': sha(png.read_bytes()),
                    'width_mm': width_mm, 'text_changes': changed,
                    'original_pdf_sha256': sha(original_pdf.read_bytes()), 'pdf_native_text_changes': len(edits),
                    'all_native_pdf_drawing_signatures_exact_equal': True,
                    'non_text_element_count': len(non_text), 'all_non_text_tags_attributes_order_exact_equal': True})
report = {'status': 'figure_source_text_polished_non_text_geometry_equal_awaiting_visual_review',
          'checked_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'figures': reports, 'old_sources_modified': False,
          'new_data_fitting_scoring_or_case_selection': False,
          'renderer': 'PyMuPDF native PDF text replacement with all original vector drawings preserved; original editable SVG text updated'}
(BASE / 'qa/writing_revision_20261004/figure_text_changes.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'figures': len(reports), 'text_changes': sum(len(r['text_changes']) for r in reports), 'all_non_text_geometry_identical': True}))
