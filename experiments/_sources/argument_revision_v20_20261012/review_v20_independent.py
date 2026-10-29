from pathlib import Path
import hashlib
import json
import re
import zipfile
from datetime import datetime
from lxml import etree

ROOT = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
OUT = ROOT / 'experiment_completion' / 'argument_revision_v20_20261012'
DOCX = ROOT / 'paper_v20_20261012_中文修订_联合模型_三线表.docx'
PRIOR = ROOT / 'paper_v19_20261012_中文修订_补充实验完成_三线表.docx'
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}


def read_docx(path):
    with zipfile.ZipFile(path) as z:
        tree = etree.fromstring(z.read('word/document.xml'))
        paragraphs = tree.findall('.//w:body/w:p', NS)
        texts = [''.join(p.xpath('.//w:t/text()', namespaces=NS)) for p in paragraphs]
        reltree = etree.fromstring(z.read('word/_rels/document.xml.rels'))
        rels = {r.get('Id'): r.get('Target') for r in reltree}
        figures = []
        for index, paragraph in enumerate(paragraphs):
            embeds = paragraph.xpath('.//a:blip/@r:embed', namespaces=NS)
            for embed in embeds:
                target = rels[embed]
                archive_path = 'word/' + target if not target.startswith('/') else target[1:]
                caption = None
                caption_index = None
                for nxt in range(index + 1, min(index + 5, len(texts))):
                    match = re.match(r'^图\s*(S?\d+)\s*[:：]', texts[nxt])
                    if match:
                        caption = match.group(1)
                        caption_index = nxt
                        break
                figures.append({'id': caption, 'paragraph': index,
                                'caption_paragraph': caption_index,
                                'image_sha256': hashlib.sha256(z.read(archive_path)).hexdigest()})
        tables = [[[''.join(cell.xpath('.//w:t/text()', namespaces=NS))
                    for cell in row.findall('w:tc', NS)]
                   for row in table.findall('w:tr', NS)]
                  for table in tree.findall('.//w:body/w:tbl', NS)]
        comments_hash = (hashlib.sha256(z.read('word/comments.xml')).hexdigest()
                         if 'word/comments.xml' in z.namelist() else None)
        fields = tree.xpath('//w:instrText/text()', namespaces=NS)
    return texts, tables, figures, comments_hash, fields


def citation_numbers(text):
    values = []
    for bracket in re.findall(r'\[([0-9][0-9,，\s–—\-]*)\]', text):
        for token in re.split(r'[,，]', bracket):
            token = token.strip()
            match = re.fullmatch(r'(\d+)\s*[–—\-]\s*(\d+)', token)
            if match:
                values.extend(range(int(match.group(1)), int(match.group(2)) + 1))
            elif token.isdigit():
                values.append(int(token))
    return values


def crossreferences(text):
    hits = []
    for match in re.finditer(r'([图表])\s*(S?\d+)(?:\s*[–—\-]\s*(S?\d+))?', text):
        kind, start, end = match.groups()
        if end:
            prefix = 'S' if start.startswith('S') else ''
            suffix = 'S' if end.startswith('S') else ''
            if prefix == suffix:
                hits.extend((kind, prefix + str(n))
                            for n in range(int(start.lstrip('S')), int(end.lstrip('S')) + 1))
            else:
                hits.extend(((kind, start), (kind, end)))
        else:
            hits.append((kind, start))
    return hits


sha_before = hashlib.sha256(DOCX.read_bytes()).hexdigest()
texts, tables, figures, comments_hash, fields = read_docx(DOCX)
prior_texts, prior_tables, prior_figures, prior_comments_hash, prior_fields = read_docx(PRIOR)
reference_paragraphs = [(i, t) for i, t in enumerate(texts) if re.match(r'^\[\d+\]\s', t)]
reference_ids = [int(re.match(r'^\[(\d+)\]', t).group(1)) for _, t in reference_paragraphs]
body_text = '\n'.join(t for t in texts if not re.match(r'^\[\d+\]\s', t))
cited = sorted(set(citation_numbers(body_text)))
captions = []
for i, text in enumerate(texts):
    match = re.match(r'^([图表])\s*(S?\d+)\s*[:：]', text)
    if match:
        captions.append({'kind': match.group(1), 'id': match.group(2), 'paragraph': i})
caption_set = {(c['kind'], c['id']) for c in captions}
missing_crossrefs = sorted(set(crossreferences(body_text)) - caption_set)
duplicate_captions = sorted({key for key in caption_set
                            if sum((c['kind'], c['id']) == key for c in captions) > 1})

new_to_old = {'1': '1', '2': '2', '3': '3', '4': '4', '5': '5',
              'S1': 'S2', 'S2': 'S3', 'S3': 'S4', 'S4': 'S5',
              'S5': 'S6', 'S6': 'S7', 'S7': 'S8', 'S8': 'S9',
              'S9': 'S12', 'S10': 'S13', 'S11': 'S14'}
old_images = {f['id']: f for f in prior_figures}
mapping = []
for figure in figures:
    old_id = new_to_old.get(figure['id'])
    prior_image = old_images.get(old_id)
    mapping.append({'new_id': figure['id'], 'old_id': old_id,
                    'image_sha256': figure['image_sha256'],
                    'same_image': bool(prior_image and prior_image['image_sha256'] == figure['image_sha256'])})

checks = {
    '31_references_present_and_all_cited': reference_ids == list(range(1, 32)) and cited == list(range(1, 32)),
    'no_missing_crossreferences': not missing_crossrefs,
    'no_duplicate_captions': not duplicate_captions,
    'figures_main_1_5_supplement_1_11': {c['id'] for c in captions if c['kind'] == '图'} == {str(n) for n in range(1, 6)} | {'S' + str(n) for n in range(1, 12)},
    'tables_main_1_2_supplement_1_11': {c['id'] for c in captions if c['kind'] == '表'} == {'1', '2'} | {'S' + str(n) for n in range(1, 12)},
    'image_mapping_matches_v19_after_deletions': len(mapping) == 16 and all(m['same_image'] for m in mapping),
    'comments_xml_unchanged_from_v19': comments_hash == prior_comments_hash,
    'system_and_component_family_protocols_retained': all(s in body_text for s in ['五折体系留出', '全部组成链', '1,662个连通家族']),
    'itrace_family_fold0_training_sources_retained': all(s in body_text for s in ['全部244条ITrace位于第0折测试集合', '该折训练集合来自ATLAS和DynaRepo']),
    'public_weights_different_training_corpora_retained': all(s in body_text for s in ['原发布权重，训练语料不同', '发布模型采用原权重']),
    'planet_native_units_dtype_temperature_structure_mode_retained': all(s in body_text for s in ['encoded结构特征', '300 K与float32', 'nm乘以10转换为Å']),
    'rmsfnet_simulated_density_relative_units_retained': all(s in body_text for s in ['4 Å分辨率、1.5 Å网格', '高斯模拟密度', '相对柔性分数', '相对分数的Å误差不适用']),
    'default_joint_method_and_small_mae_comparison_retained': all(s in body_text for s in ['联合训练作为默认方案', '0.577 Å，单目标模型为0.582 Å', '分别为0.166和0.169 nats']),
    'family_full_intra_comparison_kept_without_drop_warning': '分别降低7.5%和3.8%' in body_text and not any(s in body_text for s in ['家族划分下降', '家族评价下降', '泄漏', '强协同', '显著协同']),
    'prody_versions_use_refs_3_30': '使用ProDy 2.6.1[3,30]' in body_text,
    'no_old_s12_s14_crossreferences': not re.search(r'[图表]\s*S1[2-4](?!\d)', body_text),
    'no_unresolved_ref_fields': not any(re.search(r'\b(?:REF|PAGEREF)\b', f) for f in fields),
}

text_file = OUT / 'v20_content.txt'
text_export = text_file.read_text(encoding='utf-8')
txt_paragraphs = {}
for block in re.split(r'(?m)^P(?=\d+:)', text_export):
    match = re.match(r'(\d+):\s?(.*)', block, flags=re.S)
    if match:
        content = match.group(2)
        content = re.split(r'(?m)^T\d+', content)[0].rstrip('\r\n')
        txt_paragraphs[int(match.group(1))] = content
# Extracted text includes Office Math; compare prose paragraphs without formulas.
with zipfile.ZipFile(DOCX) as z:
    current_tree = etree.fromstring(z.read('word/document.xml'))
    current_paragraphs = current_tree.findall('.//w:body/w:p', NS)
    math_paragraph_indices = [index for index, p in enumerate(current_paragraphs)
                              if p.find('.//{http://schemas.openxmlformats.org/officeDocument/2006/math}oMath') is not None]
plain_export_mismatches = []
for index, text in enumerate(texts):
    if index in txt_paragraphs and text and index not in math_paragraph_indices:
        if re.sub(r'\s+', '', text) not in re.sub(r'\s+', '', txt_paragraphs[index]):
            plain_export_mismatches.append(index)
checks['plain_paragraph_export_matches_docx'] = not plain_export_mismatches

sha_after = hashlib.sha256(DOCX.read_bytes()).hexdigest()
checks['docx_unchanged_during_review'] = sha_before == sha_after
issues = []
for name, passed in checks.items():
    if not passed:
        issues.append({'check': name, 'severity': 'requires_inspection'})
report = {
    'review': 'Independent read-only v20 content, references, crossreferences, methods and scientific argument review',
    'timestamp_local': datetime.now().isoformat(timespec='seconds'),
    'docx': str(DOCX), 'docx_sha256': sha_after,
    'text': str(text_file), 'text_sha256': hashlib.sha256(text_file.read_bytes()).hexdigest(),
    'status': 'pass' if not issues else 'inspect', 'issues': issues,
    'counts': {'paragraphs': len(texts), 'tables': len(tables), 'images': len(figures), 'references': len(reference_ids)},
    'checks': checks, 'reference_ids': reference_ids, 'cited_ids': cited,
    'missing_reference_ids': sorted(set(cited) - set(reference_ids)),
    'uncited_reference_ids': sorted(set(reference_ids) - set(cited)),
    'captions': captions, 'missing_crossreference_captions': missing_crossrefs,
    'duplicate_captions': duplicate_captions, 'image_renumbering': mapping,
    'field_instructions': fields, 'comments_xml_sha256': comments_hash,
    'plain_export_mismatch_candidates': plain_export_mismatches,
    'office_math_paragraph_indices': math_paragraph_indices,
    'numeric_spot_checks': {
        'table1_full': {'RMSF_MAE_Angstrom': 0.577, 'RMSF_rho': 0.796, 'Hchi_MAE_nats': 0.166},
        'table1_planet_full': {'RMSF_MAE_Angstrom': 0.568, 'RMSF_rho': 0.828},
        'table1_planet_seq': {'RMSF_MAE_Angstrom': 0.714, 'RMSF_rho': 0.771},
        'table1_rmsfnet': {'RMSF_MAE_Angstrom': None, 'RMSF_rho': 0.658},
        'tableS9_joint_RMSF': {'MAE': 0.577, 'RMSE': 0.792, 'rho': 0.796},
        'tableS9_single_RMSF': {'MAE': 0.582, 'RMSE': 0.790, 'rho': 0.800},
        'tableS9_joint_Hchi': {'MAE': 0.166, 'RMSE': 0.239, 'rho': 0.941},
        'tableS9_single_Hchi': {'MAE': 0.169, 'RMSE': 0.247, 'rho': 0.938},
        'tableS10_family_full': {'RMSF_MAE': 0.616, 'Hchi_MAE': 0.173},
        'tableS10_family_intra': {'RMSF_MAE': 0.666, 'Hchi_MAE': 0.180},
        'tableS11_standardizedMAE': {'Full': 0.371, 'Intra': 0.399, 'PLANET_full': 0.364,
                                    'PLANET_seq': 0.400, 'RMSFnet': 0.551, 'ANM': 0.450, 'GNM': 0.487},
        'note': 'These displayed values agree with the previously independently audited completed result tables. They are rounded manuscript values; the 7.5% and 3.8% family gains use frozen unrounded means.'
    },
    'manual_findings': [
        '31 references are contiguous and cited; PLANET remains reference 31 and ProDy methods cite 3/30.',
        'Joint model default is justified by dual outputs and small MAE differences, without claiming large synergy or superiority for every RMSF metric.',
        'Family protocol, component-chain grouping, ITrace fold-0-only evaluation, public-weight training difference and native-output units remain explicit.',
        'Family Full/Intra MAE gains 7.5% and 3.8% agree with unchanged supplementary values; no repeated family-performance-decline warning remains.',
        'RMSF-net is described as PDB-derived simulated density with relative scores; no Angstrom MAE is assigned to its output.',
        'KIX retains single 100 ns trajectory per state, rank/localization purpose, opposite state deltas and half-window MLL-minus-cMyb sign change in supplementary tables.',
        'No new unsupported universal superiority, entropy equivalence, independent external-repeat or allosteric-response recovery claim was found.',
        'The six initial exported-text mismatch candidates were paragraphs 15/22/25/27/32/45; each contains Office Math. All non-math paragraph exports agree with DOCX.',
        'This review concerns content and numbering; Word page layout is handled separately by the root agent.'
    ]
}
(OUT / 'v20_independent_content_review.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'status': report['status'], 'docx_sha256': sha_after, 'counts': report['counts'],
                  'failed_checks': [name for name, passed in checks.items() if not passed],
                  'plain_export_mismatch_candidates': plain_export_mismatches,
                  'output': str(OUT / 'v20_independent_content_review.json')}, ensure_ascii=False, indent=2))
