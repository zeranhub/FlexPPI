from pathlib import Path
import hashlib
import json
import re
import zipfile
from datetime import datetime
from lxml import etree

ROOT = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
OUT = ROOT / 'experiment_completion/argument_revision_v23_20261012/manuscript_QA'
TXT = OUT / 'v23_content.txt'
DOCX = ROOT / 'paper_v23_20261012_中文修订_论证精简与补充整合_三线表.docx'
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
      'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math'}


def hash_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def citations(text):
    result = []
    for bracket in re.findall(r'\[([0-9][0-9,，\s–—\-]*)\]', text):
        for token in re.split(r'[,，]', bracket):
            token = token.strip()
            match = re.fullmatch(r'(\d+)\s*[–—\-]\s*(\d+)', token)
            if match:
                result.extend(range(int(match.group(1)), int(match.group(2)) + 1))
            elif token.isdigit():
                result.append(int(token))
    return result


def figure_table_refs(text):
    result = []
    for match in re.finditer(r'([图表])\s*(S?\d+)(?:\s*[–—\-]\s*(S?\d+))?', text):
        kind, start, end = match.groups()
        if end and start.startswith('S') == end.startswith('S'):
            prefix = 'S' if start.startswith('S') else ''
            result.extend((kind, prefix + str(n)) for n in range(int(start.lstrip('S')), int(end.lstrip('S')) + 1))
        else:
            result.append((kind, start))
            if end:
                result.append((kind, end))
    return result


docx_hash_before = hash_file(DOCX)
with zipfile.ZipFile(DOCX) as z:
    x = etree.fromstring(z.read('word/document.xml'))
    ps = x.findall('.//w:body/w:p', NS)
    paragraph_texts = [''.join(p.xpath('.//w:t/text()', namespaces=NS)) for p in ps]
    all_paragraph_texts = [''.join(p.xpath('.//w:t/text() | .//m:t/text()', namespaces=NS)) for p in ps]
    math_indices = [index for index, p in enumerate(ps) if p.find('.//m:oMath', NS) is not None]
    table_nodes = x.findall('.//w:body/w:tbl', NS)
    tables = [[[''.join(c.xpath('.//w:t/text()', namespaces=NS)) for c in row.findall('w:tc', NS)]
               for row in table.findall('w:tr', NS)] for table in table_nodes]
    image_count = len(x.xpath('//a:blip', namespaces=NS))
    fields = x.xpath('//w:instrText/text()', namespaces=NS)
    comment_hash = hashlib.sha256(z.read('word/comments.xml')).hexdigest() if 'word/comments.xml' in z.namelist() else None
references = [(index, text) for index, text in enumerate(paragraph_texts) if re.match(r'^\[\d+\]\s', text)]
reference_ids = [int(re.match(r'^\[(\d+)\]', text).group(1)) for _, text in references]
body = '\n'.join(text for text in paragraph_texts if not re.match(r'^\[\d+\]\s', text))
cited_ids = sorted(set(citations(body)))
captions = []
for index, text in enumerate(paragraph_texts):
    match = re.match(r'^([图表])\s*(S?\d+)\s*[:：]', text)
    if match:
        captions.append({'kind': match.group(1), 'id': match.group(2), 'paragraph': index})
caption_set = {(c['kind'], c['id']) for c in captions}
missing_crossrefs = sorted(set(figure_table_refs(body)) - caption_set)
duplicate_captions = sorted({key for key in caption_set if sum((c['kind'], c['id']) == key for c in captions) > 1})
supplement_modules = [match.group(1) for text in paragraph_texts if (match := re.match(r'^([A-Z])\s+\S', text))]
supplement_refs = sorted(set(re.findall(r'补充材料\s*([A-Z])', body)))
source_text = TXT.read_text(encoding='utf-8')
source_paragraphs = {int(m.group(1)): m.group(2) for line in source_text.splitlines()
                     if (m := re.match(r'^P(\d+):\s?(.*)', line))}
plain_export_mismatches = [index for index, text in enumerate(paragraph_texts)
                          if index not in math_indices and text
                          and re.sub(r'\s+', '', text) != re.sub(r'\s+', '', source_paragraphs.get(index, ''))]

checks = {
    'references_1_through_35_complete_and_cited': reference_ids == list(range(1, 36)) and cited_ids == list(range(1, 36)),
    'main_figures_1_through_9_and_only_supplement_S1': {c['id'] for c in captions if c['kind'] == '图'} == {str(n) for n in range(1, 10)} | {'S1'},
    'main_tables_1_2_and_supplement_S1_through_S9': {c['id'] for c in captions if c['kind'] == '表'} == {'1', '2'} | {'S' + str(n) for n in range(1, 10)},
    'no_dangling_crossrefs': not missing_crossrefs,
    'no_duplicate_captions': not duplicate_captions,
    'supplement_modules_A_through_G': supplement_modules == list('ABCDEFG'),
    'supplement_letter_references_connect': set(supplement_refs).issubset(set(supplement_modules)),
    'plain_docx_text_export_corresponds': not plain_export_mismatches,
    'released_RMSF_results_honestly_retained': all(s in body for s in ['总体成绩以BackFlip2和PLANET full较好', '0.537 Å/0.848', '0.568 Å/0.828', '0.577 Å/0.796']),
    'discordant_reference_presence_separated_from_low_recall': all(s in body for s in ['47,354个', '这些参考组合', '15.8%和16.5%', '当前模型对不同组合的恢复程度存在差异']),
    'md_Fstar_separate_from_model_scalar_outputs': '模型输出RMSF与Hχ，轨迹提供占据概率及对应的F*曲线' in body,
    'Fstar_gaps_density_and_dimensionless_scale_retained': all(s in body for s in ['频率除以箱宽', '零计数箱断线', '无量纲热能标度']),
    'main_partition_section_keeps_only_core_jointH_CA_numbers': all(s in body for s in ['1.484和0.962', '6.2%和7.8%']),
    'RMSF2_not_RMSF_is_additively_decomposed': 'RMSF²可以分解' in body,
    'sampling_change_JSD_details_moved_to_supplement_E': all(s in body for s in ['189.2/21.9', '0.050/0.067', '0.021/0.040', '0.022/0.040']),
    'KIX_delta_direction_failures_retained': '六个状态差值的预测方向均与单轨迹参考相反' in body and all(s in body for s in ['+1.017与−1.589', '+0.014与−1.263']),
    'native_coverage_values_retained': all(s in body for s in ['51.76% ± 2.71%', '85.19% ± 1.31%']),
    'seed_repeats_explicitly_original_system_fold0': '固定原体系划分的第0折训练与留出集合，以种子7、11和23分别训练完整与链内模型' in body,
    'tableS2_caption_states_fold0_and_three_seeds': '表S2：原体系划分第0折的训练随机初始化对照。相同训练与留出记录分别以三个种子训练40轮' in body,
    'single_target_models_distinguished_from_seed_fold0_scope': '单目标模型每目标训练五个' in body,
    'seed_report_describes_each_individual_holdout_model': '固定划分的随机初始化对照报告每个初始化模型的留出成绩' in body,
    'conclusion_restricts_average_gains_to_fivefold_system_holdout': any('五折体系留出评价中' in t and '完整几何环境' in t and re.search(r'(?:相同|同)训练集合', t) for t in paragraph_texts if t.startswith('FlexPPI从静态主链结构联合估计')),
    'representative_case_explanation_precedes_figure9': next(i for i, t in enumerate(paragraph_texts) if t.startswith('逐角F*进一步区分')) < next(i for i, t in enumerate(paragraph_texts) if t.startswith('图9：')),
    'family_protocol_Itrace_fold0_retained': all(s in body for s in ['全部244条ITrace位于第0折测试集合', '训练集合来自ATLAS和DynaRepo']),
    'planet_mode_dtype_unit_retained': all(s in body for s in ['encoded结构特征', '300 K与float32', 'nm乘以10转换为Å']),
    'rmsfnet_simulated_density_relative_score_retained': all(s in body for s in ['4 Å分辨率、1.5 Å网格', '高斯模拟密度', '相对柔性分数']),
    'ref16_only_future_local_loop_use': [index for index, text in enumerate(paragraph_texts) if '[16]' in text and not text.startswith('[16]')] == [118],
    'no_stale_calibration_supplement_reference': not any('校准见' in text or '补充材料H' in text or '补充材料I' in text or '补充材料K' in text or '补充材料L' in text for text in paragraph_texts),
    'no_unresolved_word_REF_fields': not any(re.search(r'\b(?:REF|PAGEREF)\b', field) for field in fields),
}

substantive_issues = []
minor_issues = []
for index, text in enumerate(paragraph_texts):
    if '两項' in text:
        minor_issues.append({'paragraph': index, 'original': '两項', 'replace': '两项', 'kind': 'mixed_script'})
    if '各自产首帧' in text:
        minor_issues.append({'paragraph': index, 'original': '各自产首帧', 'replace': '各自生产首帧', 'kind': 'missing_character'})
    if '状态間' in text:
        minor_issues.append({'paragraph': index, 'original': '状态間', 'replace': '状态间', 'kind': 'mixed_script'})
    if 'Pro及原联合标签保持原定义' in text:
        minor_issues.append({'paragraph': index, 'original': 'Pro及原联合标签保持原定义',
                             'replace': '所有原定义χ角及Pro保留在联合分区敏感性分析中，训练标签保持冻结',
                             'kind': 'possible_partition_scope_ambiguity'})
for name, passed in checks.items():
    if not passed:
        substantive_issues.append({'check': name, 'kind': 'requires_inspection'})
docx_hash_after = hash_file(DOCX)
checks['docx_unchanged_during_review'] = docx_hash_before == docx_hash_after
if not checks['docx_unchanged_during_review']:
    substantive_issues.append({'check': 'docx_unchanged_during_review', 'kind': 'candidate_changed'})

report = {
    'review': 'Independent read-only v23 logic, wording, references and crossreferences',
    'timestamp_local': datetime.now().isoformat(timespec='seconds'),
    'docx': str(DOCX), 'docx_sha256': docx_hash_after, 'text': str(TXT), 'text_sha256': hash_file(TXT),
    'status': 'requires_inspection' if substantive_issues else 'minor_edits' if minor_issues else 'pass',
    'substantive_issues': substantive_issues, 'minor_issues': minor_issues,
    'counts': {'paragraphs': len(ps), 'tables': len(tables), 'images': image_count, 'references': len(references)},
    'checks': checks, 'reference_ids': reference_ids, 'cited_ids': cited_ids,
    'uncited_reference_ids': sorted(set(reference_ids) - set(cited_ids)),
    'undefined_reference_ids': sorted(set(cited_ids) - set(reference_ids)),
    'captions': captions, 'missing_crossrefs': missing_crossrefs, 'duplicate_captions': duplicate_captions,
    'supplement_modules': supplement_modules, 'supplement_refs': supplement_refs,
    'plain_export_mismatches': plain_export_mismatches, 'field_instructions': fields,
    'comments_xml_sha256': comment_hash,
    'equation5_manual_verification': 'F_c*=-ln(p_c/pmax)=-lnp_c+lnpmax; sum_c p_cF_c*=Hchi+lnpmax, hence Hchi=sum_c p_cF_c*−lnpmax, exactly as written for one shared category distribution.',
    'manual_logic_findings': [
        'Reference discordant categories and prediction recovery are explicitly separated, and low15.8%/16.5% discordant recall remains visible.',
        'RMSF released models BF2 and PLANETfull are explicitly identified as better overall MAE/rho; no tail-RMSE reversal of this comparison is asserted.',
        'New free-energy content stays in the main text as MD-derived empirical occupancy explanation; no model-predicted density, barrier or peak-location capability is claimed.',
        'Main3.6.2 is compressed to jointH1.484/.962 and state-conditionedCA6.2%/7.8%; changes/JSD/window details are retained in supplementE.',
        'The three discussion paragraphs each perform a distinct role: geometric dual-output value, occupancy interpretation, subsequent use/research direction.',
        'Reference16 only accompanies future comparison of antibody local-loop conformations before/after binding.',
        'KIXwithin-state localization and failed between-state delta directions are both retained in their respective main/supplement roles.',
        'Final supplementB correctly states the three-seed comparison uses the original system partition fold0 training/holdout sets; single-target five-model-per-target training remains separately stated.',
        'Final methodP49 reports each individual initialization model holdout score, matching single-model rows in tableS2; conclusion explicitly scopes the gains to fivefold system holdout evaluation.',
        'The two representative-case explanation paragraphs now precede Figure9; content and reference connections remain correct after this layout reorder.',
        'No substantive defensive repetition remains; main text repeats provenance only where necessary to identify what models and trajectories provide.'
    ],
    'docx_modified': False, 'network_used': False
}
OUT.mkdir(exist_ok=True)
(OUT / 'v23_independent_content_review.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'status': report['status'], 'docx_sha256': docx_hash_after, 'counts': report['counts'],
                  'failed_checks': [name for name, passed in checks.items() if not passed],
                  'minor_issues': minor_issues,
                  'output': str(OUT / 'v23_independent_content_review.json')}, ensure_ascii=False, indent=2))
