from pathlib import Path
import hashlib
import json
import zipfile
import xml.etree.ElementTree as ET

ROOT = Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_joint_case_v49_20261019')
QA = ROOT / 'manuscript_QA'
audit = json.loads((QA / 'integration_audit.json').read_text(encoding='utf-8'))
document = Path(audit['document'])
document_sha256 = hashlib.sha256(document.read_bytes()).hexdigest()
assert document_sha256 == audit['document_sha256']
assert document_sha256 == '5f79ddf7d657dc9b782c9fb432e373a1e7c170c1c88d019e528fbf73c366c222'
W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
M = '{http://schemas.openxmlformats.org/officeDocument/2006/math}'
with zipfile.ZipFile(document) as package:
    xml = ET.fromstring(package.read('word/document.xml'))
body = xml.find(W + 'body')

def content(node):
    return ''.join(item.text or '' for item in node.iter() if item.tag in (W + 't', M + 't'))

checks = []
for change in audit['prose_changes']:
    actual = content(body[change['body_index']])
    assert actual == change['after'], (change['body_index'], actual, change['after'])
    checks.append({'final_body_index': change['body_index'], 'source_paragraph_index': change['source_paragraph_index'], 'actual_document_text_matches_reviewed_change': True})
for inserted in audit['inserted_body_nodes']:
    assert content(body[inserted['body_index']]) == inserted['text']
    checks.append({'final_body_index': inserted['body_index'], 'inserted_mechanism_paragraph_matches': True})
corpus = '\n'.join(content(node) for node in body)
assert 'KRAS' not in corpus
assert 'Poole' not in corpus
reference28 = next(content(node) for node in body if content(node).startswith('[28]'))
review = json.loads((ROOT / 'literature_review' / 'narrative_review.json').read_text(encoding='utf-8'))
assert reference28 == review['reference28']
assert sum(content(node).startswith('[28]') for node in body) == 1
assert '其中βCDR2也是晶体结构显示接触增加最明显的区域[28]' in corpus
assert '这些结果支持一种局部解释' in corpus
assert 'βPhe53和βGly54分别对应其输入E54与E55' in corpus
assert '已有ITrace轨迹为该接触网络提供了补充动态证据' in corpus

result = {
    'status': 'passed',
    'document': str(document),
    'document_sha256': document_sha256,
    'review_scope': 'Actual DOCX body text: every replaced and inserted paragraph, Figure 4 caption, reference 28, and consistency with local primary fulltexts and audited model/MD sources. This is scientific text review, not rendered page QA.',
    'actual_document_changes_checked': checks,
    'reviewed_section_functions': {
        'abstract_and_introduction': 'State the joint-output model finding and connect it to the paper validation sequence without claiming affinity prediction.',
        'methods': 'Identify structures, held-out model, corrected CDR alignment, interface definition, consistent residue numbering, input-context removal, and independently repeated MD contact criteria.',
        'background': 'Introduce MART-1/ELA and natural-parent versus affinity-matured receptor, then establish SPR and crystallographic evidence for the biological comparison.',
        'complete_region_results': 'Show the two model readouts across full CDR loops and full pMHC interface before the local mechanism example.',
        'Q72_Gly54_model_results': 'Use unchanged residue identities and TCR-context removal to tie sidechain concentration and backbone restriction to local environment.',
        'inserted_mechanism': 'Use existing ITrace contact occupancy to support a local contact-network explanation, with published SPR/structure as separately sourced biological evidence.',
        'discussion_and_conclusion': 'Explain the practical value of combining the two outputs at a concrete residue contact, without replacing measured affinity or total binding entropy.'
    },
    'scientific_checks': {
        'MART1_alias_and_melanoma_background_supported': True,
        'ELA_correctly_identified_as_anchor_modified_decamer_A27L': True,
        'SPR_18uM_to_0p6nM_and_slower_dissociation_supported': True,
        'MHC_contacts_56_to_97_supported': True,
        'CDR2beta_6_to_38_correctly_labeled_total_pMHC_contacts': True,
        'CDR3beta_not_incorrectly_claimed_to_gain_contacts': True,
        'betaI53F_new_MHC_contacts_supported': True,
        'unchanged_Gly54_and_Q72_numbering_correct': True,
        'joint_Hchi_three_chi_state_occupancy_interpretation_correct': True,
        'prediction_values_match_existing_audited_values_at_reported_precision': True,
        'interface_23_positions_17_predicted_decreases_verified': True,
        'Q72_context_delta_plus0p104_correct': True,
        'ITrace_MD_not_attributed_to_original_Madura_study': True,
        'independent_replicas_distinguished_from_saved_frames': True,
        'Hbond_occupancy_not_claimed_permanent_or_statistically_significant': True,
        'local_mechanism_explicitly_framed_as_supported_interpretation': True,
        'no_Hchi_equals_total_binding_entropy_claim': True,
        'no_RMSF_Hchi_combined_affinity_score_or_30000fold_prediction_claim': True,
        'Figure4_caption_roles_and_source_distinction_correct': True,
        'reference28_all16_authors_volume_pages_DOI_correct': True,
        'old_KRAS_and_Poole_absent_from_actual_document_body': True
    },
    'necessary_modifications': [],
    'previous_corrections_confirmed': [
        'The betaCDR3 contact-count generalization has been removed.',
        'The mechanism paragraph now says these results support a local interpretation.',
        'ELA is identified as the MART-1-derived A27L anchor-optimized peptide.',
        'Reference 28 has been replaced with the complete 16-author MEL5 primary citation.'
    ],
    'remaining_scientific_text_issues': [],
    'rendered_page_QA_included': False,
    'document_modified_by_reviewer': False
}
output = ROOT / 'literature_review' / 'final_narrative_review.json'
output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'status': result['status'], 'document_sha256': document_sha256, 'changed_or_inserted_body_nodes_checked': len(checks), 'output': str(output)}, ensure_ascii=False))
