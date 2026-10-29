"""Read-only science and figure-reference audit for compact v55 manuscript."""
from pathlib import Path
import hashlib,json,zipfile
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parent.parent
EXP=ROOT.parent
BASE=EXP.parent
OLD=EXP/'mel5_composite_v54_20261022'
DOC=BASE/'paper_v55_20261022_中文修订_MEL5精简联合图4_三线表.docx'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(DOC)=='056c2500e95d3992cf25f31980751b9936587c3dbe9faee6cf9738c88388153c'
prior=json.loads((OLD/'manuscript_QA/independent_science_review.json').read_text(encoding='utf-8'))
assert prior['status']=='passed'
oldprose=json.loads((OLD/'prose_v54.json').read_text(encoding='utf-8'))
prose=json.loads((ROOT/'prose_v55.json').read_text(encoding='utf-8'))
with zipfile.ZipFile(DOC) as archive:
    doc=ET.fromstring(archive.read('word/document.xml'))
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    paragraphs=[''.join(p.itertext()) for p in doc.findall('.//w:p',ns)]
    for value in prose.values():assert value in paragraphs
for key,value in oldprose.items():
    if key not in prose:assert value in paragraphs
assert prose['91']==oldprose['91'].replace('（图4d）','（图4d下部）')
assert prose['93']==oldprose['93'].replace('（图4e）','（图4d上部）').replace('（图4f）','（图4d下部）')
caption=prose['95']
for required in ['（a）','（b）','（c）','（d）','上部为CDR ΔRMSF','下部为pMHC界面ΔHχ','两项分别标示单位和色标','4JFF−3HG1','4JFD−4JFF','4JFE−4JFF','蓝、红表示降低、升高','灰格表示Hχ无定义或氨基酸类型不同','最近Cα距离≤8 Å','共30位','三行Hχ有效数为23、23、22','TCR采用MEL5共同编号','同一模型对已有复合物结构预测']:
    assert required in caption,required
for removed in ['（e','（f','排序','p表示肽']:
    assert removed not in caption
assert all('图4e' not in p and '图4f' not in p for p in paragraphs)
data=pd.read_csv(ROOT/'figure4/source_interface_three_contrasts.csv',float_precision='round_trip')
q=data[data.pdb_before.eq('3hg1')&data.pdb_after.eq('4jff')&data.hchi_valid&data.same_aa]
assert len(q)==23
anchor=q.loc[q.hchi_change_after_minus_before.idxmin()]
assert anchor.chain_before=='A' and int(anchor.id_before)==72
assert f'{anchor.hchi_full_before:.3f}'=='1.792' and f'{anchor.hchi_full_after:.3f}'=='1.375'
arrays=np.load(ROOT/'figure4/source_heatmap_arrays.npz',allow_pickle=False)
assert arrays['CDR_delta_RMSF'].shape==(3,45)
assert arrays['interface_delta_Hchi'].shape==(3,30)
assert np.isfinite(arrays['interface_delta_Hchi']).sum(axis=1).tolist()==[23,23,22]
figreview=json.loads((ROOT/'figure4/independent_data_review.json').read_text(encoding='utf-8'))
assert figreview['status']=='passed' and figreview['all_remaining_scientific_values_identical_to_reviewed_v54']
report={
 'status':'passed','DOCX_path':str(DOC),'DOCX_SHA256':sha(DOC),'prose_JSON_SHA256':sha(ROOT/'prose_v55.json'),
 'read_only_review':True,'all3_modified_paragraphs_exactly_present_in_DOCX':True,'other5_previously_reviewed_case_paragraphs_unchanged_and_exactly_present':True,
 'P91':{'only_change':'Figure4d -> Figure4d lower block','all23_valid_sites_retained_in_heatmap':True,'Q72_largest_Hchi_decrease_rechecked':True,'Q72_Hchi_1_792_to_1_375_nats_preserved':True,'CDR_and_Gly54_RMSF_numbers_unchanged':True},
 'P93':{'only_change':'Figure4e/f -> Figure4d upper/lower blocks','CDR2alpha_CDR3alpha_regional_response_claim_unchanged_from_passed_v54':True,'MD_occupancies_unchanged':True},
 'caption':{'four_panels_a_to_d':True,'old_ranked_bar_and_p_prefix_description_removed':True,'d_upper_CDR_deltaRMSF_lower_pMHC_deltaHchi':True,'separate_units_and_colorbars':True,'three_contrast_references_correct':True,'30_interface_sites_23_23_22_mask_rule_correct':True,'common_numbering_and_observed_input_scope_correct':True},
 'fixed_atom_distance_statement_preserved':{'WT_A':['5.85','10.58'],'mature_A':['2.94','3.12'],'metric':'Q72 OE1–beta53 backboneN and Q72 NE2–beta54 backboneO; unchanged from independently recomputed v54 science review'},
 'supporting_prior_science_review':{'path':str(OLD/'manuscript_QA/independent_science_review.json'),'DOCX_SHA256':prior['DOCX_SHA256'],'report_SHA256':sha(OLD/'manuscript_QA/independent_science_review.json')},
 'blocking_issues':[],'required_scientific_edits':[]
}
(ROOT/'manuscript_QA/independent_science_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','DOCX_SHA256':sha(DOC),'file':str(ROOT/'manuscript_QA/independent_science_review.json')},ensure_ascii=False))
