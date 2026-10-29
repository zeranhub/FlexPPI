"""Read-only, hash-bound narrative audit for dense Figure 2 revision."""
import argparse, hashlib, json
from pathlib import Path
from zipfile import ZipFile
from lxml import etree
ROOT=Path(__file__).resolve().parents[1]
V44=ROOT.parent/'predictive_figure_v44_20261018'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--new',required=True,type=Path);ap.add_argument('--figure2',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
 package=json.loads((ROOT/'manuscript_QA/independent_package_audit.json').read_text(encoding='utf8'))
 numeric=json.loads((ROOT/'text_review/independent_metrics_audit.json').read_text(encoding='utf8'))
 proposal=json.loads((ROOT/'text_review/manuscript_proposal.json').read_text(encoding='utf8'))
 summary=json.loads((V44/'analysis/summary.json').read_text(encoding='utf8'))
 with ZipFile(args.new) as z:tree=etree.fromstring(z.read('word/document.xml'))
 ps=[''.join(p.xpath('.//w:t/text()',namespaces=NS)) for p in tree.xpath('//w:body/w:p',namespaces=NS)]
 find=lambda prefix:next(p for p in ps if p.startswith(prefix))
 method=find('高动态位点筛选以每条记录MD参考');result=find('进一步将连续预测用于高动态位点筛选');within=find('在同一记录、同一种氨基酸内重新定义高熵位置后');source=find('按来源分别汇总高动态筛选结果');caption=find('图2：')
 reviewed=[p['new'] for p in proposal['replacements']]+[i['text'] for i in proposal['insertions']]+[caption]
 numerical_strings=[f"{summary['metrics'][t][m]['mean_ap']:.3f}" for t in ['rmsf','chi_entropy'] for m in ['Full','Intra','Static','AA']]
 numerical_strings += [f"{summary['metrics'][t]['Full']['mean_auroc']:.3f}" for t in ['rmsf','chi_entropy']]
 ci_strings=[f'{v:.3f}' for t in ['rmsf','chi_entropy'] for v in summary['metrics'][t]['Full']['ap_ci']]
 type_strings=[f"{summary['within_type'][m]['mean_ap']:.3f}" for m in ['Full','Intra','Static','AA']]+[f"{summary['within_type']['Full']['mean_auroc']:.3f}"]
 checks={
 'package_audit_passed_same_document':package['status']=='passed' and package['document_sha256']==sha(args.new),
 'numeric_audit_passed_and_all_output_hashes_bound':numeric['status']=='passed' and all(sha(path)==h for path,h in numeric['outputs_sha256'].items()),
 'proposal_v44_source_SHA_bound':proposal['source_sha256']=='4156791ba646d8137754386ea387b920aaaf1bdac041614414a17028915ffcab',
 'all_reviewed_text_present_once':all(ps.count(s)==1 for s in reviewed),
 'AP_q90_ties_definitions_unchanged':'第90百分位及以上为阳性' in method and '并列参考全部纳入' in method and '并列预测作为同一阈值组' in method and '召回率增量乘对应精确率求和' in method,
 'AUROC_definition_unchanged':'所有阳性–阴性位点对' in method and '并列计半分' in method,
 'record_unit_bootstrap_and_curve_AP_distinction':'对记录等权平均' in method and '来源分层重抽记录10,000次' in method and '曲线面积不代替平均AP' in method,
 'type_to_record_curve_aggregation_explicit':'指标与曲线均先对有效类型、再对记录等权平均' in method,
 'source_AP_definition_explicit':'来源分层AP在各来源内对记录等权平均' in method,
 'original_main_AP_AUC_and_CIs_preserved':all(s in result for s in numerical_strings+ci_strings),
 'within_type_AP_AUC_and_coverage_preserved':all(s in within for s in type_strings) and '3,063' in within and '42,108' in within,
 'new_type_PR_reference_correct':'图2e' in within,
 'new_source_matrix_reference_correct':'图2f' in source,
 'source_Full_values_correct':all(s in source for s in ['0.746','0.664','0.840','0.922','0.924','0.968']),
 'source_Full_advantage_supported':numeric['checks']['all_6_source_target_Full_AP_higher_than_controls'],
 'new_source_analysis_not_called_source_holdout':not any('来源留出' in s or 'leave-source-out' in s for s in reviewed),
 'removed_sparse_panels_not_referenced_as_type_gain_figures':'图2' not in find('氨基酸类型均值对照的整体Hχ相关已达') and '图2' not in find('RMSF方面，完整模型在20种氨基酸类型'),
 'type_MAE_strength_and_incomplete_Intra_gains_preserved':'18种' in find('氨基酸类型均值对照的整体Hχ相关已达') and '16种' in find('氨基酸类型均值对照的整体Hχ相关已达') and '17种' in find('氨基酸类型均值对照的整体Hχ相关已达'),
 'existing_RMSF_method_strength_preserved':'低于FlexPPI的0.577' in find('在已有RMSF方法中，BackFlip'),
 'other_images_cases_background_tables_math_comments_preserved':package['checks']['all_other_body_nodes_exactly_preserved'] and package['checks']['six_images_same_order_and_Fig5_unchanged'],
 'caption_matches_final_asset':caption==(ROOT/'figure2/caption_zh.md').read_text(encoding='utf8').strip(),
 }
 out={'status':'passed' if all(checks.values()) else 'failed','document':str(args.new),'document_sha256':sha(args.new),'figure2':str(args.figure2),'figure2_sha256':sha(args.figure2),'checks':checks,'numeric_audit_sha256':sha(ROOT/'text_review/independent_metrics_audit.json'),'package_audit_sha256':sha(ROOT/'manuscript_QA/independent_package_audit.json'),'scope':'No training, inference, MD or test-set threshold optimization. Complete Full density and original OOF rankings retained. Source-stratified comparisons are not source-held-out validation. Page rendering handled separately.'}
 args.output.parent.mkdir(exist_ok=True,parents=True);args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps({'status':out['status'],'failed':[k for k,v in checks.items() if not v]},ensure_ascii=False))
if __name__=='__main__':main()
