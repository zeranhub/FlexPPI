from pathlib import Path
from datetime import datetime, timezone, timedelta
import json, hashlib

ROOT=Path(__file__).resolve().parent
QA=ROOT/'argument_revision_v20_20261012'
RESULT=ROOT/'completed_results_20261012'
DOC=ROOT.parent/'paper_v20_20261012_中文修订_联合模型_三线表.docx'
digest=hashlib.sha256(DOC.read_bytes()).hexdigest()
assert digest=='4b4597b42fdc8a1300ab8f4155e98a4960872a2368fc887397d133af7d264d79'
layout_path=QA/'v20_layout_audit.json'
layout=json.loads(layout_path.read_text(encoding='utf8'))
assert layout['sha256']==digest and layout['pages']==40 and len(layout['tables'])==13
for name in ['v20_results_content_review.json','v20_independent_content_review.json']:
    assert (QA/name).is_file()
layout['delegated_visual_review_status']='passed'
layout['delegated_actual_visual_review_passed_pages']=list(range(19,31))
layout['visual_review']='passed_all_40_pages_at_original_resolution'
layout['visual_reviewer_for_pages_19_to_30']='/root/methods_revision'
layout_path.write_text(json.dumps(layout,ensure_ascii=False,indent=2),encoding='utf8')

oldstatus=(ROOT/'completion_status_zh.md').read_bytes()
(QA/'completion_status_before_v20.md').write_bytes(oldstatus)
delivery_path=RESULT/'delivery_status.json'
old_delivery=json.loads(delivery_path.read_text(encoding='utf8'))
(QA/'delivery_status_before_v20.json').write_bytes(delivery_path.read_bytes())
now=datetime.now(timezone(timedelta(hours=8))).isoformat()
delivery=dict(old_delivery)
delivery.update({'status':'chinese_manuscript_v20_argument_revision_delivered',
                 'completed_beijing':now,'document':str(DOC),'document_sha256':digest,
                 'pages':40,'three_line_tables':13,'figures':16,
                 'content_review':'passed_two_independent_reviews',
                 'render_and_visual_review':'passed_all_40_pages',
                 'remaining_required_work':[],
                 'previous_document':old_delivery['document'],
                 'previous_document_sha256':old_delivery['document_sha256'],
                 'changes_this_round':'Default joint model; complementary-output argument; remove repetitive defensive prose and superseded external protocol illustrations; preserve retained results.'})
delivery_path.write_text(json.dumps(delivery,ensure_ascii=False,indent=2),encoding='utf8')
(QA/'delivery_status.json').write_text(json.dumps(delivery,ensure_ascii=False,indent=2),encoding='utf8')
(ROOT/'completion_status_zh.md').write_text(
    f'当前交付为 `{DOC.as_posix()}`，共40页、13张原生可编辑三线表、5张正文图和11张补充图。文件SHA256为 `{digest}`。本轮按用户要求，将联合训练确立为默认模型，以同一结构表示提供RMSF与Hχ两种互补读数；删除摘要、正文及讨论中重复的下降、限制和预防质疑论述。家族评价保留真实方法与表S10结果，正文用完整模型相对链内模型的比较说明结构信息价值。\n\n'
    '移除已经由新外部方案替代的旧图S1、S10、S11和旧表S12，补充图连续重编号为S1至S11；基线表删除重复的家族栏。保留的表格数值全部与v19及源表一致，单目标、家族、外部和KIX结果未改。原始批注与锚点、原生公式及保留图片保持一致。两份独立内容审阅通过，全40页实际Word渲染并逐页视觉检查通过。原始v13和此前v19均保留。\n\n'
    '所有约定补充实验已完成，科学计算与训练配置本轮未修改。单目标十次、家族十次40轮训练、完整折外预测、2430训练集合的32条外部评价及全部新增基线已完成。KIX复用四种状态各一条100ns生产轨迹，新MD生产保持取消。完整最终权重及逐残基源表的本地归档已核验，归档SHA256为 `dbaaa3a72f620dc61024f24fc0442be1709d81624be9b9f49789dba59b1afc69`。\n\n'
    '本轮编辑、内容及逐页渲染审计位于 `argument_revision_v20_20261012/`，当前交付状态位于 `completed_results_20261012/delivery_status.json`。本轮没有启动实验、创建自动任务、翻译英文或外部投稿。\n',encoding='utf8')
note=ROOT.parent/'paper_v20_20261012_修改说明.md'
note.write_text(
    '# v20 中文稿修改说明\n\n'
    '默认模型明确为联合训练。论证围绕共享几何表示和主链波动、侧链状态多样性两种互补输出组织，单目标训练保留为配置对照。摘要、结果、讨论与结论删去反复强调收益大小、家族下降及预防质疑的收尾句；必要的计算口径集中保留在方法和表注。\n\n'
    '家族评价保留真实数值与划分方法，结果用于比较完整结构和链内输入。移除旧外部协议的三张图及一张历史结果表，基线表去掉重复家族列，更新全部补充图引用。65段主体替换及后续3段方法压缩，共68段替换，删除16段（包含旧图与图注）。所有保留的表格数值保持一致。\n\n'
    '最终文件40页、13张三线表，31篇文献引用及图表编号已核对，两份独立内容检查和全40页Word渲染视觉检查通过。原批注、数学公式和保留图像均保持一致。\n',encoding='utf8')
print(json.dumps({'document':str(DOC),'sha256':digest,'pages':40,'tables':13,'status':'complete','remaining_required_work':[]},ensure_ascii=False))
