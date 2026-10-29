"""Commit the completed render/source review and update the durable handoff."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json

ROOT = Path(__file__).resolve().parent
DOCX = ROOT.parent / 'paper_v18_20261010_中文修订_复用KIX单轨迹_三线表.docx'
digest = hashlib.sha256(DOCX.read_bytes()).hexdigest()
layout_path = ROOT / 'v18_layout_audit.json'
layout = json.loads(layout_path.read_text(encoding='utf8'))
review = json.loads((ROOT / 'kix_existing_reuse_20261010/v18_content_review.json').read_text(encoding='utf8'))
reuse = json.loads((ROOT / 'kix_existing_reuse_20261010/reuse_audit.json').read_text(encoding='utf8'))
assert digest == layout['sha256'] == review['reviewed_docx_sha256']
assert layout['pages'] == 40 and len(layout['tables']) == 10 and not layout['blank_pages']
assert review['issues'] == []
assert reuse['raw_file_hashes_rechecked_now'] is True
visually_inspected = [1, 2, 7, 8, 13, 14, 16, 17, 18, 39, 40]
assert layout['changed_pages_requiring_visual_review'] == visually_inspected
layout.update(visual_review='passed', inspected_changed_pages=visually_inspected,
    unchanged_page_review='All remaining page content is pixel-identical to previously rendered and visually verified v17, excluding page-number footer.',
    review_completed_utc=datetime.now(timezone.utc).isoformat(),
    findings='No blank pages, clipped text, missing glyphs or table overflow. Tables S7/S8 stay on separate complete pages; all 10 tables have three horizontal rules only.')
layout_path.write_text(json.dumps(layout, ensure_ascii=False, indent=2), encoding='utf8')
integration_path = ROOT / 'kix_single_integration_audit.json'
integration = json.loads(integration_path.read_text(encoding='utf8'))
assert integration['output_sha256'] == digest
integration.update(render_and_independent_source_audit_pending=False,
                   layout_audit=str(layout_path), content_review='passed', reuse_audit='passed')
integration_path.write_text(json.dumps(integration, ensure_ascii=False, indent=2), encoding='utf8')

status_path = ROOT / 'completion_status_zh.md'
status = status_path.read_text(encoding='utf8')
paragraphs = status.split('\n\n')
paragraphs[0] = ('当前可审阅稿件为 `E:/Project/FlexPPI/paper/revision_v14_20261010/' + DOCX.name +
    '`。稿件共40页、10张三线表，已完成Word实际渲染核查；原批注、批注锚点与既有图片保留。当前文件SHA256为 `' + digest +
    '`，排版检查记录为 `v18_layout_audit.json`。新增补充材料F、表S7/S8，使用2026年9月四状态各一条100ns的完整KIX生产轨迹；400段DCD当前哈希及全部质量分析核验通过。独立内容审查通过。稿件如实区分活动区定位与响应预测：六组L12平均RMSF状态差值均与单轨迹MD参考方向相反，且前后半程存在时段依赖，不宣称收敛、独立重复不确定性或变构机制验证。此为已完成分析的阶段性中文稿，单目标与家族隔离模型训练和新增基线全队列评价尚未完成。')
status = '\n\n'.join(paragraphs)
status = status.replace('GPU主进程PID71135', 'GPU监控父进程已按2026-10-10最新用户要求安全接管为PID76963，原训练队列PID71188保留')
status = status.replace('初始单目标任务为chi_entropy_only fold0，2026-10-10 15:39附近已达到第14轮；后续进度以服务器日志为准。',
    '北京时间2026-10-10 18:14快照中，单目标训练已完成1/10项，当前RMSF-only第0折到第38/40轮，家族隔离训练0/10项；当前日志未观察到异常，后续进度以服务器日志为准。')
status = status.replace('2026-10-10 15:41快照观察RMSF-net172条记录、PLANET full65条链、seq69条链完成，均没有已观测失败；这些计数是进度，不是完整结果。',
    '北京时间2026-10-10 18:14快照观察RMSF-net851/3,063条记录、PLANET full2,169/4,973条链、seq2,868/4,973条链完成，均没有已观测失败；这些计数是进度，不是完整结果，三项均尚无最终独立核验完成证书。')
status = status.replace('未取得全队列输出与重复轨迹时，不填写这些结果。', '未取得完整训练、预测和全队列核验输出时，不填写相应成绩；KIX单轨迹分析已经完成，不再等待新的重复轨迹。')
if 'hmu-gpu的SSH服务可达但现有公钥认证被拒绝' not in status:
    status += '\n\n本轮历史资产搜索已完成可访问范围：本地权威轨迹在 `E:/Project/FlexPPI/artifacts/external/kix_md_20260914/{apo,myb,mll,both}/production/`，各100段、100,000帧、100ns，共约8.78GB。原始指标及全帧质量分析在 `artifacts/kix_md_20260914/pilot_results_20260916/`。400段当前SHA、帧数、蛋白原子数、文件大小及五份质量分析表哈希均已核查。模型只使用 `experiments/_case_validation_sources/data/` 下修正的第1折结果；`paper/data/case_validation_20261009/` 同名旧表不用于v18。hmu-gpu的SSH服务可达但现有公钥认证被拒绝，未能搜索其目录，是否保存另一套轨迹或原溶剂体系仍为未知；访问限制不会影响本地轨迹复用。详细证明见 `kix_existing_reuse_20261010/reuse_audit.json`、`hmu_audit/connection_audit.json` 和 `v18_content_review.json`。\n'
status_path.write_text(status, encoding='utf8')

md_audit_path = ROOT / 'audit_md.md'
md_audit = md_audit_path.read_text(encoding='utf8')
note = ('> 当前执行策略更新（2026-10-10）：用户明确要求“MD跑一次也够了，先聚焦于论文”。原计划的新12条三重复生产与5项后续分析已在实际运行父队列中取消，不会启动；现有模型训练及CPU对照保留。' +
        '本地历史四状态每状态一条100ns已重新验证全部400段轨迹和质量表并整合进v18。下文三重复制备与运行方案仅为此前计划的保留记录，不能作为当前仍将运行的任务或已完成结果。' +
        '现行策略见 `pipeline_control.json`，运行时取消证明见 `pipeline_md_cancellation_audit_20261010.json`，单轨迹完整核查见 `kix_existing_reuse_20261010/reuse_audit.json`。\n\n')
if not md_audit.startswith('> 当前执行策略更新'):
    md_audit_path.write_text(note + md_audit, encoding='utf8')

print(json.dumps({'docx': str(DOCX), 'sha256': digest, 'pages': 40,
                  'tables': 10, 'visual_and_content_review': 'passed'}, ensure_ascii=False))
