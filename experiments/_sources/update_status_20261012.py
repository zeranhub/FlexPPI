from pathlib import Path
ROOT = Path(__file__).resolve().parent
path = ROOT / 'completion_status_zh.md'
paragraphs = path.read_text(encoding='utf8').split('\n\n')
for i, text in enumerate(paragraphs):
    if text.startswith('实际计算在 '):
        paragraphs[i] = ('2026-10-12北京时间早晨已核对 `hit-gpu:/root/flexppi_completion_20261010`：10次单目标训练和10次全链家族隔离训练均已完成，单目标折外预测、配对分析及家族折外评价均有完成审计。'
            '10次单目标训练结束于服务器UTC2026-10-10 23:23，家族训练及评价结束于UTC2026-10-11 15:35。'
            '队列随后在 `evaluate_external_2430` 步骤退出，日志显示两目标累计数组与三通道标签掩码的形状不一致。当前正在修复评估脚本的明确通道对应并恢复剩余步骤，已完成训练不需重跑，冻结权重和输入保持。'
            '具体恢复结果以新的服务器状态、修复审计和最终完成门控为准；不能把已完成20个检查点等同于论文全部补充实验已经完成。')
    if text.startswith('RMSF-net独立CPU进程'):
        paragraphs[i] = ('2026-10-12早晨快照显示CPU预测输出已齐：RMSF-net 3,063/3,063条记录，PLANET full与seq各4,973/4,973条链。'
            'PLANET两项队列均complete且有独立run_audit；RMSF-net全部逐项预测complete但尚无最终独立核验run_audit，剩余主队列应完成它的核验和共同残基评分。'
            '所有这些数量是文件/完成状态核查，正式共同评分仍须使用通过完整性门控的输出。PLANET保留公开float32权重、默认300K和原训练范围；RMSF-net使用PDB模拟密度适配，并非真实冷冻电镜输入的直接复现。')
path.write_text('\n\n'.join(paragraphs), encoding='utf8')
print('Updated durable status; v18 remains the latest reviewed manuscript.')
