from pathlib import Path
import json, hashlib
OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
backup = json.loads((OUT/'local_scientific_backup_audit.json').read_text(encoding='utf8'))
delivery_path = OUT/'delivery_status.json'
delivery = json.loads(delivery_path.read_text(encoding='utf8'))
assert backup['status'] == 'passed' and backup['final_weights_verified'] == 20
doc = Path(delivery['document'])
assert hashlib.sha256(doc.read_bytes()).hexdigest() == delivery['document_sha256']
delivery.update(local_final_weights_and_per_residue_backup='passed',
                scientific_sources_archive=backup['local_archive'],
                scientific_sources_archive_sha256=backup['archive_sha256'])
delivery_path.write_text(json.dumps(delivery, ensure_ascii=False, indent=2), encoding='utf8')
note = '''本版已完成原批注基础上的中文修订及约定补充实验整合。新增十次单目标和十次全链家族模型训练，RMSF-net、PLANET full/seq完整队列预测与共同残基评分；方法、摘要、结果、讨论、结论和三线表均按实际结果修改。外部验证主表与图4统一采用删除633条、保留2,430条训练记录的方案，旧2,431条结果单独保留为历史敏感性，不与新方案混用。

结果并非所有方向都支持原有强结论：联合训练的平均误差收益约0.8%和1.8%，RMSF排序略差于单目标；家族隔离后成绩下降，ITrace只有一个留出家族及来源外推；新严格外部的RMSF幅值误差略高于链内和静态对照；PLANET full整体RMSF MAE与排序优于FlexPPI，但RMSE较高且来源差异明显。稿件如实保留这些结果与比较边界。

KIX复用本地已有四状态各一条100ns轨迹。400段DCD及全部质量分析重新核验，新增12次生产模拟已取消。模型能定位活动环区，但L12六组状态差值方向与单轨迹MD不一致；窗口依赖和缺乏独立重复统计均明确写入正文及补充材料。

交付Word共45页，14张原生可编辑三线表。原6条批注、锚点、图片及原生公式保留。全部变化页完成实际Word渲染视觉检查，其余页正文像素与已核前版一致；独立数值与引用检查通过。原v13文件未覆盖。

复现资产保存在experiment_completion/completed_results_20261012。最终20个权重、四份逐残基压缩源表、配置及完成证书已备份为completed_scientific_sources_20261012.tar.gz（约519MB），本地整个压缩包及59个成员逐文件SHA256均通过核验。新图4的矢量PDF/SVG及300dpi PNG在figures目录，审计和版本对应记录随源文件保留。
'''
(ROOT.parent/'paper_v19_20261012_修改说明.md').write_text(note, encoding='utf8')
print(json.dumps({'delivery':'complete','word':str(doc),'backup_members':59,'sha256':delivery['document_sha256']}, ensure_ascii=False))
