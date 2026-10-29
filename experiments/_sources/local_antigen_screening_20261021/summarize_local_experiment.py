"""Create a readable, source-bound report for the completed local pilot."""
from pathlib import Path
import csv,json,hashlib,datetime
ROOT=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf8')
def main():
    checks={
        'experimental_labels/delivery_manifest.json':'passed',
        'experimental_labels/independent_labels_audit.json':'passed',
        'structures/static_complete_manifest.json':'passed',
        'evaluation_review/final_structure_and_U_review.json':'passed',
        'predictions/fixed_backbone/run_audit.json':'passed',
        'local_runtime/inference_final_review.json':'passed',
        'feature_extraction/feature_manifest.json':'passed',
        'evaluation_review/results_fixed_protocol/evaluation_execution_audit.json':'completed',
        'evaluation_review/independent_evaluation_review.json':'passed',
    }
    for p,prefix in checks.items():assert read(ROOT/p)['status'].startswith(prefix),(p,read(ROOT/p)['status'])
    batch=read(ROOT/'predictions/fixed_backbone/run_audit.json')
    protocol=read(ROOT/'evaluation_review/fixed_analysis_protocol.json')
    assert len(batch['candidates'])==76 and batch['device']=='cpu' and batch['remote_compute'] is False
    assert batch['checkpoint_sha256']==sha(ROOT/'model/epoch_040.pt')==protocol['base_model']['sha256']
    for p in read(ROOT/'evaluation_review/results_fixed_protocol/evaluation_execution_audit.json')['source_files']:
        assert sha(p['path'])==p['sha256'],p
    manuscript=ROOT.parent.parent/'paper_v49_20261019_中文修订_MEL5双动态机制实例_三线表.docx'
    assert sha(manuscript)=='5f79ddf7d657dc9b782c9fb432e373a1e7c170c1c88d019e528fbf73c366c222'
    with (ROOT/'evaluation_review/primary_results_7models.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
    assert len(rows)==7
    lines=['# FlexPPI本地抗原变体排序实验：首轮结果','',
      '首轮计算和独立复算已完成。此次预先固定的用法未显示RMSF与Hχ联合改善抗原候选排序的稳定收益：加入两项后，两个主体系的亲和力变化误差均增加。1G4中Hχ单项出现小幅改善，但在备用模板及逐研究核对中没有稳定复现。这里评估的是候选优选价值，不能据此否定已有残基动态预测任务，也不能据此把模型描述为已经验证的抗原筛选或定向改造工具。','',
      '实验收集了真实SPR测量的1G4–NY-ESO-1面板14种肽与A6–Tax面板16种肽，共30个生物候选、36条研究级均值。主检验使用26种单点变体，每体系13种，按各自6个突变位置留出：某位置的全部候选一起测试，拟合时只使用其他位置。两个母肽仅作差值参照；A6两个双点变体单列附加测试。五个预设模板均用于各自的全部候选，得到76份结构与模型输入；这些模板和多种拟合版本不增加独立实验样本数。','',
      '主标签遵循原始SPR研究的方法，先在同一研究内计算候选相对母肽的亲和力变化，再对重复候选的研究内变化取中位数。pKD=-log10(KD，以M计)，KD下降10倍对应pKD上升1；本报告中的MAE是这种亲和力变化预测的平均绝对误差，越小越好。Pettmann研究的主拟合按原文规则选择，同时完整保留all-constrained与all-fitted敏感性版本，未依据模型成绩选择标签。','',
      '| 方法 | 1G4误差 | 1G4排序相关 | 1G4前三重合 | A6误差 | A6排序相关 | A6前三重合 |',
      '|---|---:|---:|---:|---:|---:|---:|']
    for row in rows:
        z=row['model']=='zero_change'
        def rho(k):return '—' if not row[k] else f'{float(row[k]):.3f}'
        lines.append(f"| {row['method_zh']} | {float(row['nyeso_MAE_delta_pKD']):.4f} | {rho('nyeso_Spearman_rho')} | {'无排序' if z else row['nyeso_top3_overlap']+'/3'} | {float(row['tax_MAE_delta_pKD']):.4f} | {rho('tax_Spearman_rho')} | {'无排序' if z else row['tax_top3_overlap']+'/3'} |")
    lines += ['',
      '联合两项后，1G4误差由0.4324升至0.4488（增加3.8%），A6由0.5985升至0.6474（增加8.2%）。1G4仅加Hχ时降至0.4092（降低5.4%），前三名重合由0/3升至1/3；但这一信号未在2BNR备用模板和Aleksic单研究误差中稳定重现。A6中序列基线自身优于加入静态结构的基线，加入动态特征也未改善主模板误差。五模板乘三种SPR拟合方式的15个体系内联合比较中，13个误差增加；另外两个微小改善均在A6备用1QSE模板，且联合仍不及其中一个单项。因此没有以备用设置替换主检验。','',
      '本轮比较的是一个固定、低自由度的具体方案：序列替换程度与BLOSUM62基线、完整重原子结构的埋藏面积及接触数，再分别加入固定8 Å环境中未替换HLA/TCR残基的RMSF均值变化、Hχ均值变化。排序校准使用预先固定的岭回归，所有尺度和系数仅由训练位置拟合；基础FlexPPI权重没有更新。Hχ特征不包含改变氨基酸类型的肽位点，避免把化学类型差异当成环境效应。独立程序重建630次拟合、2037条分析预测和210组指标，最大预测数值差约4.7×10⁻¹⁵。所有结果和候选排名均保留。','',
      '所有计算在本机CPU完成。五模板的全部候选侧链重原子由PDBFixer一致补齐，原有蛋白重原子（包括全部主链）坐标均保持不变；局部新原子重建使用OpenMM内部优化，没有开展新的蛋白生产MD。正式76份神经网络推理共耗时'+f"{batch['elapsed_seconds']:.1f}"+'秒；固定10例另行复算，输出逐位一致。该实验因此只检验固定实验骨架条件下的序列响应，没有生成突变后的新主链构象。','',
      '使用的模型是本地可取得的历史strict联合模型：2,431条训练记录（ATLAS 1,866、DynaRepo 565），ITrace为0；它与论文当前2,430条记录模型的权重不同。训练清单仅多一条A00IL，目标模板没有精确或既定阈值的直接同源命中，但该记录存在间接家族组件连接，已在开始候选预测前修订协议并披露。因此本报告标为本地初步实验，不冒称当前论文模型的完整家族隔离验证。','',
      '两个面板还主要包含相对母肽的亲和力损失：1G4仅一个候选稍强于母肽，A6全部变体均弱于母肽。因此即便排序表现较好，也主要说明能否区分替换损失，尚不能证明能够设计出亲和力明显提高的新候选。当前结果最直接支持的决策是：先保留并分析本轮完整负结果，不把界面两个动态均值的简单组合用作已验证的候选筛选规则。局部相反变化被均值抵消、真实突变骨架响应和游离态信息缺失等，属于下一轮需单独设定并验证的假设。','',
      '完整结果位于evaluation_review/results_fixed_protocol；逐候选特征位于feature_extraction/candidate_features.csv；预先固定的协议、位置划分、标签来源、权重与结构输入哈希、完整敏感性及独立复算记录均保留在本目录。原论文未改动。','',
      '原始数据来源：[HERMES/PNAS 2025](https://www.pnas.org/doi/10.1073/pnas.2504783122)、[Pettmann等eLife 2021](https://elifesciences.org/articles/67092)，以及上述论文引用并在experimental_labels/raw中归档的Aleksic等2010原始补充表。']
    report=ROOT/'本地抗原筛选实验_首轮结果_20261021.md'
    report.write_text('\n'.join(lines)+'\n',encoding='utf8')
    status={'status':'completed_first_pilot','requested_execution':'local','compute_device':'CPU','remote_compute':False,'biological_candidates':30,'primary_single_mutant_candidates':26,'independent_TCR_systems':2,'templates':5,'input_structures':76,'neural_inference_elapsed_seconds':batch['elapsed_seconds'],'new_base_model_training':False,'local_ridge_fits':630,'new_production_MD':False,'joint_screening_increment_supported':False,'protocol_sha256':sha(ROOT/'evaluation_review/fixed_analysis_protocol.json'),'checkpoint_sha256':batch['checkpoint_sha256'],'same_weights_as_current_manuscript':False,'report':str(report),'report_sha256':sha(report),'checks':{p:sha(ROOT/p) for p in checks},'main_manuscript_unchanged':True,'main_manuscript_sha256':sha(manuscript),'updated_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    write(ROOT/'run_status.json',status)
    write(ROOT/'delivery_status.json',status)
    print(json.dumps(status,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
