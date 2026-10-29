"""Pre-result pilot amendment: availability changes weights, never results."""
from pathlib import Path
from datetime import datetime,timezone
import copy,hashlib,json
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parents[1]
REF=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
read=lambda p:json.loads(p.read_text(encoding='utf8'))
write=lambda p,x:p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
now=datetime.now(timezone.utc).isoformat()
historical=read(ROOT/'historical_strict_2431_checkpoint_audit.json')
assert historical['status']=='passed_read_only_historical_checkpoint_audit'
archive=ROOT/'before_2431_pilot_amendment';archive.mkdir(exist_ok=True)
for name in ['fixed_analysis_protocol.json','fixed_analysis_protocol_zh.md','protocol_freeze_binding.json','protocol_gate_audit.json',
             'training_overlap_audit.json','training_overlap_audit_zh.md','raw_template_overlap_audit.json']:
    p=ROOT/name;target=archive/name
    assert not target.exists(),'Do not silently re-run or overwrite the pre-amendment record.'
    target.write_bytes(p.read_bytes())
oldprotocol=read(ROOT/'fixed_analysis_protocol.json')
new=copy.deepcopy(oldprotocol)
oldbase=new['base_model']
model=copy.deepcopy(oldbase)
model.update(sha256=historical['actual_checkpoint']['sha256'],training_manifest_sha256=historical['strict_train_manifest']['sha256'],
    records=2431,DynaRepo=565,ATLAS=1866,ITrace=0,
    name='historical strict_full fold0 no-ITrace joint model',
    local_original_path=historical['actual_checkpoint']['path'],
    local_run_path=str(ROOT.parent/'model/epoch_040.pt'),
    training_manifest_path=historical['strict_train_manifest']['path'],
    same_weights_as_current_manuscript_2430=False,
    original_2430_weight_not_available_locally=True,
    direct_template_sequence_overlap_detected=False,
    fully_transitive_family_disjoint=False,
    transitive_component_overlap_record_count=1,
    transitive_component_overlap_records=['A00IL'],
    scientific_role='Availability-limited local pilot. It does not reproduce the current manuscript 2430 checkpoint and must be labelled as the historical 2431 model in all results.')
new['protocol_id']='local_antigen_retrospective_v1_1_local2431pilot_20261021'
new['status']='pilot_model_amended_before_new_panel_model_outputs'
new['base_model']=model
new['scope']='Local retrospective candidate-ranking pilot with the available historical2431 no-ITrace joint model; no base-model training or MD and no manuscript change.'
new['amendments'].append(dict(timestamp_utc=now,kind='pre_result_checkpoint_availability_substitution',
    old_model=oldbase,new_model=model,
    prior_protocol_sha256=sha(archive/'fixed_analysis_protocol.json'),
    prior_binding_sha256=sha(archive/'protocol_freeze_binding.json'),
    no_new_panel_model_outputs_seen=True,no_feature_split_or_metric_rules_changed=True,
    reason='Human-authorized local pilot: target2430 weights are absent from searched local worktrees; remote read-only retrieval timed out/authentication failed. A verified historical2431 no-ITrace joint checkpoint exists locally. No new remote calculations or result-dependent checkpoint selection occurred.',
    audit_path=str(ROOT/'historical_strict_2431_checkpoint_audit.json'),audit_sha256=sha(ROOT/'historical_strict_2431_checkpoint_audit.json'),
    explicit_limitation='The only extra training record A00IL belongs to the newer transitive heldout component via an indirect 30.3%-identity MMseqs edge; no direct target-template threshold overlap was found. Historical strict is its run name, not a claim of complete new-family disjointness.'))
new['leakage_gates'][0]='Downloaded/copied local checkpoint bytes must equal b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f before pilot inference.'
new['leakage_gates'][1]='Actual historical train manifest SHA must equal e0c84c566dea09b436dd7d07e9e6fd8ce1a8cf0c72be20931179e6689c8b44bd, with2431 records and0 ITrace; disclose the single indirect transitive-family member A00IL.'
new['interpretation']='A successful pilot would support feasibility of candidate-ranking features from this historical2431 model in two known TCR/HLA contexts. It would not reproduce the manuscript2430 model, establish fully family-disjoint validation, novel-antigen discovery, T-cell activation, safety, broad TCR generalization, binding free energy, or causal Hchi-to-RMSF transmission.'
write(ROOT/'fixed_analysis_protocol.json',new)

old=read(archive/'training_overlap_audit.json');updated=copy.deepcopy(old)
updated['status']='complete_local_historical2431_pilot_audit_with_disclosed_transitive_overlap'
updated['model_amendment_at_utc']=now
updated['previous_2430_audit_path']=str(archive/'training_overlap_audit.json')
updated['previous_2430_audit_sha256']=sha(archive/'training_overlap_audit.json')
updated['checkpoint_provenance']=dict(local_checkpoint_path=historical['actual_checkpoint']['path'],
    local_checkpoint_sha256=historical['actual_checkpoint']['sha256'],local_checkpoint_byte_hash_verified=True,
    local_run_path=str(ROOT.parent/'model/epoch_040.pt'),sequence_backend='learned',epoch=40,seed=7,
    same_as_current_manuscript_checkpoint=False)
updated['training_manifest']=historical['strict_train_manifest']
for r in updated['record_membership']:
    r['in_reference_2430_fold0_test']=r.pop('in_fold0_test')
    r.pop('component_members_present_in_actual_fold0_train')
    r['component_members_present_in_actual_pilot_train']=1
    r['transitive_training_component_members']=['A00IL']
updated['transitive_family_audit']=dict(fully_disjoint=False,component_id='family_9cd89939c749e38c',
    actual_pilot_training_records_in_component=['A00IL'],
    direct_target_chain_training_hits=0,
    connecting_record='A00C3',connecting_pdb='2V9T',connecting_chain='A',connecting_chain_length=107,
    added_record='A00IL',added_pdb='4WV1',added_chain='C',added_chain_length=98,
    original_MMseqs_edge_identity=.303,added_integer_span_coverage=70/98,neighbour_integer_span_coverage=77/107,
    interpretation='One training member is indirectly connected to the target component by the newer conservative family graph. It is not a direct sequence match to1G4/A6, but this pilot is not completely transitive-family-disjoint.')
updated['historical_model_audit_path']=str(ROOT/'historical_strict_2431_checkpoint_audit.json')
updated['historical_model_audit_sha256']=sha(ROOT/'historical_strict_2431_checkpoint_audit.json')
updated['direct_template_and_candidate_audits_inherit_2430_plus_explicit_added_record_check']=True
updated['conclusions_not_supported'].append('Historical strict2431 has one indirect member of the newer target family component and cannot be called fully family-disjoint or the same current-manuscript2430 weights.')
write(ROOT/'training_overlap_audit.json',updated)

raw=read(archive/'raw_template_overlap_audit.json')
raw['audited_base_2430_manifest_sha256']=old['training_manifest']['sha256']
raw['current_pilot_manifest_sha256']=historical['strict_train_manifest']['sha256']
raw['added_record_check']=dict(record='A00IL',exact_matches=0,local_30pct_70pct_both_matches=0,
    source_audit=str(ROOT/'historical_strict_2431_checkpoint_audit.json'))
raw['complete_transitive_family_separation']=False
raw['scope']='No direct exact or local-threshold target-template overlap, not complete transitive-family separation.'
write(ROOT/'raw_template_overlap_audit.json',raw)

binding=read(archive/'protocol_freeze_binding.json')
binding['protocol_id']=new['protocol_id']
binding['initial_freeze_at_utc']=binding['frozen_at_utc']
binding['pilot_checkpoint_amendment_at_utc']=now
binding['current_model_sha256']=model['sha256']
binding['files']={p.name:sha(p) for p in [ROOT/'fixed_analysis_protocol.json',ROOT/'candidate_splits.json']}
binding['amendments']=new['amendments']
binding['previous_snapshot_path']=str(archive)
binding['new_model_panel_outputs_unseen']=True
binding['no_feature_split_or_metric_rules_changed']=True
write(ROOT/'protocol_freeze_binding.json',binding)

md=(archive/'fixed_analysis_protocol_zh.md').read_text(encoding='utf8')
md=md.replace('模型固定为联合训练的Full family fold0、seed7、第40轮，权重SHA为ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5，不更新基础权重。',
    '本地pilot改用已核验的历史strict_full fold0联合模型，seed7、第40轮，权重SHA为b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f，不更新基础权重。该模型训练2,431条记录，与现稿2,430条模型仅相差A00IL，但权重不同，不能称为复现现稿模型。')
insert='\n本协议在任何新面板模型预测之前作了一次可追溯修订：目标2,430条模型未在本地找到，远程只读获取未成功，因此本地试验改用确实存在的2,431条历史联合模型；标签、候选分组、特征、回归惩罚及成绩规则完全不变。额外的A00IL是4WV1链C的98位单链记录，虽然未与本轮模板直接达到既定同源阈值，但通过一条30.3%相同率的间接连接属于新版目标家族组件。因此，这次pilot既不能冒称现稿同一权重，也不能冒称完整传递家族隔离验证。原协议、修订前文件与修订时间保存在before_2431_pilot_amendment及JSON记录中。\n'
md=md.replace('\n\n本实验要检验',insert+'\n本实验要检验',1)
(ROOT/'fixed_analysis_protocol_zh.md').write_text(md,encoding='utf8')

md=(archive/'training_overlap_audit_zh.md').read_text(encoding='utf8')
md=md.replace('训练重叠与可解释范围审查','训练重叠与可解释范围审查：历史2431本地pilot',1)
md=md.replace('\n\n已核对',
    '\n\n本轮最终执行模型已在结果产生前改为本地历史strict_full权重，SHA为b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f。直接读取CPU checkpoint配置并核对实际训练manifest，清单SHA为e0c84c566dea09b436dd7d07e9e6fd8ce1a8cf0c72be20931179e6689c8b44bd，共2,431条记录：ATLAS 1,866条、DynaRepo 565条、ITrace零条，对应3,013条链与2,589种唯一序列。公共2,430条的输入路径、标签SHA和标签协议完全相同，唯一额外记录为A00IL（4WV1链C、98位）。下面原有2,430条模型审查作为明确标注的基础证据保留，不是本轮实际模型的名称或权重证明。\n\n已核对',1)
md+='\n\n补核新增A00IL后，精确目标链及候选肽重叠仍为零；与六个原始模板及17个片段的直接30%相同率、双侧70%覆盖命中也为零。不过，A00IL通过补充MMseqs搜索中与A00C3（2V9T链A、107位）的一条边进入新版大组件：相同率0.303，跨度覆盖分别70/98与77/107。因而历史模型确实有1条训练记录属于目标的传递家族组件，不能继续写成完整组件零重叠。这一限制不会被“未训练ITrace”自动消除。本次可作为使用历史权重的本地可行性pilot开展，其结果不能冒称现稿模型的精确复现、完全家族隔离验证或随机重复。新旧审查、完整配置、差异及修订原因均已归档。\n'
(ROOT/'training_overlap_audit_zh.md').write_text(md,encoding='utf8')
print(json.dumps(dict(status='amended_before_panel_predictions',protocol_id=new['protocol_id'],new_checkpoint_sha256=model['sha256'],
    new_manifest_sha256=model['training_manifest_sha256'],transitive_component_training_records=1,
    exact_or_direct_threshold_target_overlap=0,all_features_splits_metrics_unchanged=True,amended_at_utc=now),ensure_ascii=False))
