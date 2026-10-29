"""Read-only reanalysis of archived RMSF profiles, using the user's exact score."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parents[1]/'local_antigen_screening_20261021'
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def csv(df,name):
    df.to_csv(ROOT/name,index=False,encoding='utf-8-sig',float_format='%.17g')
protocol_path=ROOT/'source_only_protocol.json'
protocol=json.loads(protocol_path.read_text(encoding='utf-8'))
for name,expected in protocol['source_files_sha256'].items():
    assert sha(Path(name))==expected,name
manifest=json.loads((OLD/'structures/fixed_backbone_manifest.json').read_text())
sources={str(protocol_path):sha(protocol_path)}
profile_cache={}
def read_peptide(cid,chain,seq):
    p=OLD/'predictions/fixed_backbone'/cid/'prediction.csv'
    sources[str(p)]=sha(p)
    raw=pd.read_csv(p)
    pep=raw.loc[raw.chain.eq(chain)].sort_values('residue_index').reset_index(drop=True)
    assert len(pep)==9,(cid,len(pep))
    assert ''.join(pep.amino_acid)==seq,(cid,seq)
    assert pep.geometry_valid.astype(str).str.lower().eq('true').all(),cid
    assert np.isfinite(pep.rmsf.to_numpy(float)).all(),cid
    assert (pep.rmsf.to_numpy(float)>=0).all(),cid
    return pep
score_rows=[]
site_rows=[]
for c in manifest['candidates']:
    cid=c['candidate_id']
    pid=c['parent_inference_candidate_id']
    pep=read_peptide(cid,c['peptide_chain'],c['peptide_sequence'])
    if pid not in profile_cache:
        profile_cache[pid]=read_peptide(pid,c['peptide_chain'],c['panel_reference_peptide'])
    parent=profile_cache[pid]
    assert list(pep.author_residue_id)==list(parent.author_residue_id),cid
    drop=parent.rmsf.to_numpy(float)-pep.rmsf.to_numpy(float)
    changed=[i+1 for i,(a,b) in enumerate(zip(c['panel_reference_peptide'],c['peptide_sequence'])) if a!=b]
    winners=np.flatnonzero(drop==drop.max())+1
    row={'inference_candidate_id':cid,'candidate_id':c['biological_candidate_id'],
      'panel_id':c['panel_id'],'tcr':c['tcr'],'template':c['template_pdb'],
      'primary_template':bool(c['primary']),'peptide_sequence':c['peptide_sequence'],
      'parent_peptide':c['panel_reference_peptide'],'mutation':c['mutations_vs_panel_reference'],
      'n_mutations':len(changed),'is_parent':len(changed)==0,
      'changed_positions':'|'.join(map(str,changed)),
      'score_max_rmsf_drop_A':float(drop.max()),
      'max_drop_positions':'|'.join(map(str,winners)),
      'max_drop_is_at_mutation':bool(set(winners)&set(changed)),
      'score_definition':'max(parent-mutant RMSF) over all 9 peptide positions; signed'}
    for i in range(9):
        row[f'parent_RMSF_P{i+1}_A']=float(parent.rmsf.iloc[i])
        row[f'mutant_RMSF_P{i+1}_A']=float(pep.rmsf.iloc[i])
        row[f'drop_RMSF_P{i+1}_A']=float(drop[i])
        site_rows.append({'inference_candidate_id':cid,'candidate_id':c['biological_candidate_id'],
           'panel_id':c['panel_id'],'template':c['template_pdb'],'primary_template':bool(c['primary']),
           'peptide_position':i+1,'author_residue_id':pep.author_residue_id.iloc[i],
           'parent_aa':c['panel_reference_peptide'][i],'mutant_aa':c['peptide_sequence'][i],
           'is_mutation':i+1 in changed,'parent_rmsf_A':float(parent.rmsf.iloc[i]),
           'mutant_rmsf_A':float(pep.rmsf.iloc[i]),'rmsf_drop_A':float(drop[i]),
           'is_maximum':i+1 in winners})
    score_rows.append(row)
scores=pd.DataFrame(score_rows)
assert len(scores)==76 and scores.candidate_id.nunique()==30
assert scores[scores.primary_template].n_mutations.value_counts().to_dict()=={1:26,0:2,2:2}
assert np.all(scores.loc[scores.is_parent,'score_max_rmsf_drop_A']==0)
csv(scores,'candidate_scores_all_templates.csv')
csv(pd.DataFrame(site_rows),'per_site_rmsf_profiles_and_drops.csv')

all_joined=[]
all_study=[]
metric_rows=[]
def metric(df,label,level,panel,template,study='candidate_median'):
    assert df.candidate_id.nunique()==len(df),(level,panel,template,study)
    x=df.score_max_rmsf_drop_A.to_numpy(float)
    y=df.delta_pKD_vs_parent.to_numpy(float)
    rho=float(spearmanr(x,y).statistic) if len(np.unique(x))>1 and len(np.unique(y))>1 else None
    n_pair=0; concordant=0; discordant=0; tied_score=0
    for i in range(len(x)):
        for j in range(i):
            if y[i]==y[j]:
                continue
            n_pair+=1
            sign=np.sign(x[i]-x[j])*np.sign(y[i]-y[j])
            concordant+=int(sign>0);discordant+=int(sign<0);tied_score+=int(sign==0)
    metric_rows.append({'SPR_variant':label,'analysis_level':level,'panel_id':panel,
      'template':template,'study':study,'n_unique_single_mutants':len(df),'spearman_rho':rho,
      'n_experimentally_ordered_pairs':n_pair,'concordant_pairs':concordant,
      'discordant_pairs':discordant,'score_tied_pairs':tied_score,
      'pair_order_agreement_with_half_credit_ties':(concordant+.5*tied_score)/n_pair if n_pair else None,
      'n_stronger_than_parent':int((y>0).sum()),'n_weaker_than_parent':int((y<0).sum()),
      'n_equal_parent':int((y==0).sum()),'score_direction':'higher drop -> higher affinity',
      'independence':'Two TCR series; templates/fits/studies are not extra biological systems'})

for variant in protocol['inherited_fixed_choices']['SPR_versions']:
    lp=OLD/f'experimental_labels/SPR_labels_candidates_{variant}.csv'
    op=OLD/f'experimental_labels/SPR_labels_observations_{variant}.csv'
    labels=pd.read_csv(lp)
    obs=pd.read_csv(op)
    assert labels.candidate_id.nunique()==30 and len(labels)==30
    assert set(labels.censoring)=={'none_reported'},set(labels.censoring)
    assert set(obs.censoring)=={'none_reported'}
    assert set(obs.KD_relation)=={'='}
    assert set(obs.assay)=={'SPR'} and set(obs.temperature_C)=={37}
    assert obs.parent_same_study.astype(str).str.lower().eq('true').all()
    assert obs.parent_same_protocol_temperature.astype(str).str.lower().eq('true').all()
    assert np.allclose(obs.pKD_M,6-np.log10(obs.KD_uM),atol=1e-14)
    assert np.allclose(obs.delta_pKD_vs_parent,obs.pKD_M-obs.parent_pKD_M,atol=1e-14)
    for cid,g in obs.groupby('candidate_id'):
        assert np.isclose(g.delta_pKD_vs_parent.median(),labels.loc[labels.candidate_id.eq(cid),'delta_pKD_vs_parent'].iloc[0],atol=1e-14)
    joined=scores.merge(labels[['candidate_id','delta_pKD_vs_parent','label_variant','n_study_observations',
       'study_ids','study_delta_pKD_json','min_study_delta_pKD','max_study_delta_pKD','censoring',
       'assay','temperature_C','parent_same_batch']],on='candidate_id',validate='many_to_one')
    joined['SPR_variant']=variant
    joined['score_rank_desc_within_template']=joined.loc[joined.n_mutations.eq(1)].groupby(['panel_id','template']).score_max_rmsf_drop_A.rank(ascending=False,method='min')
    joined['affinity_rank_desc_within_template']=joined.loc[joined.n_mutations.eq(1)].groupby(['panel_id','template']).delta_pKD_vs_parent.rank(ascending=False,method='min')
    all_joined.append(joined)
    primary=joined[joined.n_mutations.eq(1)&joined.primary_template]
    assert len(primary)==26
    metric(primary,variant,'pooled_primary_templates_descriptive','two_panels','2BNQ+1AO7')
    for (panel,template),g in joined[joined.n_mutations.eq(1)].groupby(['panel_id','template']):
        assert len(g)==13
        metric(g,variant,'candidate_within_panel',panel,template)
    # Actual source-level KD values are retained, rather than fabricate a pooled absolute KD.
    obs_cols=['observation_id','candidate_id','study_id','KD_uM','pKD_M','parent_KD_uM','parent_pKD_M',
       'delta_pKD_vs_parent','KD_relation','censoring','assay','temperature_C','error_value','error_type',
       'n','primary_doi','primary_source_url','primary_locator','selected_fit','fit_selection_rule']
    sj=scores.merge(obs[obs_cols],on='candidate_id',validate='many_to_many')
    sj['SPR_variant']=variant
    all_study.append(sj)
    for (panel,template,study),g in sj[sj.n_mutations.eq(1)].groupby(['panel_id','template','study_id']):
        metric(g,variant,'within_original_study',panel,template,study)

joined_all=pd.concat(all_joined,ignore_index=True)
study_all=pd.concat(all_study,ignore_index=True)
metrics=pd.DataFrame(metric_rows)
csv(joined_all,'scores_and_affinity_all_templates_all_SPR.csv')
csv(study_all,'scores_and_original_study_KD_all_versions.csv')
csv(metrics,'spearman_and_pair_order_all_sensitivities.csv')
primary=joined_all[joined_all.SPR_variant.eq('reported_hybrid')&joined_all.primary_template&joined_all.n_mutations.eq(1)].sort_values(['panel_id','score_rank_desc_within_template','candidate_id'])
csv(primary,'primary_26_single_mutants_ranked.csv')
csv(joined_all[joined_all.n_mutations.eq(2)],'two_double_mutants_separate.csv')
primary_metrics=metrics[metrics.SPR_variant.eq('reported_hybrid')&((metrics.analysis_level.eq('pooled_primary_templates_descriptive'))|((metrics.analysis_level.eq('candidate_within_panel'))&(metrics.template.isin(['2BNQ','1AO7']))))]

lines=['# 肽任一位点最大RMSF下降与直接SPR亲和力：探索性复核','',
 '按照用户指定的方向，对每个候选取肽全部9个位点中“母肽预测RMSF减去候选预测RMSF”的最大值，数值越大越优先。没有训练、重新推理、使用Hχ或依据结果选择阈值。主模板及SPR拟合版本沿用上一轮在结果前固定的2BNQ、1AO7和reported_hybrid。上一轮标签与结果已经知道，因此本轮是新假设的探索性复核，并非新的盲测。','',
 '主比较只使用26个单突变，每个TCR各13个；母肽只作差值参照，两个A6双突变单列。亲和力为同研究、同SPR条件37℃下的ΔpKD，越大代表结合越强；同一候选跨研究时先各自减去该研究母肽，再取中位数。所有原始KD与pKD仍按研究单独保存，未把跨研究的绝对KD合成一个假值。','',
 '| 范围 | 单突变数 | Spearman相关 | 候选两两顺序一致率 |','|---|---:|---:|---:|']
for r in primary_metrics.to_dict('records'):
    lines.append(f"| {r['panel_id']}（{r['template']}） | {r['n_unique_single_mutants']} | {r['spearman_rho']:.6f} | {r['pair_order_agreement_with_half_credit_ties']:.4f} |")
lines+=['','正相关支持“大幅下降排在前面”的方向，负相关表示相反；合并值只是这两个有限面板的描述，不能增加独立TCR体系数。所有5个模板、3种SPR版本及逐研究比较都在spearman_and_pair_order_all_sensitivities.csv中完整保留，备用设置不替代主设置。','',
 '这不是完整突变库：1G4主标签只有1个变体略强于母肽，A6的变体均弱于母肽，因此不能从这些数据估计实际筛选的命中率，也不能把排在前面解释为已经设计出高亲和力抗原。模型输入保留母本骨架；变化是固定骨架条件下的序列响应，没有显式生成突变后新构象。使用历史2431条训练版联合模型，并非论文当前2430条训练版权重；相关性也不证明RMSF改变导致亲和力改变。','',
 '完整逐候选文件包含九个位点的母肽RMSF、候选RMSF、降幅、最大降幅发生位置、全部研究亲和力差值和排序，不删除负值或不利候选。']
(ROOT/'直接亲和力复核_结果.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
audit={'status':'complete','completed_at_utc':datetime.now(timezone.utc).isoformat(),
 'n_inference_profiles':len(scores),'n_unique_candidates':scores.candidate_id.nunique(),
 'n_per_site_rows':len(site_rows),'n_primary_singles':len(primary),'n_metrics':len(metrics),
 'all_sources_uncensored_SPR_37C':True,'all_pKD_and_parent_deltas_independently_recomputed':True,
 'all_candidate_median_deltas_independently_recomputed':True,'no_inference_or_training_performed':True,
 'source_protocol_sha256':sha(protocol_path),'prediction_csv_hashes':sources,
 'output_sha256':{p.name:sha(p) for p in ROOT.iterdir() if p.is_file() and p.name not in ['analysis_audit.json']}}
(ROOT/'analysis_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(primary_metrics.to_json(orient='records',force_ascii=False,indent=2))
print('PRIMARY TOP 5 EACH:')
print(primary.groupby('panel_id').head(5)[['panel_id','mutation','score_max_rmsf_drop_A','max_drop_positions','delta_pKD_vs_parent','affinity_rank_desc_within_template']].to_string(index=False))
print('METRICS:',len(metrics),'ROWS:',len(joined_all),'STUDY_ROWS:',len(study_all))
