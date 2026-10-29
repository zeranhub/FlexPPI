"""Read frozen v22 tables and scored sources; write only new v23 evidence."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,sys
import numpy as np
import pandas as pd
from docx import Document

ROOT=Path('E:/Project/FlexPPI/paper/revision_v14_20261010')
BASE=ROOT/'experiment_completion'
OUT=BASE/'argument_revision_v23_20261012'
DOC=ROOT/'paper_v22_20261012_中文修订_自由能分析融入正文_三线表.docx'
COMP=BASE/'completed_results_20261012'
SCORED=BASE/'baselines/additional_original_oof'
PRODY=BASE/'baselines/prody_common_evaluation'
NAMES=[('flexppi','完整FlexPPI'),('intrachain','链内FlexPPI'),('planet_full','PLANET full'),('planet_seq','PLANET seq'),('rmsf_net_pdb','RMSF-net PDB'),('anm','ANM'),('gnm','GNM')]

def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):digest.update(block)
    return digest.hexdigest()

def main():
    sys.stdout.reconfigure(encoding='utf-8')
    OUT.mkdir(parents=True,exist_ok=True)
    original_doc_sha=sha(DOC)
    d=Document(DOC)
    # Map native DOCX tables by their actual preceding table caption instead
    # of assuming a fixed array number after a section/paragraph promotion.
    namespace={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    caption=None;native={};index=0
    for child in d.element.body:
        if child.tag.endswith('}p'):
            text=''.join(child.xpath('.//w:t/text()'))
            if text.startswith('表') and '：' in text:caption=text.split('：',1)[0]
        elif child.tag.endswith('}tbl'):
            native[caption]=dict(index=index,cells=[[c.text for c in r.cells] for r in d.tables[index].rows])
            index+=1
    s11,s5=native['表S11'],native['表S5']
    checks=[];frozen={}
    def check(name,condition,detail=None):
        row=dict(name=name,passed=bool(condition))
        if detail is not None:row['detail']=detail
        checks.append(row)
    def csv(path):
        frozen[str(path)]=sha(path)
        return pd.read_csv(path,float_precision='round_trip')
    def js(path):
        frozen[str(path)]=sha(path)
        return json.loads(path.read_text(encoding='utf-8'))
    raw=csv(COMP/'baseline_original_overall_values_verified.csv').set_index('method')
    source=csv(COMP/'baseline_original_source_values_verified.csv')
    scored=csv(SCORED/'per_record_metrics.csv')
    cover=csv(SCORED/'coverage_audit.csv')
    prody=csv(PRODY/'per_record_metrics.csv')
    audit=js(SCORED/'scoring_audit.json')
    check('Current v22 native TableS11 has exactly7 methods4columns',len(s11['cells'])==8 and all(len(r)==4 for r in s11['cells']))
    check('All7baseline predictions fullycovered3063/961321',audit['all_methods_complete'] and audit['records']==3063 and audit['expected_rmsf_residues']==961321 and audit['cohort_freeze_sha256']=='b58795275fa2843a502fdbb68f417f9a2b0b80f76e28d226ad5d58c9bf430a91' and audit['family_fold_sha256'] is None)
    check('Actualcoverage has21441complete records',len(cover)==3063*7 and cover.complete.all() and (cover.expected_valid==cover.covered_valid).all())
    grouped=scored.groupby('method')
    check('Each method3063uniqueIDsall961321RMSF',len(grouped)==7 and all(len(g)==g.sample_id.nunique()==3063 and int(g['count'].sum())==961321 and g.common_fraction.eq(1).all() for _,g in grouped))
    recomputed=grouped[['mae','spearman','z_mae','pearson','top10_precision']].mean()
    for method,_ in NAMES:
        for metric in ['mae','spearman','z_mae','pearson','top10_precision']:
            a,b=recomputed.loc[method,metric],raw.loc[method,metric]
            check('Primaryrecordmacro '+method+'/'+metric,(pd.isna(a) and pd.isna(b)) or abs(a-b)<1e-12)
    expected=[]
    for method,label in NAMES:
        r=raw.loc[method]
        expected.append([label,'不适用' if pd.isna(r.mae) else f'{r.mae:.3f}',f'{r.spearman:.3f}',f'{r.z_mae:.3f}'])
    check('All merged cells exactlynativeS11',s11['cells'][1:]==expected,{'table_cells':s11['cells'][1:],'expected':expected})
    primary=pd.DataFrame(s11['cells'][1:],columns=s11['cells'][0])
    primary.to_csv(OUT/'merged_baseline_table.csv',index=False,encoding='utf-8')
    precise=[]
    for method,label in NAMES:
        r=raw.loc[method]
        precise.append(dict(method=method,label_zh=label,rmsf_mae_angstrom=r.mae,spearman=r.spearman,standardized_mae=r.z_mae,absolute_angstrom=bool(pd.notna(r.mae)),records=3063,common_rmsf_positions=961321,aggregation='record_macro_original_system_OOF'))
    pd.DataFrame(precise).to_csv(OUT/'merged_baseline_table_raw_values.csv',index=False)
    # Unique D evidence is kept separately; no undocumented values are
    # inserted into the compact 7method paper table.
    archived=pd.DataFrame(s5['cells'][1:],columns=s5['cells'][0])
    archived.to_csv(OUT/'tableS5_unique_elastic_network_support.csv',index=False,encoding='utf-8')
    prody_mean=prody.groupby('method')[['top10_precision','pearson']].mean()
    support_names={'完整FlexPPI':'flexppi','链内FlexPPI':'intrachain','ANM':'anm','GNM':'gnm'}
    for row in s5['cells'][1:]:
        method=support_names[row[0]]
        check('NativeS5top10 '+method,row[4]==f'{100*raw.loc[method,"top10_precision"]:.1f}' and abs(prody_mean.loc[method,'top10_precision']-raw.loc[method,'top10_precision'])<1e-12)
        check('NativeS5Pearson '+method,row[2]==f'{raw.loc[method,"pearson"]:.3f}')
    source_recomputed=scored.groupby(['method','dataset'])[['mae','spearman','bias']].mean()
    for row in source.itertuples(index=False):
        for metric in ['mae','spearman','bias']:
            a,b=source_recomputed.loc[(row.method,row.dataset),metric],getattr(row,metric)
            check('Source recordmacro '+row.method+'/'+row.dataset+'/'+metric,(pd.isna(a) and pd.isna(b)) or abs(a-b)<1e-12)
    source.to_csv(OUT/'merged_baseline_source_values.csv',index=False)
    methods='七种方法在3,063条记录的961,321个共同RMSF有效位置上计分，先合并同一记录的各条链，再对记录等权平均。FlexPPI使用原体系五折留出预测；PLANET和RMSF-net保留公开权重。标准化MAE为每条记录的参考和预测分别减去自身均值、除以自身标准差后所得的平均绝对差，不带物理单位。ANM、GNM和RMSF-net的相对分数不拟合成Å幅值，其绝对MAE记为不适用。'
    parameters='ANM与GNM由ProDy 2.6.1对完整输入组装的Cα构建，截断15 Å、弹簧常数1，采用最低20个非零模态，极短结构以可用模态为限；相对分数为平方波动总和的平方根。PLANET full与sequence-only分别逐链推断后合并记录，使用发布模型默认300 K、float32精度，RMSF从nm乘10转为Å。RMSF-net采用PDB-only权重及4 Å分辨率、1.5 Å网格的原子序数加权高斯模拟密度。代码、权重、逐记录结果和共同位置清单随复现材料保留。'
    stratified='来源分层中，ANM与完整FlexPPI的Spearman相关在ATLAS为0.785/0.788，DynaRepo为0.635/0.774，ITrace为0.453/0.938。PLANET full在三来源的MAE/ρ为0.398 Å/0.876、0.704 Å/0.752、1.412 Å/0.732；sequence-only为0.499 Å/0.826、0.782 Å/0.668、2.171 Å/0.719。ITrace中预测减参考的平均偏差分别为+1.381 Å和+2.150 Å；来源分层帮助解释整体记录等权结果。'
    caption='共同RMSF位置上的学习模型与弹性网络基线。MAE单位为Å，ρ为Spearman相关；标准化MAE无量纲，各指标先逐记录计算，再对3,063条记录等权平均。Full/Intra为原五折体系留出预测，发布模型使用原权重；相对分数的绝对幅值误差不适用。'
    note='正文合并表推荐保留原表S11的七种方法与四列，不增加命中率列。原表S5的Pearson和前10%命中率另存tableS5_unique_elastic_network_support.csv；Full65.1%、Intra62.4%、ANM56.7%、GNM57.7%均已由实际逐记录CSV复算，不需要另在论文重复呈现。'
    (OUT/'merged_baseline_text_zh.md').write_text(methods+'\n\n'+parameters+'\n\n'+stratified+'\n\n建议表注：'+caption+'\n\n'+note+'\n',encoding='utf-8')
    frozen[str(DOC)]=original_doc_sha
    check('All sourcefiles frozen before/after',all(sha(path)==digest for path,digest in frozen.items()))
    errors=[c for c in checks if not c['passed']]
    output_files=['merged_baseline_table.csv','merged_baseline_table_raw_values.csv','tableS5_unique_elastic_network_support.csv','merged_baseline_source_values.csv','merged_baseline_text_zh.md']
    report=dict(status='passed' if not errors else 'failed',created_utc=datetime.now(timezone.utc).isoformat(),source_document=str(DOC),source_document_sha256=original_doc_sha,native_tables={'S11':s11,'S5':s5},all7methods4columns_preserved=True,primary_recordmacro_not_foldmean=True,top10_not_added_to_paper_table=True,relative_methods_no_Angstrom_calibration=True,checks=checks,failures=errors,source_sha256=frozen,outputs_sha256={name:sha(OUT/name) for name in output_files},no_existing_file_modified=True,no_Word_edit=True,no_existing_certificate_rerun=True,no_training_or_GPU_or_network=True)
    (OUT/'merged_baseline_table_evidence_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=report['status'],checks=len(checks),failures=errors,table=str(OUT/'merged_baseline_table.csv'),source_doc_sha256=original_doc_sha),ensure_ascii=False,indent=2))
    return 1 if errors else 0
if __name__=='__main__':raise SystemExit(main())
