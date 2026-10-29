"""Read-only manuscript/source audit; writes only new v23 QA JSON and Markdown."""
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import hashlib, json, re, zipfile
import numpy as np
import pandas as pd
from docx import Document
from lxml import etree

ROOT = Path('E:/Project/FlexPPI/paper/revision_v14_20261010')
BASE = ROOT / 'experiment_completion'
V23 = BASE / 'argument_revision_v23_20261012'
V22 = BASE / 'energy_landscape_v22_maintext_20261012'
V21 = BASE / 'energy_landscape_v21_20261012'
QA = V23 / 'manuscript_QA'
DOC = ROOT / 'paper_v23_20261012_中文修订_论证精简与补充整合_三线表.docx'
OLD = ROOT / 'paper_v22_20261012_中文修订_自由能分析融入正文_三线表.docx'
TXT = QA / 'v23_content.txt'
sources, checks, issues = {}, [], []

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def check(name, passed, detail=None):
    item = {'name': name, 'passed': bool(passed)}
    if detail is not None: item['detail'] = detail
    checks.append(item)

def read_json(path):
    sources[str(path)] = sha(path)
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def read_csv(path, **kwargs):
    sources[str(path)] = sha(path)
    return pd.read_csv(path, float_precision='round_trip', **kwargs)

def semantic(node):
    return (node.tag, tuple(sorted(node.attrib.items())), node.text or '',
            tuple(semantic(c) for c in node))

def main():
    import sys
    sys.stdout.reconfigure(encoding='utf-8')
    doc_hash, text_hash, old_hash = sha(DOC), sha(TXT), sha(OLD)
    sources[str(OLD)] = old_hash
    doc, old = Document(DOC), Document(OLD)
    p = [''.join(n.text or '' for n in x._p.xpath('.//w:t | .//m:t')).strip()
         for x in doc.paragraphs]
    parsed = {int(i): t for i, t in re.findall(r'^P(\d+): (.*)$', TXT.read_text(encoding='utf-8'), re.M)}
    check('Exact current DOCX-to-text including native equations',
          all(p[i] == t for i, t in parsed.items()) and len(parsed) == sum(bool(x) for x in p))
    integration = read_json(QA / 'integration_audit.json')
    check('Current candidate hash matches integration audit', doc_hash == integration['document_sha256'])
    check('Frozen v22 DOCX unchanged', old_hash == '315ca2108ab9dbdfae83774f4d9fe1556087d66f8ea09697afa9c9eb67b67bdd')
    tables = [[[c.text for c in row.cells] for row in t.rows] for t in doc.tables]
    old_tables = [[[c.text for c in row.cells] for row in t.rows] for t in old.tables]
    mapping = [('表1',0,0),('表2',1,1),('表S1',2,3),('表S2',3,4),
               ('表S3',4,11),('表S4',5,13),('表S5',6,12),('表S6',7,8),
               ('表S7',8,2),('表S8',9,9),('表S9',10,14)]
    check('Eleven native tables and all expected source tables present', len(tables)==11 and len(old_tables)==15)
    table_audit = []
    for label, current, prior in mapping:
        equal = tables[current] == old_tables[prior]
        check(label + ' every native cell exactly equals frozen v22 source', equal)
        table_audit.append({'table':label,'current_table_index':current,'source_v22_table_index':prior,
                            'rows':len(tables[current]),'columns':len(tables[current][0]),
                            'cells':sum(map(len,tables[current])),'exact':equal})
    check('Original six-state delta table removed', not any(t == old_tables[10] for t in tables))
    check('Ordinary external and posthoc calibration native tables removed',
          not any(t in [old_tables[5],old_tables[6]] for t in tables))
    check('Merged baseline exactly matches separately audited merge CSV',
          tables[5] == [list(pd.read_csv(V23/'merged_baseline_table.csv',dtype=str,keep_default_na=False).columns)]
          + pd.read_csv(V23/'merged_baseline_table.csv',dtype=str,keep_default_na=False).values.tolist())
    sources[str(V23/'merged_baseline_table.csv')] = sha(V23/'merged_baseline_table.csv')
    merge = read_json(V23/'merged_baseline_table_evidence_audit.json')
    check('Independent merged baseline evidence certificate passed', merge.get('status') in ['passed','PASS','complete'])

    s = read_csv(BASE/'completed_results_20261012/summary/single_target_overall_record_macro.csv')
    rows=[]
    for metric, variant, target, name in [('rmsf','rmsf_only','RMSF','单目标'),('rmsf','full','RMSF','联合目标'),
                                         ('chi_entropy','chi_entropy_only','Hχ','单目标'),('chi_entropy','full','Hχ','联合目标')]:
        r=s[s.metric.eq(metric)&s.variant.eq(variant)].iloc[0]
        rows.append([target,name,f'{r.mae:.3f}',f'{r.rmse:.3f}',f'{r.spearman:.3f}'])
        check('S3 same3063 and exact target residue denominator '+metric+'/'+variant,
              r.records==3063 and r.valid_residues==({'rmsf':961321,'chi_entropy':824988}[metric]))
    check('S3 values independently rounded from completed single/dual summary',tables[4][1:]==rows)
    gains=[]
    for metric, variant in [('rmsf','rmsf_only'),('chi_entropy','chi_entropy_only')]:
        a=s[s.metric.eq(metric)&s.variant.eq(variant)].iloc[0].mae
        b=s[s.metric.eq(metric)&s.variant.eq('full')].iloc[0].mae
        gains.append(format(100*(1-b/a),'.1f'))
    check('Small joint MAE gains and RMSF single-target advantage honestly retained',gains==['0.8','1.8'] and 'RMSE和排序略好' in p[171])

    family = read_csv(BASE/'completed_results_20261012/summary/family_overall_record_macro.csv')
    fs = read_csv(BASE/'completed_results_20261012/summary/family_source_record_macro.csv')
    rows=[]
    for dataset, label in [('all','整体'),('atlas','ATLAS'),('dynarepo','DynaRepo'),('itrace','ITrace')]:
        for variant, name in [('full','完整'),('intrachain','链内')]:
            selected=family if dataset=='all' else fs[fs.dataset.eq(dataset)]
            r=selected[selected.metric.eq('rmsf')&selected.variant.eq(variant)].iloc[0]
            h=selected[selected.metric.eq('chi_entropy')&selected.variant.eq(variant)].iloc[0]
            rows.append([label,name,f'{r.mae:.3f}',f'{r.spearman:.3f}',f'{h.mae:.3f}',f'{h.spearman:.3f}'])
    check('S5 eight source/overall rows independently rounded from family summaries',tables[6][1:]==rows)
    check('Family ITrace explicitly244records onefold and no fabricated SD',
          '244' in p[177] and '第0折' in p[177] and '仅第0折' in p[178]
          and not any('±' in cell for row in tables[6] for cell in row))
    check('Family gains accurately rounded',
          [format(100*(1-family[family.metric.eq(m)&family.variant.eq('full')].iloc[0].mae /
                       family[family.metric.eq(m)&family.variant.eq('intrachain')].iloc[0].mae),'.1f')
           for m in ['rmsf','chi_entropy']]==['7.5','3.8'])
    external = read_csv(BASE/'completed_results_20261012/summary/external_2430_record_macro.csv')
    extrows=[]
    for variant, name in [('full','完整模型'),('intrachain','仅链内模型'),('static_sequence_complex','静态特征模型'),('aa_only','氨基酸类型均值')]:
        r=external[external.variant.eq(variant)&external.metric.eq('rmsf')].iloc[0]
        h=external[external.variant.eq(variant)&external.metric.eq('chi_entropy')].iloc[0]
        extrows.append([name,f'{r.mae:.3f}',f'{r.spearman:.3f}',f'{h.mae:.3f}',f'{h.spearman:.3f}'])
        check('External32 common-target denominators '+variant,r.records==h.records==32 and r.valid_residues==18748 and h.valid_residues==16265)
    check('MainTable2 current2430 exact rounded summary values',tables[1][1:5]==extrows)
    check('MainTable2 preserves public BackFlip2 reference row',tables[1][5]==['BackFlip2发布模型','0.749','0.593','无','无'])
    check('External purge633/train2430/familyfold0 protocol explicit',
          all(x in p[43] for x in ['633','2,430','第0折']) and all(x in p[87] for x in ['633','2,430']))
    check('ExternalRMSF slight disadvantage rather than false two-output win explicitly stated',
          all(x in p[88] for x in ['0.474','0.471','主要幅值收益来自Hχ']) and not re.search(r'两项.*(?:均低于|均优于)',p[88]))

    state=read_csv(BASE/'kix_existing_reuse_20261010/summary/state_summary.csv')
    delta=read_csv(BASE/'kix_existing_reuse_20261010/summary/l12_delta_summary.csv')
    rows=[]
    for variant, label in [('apo','游离KIX'),('myb','cMyb'),('mll','MLL'),('both','双伙伴')]:
        r=state[state.state.eq(variant)].iloc[0]
        rows.append([label,f'{r.RMSF_MAE_angstrom:.3f}',f'{r.RMSF_spearman:.3f}',
                     str(r.top20_overlap)+'/18',f'{r.L12_MD_mean_RMSF_angstrom:.3f}',f'{r.L12_fold1_mean_RMSF_angstrom:.3f}'])
        check('KIX one independent100ns/87positions '+variant,r.independent_trajectories==1 and r.production_ns==100 and r.KIX_residues==87 and r.top20_nominated==18 and r.L12_residues==8)
    check('S8 all currentfold1 state values independently rounded',tables[9][1:]==rows)
    check('All six KIX delta directions opposite in real source',len(delta)==6 and (delta.MD_delta_RMSF_angstrom*delta.fold1_delta_RMSF_angstrom<0).all())
    for contrast, expected in [('both-apo',['+1.017','−1.589']),('mll-myb',['+0.014','−1.263'])]:
        r=delta[delta.contrast.eq(contrast)].iloc[0]
        values=[format(r.MD_delta_RMSF_angstrom,'+.3f'),format(r.fold1_delta_RMSF_angstrom,'+.3f').replace('-','−')]
        check('F retained accurate delta '+contrast,values==expected and all(x in p[188] for x in expected))
    check('KIX response mismatch retained without fictitious replicates', '六个' in p[188] and '相反' in p[188] and '一条100 ns' in p[186])

    classes=read_csv(V21/'oof_analysis/reference_class_record_macro.csv')
    c=classes[classes.dataset.eq('all')].set_index('reference_class')
    names=[('low_RMSF_low_Hchi','低RMSF／低Hχ'),('low_RMSF_high_Hchi','低RMSF／高Hχ'),
           ('high_RMSF_low_Hchi','高RMSF／低Hχ'),('high_RMSF_high_Hchi','高RMSF／高Hχ'),('unclassified','未分类')]
    rows=[]
    for key,label in names:
        r=c.loc[key]
        rows.append([label,format(int(r.common_residues),','),f'{100*r.reference_fraction_record_macro:.3f}',
                     f'{r.rmsf_mae_record_macro:.3f}',f'{r.hchi_mae_record_macro:.3f}',
                     f'{100*r.rank_class_recall_record_macro:.1f}' if key!='unclassified' else '无'])
    check('S9 allclass counts/macrofractions/errors/recalls independently rounded',tables[10][1:]==rows)
    check('Class counts preserve all824988 common positions',int(c.common_residues.sum())==824988)
    check('Main3.3 lowdiscordant recall retained',all(x in p[82] for x in ['15.8%','16.5%','24,842','22,512']))
    check('Class proportions and conditional metric macrodefinition retained',
          '未分类位置仍计入比例分母' in p[190] and '存在该类的记录' in p[191])

    d=read_csv(V21/'trajectory_analysis/distribution_record_macro.csv')
    n=read_csv(V21/'trajectory_analysis/nonpro_chi1_record_macro.csv')
    d=d[d.replica_mode.eq('replica1')];n=n[n.replica_mode.eq('replica1')]
    phase=[]
    for shift in [0,30,60,90]:
        x=d[d.shift_degrees.eq(shift)];y=n[n.shift_degrees.eq(shift)]
        check('S7 all83 records and sourcecounts at phase '+str(shift), len(x)==len(y)==83 and x.groupby('dataset').size().to_dict()=={'atlas':30,'dynarepo':23,'itrace':30})
        phase.append([str(shift)+'°',f'{x.h_joint.mean():.3f} ± {x.h_joint.std(ddof=1):.3f}',f'{y.h_chi1.mean():.3f}',
                      f'{100*y.chi1_boundary_10deg_fraction.mean():.2f}',f'{y.apparent_chi1_changes.mean():.1f}',f'{x.half_joint_jsd.mean():.3f}'])
    check('S7 all phase rows independently recomputed joint and nonPro chi1',tables[8][1:]==phase)
    check('SupplementE pooled entropy versus windowmean numbers unchanged',
          [f'{d[d.shift_degrees.eq(k)].full_minus_window_h_weighted.mean():.3f}' for k in [0,60]]==['0.135','0.170'] and all(x in p[184] for x in ['0.135','0.170','帧数加权']))
    fa=read_json(V23/'figures/fig8_joint_partition_CA.audit.json')
    for path,digest in fa['input_sha256'].items():
        check('Fig8 frozen input hash '+Path(path).name,sha(path)==digest)
        sources[path]=sha(path)
    for key,meta in fa['outputs'].items():
        check('Fig8 output hash '+key,sha(meta['path'])==meta['sha256'])
        sources[meta['path']]=sha(meta['path'])
    ca=read_csv(V22/'CA_joint_state_decomposition_per_record.csv')
    plotted=read_csv(V23/'figures/fig8_joint_partition_CA.all83_source_data.csv')
    points=read_csv(V23/'figures/fig8_joint_partition_CA.points.csv')
    stats=read_csv(V23/'figures/fig8_joint_partition_CA.statistics.csv')
    source=d[d.shift_degrees.isin([0,60])].set_index(['sample_id','shift_degrees'])
    cax=ca.set_index(['sample_id','shift_degrees'])
    px=plotted.set_index(['sample_id','shift_degrees'])
    check('NewFig8 bothpanels all83records/166 pairedsource rows',len(plotted)==166 and plotted.sample_id.nunique()==83 and set(px.index)==set(source.index)==set(cax.index))
    check('NewFig8 sourcevalues actual frozenjointH all166',np.allclose(px.h_joint.sort_index(),source.h_joint.sort_index(),rtol=0,atol=1e-12))
    check('NewFig8 CA ratios actual frozenperrecord all166',np.allclose(px.ca_between_joint_percent.sort_index(),100*cax.between_fraction_record_mean.sort_index(),rtol=0,atol=1e-12))
    for panel,column in [('a','h_joint'),('b','ca_between_joint_percent')]:
        actual=points[points.panel.eq(panel)&points.scope.eq('record')].set_index(['sample_id','shift_degrees'])
        check('Fig8 actualexport all166points '+panel,len(actual)==166 and set(actual.index)==set(px.index) and np.allclose(actual.y.sort_index(),px[column].sort_index(),atol=1e-12,rtol=0))
        for r in stats[stats.panel.eq(panel)].itertuples(index=False):
            q=plotted[plotted.shift_degrees.eq(r.shift_degrees)]
            if r.dataset!='all':q=q[q.dataset.eq(r.dataset)]
            check('Fig8 independent macro mean '+panel+'/'+r.dataset+'/'+str(r.shift_degrees),abs(q[column].mean()-r.mean)<1e-12 and len(q)==r.n_records and int(q.valid_residues.sum())==r.effective_residues)
    check('Fig8 main values round to1.484/.962 and6.2/7.8',
          [f'{source.xs(k,level=1).h_joint.mean():.3f}' for k in [0,60]]==['1.484','0.962'] and
          [f'{100*cax.xs(k,level=1).between_fraction_record_mean.mean():.1f}' for k in [0,60]]==['6.2','7.8'])
    check('Fig8 caption correcttwo panels jointH includesPro common33493 RMSFsquared',
          all(x in p[109] for x in ['（a）','（b）','Pro','33,493','RMSF²']) and '（c）' not in p[109] and '（d）' not in p[109])
    check('Freeenergy main sections and MD identity still present',
          all(x in p[56] for x in ['pmax','F']) and '−ln' in p[57] and 'pc > 0' in p[57] and
          any(x.startswith('3.6.1') for x in p) and any(x.startswith('3.6.2') for x in p) and any(x.startswith('3.6.3') for x in p))
    check('Windowmean/pooled occupancy distinction retained',
          '窗口' in p[53] and '合并' in p[53] and '不是' not in p[56])
    newmath=doc._element.xpath('.//m:oMath');oldmath=old._element.xpath('.//m:oMath')
    check('All native math semantic content preserved',len(newmath)==len(oldmath) and [semantic(x) for x in newmath]==[semantic(x) for x in oldmath])
    check('Fig8 isactual embedded PNG once, exactlytenactive images',len(doc._element.xpath('.//a:blip'))==10)
    with zipfile.ZipFile(DOC) as z,zipfile.ZipFile(OLD) as oz:
        embedded=[hashlib.sha256(z.read(x)).hexdigest() for x in z.namelist() if x.startswith('word/media/')]
        check('NewFigure8 actualPNG embedded',fa['outputs']['png']['sha256'] in embedded)
        comments=[x for x in oz.namelist() if 'comments' in x.lower() and x.endswith('.xml')]
        check('Original comment files preserved',all(x in z.namelist() and z.read(x)==oz.read(x) for x in comments))
    figcaptions=[re.match(r'^图(S?\d+)：',x).group(1) for x in p if re.match(r'^图(S?\d+)：',x)]
    tablecaptions=[re.match(r'^表(S?\d+)：',x).group(1) for x in p if re.match(r'^表(S?\d+)：',x)]
    check('Figure captions exactlymain1through9 plusS1',figcaptions==list(map(str,range(1,10)))+['S1'])
    check('Table captions exactlymain1and2 plusS1through9',tablecaptions==['1','2']+['S'+str(x) for x in range(1,10)])
    figrefs=set(re.findall(r'图(S?\d+)', '\n'.join(p)))
    tabrefs=set(re.findall(r'表(S?\d+)', '\n'.join(p)))
    check('Every remaining figure/table crossreference resolves',figrefs<=set(figcaptions) and tabrefs<=set(tablecaptions),{'figure_refs':sorted(figrefs),'table_refs':sorted(tabrefs)})
    check('Supplement headings exactlyAthroughG', [x[0] for x in p[160:] if re.match(r'^[A-Z] ',x)]==list('ABCDEFG'))
    stale=[{'paragraph':i,'text':x} for i,x in enumerate(p) if any(k in x for k in ['632','2,431','2431','五模型','事后校准','分组校准','图S2','图S3','图S4','图S5','图S6','图S7','图S8','图S9','图S10','图S11'])]
    check('No withdrawn ordinaryexternal/calibration/oldplot references or old2431protocol', not stale, stale)
    seed_precise='第0折' in p[168] and '第0折' in p[169]
    if not seed_precise:
        issues.append({'severity':'minor_protocol_clarification','paragraphs':[168,169],
                       'evidence':'Original corrected training_freeze.json seed_stability=fixed fold0 seeds7,11,23, separate from fold variation; expected_neural_fits=14 (10 originalOOF+4 addedseedfit).',
                       'problem':'固定五折体系划分下的表述没有明说随机初始化评价只运行固定第0折。',
                       'suggestion':'P168改为固定原体系划分的第0折训练与留出集合；表S2图注明确固定第0折。',
                       'changes_values':False})
    check('Seedstability scope explicitlyfixedfold0',seed_precise)
    if '均值与样本标准差' in p[49]:
        issues.append({'severity':'minor_reporting_clarification','paragraphs':[49],
                       'problem':'P49承诺随机初始化对照报告均值与样本标准差，表S2实际只报告各seed单模型分数。',
                       'suggestion':'改为报告每个随机初始化模型的留出成绩，或按真实结果补汇总。',
                       'changes_values':False})
    if '五折体系留出' not in p[120] and '相同训练集合的对照中降低了两项平均误差' in p[120]:
        issues.append({'severity':'minor_scope_clarification','paragraphs':[120],
                       'problem':'结论两项误差均降低缺少体系留出范围，外部RMSF MAE完整0.474高于链内与静态0.471。',
                       'suggestion':'在该句中明确五折体系留出评价范围；无需新增防御性段落。',
                       'changes_values':False})
    for path,digest in sources.items():
        check('Source unchanged after independent readonly review '+Path(path).name,sha(path)==digest)
    check('Candidate remained same during review',sha(DOC)==doc_hash and sha(TXT)==text_hash)
    substantive_fail=[x for x in checks if not x['passed'] and x['name']!='Seedstability scope explicitlyfixedfold0']
    report={'status':'passed' if not substantive_fail and not issues else 'passed_numeric_with_protocol_clarification' if not substantive_fail else 'failed',
            'timestamp_utc':datetime.now(timezone.utc).isoformat(),'document':str(DOC),'document_sha256':doc_hash,
            'text_sha256':text_hash,'frozen_v22_document_sha256':old_hash,
            'scope':'11actual native tables; exact v22cell equality; independent completedCSV rounding; actualFig8 all83data/points; nativeequations; remainingreferences and truthful scientific claims. No DOCX/source/oldcertificate mutations.',
            'table_audit':table_audit,'total_native_cells':sum(t['cells'] for t in table_audit),
            'checks':checks,'check_counts':{'total':len(checks),'passed':sum(x['passed'] for x in checks),'failed':sum(not x['passed'] for x in checks)},
            'issues':issues,'substantive_failures':substantive_fail,'source_sha256':sources,
            'visual_scope':'Actual Figure8PNG inspected with view_image: both panels, all paired records, legends, Chinese labels and units are legible without clipping. Whole manuscript page-layout QA belongs to root reviewer.',
            'figure8_visual_QA':{'status':'passed','image_sha256':'3d60badb79fc174351471436487d75062ebf12779329c33e2dbde63e7036208a',
                                 'actual_view':'Full exported2161x1192PNG via view_image original; model display resized2048x1130.',
                                 'details':'Two aligned panels; all record points/lines and source/overall means present; Chinese labels and nats/percent units complete; 33493positions and501frames defined; no clipping.'},
            'frozen_source_files_unchanged':all(sha(path)==digest for path,digest in sources.items())}
    (QA/'v23_results_independent_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    md='对实际DOCX SHA '+doc_hash+' 的独立检查：'+report['status']+'。11张原生表共'+str(report['total_native_cells'])+'个单元格与冻结v22原表逐项完全一致，表号重排与保留值正确；新表S3/S5及主表2另从完成源CSV重新取值并舍入，表S7重算83条记录的角区间统计、表S8重算KIX四状态、表S9重算组合计数与记录等权指标。图8的两面板分别使用83条记录、166个逐记录点，实际点文件与冻结联合Hχ和Cα状态间方差比例一致，均值1.484/0.962 nats与6.2%/7.8%正确；原生公式完整保留，剩余图表引用均可解析。正文保留不同步组合15.8%/16.5%的低召回，外部RMSF较链内与静态对照略差，以及KIX全部六对差值方向相反的真实结果。\n\n'
    if issues: md += '需要同步的文字：'+'；'.join('P'+','.join(map(str,x['paragraphs']))+' '+x['problem']+' 建议：'+x['suggestion'] for x in issues)+' 数值均无需改变。\n\n'
    if substantive_fail:md += '失败检查：'+json.dumps(substantive_fail,ensure_ascii=False)+'\n'
    md+='本轮只生成v23目录独立审阅文件，未改Word、旧CSV、图或证书；完整页面视觉排版由根代理核查。\n'
    (QA/'v23_results_independent_review.md').write_text(md,encoding='utf-8')
    print(json.dumps({'status':report['status'],'document_sha256':doc_hash,'counts':report['check_counts'],
                      'substantive_failures':substantive_fail,'issues':issues},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
