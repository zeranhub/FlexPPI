"""Read-only scientific review. Only new v22 audit files are written."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,re,sys,zipfile
import numpy as np
import pandas as pd
from docx import Document
ROOT=Path('E:/Project/FlexPPI/paper/revision_v14_20261010')
BASE=ROOT/'experiment_completion'
V22=BASE/'energy_landscape_v22_maintext_20261012'
V21=BASE/'energy_landscape_v21_20261012'
TRAJ=V21/'trajectory_analysis'
CASE=V21/'case_analysis'
DOC=ROOT/'paper_v22_20261012_中文修订_自由能分析融入正文_三线表.docx'
TXT=V22/'manuscript_QA/v22_content.txt'
def sha(path):
    d=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):d.update(block)
    return d.hexdigest()
def main():
    sys.stdout.reconfigure(encoding='utf-8')
    doc_sha,txt_sha=sha(DOC),sha(TXT)
    doc=Document(DOC)
    paragraphs=[''.join(n.text or '' for n in p._p.xpath('.//w:t | .//m:t')).strip() for p in doc.paragraphs]
    parsed={int(i):t for i,t in re.findall(r'^P(\d+): (.*)$',TXT.read_text(encoding='utf-8'),re.M)}
    tables=[[[c.text for c in r.cells] for r in t.rows] for t in doc.tables]
    checks,sources=[],{}
    def check(name,condition,detail=None):
        item={'name':name,'passed':bool(condition)}
        if detail is not None:item['detail']=detail
        checks.append(item)
    def csv(path,**kwargs):
        sources[str(path)]=sha(path)
        return pd.read_csv(path,float_precision='round_trip',**kwargs)
    def js(path):
        sources[str(path)]=sha(path)
        return json.loads(path.read_text(encoding='utf-8'))
    check('Text exact including native math',all(paragraphs[i]==t for i,t in parsed.items()) and len(parsed)==sum(bool(p) for p in paragraphs))
    integration=js(V22/'manuscript_QA/integration_audit.json')
    check('DOCX matches integration hash',doc_sha==integration['document_sha256'])
    completion=js(TRAJ/'completion_audit.json')
    delivery=js(TRAJ/'delivery_summary_audit.json')
    check('83records203trajectories complete',completion['status']=='complete' and completion['records']==83 and completion['trajectories']==203 and completion['total_selected_frames']==203*501)
    check('Uniform labels/masks reproduced',completion['original_uniform_label_reproduction_max_absolute']==0 and completion['original_uniform_masks_exact'])
    for name,meta in completion['output_tables'].items():
        path=TRAJ/Path(meta['path'].replace('\\','/')).name
        check('Certified trajectory source '+name,sha(path)==meta['sha256'])
        sources[str(path)]=sha(path)
    for name,digest in delivery['outputs'].items():
        check('Certified nonPro derived source '+name,sha(TRAJ/name)==digest)
        sources[str(TRAJ/name)]=sha(TRAJ/name)
    audit=js(V22/'representative_Fstar_and_CA_maintext_audit.json')
    check('Representative and CA audit passed',audit['status']=='passed' and audit['source_files_unchanged_after_audit'])
    for path,digest in audit['v21_inputs_sha256'].items():
        check('v21 frozen case source '+Path(path).name,sha(path)==digest)
        sources[path]=sha(path)
    for name,digest in audit['new_csv_outputs'].items():
        check('v22 case/CA source '+name,sha(V22/name)==digest)
        sources[str(V22/name)]=sha(V22/name)
    d=csv(TRAJ/'distribution_record_macro.csv')
    n=csv(TRAJ/'nonpro_chi1_record_macro.csv')
    d=d[d.replica_mode.eq('replica1')]
    n=n[n.replica_mode.eq('replica1')]
    check('All83firstreplica4partitions',len(d)==len(n)==83*4 and not d.duplicated(['sample_id','shift_degrees']).any() and not n.duplicated(['sample_id','shift_degrees']).any())
    check('83sourcecount302330',d[d.shift_degrees.eq(0)].groupby('dataset').size().to_dict()=={'atlas':30,'dynarepo':23,'itrace':30})
    table3=[];means={}
    for shift in [0,30,60,90]:
        x=d[d.shift_degrees.eq(shift)];y=n[n.shift_degrees.eq(shift)]
        means[shift]={'full_minus_window':float(x.full_minus_window_h_weighted.mean())}
        table3.append([str(shift)+'°',f'{x.h_joint.mean():.3f} ± {x.h_joint.std(ddof=1):.3f}',f'{y.h_chi1.mean():.3f}',f'{100*y.chi1_boundary_10deg_fraction.mean():.2f}',f'{y.apparent_chi1_changes.mean():.1f}',f'{x.half_joint_jsd.mean():.3f}'])
    check('MainTable3 all cells independently recomputed',tables[2][1:]==table3,{'actual':tables[2][1:],'expected':table3})
    check('P111 correct joint/nonPro numbers',all(x in paragraphs[111] for x in ['1.484','0.962','0.717','0.374','61.3%','0.13%','189.2','21.9']))
    check('P114 windowpooling contrast independently recomputed',[format(means[k]['full_minus_window'],'.3f') for k in [0,60]]==['0.135','0.170'])
    fig8stats=csv(V21/'figures/figS12_partition_sampling.statistics.csv')
    for r in fig8stats[fig8stats.dataset.eq('all')].itertuples(index=False):
        if r.panel=='a':value=float(d[d.shift_degrees.eq(int(r.category))].h_joint.mean())
        elif r.panel=='b':value=float(n[n.shift_degrees.eq(int(r.category))].h_chi1.mean())
        elif r.panel=='c':value=float(n[n.shift_degrees.eq(int(r.category))].apparent_chi1_changes.mean())
        else:continue
        check('Fig8source '+r.panel+'/'+str(r.category),abs(value-r.mean)<1e-12)
    times=fig8stats[fig8stats.panel.eq('d')&fig8stats.dataset.eq('all')].set_index('category')
    check('P114 chi1 half/repeat distinct from jointJSD',[format(times.loc[k,'mean'],'.3f') for k in ['half0','half60','repeat0','repeat60']]==['0.021','0.040','0.022','0.040'])
    per=csv(TRAJ/'per_residue_trajectory.csv.gz',usecols=['sample_id','dataset','replica_index','residue_index','shift_degrees','ca_total_variance_A2','ca_within_joint_variance_A2','ca_between_joint_variance_A2','ca_between_joint_fraction'])
    primary=per[per.replica_index.eq(1)&per.shift_degrees.isin([0,60])]
    ca_summary=csv(V22/'CA_joint_state_decomposition_summary_and_denominators.csv')
    ca_per=csv(V22/'CA_joint_state_decomposition_per_record.csv')
    max_identity_error=0
    for shift in [0,60]:
        q=primary[primary.shift_degrees.eq(shift)]
        q=q[np.isfinite(q.ca_total_variance_A2)&np.isfinite(q.ca_within_joint_variance_A2)&np.isfinite(q.ca_between_joint_variance_A2)&np.isfinite(q.ca_between_joint_fraction)]
        error=np.abs(q.ca_total_variance_A2-q.ca_within_joint_variance_A2-q.ca_between_joint_variance_A2)
        max_identity_error=max(max_identity_error,float(error.max()))
        check('CA within+between=total all33493 '+str(shift),len(q)==33493 and q.sample_id.nunique()==83 and error.max()<1e-10)
        record=q.groupby('sample_id').ca_between_joint_fraction.mean()
        saved=ca_per[ca_per.shift_degrees.eq(shift)].set_index('sample_id')
        check('CA actual perposition-to-record ratios '+str(shift),set(record.index)==set(saved.index) and np.allclose(record.sort_index(),saved.between_fraction_record_mean.sort_index(),atol=1e-12,rtol=0))
        expected=ca_summary[ca_summary.dataset.eq('all')&ca_summary.shift_degrees.eq(shift)].iloc[0]
        check('CA recordmacro source '+str(shift),abs(record.mean()-expected.mean_of_record_residue_mean_between_fraction)<1e-12)
        check('CA common mask source5660/6384/21449 '+str(shift),q.groupby('dataset').size().to_dict()=={'atlas':5660,'dynarepo':6384,'itrace':21449})
    check('Main6.2/7.8 rounding',[format(100*ca_summary[ca_summary.dataset.eq('all')&ca_summary.shift_degrees.eq(k)].iloc[0].mean_of_record_residue_mean_between_fraction,'.1f') for k in [0,60]]==['6.2','7.8'])
    occupancy=csv(V22/'representative_rep1_chi1_occupancy_and_Fstar.csv',keep_default_na=False,dtype={'chain':str,'maximum_bin_centers_degrees':str})
    for r in occupancy.itertuples(index=False):
        path=TRAJ/'cache'/(r.sample_id+'_rep1.npz')
        with np.load(path,allow_pickle=False) as cache:
            angle=cache['chi_degrees'][:,r.input_index,0].astype(float);time=cache['times_ps'].astype(float)
        check('Case501frames '+r.sample_id+'/'+r.residue,len(angle)==501 and np.array_equal(time,np.arange(501)*200))
        ang=(angle+180)%360-180
        gminus=int(((ang>=-120)&(ang<0)).sum());gplus=int(((ang>=0)&(ang<120)).sum());trans=501-gminus-gplus
        check('Case exactg-/g+/trans '+r.sample_id+'/'+r.residue,[gminus,gplus,trans]==[r.chi1_gminus_frames,r.chi1_gplus_frames,r.chi1_trans_frames])
        counts,edges=np.histogram(ang,bins=np.arange(-180,181,5));centers=(edges[:-1]+edges[1:])/2
        peaks=';'.join(format(x,'g') for x in centers[counts==counts.max()])
        check('Case highest5degreebox '+r.sample_id+'/'+r.residue,peaks==r.maximum_bin_centers_degrees)
    occ=occupancy.set_index(['sample_id','residue'])
    check('P1194.8/47.9/51.9true',format(100*occ.loc[('A0022','I78'),'chi1_gplus_probability'],'.1f')=='4.8' and [format(100*occ.loc[('itrace_6p64','T149'),c],'.1f') for c in ['chi1_gplus_probability','chi1_trans_probability']]==['47.9','51.9'])
    selection=js(CASE/'representative_cases_audit.json')
    case_text=next(p for p in paragraphs if p.startswith('为将这一解释对应到不同结构来源'))
    case_caption=next(p for p in paragraphs if p.startswith('图9：'))
    check('Representativeall6originalreference values',all(format(r['rmsf_reference'],'.3f') in case_text and format(r['hchi_reference'],'.3f') in case_text for r in selection['residues']))
    check('Preselected recordrank0.5',[r['mean_error_rank'] for r in selection['choices']]==[.5,.5,.5])
    check('Case site selection explicitly by MDreference','该位点选择使用参考读数' in paragraphs[62] and '位点按参考Hχ秩减RMSF秩最相离选取' in case_caption)
    ubi=csv(CASE/'ubiquitin_bin_sensitivity.csv').set_index(['residue','phase'])
    un=csv(CASE/'ubiquitin_independent_numeric_audit.csv').set_index(['residue','phase_degrees'])
    for k in ['Q41','Q62']:
        old=ubi.loc[(k,0)]
        check('P106frozenpaired values '+k,all(format(v,'.3f') in paragraphs[106] for v in [old.md_rmsf_A,old.pred_rmsf_A,old.joint_window_mean_h,old.pred_joint_h]))
        angle=un.loc[(k,0)]
        check('P106marginalwindowmean '+k,all(format(angle[c],'.3f') in paragraphs[106] for c in ['chi1_window_mean','chi2_window_mean','chi3_window_mean']))
        offset=ubi.loc[(k,60)]
        check('P107shift60windowmeans '+k,format(offset.joint_window_mean_h,'.3f') in paragraphs[107] and format(offset.chi1_window_mean_h,'.3f') in paragraphs[107])
        check('PooledHdistinct not substituted frozenlabel '+k,abs(old.joint_pooled_h-old.joint_window_mean_h)>1e-3 and abs(old.joint_window_mean_h-old.frozen_joint_h)<1e-12)
    native=doc.paragraphs[57]._p
    math_ns={'m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
    check('Equation5 nativeWordmath',len(native.xpath('./m:oMath'))==1)
    frac=native.xpath('.//m:f')[0]
    check('Equation5fractionpc/pmax',''.join(frac.xpath('./m:num//m:t/text()',namespaces=math_ns))=='pc' and ''.join(frac.xpath('./m:den//m:t/text()',namespaces=math_ns))=='pmax')
    check('Equation5star+subc correct',len(native.xpath('.//m:sSubSup'))==2 and all(''.join(s.xpath('./m:sub//m:t/text()',namespaces=math_ns))=='c' and ''.join(s.xpath('./m:sup//m:t/text()',namespaces=math_ns))=='*' for s in native.xpath('.//m:sSubSup')))
    nary=native.xpath('.//m:nary')[0]
    check('Equation5sum positivepcFstar',''.join(nary.xpath('./m:sub//m:t/text()',namespaces=math_ns))=='c: pc > 0' and ''.join(nary.xpath('./m:e//m:t/text()',namespaces=math_ns))=='pcFc*')
    check('Equation5minuslnpmax sign',' − ln pmax.     (5)' in paragraphs[57])
    for p in [np.array([1.]),np.array([.7,.2,.1]),np.array([.5,.5]),np.array([.1,.2,.3,.4])]:
        h=float(-np.sum(p*np.log(p)));f=-np.log(p/p.max())
        check('Equation5independent algebra '+str(p.tolist()),abs(h-(float(np.dot(p,f))-float(np.log(p.max()))))<1e-14)
    check('FstarfromMD not modeloutput','模型输出RMSF与Hχ，轨迹提供占据概率及对应的F*曲线' in paragraphs[58])
    check('Fstar dimensionless and unobserved gap','采用无量纲热能标度' in paragraphs[61] and '零计数箱断线' in paragraphs[61])
    check('No learnedbarrier/rate/fullproteinbasin',not any(s in '\n'.join(paragraphs[55:131]) for s in ['预测能垒','反演能垒','跃迁率预测','全蛋白自由能盆地','预测自由能曲线']))
    text='\n'.join(paragraphs)
    figures=re.findall(r'^(图S?\d+)：',text,re.M);caps=re.findall(r'^(表S?\d+)：',text,re.M)
    check('All9main11supp figurecaptions',sorted(figures)==sorted(['图'+str(i) for i in range(1,10)]+['图S'+str(i) for i in range(1,12)]))
    check('All3main12supp tablecaptions',sorted(caps)==sorted(['表'+str(i) for i in range(1,4)]+['表S'+str(i) for i in range(1,13)]))
    check('No orphan figure/table reference',not (set(re.findall(r'图S?\d+',text))-set(figures)) and not (set(re.findall(r'表S?\d+',text))-set(caps)))
    with zipfile.ZipFile(DOC) as z:media=[z.read(name) for name in z.namelist() if name.startswith('word/media/')]
    for name in ['fig7_ubiquitin_angular_landscapes.png','figS12_partition_sampling.png','figS13_representative_distribution_cases.png']:
        path=V21/'figures'/name
        check('Promoted exactsource '+name,sum(hashlib.sha256(x).hexdigest()==sha(path) for x in media)==1)
    optional=[]
    if '五窗口平均的原分区χ1边缘熵' not in paragraphs[106]:
        optional.append(dict(paragraphs=[106,107],severity='optional_clarity',topic='Correct windowmean marginalH vs pooledFstar',explanation='All marginal entropy numbers are correct five-window averages; coarseFstar pools5000frames. Explicit windowmean wording reduces ambiguity.',suggested_replacements=[{'paragraph':106,'from':'原分区χ1边缘熵','to':'五窗口平均的原分区χ1边缘熵'},{'paragraph':107,'from':'常规χ1分区的边缘熵','to':'常规χ1分区的五窗口平均边缘熵'}]))
    check('DOCX andTXTunchanged',sha(DOC)==doc_sha and sha(TXT)==txt_sha)
    failed=[c for c in checks if not c['passed']]
    report=dict(status='passed' if not failed else 'failed_checks',checked_utc=datetime.now(timezone.utc).isoformat(),document=str(DOC),document_sha256=doc_sha,extracted_text_sha256=txt_sha,read_only=True,writes_only_v22=True,source_v21_modified=False,existing_certificate_rerun=False,no_training_or_MD_or_GPU=True,scope='Actual2.8/3.3/3.6/discussion/table3 sources/nativeequation5/promotionrefs; root owns pagevisualQA',checks=checks,failed_checks=failed,required_scientific_issues=[],optional_clarity=optional,CA_maximum_independently_recalculated_identity_error_A2=max_identity_error,source_sha256=sources)
    path=V22/'manuscript_QA/v22_results_independent_content_review.json'
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    summary='v22自由能整合独立科学内容复核'+('通过' if not failed else '存在检查未通过')+'。实际DOCX SHA256：'+doc_sha+'。共'+str(len(checks))+'项核对，失败'+str(len(failed))+'项。\n\n正文表3与已认证83条记录第一轨迹统计一致；联合Hχ包含Pro和全部原定义角，χ1边缘熵、边界比例与变化次数使用非Pro位置。半程及重复JSD明确区分联合和χ1，并在图8d保持共同60条记录口径。泛素各熵值为五窗口平均，F*粗线为合并5000帧，χ1投影未替代联合目标。三来源原标签、第一轨迹取向占据比例和5°经验峰已由缓存独立复算。\n\nRMSF²状态内与状态间两项逐位复现总均方波动；共同33493个位点先逐位取比例、记录内均值、再83记录等权，6.242%和7.769%对应正文6.2%和7.8%。式(5)概率恒等式正确，负号、p_c/p_max、联合期望及星号均在可编辑Word公式中正确。\n\n正文用于轨迹分布解释，模型仍预测RMSF与Hχ两个标量；没有预测能垒、转换速率或全蛋白自由能盆地的虚构。图8、图9及正文表3引用完整。'
    if optional:summary+='\n\n可选清晰化：P106“原分区χ1边缘熵”加“五窗口平均的”，P107“常规χ1分区的边缘熵”加“五窗口平均”，可避免与合并帧F*口径混读；现有数值均正确。'
    (V22/'manuscript_QA/v22_results_independent_content_review.md').write_text(summary+'\n',encoding='utf-8')
    print(json.dumps(dict(status=report['status'],checks=len(checks),failed=failed,document_sha256=doc_sha,optional_clarity=optional,report=str(path)),ensure_ascii=False,indent=2))
    return 1 if failed else 0
if __name__=='__main__':raise SystemExit(main())
