from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
checks=[]
def check_file(path,digest):
    assert sha(path)==digest,str(path)
    checks.append({'file':str(path),'sha256':digest})
ma=readj(ROOT/'model/verification_audit.json');assert ma['status']=='passed' and ma['checks_passed']==73
for name,digest in ma['outputs_sha256'].items():check_file(ROOT/'model'/name,digest)
ta=readj(ROOT/'trajectory/diagnostic_audit.json');assert ta['status']=='complete'
for name,info in ta['outputs'].items():check_file(ROOT/'trajectory'/name,info['sha256'])
sa=readj(ROOT/'structure/run_audit.json');assert sa['status']=='complete'
for name,info in sa['outputs'].items():check_file(ROOT/'structure'/name,info['sha256'])
sv=readj(ROOT/'structure/short_verification_and_summary_audit.json');assert sv['status']=='passed' and sv['check_counts']=={'total':110,'passed':110}
check_file(ROOT/'structure/run_audit.json',sv['original_run_audit_sha256'])
for name,info in sv['new_outputs'].items():check_file(ROOT/'structure'/name,info['sha256'])
model=pd.read_csv(ROOT/'model/all6_Gln_frozen_model_comparison.csv')
traj=pd.read_csv(ROOT/'trajectory/six_GLN_prediction_window_summary.csv')
entropy=pd.read_csv(ROOT/'trajectory/original120_entropy_by_window_and_scope.csv')
contacts=pd.read_csv(ROOT/'structure/contact_summary_all6_by_grid_window.csv')
sasa=pd.read_csv(ROOT/'structure/static_SASA_all6.csv')
fine=pd.read_csv(ROOT/'trajectory/fine5degree_501_vs5000_summary.csv')
state_probs=pd.read_csv(ROOT/'trajectory/original120_per_chi_state_probabilities.csv')
assert list(model.residue_id)==[2,31,40,41,49,62] and set(traj.residue)=={'Q2','Q31','Q40','Q41','Q49','Q62'}
q=model[model.residue_id==40].iloc[0];qt=traj[traj.residue=='Q40'].iloc[0]
np.testing.assert_allclose(model.full_Hchi_nats-model.MD_Hchi_nats,model.full_signed_error_nats,rtol=0,atol=1e-14)
windows=entropy[(entropy.residue=='Q40')&entropy.scope.str.startswith('window')].sort_values('window').copy()
assert len(windows)==5
assert abs(windows.joint_Hchi_nats.mean()-q.MD_Hchi_nats)<1e-6
assert abs(windows.joint_Hchi_nats.std(ddof=1)-q.label_window_SD_nats)<1e-7
assert sum(windows.predicted_minus_MD_nats>0)==4 and sum(windows.predicted_minus_MD_nats<0)==1
assert qt.MD_window_sample_SD_nats==traj.MD_window_sample_SD_nats.max()
assert abs(windows.marginal_entropy_sum_nats.mean()-qt.mean_window_marginal_entropy_sum_nats)<1e-14
assert abs(windows.total_correlation_nats.mean()-qt.mean_window_total_correlation_nats)<1e-14
q40contacts=contacts[(contacts.site=='Q40')&(contacts.grid=='original5000')].iloc[0]
assert q40contacts.nonadjacent_4A_interchain_any_mean==0
q40SASA=sasa[sasa.site=='Q40'].sidechain_SASA_A2.iloc[0]
assert 70<q40SASA<71
summary=model[['residue','MD_Hchi_nats','full_Hchi_nats','full_signed_error_nats','intra_Hchi_nats','static_Hchi_nats','AA_balanced_Hchi_nats','five_window_std','five_window_min','five_window_max']].merge(sasa[['site','sidechain_SASA_A2']],left_on='residue',right_on='site',validate='one_to_one').drop(columns='site').merge(fine[['residue','MD_width_5000_four_origins_mean_degrees','MD_fine_entropy_5000_four_origins_mean_nats']],on='residue',validate='one_to_one')
summary.to_csv(ROOT/'six_Gln_integrated_comparison.csv',index=False,float_format='%.17g')
windows=windows.merge(contacts[contacts.site=='Q40'][['grid','polar3p5_samechain_any_mean','polar3p5_nonadjacent_samechain_any_mean','nonadjacent_4A_samechain_residue_count_mean']],left_on='scope',right_on='grid',validate='one_to_one').drop(columns='grid')
windows.to_csv(ROOT/'Q40_five_windows_integrated.csv',index=False,float_format='%.17g')
four=windows.query('window<=4');last=windows.query('window==5').iloc[0]
increments={k:float(last[k]-four[k].mean()) for k in ['joint_Hchi_nats','Hchi1_marginal_nats','Hchi2_marginal_nats','Hchi3_marginal_nats','marginal_entropy_sum_nats','total_correlation_nats']}
assert abs(increments['joint_Hchi_nats']-(increments['marginal_entropy_sum_nats']-increments['total_correlation_nats']))<1e-14
detail={'status':'complete','finished_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'Q40_frozen_reference_Hchi':float(q.MD_Hchi_nats),'Q40_full_prediction_Hchi':float(q.full_Hchi_nats),'Q40_error_Hchi':float(q.full_signed_error_nats),'five_windows_joint_Hchi':windows.joint_Hchi_nats.tolist(),'five_window_sample_SD':float(q.five_window_std),'last_vs_first4_entropy_increments':increments,'actual_marginal_sum_mean':float(qt.mean_window_marginal_entropy_sum_nats),'actual_total_correlation_mean':float(qt.mean_window_total_correlation_nats),'observed_marginal_sum_below_model_prediction':float(q.full_Hchi_nats-qt.mean_window_marginal_entropy_sum_nats),'observations':['All four frozen predictors overestimate Q40; Full error is largest.','Q40 has the largest five-window reference variation among six GLN sites.','Q40 is exposed in the initial structure; no 4A direct UBA sidechain contact in 5000 saved frames.','Q40 has intermittent Pro37/terminal Arg72 intrachain contacts; no unique persistent polar lock or monotone last-window contact release was identified.','Coarse-state angular dependence is not exceptional compared with Q41.','Coarse-grid sensitivity is not largest for Q40.'],'supported_hypothesis':'The learned mapping from backbone environment to sidechain state occupation overestimates the diversity of this exposed but relatively concentrated site; temporal changes in occupation make a single-structure estimate of the fixed five-window average harder. Specific molecular causation is unresolved.','causal_limit':'No same-weight partner deletion, molecular intervention, independent MD repeat, or new training was performed. Contacts are distance-based observations, not verified hydrogen bonds.','paper_modified':False,'new_MD_production':False,'new_training':False,'local_output_hashes_checked':len(checks)}
write(ROOT/'integrated_diagnosis.json',detail)
report=f'''# Q40侧链联合熵预测偏差诊断

本次复核使用A003R（2QHO）的原第0折冻结输出、原500 ns单条MD轨迹的前5,000帧及原五个连续1,000帧窗口，同时保留全部六个谷氨酰胺作为对照。没有修改论文、模型或参考标签，没有新跑MD或训练。

Q40的偏差方向明确：模型预测{q.full_Hchi_nats:.6f} nats，原五窗口参考{q.MD_Hchi_nats:.6f} nats，高估{q.full_signed_error_nats:.6f} nats。现有证据最支持的解释是：Q40初始结构较暴露，但其实际角度状态占据较集中且随时间段变化，模型从静态主链环境到侧链状态多样性的映射在这里给了过高读数。该解释包含已测事实与机制推测，尚未定位到某一条相互作用的因果作用。

## 1. 时间段之间的占据差异最突出

| 原100 ns窗口 | MD联合Hχ（nats） | 模型减该窗口参考（nats） |
|---|---:|---:|
| 0–99.9 ns | {windows.iloc[0].joint_Hchi_nats:.6f} | {windows.iloc[0].predicted_minus_MD_nats:+.6f} |
| 100–199.9 ns | {windows.iloc[1].joint_Hchi_nats:.6f} | {windows.iloc[1].predicted_minus_MD_nats:+.6f} |
| 200–299.9 ns | {windows.iloc[2].joint_Hchi_nats:.6f} | {windows.iloc[2].predicted_minus_MD_nats:+.6f} |
| 300–399.9 ns | {windows.iloc[3].joint_Hchi_nats:.6f} | {windows.iloc[3].predicted_minus_MD_nats:+.6f} |
| 400–499.9 ns | {windows.iloc[4].joint_Hchi_nats:.6f} | {windows.iloc[4].predicted_minus_MD_nats:+.6f} |

Q40的窗口样本标准差为{q.five_window_std:.6f} nats，是六个GLN中最大；其余五个位点最大的为0.118428 nats。模型值接近最后一窗，而高于前四窗；固定评价目标仍是五窗平均，不能挑选末窗替换参考，也不能据此认为网络预测了末窗状态。这里能确认的是同一条轨迹的取向占据随时间明显变化。五段连续时间窗口不是五次独立实验，标准差不作为均值的置信区间。

最后一窗相对前四窗均值，联合Hχ增加{increments['joint_Hchi_nats']:.6f} nats；三个角各自的粗状态熵增加{increments['Hchi1_marginal_nats']:.6f}、{increments['Hchi2_marginal_nats']:.6f}、{increments['Hchi3_marginal_nats']:.6f} nats，χ2增幅最大。例如果将χ2分成原三个120°区间，第4窗三个区间的占比为19.1%、80.1%、0.8%，第5窗变为28.5%、53.6%、17.9%。也就是它从主要停留在一个区间，变得更频繁地访问另外两个区间。

## 2. 单角占据集中与局部几何映射误差

在原120°状态定义下，Q40三个角的五窗平均边缘熵分别为0.811216、0.761895、0.819841 nats，和为2.392953；实际联合熵为2.186479，两者相差0.206473。这个差值表示三个角的组合具有依赖。Q41对应差值为0.199947，接近Q40，因此没有发现Q40具有特别强的角间依赖。即使保留Q40实际单角分布、暂时把三个角视为独立，其联合熵上限也只有2.392953，仍低于模型的2.554096；仅用“网络忽略角间耦合”解释误差并不充分。这里的边缘熵分解用于诊断真实占据，模型本身只输出联合Hχ，不能反推出网络给每个角预测了什么。

在初始PDB中，Q40侧链溶剂可及面积为70.636 Å²，与Q2的70.659及Q49的69.100接近，Q41为0。Q40同链非相邻Cα的8 Å邻居为6个，也比Q41的10个少。从静态主链环境看，它比Q41更像一个较暴露的位置；而原500 ns中，Q40的逐角低自由能宽度平均为32.917°/角，接近Q41的30.000，低于Q2的69.167和Q49的98.750。宽度与120°联合熵是不同读数，这些值支持“暴露程度不足以确定实际侧链状态占据”，不能把宽度换算成模型Hχ误差。

| Q40的读数 | Hχ（nats） | 相对原参考误差（nats） |
|---|---:|---:|
| 原MD参考 | {q.MD_Hchi_nats:.6f} | 0 |
| 完整模型 | {q.full_Hchi_nats:.6f} | {q.full_signed_error_nats:+.6f} |
| 链内模型 | {q.intra_Hchi_nats:.6f} | {q.intra_signed_error_nats:+.6f} |
| 静态特征对照 | {q.static_Hchi_nats:.6f} | {q.static_signed_error_nats:+.6f} |
| 平衡类型均值对照 | {q.AA_balanced_Hchi_nats:.6f} | {q.AA_balanced_signed_error_nats:+.6f} |

四种输出均高估Q40，完整模型又比链内模型高0.169652 nats。完整模型将Q40排在Q31之上约0.092 nats，原参考中Q40却比Q31低约0.252 nats，显示该处局部几何到侧链多样性的映射发生了偏差。完整模型也高于训练GLN普通均值2.370698及平衡均值2.397724，所以不是简单返回氨基酸平均值。完整与链内是分别训练的网络，其差值不能单独归因为UBA伙伴链。

## 3. 接触与分箱核查给出的边界

5000个保存帧中，Q40与UBA链B没有4 Å内的直接侧链重原子接触，最近距离仍为5.743 Å。Q40与同链Pro37及Arg72的4 Å接触占比为87.58%和74.64%；NE2到Arg72末端羧基OT1/OT2的3.5 Å候选极性近接分别为23.82%和22.32%。这些是间歇接触，并非已确认的持续氢键锁定；尤其OT1/OT2属于输入链末端羧基，不能误称Arg侧链盐桥。

第4到第5窗Hχ升高时，排除相邻残基后的极性近接比例由27.8%增到32.4%，NE2–OT1/OT2近接也增加。现有结果没有支持“某个远程极性锁在末窗解除”的简单机制。初始结构也没有Q40的OE1/NE2到其他蛋白N/O的3.5 Å近接。因此具体侧链约束如何造成预测误差仍是待验证假设。

诊断性地将原120°状态网格平移24种原点后，Q40平均Hχ的跨原点范围为0.879 nats，而预测较准的Q41、Q62分别为1.381及1.248 nats；Q40并非最敏感者，不能把异常误差特别归因于分箱。平移后的标签没有用于重新计分或改变参考。原五窗标签均值及方差已逐项复算匹配。

## 可用于解释结果的结论

Q40属于本例中“初始结构较暴露、MD占据偏集中且随时段变化”的位置。模型高估了这种位置的长期平均侧链状态多样性，完整几何模型在该位置的修正又进一步偏高。目前能定位到预测映射与时间占据统计层面的原因，不能把误差确定归因于UBA结合、深度埋藏、独有强耦合或某条持续氢键。

模型输出与训练均值口径见model/diagnosis_summary.json，模型73项核验见model/verification_audit.json；原五窗及逐角分解见trajectory，最大原窗口复算差4.44×10⁻¹⁶；结构接触110项核验见structure/short_verification_and_summary_audit.json；分箱敏感性见partition。角度3.5 Å近接未使用氢原子和角度，因此只称候选极性接触。论文v33保持原样，这份诊断不是新增正文机制论断。
'''
report_path=ROOT/'Q40预测偏差分析_20261015.md';report_path.write_text(report,encoding='utf8')
write(ROOT/'completion_audit_20261015.json',{'status':'complete','local_output_hashes_rechecked':checks,'all6_sites_retained':True,'five_windows_verified':True,'model_independent_checks':73,'structure_independent_checks':110,'model_and_structure_unresolved_causation_clearly_stated':True,'report_sha256':sha(report_path),'manuscript_modified':False})
sources=ROOT/'source';sources.mkdir(exist_ok=True)
case=ROOT.parent/'free_energy_curves_v27_20261014/case_curves'
for path in [case/'angles_5000frames_all6GLN.csv',case/'joint_entropy_original_five_windows.csv',case/'selection_inventory_with_MD_outcomes.csv',Path(r'D:\data\dynarepo_non_nucleic_full\A003R\structure.pdb')]:shutil.copy2(path,sources/path.name)
ms=readj(ROOT/'model/diagnosis_summary.json')
for path,digest in ms['sources_sha256'].items():
    p=Path(path)
    if p.name=='A003R.npz':
        check_file(p,digest)
        shutil.copy2(p,sources/(p.parent.name+'_A003R.npz'))
files=[]
for p in ROOT.rglob('*'):
    if not p.is_file() or '_staging' in p.parts or p.name.startswith('diagnose_before') or p.name in ['package_manifest_20261015.json']:continue
    files.append((p,p.relative_to(ROOT).as_posix()))
manifest={'report':str(report_path),'report_sha256':sha(report_path),'source_paths_recorded_in_per_module_audits':True,'raw_MD_not_duplicated':True,'files':{name:sha(p) for p,name in files}}
write(ROOT/'package_manifest_20261015.json',manifest);files.append((ROOT/'package_manifest_20261015.json','package_manifest_20261015.json'))
package=BASE/'Q40预测偏差分析与源数据_20261015.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path,name in files:z.write(path,name)
with zipfile.ZipFile(package) as z:assert z.testzip() is None and all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
print(json.dumps({'status':'complete','report':str(report_path),'package':str(package),'package_MB':round(package.stat().st_size/1e6,2),'outputs_verified':len(checks),'Q40_reference':float(q.MD_Hchi_nats),'Q40_prediction':float(q.full_Hchi_nats),'causation':'not_identified'},ensure_ascii=False))
