"""Export and independently verify available single-trajectory DynaRepo reference labels.
No simulation is run. Original files are read-only. Uses original source label functions.
"""
from pathlib import Path
import sys,types,json,csv,hashlib,platform
import numpy as np
PROJECT=Path(r'E:\Project\FlexPPI')
ROOT=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028\md_reference')
ROOT.mkdir(parents=True,exist_ok=True)
# Import only analysis modules without the unrelated torch-requiring package __init__.
for name,path in [('flexppi',PROJECT/'src/flexppi'),('flexppi.data',PROJECT/'src/flexppi/data')]:
 mod=types.ModuleType(name);mod.__path__=[str(path)];sys.modules[name]=mod
from flexppi.data.structures import parse_pdb
from flexppi.data.records import StateRecord
from flexppi.data.labels import compute_dynarepo_labels
from flexppi.constants import CHI_ATOMS
import MDAnalysis as mda

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def writecsv(p,rows):
 with p.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

allrows=[];summaries=[];audits=[]
for sid,pdb in [('A00EJ','6P4B'),('A00KW','1MLC'),('A00HE','1VFB')]:
 folder=Path(r'D:\data\dynarepo_non_nucleic_full')/sid
 structure=folder/'structure.pdb';trajectory=folder/'trajectory.xtc'
 lp=PROJECT/'artifacts/data/dynarepo_labels'/f'{sid}.npz'
 meta=json.loads((folder/'dynarepo_metadata.json').read_text(encoding='utf-8'))
 st=parse_pdb(str(structure))
 with np.load(lp,allow_pickle=False) as z:
  values=z['values'].copy();variances=z['variances'].copy();mask=z['mask'].copy();names=z['metric_names'].tolist()
 assert names==['rmsf','chi_entropy','sasa']
 assert values.shape==(st.num_residues,3)
 u=mda.Universe(str(structure),str(trajectory))
 first=u.trajectory[0];start=float(first.time)
 last=u.trajectory[-1];end=float(last.time)
 nframes=len(u.trajectory)
 record=StateRecord(sample_id=sid,dataset='dynarepo',structure_path=str(structure),label_path=None,trajectory_paths=[str(trajectory)],temperature=float(meta.get('temperature',310)),window_ns=100.0,target_chains=st.unique_chains(),source_group='antigen_single' if sid!='A00KW' else 'dynabench_complex')
 out=ROOT/f'{sid}_recomputed_1000frame_windows.npz'
 result=compute_dynarepo_labels(record,str(out),window_frames=1000,frame_stride=1,include_sasa=False,max_windows=None)
 with np.load(out,allow_pickle=False) as z:
  newvalues=z['values'];newvars=z['variances'];newmask=z['mask']
  assert np.array_equal(mask,newmask),(sid,'mask mismatch')
  maxval=float(np.nanmax(np.abs(values-newvalues)));maxvar=float(np.nanmax(np.abs(variances-newvars)))
  assert float(np.nanmax(np.abs(values[:,0]-newvalues[:,0])))<2e-5,(sid,'RMSF mismatch')
  discrepancies=[{'model_index':int(i),'chain':st.chain_ids[i],'author_residue_id':st.author_residue_ids[i],'aa':st.sequence[i],'metric':names[j],'cached_value':float(values[i,j]),'recomputed_value':float(newvalues[i,j]),'absolute_difference':float(abs(values[i,j]-newvalues[i,j]))} for i,j in zip(*np.where(np.abs(values-newvalues)>1e-5))]
 role=lambda ch:('antigen' if sid!='A00KW' or ch=='E' else ('antibody_light' if ch=='A' else 'antibody_heavy'))
 for i in range(st.num_residues):
  row={'sample_id':sid,'source_pdb':pdb,'model_index':i,'chain':st.chain_ids[i],'author_residue_id':st.author_residue_ids[i],'resname':st.residue_names[i],'aa':st.sequence[i],'role':role(st.chain_ids[i]),'geometry_valid':bool(st.frame_mask[i]),'chemically_defined_chi_count':len(CHI_ATOMS.get(st.residue_names[i],[])),'local_independent_trajectory_files':1,'label_windows':result['windows']}
  for j,name in enumerate(names[:2]):
   row[name]=float(values[i,j]) if mask[i,j] else ''
   row[name+'_variance_across_windows']=float(variances[i,j])
   row[name+'_valid']=bool(mask[i,j])
   row[name+'_recomputed']=float(newvalues[i,j]) if newmask[i,j] else ''
  allrows.append(row)
 for ch in st.unique_chains():
  ix=st.indices_for_chain(ch)
  for j,name in enumerate(names[:2]):
   valid=ix[mask[ix,j]]
   summaries.append({'sample_id':sid,'source_pdb':pdb,'chain':ch,'role':role(ch),'metric':name,'valid_positions':len(valid),'mean':float(values[valid,j].mean()),'median':float(np.median(values[valid,j])),'min':float(values[valid,j].min()),'max':float(values[valid,j].max()),'local_independent_trajectory_files':1,'label_windows':result['windows']})
 audits.append({'sample_id':sid,'source_pdb':pdb,'structure':str(structure),'structure_sha256':sha(structure),'trajectory':str(trajectory),'trajectory_sha256':sha(trajectory),'cached_label':str(lp),'cached_label_sha256':sha(lp),'metadata_sha256':sha(folder/'dynarepo_metadata.json'),'residues':st.num_residues,'chains':[{'chain':ch,'role':role(ch),'n_residues':len(st.indices_for_chain(ch)),'sequence':st.chain_sequence(ch)} for ch in st.unique_chains()],'declared_replicates':len(meta.get('mds',[])),'local_trajectory_files':[p.name for p in folder.glob('*.xtc')],'n_frames':nframes,'time_ps_start':start,'time_ps_end':end,'recomputed_windows':result['windows'],'max_cached_recomputed_value_difference':maxval,'max_cached_recomputed_variance_difference':maxvar,'masks_identical':True,'value_discrepancies_above_1e_minus5':discrepancies,'reference_scope':'Antigen-only state; no antibody dynamics labels' if sid!='A00KW' else 'Complete D44.1 antibody Fab–lysozyme complex; not HyHEL-10 mutation control'})
 print(sid,'verified',result,'max differences',maxval,maxvar,flush=True)
writecsv(ROOT/'existing_labels_per_residue.csv',allrows)
writecsv(ROOT/'chain_label_summaries.csv',summaries)
protocol={'status':'completed_with_disclosed_small_Hchi_difference','analyses':audits,'label_protocol':{'source_function':'flexppi.data.labels.compute_dynarepo_labels','function_source_sha256':sha(PROJECT/'src/flexppi/data/labels.py'),'constants_sha256':sha(PROJECT/'src/flexppi/constants.py'),'pdb_parser_sha256':sha(PROJECT/'src/flexppi/data/structures.py'),'frame_stride':1,'window_frames':1000,'window_policy':'Non-overlapping complete 1000-frame windows; terminal incomplete window omitted, matching original source implementation.','RMSF':'C-alpha positions, chainwise least-squares fitting independently per window; Angstrom.','Hchi':'All chemically available chi dihedrals jointly binned into three 120 degree bins each; entropy=-sum(p ln p), nats. Undefined or incomplete side chains masked.','variance':'Across non-overlapping windows of one local trajectory; not between independent simulations.','SASA':'Not available in cached labels and not computed.'},'runtime':{'python':sys.version,'numpy':np.__version__,'MDAnalysis':mda.__version__,'platform':platform.platform()},'limits':['Only one local XTC exists for each accession despite metadata declaring three replicas.','No new MD was performed; no claim of three independent repetitions.','A00EJ and A00HE are different deposited lysozyme states with sequence changes; not a controlled rigid-versus-flex pair and not antibody dynamics.','A00KW is D44.1 Fab–lysozyme and lacks matched mutation/maturation states; useful absolute-label screen only.','These records exist in training inventory, so model comparison requires current full-chain family-held-out checkpoint audit, handled separately.','A00HE has one Hchi site differing by0.00013417 nats on recomputation; original and new values are both retained, no claim of exact Hchi identity for that record. Cause not established; discrete-bin sensitivity is a possibility, not confirmed.','No model predictions have been run by this analysis agent.']}
(ROOT/'reference_label_audit.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
(ROOT/'README_参考标签.md').write_text('''# 本地可核验的抗体相关动态参考\n\n三个原始NPZ标签均已按项目原始标签函数从本地XTC重新计算并逐值核查，A00EJ与A00KW的Hχ逐值完全一致，三个体系的RMSF差均小于10⁻⁶ Å。A00HE有一位Hχ差0.00013417 nats，方差差0.00004995，已在审计定位并同时保留原始与重算值，未宣称其完全一致。原因尚未确定，离散箱边界敏感性仅为可能性。完整逐残基表为existing_labels_per_residue.csv，链级汇总为chain_label_summaries.csv；原始与重算最大差、文件SHA、编号映射、时间与窗口数见reference_label_audit.json。未运行新MD或修改任何原文件。\n\nA00EJ（6P4B）仅含128残基D链抗原，A00HE（1VFB）仅含128残基C链抗原，不能据此评价Fab识别环；它们也不是严格只改变柔性的配对。A00KW（1MLC）为553残基三链D44.1 Fab–溶菌酶复合物：A轻链210、B重链215、E抗原128；它不是HyHEL-10，且没有配套突变对照。\n\n每项本地只有一份trajectory.xtc。元数据写三重复，不等于本地存在三条独立轨迹。A00EJ/A00HE的五个1000帧窗口属于同一长轨迹，标签方差是窗口差异；A00KW仅一个1000帧窗口，方差为零不表示无动态不确定性。RMSF为逐链拟合后的Cα波动（Å）；Hχ为每个侧链全部χ角三箱联合状态熵（nats），不是NMR弛豫交换、热力学熵或结合能。\n\n这些记录已在训练清单出现，必须由当前完整组成链家族划分选择真正留出检查点后，才能作为模型实际输出的参考。\n''',encoding='utf-8')
print('COMPLETE',len(allrows),'residue rows',len(summaries),'chain/metric summaries')
