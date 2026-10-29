"""Short read-only checks/summaries of completed contacts. No trajectory loop or MD."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,sys
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
CASE=BASE/'free_energy_curves_v27_20261014/case_curves'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def js(path,value):Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def main():
    sys.stdout.reconfigure(encoding='utf-8')
    checks=[]
    def check(name,value,detail=None):
        row={'name':name,'passed':bool(value)}
        if detail is not None:row['detail']=detail
        checks.append(row)
        assert row['passed'],row
    old=json.loads((HERE/'run_audit.json').read_text(encoding='utf-8'))
    frozen={str(HERE/'run_audit.json'):sha(HERE/'run_audit.json')}
    for path,digest in old['sources_sha256'].items():
        check('SourceSHA '+Path(path).name,sha(path)==digest)
        frozen[path]=digest
    for name,meta in old['outputs'].items():
        check('ExistingcompletedoutputSHA '+name,sha(HERE/name)==meta['sha256'])
        frozen[str(HERE/name)]=meta['sha256']
    check('FrozenanalysisplanSHA',sha(HERE/'analysis_plan.json')==old['plan_sha256'])
    frozen[str(HERE/'analysis_plan.json')]=old['plan_sha256']
    frames=pd.read_csv(HERE/'trajectory_contact_timeseries_all5000.csv',float_precision='round_trip')
    primary=pd.read_csv(HERE/'trajectory_contact_timeseries_primary501.csv',float_precision='round_trip')
    summary=pd.read_csv(HERE/'contact_summary_all6_by_grid_window.csv',float_precision='round_trip')
    residue=pd.read_csv(HERE/'residue_contact_persistence.csv',float_precision='round_trip')
    polar=pd.read_csv(HERE/'polar_atom_pair_persistence.csv',float_precision='round_trip')
    sasa=pd.read_csv(HERE/'static_SASA_all6.csv',float_precision='round_trip')
    static=pd.read_csv(HERE/'static_contacts_by_residue.csv',float_precision='round_trip')
    inventory=pd.read_csv(CASE/'selection_inventory_with_MD_outcomes.csv',float_precision='round_trip')
    entropy=pd.read_csv(CASE/'joint_entropy_original_five_windows.csv',float_precision='round_trip')
    for path in [CASE/'selection_inventory_with_MD_outcomes.csv',CASE/'joint_entropy_original_five_windows.csv']:
        frozen[str(path)]=sha(path)
    sites=['Q2','Q31','Q40','Q41','Q49','Q62']
    check('All6sites5000frames',frames.groupby('site').size().to_dict()==dict.fromkeys(sites,5000))
    check('All6sitesprimary501frames',primary.groupby('site').size().to_dict()==dict.fromkeys(sites,501))
    for site,g in frames.groupby('site'):
        check('Alloriginalindices '+site,np.array_equal(g.frame_zero_based,np.arange(5000)))
        check('Actualtimestamps0to499p9 '+site,np.max(abs(g.time_ns-g.frame_zero_based/10))<1e-8)
        check('Fivecontiguous1000observations '+site,g.groupby('original_window').size().to_dict()=={i:1000 for i in range(1,6)})
    a=frames[frames.in_primary501].reset_index(drop=True)
    check('Primary501exactsubsetoffull5000',a.equals(primary))
    avgcols=[x for x in frames if x.endswith('_any') or x.endswith('_count')]
    largest_error=0.
    for r in summary.itertuples(index=False):
        part=frames[frames.site.eq(r.site)]
        if r.grid=='primary501':part=part[part.in_primary501]
        elif r.grid.startswith('window'):part=part[part.original_window.eq(int(r.grid[-1]))]
        check('Summaryrecorddenominator '+r.grid+'/'+r.site,len(part)==r.frames)
        for col in avgcols:
            err=abs(part[col].mean()-getattr(r,col+'_mean'));largest_error=max(largest_error,err)
    check('Allsummariesrecalculatedfromsavedframes',largest_error<1e-14,largest_error)
    for data,keys in [(residue,['site','partner_chain','partner_resid']),(polar,['site','source_atom_name','target_chain','target_resid','target_atom_name'])]:
        whole=data[data.grid.eq('original5000')].set_index(keys)
        windows=data[data.grid.str.startswith('window')].groupby(keys).contact_frames.sum().sort_index()
        check('Wholecountsarefivewindows '+str(len(keys)),whole.contact_frames.sort_index().equals(windows))
        check('Occupancyexactintegercount '+str(len(keys)),np.allclose(data.contact_fraction*data.frames,data.contact_frames,atol=1e-10,rtol=0))
    # Independent direct sourcePDB geometry; does not rely on the original Universe's firstframe state.
    pdb=Path(next(k for k in old['sources_sha256'] if k.endswith('structure.pdb')))
    atoms=[]
    for line in pdb.read_text().splitlines():
        if line[:6].strip() not in ['ATOM','HETATM']:continue
        element=line[76:78].strip().upper()
        atoms.append({'atom_index':len(atoms),'chain':line[21].strip(),'resid':int(line[22:26]),'resname':line[17:20].strip(),
                      'name':line[12:16].strip(),'element':element,'xyz':np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])})
    heavy=[x for x in atoms if x['element'] not in ['H','D']]
    check('DirectsourcePDBheavy935',len(atoms)==1921 and len(heavy)==935)
    atomwise=[];errors=[]
    side=('CB','CG','CD','OE1','NE2');back=('N','CA','C','O')
    for r in static.itertuples(index=False):
        own=[x for x in heavy if x['chain']=='A' and x['resid']==int(r.site[1:]) and x['name'] in (side if r.source_group=='sidechain' else back)]
        other=[x for x in heavy if x['chain']==r.partner_chain and x['resid']==r.partner_resid]
        d=np.sqrt(((np.array([x['xyz'] for x in own])[:,None,:]-np.array([x['xyz'] for x in other])[None,:,:])**2).sum(axis=2))
        errors.append(abs(float(d.min())-r.minimum_distance_angstrom))
        if d.min()<=5. and r.eligible_nonadjacent_4A:
            ii,jj=np.unravel_index(np.argmin(d),d.shape)
            atomwise.append({'site':r.site,'source_group':r.source_group,'source_atom':own[ii]['name'],
               'target_chain':r.partner_chain,'target_resid':r.partner_resid,'target_resname':r.partner_resname,'target_atom':other[jj]['name'],
               'distance_A':float(d[ii,jj])})
    check('All1416staticdistancesindependentlydirectPDBrecomputed',len(errors)==1416 and max(errors)<2e-5,{'rows':len(errors),'maximum_difference_A':max(errors)})
    pd.DataFrame(atomwise).to_csv(HERE/'static_atomwise_nearest_pairs_within5A.csv',index=False)
    base=summary[summary.grid.eq('original5000')].merge(sasa,on='site').merge(inventory,left_on='site',right_on='residue',validate='one_to_one')
    columns=['site','sidechain_SASA_A2','nonadjacent_4A_samechain_residue_count_mean','nonadjacent_4A_interchain_any_mean',
      'polar3p5_samechain_any_mean','polar3p5_nonadjacent_samechain_any_mean','polar3p5_interchain_any_mean',
      'predicted_joint_Hchi_nats','MD_joint_Hchi_mean5windows_nats','prediction_absolute_error_nats','MD_width_501frames_deg','MD_width_5000frames_deg']
    base=base[columns].set_index('site').loc[sites].reset_index()
    base.to_csv(HERE/'all6_structure_contact_comparison.csv',index=False)
    sw=summary[summary.grid.str.startswith('window')].copy()
    sw['window']=sw.grid.str[-1].astype(int)
    sw=sw.merge(entropy,left_on=['site','window'],right_on=['residue','window'],validate='one_to_one')
    check('Joinedall6all5windows',len(sw)==30)
    sw.to_csv(HERE/'all6_window_contacts_entropy.csv',index=False)
    key=['nonadjacent_4A_samechain_residue_count_mean','polar3p5_samechain_any_mean','polar3p5_nonadjacent_samechain_any_mean','Hchi_nats']
    changes=[]
    for site in sites:
        w4=sw[sw.site.eq(site)&sw.window.eq(4)].iloc[0]
        w5=sw[sw.site.eq(site)&sw.window.eq(5)].iloc[0]
        row={'site':site}
        for k in key:row[k+'_window4']=w4[k];row[k+'_window5']=w5[k];row[k+'_delta5minus4']=w5[k]-w4[k]
        changes.append(row)
    pd.DataFrame(changes).to_csv(HERE/'window4_to5_all6.csv',index=False)
    selected=polar[polar.site.eq('Q40')&(((polar.source_atom_name.eq('NE2'))&polar.target_resid.eq(72)&polar.target_atom_name.isin(['OT1','OT2']))|
      ((polar.source_atom_name.eq('OE1'))&polar.target_resid.eq(39)&polar.target_atom_name.eq('O'))|
      ((polar.source_atom_name.eq('NE2'))&polar.target_resid.eq(39)&polar.target_atom_name.eq('O')))].copy()
    selected['interpretation']=np.where(selected.target_atom_name.isin(['OT1','OT2']),'A_chain_terminal_carboxyl_oxygen_not_ARG_guanidinium',
      np.where(selected.source_atom_name.eq('OE1'),'O_to_O_proximity_not_hydrogenbond_evidence','N_to_O_candidate_only_no_geometry_test'))
    selected.to_csv(HERE/'Q40_key_polar_contacts_all_windows.csv',index=False)
    topfiles=[]
    for grid in ['original5000','primary501']:
        x=residue[residue.grid.eq(grid)&residue.contact_fraction.gt(0)].sort_values(['site','contact_fraction','partner_chain','partner_resid'],ascending=[True,False,True,True]).groupby('site',sort=False).head(10)
        name='top_residue_contacts_'+grid+'.csv';x.to_csv(HERE/name,index=False);topfiles.append(name)
        x=polar[polar.grid.eq(grid)&polar.opposite_N_O_elements&polar.contact_fraction.gt(0)].sort_values(['site','contact_fraction','source_atom_name','target_chain','target_resid','target_atom_name'],ascending=[True,False,True,True,True,True]).groupby('site',sort=False).head(10).copy()
        x['target_category']=np.where(x.target_atom_name.isin(['OT1','OT2','OXT']),'terminal_carboxyl_oxygen',np.where(x.target_atom_name.isin(['N','CA','C','O']),'backbone','sidechain'))
        name='top_opposite_NO_contacts_'+grid+'.csv';x.to_csv(HERE/name,index=False);topfiles.append(name)
    text=(
      '现有结构与500 ns单条MD轨迹不能把Q40的较大预测误差归因于某个独特、持续的极性锁定网络。Q40输入侧链SASA为70.636 Å²，与Q2的70.659及Q49的69.100接近；Q41为0，Q31为80.598，Q62为99.762。Q40初始侧链在4 Å内的非邻位伙伴主要是Pro37，最近重原子距离3.437 Å；OE1/NE2初始没有到其他蛋白N/O的3.5 Å近接。5000帧内，Q40与UBA链B的4 Å侧链接触及3.5 Å极性近接均为零，最近链B重原子距离仍为5.743 Å。因此，不能用UBA界面约束或高度埋藏解释Q40。\n\n'
      'Q40确有局部同链接触：侧链到Pro37和Arg72的4 Å接触分别出现在87.58%和74.64%的保存帧；这一判据包括CB，初始Pro37最近接触实际为Q40 CB–Pro37 O，不能据此认定末端酰胺被锁定。最常见的N/O候选近接是Q40 NE2到A链末端Arg72的OT1/OT2，分别为23.82%和22.32%。这两个原子属于输入链末端羧基，并非Arg侧链胍基，不能称为GLN–ARG盐桥或已经确认的氢键。Q40侧链末端OE1/NE2与任意同链N/O近接的帧比例为52.84%，排除直接序列邻位后为42.08%；Q41的对应非邻位比例为99.62%，NE2–Ile36 O单对为96.60%，明显更持久。六个GLN并排结果表明Q40不是最埋藏或具有最持久极性网络的位点；接触数量或单个距离阈值不足以解释其误差排序。\n\n'
      '末窗变化也不支持“关键远程极性接触解除导致Hχ升高”的简单解释。Q40第4至第5个100 ns窗口的Hχ由1.956升至2.557 nats；任意同链N/O近接由48.6%降至36.8%，但排除邻位后反而由27.8%升至32.4%。NE2–OT1/OT2分别由15.1%/13.8%升至20.3%/19.3%；非邻位4 Å伙伴数也由1.877增至1.923。下降最明显的局部项之一是Q40 OE1与邻位Asp39主链O的O–O近接，18.3%降至6.1%，两者均为氧，不能把它解释为氢键断裂。因此，接触结果最多支持Q40周围存在随时间重排的局部接触；尚不能确定这些接触导致较窄χ分布或模型高估，更不能排除角度联合占据、分区定义和单轨迹采样的影响。\n\n'
      '可向用户概括：Q40比其他GLN难预测，当前能确认的是其五个窗口的侧链状态统计变化较大，Full预测2.554 nats而五窗参考均值为2.186 nats，绝对误差0.368 nats；结构接触没有给出一个独特且与末窗变化一致的“锁”。它并非被UBA抓住或深埋在核心，较强的接触主要来自同链Pro37及链末端，且关键远程极性近接在末窗并未减少。具体误差机制仍应保持为待验证假设。\n\n'
      '计算口径：完整A链72+B链47，935个蛋白重原子；侧链使用CB/CG/CD/OE1/NE2，4 Å接触排除自身及同链直接相邻残基。极性近接仅用OE1/NE2与其他蛋白N/O的3.5 Å距离，单独标注相邻残基和同元素对，忽略氢及角度。主采样为0–100 ns、每200 ps共501帧；时间诊断另用原标签前5000帧，0–499.9 ns每100 ps，并严格分为五个1000帧窗口。接触采用实际周期盒的最小镜像；静态SASA为原输入PDB的Biopython Shrake–Rupley、1.4 Å探针、960球点，不能称为MD平均暴露。所有时间帧属于同一轨迹，不是独立重复。\n')
    (HERE/'conclusion_zh.md').write_text(text,encoding='utf-8')
    methods=(
      '复算命令（只读已有结果的短核验）：bundled python执行本目录summarize_and_verify_existing.py。原计算脚本analyze_gln_contacts.py已于2026-10-14完成，CPU总耗时21.680秒，不需重跑。它只对经过SHA验证的临时PDB/XTC副本读取5000帧，原轨迹源未改；_staging目录只作临时读取。\n\n'
      '原始完整结果：trajectory_contact_timeseries_all5000.csv（30000位点×帧记录）、trajectory_contact_timeseries_primary501.csv（3006条）、residue_contact_persistence.csv（所有合资格伙伴均保留，含零占据）、polar_atom_pair_persistence.csv（完整5000帧内曾发生过近接的全部原子对，按全长/主窗口/五窗提供占据，未发生者为零且不在此稀疏表）。contact_summary_all6_by_grid_window.csv含全部六位点七种网格。静态环境在static_contacts_by_residue.csv，SASA在static_SASA_all6.csv。\n\n'
      '本轮追加all6_structure_contact_comparison.csv、all6_window_contacts_entropy.csv、window4_to5_all6.csv以及按相同规则排序的六位点topcontact表；Q40_key_polar_contacts_all_windows.csv仅是解释所需的透明摘录，完整原子对表不变。原通用target_backbone标志未把OT1/OT2归到末端主链羧基，本轮top表新增target_category准确标注，距离和计数未改。\n\n'
      '必要核验：原输出/脚本/plan和所有源文件SHA全部重查；从30000条帧记录重算全部摘要与501子集；验证五窗接触计数相加等于全长计数；从原始PDB文本独立重算全部1416静态残基距离，而不依赖原Universe初始化时的坐标状态。所有新文件仍只写structure目录。\n')
    (HERE/'README_结果与复算.md').write_text(methods,encoding='utf-8')
    for path,digest in frozen.items():check('Unchangedaftershortreview '+Path(path).name,sha(path)==digest)
    outputs=[p for p in HERE.iterdir() if p.is_file() and p.name not in set(old['outputs'])|{'run_audit.json','analysis_plan.json','short_verification_and_summary_audit.json'}]
    result={'status':'passed','timestamp_utc':datetime.now(timezone.utc).isoformat(),'original_run_audit_sha256':sha(HERE/'run_audit.json'),
      'sources_frozen_sha256':frozen,'checks':checks,'check_counts':{'total':len(checks),'passed':sum(x['passed'] for x in checks)},
      'new_MD':False,'original_trajectory_reread':False,'original_source_files_modified':False,'all6_and_all5windows_retained':True,
      'conclusion':'Q40 has intermittent local samechain contacts, but no unique persistent polar lock or loss-of-remote-contact explanation for lastwindow Hchi increase. No causal attribution.',
      'new_outputs':{p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in outputs}}
    js(HERE/'short_verification_and_summary_audit.json',result)
    print(json.dumps({'status':'passed','checks':len(checks),'original_run_audit_sha256':result['original_run_audit_sha256'],
                     'independent_static_max_difference_A':max(errors),'output_directory':str(HERE)},ensure_ascii=False))

if __name__=='__main__':main()
