# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) True input structures -> render_SB27_JM22.py -> param inherit
# (b) Complete13 peptide RMSF -> LineTrend/plot_trend.py -> param inherit
# (c) Complete11 valid Hchi -> PairedBoxScatter/plot_PairedBoxScatter1.R -> param inherit
# (d) Q7 input control and complete Hchi diagnostic -> case source CSV -> param inherit
# RULE: all13 positions and all11 valid Hchi pairs retained; no native-run panels.
import os
os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from atlas_style import *
import pandas as pd
EXP=ROOT.parent.parent
OUT=ROOT/'SB27_SB47'
SRC=EXP/'tcr_case_v40_20261017'
source=SRC/'figure6'/'all13_model_and_MD_changes.csv'
ctxsource=SRC/'sb27_sb47'/'context_inference.csv'
repsource=SRC/'sb27_sb47'/'Q7_per_replica_physical_summary.csv'
p=pd.read_csv(source,float_precision='round_trip').sort_values('position_a')
ctx=pd.read_csv(ctxsource,float_precision='round_trip')
rep=pd.read_csv(repsource,float_precision='round_trip')
assert len(p)==13 and p.position_a.tolist()==list(range(1,14))
assert ''.join(p.aa_a)=='LPEPLPQGQLTAY' and p.same_aa.all() and p.rmsf_valid.all()
q=p[p.hchi_valid].copy()
assert len(q)==11 and p.loc[~p.hchi_valid,'position_a'].tolist()==[8,12]
for metric in ['rmsf','hchi']:
    v=p[p[f'{metric}_valid']]
    assert np.allclose(v[f'{metric}_full_b']-v[f'{metric}_full_a'],v[f'{metric}_full_delta'],atol=1e-14)
p['label']=[f'{a}{int(i)}' for a,i in zip(p.aa_a,p.position_a)]
fig=make_canvas(190)
heading(fig,'a','SB27 and SB47: distinct TCR docking on the same bulged peptide–HLA',.973)
assets=[]
for path,rect,title,col in [
    (OUT/'2ak4_binding.png',[.055,.765,.255,.18],'SB27 · 2AK4',FIRST),
    (OUT/'4jry_binding.png',[.330,.765,.255,.18],'SB47 · 4JRY',SECOND),
    (OUT/'peptide_superposition.png',[.635,.793,.32,.115],'MHC-aligned peptide',GRAY),
]:
    ax,record=add_molecule(fig,path,rect,title,col);assets.append(record)
fig.text(.795,.765,'Q7',color=HLA,ha='center',fontsize=9,weight='bold')
fig.text(.795,.744,'Cα RMSD: 0.92 Å',color=GRAY,ha='center',fontsize=9)
fig.text(.065,.725,'HLA-B*35:08 · LPEPLPQGQLTAY',fontsize=9)
fig.text(.065,.704,'Comparison: receptor identity and docking geometry',color=GRAY,fontsize=9)

heading(fig,'b','Model: backbone RMSF across the complete thirteen-residue peptide',.661)
legend(fig,['SB27','SB47'],anchor=(.965,.69))
ax=fig.add_axes([.10,.485,.856,.135])
rrecords=absolute_profile(ax,p,'position_a','rmsf_full_a','rmsf_full_b','RMSF (Å)',
    label_col='label',highlights=(7,),ylim=(0,1.8))
ax.set_xlabel('Peptide position',fontsize=9,labelpad=6)
for tick in ax.get_xticklabels():
    if tick.get_text()=='Q7':tick.set_color(HLA);tick.set_fontweight('bold')

heading(fig,'c','Model: side-chain entropy at all eleven χ-defined peptide sites',.416)
ax=fig.add_axes([.10,.291,.856,.085])
hrecords=paired_residue_values(ax,q,'hchi_full_a','hchi_full_b',
    [f'{a}{int(i)}' for a,i in zip(q.aa_a,q.position_a)],id_col='position_a',highlight_ids=(7,),
    ylim=(-.15,3.5),rotate=0)
for tick in ax.get_xticklabels():
    if tick.get_text()=='Q7':tick.set_color(HLA);tick.set_fontweight('bold')
fig.text(.528,.251,'G8 and A12 have no side-chain χ angles',ha='center',fontsize=9,color=GRAY)

heading(fig,'d','Q7 context dependence and the scope of the side-chain result',.216)
ax=fig.add_axes([.115,.061,.33,.113])
context_records=[]
for pdb,col,marker,face in [('2ak4',FIRST,'o','white'),('4jry',SECOND,'s',SECOND)]:
    t=ctx[(ctx.pdb==pdb)&(ctx.position==7)].set_index('input').loc[['full','pMHC_only']]
    assert len(t)==2
    ax.plot([0,1],t.Hchi,color=col,marker=marker,mfc=face,mec=col,lw=1.05,ms=4)
    for name,row in t.iterrows():context_records.append({'pdb':pdb,'input':name,'Q7_model_Hchi':row.Hchi})
ax.set(xlim=(-.15,1.15),ylim=(1.3,3.15),xticks=[0,1],xticklabels=['Full complex','TCR removed'],
       ylabel=r'Q7 model $H_\chi$ (nats)',yticks=[1.5,2,2.5,3])
ax.tick_params(length=2.5,pad=3)
ax.set_title('Model input control',fontsize=9,pad=5)
ax=fig.add_axes([.635,.061,.325,.113])
diag=delta_comparison(ax,q.hchi_md_delta,q.hchi_full_delta,unit='nats',title='All 11 side-chain changes')
rho=float(q.hchi_md_delta.rank().corr(q.hchi_full_delta.rank()))
diag['Spearman_rho']=rho
ax.text(.03,.94,f'ρ = {rho:.2f}',transform=ax.transAxes,fontsize=9,va='top')
for pos in [7,9]:
    row=q[q.position_a==pos].iloc[0]
    ax.annotate(f'Q{pos}',(row.hchi_md_delta,row.hchi_full_delta),xytext=(4,4),textcoords='offset points',fontsize=9)

pd.DataFrame(rrecords).to_csv(OUT/'source_absolute_RMSF.csv',index=False)
pd.DataFrame(hrecords).to_csv(OUT/'source_absolute_Hchi.csv',index=False)
pd.DataFrame(context_records).to_csv(OUT/'source_Q7_context.csv',index=False)
p.to_csv(OUT/'source_complete_peptide.csv',index=False)
rep.to_csv(OUT/'source_Q7_existing_MD_replicas.csv',index=False)
rhor=float(p.rmsf_md_delta.rank().corr(p.rmsf_full_delta.rank()))
audit=export_audited(fig,OUT/'case_SB27_SB47',extra={
    'case':'SB27_SB47','all_peptide_positions':13,'RMSF_positions':13,'Hchi_positions':11,
    'Hchi_invalid_positions':[8,12],'Hchi_invalid_reason':'Gly and Ala lack side-chain chi angles',
    'source_files':{str(v):sha256(v) for v in [source,ctxsource,repsource,OUT/'molecular_view.json']},
    'molecular_assets':assets,'MD_Hchi_diagnostic':diag,'MD_RMSF_change_Spearman_all13':rhor,
    'model_checkpoint_SHA256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',
    'model_training':'1866 ATLAS + 564 DynaRepo; zero ITrace training systems; joint Full fold0',
    'context_control':'Remove TCR atoms while preserving each own-state pMHC coordinates; model-input dependence only, not an unbinding or apo experiment.',
    'MD_sampling':'Three independent existing 200 ns trajectories per state; all saved frames.',
    'statistics':'No hypothesis test or confidence interval. Spearman over residue positions; points are not experimental replicas. Single model inference has no error bars.',
    'limitations':'Q7 is a strong local entropy/context example, but whole-peptide Hchi change agreement is modest and Q9 is retained. Different TCR identities and docking cannot isolate a single causal mutation.',
    'primary_curves':'Absolute model output only; no MD curves',
    'new_training_inference_MD':False,
})
(OUT/'caption_zh.md').write_text('SB27/SB47对照。a，两种TCR识别同一HLA-B*35:08呈递的13肽LPEPLPQGQLTAY；真实模型输入结构按MHC对齐，肽独立放大，蓝绿色球标记Q7的Cα。b，全部13个位点模型RMSF。c，全部11个有效χ位点模型Hχ，Q7以底色标明，G8/A12没有侧链χ角。d，左图保持各自pMHC坐标而移除TCR输入，展示Q7读数对结合伙伴信息的依赖；这不是解离实验。右图保留全部11个有效位点的模型变化与已有三次200 ns轨迹均值变化比较；Q9同时保留。模型主链变化的全13位点相关为0.945，而侧链全11位点相关约0.321，表明Q7局部解释不能推广为整个肽侧链变化均被准确捕获。MD只提供独立佐证，不作为主曲线。',encoding='utf8')
print(json.dumps({'output':str(OUT),'PNG_SHA256':audit['PNG_SHA256'],'rmsf_rho':rhor,'Hchi_rho':rho},ensure_ascii=False))
