# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) True input structures -> render_SB27_JM22.py -> param inherit
# (b) Complete9 absolute RMSF -> LineTrend/plot_trend.py -> param inherit
# (c) Complete7 valid Hchi -> PairedBoxScatter/plot_PairedBoxScatter1.R -> param inherit
# (d) Complete RMSF/Hchi diagnostic -> validated source CSV -> param inherit
# RULE: no native-run panels; retain all9 positions and flag F5L chemistry change.
import os
os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from atlas_style import *
import pandas as pd
EXP=ROOT.parent.parent
OUT=ROOT/'JM22_F5L'
SRC=EXP/'tcr_expanded_screen_v41_20261017'/'jm22_f5l_validation'
source=SRC/'all9_model_MD_physical_changes.csv'
p=pd.read_csv(source,float_precision='round_trip').sort_values('position')
assert len(p)==9 and p.position.tolist()==list(range(1,10))
assert ''.join(p.aa_a)=='GILGFVFTL' and ''.join(p.aa_b)=='GILGLVFTL'
assert p.rmsf_valid_a.all() and p.rmsf_valid_b.all()
q=p[p.hchi_valid_a&p.hchi_valid_b].copy()
assert len(q)==7 and p.loc[~(p.hchi_valid_a&p.hchi_valid_b),'position'].tolist()==[1,4]
assert p.loc[~p.same_aa,'position'].tolist()==[5]
q['altered_chemistry']=~q.same_aa
for metric in ['rmsf','hchi']:
    v=q if metric=='hchi' else p
    assert np.allclose(v[f'{metric}_prediction_b']-v[f'{metric}_prediction_a'],v[f'{metric}_prediction_delta'],atol=1e-14)
p['label']=[f'{a}{int(i)}' if i!=5 else 'F5→L5' for a,i in zip(p.aa_a,p.position)]
fig=make_canvas(185)
heading(fig,'a','JM22: recognition of the influenza M1 epitope and its F5L variant',.973)
assets=[]
for path,rect,title,col in [
    (OUT/'1oga_binding.png',[.055,.765,.255,.18],'WT · 1OGA',FIRST),
    (OUT/'5hhm_binding.png',[.330,.765,.255,.18],'F5L · 5HHM',SECOND),
    (OUT/'peptide_superposition.png',[.635,.793,.32,.115],'MHC-aligned peptides',GRAY),
]:
    ax,record=add_molecule(fig,path,rect,title,col);assets.append(record)
fig.text(.795,.765,'F5 / L5',color=SECOND,ha='center',fontsize=9,weight='bold')
fig.text(.795,.744,'Cα RMSD: 0.51 Å',color=GRAY,ha='center',fontsize=9)
fig.text(.065,.725,'HLA-A*02:01 · GILGFVFTL → GILGLVFTL',fontsize=9)
fig.text(.065,.704,'Natural M1 variant · Valkenburg et al., 2016 · PNAS',color=GRAY,fontsize=9)

heading(fig,'b','Model: backbone RMSF across the complete nine-residue peptide',.661)
legend(fig,['WT','F5L'],anchor=(.965,.69))
ax=fig.add_axes([.10,.485,.856,.135])
rrecords=absolute_profile(ax,p,'position','rmsf_prediction_a','rmsf_prediction_b','RMSF (Å)',
    label_col='label',highlights=(5,),ylim=(0,1.45))
ax.set_xlabel('Peptide position',fontsize=9,labelpad=6)
for tick in ax.get_xticklabels():
    if tick.get_text()=='F5→L5':tick.set_color(SECOND);tick.set_fontweight('bold')

heading(fig,'c','Model: side-chain entropy at all seven χ-defined peptide sites',.416)
ax=fig.add_axes([.10,.291,.856,.085])
hrecords=paired_residue_values(ax,q,'hchi_prediction_a','hchi_prediction_b',
    [f'{a}{int(i)}' if i!=5 else 'F5→L5' for a,i in zip(q.aa_a,q.position)],
    id_col='position',highlight_ids=(5,),mutation_col='altered_chemistry',ylim=(0,1.9),rotate=0)
for tick in ax.get_xticklabels():
    if tick.get_text()=='F5→L5':tick.set_color(SECOND);tick.set_fontweight('bold')
fig.text(.528,.251,'▼ Altered side-chain chemistry at position 5; G1 and G4 have no χ angles',
         ha='center',fontsize=9,color=GRAY)

heading(fig,'d','Existing MD comparison: backbone and side-chain changes',.216)
diagnostics={}
for metric,data,rect,unit,title in [
    ('rmsf',p,[.11,.061,.34,.116],'Å','Backbone · all 9 sites'),
    ('hchi',q,[.62,.061,.34,.116],'nats','Side chains · all 7 sites'),
]:
    ax=fig.add_axes(rect)
    diag=delta_comparison(ax,data[f'{metric}_MD_delta'],data[f'{metric}_prediction_delta'],unit=unit,title=title)
    rho=float(data[f'{metric}_MD_delta'].rank().corr(data[f'{metric}_prediction_delta'].rank()))
    common=data[data.same_aa]
    diag.update({'Spearman_rho':rho,'same_AA_n':int(len(common)),
        'same_AA_Spearman_rho':float(common[f'{metric}_MD_delta'].rank().corr(common[f'{metric}_prediction_delta'].rank())),
        'same_AA_MAE':float(np.mean(np.abs(common[f'{metric}_MD_delta']-common[f'{metric}_prediction_delta'])))})
    diagnostics[metric]=diag
    ax.text(.03,.94,f'ρ = {rho:.2f}',transform=ax.transAxes,fontsize=9,va='top')
    if metric=='hchi':
        row=data[data.position==5].iloc[0]
        ax.scatter([row.hchi_MD_delta],[row.hchi_prediction_delta],s=32,facecolors='white',
            edgecolors=SECOND,linewidths=1,marker='v',zorder=5)
        ax.annotate('F5→L5',(row.hchi_MD_delta,row.hchi_prediction_delta),xytext=(-4,-16),
            ha='right',textcoords='offset points',fontsize=9,color=SECOND)

pd.DataFrame(rrecords).to_csv(OUT/'source_absolute_RMSF.csv',index=False)
pd.DataFrame(hrecords).to_csv(OUT/'source_absolute_Hchi.csv',index=False)
p.to_csv(OUT/'source_complete_peptide.csv',index=False)
audit=export_audited(fig,OUT/'case_JM22_F5L',extra={
    'case':'JM22_F5L','all_peptide_positions':9,'RMSF_positions':9,'Hchi_positions':7,
    'Hchi_invalid_positions':[1,4],'Hchi_invalid_reason':'Gly lacks side-chain chi angles',
    'substituted_position':5,'substitution':'Phe to Leu; changes side-chain chemistry, not treated as same-AA entropy constraint',
    'source_files':{str(v):sha256(v) for v in [source,SRC/'raw_trajectory_audit.json',OUT/'molecular_view.json']},
    'molecular_assets':assets,'diagnostics':diagnostics,
    'model_checkpoint_SHA256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',
    'model_training':'1866 ATLAS + 564 DynaRepo; zero ITrace training systems; joint Full fold0',
    'MD_sampling':'Three independent existing 200 ns trajectories per state; all saved frames.',
    'statistics':'No hypothesis test or confidence interval. Spearman over residue positions; points are not experimental replicas. Single model inference has no error bars.',
    'limitations':'The large site5 Hchi difference includes a Phe-to-Leu chemistry change. Smaller changes at six unchanged chi-defined sites and backbone readouts cannot establish a directed inter-residue pathway or altered binding affinity.',
    'primary_curves':'Absolute model output only; no MD curves','new_training_inference_MD':False,
})
(OUT/'caption_zh.md').write_text('JM22识别流感M1表位及F5L天然变体。a，WT与F5L的真实模型输入结构，经MHC对齐后比较肽构象；蓝绿色球标记第5位Cα。b，完整9肽的模型RMSF绝对读数。c，全部7个χ角定义位点的模型Hχ，倒三角单列发生Phe→Leu替换的位置；G1和G4无侧链χ角。d，保留全部9个RMSF位点及7个Hχ位点的模型变化与已有三次独立200 ns轨迹均值变化比较；侧链图中倒三角表示第5位，虚线表示数值一致。全部位点的RMSF与Hχ变化相关分别为0.750和0.786，去除替换位点后分别为0.714（8位点）和0.657（6个有效侧链）。该例提示侧链读数能够提供主链之外的变化线索，但第5位熵差不能脱离氨基酸化学身份变化解释，不能单独证明跨残基因果路径。',encoding='utf8')
print(json.dumps({'output':str(OUT),'PNG_SHA256':audit['PNG_SHA256'],'diagnostics':diagnostics},ensure_ascii=False))
