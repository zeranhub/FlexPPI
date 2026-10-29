# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) True-coordinate structures -> v48 figure4 production molecular PNGs -> param inherit
# (b) All10 absolute RMSF -> LineTrend/plot_trend.py -> param inherit
# (c) All5 comparable Hchi sites -> PairedBoxScatter/plot_PairedBoxScatter1.R -> param inherit
# (d) Full-scope MD/model changes -> manuscript source delta comparison -> param inherit
# RULE: no native-run panel; all complete peptide positions are retained.
import os
os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from atlas_style import *
import pandas as pd

EXP = ROOT.parent.parent
SRC = EXP / 'kras_affinity_v48_20261019' / 'figure4'
OUT = ROOT / 'KRAS_G12D'
OUT.mkdir(exist_ok=True)
source = SRC / 'source_archive' / 'source_all10_peptide_model_and_MD_changes.csv'
p = pd.read_csv(source, float_precision='round_trip').sort_values('position_a')
assert len(p) == 10 and p.position_a.tolist() == list(range(1, 11))
assert ''.join(p.aa_a) == 'VVVGAGGVGK' and ''.join(p.aa_b) == 'VVVGADGVGK'
assert p.rmsf_valid.all() and p.loc[~p.same_aa, 'position_a'].tolist() == [6]
assert np.allclose(p.rmsf_full_b-p.rmsf_full_a, p.rmsf_full_delta, atol=1e-14)
q = p[p.hchi_valid].copy()
assert q.position_a.tolist() == [1, 2, 3, 8, 10]
assert q.same_aa.all()
p['label'] = [f'{a}{int(i)}' if i != 6 else 'G6→D6' for a, i in zip(p.aa_a, p.position_a)]
fig = make_canvas(185)
heading(fig, 'a', 'KRAS G12D: antigen substitution with the same TCR and HLA', .973)
assets = []
for path, rect, title, col in [
    (SRC/'7ow5_binding.png', [.055,.765,.255,.18], 'WT · 7OW5', FIRST),
    (SRC/'7ow6_binding.png', [.330,.765,.255,.18], 'G12D · 7OW6', SECOND),
    (SRC/'peptide_superposition.png', [.625,.787,.33,.12], 'MHC-aligned peptides', GRAY),
]:
    ax, record = add_molecule(fig, path, rect, title, col)
    assets.append(record)
fig.text(.79,.765,'G6 / D6',color=SECOND,fontsize=9,ha='center')
fig.text(.79,.746,'Peptide Cα RMSD: 0.61 Å',color=GRAY,fontsize=9,ha='center')
fig.text(.065,.725,'WT: VVVGAGGVGK',color=FIRST,fontsize=9)
fig.text(.405,.725,'G12D: VVVGADGVGK',color=SECOND,fontsize=9)
fig.text(.065,.704,r'SPR $K_D$: 3 μM → 743 pM',fontsize=9)
fig.text(.62,.704,'Poole et al., 2022 · Nat. Commun.',color=GRAY,fontsize=9)

heading(fig, 'b', 'Model: backbone RMSF across the complete ten-residue peptide', .661)
legend(fig, ['WT', 'G12D'], anchor=(.965,.690))
ax = fig.add_axes([.10,.485,.856,.135])
rmsf_records = absolute_profile(ax,p,'position_a','rmsf_full_a','rmsf_full_b','RMSF (Å)',
                                label_col='label',highlights=(6,),ylim=(0,1.45))
ax.set_xlabel('Peptide position',fontsize=9,labelpad=6)
for tick in ax.get_xticklabels():
    if tick.get_text() == 'G6→D6': tick.set_color(SECOND); tick.set_fontweight('bold')

heading(fig, 'c', 'Model: side-chain entropy at every comparable peptide site', .416)
ax = fig.add_axes([.10,.29,.856,.085])
hchi_records = paired_residue_values(ax,q,'hchi_full_a','hchi_full_b',
    [f'{a}{int(i)}' for a,i in zip(q.aa_a,q.position_a)],ylim=(0,3.5),rotate=0)
fig.text(.528,.251,'G4, A5, G6→D6, G7 and G9 have no paired χ-defined comparison',
         ha='center',fontsize=9,color=GRAY)

heading(fig, 'd', 'Existing MD comparison: complete-scope changes, without residue selection', .216)
diagnostics = {}
for metric, data, rect, unit, name in [
    ('rmsf',p,[.11,.061,.34,.116],'Å','Backbone · all 10 sites'),
    ('hchi',q,[.62,.061,.34,.116],'nats','Side chains · all 5 sites'),
]:
    ax = fig.add_axes(rect)
    diag = delta_comparison(ax,data[f'{metric}_md_delta'],data[f'{metric}_full_delta'],unit=unit,title=name)
    rho = float(data[f'{metric}_md_delta'].rank().corr(data[f'{metric}_full_delta'].rank()))
    diag['Spearman_rho'] = rho
    diagnostics[metric] = diag
    ax.text(.03,.94,f'ρ = {rho:.2f}',transform=ax.transAxes,fontsize=9,va='top')
    ax.xaxis.labelpad=3
    ax.yaxis.labelpad=3

pd.DataFrame(rmsf_records).to_csv(OUT/'source_absolute_RMSF.csv',index=False)
pd.DataFrame(hchi_records).to_csv(OUT/'source_absolute_Hchi.csv',index=False)
p.to_csv(OUT/'source_complete_peptide.csv',index=False)
audit = export_audited(fig,OUT/'case_KRAS_G12D',extra={
    'case':'KRAS_G12D', 'source_file':str(source),'source_SHA256':sha256(source),
    'molecular_assets':assets,'all_peptide_positions':10,'RMSF_positions':10,
    'paired_Hchi_positions':5,'Hchi_no_pair_positions':[4,5,6,7,9],
    'Hchi_exclusion_reason':'At least one state has no side-chain chi angle; G6→D6 changes side-chain chemistry.',
    'diagnostics':diagnostics,'primary_curves':'Absolute model predictions; no MD curves',
    'model_checkpoint_SHA256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',
    'model_training':'1866 ATLAS + 564 DynaRepo; zero ITrace training systems; joint Full fold0',
    'MD_sampling':'Three independent existing 200 ns trajectories per state; source CSV has per-replica values.',
    'statistics':'Description only. Spearman over complete predefined residues; no significance test, no confidence interval. Model curves have no error bars.',
    'limitations':'Most comparable side-chain changes are weak; the substituted G6 has no paired Hchi comparison. This case primarily supports local backbone interpretation, not a joint side-chain/backbone causal mechanism.',
    'new_training_inference_MD':False,
})
(OUT/'caption_zh.md').write_text('KRAS G12D抗原替换案例。a，WT与G12D复合物的真实结构以及文献SPR结合常数；两种肽在同一MHC参考系中叠合。b，全部10个肽位置的模型RMSF绝对预测。c，两种状态均具有χ角定义的全部5个位点的模型侧链联合熵；G6→D6等缺乏配对χ定义的位置不构造熵差。d，完整预设区域中模型变化与已有三次独立轨迹均值变化的描述性比较。每点为残基而非独立实验重复，虚线表示数值一致；ρ为Spearman相关，无显著性检验。该例以主链局部变化的解释为主，不能由此建立G6侧链熵变化与邻位RMSF之间的因果机制。',encoding='utf8')
print(json.dumps({'output':str(OUT),'SHA256':audit['PNG_SHA256'],'diagnostics':diagnostics},ensure_ascii=False))
