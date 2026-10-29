from pathlib import Path
import shutil
root=Path(__file__).resolve().parent;old=root.parent/'revision_v107_soft_alpha_overlay';out=root/'figure3';out.mkdir(exist_ok=True)
shutil.copytree(old/'figure3/source_data',out/'source_data',dirs_exist_ok=True)
shutil.copytree(old/'pymol',root/'pymol',dirs_exist_ok=True)
shutil.copy2(old/'figure3/publication_style.py',out/'publication_style.py')
shutil.copy2(old/'plot_overlay.py',root/'plot_overlay.py')
source=(old/'figure3/compose_fig3_v107.py').read_text(encoding='utf-8')
source=source.replace('# b: resolve actual bound-state flexibility, with a separate static geometry track.', '''# b/c: both metrics per case; color = metric, line style = source.
import sys
sys.path.insert(0,str(ROOT))
from case_profiles import draw_case,method_handles
# Original static context tracks and all their values are preserved.
# b: dual metrics plus original shielding track.''')
start=source.index("panel('b','CaMKI backbone and shielding'")
end=source.index('# c/d: complete sequence accuracy',start)
b='''panel('b','CaMKI dynamics and shielding',98,2)
fig.legend(handles=method_handles(),loc='center',bbox_to_anchor=(145/W,1-10/H),ncol=2,handlelength=1.8,columnspacing=1.2,frameon=False)
x=cam.CaMKI_native_position.to_numpy()
a,b,profile_values=draw_case(fig,(W,H),(111,18,60,22),'CaMKI',xlabel=False)
for (metric,series),z in profile_values.groupby(['metric','series'],sort=False):
    values('b','RMSF_A' if metric=='RMSF' else 'Hchi_nats',z.residue_id,z.value,series=series)
mae_full=float(np.mean(abs(cam.pred_rmsf-cam.target_rmsf)));mae_intra=float(np.mean(abs(cam.intra_rmsf-cam.target_rmsf)))
# Identify sites along the residue axis without implying equality across y axes.
for pos in [303,316]:
    a.axvline(pos,color=TEAL,lw=.4,ls=':',alpha=.6,zorder=0)
    a.text(pos,.98,('W' if pos==303 else 'M')+str(pos),transform=a.get_xaxis_transform(),ha='center',va='top',fontsize=9,color=TEAL)
a=ax(111,44,60,9);style_axis(a,'b_static_shielding')
shield=cam.sidechain_shielding_fraction.to_numpy()*100
a.bar(x,shield,width=.78,color=TEAL,alpha=.60,edgecolor='none')
a.set(xlim=(293.5,318.5),ylim=(0,105),ylabel='CaM\\nshielding (%)',xlabel='CaMKI residue')
a.set_yticks([0,100]);a.set_xticks([294,300,306,312,318]);a.tick_params(axis='y',length=2,pad=1)
values('b','static_sidechain_shielding_percent',x,shield,series='starting_structure')

'''
source=source[:start]+b+source[end:]
start=source.index("panel('c','Ubiquitin side chains'")
end=source.index("ordered=gln.sort_values('MD_Hchi_nats')",start)
c='''panel('c','Ubiquitin backbone and side chains',3,62)
panel('d','Same side chain, different environments',98,62)
a,b,profile_values=draw_case(fig,(W,H),(14,74,69,30),'Ubiquitin')
for (metric,series),z in profile_values.groupby(['metric','series'],sort=False):
    values('c','RMSF_A' if metric=='RMSF' else 'Hchi_nats',z.residue_id,z.value,series=series)
for pos in [41,62]:a.axvspan(pos-.55,pos+.55,color=TEAL,alpha=.065,lw=0,zorder=0)
# Positional annotations sit on the shared x axis, outside the measured y values.
a.text(41,.99,'Q41',transform=a.get_xaxis_transform(),ha='center',va='top',fontsize=9,color=TEAL)
a.text(62,.99,'Q62',transform=a.get_xaxis_transform(),ha='center',va='top',fontsize=9,color=TEAL)
'''
source=source[:start]+c+source[end:]
source=source.replace("assert df.groupby('panel').size().to_dict()=={'b':75,'c':132,'d':18,'e':302}","assert df.groupby('panel').size().to_dict()=={'b':117,'c':276,'d':18,'e':302}")
source=source.replace("'panels_a_to_e_unchanged_from_v103':True","'a_d_e_and_f_unchanged_from_v107':True,'b_c_both_metrics_added':True")
(out/'compose_fig3_v108.py').write_text(source,encoding='utf-8')
print('Prepared updated b/c dual-metric profiles; a/d/e/f unchanged.')
