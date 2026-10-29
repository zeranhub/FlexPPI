from pathlib import Path
import shutil
root=Path(__file__).resolve().parent;old=root.parent/'revision_v103_direct_RMSF_Hchi';out=root/'figure3';out.mkdir(exist_ok=True)
shutil.copytree(old/'figure3'/'source_data',out/'source_data',dirs_exist_ok=True)
shutil.copytree(old/'pymol',root/'pymol',dirs_exist_ok=True)
shutil.copy2(old/'figure3'/'publication_style.py',out/'publication_style.py')
source=(old/'figure3'/'compose_fig3_v103.py').read_text(encoding='utf-8')
prefix=source.split('# f: use the user-requested direct MD and predicted joint distributions.')[0]
prefix=prefix.replace('W,H=183.,227.','W,H=183.,230.')
footer='''# f: one-axis overlay, using the identical common MD/Full cohort.
import sys
sys.path.insert(0,str(ROOT))
from plot_overlay import draw_overlay,load_data
panel('f','RMSF–Hχ relationship',3,173)
density_audit=draw_overlay(fig,(W,H),(14,190,157,30),xmax=10,gridsize=(80,35),marginal_mm=4.)
direct_arrays,direct_stats=load_data()
df=pd.DataFrame(plotted);df.to_csv(OUT/'main_original_panels_a_to_e_values.csv',index=False)
assert df.groupby('panel').size().to_dict()=={'b':75,'c':132,'d':18,'e':302}
fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
for t in fig.findobj(mpl.text.Text):
    if t.get_visible() and t.get_text():
        bb=t.get_window_extent(renderer)
        if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'fig3_dual_case_entropy.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
qa={'size_mm':[W,H],'font_pt':9,'panels_a_to_e_unchanged_from_v103':True,'e_four_horizontal':True,
 'f':'single-axis overlay of MD and Full RMSF-Hchi densities and marginals','common_f_positions':824988,
 'f_records':3063,'density':density_audit,'outside_text':outside,'model_and_MD_rerun':False,'word_changed':False,
 'visual_review':'pending','sha256':{e:hashlib.sha256((OUT/f'fig3_dual_case_entropy.{e}').read_bytes()).hexdigest() for e in ['png','pdf','svg']}}
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
plt.close(fig)
print(json.dumps(qa,ensure_ascii=False))
'''
(out/'compose_fig3_v107.py').write_text(prefix+footer,encoding='utf-8')
print('Prepared one-overlay f main figure; a-e unchanged.')
