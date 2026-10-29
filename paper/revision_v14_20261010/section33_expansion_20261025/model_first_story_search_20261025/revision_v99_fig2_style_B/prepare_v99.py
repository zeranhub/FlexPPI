from pathlib import Path
import hashlib
import shutil
import json

root=Path(__file__).resolve().parent
old=root.parent/'revision_v98_structure_fill'
base=root.parents[2]
out=root/'figure3'
out.mkdir(parents=True,exist_ok=True)
shutil.copytree(old/'figure3'/'source_data',out/'source_data',dirs_exist_ok=True)
shutil.copytree(old/'pymol',root/'pymol',dirs_exist_ok=True)
style=base/'figure2_dense_20261024'/'source_internal'/'publication_style.py'
shutil.copy2(style,out/'publication_style.py')
archive=root/'source_archive'
archive.mkdir(exist_ok=True)
shutil.copy2(old/'figure3'/'compose_fig3_v98.py',archive/'compose_fig3_v98.py')
shutil.copy2(base/'figure2_dense_20261024'/'plot_dense.py',archive/'plot_fig2_dense.py')
prefix=(archive/'compose_fig3_v98.py').read_text(encoding='utf-8').split('# a: adjacent intact structures')[0]
prefix=prefix.replace('W,H=179.,200.','W,H=183.,202.')
prefix=prefix.replace('FS=8.75','FS=9.0')
prefix=prefix.replace("'axes.linewidth':.55","'axes.linewidth':.6")
prefix=prefix.replace("txt(x,y,label,10,fontweight='bold')","txt(x,y,label,11,fontweight='bold')")
prefix=prefix.replace("txt(x+6,y+.25,title,FS,fontweight='normal')","txt(x+6,y+.25,title,10,fontweight='normal')")
prefix=prefix.replace("# All model and MD arrays are frozen; layout, typography, and annotations only.","# All model and MD arrays are frozen; layout, typography, and annotations only.\n# Style target: verified current manuscript Fig2; layout: PaperBanana B concept.")
body=(root/'native_body.py.inc').read_text(encoding='utf-8')
(out/'compose_fig3_v99.py').write_text(prefix+body,encoding='utf-8')
manifest={'native_data_source':str(old),'figure2_style_source':str(style),
 'figure2_png':str(base/'figure2_dense_20261024'/'fig2_dense.png'),
 'figure2_png_sha256':hashlib.sha256((base/'figure2_dense_20261024'/'fig2_dense.png').read_bytes()).hexdigest(),
 'layout_reference_only':str(base/'figure3_paperbanana_20261027_从头生成方案B.png'),
 'all_data_and_structure_files':[{'relative_path':str(p.relative_to(old)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for folder in [old/'figure3'/'source_data',old/'pymol'] for p in sorted(folder.rglob('*')) if p.is_file()],
 'style':{'family':'Arial','axis_tick_legend_pt':9,'title_pt':10,'panel_letter_pt':11,'axis_width_pt':.6,'Full':'#2166AC','Intra':'#B87936','MD':'#62676B','Full_linewidth_pt':1.25,'Intra_linewidth_pt':1.1},
 'existing_manuscript_not_modified':True}
(root/'input_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'prepared':str(root),'source_numbers':853,'style_reference':manifest['figure2_png']},ensure_ascii=False))
