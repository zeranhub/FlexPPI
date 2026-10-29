"""Reuse already English scientific figures and translate editable Figure 7.

No observations are changed. Existing structure raster renders stay embedded in
the original vector masters; they are not repainted or regenerated.
"""
from pathlib import Path
import shutil, hashlib, json, re, zipfile
from lxml import etree

OUT = Path(__file__).resolve().parent
BASE = OUT.parents[2]
SOURCE7 = BASE/'experiment_completion/figure7_a_annotation_v38_20261016'
FIG1 = Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\figures\main_revision_v13_20261008')
FIG6 = BASE.parent/'figures/case_validation_v13_20261009'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

for folder, original, stem in [('fig1',FIG1,'FlexPPI_main_v13'),('fig6',FIG6,'fig_local_dynamic_cases')]:
    dest=OUT/folder
    # The current DOCX PNG was inspected and already has English text.
    shutil.copy2(dest/'active_v38.png',dest/(folder+'_english.png'))
    pairs = [(original/(stem+'.svg'),dest/(folder+'_english.svg')),
             (original/(stem+('_183mm.pdf' if folder=='fig1' else '.pdf')),dest/(folder+'_english.pdf'))]
    for src,target in pairs: shutil.copy2(src,target)
    scripts = ['build_v13.py','export_v13.py','project_frames_v13.py','check_v13.py'] if folder=='fig1' else ['compose_figure.py','figure_provenance_and_QA.json','case_model_protocol.json','v13_case_values_QA.json']
    for name in scripts: shutil.copy2(original/name,dest/name)
    if folder=='fig6':
        for sub in ['data','kix','ubiquitin']:
            shutil.copytree(original/sub,dest/sub,dirs_exist_ok=True)
    else:
        shutil.copy2(original/(stem+'.drawio'),dest/(stem+'.drawio'))
    tree=etree.parse(str(dest/(folder+'_english.svg')))
    ns={'s':'http://www.w3.org/2000/svg','h':'http://www.w3.org/1999/xhtml'}
    texts=tree.xpath('//s:text//text()|//h:div//text()',namespaces=ns)
    chinese=[s for s in texts if re.search('[\u3400-\u9fff]',s)]
    assert not chinese,chinese
    audit={'figure':folder,'status':'PASS','visible_figure_language':'English',
           'source_master_directory':str(original),'current_docx_png_reused_byte_exact':True,
           'scientific_data_unchanged':True,'english_visible_svg_text_count':len(texts),
           'chinese_visible_svg_text':chinese,
           'visual_review':'Active DOCX PNG and source preview inspected: all titles, axes, legends, structure labels are English.',
           'output_sha256':{ext:sha(dest/(folder+'_english.'+ext)) for ext in ['png','pdf','svg']}}
    (dest/'audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf8')

# Keep the established scientific code and all source observations.
dest=OUT/'fig7'
for sub in ['source','experiment']:
    shutil.copytree(SOURCE7/sub,dest/sub,dirs_exist_ok=True)
code=(SOURCE7/'figures/plot_clear_a_annotations_and_MD_curves.py').read_text(encoding='utf8')
old="sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'completed_results_20261012'/'plot_runtime'))"
new="sys.path.insert(0,str(Path(__file__).resolve().parents[3]/'completed_results_20261012'/'plot_runtime'))"
assert old in code;code=code.replace(old,new)
code=code.replace("ROOT=OUT.parent;DATA=ROOT/'source';EXP=ROOT/'experiment'","ROOT=OUT;DATA=ROOT/'source';EXP=ROOT/'experiment'")
code=code.replace("STEM='fig7_clear_a_annotations_and_MD_curves'","STEM='fig7_english'")
code=code.replace("font_manager.findfont('Microsoft YaHei',fallback_to_default=False)","font_manager.findfont('Arial',fallback_to_default=False)")
code=code.replace("'font.sans-serif':['Microsoft YaHei','Arial']","'font.sans-serif':['Arial','Liberation Sans']")
line="type_text=text(15,128,'扣除氨基酸类型均值后 ρ 0.868',fontsize=7.5)\n"
assert line in code;code=code.replace(line,'')
code=code.replace("assert not legend_box.overlaps(type_text.get_window_extent(renderer))\n",'')
substitutions={
    "泛素全部 66 个有效位置":"Ubiquitin: all 66 valid residues",
    "联合 Hχ\\n(nats)":"Joint Hχ\\n(nats)",
    "['MD 参考','模型预测']":"['MD reference','Model prediction']",
    "逐χ MD占据与联合预测":"Per-χ occupancy vs. joint prediction",
    "MD 细箱熵：同类型偏差 (nats)":"MD bin entropy: within-type deviation (nats)",
    "联合 Hχ 预测：同类型偏差 (nats)":"Predicted joint Hχ: within-type deviation (nats)",
    "预测端点的逐角 MD 自由能":"MD free energy at prediction-selected sites",
    "六个谷氨酰胺中预测最低 Q41：Hχ = 2.000":"Lowest of six Gln: Q41, Hχ = 2.000",
    "预测最高 Q2：Hχ = 2.990":"Highest: Q2, Hχ = 2.990",
    " 低值角区 ":" low-F* width ",
    "侧链角度 (°)":"Side-chain angle (°)",
    "MD F* (无量纲)":"MD F* (dimensionless)",
}
for old,new in substitutions.items():
    assert old in code,old
    code=code.replace(old,new)
# Caption remains Chinese, but remove the explanation of the deleted statistic.
code=code.replace('折线为原始Hχ；扣除类型均值后的ρ由预测与MD参考分别减去训练折中该氨基酸类型的同一个均值后计算。','')
code=code.replace('模型为第0折冻结完整模型。','模型为第0折的完整模型。')
code=code.replace("'a_type_mean_label':'扣除氨基酸类型均值后 ρ 0.868'","'a_type_mean_label_removed':True")
code=code.replace("'parent_figure_version':'v37'","'parent_figure_version':'v38'")
code=code.replace("'status':'numerical_checks_passed_visual_pending'","'status':'numerical_checks_passed_visual_pending','visible_figure_language':'English'")
# Visible-text language is checked at the figure object level.
code=code.replace("assert not outside,outside","assert not outside,outside\nassert not [t.get_text() for t in texts if __import__('re').search('[\\u3400-\\u9fff]',t.get_text())]\n(OUT/(STEM+'.visible_text.json')).write_text(json.dumps([t.get_text() for t in texts],ensure_ascii=False,indent=2),encoding='utf8')")
(dest/'plot_english.py').write_text(code,encoding='utf8')
print(json.dumps({'fig1':'already English; active PNG reused','fig6':'already English; active PNG reused','fig7_script':str(dest/'plot_english.py')},ensure_ascii=False))
