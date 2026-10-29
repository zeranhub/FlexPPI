from pathlib import Path
import hashlib,json,shutil
root=Path(__file__).resolve().parent
old=root.parent/'revision_v99_fig2_style_B'
out=root/'figure3';out.mkdir(parents=True,exist_ok=True)
shutil.copytree(old/'figure3'/'source_data',out/'source_data',dirs_exist_ok=True)
shutil.copytree(old/'pymol',root/'pymol',dirs_exist_ok=True)
shutil.copy2(old/'figure3'/'publication_style.py',out/'publication_style.py')
prefix=(old/'figure3'/'compose_fig3_v99.py').read_text(encoding='utf-8').split('# Copy-first native v98 production logic;')[0]
prefix=prefix.replace('W,H=183.,202.','W,H=183.,222.')
prefix=prefix.replace('# Style target: verified current manuscript Fig2; layout: PaperBanana B concept.','# Style target: verified current manuscript Fig2; evidence-focused redesign with horizontal per-chi row.')
(out/'compose_fig3_v100.py').write_text(prefix+(root/'native_body.py.inc').read_text(encoding='utf-8'),encoding='utf-8')
(root/'input_manifest.json').write_text(json.dumps({'source':str(old),'no_new_model_or_MD':True,'source_data_files':[{'file':str(p.relative_to(old)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for folder in [old/'figure3'/'source_data',old/'pymol'] for p in sorted(folder.rglob('*')) if p.is_file()],'changes':['b Intra moved to supplement; absolute MD/Full plus existing static CaM shielding','d all six Gln ordered pairs and separate existing SASA','e four horizontally aligned axes','f paired correlations for the same 83 records; identity comparison'],'existing_word_unchanged':True},ensure_ascii=False,indent=2),encoding='utf-8')
print('Prepared '+str(root))
