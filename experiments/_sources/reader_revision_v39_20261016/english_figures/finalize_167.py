from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parents[2]
OLD=BASE/'experiment_completion/figure7_a_annotation_v38_20261016'
NEW=ROOT/'fig7'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
comparisons=[]
for ext in ['points.csv','display_widths.csv','all151_scatter_source.csv']:
    a=OLD/'figures'/('fig7_clear_a_annotations_and_MD_curves.'+ext)
    b=NEW/('fig7_english.'+ext)
    assert a.read_bytes()==b.read_bytes(),ext
    comparisons.append({'data':ext,'sha256':sha(b),'byte_exact_to_v38':True})
for sub in ['source','experiment']:
    for p in (OLD/sub).iterdir():
        if p.is_file():assert p.read_bytes()==(NEW/sub/p.name).read_bytes(),str(p)
audit=json.loads((NEW/'fig7_english.audit.json').read_text(encoding='utf8'))
audit['status']='PASS'
audit['visual_review']='PASS: full PNG inspected; all labels English, no overlap/cropping, a legend above curve at upper right; a secondary type-mean label removed; b has no coefficient table; c retains all three MD curves.'
audit['exact_data_comparisons']=comparisons
audit['all_source_and_experiment_files_byte_exact_to_v38']=True
(NEW/'fig7_english.audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
manifest={'status':'PASS','scope':'Figures 1, 6 and 7 English visible text',
          'figures':{'1':{'replace_docx_image':False,'reason':'Active image1 already English','png':str(ROOT/'fig1/fig1_english.png'),'pdf':str(ROOT/'fig1/fig1_english.pdf'),'svg':str(ROOT/'fig1/fig1_english.svg')},
                     '6':{'replace_docx_image':False,'reason':'Active image5 already English','png':str(ROOT/'fig6/fig6_english.png'),'pdf':str(ROOT/'fig6/fig6_english.pdf'),'svg':str(ROOT/'fig6/fig6_english.svg')},
                     '7':{'replace_docx_image':True,'png':str(NEW/'fig7_english.png'),'pdf':str(NEW/'fig7_english.pdf'),'svg':str(NEW/'fig7_english.svg'),'caption_zh':str(NEW/'fig7_english.caption_zh.md'),'all715_points_byte_exact_to_v38':True}},
          'no_model_training_or_MD':'No training, inference or trajectory changes; label-only figure edit.'}
(ROOT/'manifest_167.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(manifest,ensure_ascii=False))
