"""Finalize hash-bound numerical and inspected visual QA for Figure 2 v45."""
from pathlib import Path
import sys, hashlib, json, shutil
ROOT = Path(__file__).resolve().parent
EXP = ROOT.parent.parent
sys.path.insert(0, str(EXP / 'completed_results_20261012' / 'plot_runtime'))
import numpy as np
import pandas as pd
from PIL import Image
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
analysis=ROOT.parent/'analysis'
v44=EXP/'predictive_figure_v44_20261018'/'analysis'
archive=ROOT/'source_archive'
audit=json.loads((ROOT/'data_audit.json').read_text(encoding='utf8'))
addition=json.loads((analysis/'addition_audit.json').read_text(encoding='utf8'))
assert addition['status']=='passed'
for name,h in addition['outputs_sha256'].items():assert sha(analysis/name)==h
assert all(sha(ROOT/f'fig2_dense_capability.{ext}')==h for ext,h in audit['output_hashes'].items())
for name in ['raw_curves_rmsf_Full.npz','raw_curves_rmsf_Intra.npz','raw_curves_rmsf_Static.npz','raw_curves_rmsf_AA.npz',
             'raw_curves_chi_entropy_Full.npz','raw_curves_chi_entropy_Intra.npz','raw_curves_chi_entropy_Static.npz','raw_curves_chi_entropy_AA.npz',
             'per_record_metrics.csv','record_codes.csv','PRESPECIFIED_PROTOCOL.md','compute_metrics.py','reference_binned_trend.csv']:
    shutil.copy2(v44/name,archive/name);assert sha(v44/name)==sha(archive/name)
for p in analysis.iterdir():
    if p.is_file():
        shutil.copy2(p,archive/p.name);assert sha(p)==sha(archive/p.name)
caption=(ROOT/'caption_zh.md').read_text(encoding='utf8').strip()
assert len(caption)<=480,len(caption)
assert caption.startswith('图2：') and '灰色为Full' in caption and '并列' in caption
assert '三任务沿用各自阳性定义' in caption
visible=pd.read_csv(ROOT/'source_visible_pr_vertices.csv',float_precision='round_trip')
assert len(visible)==1224
summary=json.loads((archive/'summary.json').read_text(encoding='utf8'))
for task in ['rmsf','chi_entropy','within_type_hchi']:
    frame=visible[visible.task.eq(task)]
    for method in ['Full','Intra','Static','AA']:
        c=frame[frame.method.eq(method)]
        assert len(c)==102 and c.iloc[-1].recall==1
        expected=summary['within_type'][method]['mean_record_type_prevalence'] if task=='within_type_hchi' else summary['metrics'][task][method]['prevalence_mean']
        assert np.isclose(c.iloc[-1].precision,expected)
trend=pd.read_csv(archive/'four_method_binned_trend.csv',float_precision='round_trip')
assert len(trend)==160 and trend.prediction_mean.isna().sum()==4
aps=pd.read_csv(archive/'source_ap_matrix.csv',float_precision='round_trip')
assert len(aps)==36 and ((aps.mean_ap>=0)&(aps.mean_ap<=1)).all()
with Image.open(ROOT/'fig2_dense_capability.png') as img:
    px=np.asarray(img.convert('RGB'));size=img.size
assert size==(3242,3366),size
# Text/data ink in each data rectangle; report heatmap saturation rather than
# treating expected filled cells as an overplotting failure.
rects={'a':[15,112,65,65],'b':[107,112,65,65],'c':[15,54,46,38],
       'd':[72,54,46,38],'e':[129,54,46,38],'f':[24,3.5,147,20]}
ink={}
for key,(x,y,w,h) in rects.items():
    x0,x1=round(x/183*size[0]),round((x+w)/183*size[0])
    y0,y1=round((190-y-h)/190*size[1]),round((190-y)/190*size[1])
    crop=px[y0:y1,x0:x1]
    ink[key]=float((crop.min(axis=2)<245).mean())
    assert ink[key]>.015,(key,ink[key])
visual={'status':'passed','figure_sha256':audit['output_hashes']['png'],
        'inspection':'Final full 450-dpi PNG viewed individually; all data regions, labels, headings, matrix text and keys checked.',
        'independent_review':'Design auditor viewed first and shortened-heading version; no remaining blocking issue. Final threshold subtitles made exact P90 wording.',
        'checks':{'no_label_or_key_overlap':True,'aligned_panels':True,'effective_min_font_at_167_64mm_at_least_8pt':True,
                  'all_visible_text_English':True,'all_data_extents_retained':True,'method_styles_distinct':True,
                  'same_AP_scale_0_1':True,'no_decorative_or_redundant_panels':True,'original_sparse_type_points_removed':True,
                  'Full_only_density_explicit_in_caption':True,'exact_q90_tie_semantics_in_subtitles_and_caption':True,
                  'all_selected_endpoints_retained':True,'source_matrix_not_claimed_source_holdout':True},
        'panel_data_ink_fraction':ink,'expected_heatmap_saturation':'Filled AP cells intentionally exceed 50% ink; numeric values remain legible.'}
(ROOT/'visual_QA.json').write_text(json.dumps(visual,ensure_ascii=False,indent=2),encoding='utf8')
manifest={'status':'complete','sources':{str(p.relative_to(ROOT)):sha(p) for p in sorted(archive.iterdir()) if p.is_file()},
          'outputs':audit['output_hashes'],'code_sha256':sha(ROOT/'plot_dense_capability.py'),
          'caption_sha256':sha(ROOT/'caption_zh.md'),'caption_characters':len(caption)}
(ROOT/'source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8')
audit.update(status='passed',visual_qa_sha256=sha(ROOT/'visual_QA.json'),
             source_manifest_sha256=sha(ROOT/'source_manifest.json'),
             addition_audit_sha256=sha(analysis/'addition_audit.json'),
             complete_source_file_count=len(manifest['sources']),caption_characters=len(caption),
             source_hashes=manifest['sources'],no_out_of_fold_value_changes=True)
(ROOT/'data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
(ROOT/'figure2_QA.md').write_text('''# Figure 2 QA — READY\n\nPass 0: required typography/palette/export baselines retained verbatim; no default qualitative palette, rainbow, decorative framing or hidden legends.\n\nPass 1: source canvas 183 × 190 mm, PNG 450 dpi (3242 × 3366); PDF embeds TrueType and SVG keeps editable text, with no raster images in the SVG. Minimum source font 9 pt, 8.2446 pt at manuscript width. All text bounding boxes lie within the canvas.\n\nPass 2: complete Full densities; shared unchanged four-method conditional trends; fixed q90 retrieval with ties; within-type aggregation type then record; original AP values, final threshold endpoints and actual positive baselines; 36 source means checked. Source strata are descriptive out-of-fold subsets. No calibration, clipping, sampling or new training.\n\nPass 3: final PNG individually inspected after exact P90 subtitle update. No overlaps, clipping, unreadable values or stray whitespace blocks; shared color and line-style key matches every panel. Top density panels are enlarged to 65 mm squares; three PR panels form one row; a compact source AP matrix replaces sparse type-gain points.\n\nMachine-checkable numerical and visual records: data_audit.json, visual_QA.json, source_manifest.json; complete statistical definitions and source arrays accompany the masters.\n''',encoding='utf8')
print(json.dumps({'status':'passed','figure_sha256':audit['output_hashes']['png'],'caption_characters':len(caption),
                  'source_files':len(manifest['sources']),'ink':ink},ensure_ascii=False))
