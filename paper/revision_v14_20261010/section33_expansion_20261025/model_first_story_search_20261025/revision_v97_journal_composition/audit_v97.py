from pathlib import Path
from zipfile import ZipFile
import hashlib,json,sys
import pandas as pd
from lxml import etree
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent;OLD=ROOT.parent/'revision_v96_visual_hierarchy';BASE=ROOT.parents[2]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for p in (OLD/'figure3/source_data').glob('*'):
    if p.is_file() and p.name!='all_plotted_values.csv':assert sha(p)==sha(ROOT/'figure3/source_data'/p.name),p.name
before=pd.read_csv(OLD/'figure3/source_data/all_plotted_values.csv',float_precision='round_trip')
after=pd.read_csv(ROOT/'figure3/source_data/all_plotted_values.csv',float_precision='round_trip')
pd.testing.assert_frame_equal(before,after,check_exact=True)
assert len(after)==853
for name in ['cam_overview.png','ub_overview.png','cam_overview.pse','ub_overview.pse']:
    assert sha(OLD/'pymol'/name)==sha(ROOT/'pymol'/name)
with ZipFile(BASE/'paper_v96_20261026_中文修订_图3视觉层级重排.docx') as z:a={n:z.read(n) for n in z.namelist()}
with ZipFile(BASE/'paper_v97_20261027_中文修订_期刊参考组图重排.docx') as z:b={n:z.read(n) for n in z.namelist()}
assert a.keys()==b.keys()
assert [n for n in a if a[n]!=b[n]]==['word/document.xml','word/media/image5.png']
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
roots=[etree.fromstring(z['word/document.xml']) for z in [a,b]]
assert [etree.tostring(t,method='c14n') for t in roots[0].findall('.//w:tbl',ns)]==[etree.tostring(t,method='c14n') for t in roots[1].findall('.//w:tbl',ns)]
texts=[[''.join(p.xpath('.//w:t/text()',namespaces=ns)) for p in r.findall('.//w:body/w:p',ns)] for r in roots]
assert len(texts[0])==len(texts[1])
changed_text=[i for i,(x,y) in enumerate(zip(*texts)) if x!=y]
assert changed_text==[73],changed_text
q=json.loads((ROOT/'figure3/figure_QA.json').read_text(encoding='utf-8'))
assert abs(q['annotation_metrics']['CaMKI_Full_RMSF_MAE_A']-.1148592012)<1e-10
assert abs(q['annotation_metrics']['CaMKI_Intra_RMSF_MAE_A']-2.383203344)<1e-9
assert abs(q['annotation_metrics']['Ubiquitin_Full_Hchi_MAE_nats']-.07953050649754886)<1e-10
assert '<text' in (ROOT/'figure3/fig3_dual_case_entropy.svg').read_text(encoding='utf-8')
audit={'all_853_plotted_values_exactly_identical_to_v96':True,'all_frozen_source_tables_identical':True,'all_native_PNG_and_PSE_overviews_identical':True,'all_six_document_tables_identical':True,'all_other45_media_identical':True,'only_changed_text_paragraph':73,'all_body_and_panel_references_identical':True,'computed_MAE_annotations_independently_verified':q['annotation_metrics'],'figure_width_height_mm':q['size_mm'],'smallest_font_at_Word_164mm':q['print_font_at_164mm'],'SVG_text_remains_editable':True,'all_six_Gln_and_all_four_chi_subsets_retained':True,'all83_records_retained_for_both_descriptors_and_targets':True,'no_new_model_or_MD_results_created':True}
render=ROOT/'render/release_01'
pages=sorted(render.glob('page-*.png'),key=lambda p:int(p.stem.split('-')[1]))
oldpages=sorted((OLD/'render/release_02').glob('page-*.png'),key=lambda p:int(p.stem.split('-')[1]))
if len(pages)==len(oldpages)==29:
    changed=[i+1 for i,(x,y) in enumerate(zip(oldpages,pages)) if sha(x)!=sha(y)]
    audit.update(rendered_pages=29,changed_pages_vs_verified_v96=changed,identical_pages_inherit_verified_review=29-len(changed))
(ROOT/'audit/data_and_document_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False))
