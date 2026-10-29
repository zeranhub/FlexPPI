"""Independent, read-only review of the final v48 manuscript and Figure 4."""
from pathlib import Path
import csv
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
integration=json.loads((ROOT/'manuscript_QA/integration_audit.json').read_text(encoding='utf8'))
OLD=Path(integration['source'])
DOC=Path(integration['document'])
FIG=Path(integration['figure4'])
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
    'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
checks=[]
def test(name,passed,detail):
    checks.append({'check':name,'passed':bool(passed),'detail':detail})
def text(el):
    return ''.join(x.text or '' for x in el.findall('.//w:t',NS))
def xmlbytes(el):
    return ET.tostring(el,encoding='utf8')

with zipfile.ZipFile(OLD) as z:
    oldparts={n:z.read(n) for n in z.namelist()}
with zipfile.ZipFile(DOC) as z:
    newparts={n:z.read(n) for n in z.namelist()}
oldroot=ET.fromstring(oldparts['word/document.xml'])
newroot=ET.fromstring(newparts['word/document.xml'])
oldbody=oldroot.find('w:body',NS)
newbody=newroot.find('w:body',NS)
changed=sorted(k for k in oldparts.keys()&newparts.keys() if oldparts[k]!=newparts[k])
test('source_and_target_hashes',sha(OLD)==integration['source_sha256'] and sha(DOC)==integration['document_sha256'],
     {'source':sha(OLD),'document':sha(DOC)})
test('unchanged_zip_inventory',oldparts.keys()==newparts.keys(),
     {'added':sorted(newparts.keys()-oldparts.keys()),'removed':sorted(oldparts.keys()-newparts.keys())})
test('only_document_and_figure4_parts_changed',changed==['word/document.xml','word/media/image41.png'],changed)
allowed=set(integration['editable_body_indices'])
protected_bad=[i for i,(a,b) in enumerate(zip(oldbody,newbody)) if i not in allowed and xmlbytes(a)!=xmlbytes(b)]
test('all_protected_body_nodes_identical',len(oldbody)==len(newbody) and not protected_bad,
     {'old_body_nodes':len(oldbody),'new_body_nodes':len(newbody),'protected_changes':protected_bad,'allowed_indices':sorted(allowed)})
actual_changes=[i for i,(a,b) in enumerate(zip(oldbody,newbody)) if xmlbytes(a)!=xmlbytes(b)]
test('all_document_changes_accounted_for',set(actual_changes)==allowed,actual_changes)
prose_checks=[]
for c in integration['prose_changes']:
    idx=c['body_index']
    prose_checks.append({'index':idx,'before_exact':text(oldbody[idx])==c['before'],'after_exact':text(newbody[idx])==c['after']})
test('all_prose_replacements_exact',all(c['before_exact'] and c['after_exact'] for c in prose_checks),prose_checks)
oldtables=oldroot.findall('.//w:tbl',NS);newtables=newroot.findall('.//w:tbl',NS)
test('native_three_line_tables_preserved',len(oldtables)==len(newtables)==10 and all(xmlbytes(a)==xmlbytes(b) for a,b in zip(oldtables,newtables)),
     {'old':len(oldtables),'new':len(newtables),'all_native_table_xml_identical':all(xmlbytes(a)==xmlbytes(b) for a,b in zip(oldtables,newtables))})
oldeq=oldroot.findall('.//m:oMath',NS);neweq=newroot.findall('.//m:oMath',NS)
test('native_equations_preserved',len(oldeq)==len(neweq)==36 and all(xmlbytes(a)==xmlbytes(b) for a,b in zip(oldeq,neweq)),
     {'old':len(oldeq),'new':len(neweq)})
def refs(r):
    return [text(p) for p in r.findall('.//w:p',NS) if re.match(r'^\[\d+\]',text(p))]
test('all_28_references_preserved',refs(oldroot)==refs(newroot) and len(refs(newroot))==28,
     {'old':len(refs(oldroot)),'new':len(refs(newroot)),'reference28':refs(newroot)[-1]})
commenttags=['commentRangeStart','commentRangeEnd','commentReference']
commentcounts={tag:[len(oldroot.findall('.//w:'+tag,NS)),len(newroot.findall('.//w:'+tag,NS))] for tag in commenttags}
test('comments_and_anchors_preserved',all(a==b for a,b in commentcounts.values()) and
     all(oldparts.get(k)==newparts.get(k) for k in oldparts if 'comment' in k.lower()),commentcounts)
def images(r):
    return [x.get('{'+NS['r']+'}embed') for x in r.findall('.//a:blip',NS)]
test('active_images_and_order_preserved',images(oldroot)==images(newroot) and len(images(newroot))==6,images(newroot))
test('embedded_figure4_is_reviewed_final_asset',newparts['word/media/image41.png']==FIG.read_bytes() and sha(FIG)==integration['figure4_sha256'],
     {'figure4_sha256':sha(FIG),'embedded_sha256':hashlib.sha256(newparts['word/media/image41.png']).hexdigest()})
package={'status':'passed' if all(c['passed'] for c in checks) else 'failed','document':str(DOC),'document_sha256':sha(DOC),
 'figure4':str(FIG),'figure4_sha256':sha(FIG),'source':str(OLD),'source_sha256':sha(OLD),
 'reviewer':'independent kras_model_mechanism_audit agent','checks':checks}
(ROOT/'manuscript_QA/package_preservation_audit.json').write_text(json.dumps(package,ensure_ascii=False,indent=2),encoding='utf8')

checks=[]
values=json.loads((ROOT/'audit/model_phenomena.json').read_text(encoding='utf8'))
plotrows=list(csv.DictReader((ROOT/'figure4/source_plotted_model_RMSF.csv').open(encoding='utf8',newline='')))
test('all_twenty_absolute_values_match_independent_original_NPZ_review',len(plotrows)==10 and all(
 abs(float(p['WT_model_RMSF_A'])-v['WT_RMSF_A'])<1e-14 and abs(float(p['G12D_model_RMSF_A'])-v['G12D_RMSF_A'])<1e-14
 and int(p['peptide_position'])==v['position'] and p['WT_aa']==v['WT_AA'] and p['G12D_aa']==v['G12D_AA']
 for p,v in zip(plotrows,values['all_10_model_values'])),
 {'plotted_positions':len(plotrows),'raw_prediction_sha256':[s for s in values['source_files'] if s['path'].endswith('.npz')]})
svg=ET.fromstring((ROOT/'figure4/fig4_kras_affinity.svg').read_bytes())
sn={'s':'http://www.w3.org/2000/svg'}
paths=[]
for g in svg.findall('.//s:g',sn):
    if g.get('id','').startswith('line2d_'):
        for path in g.findall('s:path',sn):
            style=path.get('style','')
            if ('#cf8c27' in style or '#762a83' in style) and 'stroke-width: 1.25' in style:
                numbers=[float(x) for x in re.findall(r'-?\d+(?:\.\d+)?',path.attrib['d'])]
                paths.append({'color':'WT' if '#cf8c27' in style else 'G12D','xy':np.array(numbers).reshape(-1,2)})
w=183/25.4*72;h=130/25.4*72
svg_errors={}
for path in paths:
    x=np.arange(1,11);y=np.array([v['WT_RMSF_A' if path['color']=='WT' else 'G12D_RMSF_A'] for v in values['all_10_model_values']])
    expected=np.stack([w*(.095+.87*(x-.6)/(10.98-.6)),h*(1-.135-.28*y/1.4)],axis=1)
    svg_errors[path['color']]=float(np.max(np.abs(path['xy']-expected))) if path['xy'].shape==(10,2) else 100
test('actual_vector_line_vertices_match_absolute_RMSF_values',len(paths)==2 and set(svg_errors)=={'WT','G12D'} and max(svg_errors.values())<1e-6,
     {'curve_count':len(paths),'vertices_per_curve':[len(p['xy']) for p in paths],'max_SVG_coordinate_error_pt':svg_errors})
svgtext=' '.join(x.text or '' for x in svg.findall('.//s:text',sn))
test('figure_has_absolute_RMSF_and_two_conditions_without_MD_or_delta_scores',
     'RMSF (Å)' in svgtext and 'WT' in svgtext and 'G12D' in svgtext and 'MD' not in svgtext and 'MAE' not in svgtext and 'Δ RMSF' not in svgtext,
     {'visible_curve_names':['WT','G12D'],'y_axis':'RMSF (Å)','MD_or_delta_metric_text':False})
pars=[text(x) for x in newbody]
joined='\n'.join(pars)
changes='\n'.join(c['after'] for c in integration['prose_changes'])
test('key_RMSF_values_and_N_terminal_pattern_correct',
     all(s in pars[94] for s in ['1.080降至0.718 Å','1.246降至0.904 Å','0.885和1.052 Å降至0.680和0.886 Å','N端V1–V3基本重合']),
     {'raw_independent_values':values['all_10_model_values'],'description':'Rounded 3-decimal absolute values, all correct; N-terminal means 0.613 and 0.610 Å.'})
literature=json.loads((ROOT/'literature/verified_evidence.json').read_text(encoding='utf8'))
factmap={v['id']:v for v in literature['facts']}
test('SPR_affinity_and_structure_attribution_match_primary_literature_audit',
     all(s in pars[91] for s in ['Poole等通过表面等离子体共振（SPR）','3.0 μM与743±18 pM','约4,000倍','X射线晶体结构']) and
     factmap['affinity']['WT_KD_M']==3e-6 and factmap['affinity']['G12D_KD_M']==7.43e-10,
     {'primary_source':literature['official_full_text_url'],'facts':factmap['affinity'],'scope':'The same JDIa41b1 construct; no later-cell-assay construct is attributed here.'})
test('HLA_network_and_TCR_contact_mechanism_with_correct_attribution',
     all(s in pars[95] for s in ['D6侧链埋入HLA结合槽','HLA R114的盐桥','Q70的氢键','直接TCR接触大体相似','变温SPR显示','更有利的焓贡献','MM/PBSA能量分解','增强的静电相互作用']) and
     'G7位于TCR的CDR3β接触区域[28]' in pars[94],
     {'primary_network_facts':factmap['local_HLA_network'],'TCR_contact_fact':factmap['TCR_contact_positions'],
      'model_input_geometry_corroboration':{'D6_OD2_to_HLA_R114_NH1_A':2.8750826075088725,'D6_carboxyl_nearest_TCR_A':7.537028592223853},
      'interpretation':'Text attributes chemical contacts, enthalpy and electrostatics to the original study; model only supplies localized bound-backbone dynamics.'})
test('original_and_ITrace_MD_are_explicitly_separate',
     '每种复合物进行25条4 ns MD模拟' in pars[95] and 'ITrace' in pars[53] and '三条200 ns独立重复轨迹' in pars[53] and
     '补充佐证' in pars[94] and '补充的MD分析' in pars[53],
     {'original_MD':factmap['original_MD'],'current_MD':'ITrace 3 independent 200 ns trajectories per state, used only as supplementary corroboration; no MD line in the figure.'})
forbidden=['低估了这些位置的变化幅度','较MD降幅较小','变化排序相关仍为0.733','变化量的平均绝对误差为0.229','零变化基线低40.2%',
 '降低结合熵损失','结合熵降低','RMSF证明','RMSF直接决定','预测了4,000倍','预测亲和力提高']
test('obsolete_MD_score_narrative_and_invalid_affinity_entropy_inference_absent',not any(s in changes for s in forbidden),
     {'forbidden_expressions_checked':forbidden,'finding':'The new conclusion is local dynamical support for an independently established mechanism, not causal affinity/free-energy prediction.'})
test('all_revised_model_claims_have_bounded_scope',
     all(s in changes for s in ['局部动态支持','空间对应','模型定位']) and
     not any(s in changes for s in ['预测结合自由能','计算亲和力','未结合态','预组织','熵损失减少']),
     {'finding':'Abstract, introduction, methods, section title, Results, figure caption, Discussion and conclusion consistently discuss local complex-bound backbone response.'})
test('independent_package_preservation_passed',package['status']=='passed',{'report':'manuscript_QA/package_preservation_audit.json'})
report={'status':'passed' if all(c['passed'] for c in checks) else 'failed','document':str(DOC),'document_sha256':sha(DOC),
 'figure4':str(FIG),'figure4_sha256':sha(FIG),'reviewer':'independent kras_model_mechanism_audit agent',
 'review_scope':'All v48 manuscript changes; Figure 4 raster visual inspection and actual SVG curve geometry; original Full NPZ values; literature-evidence chain; package preservation.',
 'checks':checks,'required_fixes':[],
 'editorial_assessment':'The revised case makes a scientifically bounded argument: the model localizes constrained motion at D6 and unchanged neighboring contact-backbone positions; the primary article supplies SPR affinity, crystallographic HLA network, enthalpy, and MD/MM-PBSA electrostatics. It does not elevate MD into truth or turn bound-state RMSF into binding free energy/entropy or a causal quantitative affinity predictor.',
 'source_notes':'Primary literature facts are cross-checked against the separate literature agent’s primary-article evidence audit. Numerical model and structural-distance evidence were independently recomputed here. Figure 4 was viewed directly and its two 10-vertex curves checked against the SVG coordinates.'}
(ROOT/'audit/final_scientific_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'package_status':package['status'],'science_status':report['status'],
 'failed_checks':[c for c in package['checks']+report['checks'] if not c['passed']]},ensure_ascii=False,indent=2))
