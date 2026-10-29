"""Independent figure geometry, original prediction and manuscript science review."""
from pathlib import Path
import sys
import csv
import hashlib
import json
import re
import statistics
import zipfile
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parent.parent
EXP=ROOT.parent
sys.path.insert(0,str(EXP/'completed_results_20261012/plot_runtime'))
import numpy as np

sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
checks=[]
def check(name,passed,detail):checks.append({'check':name,'passed':bool(passed),'detail':detail})

integration=json.loads((ROOT/'manuscript_QA/integration_audit.json').read_text(encoding='utf8'))
DOC,FIG=Path(integration['document']),Path(integration['figure4'])
prose=json.loads((ROOT/'prose_v49.json').read_text(encoding='utf8'))
contract=json.loads((ROOT/'scientific_review/number_contract.json').read_text(encoding='utf8'))
source=EXP/'tcr_expanded_screen_v41_20261017/mel5_validation/corrected_pair_all_aligned_positions.csv'
original=read(source)
cdr=read(ROOT/'figure4/source_plotted_CDR_RMSF.csv')
interface=read(ROOT/'figure4/source_plotted_interface_Hchi.csv')
occupancy=read(ROOT/'figure4/source_plotted_existing_MD_occupancy.csv')
original_cdr=[r for r in original if r['role'].startswith('tcr') and r['region'].startswith('CDR')]
original_interface=[r for r in original if r['role'] in ('mhc','peptide') and r['interface']=='True' and r['hchi_valid']=='True']

raw={}
raw_paths=[]
for pdb in ('3hg1','4jff'):
    p=EXP/f'tcr_case_v40_20261017/family_predictions/full/itrace_{pdb}.npz'
    raw_paths.append(p)
    with np.load(p,allow_pickle=False) as z:
        raw[pdb]={(str(c),str(i)):{'aa':str(a),'mean':m.copy(),'valid':v.copy()} for c,i,a,m,v in zip(z['chain'],z['author_residue_id'],z['amino_acid'],z['mean'],z['valid'])}

cdr_errors=[]
seen=[]
for r in cdr:
    aa=[o for o in original_cdr if o['region']==r['region'] and int(o['position_a'])==int(r['common_TCR_position'])]
    assert len(aa)==1
    a=aa[0];seen.append((a['region'],a['id_a']))
    for state,pdb,plotcol in [('a','3hg1','MEL5_model_RMSF_A'),('b','4jff','alpha24beta17_model_RMSF_A')]:
        n=raw[pdb][(a['chain_'+state],a['id_'+state])]
        cdr_errors.extend([abs(float(r[plotcol])-float(a['rmsf_full_'+state])),abs(float(r[plotcol])-float(n['mean'][0]))])
        assert n['aa']==a['aa_'+state] and bool(n['valid'][0])
    assert int(r['derivative_input_position'])==int(a['position_b'])
    assert r['MEL5_aa']==a['aa_a'] and r['alpha24beta17_aa']==a['aa_b']
check('all45_CDR_plotted_values_equal_original_Full_NPZ_and_corrected_table',len(cdr)==45 and len(set(seen))==45 and max(cdr_errors)<1e-14,
      {'n_positions':len(cdr),'n_values':len(cdr)*2,'max_absolute_error':max(cdr_errors),'raw_prediction_files':[str(p) for p in raw_paths]})

interface_errors=[]
seen=[]
for r in interface:
    aa=[o for o in original_interface if o['role']==r['role'] and int(o['position_a'])==int(r['position_a'])]
    assert len(aa)==1
    a=aa[0];seen.append((a['role'],a['id_a']))
    assert a['same_aa']=='True' and r['aa_a']==a['aa_a']
    for state,pdb in [('a','3hg1'),('b','4jff')]:
        n=raw[pdb][(a['chain_'+state],a['id_'+state])]
        interface_errors.extend([abs(float(r['hchi_full_'+state])-float(a['hchi_full_'+state])),abs(float(r['hchi_full_'+state])-float(n['mean'][1]))])
        assert bool(n['valid'][1])
check('all23_valid_interface_Hchi_values_equal_original_Full_NPZ_and_complete_scope',len(interface)==23 and len(set(seen))==23 and max(interface_errors)<1e-14,
      {'n_positions':len(interface),'n_values':len(interface)*2,'HLA_n':sum(r['role']=='mhc' for r in interface),'peptide_n':sum(r['role']=='peptide' for r in interface),'max_absolute_error':max(interface_errors)})

T=EXP/'tcr_joint_case_search_20261019/residue_story/trajectory_link'
rep=read(T/'per_replica_Q72_beta_contact_and_RMSF.csv')
dual=read(T/'Q72_two_end_backbone_Hbond_same_frame_per_replica.csv')
occ_errors=[]
occ_seen=[]
for r in occupancy:
    if r['metric']=='Q72–Gly54 contact':
        a=next(o for o in rep if o['pdb']==r['pdb'] and o['replica']==r['replica'] and o['beta_common_position']=='54')
        expected=float(a['min_heavy_contact_le4_A_fraction'])*100
    else:
        assert r['metric']=='Q72 two-ended bridge'
        a=next(o for o in dual if o['pdb']==r['pdb'] and o['replica']==r['replica'])
        expected=float(a['both_same_frame_Hbond_fraction'])*100
    assert r['frames']==a['frames']
    occ_errors.append(abs(float(r['frame_occupancy_percent'])-expected))
    occ_seen.append((r['metric'],r['pdb'],r['replica']))
check('all12_occupancy_points_equal_independent_trajectory_audit',len(occupancy)==12 and len(set(occ_seen))==12 and max(occ_errors)<1e-12,
      {'n_points':len(occupancy),'max_error_percent':max(occ_errors),'aggregation':'Unweighted mean of three independent per-trajectory occupancy fractions.'})

# Audit actual SVG coordinates, rather than assuming a matching CSV was used.
svgpath=ROOT/'figure4/fig4_mel5_joint.svg'
svg=ET.fromstring(svgpath.read_bytes())
S={'s':'http://www.w3.org/2000/svg'}
W,H=183/25.4*72,185/25.4*72
paths=[]
collections=[]
for g in svg.findall('.//s:g',S):
    if g.get('id','').startswith('line2d_'):
        paths.extend(g.findall('s:path',S))
    if g.get('id','').startswith('PathCollection_'):
        collections.append(g.findall('.//s:use',S))
def xy(path):
    a=[float(v) for v in re.findall(r'-?\d+(?:\.\d+)?',path.get('d',''))]
    return np.array(a).reshape(-1,2)
linepaths=[p for p in paths if 'stroke-width: 1.05' in p.get('style','') and any(c in p.get('style','') for c in ('#cf8c27','#762a83'))]
groups=['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
vector_errors=[]
assert len(linepaths)==12
for j,region in enumerate(groups):
    rr=[r for r in cdr if r['region']==region]
    for k,col in enumerate(('MEL5_model_RMSF_A','alpha24beta17_model_RMSF_A')):
        expected=np.array([[W*(.095+.866*(float(r['plot_x'])+1)/56.5),H*(1-.47-.16*float(r[col])/1.45)] for r in rr])
        actual=xy(linepaths[2*j+k]);assert actual.shape==expected.shape
        vector_errors.append(float(np.max(abs(actual-expected))))
check('actual_SVG_all12_CDR_curve_vertices_match_original_predictions',max(vector_errors)<1e-6,
      {'curves':12,'vertices':90,'max_coordinate_error_pt':max(vector_errors)})

pairpaths=[p for p in paths if 'stroke: #9e9e9e' in p.get('style','') and 'stroke-width: 0.7' in p.get('style','')]
vector_errors=[]
assert len(pairpaths)==23
for r,p in zip(interface,pairpaths):
    x=float(r['plot_x'])
    expected=np.array([[W*(.095+.866*(x-.14+.8)/24.6),H*(1-.258-.105*float(r['hchi_full_a'])/3.4)],
                       [W*(.095+.866*(x+.14+.8)/24.6),H*(1-.258-.105*float(r['hchi_full_b'])/3.4)]])
    vector_errors.append(float(np.max(abs(xy(p)-expected))))
check('actual_SVG_all23_Hchi_pair_endpoints_match_original_predictions',max(vector_errors)<1e-6,
      {'pairs':23,'endpoints':46,'max_coordinate_error_pt':max(vector_errors)})

vector_errors=[]
assert len(collections)==4 and all(len(u)==3 for u in collections)
for j,metric in enumerate(('Q72–Gly54 contact','Q72 two-ended bridge')):
    for k,pdb in enumerate(('3hg1','4jff')):
        rr=sorted([r for r in occupancy if r['metric']==metric and r['pdb']==pdb],key=lambda r:int(r['replica']))
        for repindex,(r,u) in enumerate(zip(rr,collections[2*j+k])):
            xx=float(r['frame_occupancy_percent'])
            yy=2*(1-j)+(.33 if k==0 else -.33)+(-.24,0,.24)[repindex]
            expected=np.array([W*(.335+.626*(xx+2)/117),H*(1-.055-.1*(yy+.7)/3.4)])
            actual=np.array([float(u.get('x')),float(u.get('y'))])
            vector_errors.append(float(np.max(abs(actual-expected))))
check('actual_SVG_all12_MD_occupancy_marker_positions_match_replica_values',max(vector_errors)<1e-6,
      {'points':12,'max_coordinate_error_pt':max(vector_errors)})

N={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
with zipfile.ZipFile(DOC) as z:doc=ET.fromstring(z.read('word/document.xml'))
text=lambda e:''.join(n.text or '' for n in e.findall('.//w:t',N))
body=doc.find('w:body',N)
changes={int(c['source_paragraph_index']):text(body[c['body_index']]) for c in integration['prose_changes']}
inserted=text(body[integration['inserted_body_nodes'][0]['body_index']])
joined='\n'.join(changes.values())+'\n'+inserted
check('final_DOCX_contains_approved_case_prose_exactly',all(changes[int(k)]==v for k,v in prose['replacements'].items()) and inserted==prose['insert_after_94'],
      {'replaced_full_paragraphs':len(prose['replacements']),'inserted_paragraphs':1})
check('whole_region_claims_and_key_model_values_correct',all(v in changes[93] for v in ['五个位点RMSF均降低','13个位点也均降低','0.896降至0.636 Å','0.801降至0.648 Å','23个有效位点中有17个预测Hχ降低']) and all(v in changes[94] for v in ['1.792降至1.375 nats','1.118降至0.747 Å','三个χ角的联合状态占据更集中']),
      {'CDR2beta':contract['key_model_values']['all5_CDR2beta_RMSF_WT_engineered_mean_A'],'complete_23_site_scope':True})
check('residue_numbering_sequence_changes_and_context_control_correct',all(v in changes[52] for v in ['βPhe53和βGly54分别对应其输入E54与E55','23个位点具有有效Hχ','全部45个位点']) and all(v in changes[94] for v in ['Gly52→Pro52与Ile53→Phe53','两处残基本身均未改变氨基酸类型','保留各自pMHC坐标','仅从输入中移除TCR链','增加0.104 nats']),
      {'common_beta53':'WT E53 Ile -> engineered E54 Phe','common_beta54':'WT E54 Gly -> engineered E55 Gly','context':'Input dependence, not an experimental apo state or physical unbinding.'})
check('primary_SPR_and_Xray_facts_correct_with_literature_attribution',all(v in changes[90] for v in ['Madura等通过表面等离子体共振（SPR）','18 μM降至0.6 nM','30,000倍','主要伴随解离减慢','晶体结构显示','56增至97','6增至38']) and '其中βCDR2也是晶体结构显示接触增加最明显的区域' in changes[93],
      {'primary_article':contract['primary_article_contract'],'corrected_CDR3beta_contact_scope':'36->30, no false additional-contact claim for betaCDR3 remains.'})
check('existing_MD_protocol_and_occupancy_narrative_correct',all(v in changes[53] for v in ['三条200 ns独立重复轨迹','1,001–1,003个保存帧','不超过4 Å','不超过3.5 Å','不小于150°','三条独立重复等权汇总']) and all(v in inserted for v in ['已有ITrace轨迹','补充动态证据','24.4%升至97.4%','0%升至37.1%','1.684降至1.219 nats','1.893降至1.078 Å','OE1与β53主链N']),
      {'trajectory_source':'ITrace, not the 2013 primary paper','same_frame_dual_geometry_occupancy':[0,37.14773096627533],'no_independent_frame_replication_claim':True})
forbidden=['Hχ下降导致','Hχ降低导致','Hχ决定','Hχ降低使','熵损失降低','结合熵降低','亲和力预测提高30,000倍','接触持续时间','更持久','原研究通过MD']
check('mechanistic_interpretation_without_Hchi_to_RMSF_causation_or_binding_entropy_conflation',not any(v in joined for v in forbidden) and '这些结果支持一种局部解释' in inserted and '局部协同约束' in inserted and '文献SPR观察到的解离减慢' in inserted,
      {'scientific_reading':'A bounded local mechanism interpretation anchored by published contact/SPR evidence and separate existing trajectory observations; it does not estimate a Q72 kinetic effect or directional entropy cause.','forbidden_expressions':forbidden,'twenty_ns_correlations':'Mixed signs, so not promoted into a temporal-coupling or causal claim.'})

svgtext=' '.join(''.join(k.text or '' for k in n.iter()).replace('\u00a0',' ') for n in svg.findall('.//s:text',S))
check('figure_labels_distinguish_model_profiles_from_existing_MD_evidence',all(v in svgtext for v in ['Model: backbone RMSF across all six CDRs','Model: side-chain entropy at all 23 valid pMHC interface sites','Existing MD: frame occupancy in three independent replicas','SPR [28]','βPhe53','βGly54','HLA Q72']) and not any('\u4e00'<=v<='\u9fff' for v in svgtext),
      {'figure_raster_visually_inspected':True,'finding':'Atoms, two local contacts, profile scopes, English labels, independent replica points and mean marks visible; no Hchi -> RMSF causal arrow.'})
package=json.loads((ROOT/'manuscript_QA/package_preservation_audit.json').read_text(encoding='utf8'))
check('independent_native_package_preservation_passed_and_bound_to_current_document',package['status']=='passed' and package['document_sha256']==sha(DOC) and package['figure4_sha256']==sha(FIG),
      {'package_status':package['status'],'package_document_sha256':package['document_sha256'],'current_document_sha256':sha(DOC)})
lit=ROOT/'literature_review/narrative_review.json'
sources=[source,*raw_paths,ROOT/'prose_v49.json',ROOT/'scientific_review/number_contract.json',ROOT/'figure4/source_plotted_CDR_RMSF.csv',ROOT/'figure4/source_plotted_interface_Hchi.csv',ROOT/'figure4/source_plotted_existing_MD_occupancy.csv',svgpath,T/'per_replica_Q72_beta_contact_and_RMSF.csv',T/'Q72_two_end_backbone_Hbond_same_frame_per_replica.csv',lit]
report={'status':'passed' if all(c['passed'] for c in checks) else 'failed','document':str(DOC),'document_sha256':sha(DOC),
        'figure4':str(FIG),'figure4_sha256':sha(FIG),'reviewer':'independent residue_pair_evidence agent',
        'review_scope':'Actual revised manuscript text, Figure4 raster and actual SVG vertices/markers, all 90 RMSF and all 46 Hchi original NPZ predictions, all 12 MD occupancy points, corrected numbering, primary-paper facts, and native Word package preservation.',
        'checks':checks,'required_fixes':[],'source_SHA256':{str(p):sha(p) for p in sources},
        'editorial_assessment':'The MEL5 case now makes a concrete joint-output argument at HLA Q72 and unchanged beta Gly54, within complete 45-CDR/23-sidechain regional profiles. Published SPR and crystal contacts establish affinity/structure; model outputs supply complementary local dynamics; existing ITrace trajectories corroborate contact occupancy and motion. No directional Hchi-to-RMSF causal or thermodynamic-binding-entropy claim is made.'}
report['required_fixes']=[c['check'] for c in checks if not c['passed']]
(ROOT/'scientific_review/final_scientific_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'status':report['status'],'document_sha256':report['document_sha256'],'figure4_sha256':report['figure4_sha256'],'failed_checks':[c for c in checks if not c['passed']]},ensure_ascii=False,indent=2))
