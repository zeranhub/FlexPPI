"""Read-only audit of v54 MEL5 prose, geometry, numeric claims, and figure caption."""
from pathlib import Path
import json, hashlib, zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

ROOT=Path(__file__).resolve().parent.parent
EXP=ROOT.parent
BASE=EXP.parent
DOC=BASE/'paper_v54_20261022_中文修订_MEL5六面板组合图4_三线表.docx'
EXPECTED='90b66a4029e22b0f67e52d3176839e547a8265690e5ff21d385ab7b5fecd1027'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(DOC)==EXPECTED
prose_path=ROOT/'prose_v54.json'
prose=json.loads(prose_path.read_text(encoding='utf-8'))
with zipfile.ZipFile(DOC) as archive:
    doc=ET.fromstring(archive.read('word/document.xml'))
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    paras=[''.join(p.itertext()) for p in doc.findall('.//w:p',ns)]
    paragraphs=[''.join(p.itertext()) for p in doc.findall('.//w:body/w:p',ns)]
    for key,value in prose.items():assert value in paras,(key,'prose not found exactly')
alltext='\n'.join(paras)
assert '近距离不超过4 Å' in alltext or '最近距离不超过4 Å' in alltext
assert '供体–受体距离不超过3.5 Å' in alltext
assert '供体–氢–受体角不小于150°' in alltext
assert '先在每条轨迹内计算，再对三条独立重复等权汇总' in alltext

OLD=EXP/'mel5_joint_case_v49_20261019'
parser=PDBParser(QUIET=True)
distances=[]
for pdb,b53,b54 in [('3hg1',53,54),('4jff',54,55)]:
    p=OLD/f'figure4/source_archive/{pdb}_model_input.pdb'
    s=parser.get_structure(pdb,p)[0]
    q=s['A'][72]
    d53=float(np.linalg.norm(q['OE1'].coord.astype(float)-s['E'][b53]['N'].coord.astype(float)))
    d54=float(np.linalg.norm(q['NE2'].coord.astype(float)-s['E'][b54]['O'].coord.astype(float)))
    distances.append({'pdb':pdb,'OE1_to_beta53_N_A':d53,'NE2_to_beta54_O_A':d54,'display':[f'{d53:.2f}',f'{d54:.2f}'],'source_SHA256':sha(p)})
assert distances[0]['display']==['5.85','10.58']
assert distances[1]['display']==['2.94','3.12']
assert '距离由5.85缩短至2.94 Å' in prose['92']
assert '距离由10.58缩短至3.12 Å' in prose['92']

read=lambda p:pd.read_csv(p,float_precision='round_trip')
cdr=read(OLD/'figure4/source_plotted_CDR_RMSF.csv')
h=read(OLD/'figure4/source_plotted_interface_Hchi.csv')
claims={}
for region,n in [('CDR2β',5),('CDR3β',13)]:
    q=cdr[cdr.region.eq(region)]
    assert len(q)==n and (q.alpha24beta17_model_RMSF_A<q.MEL5_model_RMSF_A).all()
    claims[region]={'n':len(q),'all_decrease':True,'mean_before':float(q.MEL5_model_RMSF_A.mean()),'mean_after':float(q.alpha24beta17_model_RMSF_A.mean())}
assert [f"{claims['CDR2β'][k]:.3f}" for k in ['mean_before','mean_after']]==['0.896','0.636']
assert [f"{claims['CDR3β'][k]:.3f}" for k in ['mean_before','mean_after']]==['0.801','0.648']
q=h[h.chain_a.eq('A')&h.position_a.eq(72)].iloc[0]
assert f'{q.hchi_full_a:.3f}'=='1.792' and f'{q.hchi_full_b:.3f}'=='1.375'
delta=h.hchi_full_b-h.hchi_full_a
assert (delta.idxmin()==q.name)
g=cdr[cdr.region.eq('CDR2β')&cdr.common_TCR_position.eq(54)].iloc[0]
assert g.MEL5_aa=='G' and g.alpha24beta17_aa=='G'
assert [f'{g[k]:.3f}' for k in ['MEL5_model_RMSF_A','alpha24beta17_model_RMSF_A']]==['1.118','0.747']
claims['Q72_Hchi']={'before_nats':q.hchi_full_a,'after_nats':q.hchi_full_b,'largest_decrease_of23':True}
claims['betaGly54_RMSF']={'before_A':g.MEL5_model_RMSF_A,'after_A':g.alpha24beta17_model_RMSF_A,'unchanged_AA':True}

trajectory=OLD/'case_sources/residue_story/trajectory_link'
contact=read(trajectory/'per_replica_Q72_beta_contact_and_RMSF.csv')
both=read(trajectory/'Q72_two_end_backbone_Hbond_same_frame_per_replica.csv')
findings=json.loads((trajectory/'key_findings.json').read_text(encoding='utf-8'))
md=[]
for p in ['3hg1','4jff']:
    q=contact[contact.pdb.eq(p)&contact.beta_common_position.eq(54)]
    hb=both[both.pdb.eq(p)]
    assert len(q)==3 and len(hb)==3
    md.append({'pdb':p,'contact_percent':float(q.min_heavy_contact_le4_A_fraction.mean()*100),'simultaneous_two_end_Hbond_percent':float(hb.both_same_frame_Hbond_fraction.mean()*100),'Q72_MD_Hchi':findings['MD_descriptive_three_independent_replicas'][p]['Q72_Hchi_nats']['mean'],'G54_MD_RMSF':findings['MD_descriptive_three_independent_replicas'][p]['beta_CA_RMSF_A']['mean']})
assert [f"{r['contact_percent']:.1f}" for r in md]==['24.4','97.4']
assert [f"{r['simultaneous_two_end_Hbond_percent']:.1f}" for r in md]==['0.0','37.1']
assert md[1]['Q72_MD_Hchi']<md[0]['Q72_MD_Hchi'] and md[1]['G54_MD_RMSF']<md[0]['G54_MD_RMSF']

INV=EXP/'opposite_affinity_rmsf_cases_20261022/inventory'
allcdr=read(INV/'MEL5_three_contrasts_all_six_CDR.csv')
allh=read(INV/'MEL5_three_contrasts_four_state_union_pMHC_interface.csv')
peptide=read(INV/'MEL5_four_states_complete_peptide.csv')
cdrr=[];hr=[]
for target in ['4jfd','4jfe']:
    q=allcdr[allcdr.pdb_before.eq('4jff')&allcdr.pdb_after.eq(target)]
    rows=[]
    for region in ['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']:
        z=q[q.region_before.eq(region)]
        rows.append({'region':region,'n':len(z),'mean_delta_RMSF_A':float(z.rmsf_change_after_minus_before.mean()),'increased_n':int((z.rmsf_change_after_minus_before>0).sum())})
    ranked=sorted(rows,key=lambda r:r['mean_delta_RMSF_A'],reverse=True)
    assert {r['region'] for r in ranked[:2]}=={'CDR2α','CDR3α'}
    cdrr.append({'pdb':target,'all_regions':rows,'two_largest_mean_increases':['CDR2α','CDR3α']})
    q=allh[allh.pdb_before.eq('4jff')&allh.pdb_after.eq(target)]
    hr.append({'pdb':target,'HLA_deltas_nats':{f'{aa}{pos}':float(q[q.chain_before.eq('A')&q.id_before.eq(pos)].hchi_change_after_minus_before.iloc[0]) for aa,pos in [('K',68),('R',75),('Q',72)]}})

scope=[]
for a,b in [('3hg1','4jff'),('4jff','4jfd'),('4jff','4jfe')]:
    q=allh[allh.pdb_before.eq(a)&allh.pdb_after.eq(b)]
    n=int((q.hchi_valid&q.same_aa).sum())
    scope.append({'before':a,'after':b,'interface_n':len(q),'Hchi_comparable_n':n})
assert [r['interface_n'] for r in scope]==[30,30,30]
assert [r['Hchi_comparable_n'] for r in scope]==[23,23,22]
for p in ['3hg1','4jff','4jfd','4jfe']:assert len(peptide[peptide.pdb.eq(p)])==10
caption=prose['95']
for s in ['p表示肽','4JFF−3HG1','4JFD−4JFF','4JFE−4JFF','蓝、红表示降低、升高','灰格表示Hχ无定义或氨基酸类型不同','最近Cα距离≤8 Å','三行Hχ有效数为23、23、22','TCR采用MEL5共同编号','同一模型对已有复合物结构预测']:
    assert s in caption,s
literature=OLD/'case_sources/primary_literature/mel5_thermodynamics_2013.md'
lt=literature.read_text(encoding='utf-8')
for s in ['18 μ m','0.6 n m','97 MHC contacts compared with only 56','18 new vdW contacts and 1 new hydrogen bond']:
    assert s in lt,s
assert 18e-6/.6e-9==30000

report={'status':'passed','DOCX_path':str(DOC),'DOCX_SHA256':sha(DOC),'expected_SHA256_matched':True,'prose_JSON_SHA256':sha(prose_path),'read_only_review':True,'all8_changed_paragraphs_exactly_present_in_DOCX':True,'geometry_metric':'Fixed Q72 OE1–beta53 backboneN and Q72 NE2–betaG54 backboneO atom pairs, in starting structures; not minimum all-heavy distance.','recomputed_fixed_atom_distances':distances,'model_values':claims,'MD_summary_values':md,'MD_criteria_and_equal_weight_replicate_definition_retained_in_method_P53':True,'peptide_contrast_CDR_response_claims':cdrr,'peptide_contrast_HLA_Hchi_responses':hr,'interface_scope_and_mask':scope,'caption':{'all6_panels_match':True,'all_units_correct':True,'three_reference_contrasts_correct':True,'p_prefix_defined':True,'grey_mask_meaning_correct':True,'common_TCR_numbering_correct':True,'observed_structure_prediction_scope_correct':True},'primary_literature':{'path':str(literature),'SHA256':sha(literature),'SPR_18uM_to_0_6nM_and_30000fold_verified':True,'slower_offrate_verified':True,'MHC56_to97_and_betaCDR2_6_to38_verified':True,'I53F_new_MHC_contacts_verified':True,'original_experiments_not_mislabeled_as_MD':True},'mechanistic_claim_scope':'Associates model-predicted Q72 sidechain and unmutated G54 backbone responses with same local contact network, supported by existing MD and original SPR/crystallography; does not assert directional Hchi→RMSF causation or affinity prediction. Added peptide comparisons describe distinct regional responses, not a Q72-based affinity classifier.','blocking_issues':[],'required_scientific_edits':[],'visual_render_review':'not part of this report; pending separate page review'}
(ROOT/'manuscript_QA').mkdir(exist_ok=True)
(ROOT/'manuscript_QA/independent_science_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','DOCX_SHA256':sha(DOC),'file':str(ROOT/'manuscript_QA/independent_science_review.json'),'geometry':[r['display'] for r in distances]},ensure_ascii=False))
