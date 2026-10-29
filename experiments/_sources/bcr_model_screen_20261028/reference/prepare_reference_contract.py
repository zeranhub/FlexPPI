from pathlib import Path
from collections import OrderedDict
from datetime import datetime, timezone
import xml.etree.ElementTree as ET
import json,csv,hashlib
ROOT=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028')
OUT=ROOT/'reference';OUT.mkdir(exist_ok=True)
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_bcr_case_search_20261028')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse_pdb(code):
 d={}
 for l in (ROOT/'structures'/f'{code.lower()}_raw.pdb').read_text().splitlines():
  if l.startswith('ENDMDL'):break
  if l.startswith('ATOM  ') and l[16] in ' A':d.setdefault(l[21],OrderedDict()).setdefault(l[22:27].strip(),l[17:20])
 return d
pdbs={p:parse_pdb(p) for p in ['1C08','5AYU','6P4A','6P4B','6P4C','6P4D']}
AA3={'Q':'GLN','V':'VAL','T':'THR','S':'SER','W':'TRP','D':'ASP','H':'HIS','G':'GLY','F':'PHE','E':'GLU','R':'ARG','L':'LEU','I':'ILE','N':'ASN','Y':'TYR','C':'CYS'}
sites=[]
def addsite(contrast,role,number,aa,effect,source):
 ch={'H':'B','L':'A'}[role]
 actual=pdbs['1C08'][ch].get(str(number))
 sites.append({'case':'HyHEL10_2017','contrast':contrast,'paper_chain':role,'paper_number':number,'paper_aa':aa,'PDB_1C08_chain':ch,'PDB_1C08_author_id':str(number),'PDB_1C08_aa':actual,'identity_matches':actual==AA3[aa],'unchanged_in_contrast':True,'reported_effect':effect,'source_location':source,'can_validate_actual_model_label':False,'reference_scope':'Free-state millisecond NMR; qualitative cross-timescale hypothesis only, not RMSF/Hchi numerical truth.'})
for n,aa in [(3,'Q'),(24,'V'),(25,'T'),(28,'S'),(31,'S'),(36,'W'),(56,'T'),(68,'S'),(72,'D')]:
 addsite('free_R71A_minus_free_WT','H',n,aa,'NMR relaxation dispersion became flat; exchange contribution not detected under reported conditions.','PMC5577056 Results Par15 / Figure4e')
for n,aa in [(33,'H'),(84,'G'),(87,'F'),(98,'F')]:
 addsite('free_R45A_minus_free_WT','L',n,aa,'NMR relaxation dispersion no longer observed.','PMC5577056 Results Par16 / Figure4f')
for n,aa in [(46,'E'),(24,'V'),(25,'T'),(28,'S'),(31,'S'),(36,'W'),(72,'D')]:
 addsite('free_R45A_minus_free_WT','H',n,aa,'NMR relaxation dispersion no longer observed.','PMC5577056 Results Par16 / Figure4f')
with (OUT/'all_literature_named_NMR_sites.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(sites[0]));w.writeheader();w.writerows(sites)
# XML paragraph IDs are retrieved deterministically, avoiding inferred locators.
source_pars={}
for pid in ['PMC5577056','PMC7486785']:
 doc=ET.fromstring((TMP/(pid+'.xml')).read_bytes())
 hits=[]
 for e in doc.findall('.//body//p'):
  t=' '.join(''.join(e.itertext()).split())
  if any(k in t for k in ['Gln3, Val24','His33, Gly84','H:Q3A','L:H34A','0.38','1.1 Å','73 to 78','W34','Y53','three independent simulations','L4F','1,900-fold']):hits.append({'id':e.get('id'),'text':t})
 source_pars[pid]=hits
(OUT/'source_paragraphs.json').write_text(json.dumps(source_pars,ensure_ascii=False,indent=2),encoding='utf-8')
contract={
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'analyst_has_seen_prediction_outputs':False,
 'status':'Frozen literature-derived scope, not an outcome-selected case selection or a claim that the model will recover it.',
 'scope':'Only explicit names/regions and effect types in the two supplied primary papers. Figure color maps and the 32 free-WT NMR sites are not all numerically tabulated; no missing numeric site list invented.',
 'general_rules':[
  'Retain complete shared variable-region chains, all six CDRs under a separately verified numbering scheme, and full antigen residue profiles; no best-site filtering.',
  'For change Hchi comparisons require same amino-acid identity and chemically defined chi on both sides. Changed sites, Ala/Gly absence and missing coordinates remain in raw output but cannot be called entropy recovery.',
  'A static rotamer change, shorter C-alpha distance or crystal B factor is not Hchi or RMSF validation. Published free-state NMR millisecond exchange is not the model nanosecond-MD target.',
  'Only genuine independent MD/reference RMSF and Hchi on matched definitions can evaluate numeric prediction accuracy. Current paper-derived checks below are qualitative or structural, except the published FlexHEL antigen RMSF trend.',
  'If mutant input reuses WT backbone and only changes amino-acid identity, label it fixed-backbone sequence-only sensitivity. It is not an independent mutant crystal or the paper equilibrated MD state.',
  'Checkpoint/family and state selection must be based on structural/sequence provenance, not prediction outcome. Crystal copies and checkpoints are not biological repeats.'
 ],
 'structure_map':{
  'HyHEL10_2017':{'WT_bound':'1C08: A=VL, B=VH, C=HEL','WT_free':'5AYU: L=VL, H=VH','mutant_crystals':'None reported in the paper; R71A/R45A MD starts from PyMOL mutations of1C08.','PDB_mutation_checks':{'H:R71':'1C08 B71 ARG; 5AYU H71 ARG','L:R45':'1C08 A45 ARG; 5AYU L45 ARG','H:Q3':'1C08 B3 GLN; 5AYU H3 GLN','L:H34':'1C08 A34 HIS; 5AYU L34 HIS'}},
  'FlexHEL_2020':{'6P4A':'WT Hy10 Fab L/H + rigid HEL2x chain C (R21Q/R73E; D101 remains Asp)','6P4B':'Hy10-4x Fab L/H + flex antigen C, plus crystal copy A/B + D. Heavy changesL4F/Y33H/S56N/Y58F; antigenC76S/C94S relative rigid.','6P4C':'Free Hy10-4x Fab L/H','6P4D':'Free HEL3x A: R21Q/R73E/D101R; NOT sequence-identical free counterpart of the D101Asp antigens in6P4A/B.','paper_Fig5_ID_warning':'Original Fig5 PDB text is inconsistent with current deposited identities; use sequence and current RCSB metadata, not its ID labels alone.'}
 },
 'contrasts':[
  {'case':'HyHEL10_2017','state_comparison':'free H:R71A minus free WT','reported_effect':'Reduction of detectable millisecond backbone exchange on the complete reported unchanged VH site set; VL reported little affected.','unchanged_residue_eligibility':{'all_reported':[f"{x['paper_chain']}:{x['paper_aa']}{x['paper_number']}" for x in sites if 'R71A' in x['contrast']],'exact_identity_evaluable':[f"{x['paper_chain']}:{x['paper_aa']}{x['paper_number']}" for x in sites if 'R71A' in x['contrast'] and x['identity_matches']],'unresolved':['H:Thr56: 1C08 B56 isSer, B57 isThr. Do not silently shift mapping.'],'changed_site':'H71 Arg→Ala excluded from same-type Hchi comparison.'},'source_location':'PMC5577056 Results Par15 / Figure4c,e','can_validate_actual_model_label':False,'validation_kind':'Qualitative cross-timescale; no expected Hchi direction fixed. RMSF decrease would be a hypothesis, not exact NMR-label validation.'},
  {'case':'HyHEL10_2017','state_comparison':'free L:R45A minus free WT','reported_effect':'Detectable exchange lost at VL–VH interface sites and additional VH outer-loop/CDR sites.','unchanged_residue_eligibility':{'all_reported':[f"{x['paper_chain']}:{x['paper_aa']}{x['paper_number']}" for x in sites if 'R45A' in x['contrast']],'exact_identity_evaluable':[f"{x['paper_chain']}:{x['paper_aa']}{x['paper_number']}" for x in sites if 'R45A' in x['contrast'] and x['identity_matches']],'unresolved':['L:His33: 1C08 A33 isLeu, A34 isHis. Do not silently shift or equate it to theL:H34A affinity-control site.'],'changed_site':'L45 Arg→Ala excluded from same-type Hchi comparison.'},'source_location':'PMC5577056 Results Par16 / Figure4d,f','can_validate_actual_model_label':False,'validation_kind':'Qualitative cross-timescale; Gly84 is valid backbone site but has no side-chain chi.'},
  {'case':'HyHEL10_2017','state_comparison':'free WT versus bound WT','reported_effect':'Most of32 free-WT exchanging residues lose observable relaxation dispersion after binding; free clusters kex840±15 and2033±50 s^-1.','unchanged_residue_eligibility':'Full shared VH/VL and all six CDRs; named H71 andL45 are examples, not a complete32-site list. VL–VH remain associated in Fv, so this is internal dynamics.','source_location':'PMC5577056 Results Par9,Figure3','can_validate_actual_model_label':False,'validation_kind':'No numeric32-site table extracted or invented; qualitative free-to-bound consistency only.'},
  {'case':'HyHEL10_2017','state_comparison':'bound H:R71A minus bound WT','reported_effect':'130ns MD median C-alpha distance H:I29–H:S74 is0.38Angstrom shorter; outer-loop/CDR packing closer.','unchanged_residue_eligibility':['H:I29','H:S74'],'PDB_pairs':[['1C08:B:29:CA','1C08:B:74:CA']],'source_location':'PMC5577056 Results Par18 andFigure5a,b; Methods Sec15','can_validate_actual_model_label':False,'validation_kind':'Distance/packing reference, not RMSF or Hchi. No per-residue MD-RMSF direction inferred from reduced distance.'},
  {'case':'HyHEL10_2017','state_comparison':'bound L:R45A minus bound WT','reported_effect':'130ns MD median C-alpha distance L:H41–H:G42 is1.1Angstrom shorter; VL–VH packing closer. WT L:Arg45 side chain described as flipping between domains.','unchanged_residue_eligibility':['L:H41','H:G42'],'PDB_pairs':[['1C08:A:41:CA','1C08:B:42:CA']],'source_location':'PMC5577056 Results Par18,Figure5c,d; Methods Sec15','can_validate_actual_model_label':False,'validation_kind':'Distance/packing only. Changed L45 Arg→Ala cannot provide same-type Hchi; H:G42 has no chi.'},
  {'case':'HyHEL10_2017','state_comparison':'WT, R71A, R45A, Q3A, H34A affinity comparison','reported_effect':{'WT_KA_M_inverse':1.8e8,'H_R71A_KA_M_inverse':2.0e9,'L_R45A_KA_M_inverse':2.0e9,'H_Q3A_KA_M_inverse':2.2e8,'L_H34A_KA_M_inverse':2.8e8},'unchanged_residue_eligibility':'Compare full shared unchanged VH/VL profiles; never use Ala changed-site Hchi absence as positive signal.','source_location':'PMC5577056 Results Par11–12; SupplementaryFigure2','can_validate_actual_model_label':False,'validation_kind':'Affinity is contextual, not FlexPPI target. Q3A andH34A are approximately unchanged-affinity controls, not experimentally established no-dynamic-change controls. Q3A unstable for further fluctuation analysis; H:R44A aggregates and cannot be an affinity-performance negative.'},
  {'case':'FlexHEL_2020','state_comparison':'free matched rigid antigen versus free matched flex antigen','reported_effect':'Three independent simulations report greater per-residueRMSF in antigen73–78 for flex, with overall fold stable.','unchanged_residue_eligibility':{'primary_local_loop':[73,74,75,76,77,78],'mutated_in_loop':[76],'all_antigen_profile':'all shared protein residues; retain all129 positions, with same-typeHchi limited by identity/chemistry','isogenic_sequence_rule':'Use identical R21Q/R73E andD101Asp backgrounds; onlyC76/C94→Ser differ.6P4D hasD101Arg and is not an isogenic rigid input.'},'source_location':'PMC7486785 Figure3A andSI Appendix FigureS1,MD methods inSI','can_validate_actual_model_label':'Qualitative regionalRMSF trend only until raw trajectories and model-specific alignment definitions are verified.','validation_kind':'Published MD supports antigenRMSF73–78; not an antibody/Fab or complex Hchi truth; no side-chain entropy direction reported.'},
  {'case':'FlexHEL_2020','state_comparison':'6P4A WT-Fab/rigid versus6P4B 4x-Fab/flex','reported_effect':'Antigen73–78 loop and antibodyCDRH1 reorient; remote framework/side-chain network linksH4/H33 throughH34 toH29, andH58/H56 throughAg101/H53/Ag62 to the loop.','unchanged_residue_eligibility':{'primary_unmutated_sidechain_network':['FabH:W34','FabH:I29','FabH:Y53','antigen:D101','antigen:W62'],'additional_reported_contact':['antigen:K97','FabH:53 backbone'],'changed_sites_do_not_compare_same_typeHchi':['FabH:L4F','FabH:Y33H','FabH:S56N','FabH:Y58F','antigen:C76S','antigen:C94S'],'backbone_regions':'All sixCDRs under a separately verified scheme, entire sharedFvH/L and full antigen, with antigen73–78 prespecified local region.'},'source_location':'PMC7486785 Figure5B–D and Results structural paragraphs','can_validate_actual_model_label':False,'validation_kind':'Crystal rotamer/backbone conformation evidence only. Both antigen and antibody sequences change, so not a one-site causal control. No Hchi direction inferred from staticrotamer differences.'},
  {'case':'FlexHEL_2020','state_comparison':'6P4C free4x-Fab versus6P4B bound4x-Fab','reported_effect':'Free4xCDRH1 broadly resembles WT rigid-complex conformation rather than the flex-bound arrangement; localH4/H34 adjustments retained.','unchanged_residue_eligibility':'Same4x antibody sequence on full sharedFv domains and all sixCDRs; match all common exact amino-acid identities.','source_location':'PMC7486785 Results final structural paragraph immediately beforeDiscussion,SI AppendixFigureS8','can_validate_actual_model_label':False,'validation_kind':'Independent free/bound crystals support structural context, not directly ensemble entropy; absent matched apo/bound rawMD cannot assess RMSF/Hchi accuracy.'}
 ],
 'fixed_analysis_reporting':{
  'HyHEL10_main_mutations':['H:R71A','L:R45A'],'HyHEL10_affinity_controls':['H:Q3A','L:H34A'],
  'full_scope':'Retain every common residue output; report full sharedVH/VL, all sixCDRs using verified numbering and all antigen residues before selected local discussion.',
  'direction_reporting':'Report raw deltas on all predefined sets without scoring millisecondNMR or staticpacking as quantitative MD-targetaccuracy. No sign assigned forHchi when literature lacks it.',
  'chemical_exclusion':'Mutation-created Ala/Gly or changed amino-acid type must not create apparent same-typeHchi agreement.',
  'technical_copy_rules':'6P4B two crystallographic copies may be sensitivity inputs only; do not count them as independent simulation/biological replicates.',
  'missing_reference':'If raw MD orNMR per-site numeric values unavailable, report that absence, not fabricated reference data or agreement.'
 },
 'source_hashes':[{'path':str(p),'sha256':sha(p)} for p in [TMP/'PMC5577056.xml',TMP/'PMC7486785.xml']]+[{'path':str(ROOT/'structures'/f'{p.lower()}_raw.pdb'),'sha256':sha(ROOT/'structures'/f'{p.lower()}_raw.pdb')} for p in pdbs]
}
(OUT/'literature_validation_contract.json').write_text(json.dumps(contract,ensure_ascii=False,indent=2),encoding='utf-8')
md='''# BCR候选的预定义验证范围\n\n本文件在分析者未查看模型输出时固定文献范围，不预设候选会成功。全部共享残基、完整Fv重轻链、六个经独立编号的CDR和完整抗原输出均须保留；文献局部位点只作预先指定核查，不取代完整结果。数值与完整范围见literature_validation_contract.json，逐位点核查见all_literature_named_NMR_sites.csv。\n\nHyHEL10的主要单点为H:R71A和L:R45A，近似无亲和变化的对照为H:Q3A和L:H34A。后两者不是已经证明“无动态变化”的对照；Q3A由于稳定性不足未继续做动态分析，H:R44A则聚集而不能作真实亲和力阴性。上述变体均没有该论文报告的独立突变晶体；若复用1C08或5AYU骨架只改序列，只能称固定骨架敏感性。\n\n2017论文中游离态NMR明确点名的R71A未替换响应位为H:Q3/V24/T25/S28/S31/W36/T56/S68/D72，共9位；R45A响应位为L:H33/G84/F87/F98及H:E46/V24/T25/S28/S31/W36/D72，共11位。这里H:T56在1C08和5AYU同号为Ser，L:H33同号为Leu，存在AA/编号不一致，两项应标为映射未解决，不能静默改成H57或L34。另有被改变的H71、L45本身，不能用于同型侧链熵比较。全32个游离WT交换位点未从正文得到数值表，不能捏造完整逐位真值。\n\nNMR是游离态毫秒交换，与模型基于MD的RMSF/Hχ并非同一标量、同一时间尺度。论文结合态130nsMD明确给出的主要局部证据是H:I29–H:S74的Cα距离收缩0.38Å，以及L:H41–H:G42距离收缩1.1Å；更紧的堆积本身不等于RMSF或Hχ下降。没有原始轨迹时，只能将模型输出作为定性假设核查，不能宣称已测得定量预测准确率。\n\nFlexHEL的当前PDB对应为6P4A=WT Hy10+rigid抗原、6P4B=4xHy10+flex抗原、6P4C=游离4xFab、6P4D=游离HEL3x。原文Fig5的ID对应与当前结构不符，需用当前序列核实。6P4D的101位是Arg，6P4A/B抗原101位是Asp，不能把6P4D直接当作序列相同的rigid apo。6P4B含两份晶体copy，可作敏感性但不能算独立重复。\n\nFlexHEL发表的三重复MD支持游离抗原73–78环在flex中RMSF增大，不是抗体或复合物的MD真值。侧链网络预定义为未替换的抗体H:W34/I29/Y53、抗原D101/W62，另保留抗原K97与抗体H53骨架联系；四个抗体突变及C76S/C94S不作同型Hχ比较。这些网络来自晶体构象重排，没有侧链熵方向先验。6P4A→6P4B同时改变抗体和抗原，属于状态比较；6P4C→6P4B则可比较同一4x抗体的游离/结合结构。\n'''
(OUT/'README_验证范围.md').write_text(md,encoding='utf-8')
print(json.dumps({'contract':str(OUT/'literature_validation_contract.json'),'named_NMR_sites':len(sites),'identity_matched':sum(s['identity_matches'] for s in sites),'unresolved':[f"{s['paper_chain']}{s['paper_number']} {s['paper_aa']} vs {s['PDB_1C08_aa']}" for s in sites if not s['identity_matches']],'contract_sha256':sha(OUT/'literature_validation_contract.json')},ensure_ascii=False,indent=2))
