from pathlib import Path
import json,pandas as pd,numpy as np,hashlib
R=Path(__file__).resolve().parent
O=R/'biological_annotation'
c=pd.read_csv(R/'shared/candidates_frozen56.csv')
abstracts=json.loads((O/'primary_abstracts_by_pdb.json').read_text(encoding='utf-8'))
official=json.loads((O/'official_mutation_metadata.json').read_text(encoding='utf-8'))
entries={e['rcsb_id']:e for e in official['response']['data']['entries']}
# Biological categories fixed from existing official metadata/primary abstract evidence, not per-member prediction quality.
CAT={
'1CLV':('folded_miniprotein',False,'原摘要明确为knottin样折叠、三股β片和二硫键，非线性识别片段。'),
'1DP5':('natural_inhibitory_fragment',True,'原摘要明确IA3的N端2–34抑制片段，游离态未检测到二级结构，结合态形成螺旋。'),
'1GL1':('folded_miniprotein',False,'原摘要明确三股反平行β片、三个二硫键的天然蛋白抑制剂。'),
'1MXE':('regulatory_recognition_fragment',True,'原摘要明确CaMKI 294–318调控区靶肽与CaM识别；不称完整激酶。'),
'1PPE':('folded_miniprotein',False,'原摘要明确29残基天然CMTI-I及三对二硫键、独立折叠。'),
'1SYQ':('regulatory_recognition_fragment',True,'原摘要明确talin VBS结合位点与vinculin头域相互作用。'),
'1T63':('motif_recognition',True,'原摘要明确核受体共激活因子的LXXLL/FXXLF招募基序识别。'),
'2D10':('motif_recognition',True,'原摘要明确NHERF末端13残基Motif-2与FERM结合。'),
'2I32':('functional_binding_fragment',True,'原摘要明确HIRA B区域与ASF1a特异结合，局部β发夹识别片段；不据此声称天然无序。'),
'2JBY':('motif_recognition',True,'原摘要将短链标识为与结合伙伴相互作用的BH3区域，官方短链身份对应；仅纳入既有结构评价。'),
'2OIN':('cofactor_binding_fragment',True,'原摘要明确既有结构中存在NS4A辅因子肽，官方短链身份一致；标记受体突变，非野生型泛化证据。'),
'2QXV':('functional_binding_fragment',True,'原摘要明确EZH2的30残基识别肽与EED结合。'),
'2XZE':('motif_recognition',True,'原摘要明确CHMP3 C端MIM与AMSH结合；该片段游离/结合状态证据来自原文而非本模型。'),
'3C5T':('natural_peptide_hormone_fragment',True,'原摘要明确Exendin-4(9–39)天然肽衍生的截短拮抗片段，不等同于完整天然激动肽。'),
'3C66':('functional_binding_fragment',True,'原摘要明确Fip1 80–105与poly(A) polymerase的结合片段。'),
'3DXC':('motif_recognition',True,'原摘要明确APP胞内尾与Fe65-PTB2的NPxY和相邻螺旋识别。'),
'3E2U':('folded_small_domain',False,'官方聚合物片段为Zn-knuckle 2，原摘要以C端锌指结构域为对象；观察到20残基不能视作天然短肽。'),
'3EJH':('functional_binding_fragment',True,'原摘要明确I型胶原α1来源肽与fibronectin结合。'),
'3G9W':('functional_binding_fragment',True,'原摘要明确integrin β尾与talin结合；天然长蛋白功能尾片段。'),
'3GJ8':('folded_small_domain',False,'原摘要明确Nup153的Ran结合锌指模块，非线性识别肽。'),
'3HTU':('functional_binding_fragment',True,'原摘要将VPS20首个螺旋映射为与VPS25相互作用的功能区域。'),
'3K6G':('functional_binding_fragment',True,'原摘要明确RAP1 RCT与TRF2结合；官方记录短链为TRF2功能片段。'),
'3OAK':('functional_binding_fragment',True,'原摘要明确Spn1与Spt6 binding determinant的复合结构。'),
'3ONW':('motif_recognition',True,'原摘要明确RGS14 GoLoco基序；受体存在已记录突变，作为突变敏感性单列。'),
'3SYN':('assembly_fragment_sensitivity',False,'官方短链为YlxH N端1–23；摘要讨论效应蛋白但未明确该截片能否代表独立天然识别肽，保守列敏感性。'),
'3TZ1':('regulatory_recognition_fragment',True,'原摘要明确短troponin I片段与troponin C末端域识别。'),
'3VUX':('folded_small_domain',False,'原摘要明确A20 ZF7识别diubiquitin，短链属于锌指结构域。'),
'3WA0':('motif_recognition',True,'原摘要明确DCAF1 C端IILXLN尾基序与merlin FERM识别。'),
'3WN7':('motif_recognition',True,'原摘要明确Nrf2 DLGex degron与Keap1识别。'),
'3ZRZ':('functional_binding_fragment',True,'原摘要明确来源于FnBP的Fn结合重复肽；仅纳入既有结构的数值评价。'),
'3ZWL':('functional_binding_fragment',True,'原摘要明确eIF3b C端α螺旋与eIF3i的保守功能性结合片段。'),
'4AOQ':('folded_miniprotein',False,'原摘要明确SOTI-III属于三二硫键knottin小蛋白，即使样品由合成获得也不改其天然折叠类别。'),
'4BQD':('designed_peptide',False,'原摘要明确经文库筛选和优化的compound 3设计肽；独立合成肽敏感性队列。'),
'4BRU':('motif_recognition',True,'官方短链Edc3与原摘要FDF结合位点相符；标题亦讨论Pat1但不将短链误注为Pat1。'),
'4BWQ':('motif_recognition',True,'原摘要明确PQBP1 C端YxxPxxVL结合基序。'),
'4J2L':('functional_binding_fragment',True,'原摘要明确CIC与ATXN1 AXH结合及功能性相互作用。'),
'4LZX':('motif_recognition',True,'原摘要明确IQCG的IQ基序与CaM识别。'),
'4MP0':('regulatory_recognition_fragment',True,'原摘要明确PNUTS的PP1结合片段与识别基序；重复记录先按PDB合并。'),
'4NIQ':('motif_recognition',True,'原摘要明确Vfa1末端17残基MIT-interacting motif 2。'),
'4REY':('functional_binding_fragment',True,'原摘要明确GM130 C端肽与GRASP65 PDZ结构域结合。'),
'4Z80':('annotation_uncertain',False,'官方聚合物采用Cytoadherence-linked asexual protein/EGF名称，原摘要采用另一组生物命名；未独立确认名称对应前仅宽队列保留。'),
'4ZRK':('functional_binding_fragment',True,'原摘要明确Lats1/2与merlin FERM的结合，短链为相应功能片段。'),
'4ZW2':('regulatory_recognition_fragment',True,'原摘要明确α1s I-II loop AID识别肽与DHPR β域相互作用。'),
'5B78':('chemically_modified_peptide',False,'原摘要明确组蛋白酰化修饰肽；规范化标准AA输入可能不完整表示化学状态，核心保守排除、宽队列保留。'),
'5ET1':('motif_recognition',True,'原摘要明确Myo3 THDI重复结合序列与Espin1识别。'),
'5F5O':('cofactor_binding_fragment',True,'官方身份及原摘要明确天然来源的辅因子结合肽；仅保留既有结构预测评价。'),
'5FZT':('motif_recognition',True,'原摘要明确DLC1 LD基序与talin R8识别。'),
'5KY5':('folded_small_domain',False,'原摘要明确EGF-like结构域是POFUT1底物，40残基不等于线性识别肽。'),
'5TGH':('functional_binding_fragment',True,'官方短链身份与原始复合物摘要匹配，作为既有功能结合片段，不展开调控或工程分析。'),
'5WUJ':('assembly_fragment_sensitivity',False,'原摘要以FliF胞质域与FliG N域装配为对象，短截片与线性识别基序不能直接等同。'),
'5YIR':('motif_recognition',True,'原摘要明确ankyrin来源Atg8结合序列；本文不能据此将所有构建体叫野生型，另查官方突变记录。'),
'6H9N':('functional_binding_fragment',True,'原摘要明确FtsB 64–87的24残基肽可重现对应相互作用，满足最小功能结合片段规则。'),
'6J4S':('observed_fragment_sensitivity',False,'官方SVBP沉积66残基而输入仅观察32；原摘要描述整段辅因子，核心不将部分可见片段当天然短肽。'),
'6JWJ':('motif_recognition',True,'原摘要明确Ufd1的Npl4结合基序。'),
}
assert set(c.pdb_id)==set(CAT)
rows=[]
for q in c.itertuples():
    cat,core,why=CAT[q.pdb_id]
    entry=entries[q.pdb_id]
    peptides=entry['polymer_entities']
    # Exact official description match; handle trailing whitespace only, not inferred biological aliases.
    exact=[p for p in peptides if p['rcsb_polymer_entity']['pdbx_description'].strip().casefold()==q.short_description.strip().casefold()]
    if not exact:
        # Candidate identities contain UniProt when available.
        uid=str(q.short_identity).removeprefix('UniProt:')
        exact=[p for p in peptides if uid in (p['rcsb_polymer_entity_container_identifiers'].get('uniprot_ids') or [])]
    match_rule='exact description or UniProt identity'
    if len(exact)!=1:
        exact=[p for p in peptides if q.short_sequence in p['entity_poly']['pdbx_seq_one_letter_code_can'].replace('\n','').replace(' ','')]
        match_rule='unique exact observed sequence substring in official polymer'
    ent=exact[0] if len(exact)==1 else None
    info=ent['rcsb_polymer_entity'] if ent else {}
    mut=[p['rcsb_polymer_entity'].get('pdbx_mutation') for p in peptides]
    mut=[x for x in mut if x and str(x).lower() not in ['none','no','?']]
    a=abstracts.get(q.pdb_id,{})
    seq=ent['entity_poly']['pdbx_seq_one_letter_code_can'].replace('\n','').replace(' ','') if ent else ''
    qseq=q.short_sequence
    # Sequence containment is a cross-check, not sequence mapping or label alignment.
    contain= qseq in seq if ent else False
    rows.append(dict(sample_id=q.sample_id,pdb_id=q.pdb_id,family_id=q.family_id,short_chain=q.short_chain,short_length=q.short_length,short_description=q.short_description,receptor_description=q.receptor_descriptions,category=cat,core_functional_peptide=core,reason=why,annotation_level='root manual review of official metadata plus primary abstract',primary_title=a.get('title',''),doi=a.get('doi',''),source_url=a.get('pubmed_url',''),rcsb_url='https://www.rcsb.org/structure/'+q.pdb_id,official_entity_matched=bool(ent),official_entity_id=ent.get('rcsb_id','') if ent else '',deposited_peptide_length=len(seq) if ent else None,observed_short_chain_is_substring_of_canonical_deposit=contain,official_fragment=info.get('pdbx_fragment'),official_peptide_mutation=info.get('pdbx_mutation'),any_official_polymer_mutation=bool(mut),all_official_polymer_mutations=' | '.join(mut),claim_native_disorder=False))
ann=pd.DataFrame(rows)
ann.to_csv(O/'all56_biological_annotations.csv',index=False)
ann.drop_duplicates('pdb_id').to_csv(O/'all54_PDB_biological_annotations.csv',index=False)
core=ann[ann.core_functional_peptide]
summary=dict(n_records=len(ann),n_PDBs=ann.pdb_id.nunique(),n_core_records=len(core),n_core_PDBs=core.pdb_id.nunique(),n_core_split_groups=core.family_id.nunique(),categories=ann.drop_duplicates('pdb_id').category.value_counts().to_dict(),n_official_entities_matched=int(ann.official_entity_matched.sum()),classification_independent_of_per_member_prediction=True,core_definition='Primary abstract/official polymer identifies a functional binding motif or recognized peptide fragment, natural peptide or inhibitory peptide fragment. No claim that the whole cohort is natively disordered.',exclusions_are_scope_not_model_failure=True,all_candidates_remain_in_wide_sensitivity=True,metadata_file_sha256=hashlib.sha256((O/'official_mutation_metadata.json').read_bytes()).hexdigest(),abstract_file_sha256=hashlib.sha256((O/'primary_abstracts_by_pdb.json').read_bytes()).hexdigest(),caveats=['Abstract-level classification; no claim full-text methods or all constructs verified.','Mutated receptors remain flagged; sensitivity excludes all official mutation-positive entries. Absence of mutation annotation does not prove wild type.','Observed chain length differs from deposited construct length; flags retained.','Public PDB auth/asym chain IDs are not assumed to match normalized A/B without evidence. Actual target labels use local exact chain mapping.'])
(O/'annotation_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
