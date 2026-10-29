import csv
import json
from pathlib import Path
from datetime import datetime

BASE = Path(__file__).resolve().parent
META = BASE.parent / 'all244_ITrace_metadata_and_sequences.csv'
metadata = {row['pdb'].lower(): row for row in csv.DictReader(META.open(encoding='utf-8-sig'))}

# Biology-first order was established without opening prediction or scoring files.
cases = [
    dict(case_id='B01', priority='A1', system='MAGE-A10: parental c728 to c796',
         pdbs='7pdw;7pbc', contrasts='7pdw>7pbc',
         controlled='Same HLA-A*02:01 and GLYDGMEHL peptide; parental versus direct derivative TCR',
         biological_change='Published beta-CDR1 E31D substitution; about sixfold affinity increase with retained target selectivity',
         evidence='SPR and previously published functional X-scan selectivity; 2022 structural explanation',
         source_doi='10.1136/jitc-2022-004600', source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC9295655/',
         source_local='sources/magea10_2022_exa.txt',
         prespecified_regions='Whole peptide; all six CDRs separately and pooled; unchanged-AA pMHC interface',
         caveat='7QPJ c756 is from a different parent, not a second derivative of c728. N-terminal tags and chain naming differ. Apo c728 crystal 7PDX exists in the source paper but is not in ITrace 244; apo MD was not identified. Solvent exclusion is not residue Hchi.',
         why_suitable='Single TCR change perturbs the same antigen environment; complete biological contrast with modest structural rearrangement'),
    dict(case_id='B02', priority='A2', system='AHIII12.2: single MHC substitutions with thermodynamic and cell-response controls',
         pdbs='1lp9;2j8u;2jcc;2uwe', contrasts='1lp9>2j8u;1lp9>2jcc;1lp9>2uwe',
         controlled='Same AHIII12.2 TCR and ALWGFFPVL peptide; HLA-A2 WT versus K66A/W167A/T163A',
         biological_change='SPR equilibrium Kd WT 9.3, K66A 31.8, W167A 15.4, T163A 4.7 uM; K66A has major CDR3alpha rearrangement and lowest cytolysis',
         evidence='SPR all four conditions; cytolytic response panel; direct ITC WT versus K66A',
         source_doi='10.1016/j.jmb.2007.07.028', source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC2065754/',
         source_local='sources/ahiii_mhc_2007_exa.txt',
         prespecified_regions='Whole peptide; all six CDRs separately and pooled; unchanged-AA pMHC interface; all three variants retained',
         caveat='2J8U=K66A; 2JCC=W167A; 2UWE=T163A. T163A and W167A provide essential controls. No apo trajectory identified in ITrace inventory. Bound Hchi does not equal binding entropy.',
         why_suitable='Direct MHC-environment change at fixed TCR/peptide, with a strong alteration and milder controls rather than one favorable pair'),
    dict(case_id='B03', priority='A3', system='1G4: NY-ESO-1 peptide C9V anchor substitution',
         pdbs='2bnr;2bnq', contrasts='2bnr>2bnq',
         controlled='Same 1G4 TCR and HLA-A2; peptide SLLMWITQC versus SLLMWITQV',
         biological_change='TCR Kd 13.3 to 5 uM at 25 C; analogue also has stronger HLA binding and faster lytic-granule polarization',
         evidence='SPR, peptide-HLA binding, cellular synapse/granule and cytotoxicity observations',
         source_doi='10.1084/jem.20042323', source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC2213140/',
         source_local='../../tcr_expanded_screen_v41_20261017/biology_sources/nyeso_vaccine_2005.md',
         prespecified_regions='Whole peptide including separate C9V position; unchanged peptide positions; all six CDRs; MHC interface',
         caveat='P9 chemistry changes directly; main evidence for environmental prediction must exclude P9 Hchi. Stronger binding does not imply every bound residue becomes rigid.',
         why_suitable='Single cancer-antigen anchor change with clear published functional and thermodynamic consequences'),
    dict(case_id='B04', priority='A4', system='N17.1.2: NRAS Q61K versus Q61R',
         pdbs='8yiv;8yj2', contrasts='8yiv>8yj2',
         controlled='Same N17.1.2 TCR and HLA-A*01:01; ILDTAGKEEY versus ILDTAGREEY, peptide P7 K/R',
         biological_change='SPR Kd 1.2 +/- 0.1 versus 3.4 +/- 0.2 uM; TCR cross-recognizes both neoantigens but not WT NRAS',
         evidence='SPR and prior peptide-dose functional reactivity; co-crystal structures',
         source_doi='10.1126/sciadv.adq6150', source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11584006/',
         source_local='../../tcr_expanded_screen_v41_20261017/biology_sources/nras_neoantigen_2024.md',
         prespecified_regions='All nine unchanged peptide positions; substituted P7 separately; all six CDRs; unchanged pMHC interface',
         caveat='This pair is two recognized neoantigens, not WT-versus-mutant escape. P7 K/R entropy cannot alone establish contextual prediction.',
         why_suitable='Same TCR and MHC with one side-chain replacement and quantitative modest biological difference'),
    dict(case_id='B05', priority='A5', system='PRAME: S2 TCR across HLA-A3 and HLA-A11',
         pdbs='8rym;8ryn', contrasts='8rym>8ryn',
         controlled='Same S2 WT TCR and ELFSYLIEK PRAME peptide; HLA-A*03:01 versus HLA-A*11:01',
         biological_change='SPR Kd 21 versus 72 uM; HLA position 163 changes T/R and is the directly contacted polymorphic determinant',
         evidence='SPR and structural comparison; published HLA-mutant mechanistic controls',
         source_doi='10.1038/s41598-024-75140-7', source_url='https://www.nature.com/articles/s41598-024-75140-7',
         source_local='sources/alloselectivity_2024_exa.txt',
         prespecified_regions='Whole identical peptide and TCR variable domains; all six CDRs; unchanged MHC-interface residues',
         caveat='HLA alleles differ at multiple sites, not only 163. Do not treat 8RYN as an engineered S2 TCR.',
         why_suitable='Exactly fixed antigen/TCR sequence isolates changed presenting environment, with quantified alloselectivity'),
    dict(case_id='B06', priority='A6', system='4TCR2: KRAS G12V affinity engineering',
         pdbs='8wte;8wul', contrasts='8wte>8wul',
         controlled='Same HLA-A*11:01 and VVGAVGVGK peptide; 4TCR2 parent versus published MH derivative',
         biological_change='Published two TCR-beta substitutions yield about 15-fold affinity gain; MH Kd about 2 uM and retains discrimination from WT peptide',
         evidence='SPR and co-crystal structures',
         source_doi='10.1038/s42003-024-06209-2', source_url='https://www.nature.com/articles/s42003-024-06209-2',
         source_local='../../tcr_expanded_screen_v41_20261017/biology_sources/kras_g12v_enhance_2024.md',
         prespecified_regions='Whole identical peptide; all six CDRs; unchanged pMHC interface',
         caveat='The current sequence inventory also has constant-domain truncation/missing-residue differences. Align variable domains carefully and do not interpret missing residues as engineered changes.',
         why_suitable='Cancer neoantigen with a true parent-to-derivative contrast and measured affinity effect'),
    dict(case_id='B07', priority='A7', system='KRAS G12D: engineered TCR recognizes WT versus mutant peptide',
         pdbs='7ow5;7ow6', contrasts='7ow5>7ow6',
         controlled='Same engineered TCR and HLA-A*11:01; VVVGAGGVGK versus VVVGADGVGK at P6',
         biological_change='JDIa41b1 Kd about 3.0 uM on WT versus 743 +/- 18 pM on G12D; original study reports thermodynamics and MD context',
         evidence='Biochemical affinity/specificity and cellular activity in primary 2022 study',
         source_doi='10.1038/s41467-022-32811-1', source_url='https://www.nature.com/articles/s41467-022-32811-1',
         source_local='sources/kras_g12d_2022_exa.txt',
         prespecified_regions='All nine unchanged peptide positions; P6 separately; all six CDRs; MHC interface',
         caveat='Structure is the engineered TCR bound to both peptides, not the parent TCR. 7PB2 is parent-JDI on mutant peptide and is a separate engineering contrast.',
         why_suitable='Single tumour mutation at identical TCR/MHC is a direct environment-prediction question'),
    dict(case_id='B08', priority='B1', system='PRAME: S2 WT versus S2-198 on HLA-A3',
         pdbs='8rym;8ryo', contrasts='8rym>8ryo',
         controlled='Same HLA-A*03:01 and ELFSYLIEK peptide; S2 WT versus published CDR1alpha-engineered S2-198',
         biological_change='Kd 21 uM to 76 nM for A3; A3/A11 affinity preference widens, contrary to engineering objective',
         evidence='SPR, HLA-mutant controls, co-crystal structure',
         source_doi='10.1038/s41598-024-75140-7', source_url='https://www.nature.com/articles/s41598-024-75140-7',
         source_local='sources/alloselectivity_2024_exa.txt',
         prespecified_regions='Whole peptide; all six CDRs; unchanged pMHC interface',
         caveat='8RYP versus 8RYQ changes both S8 TCR and HLA allele; it is not a clean S8 parent/engineered pair on fixed HLA and is lower priority.',
         why_suitable='Same antigen environment with defined engineering and a biologically informative specificity outcome'),
    dict(case_id='B09', priority='B2', system='MEL5: native MART-1 decamer versus nonamer',
         pdbs='4qok;6eqa', contrasts='4qok>6eqa',
         controlled='Same MEL5 TCR and HLA-A*02:01; native EAAGIGILTV versus AAGIGILTV',
         biological_change='MEL5 can induce alternate peptide anchoring in AAG; NMR and thermodynamics support two N-terminal configurations in bound AAG',
         evidence='Crystallography, peptide thermodynamic comparison, NMR and rigidity analysis in 2019 study',
         source_doi='10.1002/eji.201948085', source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC6618058/',
         source_local='sources/mart1_anchor_2019_exa.txt',
         prespecified_regions='Entire peptide with explicit sequence alignment; all six CDRs; unchanged MHC interface',
         caveat='Peptide length differs; register must be explicitly aligned. The static model input may contain only one of alternate crystallographic occupancies. 6EQB is engineered alpha24beta17 bound to AAG, not an apo state.',
         why_suitable='Existing orthogonal evidence explicitly concerns molecular plasticity, close to the model outputs'),
    dict(case_id='B10', priority='B3', system='DMF5: MART-1 nonamer versus heteroclitic decamer, with DMF4 control',
         pdbs='3qdg;3qdj;3qdm;3qeq', contrasts='3qdg>3qdj;3qdm>3qeq',
         controlled='Within each TCR, same HLA-A2; ELAGIGILTV versus AAGIGILTV',
         biological_change='DMF5 Kd 6 versus 40 uM; DMF4 29 versus 170 uM. DMF5 engages both with a shared open architecture; DMF4 changes docking orientation.',
         evidence='SPR and paired co-crystal structures; free DMF5 structure in source study',
         source_doi='10.4049/jimmunol.1101268', source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC3166883/',
         source_local='../../tcr_expanded_screen_v41_20261017/biology_sources/dmf5_crossreactivity_2011.md',
         prespecified_regions='Whole aligned peptide and six CDRs for both receptor series; unchanged MHC interface',
         caveat='Length and anchor identity differ together. Avoid comparing only DMF4 versus DMF5 and calling the difference caused by one molecular change.',
         why_suitable='A mechanistic contrast can be checked across two separately controlled receptor series'),
    dict(case_id='B11', priority='B4', system='MEL5: recognition of three endogenous tumor antigens',
         pdbs='4qok;7q99;7q9a', contrasts='4qok>7q99;4qok>7q9a;7q99>7q9a',
         controlled='Same MEL5 TCR and HLA-A2; Melan-A EAAGIGILTV, IMP2 NLSALGIFST, BST2 LLLGIGILVL',
         biological_change='MEL5 recognizes all three epitopes, with functional activation hierarchy Melan-A > BST2 > IMP2; shared central recognition motif despite different peptide sequences',
         evidence='Cell 2023 biochemical/functional study plus three TCR-bound antigen structures',
         source_doi='10.1016/j.cell.2023.06.020', source_url='https://pubmed.ncbi.nlm.nih.gov/37490916/',
         source_local='sources/multi_tumor_2023_exa.txt;sources/mel5_multi_affinity_exa.txt',
         prespecified_regions='Entire peptide with changed and unchanged positions separated; all six CDRs; unchanged MHC interface',
         caveat='Multiple peptide substitutions. 7Q9B is MEL8, not MEL5. Strong antigen-knockout tumor-cell evidence in this article is for MEL8 and must not automatically be attributed to MEL5.',
         why_suitable='Biologically substantive cross-recognition, but chemical differences make it less clean than single-change pairs'),
    dict(case_id='B12', priority='B5', system='1G4: CDR2 affinity maturation',
         pdbs='2f54;2f53', contrasts='2f54>2f53',
         controlled='Published WT 1G4 versus c49c50 CDR2-alpha/beta derivative; same HLA-A2 and SLLMWITQC peptide',
         biological_change='c49c50 soluble affinity about 1 nM with retained peptide specificity in tested panels',
         evidence='SPR and cell-binding specificity assays; paired structures',
         source_doi='10.1110/ps.051936406', source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC2242494/',
         source_local='../../tcr_expanded_screen_v41_20261017/biology_sources/1g4_cdr2_affinity_2006.md',
         prespecified_regions='Whole peptide; all six CDRs, not only edited CDR2; unchanged MHC interface',
         caveat='Current ITrace alpha sequences contain missing-residue/construct differences outside engineered CDR2. Validate correspondence against SEQRES/author numbering; consider 2BNR as secondary WT construct control. Affinity is not cellular killing.',
         why_suitable='Known parent-to-derivative relationship, but sequence provenance is less clean in current inventory'),
]

for case in cases:
    pdbs = case['pdbs'].split(';')
    assert all(pdb in metadata for pdb in pdbs), case['case_id']
    case['in_ITrace_244'] = 'all_present'
    case['model_scores_used_to_select'] = 'no'
    case['peptide_sequences_in_inventory'] = ';'.join(f'{p}:{metadata[p]["peptide_sequence"]}' for p in pdbs)
    case['PDB_primary_citations'] = ';'.join(f'{p}:{metadata[p]["DOI"]}' for p in pdbs)

fields = list(cases[0])
with (BASE / 'biology_first_ITrace_candidate_series.csv').open('w', encoding='utf-8-sig', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(cases)

audit = {'created': datetime.now().isoformat(timespec='seconds'),
         'selection_basis': 'Primary published biology, controlled variables and availability in ITrace metadata, without reading model predictions.',
         'source_metadata': str(META),
         'candidate_series_count': len(cases),
         'unique_pdbs': sorted({p for case in cases for p in case['pdbs'].split(';')}),
         'manuscript_modified': False,
         'new_variants_designed': False,
         'new_MD_or_training_started': False,
         'Jina_status': 'Direct Jina connection unavailable; agent-reach Exa returned primary-paper excerpts, supplemented by primary-paper web reads and existing local full text.'}
(BASE / 'selection_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')

report = '''# ITrace 生物学优先候选盘点

本轮仅按已发表生物学差异、变量是否清楚、ITrace 是否已有结构和轨迹建立候选顺序，没有读取模型预测值或按某个环的相关系数挑选。12 组系列涉及的 PDB 均在当前 244 个 ITrace 条目中。这个顺序表示优先复算价值，不表示模型已经预测成功。训练完全不含 ITrace 的既有模型应在候选与统计区域确定后统一测试；论文不作改动。

最优先的是 **MAGE-A10 的 7PDW→7PBC**。它比较同一亲本 c728 与直接衍生的 c796，对应已发表的 βCDR1 E31D 单点变化；pMHC 与肽 GLYDGMEHL 保持相同，亲和力约增强 6 倍。这样可以具体检验：仅改变受体的一处已有变体，同一抗原周围未改变残基的主链和侧链状态能否被模型辨认。[原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9295655/)。**7QPJ 的 c756 来自另一个亲本**，不能把它写成 c728 的第二代双点改造。c756 与 c796 的 X-scan 选择性差别可作为额外比较，但必须保留不同背景的事实。

第二优先的是 **AHIII12.2 四状态系列：1LP9 为 WT，2J8U 为 MHC K66A，2JCC 为 W167A，2UWE 为 T163A**。同一 TCR 和 ALWGFFPVL 肽保持不变。SPR 的平衡解离常数依次为 9.3、31.8、15.4、4.7 μM；K66A 具有显著 CDR3α 构象改变且细胞杀伤最低，而 T163A、W167A 提供不同程度的对照。WT/K66A 另有直接 ITC 比较，尤其适合联系结构环境变化与识别环运动。[原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2065754/)。当前 ITrace 盘点未发现相应游离 TCR 状态，不能拿两个结合态的侧链熵直接解释结合熵。原文关于 K66A 热效应方向的正文和符号表述需在引用具体热值前另行核查，现阶段只采用“焓贡献显著减小”的定性结果。

另一类干净的对照是同 TCR、同 MHC 的单肽位点改变。**1G4 的 2BNR→2BNQ** 为 NY-ESO-1 肽 C9V，实测 TCR 解离常数从 13.3 降至 5 μM，并有杀伤颗粒极化等功能证据；锚定位点变化还能传递到未改变的识别界面，不能只看 C9/V9 本身的熵。[原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2213140/)。**NRAS 的 8YIV→8YJ2** 为相同 N17.1.2 受体识别 P7 Lys/Arg 两种新抗原，解离常数为 1.2 与 3.4 μM；这是两种均可识别抗原的差别，并非“野生型不识别”的结构成对证据。[原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11584006/)。**KRAS 的 7OW5→7OW6** 也有固定 TCR/MHC 的 WT 肽与 G12D 肽对照，应与亲本 JDI 的 7PB2 分开。

本轮新增加的 **PRAME 8RYM/8RYN/8RYO** 很值得测。8RYM→8RYN 保持 S2 WT 受体及 ELFSYLIEK 肽相同，只更换 HLA-A3/A11，实测亲和力为 21/72 μM；8RYM→8RYO 则固定 A3，比较 WT 与 S2-198，后者亲和力达 76 nM。原文还指出工程后对 A3 的偏好反而加大，给出了比“增强亲和力”更明确的机制问题。[原始论文](https://www.nature.com/articles/s41598-024-75140-7)。8RYP→8RYQ 同时改变 S8 受体和 HLA，不能当成固定 HLA 的亲本改造对照。

其余已按生物学关系登记的系列包括 8WTE→8WUL（KRAS G12V 的 4TCR2 亲本/已发表 MH 衍生体）、4QOK→6EQA（MEL5 识别 MART-1 十肽/九肽与锚定改变）、3QDG→3QDJ 并行 3QDM→3QEQ（DMF5/DMF4 对同一组十肽/九肽的不同识别方式）、4QOK/7Q99/7Q9A（MEL5 对 Melan-A、IMP2、BST2 的多肿瘤抗原识别）以及 2F54→2F53（1G4 的 CDR2 亲和成熟）。7Q9B 是 MEL8，不能并入 MEL5；Cell 2023 中的抗原敲除细胞杀伤主要为 MEL8，也不能移到 MEL5。2F54 与 8WUL 的当前 ITrace 序列还有缺失残基或构建差别，先确认相应残基再计算差值。

这些系列统一预先保留完整肽段、六条 CDR（分别及合并）和两状态界面残基的并集，界面距离阈值应与现有分析保持一致；氨基酸被替换的位置另外列出。比较主链摆动幅度与侧链联合熵的变化，同时报告模型是否优于“所有位置都预测不变”，再用独立重复 MD 检查改变是否稳定。不能只挑某一条 CDR、某一个残基或只保留高相关结果。原文亲和力与功能结果提供要解释的生物学问题，模型是否有解释力仍须由这些固定区域的 MD 对比证明。

详细候选、PDB 对、控制变量、原始 DOI、适合测试的区域及已发现的对应关系问题，见同目录 `biology_first_ITrace_candidate_series.csv`。本目录不包含论文改稿，也不包含新设计变体或新启动的 MD/训练。
'''
(BASE / 'ITrace生物学优先候选说明.md').write_text(report, encoding='utf-8')
print(json.dumps({'series': len(cases), 'unique_pdbs': len(audit['unique_pdbs']), 'csv': str(BASE / 'biology_first_ITrace_candidate_series.csv')}, ensure_ascii=False))
