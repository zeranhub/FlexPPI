from pathlib import Path
import csv
import json
import math
import hashlib

ROOT = Path(r'E:\Project\FlexPPI\experiments\_sources')
OLD = ROOT / 'tcr_expanded_screen_v41_20261017'
OUT = ROOT / 'tcr_joint_case_search_20261019' / 'residue_story' / 'independent_review'
OUT.mkdir(parents=True, exist_ok=True)
SOURCE = OLD / 'mel5_validation' / 'corrected_pair_all_aligned_positions.csv'
with SOURCE.open(encoding='utf-8-sig', newline='') as f:
    rows = list(csv.DictReader(f))
with (OLD / 'mel5_context' / 'all_natural_pMHC_context_comparison.csv').open(encoding='utf-8-sig', newline='') as f:
    context = list(csv.DictReader(f))

def atoms(path):
    result = []
    seen = set()
    for line in path.read_text().splitlines():
        if not line.startswith('ATOM  ') or line[16] not in (' ', 'A'):
            continue
        name = line[12:16].strip()
        if name.startswith('H') or (len(name) > 1 and name[0].isdigit() and name[1] == 'H'):
            continue
        key = (line[21], int(line[22:26]), line[26], name)
        if key in seen:
            continue
        seen.add(key)
        result.append({'chain': line[21], 'id': int(line[22:26]), 'res': line[17:20], 'atom': name,
                       'xyz': tuple(float(line[i:i+8]) for i in (30, 38, 46))})
    return result

structures = {p: atoms(OLD / 'mel5_validation' / f'{p}_model_input.pdb') for p in ('3hg1', '4jff')}

def distance_record(pdb, mhc_id, chain, residue, sidechain_only=False):
    aa = [a for a in structures[pdb] if a['chain'] == 'A' and a['id'] == mhc_id]
    if sidechain_only:
        aa = [a for a in aa if a['atom'] not in ('N', 'CA', 'C', 'O', 'OXT')]
    bb = [a for a in structures[pdb] if a['chain'] == chain and a['id'] == residue]
    distances = [(math.dist(a['xyz'], b['xyz']), a['atom'], b['atom']) for a in aa for b in bb]
    best = min(distances)
    return {'minimum_distance_A': best[0], 'MHC_atom': best[1], 'TCR_atom': best[2],
            'atom_pairs_within_4A': sum(d[0] <= 4 for d in distances),
            'MHC_atom_scope': 'side chain only' if sidechain_only else 'all heavy atoms'}

def summary(row):
    fields = ['role', 'region', 'id_a', 'id_b', 'aa_a', 'aa_b', 'same_aa']
    result = {k: row[k] for k in fields}
    for m in ('rmsf', 'hchi'):
        result[m+'_valid'] = row[m+'_valid'] == 'True'
        for source in ('full', 'md'):
            for state in ('a', 'b'):
                val = row[f'{m}_{source}_{state}']
                result[f'{m}_{source}_{state}'] = float(val) if val else None
    return result

beta = []
for row in rows:
    if row['region'] != 'CDR2β':
        continue
    z = summary(row)
    z['publication_beta_position'] = int(row['id_a'])
    z['publication_numbering_basis'] = ('53 is explicit in the primary article; contiguous 50–54 labels are projected from the five-position sequence alignment, not relabelled input IDs.')
    for state, pdb in (('a', '3hg1'), ('b', '4jff')):
        z['Q72_contacts_'+state] = distance_record(pdb, 72, 'E', int(row['id_'+state]))
        z['Q72_sidechain_contacts_'+state] = distance_record(pdb, 72, 'E', int(row['id_'+state]), True)
    beta.append(z)

R65_partners = []
for row in rows:
    if row['role'] not in ('tcr_alpha', 'tcr_beta'):
        continue
    ds = [distance_record(p, 65, row['chain_'+s], int(row['id_'+s])) for s,p in (('a','3hg1'),('b','4jff'))]
    if not any(d['minimum_distance_A'] <= 4 for d in ds):
        continue
    z = summary(row)
    z['R65_contacts_a'], z['R65_contacts_b'] = ds
    R65_partners.append(z)

mhc_keys = []
for position in (65,68,72,75,76):
    row = next(r for r in context if r['role'] == 'mhc' and r['id_a'] == str(position))
    z = {k: row[k] for k in ('id_a','aa_a')}
    for k in ('hchi_full_a','hchi_full_b','hchi_full_delta','hchi_pMHC_only_delta','hchi_md_a','hchi_md_b'):
        z[k] = float(row[k])
    mhc_keys.append(z)

assessment = {
    'status':'complete', 'method':'Read-only reanalysis of corrected residue-pair tables and existing model-input snapshots; no new prediction, training, or MD.',
    'source':str(SOURCE), 'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'geometry_scope':'Simulation initial snapshots used as model input. A 4 Å heavy-atom cutoff defines proximity, not occupancy, binding energy, nor a chemically assigned hydrogen bond.',
    'state_a':'3HG1 / MEL5', 'state_b':'4JFF / alpha24beta17',
    'all_five_beta_CDR2_positions':beta, 'R65_all_TCR_contact_partners_union_4A':R65_partners,
    'all_MHC_CDR2beta_contact_positions_plus_R65':mhc_keys,
    'preferred_pair':{'Hchi_site':'HLA Q72', 'RMSF_site':'TCR beta Gly54 in publication-compatible numbering; input E54 in 3HG1 and E55 in 4JFF',
        'Hchi_prediction_nats':[1.792449712753296,1.3750494718551636],
        'RMSF_prediction_A':[1.1178646087646484,0.7468647956848145],
        'mechanistic_bridge':'Q72 side-chain proximity expands from beta Val51 in MEL5 to the beta Gly51–Pro52–Phe53–Gly54 patch in alpha24beta17. Q72 sidechain state restriction and lower motion of an unchanged adjacent Gly54 are consistent consequences of the redesigned local interface.',
        'causal_status':'Compatible local coupling hypothesis, not demonstrated Q72 Hchi -> Gly54 RMSF causation.',
        'sequence_confounder':'Gly52 -> Pro52 changes intrinsic loop-backbone preferences; Ile53 -> Phe53 and the other CDR2 substitutions also contribute.'},
    'R65_alternative':'R65 Hchi decreases, but its unchanged alpha Gly aligned WT94/engineered95 exits the 4 Å contact set in the engineered input (3.567 to 4.957 Å), and model RMSF is almost unchanged (0.856 to 0.851 Å). Its direct engineered CDR3alpha partners are mutated residues, so this is weaker than Q72–Gly54 as a two-output residue story.',
    'complete_patch_check':'All five CDR2beta positions show lower model RMSF. All four contacting MHC sidechains are retained: K68 and Q72 Hchi decrease; R75 is almost unchanged in prediction and increases in MD; V76 slightly increases. Do not turn extra contacts into a universal entropy-decrease claim.',
    'recommended_existing_trajectory_analysis':['For all three replicas in each state, report Q72 chi1–chi3 joint-state occupancy and lifetimes.', 'Within each replica, measure Q72-sidechain contact occupancy to beta Phe53 and unchanged Gly54, rather than relying on a single snapshot.', 'Compare Gly54 backbone fluctuations in fixed-duration trajectory windows classified by Q72 state and by contact formation. Align against a consistent rigid reference; windows require sufficient state counts and repeated support across replicas.', 'Measure time-lagged relationships and residence times as descriptive coupling support, accounting for autocorrelation; no directional causality claim follows merely from a correlation.'],
    'primary_article':'https://pmc.ncbi.nlm.nih.gov/articles/PMC3696650/',
    'article_fact':'SPR measures affinity 18 microM to 0.6 nM and slower dissociation; crystallography assigns new MHC contacts, including Ile53 to Phe53. The article does not perform MD. Its total binding entropy and solvent contribution are not bound-state Hchi.'
}
(OUT/'paired_residue_assessment.json').write_text(json.dumps(assessment,ensure_ascii=False,indent=2),encoding='utf8')

lines = ['# MEL5 双输出的逐残基故事：独立复核','',
'最值得展开的是 HLA Q72 与 TCR β链 Gly54 所在的 CDR2β 接触片区。Gly54 在两种受体中没有发生替换，模型 RMSF 从 1.118 降至 0.747 Å；其对侧 Q72 的模型 Hχ 从 1.792 降至 1.375 nats。更准确的机制是亲和成熟改变了两者之间的接触网络，Q72 取向状态和 Gly54 主链运动共同变得受限。不能把 Hχ 本身当作作用力，写成“Hχ 下降导致 RMSF 下降”。', '',
'## 全部五个 CDR2β 位置（不按预测效果筛点）','',
'原文明确使用 β53 Ile→Phe；下表其余出版兼容标号按相邻连续配准投影。模型输入中工程受体因 N 端多一个残基，对应编号整体加 1。正文和图中若采用原文编号，Gly54 对应两个输入 E54 与 E55，不应把后者误认成 Gly55。','',
'| 原文兼容位置 | 输入 WT→工程 | 氨基酸 | 模型 RMSF / Å | MD RMSF / Å | Q72 最小重原子距离 WT→工程 / Å |','|---|---|---|---|---|---|']
for z in beta:
    lines.append(f"| β{z['publication_beta_position']} | E{z['id_a']}→E{z['id_b']} | {z['aa_a']}→{z['aa_b']} | {z['rmsf_full_a']:.3f}→{z['rmsf_full_b']:.3f} | {z['rmsf_md_a']:.3f}→{z['rmsf_md_b']:.3f} | {z['Q72_contacts_a']['minimum_distance_A']:.3f}→{z['Q72_contacts_b']['minimum_distance_A']:.3f} |")
lines += ['', 'Q72 在 WT 输入中以 ≤4 Å 接近 βVal51；工程输入中接近 βGly51、Pro52、Phe53、Gly54。工程输入 Q72 OE1—βPhe53 N 的距离为 2.943 Å。尤其 Q72 侧链到未替换 Gly54 主链的最近距离从 WT 的 8.265 Å 变为工程态 Q72 NE2—Gly54 O 的 3.120 Å。原始输入是模拟起点，距离只用于指明空间邻接，不替代晶体结构化学键判读或轨迹接触占据率。', '',
'## 完整对侧片区与第二候选','',
'四个与 CDR2β 相接近的 HLA 侧链是 K68、Q72、R75、V76。模型 Hχ 依次为 3.017→2.769、1.792→1.375、2.489→2.444、0.558→0.633 nats。这说明局部取向重分配有方向差异，不能说所有新增接触都使 Hχ 降低。Q72 是空间接触和双输出均较清楚的局部放大例子，应和完整五位 CDR2β、完整四位 MHC 片区共同呈现。', '',
'R65 的模型 Hχ 为 1.794→1.412 nats，但其工程输入直接接触的是 αAsp93、Gly94、Arg96 和 βIle57（此处均为输入编号）。其中 αGly94 来自 WT 的 Ala93 替换，不能当作未突变 Gly 对照；真正保持 Gly 的 WT α94→工程 α95 从 WT 的 3.567 Å 接近变为工程态 4.957 Å，退出 4 Å 接触集合，其预测 RMSF 仅 0.856→0.851 Å。R65 可补充接触网络重排，作为“Hχ 变化对应另一个未变残基 RMSF”的故事弱于 Q72—βGly54。', '',
'## 可用的准确中文叙述','',
'MEL5 的亲和成熟使 TCR β链第二识别环与 HLA 的接触片区发生重排。在原研究重点分析的 βIle53→Phe53 位点附近，HLA Q72 与第二识别环形成更广泛的空间接近；在用于预测的结构中，Q72 侧链与相邻且未替换的 βGly54 主链由相距 8.27 Å 变为 3.12 Å 的近接触。模型同时预测 Q72 的侧链联合状态熵由 1.792 降至 1.375 nats，βGly54 的主链 RMSF 由 1.118 降至 0.747 Å。两项输出分别把这一片区的侧链取向集中与识别环摆动收紧定位到具体位置，使新增的 TCR–MHC 接触获得了局部动态解释。原研究通过 SPR 测得亲和力增强并主要归因于解离减慢，通过晶体结构观察到新增 MHC 接触；模型预测进一步指出，这种界面重排伴随 Q72 取向状态与 βGly54 骨架运动的共同约束。', '',
'这段叙述的因果起点是已发表的受体替换及接触网络变化，而非“熵数值变小”。同一片区还存在 Gly52→Pro52 等替换，可直接改变识别环骨架偏好；因此目前证据支持局部共同响应，不能单独归因于 Q72 对 Gly54 的单向作用。', '',
'## 最有价值的补强','',
'直接复用六条已有 MD 轨迹，按重复分别给出 Q72 χ1–χ3 状态占据、Q72—βPhe53/βGly54 接触占据与持续时间，再在固定时长窗口中检验：Q72 处于主要取向状态且保持对应接触时，βGly54 的骨架波动是否更小。若这一关系在各条重复中出现，才能把两个预测读数之间的“空间相邻”进一步推进为“同一接触状态中的动态关联”。这是有价值的证据，不需要重新生产 MD；即使结果成立也不应宣称单向因果。', '',
'原研究：[T-cell Receptor Specificity Maintained by Altered Thermodynamics](https://pmc.ncbi.nlm.nih.gov/articles/PMC3696650/)。其方法包括 SPR、晶体结构、肽替换与温度依赖热力学分析，未报告 MD；总结合熵中的溶剂贡献不能用结合态 Hχ 取代。', '',
'本次仅只读复核，没有修改论文，没有新增训练、推理或 MD。完整数值、所有直接接触伙伴和方法边界保存在同名 JSON。']
(OUT/'paired_residue_assessment.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
print(json.dumps({'beta_CDR2':beta,'R65_partners':R65_partners,'outputs':[str(OUT/'paired_residue_assessment.md'),str(OUT/'paired_residue_assessment.json')]},ensure_ascii=False,indent=2))
