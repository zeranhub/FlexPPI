"""Independent biological boundary audit. Does not open predictions or membership."""
from pathlib import Path
import csv
import datetime
import hashlib
import json

HERE = Path(__file__).resolve().parent
GROUPING = HERE.parent
ENTITIES = json.loads((GROUPING / 'official_entities_normalized.json').read_text(encoding='utf-8'))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

cases = [
    ('1E3U', 'protein', 'same_mature_protein_variant',
     '两个实体均为 OXA-10，均映射 P14489 的 21–266；长度均 246 aa，仅一位不同。实体数不等于不同蛋白身份。',
     '同一种成熟酶的同源/近同源变体；单点差异不产生异源蛋白伙伴。', []),
    ('3F4Y', 'protein', 'same_mature_protein_derived_domain_design',
     'NHR 与 CHR 均来源于成熟 gp41；NHR 对应 P04580 545–580，gp41 成熟链为 511–855。原论文明确 3F4Y 是 NHR peptide 1 与 CHR alpha-peptide 3 的六螺旋束；peptide 3 即 T-2635，较 native CHR 约 50% 位点修改，并具有独立融合抑制剂用途。',
     'root 于本次性能分组前裁定：同一成熟 gp41 的拆域及同源派生设计变体仍归 protein，保留人工边界说明。不能称为未经修改的天然片段。',
     ['https://www.pnas.org/doi/10.1073/pnas.0902663106', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC2736470/', 'https://rest.uniprot.org/uniprotkb/P04580.json']),
    ('3P9W', 'protein_complex', 'antigen_engineered_antibody_domain',
     'entry 标题和原论文均明确为 engineered autonomous human VH domain 与 VEGF。entity 1 为 P15692 35–138 的 VEGF；entity 2 为 123 aa VH 序列，以 EVQLVESGGGLVQPGGSLRLSCAAS 起始。entity 2 的 human VEGF 描述与标题、序列及原论文相矛盾。',
     '以原论文身份、entry 标题和 VH 序列共同判定为抗原–抗体域复合物；不机械使用错误实体描述。',
     ['https://doi.org/10.1016/j.jmb.2013.03.020']),
    ('4LLO', 'protein', 'same_mature_protein_split_domains',
     '两个实体均映射同一个 mouse EAG1/KCNH1 蛋白 Q60603：eag domain 6–136，CNBHD 552–724。',
     '同一成熟通道蛋白的两个实验拆分结构域，按用户规则归 protein；标题中的 complex 本身不足以判定异源复合物。',
     ['https://doi.org/10.1038/nature12487']),
    ('5JJD', 'protein', 'same_mature_protein_split_domains',
     'PH 与 START 两实体均映射同一个 CERT/COL4A3BP 蛋白 Q9Y5P4。当前 RCSB reference alignment 为 20–122 与 385–624；depositor fragment 字段写 148–250 与 359–598，坐标采用的 isoform/编号不一致，因此不把 fragment 文本编号当作另一蛋白证据。',
     '同一成熟神经酰胺转运蛋白的 PH/START 域互作，归 protein。域编号差异不影响身份结论。',
     ['https://doi.org/10.1074/jbc.M117.780007']),
    ('3VV2', 'protein_complex', 'mature_enzyme_separate_inhibitory_propeptide',
     'P58502 同一前体的两个不同功能产物：RCSB 将 enzyme 映射 94–422，propeptide 映射 25–93。原论文明确先发生自切，留下抑制性 enzyme:propeptide 复合物；前肽随后解离并被降解才释放活性成熟酶。3VV2 是 L69P-propeptide 与 S324A-subtilisin 的共晶。',
     'root 裁定独立前肽伙伴属于成熟功能酶之外的 inhibitor/chaperone，不是同一活性成熟蛋白保留的拆分链，归 protein_complex。共享 precursor accession 不充分。',
     ['https://doi.org/10.1111/febs.12091', 'https://rest.uniprot.org/uniprotkb/P58502.json']),
    ('4ZRJ', 'protein', 'same_mature_protein_split_domains',
     '两个实体均为 Merlin/NF2 P35240：FERM 1–320，CTD 506–595；论文描述 Merlin 自抑制尾部对 FERM 的调控。',
     '同一成熟蛋白头尾的实验拆域互作，归 protein。',
     ['https://doi.org/10.1038/cr.2015.69']),
    ('1RKE', 'protein', 'same_mature_protein_split_domains',
     'Vinculin head 与 tail 均映射 P18206；head reference alignment 1–258，tail 当前 reference alignment 947–1131；depositor/title 使用 879–1066 的另一个编号，编号不完全一致但蛋白身份一致。',
     '同一成熟 vinculin 头尾实验拆域，归 protein；不因为两条链或标题用 complex 而判定为不同蛋白。',
     ['https://doi.org/10.1038/nature02281']),
    ('4BQD', 'protein_complex', 'de_novo_selected_peptide_antagonist',
     '原论文 Results 的 Identification and Optimization of Peptides That Neutralize TFPI 节说明：27-mer mRNA display 文库筛出 hit，截短/替换后得到 compound 3。其化学序列为 Ac-FQSKpNVHVDGYFERL-Aib-AKL-NH2，是 20-mer TFPI 拮抗剂；RCSB 的 21 个位置包含端基/非标准残基的规范化表示。官方 Q59EE5 与 canonical TFPI P10646 均不含 FQSKPNVHVDGYFERLAAKL 这一连续序列。',
     '其来源是独立文库筛选与优化，不是 TFPI 自身成熟蛋白的片段；归 protein_complex。判断依据不是 synthetic 标签。',
     ['https://pmc.ncbi.nlm.nih.gov/articles/PMC3894350/', 'https://doi.org/10.1074/jbc.M113.533836', 'https://rest.uniprot.org/uniprotkb/Q59EE5.json', 'https://rest.uniprot.org/uniprotkb/P10646.json']),
    ('4JO6', 'protein_complex', 'streptavidin_independent_binding_tag',
     '主体为 streptavidin P22629 25–183；伙伴是独立 38 aa SBP-Tag。原论文及 RCSB 摘要说明该 tag 跨接两个 streptavidin 亚基的 biotin-binding pockets。',
     'streptavidin 同源多聚体本身可归 protein，但存在独立 SBP-Tag 结合伙伴，完整体系归 protein_complex。',
     ['https://doi.org/10.1107/S0907444913002576']),
]

rows = []
for pdb, decision, rule, facts, inference, extra in cases:
    ent = ENTITIES[pdb]
    payloads = {e['payload_path']: e['payload_sha256'] for e in ent}
    for path, expected in payloads.items():
        assert sha(GROUPING / path) == expected, (pdb, path)
    entry_path = HERE / f'{pdb}_entry.json'
    entry = json.loads(entry_path.read_text(encoding='utf-8'))
    sources = [f'https://www.rcsb.org/structure/{pdb}', f'https://data.rcsb.org/rest/v1/core/entry/{pdb}'] + extra
    rows.append({
        'pdb_id': pdb, 'decision_approved_by_root': decision, 'rule': rule,
        'facts_zh': facts, 'classification_reason_zh': inference,
        'factual_identity_resolved': True,
        'remaining_unresolved_membership': False,
        'operational_boundary_note': pdb in {'3F4Y', '3VV2'},
        'entry_title': entry['struct']['title'],
        'rcsb_graphql_payloads_sha256': payloads,
        'rcsb_entry_snapshot_path': str(entry_path),
        'rcsb_entry_snapshot_sha256': sha(entry_path),
        'sources': sources,
        'primary_citation': entry.get('rcsb_primary_citation'),
        'entities': [{k: v for k, v in e.items() if k != 'sequence'} for e in ent],
    })

warnings = [
    'P58502 当前 UniProt 自动预测 Propeptide 25–106 / Chain 107–422 与该实验原论文及 RCSB 25–93 / 94–422 不一致。原文以去信号肽构建体编号 Gly1–Leu69 / Gly70–Gly398。该处实验功能身份以原论文为据，不把当前预测 cleavage site 当作已验证实验边界。',
    '3P9W entity 2 的描述 human VEGF 为误名，不能据此把两实体并为一个 VEGF 蛋白。',
    '3F4Y 属人为定义的边界：同一成熟 gp41 派生的设计 CHR 抑制肽仍按源蛋白/拆域规则归 protein。应保留该定义；若后续做边界敏感性分析，应一次性预先指定移动该整个案例，不能按模型成绩挑选。当前未读取任何性能、未运行任何预测或性能敏感性计算。',
    '4BQD 化学 compound 3 是 20-mer；RCSB entity canonical sequence 为 21 个符号，并含 X 与非标准化学残基。不能直接用规范化序列长度断言为天然 21 aa 蛋白。',
]
result = {
    'generated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'scope': 'independent read-only biological boundary audit; no performance or model output read',
    'membership_files_modified': False, 'model_inference_run': False,
    'approved_rule': 'protein includes homooligomers and processed/split chains/domains of one mature functional protein, including source-derived domain design variants; independently functioning partner products/selected binders are protein_complex',
    'normalized_entities_sha256': sha(GROUPING / 'official_entities_normalized.json'),
    'warnings': warnings, 'records': rows,
}
(HERE / 'review.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
fields = ['pdb_id', 'decision_approved_by_root', 'rule', 'facts_zh', 'classification_reason_zh', 'factual_identity_resolved', 'remaining_unresolved_membership', 'operational_boundary_note', 'rcsb_entry_snapshot_path', 'rcsb_entry_snapshot_sha256', 'sources']
with (HERE / 'review.csv').open('w', encoding='utf-8-sig', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for r in rows:
        row = {k: r[k] for k in fields}
        row['sources'] = ' | '.join(row['sources'])
        writer.writerow(row)

md = ['# 蛋白 / 蛋白复合物二分类的独立边界复核', '',
      '本复核只使用官方结构实体、官方 UniProt 序列及原始论文，未读取模型预测或性能。root 已确认以下十例分类；本目录不改 membership。蛋白组包含同源多聚体以及同一成熟功能蛋白的加工/拆分链或拆域；复合物组包含成熟功能蛋白之外的独立结合伙伴。共享 accession、链数、entity 数量或标题单独出现 complex 都不是充分判据。', '',
      '| PDB | 最终分组 | 身份与边界依据 |', '|---|---|---|']
for r in rows:
    group = '蛋白' if r['decision_approved_by_root'] == 'protein' else '蛋白复合物'
    md.append(f"| [{r['pdb_id']}](https://www.rcsb.org/structure/{r['pdb_id']}) | {group} | {r['facts_zh']} {r['classification_reason_zh']} |")
md.extend(['', '这些案例中没有身份来源仍未查清的记录。3F4Y 与 3VV2 的事实已明确；它们的分组取决于上面事先明确的成熟功能蛋白边界，root 已作出裁定。', '', '## 原始证据定位', ''])
for r in rows:
    md.append(f"**{r['pdb_id']}**：" + '；'.join(f'[{u}]({u})' for u in r['sources']))
    for path, digest in r['rcsb_graphql_payloads_sha256'].items():
        md.append(f"原 GraphQL payload：`{path}`，SHA-256 `{digest}`；独立 entry REST 快照 `{r['rcsb_entry_snapshot_path']}`，SHA-256 `{r['rcsb_entry_snapshot_sha256']}`。")
    md.append('')
md.extend(['## 必须保留的限定', ''] + [f'- {w}' for w in warnings])
md.extend(['', '下载记录见 fetch_manifest.json。Jina Reader 连接失败后，采用官方 REST 以及 web 工具的原论文读取作为回退。4BQD 原文完整 HTML 已保存并核对有 27-mer mRNA display 与 compound 3 段落。3F4Y 的 PNAS 正文通过 web 工具已读到，直接 HTTP 下载被 403 阻止；PMC 的 HTTP 200 返回的是挑战页面，已标为非文章，不冒充正文快照。'])
(HERE / 'review_zh.md').write_text('\n'.join(md) + '\n', encoding='utf-8')

# Mark the HTTP-200 challenge honestly, without deleting the evidence of failure.
manifest_path = HERE / 'fetch_manifest.json'
manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
for item in manifest['files']:
    if item['file'] == '3F4Y_primary_pmc.html':
        item['content_validation'] = 'blocked_challenge_page_not_article; not used as primary text evidence'
    elif item['file'] == '4BQD_primary_fulltext.html':
        item['content_validation'] = 'article text contains Identification and Optimization, 27-mer mRNA display, and compound 3 sequence'
manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')

out_manifest = []
for name in ['build_review.py', 'review.json', 'review.csv', 'review_zh.md', 'fetch_manifest.json']:
    out_manifest.append({'path': str(HERE / name), 'sha256': sha(HERE / name)})
(HERE / 'delivery_hashes.json').write_text(json.dumps(out_manifest, indent=2), encoding='utf-8')
assert len(rows) == 10
assert sum(r['decision_approved_by_root'] == 'protein' for r in rows) == 6
print(json.dumps({'records': len(rows), 'protein': 6, 'protein_complex': 4, 'files': out_manifest}, ensure_ascii=True))
