# 蛋白 / 蛋白复合物二分类的独立边界复核

本复核只使用官方结构实体、官方 UniProt 序列及原始论文，未读取模型预测或性能。root 已确认以下十例分类；本目录不改 membership。蛋白组包含同源多聚体以及同一成熟功能蛋白的加工/拆分链或拆域；复合物组包含成熟功能蛋白之外的独立结合伙伴。共享 accession、链数、entity 数量或标题单独出现 complex 都不是充分判据。

| PDB | 最终分组 | 身份与边界依据 |
|---|---|---|
| [1E3U](https://www.rcsb.org/structure/1E3U) | 蛋白 | 两个实体均为 OXA-10，均映射 P14489 的 21–266；长度均 246 aa，仅一位不同。实体数不等于不同蛋白身份。 同一种成熟酶的同源/近同源变体；单点差异不产生异源蛋白伙伴。 |
| [3F4Y](https://www.rcsb.org/structure/3F4Y) | 蛋白 | NHR 与 CHR 均来源于成熟 gp41；NHR 对应 P04580 545–580，gp41 成熟链为 511–855。原论文明确 3F4Y 是 NHR peptide 1 与 CHR alpha-peptide 3 的六螺旋束；peptide 3 即 T-2635，较 native CHR 约 50% 位点修改，并具有独立融合抑制剂用途。 root 于本次性能分组前裁定：同一成熟 gp41 的拆域及同源派生设计变体仍归 protein，保留人工边界说明。不能称为未经修改的天然片段。 |
| [3P9W](https://www.rcsb.org/structure/3P9W) | 蛋白复合物 | entry 标题和原论文均明确为 engineered autonomous human VH domain 与 VEGF。entity 1 为 P15692 35–138 的 VEGF；entity 2 为 123 aa VH 序列，以 EVQLVESGGGLVQPGGSLRLSCAAS 起始。entity 2 的 human VEGF 描述与标题、序列及原论文相矛盾。 以原论文身份、entry 标题和 VH 序列共同判定为抗原–抗体域复合物；不机械使用错误实体描述。 |
| [4LLO](https://www.rcsb.org/structure/4LLO) | 蛋白 | 两个实体均映射同一个 mouse EAG1/KCNH1 蛋白 Q60603：eag domain 6–136，CNBHD 552–724。 同一成熟通道蛋白的两个实验拆分结构域，按用户规则归 protein；标题中的 complex 本身不足以判定异源复合物。 |
| [5JJD](https://www.rcsb.org/structure/5JJD) | 蛋白 | PH 与 START 两实体均映射同一个 CERT/COL4A3BP 蛋白 Q9Y5P4。当前 RCSB reference alignment 为 20–122 与 385–624；depositor fragment 字段写 148–250 与 359–598，坐标采用的 isoform/编号不一致，因此不把 fragment 文本编号当作另一蛋白证据。 同一成熟神经酰胺转运蛋白的 PH/START 域互作，归 protein。域编号差异不影响身份结论。 |
| [3VV2](https://www.rcsb.org/structure/3VV2) | 蛋白复合物 | P58502 同一前体的两个不同功能产物：RCSB 将 enzyme 映射 94–422，propeptide 映射 25–93。原论文明确先发生自切，留下抑制性 enzyme:propeptide 复合物；前肽随后解离并被降解才释放活性成熟酶。3VV2 是 L69P-propeptide 与 S324A-subtilisin 的共晶。 root 裁定独立前肽伙伴属于成熟功能酶之外的 inhibitor/chaperone，不是同一活性成熟蛋白保留的拆分链，归 protein_complex。共享 precursor accession 不充分。 |
| [4ZRJ](https://www.rcsb.org/structure/4ZRJ) | 蛋白 | 两个实体均为 Merlin/NF2 P35240：FERM 1–320，CTD 506–595；论文描述 Merlin 自抑制尾部对 FERM 的调控。 同一成熟蛋白头尾的实验拆域互作，归 protein。 |
| [1RKE](https://www.rcsb.org/structure/1RKE) | 蛋白 | Vinculin head 与 tail 均映射 P18206；head reference alignment 1–258，tail 当前 reference alignment 947–1131；depositor/title 使用 879–1066 的另一个编号，编号不完全一致但蛋白身份一致。 同一成熟 vinculin 头尾实验拆域，归 protein；不因为两条链或标题用 complex 而判定为不同蛋白。 |
| [4BQD](https://www.rcsb.org/structure/4BQD) | 蛋白复合物 | 原论文 Results 的 Identification and Optimization of Peptides That Neutralize TFPI 节说明：27-mer mRNA display 文库筛出 hit，截短/替换后得到 compound 3。其化学序列为 Ac-FQSKpNVHVDGYFERL-Aib-AKL-NH2，是 20-mer TFPI 拮抗剂；RCSB 的 21 个位置包含端基/非标准残基的规范化表示。官方 Q59EE5 与 canonical TFPI P10646 均不含 FQSKPNVHVDGYFERLAAKL 这一连续序列。 其来源是独立文库筛选与优化，不是 TFPI 自身成熟蛋白的片段；归 protein_complex。判断依据不是 synthetic 标签。 |
| [4JO6](https://www.rcsb.org/structure/4JO6) | 蛋白复合物 | 主体为 streptavidin P22629 25–183；伙伴是独立 38 aa SBP-Tag。原论文及 RCSB 摘要说明该 tag 跨接两个 streptavidin 亚基的 biotin-binding pockets。 streptavidin 同源多聚体本身可归 protein，但存在独立 SBP-Tag 结合伙伴，完整体系归 protein_complex。 |

这些案例中没有身份来源仍未查清的记录。3F4Y 与 3VV2 的事实已明确；它们的分组取决于上面事先明确的成熟功能蛋白边界，root 已作出裁定。

## 原始证据定位

**1E3U**：[https://www.rcsb.org/structure/1E3U](https://www.rcsb.org/structure/1E3U)；[https://data.rcsb.org/rest/v1/core/entry/1E3U](https://data.rcsb.org/rest/v1/core/entry/1E3U)
原 GraphQL payload：`official_rcsb_payloads\batch_000.json`，SHA-256 `4316ed5146520776fad160aaf3dfa3568774abfd30e8d97aab93c52a95537eb3`；独立 entry REST 快照 `E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\grouping_protein_complex_20261024\boundary_review\1E3U_entry.json`，SHA-256 `8bbd91e79a8c48c68bd277b0b8dbd9b4c8b4a2b8fe9b956870d1012fdf85a87e`。

**3F4Y**：[https://www.rcsb.org/structure/3F4Y](https://www.rcsb.org/structure/3F4Y)；[https://data.rcsb.org/rest/v1/core/entry/3F4Y](https://data.rcsb.org/rest/v1/core/entry/3F4Y)；[https://www.pnas.org/doi/10.1073/pnas.0902663106](https://www.pnas.org/doi/10.1073/pnas.0902663106)；[https://pmc.ncbi.nlm.nih.gov/articles/PMC2736470/](https://pmc.ncbi.nlm.nih.gov/articles/PMC2736470/)；[https://rest.uniprot.org/uniprotkb/P04580.json](https://rest.uniprot.org/uniprotkb/P04580.json)
原 GraphQL payload：`official_rcsb_payloads\batch_011.json`，SHA-256 `0fd472224efc7d9e6a8fdeecfdd3ebd97349264b67a1a6eecd85b5d322c69378`；独立 entry REST 快照 `E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\grouping_protein_complex_20261024\boundary_review\3F4Y_entry.json`，SHA-256 `86b8ff534b06b8523f1c93bfd9d9e5199c57f0acd67bd5baf3f37c9251e69f32`。

**3P9W**：[https://www.rcsb.org/structure/3P9W](https://www.rcsb.org/structure/3P9W)；[https://data.rcsb.org/rest/v1/core/entry/3P9W](https://data.rcsb.org/rest/v1/core/entry/3P9W)；[https://doi.org/10.1016/j.jmb.2013.03.020](https://doi.org/10.1016/j.jmb.2013.03.020)
原 GraphQL payload：`official_rcsb_payloads\batch_012.json`，SHA-256 `d738db97c79e71e6eed28545dd528feab7287ee03530961ba5ecc6c0426ce30b`；独立 entry REST 快照 `E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\grouping_protein_complex_20261024\boundary_review\3P9W_entry.json`，SHA-256 `bd513106029c5e93b1b4aebd9acb7d261cb45c7f37e84e7b44563c541ade30b1`。

**4LLO**：[https://www.rcsb.org/structure/4LLO](https://www.rcsb.org/structure/4LLO)；[https://data.rcsb.org/rest/v1/core/entry/4LLO](https://data.rcsb.org/rest/v1/core/entry/4LLO)；[https://doi.org/10.1038/nature12487](https://doi.org/10.1038/nature12487)
原 GraphQL payload：`official_rcsb_payloads\batch_016.json`，SHA-256 `4e22dddd99c77d640dfc75ed14e421e48a0dd594248c5351c37817944296ab5d`；独立 entry REST 快照 `E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\grouping_protein_complex_20261024\boundary_review\4LLO_entry.json`，SHA-256 `c44fab599c2a9527a94244f220d54643cfc084b3dbb55da5f1337e127800d54c`。

**5JJD**：[https://www.rcsb.org/structure/5JJD](https://www.rcsb.org/structure/5JJD)；[https://data.rcsb.org/rest/v1/core/entry/5JJD](https://data.rcsb.org/rest/v1/core/entry/5JJD)；[https://doi.org/10.1074/jbc.M117.780007](https://doi.org/10.1074/jbc.M117.780007)
原 GraphQL payload：`official_rcsb_payloads\batch_019.json`，SHA-256 `8026e502e018a93f019e23e8fa77cf9d64a62fd9887174e99386dd2b228b49dc`；独立 entry REST 快照 `E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\grouping_protein_complex_20261024\boundary_review\5JJD_entry.json`，SHA-256 `680e06f89f6ebe7d4e2b74146dba5f977aa1eaa9a9c97e678e9cb67aaf193599`。

**3VV2**：[https://www.rcsb.org/structure/3VV2](https://www.rcsb.org/structure/3VV2)；[https://data.rcsb.org/rest/v1/core/entry/3VV2](https://data.rcsb.org/rest/v1/core/entry/3VV2)；[https://doi.org/10.1111/febs.12091](https://doi.org/10.1111/febs.12091)；[https://rest.uniprot.org/uniprotkb/P58502.json](https://rest.uniprot.org/uniprotkb/P58502.json)
原 GraphQL payload：`official_rcsb_payloads\batch_014.json`，SHA-256 `30f662b873ac1bc366301685e1025156249008d913f79d8ded9c07782c68fd70`；独立 entry REST 快照 `E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\grouping_protein_complex_20261024\boundary_review\3VV2_entry.json`，SHA-256 `497955b0041b0e8a5f1f70ba2267c0f3b6fcd92749475f8d34e766b4e430535a`。

**4ZRJ**：[https://www.rcsb.org/structure/4ZRJ](https://www.rcsb.org/structure/4ZRJ)；[https://data.rcsb.org/rest/v1/core/entry/4ZRJ](https://data.rcsb.org/rest/v1/core/entry/4ZRJ)；[https://doi.org/10.1038/cr.2015.69](https://doi.org/10.1038/cr.2015.69)
原 GraphQL payload：`official_rcsb_payloads\batch_018.json`，SHA-256 `eef5a3f95ec868db589cd71ab7e5983adfc92e966dad120aae0e3c66f44f8ceb`；独立 entry REST 快照 `E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\grouping_protein_complex_20261024\boundary_review\4ZRJ_entry.json`，SHA-256 `5a5b332a8254f3c3b85caf4a0a8c787644d20004df1569e273a095ca525748be`。

**1RKE**：[https://www.rcsb.org/structure/1RKE](https://www.rcsb.org/structure/1RKE)；[https://data.rcsb.org/rest/v1/core/entry/1RKE](https://data.rcsb.org/rest/v1/core/entry/1RKE)；[https://doi.org/10.1038/nature02281](https://doi.org/10.1038/nature02281)
原 GraphQL payload：`official_rcsb_payloads\batch_004.json`，SHA-256 `01e621d94012ff052c36c8100b89de53df8ba571d226124104858ff86befd70b`；独立 entry REST 快照 `E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\grouping_protein_complex_20261024\boundary_review\1RKE_entry.json`，SHA-256 `ad012f00a1e585d3723697032a527e357319be2e327fa07ec28b95742fc63ad5`。

**4BQD**：[https://www.rcsb.org/structure/4BQD](https://www.rcsb.org/structure/4BQD)；[https://data.rcsb.org/rest/v1/core/entry/4BQD](https://data.rcsb.org/rest/v1/core/entry/4BQD)；[https://pmc.ncbi.nlm.nih.gov/articles/PMC3894350/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3894350/)；[https://doi.org/10.1074/jbc.M113.533836](https://doi.org/10.1074/jbc.M113.533836)；[https://rest.uniprot.org/uniprotkb/Q59EE5.json](https://rest.uniprot.org/uniprotkb/Q59EE5.json)；[https://rest.uniprot.org/uniprotkb/P10646.json](https://rest.uniprot.org/uniprotkb/P10646.json)
原 GraphQL payload：`official_rcsb_payloads\batch_014.json`，SHA-256 `30f662b873ac1bc366301685e1025156249008d913f79d8ded9c07782c68fd70`；独立 entry REST 快照 `E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\grouping_protein_complex_20261024\boundary_review\4BQD_entry.json`，SHA-256 `928b93923a3d20f009db0e5027f97b193c3f4da49195b8679d7e5e32b2018420`。

**4JO6**：[https://www.rcsb.org/structure/4JO6](https://www.rcsb.org/structure/4JO6)；[https://data.rcsb.org/rest/v1/core/entry/4JO6](https://data.rcsb.org/rest/v1/core/entry/4JO6)；[https://doi.org/10.1107/S0907444913002576](https://doi.org/10.1107/S0907444913002576)
原 GraphQL payload：`official_rcsb_payloads\batch_016.json`，SHA-256 `4e22dddd99c77d640dfc75ed14e421e48a0dd594248c5351c37817944296ab5d`；独立 entry REST 快照 `E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\grouping_protein_complex_20261024\boundary_review\4JO6_entry.json`，SHA-256 `7d4a8df6427f4833a8d740d50bd0c953c03545d2adeee4ba81d1ed1b3067b65b`。

## 必须保留的限定

- P58502 当前 UniProt 自动预测 Propeptide 25–106 / Chain 107–422 与该实验原论文及 RCSB 25–93 / 94–422 不一致。原文以去信号肽构建体编号 Gly1–Leu69 / Gly70–Gly398。该处实验功能身份以原论文为据，不把当前预测 cleavage site 当作已验证实验边界。
- 3P9W entity 2 的描述 human VEGF 为误名，不能据此把两实体并为一个 VEGF 蛋白。
- 3F4Y 属人为定义的边界：同一成熟 gp41 派生的设计 CHR 抑制肽仍按源蛋白/拆域规则归 protein。应保留该定义；若后续做边界敏感性分析，应一次性预先指定移动该整个案例，不能按模型成绩挑选。当前未读取任何性能、未运行任何预测或性能敏感性计算。
- 4BQD 化学 compound 3 是 20-mer；RCSB entity canonical sequence 为 21 个符号，并含 X 与非标准化学残基。不能直接用规范化序列长度断言为天然 21 aa 蛋白。

下载记录见 fetch_manifest.json。Jina Reader 连接失败后，采用官方 REST 以及 web 工具的原论文读取作为回退。4BQD 原文完整 HTML 已保存并核对有 27-mer mRNA display 与 compound 3 段落。3F4Y 的 PNAS 正文通过 web 工具已读到，直接 HTTP 下载被 403 阻止；PMC 的 HTTP 200 返回的是挑战页面，已标为非文章，不冒充正文快照。
