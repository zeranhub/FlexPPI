# 局部蛋白结构图的本地文献参照

本次按用户提供的 E:/paper/papers_MinerU 检索并实际查看了以下论文原图。这里只借鉴表示层次、连接显示及整体与局部之间的关系；本文全部蛋白原子、MD 代表帧、指标和 PyMOL 渲染仍来自 FlexPPI 自身案例数据。没有将参考文献中的结构图片放入论文。

1. Xu 等，2021，OPUS-Rota4，Figure 4。原文 [Xu_2021_bbab529.md](E:/paper/papers_MinerU/Briefings_in_Bioinformatics/2021/Xu_2021_bbab529.md:153)，原图 [Figure 4](E:/paper/papers_MinerU/Briefings_in_Bioinformatics/2021/Xu_2021_bbab529.assets/images/dca1631a3835348185aefd4850ee2ee5e228fe8f1abf1d84102819906880edcc.jpg)。观察：骨架保留完整 cartoon，侧链以更细的原子键层叠加，结构对照不丢失各自骨架。应用：KIX 保留三条完整链，L12 用真实 N–CA–C 主链细棒承托侧链，八残基的七条肽键与两处脯氨酸环完整显示。

2. Liang 等，2024，IsAb2.0，Figure 3。原文 [Liang_2024_bbae445.md](E:/paper/papers_MinerU/Briefings_in_Bioinformatics/2024/Liang_2024_bbae445.md:106)，原图 [Figure 3](E:/paper/papers_MinerU/Briefings_in_Bioinformatics/2024/Liang_2024_bbae445.assets/images/235487e2b07fdc5b17bd9cc3fbde65e74c91d9494e14d81371d9d09226081b61.jpg)。观察：浅色 cartoon 提供局部背景，重点残基棒状结构能够追溯到各自骨架，残基标签与颜色支持定位。应用：泛素 Q41/Q62 灰色邻居保留完整重原子键，真实 MD 侧链加粗、自己的 N/CA/C/O 细棒与整条泛素 Cα 迹线显示其真实归属。本文没有根据参考图加入未经计算的氢键虚线。

3. Sequeiros-Borja 等，2020，Recent advances in user-friendly computational tools to engineer protein function，Figure 3A。原文 [Sequeiros-Borja_2020_bbaa150.md](E:/paper/papers_MinerU/Briefings_in_Bioinformatics/2020/Sequeiros-Borja_2020_bbaa150.md:67)，原图 [Figure 3](E:/paper/papers_MinerU/Briefings_in_Bioinformatics/2020/Sequeiros-Borja_2020_bbaa150.assets/images/383c1b61dbb1d83a393eee63d756cdac16047f8a601835bd94909cec53b9d600.jpg)。观察：整体复合物给出位置，局部放大保留该位置的真实连接与附近二级结构。应用：每个局部图先在 PyMOL 中表示完整复合物、渲染同朝向的全景原图，再按记录像素框直接裁切；PSE 保存全景与局部相机，两种场景均保留完整蛋白。

4. Rodrigues 等，2022，Figure 2，作为补充比较。原文 [Rodrigues_2022_bbac165.md](E:/paper/papers_MinerU/Briefings_in_Bioinformatics/2022/Rodrigues_2022_bbac165.md:54)，原图 [Figure 2](E:/paper/papers_MinerU/Briefings_in_Bioinformatics/2022/Rodrigues_2022_bbac165.assets/images/dccff9d6213d2e95fba3e09be49c2438bb429dc71425693b39f0bbd3287a1c2d.jpg)。这张图采用蛋白表面与界面残基强调凹腔形状。当前问题是侧链连接和 MD 朝向，因此以 cartoon、真实键及细骨架迹线表达，不额外叠加会遮挡连接的表面层。

## 当前交付画法与检查

主背景分别保留 KIX 三元复合物的 131 个 Cα cartoon 与泛素–UBA 的 119 个 Cα cartoon。KIX L12 的完整 56 个重原子、57 条真实键为一个连通片段。泛素邻近残基完整显示 N/CA/C/O 与侧链；每个实际 MD 谷氨酰胺的 9 个重原子、8 条真实键为一个连通片段，并叠加该帧自己的整条泛素 72 个 Cα 细迹线。彩色根部小球标出该帧的真实 Cα，侧链和骨架均不平移或嫁接到静态输入上。cartoon 关闭环区平滑和片层压平，所有原子坐标保持原值。连接拓扑、源坐标、全景可见范围和精确像素裁切由 PyMOL 和 Pillow 实测核对。
