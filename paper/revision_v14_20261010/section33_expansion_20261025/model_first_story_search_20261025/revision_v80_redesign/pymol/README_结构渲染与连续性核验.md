# PyMOL结构素材

这些结构全部直接取自已用于模型推理及MD标签匹配的输入PDB，不是生成图，也不是从其他复合物替换的示意结构。CaM沿用此前核实的成熟编号4–148，只改变链号及残基标识，不改变输入原子坐标；PyMOL中的拟合仅对显示对象作刚体旋转和平移。没有添加钙离子、补全侧链或虚构缺失原子。图中的结构快照只负责定位侧链和结构环境，不代表MD分布或运动幅度。

当前正文配色经根任务从原DOCX像素核实：CaM N叶蓝色`#2166AC`、C叶浅蓝`#85ABD1`、连接区灰色`#A6AEB4`、NHE1或IQCG伙伴暖赭`#B87936`、M72/M109重点侧链青绿`#3A8D8C`。硫原子使用暗金色，避免亮黄色抢占视觉。ribbon颜色独立于stick颜色，因此强调侧链不会在螺旋上制造一条绿色色带。

优先使用`whole_NHE1_trim.png`与`whole_IQCG_trim.png`展示全复合物。两态以CaM N叶4–75的CA原子拟合，保持同一相机和同一显示比例，紧裁也使用同一裁切框；因此IQCG状态在画面中更紧凑是实际结构几何差异。N叶全范围的拟合RMSD为5.401 Å，不能将两态称为“完全相同的N叶构象”。`render_manifest.json`记录相机、颜色、原文件尺寸、裁切框和精确像素标注坐标，紧裁后的坐标在`assets.<素材名>.paired_trim.landmarks`中。

正文小面板优先使用`detail_focus_M72_NHE1_trim.png`、`detail_focus_M72_IQCG_trim.png`及相应的M109两图。这组按约15 Å局部视野取景，保留连续ribbon，侧链以更深的青绿`#276F6F`、不透明0.24 Å半径sticks突出；CaM背景cartoon透明度0.25、伙伴0.42，邻近侧链灰蓝或浅赭弱化，保证25 mm左右的最终面板仍可读。`detail_M72_*`/`detail_M109_*`为更宽的完整叶片背景备选，颜色较淡。旧图曾以9 Å空间球体选择cartoon，这会保留多个不连续短片段，形成破碎印象。本次局部图始终保留完整链的cartoon选择，仅近邻stick采用4.3 Å距离选择。某些链会自然延伸到近景边界之外，这是取景边界，不是空间球体切碎cartoon。局部两态按位点附近±5残基的N/CA/C原子拟合后使用同一相机及显示比例。相机使侧链长轴处于画面平面内，以便看清真实原子连接；没有改变原子构象。`M72_*`/`M109_*`不带detail的图保留全链在视野内，可作为宽视角备选，但重点侧链较小。

`CaM_continuous_context_v80.pse`为可编辑场景。运行顺序是`render_structures_v80.py`、`render_local_details.py`、`render_tight_focus.py`；环境为`D:\Anaconda\envs\pvb_cpu\python.exe`，脚本已经配置其Library/bin DLL目录。后两个脚本依赖第一个脚本生成的会话。原始输入副本保存在`inputs/`并附SHA256。

泛素素材`ubiquitin_Q41_Q2_locator_trim.png`和同名PSE取自真实A003R模型输入，Ub链A为72残基、UBA链B为47残基；Ub为蓝色、UBA为灰色，Q41暖赭、Q2青绿。这是Q41/Q2真实静态定位图，没有把旧Q62的MD帧冒充Q2。其相机及紧裁标注坐标在`ubiquitin_Q41_Q2_locator.json`。相应复现脚本为`render_ubiquitin_locator.py`。

`audit_continuity.py`逐个检查所有输入链的残基编号、N/CA/C骨架原子和相邻残基距离。结果`backbone_continuity_QA.json`表明，CaM两态的145残基、NHE1的66残基、IQCG的37残基、Ub的72残基及UBA的47残基均没有编号缺口或骨架原子缺失；相邻CA没有超过4.5 Å的间断，肽键C–N没有超过1.8 Å的间断。没有通过强行加键掩盖真实缺口。直接查看透明PNG时部分查看器会忽略或错误显示alpha，请以白底合成后的QA预览或最终Matplotlib组图判断透明度。

注意：伙伴链沿用模型输入的全局残基编号，IQCG如GLN163不等于原始生物学IQCG编号。没有另行核对作者编号时，只标注伙伴名称和成熟CaM位点，避免把全局输入编号误写成生物学位点。
