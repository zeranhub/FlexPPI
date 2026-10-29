# 1A06 / 1MXE：自抑制结构参照与CaM结合短肽

本目录专用于图3f的结构背景补证，不含新模型预测或新MD结果。官方1A06坐标来自 https://files.rcsb.org/download/1A06.pdb ，原始研究为 Goldberg、Nairn、Kuriyan，Cell 1996，DOI:10.1016/S0092-8674(00)81066-1。1MXE一侧使用实际家族留出模型输入的同一坐标，复制自v85，不是完整激酶的激活态结构。

编号核查：1A06链A为原生CaMKI编号，实际建模279个残基，观察范围10–316但中间有实验缺失区域。与1MXE肽294–318对应的294–316共23个残基均存在且序列完全相同；末尾317R、318K在1A06的REMARK465中明确标为实验未解析，没有补造。W303、M316、N297、N308的主链和所有侧链重原子完整。具体证据见native_numbering_audit.json及原始PDB。

图形：完整CaMKI或CaM作为灰色背景，实际存在的294–316/294–318调控段为紫色，W303和M316为青色；1A06的L45、S99侧链保留淡灰棒状，供标注其与M316的实际接近。1A06背景cartoon为0.43透明以显示被前方激酶骨架遮住的调控片段；紫色调控段不透明。两张结构的共同23个CA进行一次全对象刚体拟合，RMSD5.086 Å反映两个真实片段构象不同；没有局部扭转、插值或强行变形。原子对距离最大变化1.53e-5 Å，为PyMOL浮点误差。1MXE坐标完全未动。1A06全部原子统一刚体变换后的PDB另存，但原始与原生输入均保留。

两状态使用相同相机旋转方向，各自缩放以容纳完整体系，因此并不声明其显示比例完全相同。局部自抑制接触残基没有从整体结构中移出。PNG为透明背景。原始画布1500×700，trim是alpha边界向外15像素的纯裁切。render_manifest.json的assets/autoinhibited与assets/CaM_bound各含file、width、height、view、alpha_bbox、landmarks；assets/*/landmarks是原始PNG左上角像素。assets/*/trim另含file、crop_box、width、height、landmarks；trim/landmarks已减去裁切原点，是裁图后像素。W303指向CZ2原子，M316指向CE原子，L45/S99为侧链重原子质心。

可复算结构事实：1A06的W303侧链与294–316之外蛋白背景最近距离5.188 Å，在4.5 Å阈值下无接触残基。M316则与L45、S99接近：M316 SD–L45 CD1为4.143 Å，M316 CE–S99 OG为3.161 Å。这是空间接近证据，不等同于氢键、相互作用能或动力学因果。可在f标W303 exposed及M316–L45/S99，而不应把M316写成深埋。

辅助SASA审计保留在render_manifest.structural_descriptors中。PyMOL设置dot_solvent=1、solvent_radius=1.4 Å、dot_density=4，按完整对象计算每原子SASA，然后汇总关键残基侧链（排除N/CA/C/O）；以同状态实际调控片段单独提取的SASA作暴露比例分母。1A06的W303/M316侧链SASA为167.46/73.13 Å²，相对孤立调控片段保留97.4%/68.6%暴露。1MXE对应2.57/3.02 Å²，保留1.37%/2.74%暴露。这些静态量受可见坐标与当前构象影响，只作为结构环境核查，不把SASA直接解释为自由能，也不据此补出真实转变路径。

复现：D:\Anaconda\envs\pvb_cpu\python.exe运行render_structural_contexts.py，采用PyMOL3.0.0。audit_numbering.py核对编号/缺失，audit_specific_contacts.py独立用原始坐标复算具体距离。两份同名原生PSE均保留完整可编辑对象和相机，不依赖PNG截图。
最终主图f不显示任何SASA数值；只呈现两个真实结构背景与1A06中M316接近的L45/S99对象。辅助SASA仅保留用于溯源核查。render_structural_contexts.py已整合透明度、L45/S99侧链棒、完整PSE、trim和全部投影标记的最终代码，一次运行即可生成本目录最终资产，无须任何patch脚本。最后一次完整运行后，verify_context_assets.py已重新打开两份PSE并验证全部原子、相机、B因子保留与透明图尺寸。此后PNG冻结。
