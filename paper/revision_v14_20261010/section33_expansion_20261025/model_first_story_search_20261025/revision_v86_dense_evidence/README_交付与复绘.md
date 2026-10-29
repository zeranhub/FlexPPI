# 1MXE图3完整接触与侧链环境

最终主图为figure3/fig3_1MXE_dense_evidence.pdf（183×225 mm，Arial嵌入），同名SVG保留可编辑文本，PNG为450 dpi。图3已替换进paper_v86_20261026_中文修订_图3完整接触与侧链环境.docx。图注在figure3/caption_zh.md，正文在narrative/3.3_定稿.md。完整排版29页，经渲染核查；内部页图不作为交付图稿。

本版d为全部25×144的起始结构接触图：≤4 Å深青，4–4.5 Å浅青，>4.5 Å白；53+26+3521格，无位置抽样。e为全部21个有定义Hχ的侧链，按类型分组，以蓝圆/暖赭方点对比Full/Intra，右侧青条表示静态CaM遮蔽比例。f为1A06真实自抑制CaMKI与1MXE结合短肽结构，对比共有294–316片段，1A06缺失317–318不补造；标出M316邻近L45/S99。虚线问号为机制假说。c是完整25/21位置的绝对预测曲线，使用结构锚点定义的区域带与同坐标序列，未按预测高低筛区。

模型与MD结果沿用已完成、经审计的全组成链家族留出第1折预测；本轮没有重新训练或声称对1A06新增推理、MD或激活实验。Full与Intra是同一结合态输入的独立训练模型，不能解释为游离/结合实验状态。1MXE不含激酶催化域，f不表示模拟观察到激活或重建抑制接触。S3和其他原图保留，图3在主文中不绘MD参考。

复绘主图：用包含numpy/pandas/matplotlib/Pillow的Python运行figure3/compose_fig3_1MXE.py，脚本使用包内相对路径。当前可用解释器为E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe。结构重绘运行pymol/render_1MXE_v85.py与pymol_context/render_structural_contexts.py；后者一次生成最终透明度、接触侧链棒、PSE、裁边图和landmarks，运行环境D:/Anaconda/envs/pvb_cpu/python.exe。两结构全部原子只作刚体变换，不扭曲片段；PSE验证见pymol_context/structure_QA.json。

data_audit保留完整模型/MD与结构来源，data_design保留完整距离与静态环境字段；figure3/source_data是主图的冻结直接输入。narrative/paragraph_replacements.json是Word构建的最终文本来源；manuscript/v86_text_plan.json与v86_build_audit.json记录最终内容及最小替换审计。源论文保持不变。

验证通过：主图全部3600格、21侧链、92曲线值逐项精确；所有柔性坐标从零开始；字体嵌入；PNG无越界文字；Word图文同页，无新增孤页。29页中25页与已审定v85逐像素相同，4页仅为授权图3、正文及方法变化。自动样式检查的三个模式识别误报及人工核查理由见figure3/figure_QA.json。

参考图与具体借鉴见图3期刊参考与改图依据.md。
