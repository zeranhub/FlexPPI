# 主链与侧链柔性的互补解析

本版最终Word为paper_v87_20261026_中文修订_互补柔性背景与图3凝练.docx。3.3标题改为“主链与侧链柔性的互补解析”，增加一般问题背景后引入CaM–CaMKI作为实例，并解释CaM双叶瓣、CaMKI自抑制及已有叶瓣功能实验。最终正文在narrative/3.3_定稿.md，唯一构建文本来源为narrative/paragraph_replacements.json；manuscript/v87_text_plan.json及v87_build_audit.json记录最小替换与新增段落。

图3最终为figure3/fig3_joint_flexibility_compact.pdf/svg/png，尺寸183×190 mm，较v86的183×225 mm减高15.6%，保持原有字号与全部25个主链、21个有效侧链位置。凝练来自信息整合及布局收紧：d使用完整3600个距离统计25位置的N/C/连接段接触残基数，含零位置；完整25×144接触矩阵移入figure3/source_contact_map，不在主图重复。e保留21侧链按类型配对的Full/Intra Hχ与静态遮蔽。f删除重复1MXE小全景，a/b已给出其结合环境；改为1A06整体及M316–L45/S99真实局部放大，整体/局部完全同旋转方向。虚线为重原子距离，不称氢键。图3不再用末端问号概括机制，机制假说由正文结合功能实验提出。

主图复绘：运行figure3/compose_fig3_1MXE.py，脚本使用包内相对路径，并同时生成完整接触源图。可用数值环境为E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe。PyMOL局部由pymol_context_local/render_M316_local.py复绘，整体由pymol及pymol_context主脚本复绘；环境D:/Anaconda/envs/pvb_cpu/python.exe。所有结构均为真实实验坐标，只平移/旋转/缩放相机，不造缺失残基、不扭曲构象。

全部模型与MD读数沿用已完成的全组成链家族留出第1折结果，本轮未新训模型或对1A06增加预测/MD。Full与Intra是同一结合态输入的独立模型，不能视作实验游离/结合状态。1A06未解析317–318，未补建；1MXE只含CaM与短肽，不含完整激酶。M316相关构型维持为可检验解释，尚未直接观察回接触或活性。S3图像与其他原图不变，仅同步其图注对3.3假说的指向。

源数据在figure3/source_data及data_audit，原生编号、全部距离/计数、模型预测和方法均可追溯。独立数值与视觉审查分别记录于final_figure_review.md、narrative/final_text_and_pages10_18_review.md及pymol_context_local/final_supp_pages_visual_QA.md。PDF嵌入Arial字体，PNG为450dpi；无表格、无灰色解释小字、主图无MD参考曲线。Word仅替换授权段落与图3，六张表和其他图像的文件字节均保持不变。
