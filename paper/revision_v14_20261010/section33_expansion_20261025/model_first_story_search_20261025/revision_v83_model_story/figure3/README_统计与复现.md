# 图3 v83：以模型预测解释CaM状态变化

主图所有定量读数均来自同一冻结第0折Full预测，没有重新训练或重新推理。c只保留预测RMSF和Hχ两条全序列双态曲线；d、e、f位于c上方，依次为完整四接触组、全部九个Met、M72/M109双读数配对。原MD参考和准确度散点已移入补充图S3，避免混淆案例叙述与准确度验证。图内没有数值表或灰色说明脚注。

配色基于论文实际图1/2/4/5及publication_style.py：NHE1暖赭#B87936圆点，IQCG蓝#2166AC方点；结构CaM中性灰、伙伴紫#8B72A5、关键Met青绿#3A8D8C。组分颜色与状态颜色分开，状态标题在结构上使用黑色。色板来源和对照记录见../palette_audit/palette_decision.json。

所有145个RMSF与122个有效Hχ位置完整保留；无效Hχ位置断线。三个框固定为预先展示位点M72/M109/M144的±4残基邻域：68–76、105–113、140–148。全模型原值记录在source_data/all_plotted_model_story_quantities.csv，共568条绘图值；key_neighbourhood_model_summary.csv给出固定窗口统计，model_story_summary.json给出全局统计。d为四组均值，组内位置数8/21/3/90；e保留全部9个同类型侧链；f数字为原值四舍五入至两位，原始精度保留源表。图中没有置信区间或显著性检验。

M72的预测主链变化小而Hχ升高；M109两项降低；M144主链降低而Hχ升高。模型的135/10个RMSF降低/升高与60/62个Hχ升高/降低，不能混用MD的增减计数。Hχ为固定角区下多角联合占据熵，不等同于转动速率。RMSF沿整链拟合标签口径，包含叶间运动。两记录均在训练中，属于结构解释案例；伙伴、构型及原始钙条件同时不同，不能做单因素因果归因。

结构由PyMOL真实重渲，不使用图片换色。完整复合物只作一次刚体变换，同一参考视角和共同裁边；局部用原保存场景。所有相机与坐标核验见../pymol/render_manifest.json，青绿版两套PSE及复现脚本位于该目录。

运行python compose_fig3.py，依赖numpy/pandas/matplotlib/Pillow。输出183×190mm、450dpi PNG、嵌入TrueType字体的PDF、可编辑文字SVG。补充图复现见../supplement_cam。
