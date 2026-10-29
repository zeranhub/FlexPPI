# 图5重组、统计与复现说明

本次直接继承已发表稿候选图5的冻结数值源，重组为183×230 mm的六个panel。配色统一为Full蓝#2166AC、Intrachain暖赭#B87936、辅助青绿#3A8D8C。Q41暖赭菱形、Q2青绿方形并辅以曲线实虚线，不依赖颜色单独判读。英文图内文字最小9 pt；在Word以167.64 mm宽插入后约8.24 pt。

变化集中在d–f：d把真实PyMOL完整泛素–UBA输入结构定位与全部66位点预测/MD曲线放在同一panel；e把四种侧链角按χ1至χ4排成同一横排，共享相同轴范围；f保留全部原角度箱、零计数断线与低F*角区宽度。原Q41/Q62案例不混入本图。结构是静态模型输入，不代表MD帧；f是MD经验投影自由能，不是模型预测的自由能。

## 各panel的n、统计量与来源

(a) n=3,063条复合物/来源记录，MD和预测各有3,063个完整ECDF顶点；每个值为该记录内主链RMSF和侧链Hχ经同氨基酸类型比较后的Spearman相关。图中标注均值MD=0.391、模型=0.517；无误差条、无本次新增假设检验。源：source/original_fig4/fig4_joint_complementarity_all_ecdf_points.csv，x/y分别是相关与ECDF分位。

(b) 每种物理读数n=83条配对记录，合计33,142个同类型有效残基位置；每点为同一记录中Hχ与RMSF预测对MD细箱熵或低F*角区宽度的相关之差。大点为记录均值，误差条沿用原来源分层重采样95%区间；没有新增检验或多重比较。源：source/entropy、source/width下的replica1_per_record_correlations.csv及paired_correlation_difference_summary.csv，筛选threshold=1、scope=all（熵另筛fine_entropy_mean_per_chi）。全部83对均显示，包括负差。此差值是描述性比较，不证明条件独立贡献。

(c) 每种预测器n=83条记录；高低半组只依据预测在相同氨基酸类型内划分，纵向每个点表示对应真实MD角区总宽度的高组减低组。大点/误差条为原均值/95%来源分层重采样区间。源：source/width/replica1_high_low_per_record.csv及replica1_high_low_summary.csv；Full Hχ、Intra Hχ、Full RMSF共249点，原轴范围−1.5至17.5°/χ保留。

(d) n=66个泛素Hχ有效位置，各有冻结MD参考与Full模型预测，共132值，位置范围1–72。10、28、35、46、47、53处为原Gly/Ala缺口，不插值也不连线。曲线参考为原五窗口联合熵均值；MAE=0.0795305 nats（显示0.080），Spearman ρ=0.979375（显示0.979）。源：source/original_fig7/source/all66_residue_summary.csv与ubiquitin_all66_chi_501frames.npz冻结预测/参考字段。结构来自pymol目录的真实A003R输入，Ub72与UBA47均完整显示，PSE可编辑。

(e) χ1/χ2/χ3/χ4的n依次66/52/23/10，共151个已定义残基–角度观察；全部保留。横轴是同一个联合Hχ模型输出的氨基酸类型去均值，纵轴为MD单角细箱熵的类型去均值，不能解读成模型具有四个逐角预测头。相关沿用原Spearman值0.50/0.45/0.52/0.13；无新增显著性检验。所有小图都保留原x范围[−0.70,0.50]、y范围[−1.0,0.75]。源：source/original_fig7/experiment/defined_angles_long.csv与per_chi_correlations.csv。

(f) 由六个Gln预测极端端点选Q41/Q2，预测值2.000081/2.989756 nats。每个位点3个χ角、每角72个5°箱，共432箱，来源为同一0–100 ns共同窗口501帧。F*=−ln(p̂/p̂max)为无量纲经验投影；211个原零计数箱保留NaN断线，不加伪计数、不平滑、不替换。原点0°的F*≤1总宽度分别为Q41的20/25/35°和Q2的30/55/100°。图中是这个固定原点的宽度，不混用四原点平均27.50/64.17°/χ。曲线全保留原范围；数值源histograms_501frames.csv。

本次重绘没有重新训练、重新预测、改变种子/划分、抽样或重新挑选轨迹。使用原论文冻结模型评估结果；训练/验证安排仍按论文方法。残基、轨迹帧和记录的n含义不同，不能混为独立生物重复。报告保留原来源定义，不能把一个案例的66个位点解释为66次独立实验。

## 核验与复现

compose_fig5.py由original_production.py的全部数据加载、数值校验与三段逐字风格基线继承，后半绘图布局明确重写；compose_body.py保留新布局。source/是31个原数据与源文件的逐字节副本，production_asset_sources/保留原风格参数来源；原通用LineTrend与PairedBoxScatter示例数据不匹配本科学数据，继承已有科学生产脚本和样式参数，没有使用示例数据。此差异已在原figure_contract_and_assets.json和脚本头中记录。

verify_figure5.py逐列逐行比较新旧all_plotted_points.csv，全部7,256行的坐标、数值、NaN与标识精确一致，不只是样本计数一致。确认PDF物理页183×230 mm、Arial/Arial Bold均TrueType嵌入；SVG只有一个image元素，即真实PyMOL结构，其余全部曲线、散点、字与标注为可编辑矢量。450 dpi PNG为3242×4074。最终PNG已实际view_image查看，修复结构标签引线及f数值与曲线的干扰，无图内文字出界。

运行：使用本项目runtime_py312/Scripts/python.exe依次执行compose_fig5.py和verify_figure5.py；输出fig5_joint_readouts_conformations.pdf/.svg/.png，中文图注caption_zh.md，完整坐标all_plotted_points.csv，统计/数值核验data_audit.json与verification_summary.json。结构像素坐标每次从当前PyMOL JSON读取，避免分辨率升级后标注错位。

## 四轮QA结论

Pass 0：三段风格基线保留；无默认/彩虹色、无四边框、无装饰性格线。d内部图例位于右上空白，不覆盖数据；a/f的直接标注使用空白区。Pass 1：9 pt字体、183×230 mm、450 dpi、TrueType嵌入和矢量输出均通过。Pass 2：从跨记录互补性到完整Ub验证、逐角物理含义和预测选点的MD占据逐层递进；全部数据保留、轴范围未裁点、统计定义与局限明确。Pass 3：实际渲染图检查、文本边界检查、PNG尺寸/PDF尺寸/字体与唯一raster结构检查通过。最终主文页面中的排版仍须由根任务执行整篇DOCX渲染检查。
