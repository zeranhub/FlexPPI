# v109：双指标曲线下增加对齐的溶剂暴露柱图

主图b、c保留同残基横轴的双纵轴曲线，并在正下方增加青灰色侧链SASA柱图。原b的CaM遮蔽百分比已换成侧链溶剂可及面积，c新增同口径柱图。上下层严格对齐，仅底层显示残基编号。SASA两例统一0–160Å²；RMSF0–2.4Å/Hχ0–4.5nats仍不变。主曲线368有效读数、mask、合法零和孤立点均不变。a/d/e/f数据保持；e/f向下移8mm为SASA轨道留空间。Word未修改。

SASA来源及统一：泛素all119_input_environment.csv的A链72位与d的六Gln同源；原CaM SASA由MDTraj计算，因此本版用PyMOL3.0.0在精确原始CaM复合物重算25位，采用与泛素一致的1.4Å探针/dot_density4，复合物全部重原子参与遮挡，侧链排除N CA C O OXT OT1 OT2。新几何表与输入哈希、PyMOL版本、参数及脚本位于sasa_audit/。只做静态几何，不是模型推理/MD；SASA不受Hχ mask限制，97位置全部绘出。泛素G10/G35/G47/G53无侧链重原子，侧链SASA0表示几何定义，不是深埋；A28/A46仍有SASA柱。case_SASA_all97_starting_complex.csv保存全部条值。不要将本版b下方继续写成CaM遮蔽比例，不是同坐标删除CaM的ΔSASA。

来源：source_inputs/cam_all25.csv四列与原家族fold1源一致，25RMSF/21Hχ有效，原生294–318，A299/306/309/311 Hχ无定义。泛素源source_inputs/original_frozen_Full_A003R.npz与label NPZ，chain A72位，target[:,0/1]和mean[:,0/1]，valid[:,0/1]原mask；72RMSF/66Hχ有效，G10/A28/G35/A46/G47/G53 Hχ留空，P19/P37/P38的MD合法零及Full近零原值保留。case_profiles_all97_with_masks.csv保存97位置四项原读数及mask，standalone_all368_plotted_values.csv保存全部368有效曲线读数。模型来源为现成对应留出折单模型，CaMKI familyfold1与Ub原fold0不是同一套划分；详细哈希和MAE在case_profile_manifest.json。无新模型或MD计算。

准确度：CaMKI RMSF MAE0.1148592012Å、Hχ0.254417929857nats；泛素RMSF0.039814019369Å、Hχ0.079530506498nats。原v108主图713值中25个旧遮蔽值退出本版主图，余688值逐项保留；97个SASA值加入，总785值。f脚本、NPZ和计数CSV均与v108逐字节一致，保持无轮廓浅色透明叠加。

复现：运行case_profiles.py制作独立双案例曲线及SASA轨道，assemble_main.py生成figure3/compose_fig3_v109.py，后者绘制整图。所有冻结源和新几何输入/脚本包含于包内；已有缓存可直接画，不需重跑几何、模型或MD。PDF/SVG原生矢量，Arial嵌入，PNG450dpi。release_and_audit.py检查368曲线值、688旧值保留、97SASA条值、源哈希、f不变及字体并打包。
