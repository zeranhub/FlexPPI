# v108：每个案例同时展示主链RMSF与侧链Hχ

主图b、c改为同残基横轴的双纵轴曲线：蓝色左轴RMSF（Å）、赭色右轴Hχ（nats），与图2a指标色一致；Full实线、MD虚线。两例都用RMSF0–2.4Å/Hχ0–4.5nats，从0起，不通过单独缩放制造对应趋势。Hχ小点保留无定义位置两侧的孤立有效点，不跨缺口连线。两指标单位不同，曲线交叉或峰谷接近不表示相等、因果耦合或数值可互换。

来源：source_inputs/cam_all25.csv四列与原家族fold1源一致，25RMSF/21Hχ有效，原生294–318，A299/306/309/311 Hχ无定义。泛素源source_inputs/original_frozen_Full_A003R.npz与label NPZ，chain A72位，target[:,0/1]和mean[:,0/1]，valid[:,0/1]原mask；72RMSF/66Hχ有效，G10/A28/G35/A46/G47/G53 Hχ留空，P19/P37/P38的MD合法零及Full近零原值保留。case_profiles_all97_with_masks.csv保存97位置四项原读数及mask，standalone_all368_plotted_values.csv保存全部368有效曲线读数。模型来源为现成对应留出折单模型，CaMKI familyfold1与Ub原fold0不是同一套划分；详细哈希和MAE在case_profile_manifest.json。无新模型或MD计算。

准确度：CaMKI RMSF MAE0.1148592012Å、Hχ0.254417929857nats；泛素RMSF0.039814019369Å、Hχ0.079530506498nats。CaMKI Hχ的准确度应使用这里21有效位置的值，不能沿用不同范围旧值。主图原527项科学读数保留，新增42个CaMKI Hχ值及144个泛素RMSF值，共713项；旧导出CSV在少数泛素Hχ末位存在最大约4.44e-16读写舍入，本次严格恢复原NPZ精度。a分子结构、d同类侧链环境、e四χ横排、f无轮廓浅色叠加均保留。f脚本、NPZ和计数CSV均与v107逐字节一致。

复现：运行case_profiles.py制作独立双案例曲线，assemble_main.py生成figure3/compose_fig3_v108.py，后者绘制整图。冻结源均包含于包内；已有case_profiles_all97缓存可直接重画，f已有全量NPZ可复现。PDF/SVG原生矢量，Arial嵌入，PNG450dpi。Word未改。release_and_audit.py检查全部368曲线值、527旧值保留、源哈希、Pro/mask、f不变和字体，并打包。
