# v110：删旧d、缩窄原f，压成五面板

原d的六Gln环境对比从主图删除，原f移入右侧空位并顺序标为d；密度主轴从157mm宽缩为60mm宽、高31mm。a/b/c/e位置关系与数据保留，e仍四张横排。整图从183×238mm缩为183×181mm。没有用截图裁切替代重排，所有曲线/坐标/文本为原生矢量。

数据：主图原785项非密度读数只移除旧d的18项（12个Hχ、6个SASA），其余767逐项严格相同。移除值存removed_old_d_18_values.csv，原source_data完整保留以追溯。368曲线值、97SASA条值、151个χ观察和全量824988密度原坐标/mask完全未变，源缓存与v109哈希相同。新d双方共用50×35网格，密度统计仍视野内820754/824082、全量rho.208615155/.261904889；全部非零格画出，零格透明，无等值轮廓。网格较旧底行80×35改变，因此每格峰值不同属于显示分箱，不是预测变化。a分子图未重绘，b/c静态SASA未重算，Word未改。

复现：assemble_main.py从v109生成figure3/compose_fig3_v110.py；后者直接用包内冻结输入、case_profiles.py和plot_overlay.py重画。源inputs、PyMOL静态几何表/脚本和完整NPZ保留于包内。release_and_audit.py核对767值完全保留、全部非零密度格、字体嵌入及Word哈希，导出PNG450dpi/PDF/SVG并打包。图注见caption_layout_revision_zh.md。
