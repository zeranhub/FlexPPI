# v89 Figure 3 independent QA

审查对象：实际查看首版fig3_mechanism_story.png（SHA256 30bacaac8c86b00df425ca04d264e6751deb075445c65fb7e2c53d5bc83ba561），并读取当前composer、render_manifest、figure_QA和全部绘图点。注意当前composer已补入局部接触线、叶瓣功能与N297/N308标签，首版PNG尚未包含这些新改动；下面视觉建议针对首版，数据核验针对当前导出CSV。

## 数据核验：通过

source_data/all25_dynamics_and_static_lobe_contacts.csv与v88冻结源文件逐字节相同，SHA256均为ebc5700ffcb1708f94d632816e854ea38e5397e7edba35d2d27b59503b36c473。独立脚本重新按原生位点、模型和指标构建预期92个数值，与all_plotted_dynamics.csv一对一匹配，最大绝对差为0，没有重复、丢点或替换。25主链位置各有Full/Intra，共50值；21有效Hχ位置各有Full/Intra，共42值；仅排除没有有效侧链χ的A299/A306/A309/A311。所有数值非负，RMSF轴0–8.2 Å与Hχ轴0–4.3 nats未裁掉数据。源数据为家族留出fold 1，未变更为游离态或另一折。

独立验证结果写入independent_data_QA.json，脚本independent_data_QA.py可复算。图内高亮圆点/注释重复绘制同一位置只是强调，不改变完整92个独立模型读数。

## 标签与颜色：无实质编号错误

首版真实结构图中的K296/F298、W303/M316标签与原生编号一致，1A06明确为Autoinhibited CaMKI，1MXE明确为CaM-bound regulatory peptide，没有把1MXE短肽复合物冒称完整活化激酶。F298口袋的I210/G214指向真实局部结构，未把K296强行画成同一口袋的锁扣。首位曲线为294，与实际I294相符，图中没有错写K294。当前composer准备增加的N297/N308也都是实际Asn且χ有效。

结构用灰蛋白、紫调控肽、青功能位点，模型用Full蓝圆实线和Intra赭方虚线，两种语义未重用。黑轮廓只强调Full选定位点，不是第三个模型。标题和主要图注字号可读，没有数据表或灰色解释小字。两条绝对指标曲线仍保持全25/21位置，符合既有要求。

## 需要修改的具体点

1. 首版a标题“The regulatory overlap binds two different partners”容易把同一酶中的催化核心叫成独立partner。当前composer已改“From inhibitory contacts to CaM capture”，更准确，应保留此类状态/接触关系表述。两组结构之间不要加带时间步骤的箭头。

2. b标题“Partner context restrains the dual-anchor core”容易读成已经通过物理扰动证明伙伴导致约束；Full/Intra实际是同一bound输入上的两模型。建议改为“Backbone restraint with heterogeneous side-chain states”，或“Joint flexibility resolves the regulatory segment”。Full/Intra含义仍需图注明确，避免读者将赭线看成真实游离态。

3. 首版c把K296/F298画成CaM两叶之外的完全自由肽段，容易与真实1MXE中它们也接触CaM的事实冲突。root准备补的重叠区—C叶瓣接触线必要，但它只能表示有接触，不应画成第三个已量化锚点或用相同圆点把它升级为W303/M316等价锚点。标C叶“Tight binding”、N叶“Activation required”能直接解释为什么总体结合不等于激活。

4. c首版指向催化槽的虚线箭头视觉上像“正发生回接”，而下方“Less re-contact?”悬在图底，二者关系不够清楚。可将路径标为“Possible inhibitory re-contact”，让总标题的“may limit”承担减少回接的假说；也可改成无方向虚线联系，不把假说画成已观测移动。不要用实线禁入符号装成回接已经被阻断。

5. c测试框的“Matched CaM occupancy”很关键：它指实际结合比例相近，而不是相同CaM浓度。“Core restraint ↓”与“Kinase activity ↓”均为预期，应保留Testable prediction/Predicted总标签，不能给后者加已验证的肯定符号。单有活性下降不足以验证运动约束这个中介，正文应要求同步测局部约束或抑制接触重建。

6. b的RMSF为了完整展示Intra最大7.606 Å，Full核心0.45–0.69 Å在纵轴上必然较平。这是完整尺度的代价；不建议再次堆新小图或截断坐标轴。正文应给出中心与侧翼的清楚绝对读数，机制主句只承担“连续受到约束”，不要依赖肉眼从这条较平曲线分辨F298的强弱功能角色。

## 是否真正比旧d/e/f强

是。新c有清楚的必要性：a提供调控序列既与催化核心形成抑制接触、又可被CaM捕获的结构背景，b显示捕获结构上局部主链/侧链的联合读数，c提出“总体结合保留而局部抓持变弱，是否仍降低激活”的区分实验。它没有再把同一批点算成多个摘要图。尚未被证明的是“低RMSF降低回接概率”，因此Hypothesis与may不可删除，预期的活性箭头不可当成模型已经预测了活性。本图读者应获得一个可检验的功能解释，而不是看到已确定的激活路径。

本次只写独立QA及独立数据核验文件；未修改root composer、PNG、Word或冻结数据。