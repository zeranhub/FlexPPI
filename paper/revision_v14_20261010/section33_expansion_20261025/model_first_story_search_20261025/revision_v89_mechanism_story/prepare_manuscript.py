from pathlib import Path
import json
S=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025')
V=S/'revision_v89_mechanism_story';V88=S/'revision_v88_two_panel_story'
draft=(V/'narrative/proposed_3.3_and_caption.md').read_text(encoding='utf-8')
main=draft.split('# 图3图注')[0].strip().split('\n\n')
assert len(main)==7
repl={str(idx):text for idx,text in zip([70,71,74,75,76],main[2:])}
# Preserve the existing general heading and broad background (paragraphs 68/69).
repl['70']='以钙调蛋白（CaM）识别CaM依赖性蛋白激酶I（CaMKI）的调控片段为例。自抑制结构1A06显示，调控区跨越催化核心并干扰底物及ATP结合，其K296/F298所在序列与CaM结合区重叠[29,30]。孤立调控肽还能抑制去除调控区的CaMKI293，表明该序列本身保留了对催化核心的抑制作用[29]。在1MXE中，同一序列改由CaM识别，W303与M316分别接触C端和N端叶瓣，形成1–14双锚点构型，K296/F298也直接接触CaM[29]（图3a）。已有实验发现，C端叶瓣能够紧密结合调控序列，激酶激活仍需要N端叶瓣[29]。因此，本例进一步考察：双叶瓣抓持产生的局部约束，能否为维持解除自抑制后的调控构型提供解释？'
repl['71']='主链层面，Full预测的低波动区域覆盖连续的识别核心（图3b）。以W303–M316双锚点定义核心，14个位置的Full平均RMSF为0.513 Å，核心外11个位置为0.919 Å；Intra分别为2.666和3.750 Å。结合态MD的核心及外围均值为0.553和1.089 Å，支持Full恢复的低波动核心及其与外围的差异（图S3）。这一模式跨越两个锚点，体现连续调控片段的主链约束。两模型均采用全组成链家族划分第1折留出权重，CaM相关家族未参与训练；Full与Intra分别利用完整复合物和链内编码预测同一结合态输入，不代表实测结合态与游离态。'
repl['74']='侧链层面，受到约束的主链仍容纳不同的位置状态（图3b）。该调控肽中全部两个Asn构成同类型对照：N297与N308的Full RMSF分别为0.809和0.448 Å，Hχ为1.587和0.601 nats。两者残基种类及χ角数相同，但N297侧链未被CaM遮蔽，N308遮蔽84.3%；MD的Hχ分别为1.656和0.912 nats，保留相同高低关系（图S3）。由此，连续主链的低波动与侧链状态多样性并非同一层面的信息，局部环境仍可进一步区分具体位置。静态埋藏与Hχ并非普遍单调，不能仅由埋藏程度推断侧链柔性；Hχ反映固定角区下的联合占据多样性，不表示转动速度。'
repl['75']='上述柔性特征可与调控片段的功能角色联系起来。1A06中F298等调控位点与催化核心接触；1MXE中，W303–M316则跨越CaM双叶瓣，调控片段由另一组相互作用稳定[29,30]（图3a）。结合C端叶瓣紧密结合而N端叶瓣对激活仍不可缺少的实验事实，可以提出一种机制假说：双锚点抓持通过约束调控片段，可能帮助维持其脱离催化核心后的构型，减少自抑制接触重新形成的机会；N端叶瓣对M316一侧的约束可能参与这一过程（图3c）。这为结合与激活的功能差别提供了可检验的动态解释，并不规定两个叶瓣结合的时间顺序。'
repl['76']='该假说预期：在CaM结合占据率相近时，削弱N端叶瓣–M316一侧的相互作用，会减弱局部构象约束并降低CaMKI活化（图3c）。需共同测量结合比例、局部运动与激酶活性；若约束改变但活化不受影响，或活化下降完全由结合减少解释，则不能支持这条联系。1MXE输入仅含CaM与25残基靶肽，未包含催化域，也没有匹配的游离态轨迹；1A06仅提供静态背景，未作模型预测或MD。因此，调控片段回接触及其对激活的影响仍是待检验假说。'
repl['73']='图3：主链与侧链柔性的互补解析及调控机制假说。（a）自抑制CaMKI（1A06）与CaM结合调控肽（1MXE）的真实结构。灰色表示激酶或CaM，紫色表示所示调控片段，青色标示功能位点。F298–激酶接触在全景右侧放大；K296/F298所在序列同时参与CaM识别。1MXE中W303、M316分别接触CaM的C端、N端叶瓣。两结构不是同一轨迹的连续帧，1MXE不含催化域。（b）同一1MXE结合态输入的Full与Intra预测，蓝圆实线和暖赭方点虚线分别表示两模型；保留25个RMSF位置与21个有效Hχ位置，四个无χ角的Ala位置留空。浅底色及细条标出W303–M316核心，标示K296/F298与两个Asn。两模型不代表实测结合态与游离态。（c）双锚点约束可能减少自抑制接触重建的机制示意，形状不表示预测的完整激酶构象。叶瓣功能分工来自文献[29]；虚线回接触及右侧扰动结果均为待检验假说与预期。在结合占据率相近时，减弱N端叶瓣–M316接触，预期局部约束与激酶活化同时降低。结合态MD及误差核查见图S3，结构与功能背景见文献[29,30]。'
old=json.loads((V88/'narrative/paragraph_replacements.json').read_text(encoding='utf-8'))['paragraph_replacements']
supp=old['163']
supp=supp.replace('其M316接触按294–316调控段侧链与该段以外激酶重原子间的最短距离判定，4.5 Å内对象为L45/S99；两结构只作全体原子的刚体变换，不扭曲共有片段。','静态接触采用所选调控位点侧链重原子至背景残基重原子的最短距离，4.5 Å内保留全部接触对象。1A06中F298邻近G214、I210、F104、Y184、P216及Y215，最近为G214（2.975 Å）；K296最近的E292位于相邻调控序列，不能视作远端催化核心接触。M316的片段外接触为L45/S99。1MXE中K296最接近CaM E127（2.744 Å），F298在4.5 Å内接触L116/E120/M124。图3a两结构采用同一相机旋转，按版面分别缩放，F298接触区另作局部放大；全部原子仅作刚体变换，不扭曲共有片段，也不进行构象插值。')
repl['163']=supp
(V/'narrative/paragraph_replacements.json').write_text(json.dumps({'source':'Actual v88 source 186 body paragraphs','paragraph_replacements':repl},ensure_ascii=False,indent=2),encoding='utf-8')
oldbuild=(V88/'build_v88_manuscript.py').read_text(encoding='utf-8')
new=oldbuild.replace("revision_v88_two_panel_story';M=","revision_v89_mechanism_story';M=")
new=new.replace("paper_v87_20261026_中文修订_互补柔性背景与图3凝练.docx'","paper_v88_20261026_中文修订_图3双面板主线.docx'")
new=new.replace("paper_v88_20261026_中文修订_图3双面板主线.docx'\nNS=","paper_v89_20261026_中文修订_图3调控机制假说.docx'\nNS=")
new=new.replace('fig3_two_panel_story.png','fig3_mechanism_story.png').replace('Two panels: full structural location and complete absolute flexibility profiles','Observed inhibitory and CaM-bound structures, complete model profiles, and a testable regulatory hypothesis').replace('v88_text_plan','v89_text_plan').replace('v88_final_text','v89_final_text').replace('v88_build_audit','v89_build_audit').replace('Actual v87 source paragraph numbers','Actual v88 source paragraph numbers')
new=new.replace("V.parent/'revision_v87_general_story_compact/render'","V.parent/'revision_v88_two_panel_story/render'")
new=new.replace("print(json.dumps(audit,ensure_ascii=False))","wrapper=render/'render_docx_word_wrapper.py'\nwrapper.write_text(wrapper.read_text(encoding='utf-8').replace('26.1024.11018','26.1025.11240'),encoding='utf-8')\nprint(json.dumps(audit,ensure_ascii=False))")
(V/'build_v89_manuscript.py').write_text(new,encoding='utf-8')
print(json.dumps({'prepared_paragraphs':list(repl),'main_chars':sum(len(repl[str(i)]) for i in [70,71,74,75,76]),'caption_chars':len(repl['73'])},ensure_ascii=False))