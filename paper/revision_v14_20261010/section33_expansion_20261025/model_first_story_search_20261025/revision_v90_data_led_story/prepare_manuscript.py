from pathlib import Path
import json,shutil
S=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');O=S/'revision_v89_mechanism_story';V=S/'revision_v90_data_led_story'
old=json.loads((O/'narrative/paragraph_replacements.json').read_text(encoding='utf-8'))['paragraph_replacements']
r={}
r['71']='主链层面，Full预测的低波动区域覆盖连续的识别核心（图3b）。以W303–M316双锚点定义核心，全部14个位置的Full RMSF均为0.447–0.689 Å，平均0.513 Å；核心外11个位置平均为0.919 Å。Intra的核心及外围均值分别为2.666和3.750 Å。结合态MD中，核心14个位置均为0.453–0.624 Å，核心及外围均值为0.553和1.089 Å，支持连续低波动核心及其与外围的差异（图S3）。图3b另以较小纵轴范围重绘完整25位Full RMSF，便于读取原双曲线中被高值量程压缩的位点差异；两处为同一组数据。两模型均采用全组成链家族划分第1折留出权重，CaM相关家族未参与训练；Full与Intra分别利用完整复合物和链内编码预测同一结合态输入，不代表实测结合态与游离态。'
r['75']='这些数据提供的功能线索是，低波动区域连续跨越W303与M316两个锚点，并非只出现在单个接触位点；两个Asn的结果则表明，主链约束仍可伴随侧链状态的位置差异。结构将该模式定位到具有自抑制作用的调控片段：1A06中F298等位置与催化核心接触，1MXE中该片段转而由CaM识别（图3a）[29,30]。结合C端叶瓣能够紧密结合、N端叶瓣对激活仍不可缺少的实验结果[29]，可提出一种机制假说：双叶瓣对调控片段的持续约束，可能有助于维持其脱离催化核心后的构型，减少抑制接触重新形成的机会；N端叶瓣对M316一侧的抓持可能参与这一过程。这里的推断来自连续柔性模式与既有功能背景的对应，侧链位置差异本身并不证明激酶激活。'
r['76']='该假说可在CaM结合占据率相近的条件下，通过扰动N端叶瓣–M316一侧并联测局部运动与激酶活性来检验。1MXE输入仅含CaM与25残基靶肽，未包含催化域，也没有匹配的游离态轨迹；1A06仅提供静态背景，未作模型预测或MD。因此，现有数据尚未直接观测调控片段回接触催化核心或活化过程。'
r['73']='图3：主链约束与侧链位置差异的联合预测。（a）自抑制CaMKI（1A06）与CaM结合调控肽（1MXE）的真实结构。灰色表示激酶或CaM，紫色表示所示调控片段，青色标示功能位点。F298–激酶接触在全景右侧放大；K296/F298所在序列同时参与CaM识别。1MXE中W303、M316分别接触CaM的C端、N端叶瓣。两结构不是同一轨迹的连续帧，1MXE不含催化域。（b）同一1MXE结合态输入的Full与Intra预测，分别为蓝圆实线与暖赭方点虚线。左上展示完整25位RMSF双曲线；左下仅将同一组Full的全部25位RMSF在0–2.3 Å纵轴范围重绘，以提高低值区域的可读性，不是另一组数据。右侧保留21个有效Hχ位置，四个无χ角的Ala位置留空。浅底色及细条标出W303–M316核心，并标示K296/F298、两个锚点及两个Asn。Full/Intra不代表实测结合态与游离态；结合态MD及误差核查见图S3，结构与功能背景见文献[29,30]。'
(V/'narrative/paragraph_replacements.json').write_text(json.dumps({'source':'Actual v89 source 186 body paragraphs','paragraph_replacements':r},ensure_ascii=False,indent=2),encoding='utf-8')
b=(O/'build_v89_manuscript.py').read_text(encoding='utf-8').replace('revision_v89_mechanism_story','revision_v90_data_led_story')
b=b.replace("SOURCE=R/'paper_v88_20261026_中文修订_图3双面板主线.docx';OUTPUT=R/'paper_v89_20261026_中文修订_图3调控机制假说.docx'","SOURCE=R/'paper_v89_20261026_中文修订_图3调控机制假说.docx';OUTPUT=R/'paper_v90_20261026_中文修订_图3数据呈现与机制推断.docx'")
b=b.replace('fig3_mechanism_story.png','fig3_data_led_story.png').replace('Observed inhibitory and CaM-bound structures, complete model profiles, and a testable regulatory hypothesis','Observed structural context, complete RMSF and Hchi model profiles, and repeated Full RMSF on an expanded y scale')
b=b.replace('v89_text_plan','v90_text_plan').replace('v89_final_text','v90_final_text').replace('v89_build_audit','v90_build_audit').replace('Actual v88 source paragraph numbers','Actual v89 source paragraph numbers')
b=b.replace("V.parent/'revision_v88_two_panel_story/render'","V.parent/'revision_v89_mechanism_story/render'")
(V/'build_v90_manuscript.py').write_text(b,encoding='utf-8')
# Save prepare source inside the output, not the repository root.
shutil.copy2(Path(r'C:\Users\DELL\.codex\worktrees\5752\FlexPPI\tmp_v90_prepare.py'),V/'prepare_figure.py')
print(json.dumps({'prepared_paragraphs':list(r),'mechanism_schematic_removed':True},ensure_ascii=False))