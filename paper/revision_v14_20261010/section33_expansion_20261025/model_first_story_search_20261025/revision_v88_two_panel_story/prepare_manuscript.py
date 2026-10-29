from pathlib import Path
import json
s=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');v=s/'revision_v88_two_panel_story';proto=s/'revision_v88_three_panel_story';old=s/'revision_v87_general_story_compact'
(v/'narrative').mkdir(exist_ok=True)
j=json.loads((proto/'narrative/paragraph_replacements.json').read_text(encoding='utf-8'));repl=j['paragraph_replacements']
for key in ['71','74']:repl[key]=repl[key].replace('图3c','图3b').replace('图3a、c','图3a、b')
repl['70']=json.loads((old/'manuscript/v87_text_plan.json').read_text(encoding='utf-8'))['paragraph_replacements']['69'].replace('（图3a、b）','（图3a）')
repl['75']=repl['75'].replace('在1MXE中，它则被CaM的N端叶瓣识别（图3b）。','在1MXE中，它则被CaM的N端叶瓣识别[29,30]。')
repl['73']='图3：主链约束与侧链状态差异的联合预测。（a）CaM–CaMKI复合物1MXE的整体结构，CaM的N端和C端叶瓣为深浅灰色，调控肽为紫色，关键侧链为青色；定位W303、M316与两个Asn位置。（b）完整CaMKI294–318的25个主链位置及21个有效侧链位置的RMSF、Hχ绝对预测值；Full和Intra分别用蓝圆实线与暖赭方点虚线表示，均采用全组成链家族留出第1折权重，对同一结合态输入进行预测。序列与两条曲线共享残基坐标；浅底色和细条标出结构定义的W303–M316核心，标注N297与N308以对应同类型侧链的位置差异，四个无χ角的Ala位置留空。Full/Intra差异不代表实测结合态与游离态转变；MD核查见图S3，方法与静态结构背景见补充材料D。'
repl['163']=repl['163'].replace('静态位置对照（图3b）','静态位置对照')
j['paragraph_replacements']=repl;j['main_figure_panels']=['a structural location','b complete model profiles'];j['retired_main_panels']=['contact counts','21-sidechain paired display','autoinhibitory structural context']
(v/'narrative/paragraph_replacements.json').write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8')
code=(proto/'build_v88_manuscript.py').read_text(encoding='utf-8').replace('revision_v88_three_panel_story','revision_v88_two_panel_story').replace('paper_v88_20261026_中文修订_图3三面板主线.docx','paper_v88_20261026_中文修订_图3双面板主线.docx').replace('fig3_three_panel_story.png','fig3_two_panel_story.png').replace('Three panels: full structural location, M316 static contact contexts, complete absolute flexibility profiles','Two panels: full structural location and complete absolute flexibility profiles')
(v/'build_v88_manuscript.py').write_text(code,encoding='utf-8')
print('Prepared two-panel Word build.')