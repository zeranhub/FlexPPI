from pathlib import Path
import json
s=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');v=s/'revision_v88_three_panel_story';old=s/'revision_v87_general_story_compact'
j=json.loads((v/'narrative/v88_minimal_story_suggested_replacements.json').read_text(encoding='utf-8'));repl=j['paragraph_replacements']
repl['74']=repl['74'].replace('这一静态对照本身不揭示状态转换或柔性变化，但将本例的低波动核心定位到可参与自抑制接触的调控片段。','这一对照将本例的低波动核心对应到可参与自抑制接触的调控片段。')
repl['72']=repl['72'].replace('以共同调控片段确定观察方向的M316局部结构环境','M316的两种静态局部环境，采用相同观察方向与空间尺度').replace('两结构只作刚体变换，不补建缺失残基','两结构只作刚体变换；1A06未解析317–318，不补建')
repl['75']=repl['75'].replace('未新增预测或MD','未作模型预测或MD')
repl['162']=repl['162'].replace('也未对1A06新增模型预测或动力学标签','也未对1A06进行模型预测或生成动力学标签')
# Suggested keys refer to the original template. v87 inserted background at69.
actual={str(int(k)+1 if int(k)>=69 else int(k)):val for k,val in repl.items()}
(v/'narrative/paragraph_replacements.json').write_text(json.dumps({'index_basis':'Actual v87 source paragraphs, including its inserted background at69','paragraph_replacements':actual,'original_template_to_source_index':{k:int(k)+1 if int(k)>=69 else int(k) for k in repl}},ensure_ascii=False,indent=2),encoding='utf-8')
code=(old/'build_v87_manuscript.py').read_text(encoding='utf-8')
code=code.replace("revision_v87_general_story_compact'","revision_v88_three_panel_story'")
code=code.replace("SOURCE=R/'paper_v86_20261026_中文修订_图3完整接触与侧链环境.docx';OUTPUT=R/'paper_v87_20261026_中文修订_互补柔性背景与图3凝练.docx'","SOURCE=R/'paper_v87_20261026_中文修订_互补柔性背景与图3凝练.docx';OUTPUT=R/'paper_v88_20261026_中文修订_图3三面板主线.docx'")
code=code.replace('assert len(pars)==185','assert len(pars)==186').replace('p=pars[71]','p=pars[72]').replace("fig3_joint_flexibility_compact.png","fig3_three_panel_story.png")
code=code.replace('changed=set(map(int,repl))|{71}','changed=set(map(int,repl))|{72}')
code=code.replace('True 1MXE and 1A06 structural contacts, all25 lobe counts, all21 sidechains, absolute profiles','Three panels: full structural location, M316 static contact contexts, complete absolute flexibility profiles')
a=code.index('background=json.loads(');b=code.index("parts['word/document.xml']=",a);code=code[:a]+code[b:]
code=code.replace("'paragraph_replacements':repl,'insert_before_original_paragraph_69':background","'paragraph_replacements':repl,'index_basis':'Actual v87 source paragraph numbers'")
code=code.replace("'\\n'.join(finaltext[str(i)]+('\\n'+background if i==68 else '') for i in range(len(pars)))","'\\n'.join(finaltext.values())")
code=code.replace("'\\n\\n'.join([finaltext['68'],background]+[finaltext[str(i)] for i in [69,70,73,74,75]])","'\\n\\n'.join(finaltext[str(i)] for i in [68,69,70,71,74,75,76])")
code=code.replace("repl['72']+'\\n'","repl['73']+'\\n'")
code=code.replace(',inserted_background_before_original_paragraph=69','')
code=code.replace('v87_text_plan.json','v88_text_plan.json').replace('v87_final_text.txt','v88_final_text.txt').replace('v87_build_audit.json','v88_build_audit.json').replace('revision_v86_dense_evidence/render','revision_v87_general_story_compact/render')
(v/'build_v88_manuscript.py').write_text(code,encoding='utf-8')
print('Prepared paragraph mapping and final three-panel caption.')