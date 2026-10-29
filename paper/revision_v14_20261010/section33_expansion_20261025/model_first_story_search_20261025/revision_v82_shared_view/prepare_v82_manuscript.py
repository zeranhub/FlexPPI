from pathlib import Path
import json
root=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');old=root/'revision_v81_absolute_profiles';v=root/'revision_v82_shared_view';m=v/'manuscript';f=v/'figure3'
caption='图3：CaM两种识别构型的主链与侧链柔性。（a）两态输入结构，以N叶瓣螺旋I/IV（7–19、65–75号Cα）对整个复合物作刚体叠合，采用相同观察方向、比例及裁框，保留CaM叶间构型差异；结构颜色见图例。（b）紧邻整体结构的M72、M109局部环境及两态绝对读数；各位点局部分别按前后五个残基的主链叠合，同一位点两态使用相同视角与比例。（c）145个共同主链位置的RMSF及122个有效侧链位置的Hχ，暖赭、蓝色分别表示NHE1、IQCG态；MD与Full上下对应，同一指标共用刻度。浅带为EF手环区，底部标记跨链Cα距离≤8 Å的接触；无有效Hχ处留空。（d）完整四组接触侧链的两态平均Hχ。（e）全部九个Met的Hχ。（f）各态全部122个侧链的MD与预测Hχ对应。RMSF单位Å，Hχ单位nats；d、e细线为同一轨迹五个连续100 ns段的最小至最大值，并非独立重复的置信区间。两条记录均用于固定第0折训练，本图用于结构解释。'
(f/'caption_zh.md').write_text(caption+'\n',encoding='utf-8')
plan=json.loads((old/'manuscript/v81_text_plan.json').read_text(encoding='utf-8'));plan['output_docx']='E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v82_20261025_中文修订_图3统一视角与局部并排.docx'
plan['replacements']['69']=plan['replacements']['69'].replace('两条记录，将成熟CaM序列4–148号的145个共同位置逐位对应','两条记录。两态的CaM序列相同，但叶瓣及伙伴的空间排列均有差异；将成熟CaM序列4–148号的145个共同位置逐位对应')
plan['replacements']['70']=plan['replacements']['70'].replace('（图3b）','（图3c）')
plan['replacements']['73']=plan['replacements']['73'].replace('（图3c）','（图3d）')
plan['replacements']['74']=plan['replacements']['74'].replace('（图3d、e）','（图3b、e）')
plan['notes']=['在v81基础上只更新图3结构参照方向、局部布局与panel引用，定量原值不变。','整体图使用稳定N叶螺旋的一个全复合物刚体变换，不能将同序列CaM误述为同构象。','原e的局部放大移到整体a右边并重编号b；原a散点移为f，全部点保留。']
(m/'v82_text_plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
s=(old/'manuscript/build_v81_manuscript.py').read_text(encoding='utf-8').replace('v81_text_plan.json','v82_text_plan.json').replace('v81_source_audit.json','v82_source_audit.json').replace('v81_build_audit.json','v82_build_audit.json').replace('v81_final_text.txt','v82_final_text.txt')
s=s.replace('CaM NHE1 and IQCG training-set structural interpretation with complete residue profiles','CaM matched rigid-reference orientation with right-hand local views and absolute two-state profiles')
(m/'build_v82_manuscript.py').write_text(s,encoding='utf-8')
# Widen the four-category axes to avoid touching category labels.
for name in ['shared_view_body.py','compose_fig3.py']:
    p=f/name;s=p.read_text(encoding='utf-8').replace("[('md',13),('full0',47)]","[('md',13),('full0',49)]").replace('text(x+12,197,','text(x+13.5,197,').replace('axis(x,202,24,24)','axis(x,202,27,24)').replace("[('md',85),('full0',116)]","[('md',88),('full0',118)]").replace('axis(x,202,23,24)','axis(x,202,22,24)').replace('text(x+11.5,197,','text(x+11,197,');p.write_text(s,encoding='utf-8')
print(json.dumps({'caption_characters':len(caption),'output':plan['output_docx']},ensure_ascii=False))