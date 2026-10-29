from pathlib import Path
import json
v=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v82_shared_view');f=v/'figure3';m=v/'manuscript'
body=(f/'shared_view_body.py').read_text(encoding='utf-8')
start=body.index("panel('b','Local"); end=body.index("panel('c','Absolute",start)
newlocal='''panel('b','Local environments',105,2)
for pos,top in [(72,13),(109,45)]:
    text(106,top,'M'+str(pos),7.7,fontweight='bold',color=TEAL)
    for j,(state,_,c,_) in enumerate(STATES):
        snap('detail_focus_M'+str(pos)+'_'+state,106+j*37,top+5,35,24)
        text(123.5+j*37,top,state,7.3,ha='center',color=c)

'''
body=body[:start]+newlocal+body[end:]
for line in ["    text(x+24.5,63,'2YGG-based' if i==0 else '4LZX-based',7,ha='center',color=MD)\n", "text(50,67.5,'N-lobe helices I/IV aligned · identical view and scale',6.5,ha='center',color=MD)\n", "text(1,73,'Structure:',6.5,color=MD);text(13,73,'N-lobe',6.5,color=BLUE);text(29,73,'C-lobe',6.5,color='#85ABD1');text(45,73,'Partner',6.5,color=WARM);text(61,73,'M72 / M109',6.5,color=TEAL)\n", "text(20,180,'145 backbone / 122 side-chain positions · Pale bands: EF-hand loops · Tracks: interchain Cα ≤ 8 Å',6.5,color=MD)\n", "text(167,197,'122 per state',6.7,ha='center',color=MD)\n", "text(181,242,'Fixed fold 0, training-set illustration · Hχ: nats; RMSF: Å · MD ranges: five consecutive 100 ns windows.',6.3,ha='right',color=MD)\n"]:
    assert line in body,line;body=body.replace(line,'')
body=body.replace("panel('c','Absolute backbone and side-chain profiles'","panel('c','Backbone and side-chain flexibility'")
body=body.replace("panel('a','Matched reference view'","panel('a','Recognition structures'")
body=body.replace("panel('f','Absolute Hχ'","panel('f','Hχ prediction'")
body=body.replace("W,H=183.,245.","W,H=183.,235.")
body=body.replace("panel('d','Contact groups',1,189);panel('e','All nine Met',78,189);panel('f','Hχ prediction',143,189)","panel('d','Contact groups',1,185);panel('e','All nine Met',78,185);panel('f','Hχ prediction',143,185)")
body=body.replace('text(x+13.5,197,','text(x+13.5,193,').replace('text(x+11,197,','text(x+11,193,').replace('axis(x,202,27,24)','axis(x,198,27,24)').replace('axis(x,202,22,24)','axis(x,198,22,24)').replace('axis(153,202,27,24)','axis(153,198,27,24)')
body=body.replace("text(111,237,'Mature CaM Met position',6.4,color=MD,ha='center')","text(112,232,'Mature CaM Met position',6.8,color=BLACK,ha='center')")
body=body.replace("6.9,color=MD)","7.2,color=BLACK)").replace("7.1,ha='center',color=MD)","7.1,ha='center',color=BLACK)")
body=body.replace("local_views_right_of_whole=True,","local_views_right_of_whole=True,no_numeric_tables=True,no_gray_note_text=True,")
body=body.replace("for ext,dpi in", "assert not any(t.get_color()==MD for t in visible), 'Gray explanatory text must be in caption, not figure'\nfor ext,dpi in")
# Keep a compact, necessary legend with the whole structures.
anchor="panel('b','Local environments',105,2)"
body=body.replace(anchor,"text(10,68,'CaM N-lobe',7.4,color=BLUE);text(43,68,'C-lobe',7.4,color='#85ABD1');text(67,68,'Partner',7.4,color=WARM)\n\n"+anchor)
(f/'shared_view_body.py').write_text(body,encoding='utf-8')
s=(f/'compose_fig3.py').read_text(encoding='utf-8');header=s[:s.index('# --- panel code ---')].replace('existing renderer + absolute values -> native run / param inherit','existing renderer -> native run')
(f/'compose_fig3.py').write_text(header+'# --- panel code ---\n'+body,encoding='utf-8')
caption=(f/'caption_zh.md').read_text(encoding='utf-8').replace('M72、M109局部环境及两态绝对读数','M72、M109局部环境').replace('两态输入结构，显示时','NHE1与IQCG态分别源自2YGG、4LZX；显示时')
(f/'caption_zh.md').write_text(caption,encoding='utf-8')
plan=json.loads((m/'v82_text_plan.json').read_text(encoding='utf-8'))
old=plan['replacements']['74']
new='全部九个Met提供了相同化学类型的具体对照（图3b、e）。M72的MD RMSF由1.170降至0.759 Å，Hχ由1.644升至2.018 nats；Full分别为0.877→0.810 Å与1.862→2.112 nats。M109的MD RMSF同样由2.360降至0.724 Å，Hχ却由1.627降至1.466 nats；Full分别为1.396→0.715 Å与1.971→1.759 nats。两个位置呈现主链波动下降而侧链状态变化方向不同的对照。全部九个Met的侧链变化排序相关为0.883，五个绝对变化至少0.1 nats的位置均判对方向；方向或幅度不准确的M36、M124也保留展示。这里Hχ描述固定角区下多角联合状态的占据多样性，不能直接解释为转动更快或连续角度范围必然更宽。'
plan['replacements']['74']=new
plan['notes'].append('按用户要求删除图内数值表及全部灰色说明小字，关键数值移至3.3正文，方法与统计范围移图注。')
(m/'v82_text_plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'old_paragraph74_length':len(old),'new_paragraph74_length':len(new),'caption_characters':len(caption.strip()),'no_tables':True,'no_gray_notes':True},ensure_ascii=False))