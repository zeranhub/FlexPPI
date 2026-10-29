from pathlib import Path
s=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025')
p=s/'revision_v88_three_panel_story/figure3/compose_fig3_1MXE.py';code=p.read_text(encoding='utf-8')
code=code.replace('# Three panels only: structural location, M316 static contexts, model profiles.','# Two panels only: structural location and complete model profiles.')
code=code.replace(";B=ROOT/'pymol_b_context'",'').replace('W,H=183.,155.','W,H=183.,100.')
a=code.index(";ctx=json.loads(");b=code.index('fig=plt.figure(',a);code=code[:a]+'\n'+code[b:]
code=code.replace("panel('a','CaM–CaMKI complex',1,2);panel('b','M316 contact environments',87,2)","panel('a','Structural location',1,2);panel('b','Predicted flexibility',76,2)")
code=code.replace('a=ax(1,9,80,49)','a=ax(1,19,67,54)').replace("('W303',(.01,.67))","('W303',(.08,.67))")
a=code.index('for x,label,c in [');b=code.index('# One global numerical legend;',a)
code=code[:a]+"""for x,label,c in [(6,'CaM N-lobe',NGRAY),(34,'C-lobe',CGRAY),(57,'CaMKI',PARTNER)]:
 fig.add_artist(plt.Line2D([(x-4)/W,(x-1)/W],[1-82/H]*2,transform=fig.transFigure,color=c,lw=3));txt(x,80.6,label,7.2)
"""+code[b:]
code=code.replace('1-75/H','1-4/H').replace('73.5,name','2.5,name')
code=code.replace("panel('c','Predicted backbone and sidechain flexibility',1,73)\n",'')
code=code.replace('x0,pw=19.,158.','x0,pw=84.,94.')
code=code.replace('1-86/H','1-20/H').replace('81.5,','15.5,').replace('txt(x,89,','txt(x,23,')
code=code.replace("('rmsf',96,8.2","('rmsf',31,8.2").replace("('chi_entropy',124,4.3","('chi_entropy',64,4.3").replace('pw,20)','pw,24)')
code=code.replace("(297,(298.2,2.65))","(297,(299.0,2.65))").replace("(308,(308.8,1.80))","(308,(310.0,1.85))")
code=code.replace("panel='c'","panel='b'").replace('fig3_three_panel_story.','fig3_two_panel_story.')
code=code.replace("panels=['a','b','c']","panels=['a','b']").replace("removed_independent_panels=['d','e','f']","removed_independent_panels=['b_static_context','d','e','f']")
code=code.replace("structural_both_rotation_and_scale_matched=True,b_is_static_functional_context_not_model_validation=True,","static_context_retired_to_source_files=True,")
(s/'revision_v88_two_panel_story/figure3/compose_fig3_1MXE.py').write_text(code,encoding='utf-8')