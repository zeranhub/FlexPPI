from pathlib import Path
import json
root=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');old=root/'revision_v81_absolute_profiles';out=root/'revision_v82_shared_view/figure3'
src=(old/'figure3/absolute_body.py').read_text(encoding='utf-8'); prefix=src[:src.index("panel('a'")]
newtop='''panel('a','Matched reference view',1,2)
for i,(state,_,_,_) in enumerate(STATES):
    x=1+i*47;text(x+22,10,state+'–CaM',8.2,ha='center',fontweight='bold')
    snap('whole_'+state,x,14,44,47)
    text(x+22,63,'2YGG-based' if i==0 else '4LZX-based',7,ha='center',color=MD)
text(46,67.5,'N-lobe helices I/IV aligned · identical view and scale',6.5,ha='center',color=MD)
text(1,73,'Structure:',6.5,color=MD);text(13,73,'N-lobe',6.5,color=BLUE);text(29,73,'C-lobe',6.5,color='#85ABD1');text(45,73,'Partner',6.5,color=WARM);text(61,73,'M72 / M109',6.5,color=TEAL)

panel('b','Local views and absolute readouts',98,2)
for pos,top in [(72,10),(109,42)]:
    text(99,top,'M'+str(pos)+(' · outside both contact sets' if pos==72 else ' · IQCG contact'),7.2,fontweight='bold',color=TEAL)
    for j,(state,_,c,_) in enumerate(STATES):
        snap('detail_focus_M'+str(pos)+'_'+state,98+j*21,top+10,20,18)
        text(108+j*21,top+6,state,6.6,ha='center',color=c)
    text(164,top+6,'NHE1',6.5,ha='center',color=WARM);text(177,top+6,'IQCG',6.5,ha='center',color=BLUE)
    q=chi[chi.mature_cam_number==pos].iloc[0];b=back[back.mature_cam_number==pos].iloc[0]
    for k,(label,row,kind) in enumerate([('Hχ MD',q,'md'),('Hχ Full',q,'full0'),('RMSF MD',b,'md'),('RMSF Full',b,'full0')]):
        y=top+12+k*4;text(141,y,label,6.4)
        for state,s,c,marker in STATES:
            text(164 if s=='a' else 177,y,f'{row[kind+"_"+s]:.3f}',6.4,ha='center')
            points.append(dict(panel='b',position=pos,state=state,metric='hchi' if label.startswith('Hχ') else 'rmsf',kind=kind,value=float(row[kind+'_'+s])))
text(140.5,73,'Hχ: nats; RMSF: Å',6.5,color=MD,ha='center')

'''
middle=src[src.index("panel('b','Absolute"):src.index("panel('c','Complete")]
middle=middle.replace("panel('b','Absolute backbone and side-chain profiles',1,65)","panel('c','Absolute backbone and side-chain profiles',1,80)").replace("text(127,65,","text(127,80,").replace("text(153,65,","text(153,80,")
middle=middle.replace('axis(20,74,157,3)','axis(20,88,157,3)').replace("('rmsf','md',81),('rmsf','full0',97),('hchi','md',116),('hchi','full0',132)","('rmsf','md',95),('rmsf','full0',111),('hchi','md',130),('hchi','full0',146)").replace("panel='b'","panel='c'").replace('axis(20,148,157,5)','axis(20,162,157,5)').replace('text(20,161,','text(20,180,')
bottom=src[src.index("panel('c','Complete"):src.index("panel('e','Local")]
bottom=bottom.replace("panel('c','Complete contact groups',1,169);panel('d','All nine methionines',96,169)","panel('d','Contact groups',1,189);panel('e','All nine Met',78,189);panel('f','Absolute Hχ',143,189)")
bottom=bottom.replace("[('md',15),('full0',56)]","[('md',13),('full0',47)]").replace("text(x+15.5,176,","text(x+12,197,").replace('axis(x,181,31,20)','axis(x,202,24,24)').replace("panel='c'","panel='d'")
bottom=bottom.replace("[('md',108),('full0',149)]","[('md',85),('full0',116)]").replace('text(x+14.5,176,','text(x+11.5,197,').replace('axis(x,181,29,20)','axis(x,202,23,24)').replace("panel='d',position","panel='e',position").replace('fontsize=6.1,rotation=60','fontsize=6.1,rotation=60').replace("text(106,211,'Mature CaM Met position',6.7,color=MD)","text(111,237,'Mature CaM Met position',6.4,color=MD,ha='center')")
scatter='''
a=axis(153,202,27,24);a.plot([0,4.1],[0,4.1],color='#B3B8BC',ls='--',lw=.6)
for state,s,c,marker in STATES:
    a.scatter(chi['md_'+s],chi['full0_'+s],s=6,marker=marker,color=c,alpha=.58,edgecolors='white',linewidths=.18)
    points.extend(dict(panel='f',position=int(q.mature_cam_number),state=state,MD=q['md_'+s],Full=q['full0_'+s]) for _,q in chi.iterrows())
a.set(xlim=(0,4.1),ylim=(0,4.1),xlabel='MD Hχ (nats)',ylabel='Predicted Hχ (nats)');a.set_xticks([0,2,4]);a.set_yticks([0,2,4]);a.tick_params(pad=2,length=2,labelsize=6.5);a.set_xlabel('MD Hχ (nats)',fontsize=6.6,labelpad=2);a.set_ylabel('Predicted Hχ (nats)',fontsize=6.6,labelpad=2)
text(167,197,'122 per state',6.7,ha='center',color=MD)
text(181,242,'Fixed fold 0, training-set illustration · Hχ: nats; RMSF: Å · MD ranges: five consecutive 100 ns windows.',6.3,ha='right',color=MD)
'''
end=src[src.index('fig.canvas.draw();'):].replace('fig3_CaM_absolute_profiles.','fig3_CaM_shared_view.')
end=end.replace("two_state_profiles=True,","two_state_profiles=True,local_views_right_of_whole=True,structural_alignment='N-lobe helices I/IV; one whole-complex rigid transform',")
body=prefix+newtop+middle+bottom+scatter+end
(out/'shared_view_body.py').write_text(body,encoding='utf-8')
oldcompose=(old/'figure3/compose_fig3.py').read_text(encoding='utf-8');header=oldcompose[:oldcompose.index('# --- panel code ---')]
header=header[:header.index('# Academic Figure Skill Typography')]+header[header.index('# Academic Figure Skill Typography'):]
header=header.replace('# (a) real PyMOL inputs + 244 absolute-state points -> existing real-structure renderer and correlation style -> param inherit','# (a) true coordinates, shared rigid-reference PyMOL view -> existing renderer -> native run')
header=header.replace('# (b) shared-sequence profiles + contact tracks -> LineTrend/plot_trend.py and heatmap/plot_comparison.py -> param inherit','# (b) paired local PyMOL views next to whole structures -> existing renderer + absolute values -> native run / param inherit')
header=header.replace('# (c) group means + all5 continuous windows -> GroupedBarChart/plot_GroupedBarChartv1.py -> param inherit','# (c) shared-sequence profiles + contact tracks -> LineTrend/plot_trend.py -> param inherit')
header=header.replace('# (d) all9 Met absolute states + all5 windows -> same grouped-point visual logic -> param inherit','# (d) group means + all5 continuous windows -> existing paired-point logic -> param inherit')
header=header.replace('# (e) actual complete-chain PyMOL local contexts + paired absolute readouts -> image_plate + quant -> param inherit','# (e) all9 Met absolute states + all5 windows -> existing paired-line logic -> param inherit\n# (f) all244 absolute-state MD / Full points -> existing correlation style -> param inherit')
(out/'compose_fig3.py').write_text(header+'# --- panel code ---\n'+body,encoding='utf-8')
print('Prepared shared-view figure composer; using old structures provisionally until final new PyMOL renders arrive.')