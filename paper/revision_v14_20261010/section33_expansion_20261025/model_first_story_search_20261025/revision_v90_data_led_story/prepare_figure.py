from pathlib import Path
import shutil,json
S=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');O=S/'revision_v89_mechanism_story';V=S/'revision_v90_data_led_story'
(V/'figure3/source_data').mkdir(parents=True,exist_ok=True);(V/'narrative').mkdir(exist_ok=True)
shutil.copy2(O/'figure3/source_data/all25_dynamics_and_static_lobe_contacts.csv',V/'figure3/source_data')
shutil.copytree(O/'pymol_mechanism',V/'pymol',dirs_exist_ok=True)
old=(O/'figure3/compose_fig3_1MXE.py').read_text(encoding='utf-8')
head=old.split('# Inherited audited full-profile plotting;')[0]
body=old.split('# Inherited audited full-profile plotting;')[1].split('# A functional, falsifiable interpretation.')[0]
body=body.replace('native vector mechanism schematic.','data-led layout with no mechanism schematic.').replace("R=ROOT/'pymol_mechanism'","R=ROOT/'pymol'").replace('W,H=183.,178.','W,H=183.,154.').replace("'From inhibitory contacts to CaM capture'","'Structural context of the regulatory segment'")
body=body.replace("from matplotlib.patches import Rectangle,Ellipse,FancyArrowPatch,PathPatch,Circle","from matplotlib.patches import Rectangle").replace('from matplotlib.path import Path as MPath\n','')
start=body.index('def arr(');end=body.index("panel('a'",start);body=body[:start]+body[end:]
body=body.replace("a=ax(x0,93,pw,29.5)","a=ax(x0,93,pw,16.5 if metric=='rmsf' else 50.5)")
body=body.replace("a.set_xticks([294,298,303,308,313,318]);a.set_xlabel('CaMKI residue number',labelpad=2)","a.set_xticks([294,298,303,308,313,318])\n if metric=='rmsf':a.tick_params(axis='x',labelbottom=False);a.set_yticks([0,4,8])\n else:a.set_xlabel('CaMKI residue number',labelpad=2)")
extra=r'''
# Replot all25 Full RMSF values on a magnified absolute y scale; not independent evidence.
txt(17,113,'Full: expanded y scale',7.9,color=BLUE)
a=ax(17,119,69,24.5);a.axvspan(302.5,316.5,color='#EEF4F3',ec='none',zorder=0)
a.plot(positions,df.pred_rmsf,c=BLUE,lw=1.1,marker='o',ms=2.7,zorder=2)
for pos,loc in [(296,(298.4,1.75)),(303,(303.7,1.05)),(316,(313.3,1.3))]:
 q=df[df.CaMKI_native_position==pos].iloc[0];a.plot(pos,q.pred_rmsf,'o',ms=4.2,mfc=BLUE,mec=INK,mew=.5,zorder=4)
 a.annotate(q.amino_acid+str(pos),xy=(pos,q.pred_rmsf),xytext=loc,fontsize=7.9,color=KEY,ha='center',arrowprops=dict(arrowstyle='-',lw=.5,color=INK))
a.set(xlim=(lo,hi),ylim=(0,2.3),ylabel='RMSF (Å)');a.set_yticks([0,1,2]);a.set_xticks([294,298,303,308,313,318]);a.set_xlabel('CaMKI residue number',labelpad=2);a.tick_params(length=2,pad=2)
# Frozen original values remain unchanged, including the repeated Full track.
pl=pd.DataFrame(points);assert len(pl)==92;pl.to_csv(P/'all_plotted_dynamics.csv',index=False)
zoom=pd.DataFrame({'position':positions,'model':'Full','metric':'rmsf','value':df.pred_rmsf,'display':'expanded_y_scale'});zoom.to_csv(P/'repeated_Full_RMSF_expanded.csv',index=False);assert len(zoom)==25
for text_obj in fig.findobj(mpl.text.Text):
 if text_obj.get_visible() and text_obj.get_text() and text_obj.get_fontsize()<7.9:text_obj.set_fontsize(7.9)
fig.canvas.draw();ren=fig.canvas.get_renderer();vis=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in vis if t.get_window_extent(ren).x0<-.6 or t.get_window_extent(ren).y0<-.6 or t.get_window_extent(ren).x1>fig.bbox.width+.6 or t.get_window_extent(ren).y1>fig.bbox.height+.6];assert not outside,outside
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig3_data_led_story.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
qa=dict(canvas_mm=[W,H],source='family OOF fold1',panels=['a','b'],all25_RMSF_and21_Hchi_positions=True,unique_plotted_profile_values=len(pl),repeated_Full_RMSF_values=len(zoom),repeated_data_not_independent_evidence=True,all_quantitative_axes_start_at_zero=True,original_axes=checks,expanded_RMSF_ylim=[0,2.3],no_MD_in_main=True,no_numeric_tables=True,no_gray_explanatory_notes=True,no_mechanism_schematic=True,activity_not_measured=True,Intra_not_unbound=True,outside_text=outside,min_font=min(t.get_fontsize() for t in vis),visual_review='pending',sha256={ext:hashlib.sha256((OUT/('fig3_data_led_story.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False));plt.close(fig)
'''
(V/'figure3/compose_fig3_1MXE.py').write_text(head+'# Inherited audited full-profile plotting;'+body+extra,encoding='utf-8')
(V/'figure_contract.md').write_text('''# V90 图3

撤去V89的c机制示意。a保留真实1A06/1MXE结构与F298局部接触；b保留全部Full/Intra绝对曲线，并在左下仅重绘同一组Full的全部25位RMSF、纵轴0–2.3Å，以看清低值区域的相对变化。这不是第三份证据。所有原始92个唯一读数与冻结CSV相同，另重复25个Full读数。机制假说和检验建议仅在正文叙述。没有新推理/训练/MD，也不将Intra画成游离肽。
''',encoding='utf-8')
print(V)