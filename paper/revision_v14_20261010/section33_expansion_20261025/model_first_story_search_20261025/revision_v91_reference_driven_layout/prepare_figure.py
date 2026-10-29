from pathlib import Path
import shutil
S=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');V=S/'revision_v91_reference_driven_layout';O=S/'revision_v90_data_led_story'
(V/'figure3/source_data').mkdir(parents=True,exist_ok=True)
shutil.copy2(O/'figure3/source_data/all25_dynamics_and_static_lobe_contacts.csv',V/'figure3/source_data')
head=(O/'figure3/compose_fig3_1MXE.py').read_text(encoding='utf-8').split('# Inherited audited full-profile plotting;')[0]
body=r'''
# Reference-driven: true peptide geometry carries Full predictions, aligned complete profiles remain quantitative evidence.
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;P=OUT/'source_data';R=ROOT/'pymol_quantitative'
W,H=183.,146.;BLUE='#2166AC';WARM='#B87936';KEY='#3A8D8C';INK='#222222'
MODELS=[('Full','pred',BLUE,'o','-'),('Intra','intra',WARM,'s','--')]
mpl.rcParams.update({'savefig.bbox':None,'axes.unicode_minus':False,'axes.labelsize':8.5,'xtick.labelsize':8.2,'ytick.labelsize':8.2})
df=pd.read_csv(P/'all25_dynamics_and_static_lobe_contacts.csv').sort_values('CaMKI_native_position');assert len(df)==25 and df.valid_chi_entropy.sum()==21
manifest=json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');points=[];checks=[]
def txt(x,y,s,fs=8.5,**kw):return fig.text(x/W,1-y/H,s,fontsize=fs,va=kw.pop('va','top'),color=kw.pop('color',INK),**kw)
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def panel(lab,title,x,y):txt(x,y,lab,11,fontweight='bold');txt(x+6,y,title,9,fontweight='bold')
def line(x1,y1,x2,y2,**kw):fig.add_artist(plt.Line2D([x1/W,x2/W],[1-y1/H,1-y2/H],transform=fig.transFigure,**kw))
panel('a','Backbone RMSF · Full',1,2);panel('b','Side-chain Hχ · Full',96,2)
# Both maps show identical coordinates, camera, optical scale and crop. No predicted motion/ensemble is invented.
for metric,x0,title,vmax,unit in [('RMSF',3,'RMSF',2.3,'Å'),('Hchi',98,'Hχ',4.3,'nats')]:
 info=manifest['delivery_assets'][metric];a=ax(x0,11,81,41.6);a.imshow(Image.open(R/info['focus_file']).convert('RGBA'));a.axis('off')
 labels=[('W303',(.34,.97)),('M316',(.89,.12))] if metric=='RMSF' else [('N297',(.09,.18)),('N308',(.61,.87))]
 for key,loc in labels:
  xy=info['focus_named_landmarks'][key]['CA' if metric=='RMSF' else 'sidechain_center']
  a.annotate(key,xy=xy,xytext=loc,textcoords='axes fraction',ha='center',va='center',fontsize=8.3,color=INK,arrowprops=dict(arrowstyle='-',lw=.55,color=INK),bbox=dict(fc='white',ec='none',alpha=.95,pad=.3))
 cbax=ax(x0+18,57,45,1.9);cb=fig.colorbar(ScalarMappable(norm=Normalize(0,vmax),cmap='cividis'),cax=cbax,orientation='horizontal',ticks=[0,1,2] if metric=='RMSF' else [0,2,4]);cb.outline.set_visible(False);cb.ax.tick_params(length=2,pad=2,labelsize=8.2)
 txt(x0+66,56.4,unit,8.2)
# The lower panel uses one residue coordinate for sequence, structural anchor span and both original model readouts.
panel('c','Complete residue profiles',1,70)
for x,(name,prefix,color,marker,style) in zip([139,162],MODELS):
 line(x-4,72,x+1,72,color=color,marker=marker,ms=3.5,lw=1.1,ls=style);txt(x+3,70.4,name,8.5)
positions=df.CaMKI_native_position.to_numpy();lo,hi=293.5,318.5;x0,pw=18.,160.
coreleft=x0+(302.5-lo)/(hi-lo)*pw;corewidth=14/(hi-lo)*pw
line(coreleft,79,coreleft+corewidth,79,color=KEY,lw=1.8)
txt(coreleft+corewidth/2,75.6,'W303–M316 core',8.2,ha='center',color=KEY)
# Complete sequence, once, shared by both quantitative profiles. A grey Ala point is never substituted for undefined entropy.
for _,q in df.iterrows():
 xx=x0+(q.CaMKI_native_position-lo)/(hi-lo)*pw
 txt(xx,83,q.amino_acid,8.3,ha='center',fontweight='bold' if q.CaMKI_native_position in [297,303,308,316] else 'normal')
for metric,y,ph,top,ticks in [('rmsf',90,17,8.2,[0,4,8]),('chi_entropy',116,21,4.3,[0,2,4])]:
 a=ax(x0,y,pw,ph);a.axvspan(302.5,316.5,color='#F0F4F3',ec='none',zorder=0);ys=[]
 for name,prefix,color,marker,style in MODELS:
  values=df[prefix+'_'+metric].where(df['valid_'+metric]).to_numpy();ys.append(values)
  a.plot(positions,values,c=color,lw=1.15,ls=style,marker=marker,ms=3.1,zorder=2)
  for pos in ([303,316] if metric=='rmsf' else [297,308]):
   q=df[df.CaMKI_native_position==pos].iloc[0]
   if name=='Full':a.plot(pos,q[prefix+'_'+metric],marker=marker,ms=4.8,mfc=color,mec=INK,mew=.6,zorder=4)
  points.extend(dict(panel='c',metric=metric,position=int(pos),model=name,value=float(v)) for pos,v in zip(positions,values) if np.isfinite(v))
 a.set(xlim=(lo,hi),ylim=(0,top),ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)');a.set_yticks(ticks);a.tick_params(length=2,pad=2)
 a.set_xticks([294,297,303,308,312,316,318])
 if metric=='rmsf':
  a.tick_params(axis='x',labelbottom=False)
  for pos,loc in [(303,(304.0,1.9)),(316,(315.0,2.5))]:
   q=df[df.CaMKI_native_position==pos].iloc[0];a.annotate(q.amino_acid+str(pos),xy=(pos,q.pred_rmsf),xytext=loc,fontsize=8.2,color=INK,ha='center',arrowprops=dict(arrowstyle='-',lw=.5,color=INK))
 else:
  a.set_xlabel('CaMKI residue number',labelpad=2)
  for pos,loc in [(297,(299.2,2.0)),(308,(309.7,1.65))]:
   q=df[df.CaMKI_native_position==pos].iloc[0];a.annotate(q.amino_acid+str(pos),xy=(pos,q.pred_chi_entropy),xytext=loc,fontsize=8.2,color=INK,ha='center',arrowprops=dict(arrowstyle='-',lw=.5,color=INK))
 z=np.concatenate([v[np.isfinite(v)] for v in ys]);assert min(z)>=0 and max(z)<=top;checks.append(dict(metric=metric,ylim=[0,top],range=[float(z.min()),float(z.max())]))
pl=pd.DataFrame(points);assert len(pl)==92;pl.to_csv(P/'all_plotted_dynamics.csv',index=False)
for text_obj in fig.findobj(mpl.text.Text):
 if text_obj.get_visible() and text_obj.get_text() and text_obj.get_fontsize()<8.2:text_obj.set_fontsize(8.2)
fig.canvas.draw();ren=fig.canvas.get_renderer();vis=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in vis if t.get_window_extent(ren).x0<-.6 or t.get_window_extent(ren).y0<-.6 or t.get_window_extent(ren).x1>fig.bbox.width+.6 or t.get_window_extent(ren).y1>fig.bbox.height+.6];assert not outside,outside
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig3_spatial_sequence_study.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
qa=dict(canvas_mm=[W,H],source='family OOF fold1',panels=['a Full RMSF in actual geometry','b Full Hchi in identical geometry','c complete absolute profiles'],all25_RMSF_and21_Hchi_positions=True,unique_plotted_profile_values=len(pl),repeated_zoom_removed=True,all_quantitative_axes_start_at_zero=True,axes=checks,no_MD_in_main=True,no_numeric_tables=True,no_gray_explanatory_notes=True,no_mechanism_schematic=True,activity_not_measured=True,Intra_not_unbound=True,structures_are_original_coordinates=True,structure_metric_colorbar_units_separate=True,caM_context_cropped_not_full_overview=True,outside_text=outside,min_font=min(t.get_fontsize() for t in vis),visual_review='pending',sha256={ext:hashlib.sha256((OUT/('fig3_spatial_sequence_study.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False));plt.close(fig)
'''
(V/'figure3/compose_fig3_1MXE.py').write_text(head+body,encoding='utf-8')