# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)


# Full 25 x 144 static contact matrix; all 21 valid sidechains; complete profiles.
# Production assets inspected: heatmap/plot_composition.py and PairedBoxScatter1.R.
# They encode category counts and subject before/after inference, respectively,
# so no example data, statistics or inferred biological state are reused.
# Scientific residue tracks inherit the audited v85 native-coordinate composer.
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,FancyArrowPatch
from matplotlib.colors import ListedColormap,BoundaryNorm
from PIL import Image
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;P=OUT/'source_data';R=ROOT/'pymol';F=ROOT/'pymol_context';L=ROOT/'pymol_context_local'
W,H=183.,190.;BLUE='#2166AC';WARM='#B87936';KEY='#3A8D8C';PALE='#B7D7D6';PARTNER='#8B72A5';NGRAY='#8A9299';CGRAY='#BAC1C7'
MODELS=[('Full','pred',BLUE,'o','-'),('Intra','intra',WARM,'s','--')]
mpl.rcParams.update({'savefig.bbox':None,'axes.unicode_minus':False,'axes.labelsize':7.5,'xtick.labelsize':6.8,'ytick.labelsize':6.8})
df=pd.read_csv(P/'all25_dynamics_and_static_lobe_contacts.csv').sort_values('CaMKI_native_position')
ef=pd.read_csv(P/'panel_e_all21_Hchi_grouped_by_type.csv')
cf=pd.read_csv(P/'all3600_CaMKI_CaM_distances.csv')
assert len(df)==25 and df.valid_chi_entropy.sum()==21 and len(ef)==21 and len(cf)==3600
assert cf[['CaMKI_native_position','CaM_native_id']].drop_duplicates().shape[0]==3600
manifest=json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');axes=[];points=[];checks=[]
def txt(x,y,s,fs=8,**kw):return fig.text(x/W,1-y/H,s,fontsize=fs,va=kw.pop('va','top'),color=kw.pop('color',BLACK),**kw)
def ax(x,y,w,h):
 a=fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);axes.append(a);return a
def panel(lab,title,x,y):txt(x,y,lab,11,fontweight='bold');txt(x+6,y,title,8.3,fontweight='bold')
def photo(key,x,y,w,h,trim=False):
 info=manifest['assets'][key];info=info['trim'] if trim else info
 a=ax(x,y,w,h);a.imshow(Image.open(R/info['file']).convert('RGBA'));a.axis('off');return a,info
panel('a','CaM–CaMKI complex',1,2);panel('b','Recognition anchors',94,2)
a,info=photo('overview',1,8,88,39.5,True)
for key,loc in [('N297',(.47,1.03)),('W303',(.00,.67)),('N308',(.91,.83)),('M316',(.97,.20))]:
 a.annotate(key,xy=info['landmarks'][key],xytext=loc,textcoords='axes fraction',ha='center',va='center',fontsize=7.5,color=KEY,arrowprops=dict(arrowstyle='-',color=BLACK,lw=.55,shrinkA=2,shrinkB=1),bbox=dict(fc='white',ec='none',alpha=.85,pad=.5))
for x,label,c in [(6,'N-lobe',NGRAY),(27,'C-lobe',CGRAY),(49,'Linker','#D9DDE0'),(70,'CaMKI',PARTNER)]:
 fig.add_artist(plt.Line2D([(x-4)/W,(x-1)/W],[1-51/H]*2,transform=fig.transFigure,color=c,lw=3));txt(x,49.7,label,7.0)
for key,x,label,pos in [('local_303',94,'W303 · C-lobe',303),('local_316',139,'M316 · N-lobe',316)]:
 txt(x+21,9,label,7.8,ha='center',fontweight='bold');a,inf=photo(key,x,13,42,36)
 landmark=inf['landmarks']['W303' if pos==303 else 'M316'];a.plot(*landmark,marker='o',mfc='none',mec=KEY,ms=15,mew=.7)
for x,(name,prefix,color,marker,style) in zip([132,157],MODELS):
 fig.add_artist(plt.Line2D([(x-4)/W,(x+1)/W],[1-51/H]*2,transform=fig.transFigure,color=color,marker=marker,ms=3,lw=.9,ls=style));txt(x+3,49.5,name,7.7)
panel('d','Lobe contacts',1,57);panel('e','Sidechain environments',65,57);panel('f','Autoinhibitory contacts',123,57)
# d: aggregate every one of the 3600 distances, retaining all 25 peptide rows.
mat=cf.pivot(index='CaMKI_native_position',columns='CaM_native_id',values='minimum_heavy_atom_distance_A').sort_index().sort_index(axis=1)
assert mat.shape==(25,144) and mat.notna().all().all()
D=mat.to_numpy();levels=np.where(D<=4.,2,np.where(D<=4.5,1,0))
counts=pd.DataFrame(index=df.CaMKI_native_position.to_numpy())
for group,col in [('N_lobe','N'),('C_lobe','C'),('linker','Linker')]:
 counts[col]=cf[cf.CaM_lobe==group].groupby('CaMKI_native_position').minimum_heavy_atom_distance_A.apply(lambda x:int((x<=4.5).sum())).reindex(counts.index,fill_value=0)
assert int(counts.to_numpy().sum())==int((D<=4.5).sum())==79
for group,col in [('N_lobe','N'),('C_lobe','C'),('linker','Linker')]:assert np.array_equal(counts[col].to_numpy(),df[group+'_count_4p5A'].to_numpy())
a=ax(14,69,45,52);left=np.zeros(25)
for name,color in [('N',NGRAY),('C',CGRAY),('Linker','#D9DDE0')]:
 a.barh(range(25),counts[name],left=left,height=.65,color=color,edgecolor='#777777',linewidth=.25);left+=counts[name].to_numpy()
a.set(ylim=(24.75,-.75),xlim=(0,10));a.set_xticks([0,3,6,9]);a.set_xlabel('CaM contacts (≤4.5 Å)',labelpad=2,fontsize=7)
a.set_yticks(range(25),[r.amino_acid+str(r.CaMKI_native_position) for _,r in df.iterrows()],fontsize=6.5);a.tick_params(length=2,pad=1.4)
for t in a.get_yticklabels():
 if t.get_text() in ['W303','M316','N297','N308']:t.set_color(KEY);t.set_fontweight('bold')
for x,name,color in [(14,'N',NGRAY),(29,'C',CGRAY),(44,'Linker','#D9DDE0')]:
 fig.add_artist(Rectangle((x/W,1-66/H),2.2/W,1.4/H,transform=fig.transFigure,fc=color,ec='#777777',lw=.25));txt(x+3,64.5,name,6.8)
counts.index.name='CaMKI_native_position';counts.to_csv(P/'all25_lobe_contact_counts_4p5A.csv')
# e: all 21 valid sidechains, grouped only by residue type.
a=ax(79,69,27,52);b=ax(112,69,9,52)
for k,q in ef.iterrows():
 y=k;vals=[q.Full_Hchi,q.Intra_Hchi]
 if k and q.group_order!=ef.iloc[k-1].group_order:
  for z in [a,b]:z.axhline(k-.5,color='#D9DDE0',lw=.45,zorder=0)
 a.plot(vals,[y,y],c='#777777',lw=.6,zorder=1)
 for val,(name,_,color,marker,_) in zip(vals,MODELS):
  a.plot(val,y,marker=marker,c=color,ms=2.8,mew=.45,zorder=3)
  points.append(dict(panel='e',metric='chi_entropy',position=int(q.CaMKI_native_position),model=name,value=float(val)))
 b.barh(y,q.shielded_percent,height=.50,color=KEY,alpha=.75)
a.set(xlim=(0,4.2),ylim=(20.8,-.8));a.set_xticks([0,2,4]);a.set_xlabel('Hχ (nats)',labelpad=2,fontsize=7)
a.set_yticks(range(21),ef.site_label,fontsize=6.7);a.tick_params(length=2,pad=2)
for t in a.get_yticklabels():
 if t.get_text() in ['N297','N308','W303','M316']:t.set_color(KEY);t.set_fontweight('bold')
b.set(xlim=(0,100),ylim=(20.8,-.8));b.set_xticks([0,100],['0','100']);b.tick_params(axis='y',left=False,labelleft=False);b.tick_params(axis='x',length=2,pad=2,labelsize=6.2)
txt(116.5,64,'Shielding',6.7,ha='center');b.set_xlabel('(%)',labelpad=2,fontsize=6.7)
# f: the autoinhibited structure with its adjacent local magnification.
# The 1MXE comparison is already in a/b and is not repeated as another overview.
ctx=json.loads((F/'render_manifest.json').read_text(encoding='utf-8'))
inf=ctx['assets']['autoinhibited']['trim']
txt(152,63,'Autoinhibited CaMKI (1A06)',6.9,ha='center',fontweight='bold')
a=ax(125,67,55,24);a.imshow(Image.open(F/inf['file']).convert('RGBA'));a.axis('off')
for site,loc in [('W303',(.04,.86)),('M316',(.93,.45))]:
 a.annotate(site,xy=inf['landmarks'][site],xytext=loc,textcoords='axes fraction',color=KEY,fontsize=6.5,ha='center',arrowprops=dict(arrowstyle='-',lw=.45,color=BLACK))
# A short connector means magnification, not a biological transition.
fig.add_artist(plt.Line2D([151/W,151/W],[1-91/H,1-94/H],transform=fig.transFigure,color='#777777',lw=.45))
locman=json.loads((L/'render_manifest.json').read_text(encoding='utf-8'));inf=locman['assets']['local']['trim']
a=ax(129,94,47,27);a.imshow(Image.open(L/inf['file']).convert('RGBA'));a.axis('off')
for site,loc in [('L45',(1.09,.79)),('S99',(-.22,.43)),('M316',(1.04,.12))]:
 a.annotate(site,xy=inf['landmarks'][site],xytext=loc,textcoords='axes fraction',color=KEY if site=='M316' else BLACK,fontsize=6.7,ha='center',arrowprops=dict(arrowstyle='-',lw=.45,color=BLACK))
a.text(.83,.53,'4.14 Å',transform=a.transAxes,fontsize=6.5,ha='left',va='center')
a.text(-.22,.20,'3.16 Å',transform=a.transAxes,fontsize=6.5,ha='center',va='center')
# c: shared sequence coordinate, soft structural core, no full-height outline boxes.
panel('c','Complementary flexibility profiles',1,131)
positions=df.CaMKI_native_position.to_numpy();lo,hi=293.5,318.5
x0,pw=19.,158.
coreleft=x0+(302.5-lo)/(hi-lo)*pw;corewidth=14/(hi-lo)*pw
fig.add_artist(Rectangle((coreleft/W,1-141/H),corewidth/W,.8/H,transform=fig.transFigure,fc=KEY,ec='none'))
txt(coreleft+corewidth/2,136.5,'W303–M316 recognition core',7.2,ha='center')
for _,q in df.iterrows():
 x=x0+(q.CaMKI_native_position-lo)/(hi-lo)*pw
 txt(x,143,q.amino_acid,7.0,ha='center',fontweight='bold' if q.CaMKI_native_position in [297,303,308,316] else 'normal',color=KEY if q.CaMKI_native_position in [297,303,308,316] else BLACK)
for metric,y,top,ticks in [('rmsf',149,8.2,[0,4,8]),('chi_entropy',169,4.3,[0,2,4])]:
 a=ax(x0,y,pw,13)
 a.axvspan(302.5,316.5,color='#F0F2F2',ec='none',zorder=0)
 ys=[]
 for name,prefix,color,marker,style in MODELS:
  values=df[prefix+'_'+metric].where(df['valid_'+metric]).to_numpy();ys.append(values)
  a.plot(positions,values,c=color,lw=1,ls=style,marker=marker,ms=2.3,zorder=2)
  for pos in [297,303,308,316]:
   q=df[df.CaMKI_native_position==pos].iloc[0]
   if q['valid_'+metric] and name=='Full':a.plot(pos,q[prefix+'_'+metric],marker=marker,ms=3.5,mfc=color,mec=BLACK,mew=.5,zorder=4)
  points.extend(dict(panel='c',metric=metric,position=int(pos),model=name,value=float(v)) for pos,v in zip(positions,values) if np.isfinite(v))
 a.set(xlim=(lo,hi),ylim=(0,top),ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)');a.set_yticks(ticks);a.tick_params(length=2,pad=2)
 if metric=='rmsf':a.tick_params(axis='x',labelbottom=False)
 else:a.set_xticks([294,297,303,308,312,316,318]);a.set_xlabel('CaMKI residue number',labelpad=2)
 z=np.concatenate([v[np.isfinite(v)] for v in ys]);assert min(z)>=0 and max(z)<=top
 checks.append(dict(metric=metric,ylim=[0,top],range=[float(z.min()),float(z.max())]))
pl=pd.DataFrame(points);assert len(pl[pl.panel=='c'])==92 and len(pl[pl.panel=='e'])==42
pl.to_csv(P/'all_plotted_dynamics.csv',index=False)
pd.DataFrame(levels,index=mat.index,columns=mat.columns).to_csv(P/'contact_class_full25x144.csv')
fig.canvas.draw();ren=fig.canvas.get_renderer();vis=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in vis if t.get_window_extent(ren).x0<-.6 or t.get_window_extent(ren).y0<-.6 or t.get_window_extent(ren).x1>fig.bbox.width+.6 or t.get_window_extent(ren).y1>fig.bbox.height+.6]
assert not outside,outside
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig3_joint_flexibility_compact.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
qa=dict(canvas_mm=[W,H],source='family OOF fold1',panel_layout='a,b top; d,e,f middle; c bottom',all3600_distances_used_in_counts=True,contact_counts={'le4':int((D<=4).sum()),'le4p5':int((D<=4.5).sum())},all21_sidechains_displayed=True,no_MD_drawn_in_main=True,no_numeric_tables=True,no_gray_notes=True,all_25_RMSF_and_21_Hchi_positions=True,core_defined_by_structural_anchors=[303,316],contacts_aggregated_by_lobe=True,all25_contact_rows_retained=True,main_figure_height_reduced_percent=(225-H)/225*100,full_distance_matrix_in_source_data=True,plotted_counts=pl.groupby('panel').size().to_dict(),all_quantitative_axes_start_at_zero=True,axes=checks,full_intra_are_model_variants_not_biological_states=True,model_palette={'Full':BLUE,'Intra':WARM},outside_text=outside,min_font=min(t.get_fontsize() for t in vis),context_render_present=True,visual_review='pending',sha256={ext:hashlib.sha256((OUT/('fig3_joint_flexibility_compact.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False));plt.close(fig)

# Independent complete contact matrix retained as a source figure, not repeated
# in the main plate. All 3600 pairs and both original distance levels are visible.
src=OUT/'source_contact_map';src.mkdir(exist_ok=True)
f2=plt.figure(figsize=(183/25.4,100/25.4),facecolor='white')
a2=f2.add_axes([.085,.16,.89,.72])
a2.imshow(levels,aspect='auto',origin='upper',interpolation='nearest',cmap=ListedColormap(['#FFFFFF',PALE,KEY]),vmin=0,vmax=2,extent=[3.5,147.5,318.5,293.5])
a2.set_yticks(df.CaMKI_native_position,[q.amino_acid+str(q.CaMKI_native_position) for _,q in df.iterrows()],fontsize=7)
a2.set_xticks([4,15,32,51,55,63,71,75,82,92,100,105,116,124,128,136,141,147]);a2.set_xlabel('CaM residue number');a2.set_ylabel('CaMKI residue')
a2.tick_params(length=2,pad=2)
f2.text(.085,.95,'Complete 1MXE residue contact map',fontsize=9,fontweight='bold')
f2.text(.60,.94,'Dark teal: ≤4 Å   Pale teal: 4–4.5 Å',fontsize=7)
for lo2,hi2,color,label in [(3.5,75.5,NGRAY,'N-lobe'),(75.5,81.5,'#D9DDE0',''),(81.5,147.5,CGRAY,'C-lobe')]:
 x=.085+(lo2-3.5)/144*.89;ww=(hi2-lo2)/144*.89
 f2.add_artist(Rectangle((x,.895),ww,.008,transform=f2.transFigure,fc=color,ec='none'))
 if label:f2.text(x+ww/2,.91,label,fontsize=7,ha='center')
for ext,dpi in [('pdf',600),('svg',600),('png',350)]:f2.savefig(src/('all25x144_contacts.'+ext),dpi=dpi,bbox_inches=f2.bbox_inches)
plt.close(f2)
