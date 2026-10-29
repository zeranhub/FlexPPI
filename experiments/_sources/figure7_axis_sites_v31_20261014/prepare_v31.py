from pathlib import Path
ROOT=Path(__file__).resolve().parent;FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
old=ROOT.parent/'rich_angle_figure_v30_20261014/figures/plot_full_residue_and_chi_occupancy.py'
s=old.read_text(encoding='utf8')
s=s.replace('# (a) sequence profiles + aligned errors/angle widths -> LineTrend/plot_trend.py and heatmap/plot_composition.py -> param inherit','# (a) complete66 Hchi sequence profile + six GLN axis labels -> LineTrend/plot_trend.py -> param inherit')
s=s.replace('# Retain all66 valid positions and all151 chemically defined angles; no smoothing, sampling or imputed missing angles.','# Retain all66 valid Hchi positions and all6GLN18angles; no smoothing, sampling or imputed missing angles.')
s=s.replace("OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;DATA=ROOT/'experiment'","OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;DATA=ROOT.parent/'rich_angle_figure_v30_20261014/experiment'\nif not DATA.exists():DATA=ROOT/'source'")
s=s.replace("metrics=pd.read_csv(DATA/'all66_per_angle_metrics.csv')\n",'')
s=s.replace("STEM='fig7_full_residue_and_chi_occupancy'","STEM='fig7_Hchi_profile_six_Gln_axis_and_occupancy'")
s=s.replace('W,H=183.,214.','W,H=183.,184.')
start=s.index("text(3,207,'a'");end=s.index("text(3,116,'b'")
s=s[:start]+r'''text(3,177,'a',fontsize=11,weight='bold')
text(20,177,'完整泛素：联合熵的模型预测与 MD 参考',fontsize=8.5,weight='bold')
text(20,169,'全部 66 个有效位置   MAE 0.080 nats   ρ 0.979   类型均值校正 ρ 0.868',fontsize=7.5)
pa=axis(20,135,150,27)
ids=np.arange(1,73);model=np.full(72,np.nan);reference=model.copy()
for i,r in rows.iterrows():
    pos=int(r.author_residue_id)-1;model[pos]=r.predicted_Hchi_nats;reference[pos]=r.frozen_MD_joint_Hchi_nats
    for kind,y in [('model_prediction',model[pos]),('MD_reference',reference[pos])]:points.append(dict(panel='a_Hchi',residue_id=pos+1,kind=kind,x=pos+1,y=y))
pa.set_xlim(.5,72.5)
for site,c in [(41,ORANGE),(2,PURPLE)]:pa.axvspan(site-.5,site+.5,color=c,alpha=.16,lw=0,zorder=0)
pa.plot(ids,reference,color=BLACK,lw=.7,ls=':',marker='o',markersize=2.4,markerfacecolor='white',markeredgewidth=.55,label='MD 参考',zorder=2)
pa.plot(ids,model,color=BLUE,lw=.9,marker='s',markersize=2.2,markeredgewidth=0,label='模型预测',zorder=3)
pa.set_ylim(-.15,4.25);pa.set_yticks([0,2,4]);pa.set_ylabel('联合 Hχ\n(nats)')
pa.set_xticks([1,10,20,30,50,60,72])
pa.legend(loc='upper right',bbox_to_anchor=(1.0,1.34),ncol=2,fontsize=7.5,handlelength=1.8,columnspacing=1.4)
axis_labels=[]
for site in [2,31,40,41,49,62]:
    c=identity(site)
    pa.plot([site,site],[0,-.028],transform=pa.get_xaxis_transform(),color=c,lw=1.0,clip_on=False,zorder=5)
    y=-.48 if site==41 else -.33
    t=pa.text(site,y,'Q'+str(site),transform=pa.get_xaxis_transform(),ha='center',va='center',fontsize=7.5,color=c,clip_on=False)
    axis_labels.append(dict(residue_id=site,label='Q'+str(site),x_data=float(site),y_axes=y))
text(20,121,'泛素残基编号',fontsize=7.5)

'''+s[end:]
s=s.replace('assert len(points)==2202','assert len(points)==1872')
s=s.replace('（a）泛素全部66个Hχ有效残基，按编号对齐联合熵的预测与MD参考、预测减参考的误差，以及χ1–χ4满足F*≤1的MD角区总宽度。宽度为四个5°箱格原点的均值；灰格表示无有效侧链角。','（a）泛素全部66个Hχ有效残基的联合熵预测与MD参考，按残基编号展开，横轴标出全部六个谷氨酰胺位点。')
s=s.replace("'all151_defined_angles_retained':True,'missing113_angle_slots_nan_not_zero':True,","'a_deleted_error_and_width_tracks':True,'all6_Gln_axis_labels':axis_labels,")
s=s.replace(",'MD_width_color_scale':[0,vmax]",'')
s=s.replace(",DATA/'all66_per_angle_metrics.csv'",'')
(FIG/'plot_Hchi_profile_six_Gln_axis_and_occupancy.py').write_text(s,encoding='utf8')
(ROOT/'manuscript_QA').mkdir(exist_ok=True)
print(FIG/'plot_Hchi_profile_six_Gln_axis_and_occupancy.py')
