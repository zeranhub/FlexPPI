# a: one structural-context panel, with one common component key.
txt(1,2,'a',12,fontweight='bold')
txt(17,2,'CaM–CaMKI',10,fontweight='bold')
txt(110,2,'Ubiquitin–UBA',10,fontweight='bold')
molecule('cam_overview',16,7,44,38.5,{'W303':(.14,.61),'M316':(.99,.28)})
molecule('ub_overview',109,7,44,38.5,{'Q41':(.56,.94),'Q62':(.96,.09)})
line(53,48.3,58,48.3,color=PURPLE,lw=2.6);txt(60,46.7,'Analysed chain')
line(109,48.3,114,48.3,color=GRAY,lw=2.6);txt(116,46.7,'Partner')

# b: one shared model/reference legend for b–d; sequence letters are omitted.
panel('b','Backbone RMSF',1,56)
for xx,label,col,mark,ls in [(116,'MD',INK,'o',':'),(137,'Full',BLUE,'o','-'),(158,'Intra',WARM,None,'--')]:
    line(xx,58.1,xx+5,58.1,color=col,lw=1.15,marker=mark,ms=2.7,markevery=[0],mfc='white' if label=='MD' else col,mec=col,ls=ls)
    txt(xx+6.5,56.5,label)
x=cam.CaMKI_native_position.to_numpy();lo,hi=293.5,318.5
a=ax(16,67,164,24)
a.axvspan(302.5,316.5,color='#F1F6F5',zorder=0)
for label,col,mark,ls,field in [('MD',INK,'o',':','target_rmsf'),('Full',BLUE,None,'-','pred_rmsf'),('Intra',WARM,None,'--','intra_rmsf')]:
    vals=cam[field].to_numpy()
    a.plot(x,vals,color=col,lw=1.15 if label=='Full' else 1,marker=mark,ms=2.1,ls=ls,mfc='white',mew=.6,zorder=4 if label=='Full' else 3)
    add_values('b','RMSF_A',x,vals,series=label)
a.set(xlim=(lo,hi),ylim=(0,8.2),ylabel='RMSF (Å)',xlabel='CaMKI residue')
a.set_yticks([0,4,8]);a.set_xticks([294,300,306,312,318]);a.tick_params(length=2,pad=2);a.xaxis.labelpad=2
for pos,label,loc in [(303,'W303',(303.1,1.8)),(316,'M316',(315.7,2.1))]:
    q=cam[cam.CaMKI_native_position.eq(pos)].iloc[0]
    a.annotate(label,xy=(pos,q.pred_rmsf),xytext=loc,fontsize=9.2,ha='center',color=TEAL,arrowprops=dict(arrowstyle='-',color=TEAL,lw=.6))
# The 25 static contact counts remain in cam_all25.csv and the manuscript.

# c/d: identical column boundaries to e/f; no repeated MD/Full legend.
panel('c','Side-chain Hχ',1,105)
panel('d','Six Gln positions',120,105)
a=ax(16,116,96,23)
for pos in [41,62]:a.axvspan(pos-.55,pos+.55,color='#EDF4F2',lw=0,zorder=0)
for label,col,mark,ls,field in [('MD',INK,'o',':','MD_Hchi_nats'),('Full',BLUE,None,'-','Full_Hchi_nats')]:
    vals=ub[field].where(ub.valid_hchi).to_numpy()
    a.plot(ub.residue_id,vals,color=col,lw=1.05 if label=='Full' else .85,marker=mark,ms=1.9,ls=ls,mfc='white',mew=.5,zorder=3 if label=='Full' else 2)
    add_values('c','Hchi_nats',ub.residue_id,vals,series=label)
for pos,loc in [(41,(44,3.1)),(62,(64,3.7))]:
    q=ub[ub.residue_id.eq(pos)].iloc[0]
    a.annotate('Q'+str(pos),xy=(pos,q.Full_Hchi_nats),xytext=loc,fontsize=9.2,color=TEAL,arrowprops=dict(arrowstyle='-',color=TEAL,lw=.55))
a.set(xlim=(.5,72.5),ylim=(0,4.3),ylabel='Hχ (nats)',xlabel='Ubiquitin residue')
a.set_xticks([1,20,40,60,72]);a.set_yticks([0,2,4]);a.tick_params(length=2,pad=2);a.xaxis.labelpad=2

a=ax(136,116,44,23)
for row in gln.itertuples():
    sx=row.sidechain_SASA_A2;md=row.MD_Hchi_nats;pred=row.Full_Hchi_nats
    a.plot([sx,sx],[md,pred],c='#AAB3B8',lw=.7,zorder=1)
    a.scatter(sx,md,s=19,facecolor='white',edgecolor=INK,linewidth=.8,zorder=3)
    a.scatter(sx,pred,s=15,facecolor=BLUE,edgecolor=BLUE,linewidth=.5,zorder=4)
    if row.residue in ['Q41','Q62']:
        loc={'Q41':(15,1.3),'Q62':(99,3.2)}[row.residue]
        a.annotate(row.residue,xy=(sx,(md+pred)/2),xytext=loc,fontsize=9.2,color=TEAL,ha='center',arrowprops=dict(arrowstyle='-',color=TEAL,lw=.55))
    add_values('d','Hchi_nats',[row.residue],[md],series='MD',SASA_A2=float(sx))
    add_values('d','Hchi_nats',[row.residue],[pred],series='Full',SASA_A2=float(sx))
a.set(xlim=(-9,114),ylim=(0,3.65),ylabel='Hχ (nats)',xlabel='Side-chain SASA (Å²)')
a.set_xticks([0,50,100]);a.set_yticks([0,1,2,3]);a.tick_params(length=2,pad=2);a.xaxis.labelpad=2

# e: all four angular subsets; stronger hierarchy with one shared x/y label.
defined=pd.read_csv(P/'defined_angles_long.csv',float_precision='round_trip')
chi_stats=pd.read_csv(P/'per_chi_correlations.csv',float_precision='round_trip')
chi_stats=chi_stats[chi_stats.predictor.eq('predicted_joint_Hchi')&chi_stats.metric.eq('fine_entropy')]
assert defined.groupby('chi').size().to_dict()=={'chi1':66,'chi2':52,'chi3':23,'chi4':10}
panel('e','Per-χ association',1,151)
for xx,label,marker,color in [(17,'Other','o',BLUE),(52,'Q41','D',TEAL),(84,'Q2','s',TEAL)]:
    line(xx,160.3,xx+2,160.3,color=color,marker=marker,ms=3.6,markevery=[0],lw=0)
    txt(xx+4,158.6,label)
for chi,left,top in [(1,16,171),(2,73,171),(3,16,200),(4,73,200)]:
    a=ax(left,top,39,18)
    z=defined[defined.chi.eq('chi'+str(chi))];st=chi_stats[chi_stats.chi.eq('chi'+str(chi))].iloc[0]
    txt(left,top-4.8,'χ'+str(chi)+'  n='+str(len(z))+'  ρ='+format(st.rho_AA_centered,'.2f'))
    for name,marker,col in [('Other','o',BLUE),('Q41','D',TEAL),('Q2','s',TEAL)]:
        sub=z[~z.author_residue_id.isin([2,41])] if name=='Other' else z[z.author_residue_id.eq(int(name[1:]))]
        a.scatter(sub.predicted_Hchi_nats_AA_centered,sub.fine_entropy_four_origins_mean_nats_AA_centered,s=10 if name=='Other' else 20,marker=marker,color=col,alpha=.5 if name=='Other' else 1,linewidths=.25,edgecolors='white' if name=='Other' else col,zorder=3 if name=='Other' else 4)
    for row in z.itertuples():
        plotted.extend([dict(panel='e',kind='predicted_joint_Hchi_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='Full',value=float(row.predicted_Hchi_nats_AA_centered)),dict(panel='e',kind='MD_single_chi_fine_entropy_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='MD',value=float(row.fine_entropy_four_origins_mean_nats_AA_centered))])
    a.axhline(0,color='#CDD4D8',lw=.45,ls=':');a.axvline(0,color='#CDD4D8',lw=.45,ls=':')
    a.set(xlim=(-.7,.5),ylim=(-1,.75));a.set_xticks([-.5,0,.5]);a.set_yticks([-.5,0,.5]);a.tick_params(length=2,pad=1.4)
    if chi in [2,4]:a.tick_params(labelleft=False)
    if chi in [1,2]:a.tick_params(labelbottom=False)
txt(3,194.5,'MD entropy deviation (nats)',9.2,rotation=90,ha='center',va='center')
txt(64,224.5,'Predicted joint Hχ deviation (nats)',9.2,ha='center')

# f: same cohort values; common right boundary with d, shorter direct title.
panel('f','Across 83 records',120,151)
for xx,label,mark,filled in [(130,'Hχ','o',True),(157,'RMSF','^',False)]:
    line(xx,160.3,xx+3,160.3,color=BLUE,lw=0,marker=mark,ms=4,markevery=[0],mfc=BLUE if filled else 'white',mec=BLUE)
    txt(xx+5,158.6,label)
a=ax(140,171,40,47)
for metric,base in [('fine_entropy',1.),('lowF_width',0.)]:
    for predictor,offset,marker,filled in [('full_hchi',.16,'o',True),('full_rmsf',-.16,'^',False)]:
        z=corr[corr.metric_key.eq(metric)&corr.predictor.eq(predictor)].sort_values('sample_id');assert len(z)==83
        jitter=np.linspace(-.085,.085,83)
        a.scatter(z.rho_type_centered,base+offset+jitter,s=6.5,marker=marker,facecolors=BLUE if filled else 'none',edgecolors=BLUE,linewidths=.3,alpha=.22,zorder=1)
        q=stats[stats.panel_kind.eq('correlation')&stats.metric_key.eq(metric)&stats.predictor.eq(predictor)].iloc[0]
        a.errorbar(q['mean'],base+offset,xerr=[[q['mean']-q.bootstrap_95_low],[q.bootstrap_95_high-q['mean']]],fmt=marker,ms=4.4,color=BLUE,mfc=BLUE if filled else 'white',mew=.8,elinewidth=1.15,capsize=2.1,zorder=5)
        add_values('f','within_AA_Spearman_rho',z.sample_id,z.rho_type_centered,series=predictor,metric=metric)
a.axvline(0,color='#CDD4D8',lw=.5,ls=':')
a.set(xlim=(-.05,.85),ylim=(-.42,1.42))
a.set_yticks([1,0],['Fine-bin\nentropy','Low-F*\nwidth']);a.set_xticks([0,.4,.8]);a.tick_params(length=2,pad=2)
txt(155,224.5,'Within-type ρ',9.2,ha='center')

df=pd.DataFrame(plotted);df.to_csv(P/'all_plotted_values.csv',index=False)
assert df.groupby('panel').size().to_dict()=={'b':75,'c':132,'d':12,'e':302,'f':332}
fig.canvas.draw();renderer=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in texts if t.get_window_extent(renderer).x0<-.5 or t.get_window_extent(renderer).y0<-.5 or t.get_window_extent(renderer).x1>fig.bbox.width+.5 or t.get_window_extent(renderer).y1>fig.bbox.height+.5]
assert not outside,outside
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'fig3_dual_case_entropy.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
qa=dict(size_mm=[W,H],min_font=min(t.get_fontsize() for t in texts),print_font_at_164mm=min(t.get_fontsize() for t in texts)*164/183,main_plot_counts=df.groupby('panel').size().to_dict(),all_original_mask_positions=True,all_six_Gln=True,all83_records_in_both_main_diagnostics=True,no_SASA_regression=True,no_tables_or_gray_notes=True,axes_with_absolute_amplitudes_start_at_zero=True,no_binding_free_energy_claim=True,not_apo_holo=True,structure_panels_grouped=True,local_zoom_panels_removed=True,common_columns_for_c_d_and_e_f=True,shared_MD_Full_Intra_legend=True,sequence_letters_not_displayed_but_all25_positions_retained=True,static_contact_counts_not_in_main_but_all25_retained_in_source=True,nonfocus_Gln_labels_removed_but_all6_points_retained=True,outside_text=outside,visual_review='pending',sha256={e:hashlib.sha256((OUT/f'fig3_dual_case_entropy.{e}').read_bytes()).hexdigest() for e in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(qa,ensure_ascii=False));plt.close(fig)
