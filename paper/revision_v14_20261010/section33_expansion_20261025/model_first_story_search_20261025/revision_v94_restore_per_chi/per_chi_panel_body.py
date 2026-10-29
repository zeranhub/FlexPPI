# f: restore the requested two-by-two direct model-to-MD per-chi association.
defined=pd.read_csv(P/'defined_angles_long.csv',float_precision='round_trip')
chi_stats=pd.read_csv(P/'per_chi_correlations.csv',float_precision='round_trip')
chi_stats=chi_stats[chi_stats.predictor.eq('predicted_joint_Hchi')&chi_stats.metric.eq('fine_entropy')]
assert defined.groupby('chi').size().to_dict()=={'chi1':66,'chi2':52,'chi3':23,'chi4':10}
panel('f','Per-χ occupancy association',1,181)
for xx,label,marker,color in [(12,'Other','o',BLUE),(40,'Q41','D',TEAL),(64,'Q2','s',TEAL)]:
    line(xx,189.3,xx+2,189.3,color=color,marker=marker,ms=3.6,lw=0);txt(xx+4,187.6,label,9)
for chi,left,top in [(1,19,198),(2,63,198),(3,19,224),(4,63,224)]:
    a=ax(left,top,27,14)
    z=defined[defined.chi.eq('chi'+str(chi))];st=chi_stats[chi_stats.chi.eq('chi'+str(chi))].iloc[0]
    txt(left,top-4.6,'χ'+str(chi)+' n='+str(len(z))+' ρ='+format(st.rho_AA_centered,'.2f'),9)
    for name,marker,col in [('Other','o',BLUE),('Q41','D',TEAL),('Q2','s',TEAL)]:
        sub=z[~z.author_residue_id.isin([2,41])] if name=='Other' else z[z.author_residue_id.eq(int(name[1:]))]
        a.scatter(sub.predicted_Hchi_nats_AA_centered,sub.fine_entropy_four_origins_mean_nats_AA_centered,s=10 if name=='Other' else 19,marker=marker,color=col,alpha=.55 if name=='Other' else 1,linewidths=.25,edgecolors='white' if name=='Other' else col,zorder=3 if name=='Other' else 4)
    for row in z.itertuples():
        plotted.extend([dict(panel='f',kind='predicted_joint_Hchi_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='Full',value=float(row.predicted_Hchi_nats_AA_centered)),dict(panel='f',kind='MD_single_chi_fine_entropy_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='MD',value=float(row.fine_entropy_four_origins_mean_nats_AA_centered))])
    a.axhline(0,color='#AAB3B8',lw=.5,ls=':');a.axvline(0,color='#AAB3B8',lw=.5,ls=':')
    a.set(xlim=(-.7,.5),ylim=(-1,.75));a.set_xticks([-.5,0,.5]);a.set_yticks([-.5,0,.5]);a.tick_params(length=2,pad=1.4)
    if chi in [2,4]:a.tick_params(labelleft=False)
    if chi in [1,2]:a.tick_params(labelbottom=False)
txt(4.8,218,'MD entropy deviation (nats)',9,rotation=90,ha='center',va='center')
txt(55,242.5,'Predicted joint Hχ deviation (nats)',9,ha='center')

# g: retain both cohort descriptors, rather than only plotting isolated MD curves.
panel('g','Association across 83 records',99,181)
for xx,label,mark,filled in [(110,'Hχ','o',True),(142,'RMSF','^',False)]:
    line(xx,190,xx+3,190,color=BLUE,lw=0,marker=mark,ms=4,mfc=BLUE if filled else 'white',mec=BLUE);txt(xx+5,188.3,label,9)
a=ax(124,199,56,33)
for metric,base in [('fine_entropy',1.),('lowF_width',0.)]:
    for predictor,offset,marker,filled in [('full_hchi',.16,'o',True),('full_rmsf',-.16,'^',False)]:
        z=corr[corr.metric_key.eq(metric)&corr.predictor.eq(predictor)].sort_values('sample_id');assert len(z)==83
        jitter=np.linspace(-.085,.085,83)
        a.scatter(z.rho_type_centered,base+offset+jitter,s=7,marker=marker,facecolors=BLUE if filled else 'none',edgecolors=BLUE,linewidths=.35,alpha=.22,zorder=1)
        q=stats[stats.panel_kind.eq('correlation')&stats.metric_key.eq(metric)&stats.predictor.eq(predictor)].iloc[0]
        a.errorbar(q['mean'],base+offset,xerr=[[q['mean']-q.bootstrap_95_low],[q.bootstrap_95_high-q['mean']]],fmt=marker,ms=4.6,color=BLUE,mfc=BLUE if filled else 'white',mew=.8,elinewidth=1.2,capsize=2.2,zorder=5)
        add_values('g','within_AA_Spearman_rho',z.sample_id,z.rho_type_centered,series=predictor,metric=metric)
a.axvline(0,color='#BFC5C9',lw=.6,ls=':');a.set(xlim=(-.05,.85),ylim=(-.42,1.42),xlabel='Within-type Spearman ρ');a.set_yticks([1,0],['Fine-bin\nentropy','Low-F*\nwidth']);a.set_xticks([0,.4,.8]);a.tick_params(length=2,pad=2);a.xaxis.labelpad=2
