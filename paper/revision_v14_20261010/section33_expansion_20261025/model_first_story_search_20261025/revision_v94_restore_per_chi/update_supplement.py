from pathlib import Path
import shutil,subprocess,sys
R=Path(__file__).resolve().parent;OLD=R.parent/'revision_v93_dual_case_entropy'
O=R/'supplement4';O.mkdir(exist_ok=True)
shutil.copytree(OLD/'entropy_audit/supplement_source',O/'supplement_source',dirs_exist_ok=True)
for name in ['retained_prediction_groups_all498.csv','retained_panel_statistics.csv']:
    shutil.copy2(OLD/'entropy_audit'/name,O/name)
source=(OLD/'entropy_audit/plot_supplement_S4.py').read_text(encoding='utf-8')
head=source.split('# This module body is appended')[0]
body=(OLD/'entropy_audit/supplement_body.py').read_text(encoding='utf-8')
start=body.index('for path in source_paths:');end=body.index('\necdf = read',start)
body=body[:start]+"# Frozen input copies are bundled in supplement_source.\n"+body[end:]
start=body.index('# c: complete 151');end=body.index('# d: all 432',start)
replacement='''# c: prediction-selected actual MD metrics moved from the old main figure.
panel('c',2,100,'Actual MD occupancy in prediction-selected halves')
selection=read(OUT/'retained_prediction_groups_all498.csv')
selection_stats=read(OUT/'retained_panel_statistics.csv')
for metric,left,title,lim,ticks in [('fine_entropy',30,'Fine-bin entropy (nats/χ)',(0,3.4),[0,1.5,3]),('lowF_width',116,'Low-F* width (°/χ)',(0,60),[0,30,60])]:
    ax=axis(left,64,46,27);register('c_'+metric,ax)
    part=selection[selection.metric_key.eq(metric)&selection.predictor.eq('full_hchi')].sort_values('sample_id').reset_index(drop=True)
    assert len(part)==83
    for j,r in part.iterrows():
        jitter=((j*37)%83/82-.5)*.18
        ax.plot([0+jitter,1+jitter],[r.low_MD_mean,r.high_MD_mean],color='#AEB7BC',lw=.45,alpha=.5,zorder=1)
    for xi,field in enumerate(['low_MD_mean','high_MD_mean']):
        jitter=((np.arange(83)*37)%83/82-.5)*.18
        ax.scatter(xi+jitter,part[field],s=6,color=BLACK,alpha=.3,linewidths=0,zorder=2)
        st=selection_stats[selection_stats.panel_kind.eq('prediction_selection')&selection_stats.metric_key.eq(metric)&selection_stats.predictor.eq('full_hchi')&selection_stats.value_key.eq(field)].iloc[0]
        ax.errorbar(xi,st['mean'],yerr=[[st['mean']-st.bootstrap_95_low],[st.bootstrap_95_high-st['mean']]],fmt='D',ms=4.5,mfc='white',mec=BLACK,mew=.8,ecolor=BLACK,elinewidth=1,capsize=2,zorder=5)
        for j,r in part.iterrows():all_points.append(dict(panel='c',metric=metric,sample_id=r.sample_id,kind=field,x=float(xi+jitter[j]),y=r[field]))
    ax.set(xlim=(-.3,1.3),ylim=lim);ax.set_xticks([0,1],['Low','High']);ax.set_yticks(ticks)
    ax.set_title(title,fontsize=8.5,loc='left',pad=4)
text(99,55,'Groups selected by predicted joint Hχ',fontsize=8.5,ha='center')

'''
body=body[:start]+replacement+body[end:]
body=body.replace("{'a':6126,'b':498,'c':151,'d':432}","{'a':6126,'b':498,'c':332,'d':432}")
(O/'plot_supplement_S4_v94.py').write_text(head+body,encoding='utf-8')
subprocess.run([sys.executable,str(O/'plot_supplement_S4_v94.py')],check=True)
