from pathlib import Path

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True)
source=ROOT.parent.parent/'figures/redraw_quantitative.py'
original=source.read_text(encoding='utf8')
new=original[:original.index('def main():')]
new=new.replace("OUT = Path(__file__).resolve().parent", "OUT = Path(__file__).resolve().parent\nENERGY = OUT.parent")
new += r'''

def save_new(fig,name,caption,statistics):
    save_cns_figure(fig,OUT/name)
    fig.savefig(OUT/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)
    (OUT/(name+'_caption_zh.md')).write_text(caption,encoding='utf8')
    (OUT/(name+'_statistics.json')).write_text(json.dumps(statistics,ensure_ascii=False,indent=2),encoding='utf8')
    MANIFEST['figures'][name]={'pdf':sha(OUT/(name+'.pdf')),'svg':sha(OUT/(name+'.svg')),
                              'png':sha(OUT/(name+'.png')),'all_input_rows_used':True}

def ubiquitin_distributions():
    hist=read(ENERGY/'case_analysis/ubiquitin_chi1_histograms_and_Fstar.csv')
    sens=read(ENERGY/'case_analysis/ubiquitin_bin_sensitivity.csv')
    angles={site:read(ENERGY/('case_analysis/'+site+'_angles_used.csv')) for site in ['Q41','Q62']}
    assert len(hist)==2*6*3*72
    assert all(len(d)==5000 for d in angles.values())
    colors={'Q41':CATEGORICAL[3],'Q62':CATEGORICAL[4]}
    fig,axes=plt.subplots(3,3,figsize=(183/25.4,190/25.4))
    fig.subplots_adjust(left=.10,right=.98,bottom=.065,top=.875,wspace=.35,hspace=.40)
    for col,chi in enumerate([1,2,3]):
        ax=axes[0,col];label(ax,'abc'[col],f'χ{chi}：角度时间序列')
        for site,d in angles.items():
            ax.scatter(d.time_ns,d[f'chi{chi}_degrees'],s=.45,alpha=.30,
                       color=colors[site],edgecolors='none',rasterized=True)
        ax.set(xlim=(0,500),ylim=(-180,180),xlabel='时间（ns）',ylabel='角度（°）')
        ax.set_yticks([-180,-60,60,180]);ax.set_xticks([0,250,500])
        ax=axes[1,col];label(ax,'def'[col],f'χ{chi}：经验概率密度')
        for site in ['Q41','Q62']:
            d=hist[(hist.residue==site)&(hist.window==0)&(hist.chi==chi)]
            assert d['count'].sum()==5000
            ax.plot(d.angle_center,d.density_per_degree,color=colors[site],lw=1.1,label=site)
        if chi==1:
            for b in [-60,60]:ax.axvline(b,color=GREY,lw=.7,ls='--',zorder=0)
            for b in [-120,0,120]:ax.axvline(b,color=BLACK,lw=.7,ls=':',alpha=.6,zorder=0)
        ax.set(xlim=(-180,180),xlabel='角度（°）',ylabel='密度（每度）')
        ax.set_xticks([-180,-60,60,180])
        ax=axes[2,col];label(ax,'ghi'[col],f'χ{chi}：经验投影自由能')
        for site in ['Q41','Q62']:
            for window in range(1,6):
                d=hist[(hist.residue==site)&(hist.window==window)&(hist.chi==chi)]
                ax.plot(d.angle_center,d.F_star,color=colors[site],lw=.45,alpha=.18)
            d=hist[(hist.residue==site)&(hist.window==0)&(hist.chi==chi)]
            ax.plot(d.angle_center,d.F_star,color=colors[site],lw=1.1)
        ax.set(xlim=(-180,180),xlabel='角度（°）',ylabel='F*（无量纲）')
        ax.set_xticks([-180,-60,60,180])
    fig.legend(handles=[Line2D([0],[0],color=colors[k],lw=1.5,label=k) for k in ['Q41','Q62']],
               loc='upper center',bbox_to_anchor=(.54,.998),ncol=2)
    text=[]
    for site in ['Q41','Q62']:
        r=sens[(sens.residue==site)&(sens.phase==0)].iloc[0]
        text.append(f"{site}  RMSF：MD / 模型 = {r.md_rmsf_A:.3f} / {r.pred_rmsf_A:.3f} Å；Hχ：{r.frozen_joint_h:.3f} / {r.pred_joint_h:.3f} nats")
    fig.text(.10,.948,text[0],fontsize=8,color=colors['Q41'])
    fig.text(.10,.923,text[1],fontsize=8,color=colors['Q62'])
    caption='图7：泛素Q41与Q62的连续侧链角度分布及经验投影自由能。两位点沿用图6的静态输入、标签与折外预测。（a–c）三个χ角在前5,000个保存帧中的时间序列，0–499.9 ns、间隔100 ps，每个位点全部5,000帧均绘出。（d–f）以5°等宽箱估计经验密度；χ1面板灰色虚线为原分区的−60°/+60°边界，黑色点线为常规χ1中心分区的−120°/0°/+120°边界，±180°为圆周接缝。（g–i）F*=−ln(p̂/p̂max)，最低点为零，粗线使用全部5,000帧，淡线为五个连续1,000帧窗口；零计数箱保留断线，无平滑或伪计数。各χ角投影由轨迹计算，图顶端RMSF及完整联合Hχ读数分别来自冻结标签与模型。'
    save_new(fig,'fig7_ubiquitin_angular_landscapes',caption,
             {'question':'How do paired RMSF/jointHchi outputs correspond to continuous marginal angular distributions?',
              'sites':2,'frames_per_site':5000,'angles':3,'density_bins_per_angle':72,
              'all_time_points':30000,'all_histogram_bins':len(hist),'source_sha256':dict(EVIDENCE),
              'Fstar_smoothing':False,'Fstar_empty_bins':'missing','Fstar_physical_units':'dimensionless',
              'single_angle_entropy_different_from_joint_target':True})

def representative_distributions():
    profile=read(ENERGY/'case_analysis/representative_full_profiles.csv')
    selected=read(ENERGY/'case_analysis/representative_selected_residues.csv')
    time=read(ENERGY/'case_analysis/representative_chi1_timeseries.csv')
    hist=read(ENERGY/'case_analysis/representative_chi1_Fstar.csv')
    audit=json.loads((ENERGY/'case_analysis/representative_cases_audit.json').read_text(encoding='utf8'))
    choices=audit['choices']
    fig,axes=plt.subplots(3,4,figsize=(183/25.4,170/25.4))
    fig.subplots_adjust(left=.09,right=.985,bottom=.065,top=.875,wspace=.50,hspace=.54)
    point_count=0
    for row,choice in enumerate(choices):
        sample=choice['sample_id'];d=profile[profile.sample_id==sample]
        chosen=selected[selected.sample_id==sample]
        letter='abcdefghijkl'[row*4:row*4+4]
        for col,(metric,unit) in enumerate([('rmsf','Å'),('hchi','nats')]):
            ax=axes[row,col]
            label(ax,letter[col],f"{sample}\n{['RMSF','Hχ'][col]}")
            ax.plot(d.input_index+1,d[metric+'_reference'],color=CATEGORICAL[5],lw=.55)
            ax.plot(d.input_index+1,d[metric+'_prediction'],color=CATEGORICAL[0],lw=.55)
            for i,site in enumerate(chosen.itertuples()):
                ax.axvline(site.input_index+1,color=CATEGORICAL[3+i],lw=.55,ls=':')
            for i in np.flatnonzero(d.chain.to_numpy()[1:]!=d.chain.to_numpy()[:-1]):
                ax.axvline(i+1.5,color=GREY,lw=.35,alpha=.5)
            ax.set(xlabel='输入残基顺序',ylabel=unit)
        ax=axes[row,2];label(ax,letter[2],'χ1：保存时间序列')
        ah=axes[row,3];label(ah,letter[3],'χ1：经验投影F*')
        handles=[]
        for i,site in enumerate(chosen.itertuples()):
            color=CATEGORICAL[3+i]
            v=time[(time.sample_id==sample)&(time.input_index==site.input_index)&(time.replica==1)]
            assert len(v)==501
            point_count+=len(v)
            ax.scatter(v.time_ns,v.chi1_degrees,s=.65,alpha=.5,color=color,edgecolors='none',rasterized=True)
            for rep in sorted(hist[hist.sample_id==sample].replica.unique()):
                h=hist[(hist.sample_id==sample)&(hist.input_index==site.input_index)&(hist.replica==rep)]
                assert len(h)==72 and h['count'].sum()==501
                ah.plot(h.chi1_degrees,h.F_star,color=color,lw=.9 if rep==1 else .5,
                        alpha=1 if rep==1 else .28)
            site_label=site.amino_acid+str(site.author_residue_id)
            if site.chain!='_':site_label=str(site.chain)+':'+site_label
            handles.append(Line2D([0],[0],color=color,lw=1,label=site_label))
        ax.set(xlim=(0,100),ylim=(-180,180),xlabel='时间（ns）',ylabel='角度（°）')
        ax.set_xticks([0,50,100]);ax.set_yticks([-180,0,180])
        ah.set(xlim=(-180,180),xlabel='角度（°）',ylabel='F*（无量纲）')
        ah.set_xticks([-180,0,180]);ah.legend(handles=handles,fontsize=6.5,loc='upper right')
    fig.legend(handles=[Line2D([0],[0],color=CATEGORICAL[5],label='MD参考标签'),
                        Line2D([0],[0],color=CATEGORICAL[0],label='折外预测')],
               loc='upper center',ncol=2,bbox_to_anchor=(.43,1.0))
    caption='图S13：三个来源的代表蛋白及局部角度分布。（a–d）ATLAS 2py5_A，（e–h）DynaRepo A0022，（i–l）ITrace 6p64。记录在各来源的83条轨迹子集内，按原折外两项目标MAE的来源内百分位秩均值最接近0.5选取；三条记录的均值秩均为0.5。每行前两面板使用完整输入残基顺序展示原MD标签及模型读数，灰色竖线为链边界。每条记录在相同氨基酸类型内，选取Hχ秩减RMSF秩最相离的一对非Pro位点，不使用预测误差选位点。后两面板给出该对位点第1条轨迹的全部501帧χ1及5°箱宽的F*，0–100 ns、间隔200 ps；粗线为第1条轨迹，淡线为其他已保存重复，零计数箱断开，无平滑。轨迹窗口用于连续角度分布，前两面板保留原标签的观测汇总口径。位点标签为模型输入PDB中的链与残基编号。'
    save_new(fig,'figS13_representative_distribution_cases',caption,
             {'proteins':3,'selected_sites':6,'first_replica_time_points':point_count,
              'all_profile_rows_used':len(profile),'all_histogram_rows_used':len(hist),
              'replica_curve_count':int(hist[['sample_id','replica','input_index']].drop_duplicates().shape[0]),
              'selection':audit['choices'],'source_sha256':dict(EVIDENCE),
              'Fstar_empty_bins':'missing','smoothing':False,'all_raw_source_caches_retained':True})

def main():
    OUT.mkdir(exist_ok=True)
    ubiquitin_distributions()
    representative_distributions()
    MANIFEST['sources']=dict(EVIDENCE)
    MANIFEST['visual_review']='pending'
    (OUT/'case_figure_manifest.json').write_text(json.dumps(MANIFEST,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({'figures':list(MANIFEST['figures'])},ensure_ascii=False))

if __name__=='__main__':main()
'''
(OUT/'plot_case_landscapes.py').write_text(new,encoding='utf8')
print(OUT/'plot_case_landscapes.py')
