from pathlib import Path
s=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025')
old=(s/'revision_v86_dense_evidence/figure3/compose_fig3_1MXE.py').read_text(encoding='utf-8')
code=old.replace('W,H=183.,225.','W,H=183.,190.').replace("R=ROOT/'pymol';F=ROOT/'pymol_context'","R=ROOT/'pymol';F=ROOT/'pymol_context';L=ROOT/'pymol_context_local'")
a=code.index("panel('a'");b=code.index('# e:')
new=r'''panel('a','CaM–CaMKI complex',1,2);panel('b','Recognition anchors',94,2)
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
a.set(ylim=(24.75,-.75),xlim=(0,12));a.set_xticks([0,3,6,9,12]);a.set_xlabel('Contacting CaM residues',labelpad=2,fontsize=7)
a.set_yticks(range(25),[r.amino_acid+str(r.CaMKI_native_position) for _,r in df.iterrows()],fontsize=6.5);a.tick_params(length=2,pad=1.4)
for t in a.get_yticklabels():
 if t.get_text() in ['W303','M316','N297','N308']:t.set_color(KEY);t.set_fontweight('bold')
for x,name,color in [(14,'N',NGRAY),(29,'C',CGRAY),(44,'Linker','#D9DDE0')]:
 fig.add_artist(Rectangle((x/W,1-66/H),2.2/W,1.4/H,transform=fig.transFigure,fc=color,ec='#777777',lw=.25));txt(x+3,64.5,name,6.8)
counts.index.name='CaMKI_native_position';counts.to_csv(P/'all25_lobe_contact_counts_4p5A.csv')
'''
code=code[:a]+new+code[b:]
code=code.replace('a=ax(79,84,27,58);b=ax(112,84,9,58)','a=ax(79,69,27,52);b=ax(112,69,9,52)').replace("txt(116.5,79,'Shielding'","txt(116.5,64,'Shielding'")
a=code.index('# f:');b=code.index('# c:')
new=r'''# f: the autoinhibited structure with its adjacent local magnification.
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
'''
code=code[:a]+new+code[b:]
for a,b in [("panel('c','Predicted CaMKI profiles',1,161)","panel('c','Complementary flexibility profiles',1,131)"),('1-173/H','1-141/H'),('168.5','136.5'),("txt(x,175,","txt(x,143,"),("('rmsf',181,8.2","('rmsf',149,8.2"),("('chi_entropy',203,4.3","('chi_entropy',171,4.3"),('pw,15)','pw,13)'),("fig3_1MXE_dense_evidence.","fig3_joint_flexibility_compact."),('all3600_contacts_displayed=True','all3600_distances_used_in_counts=True'),("context_render_present=context_path.exists()","context_render_present=True")]:code=code.replace(a,b)
code=code.replace("core_defined_by_structural_anchors=[303,316]","core_defined_by_structural_anchors=[303,316],contacts_aggregated_by_lobe=True,all25_contact_rows_retained=True,main_figure_height_reduced_percent=(225-H)/225*100,full_distance_matrix_in_source_data=True")
(s/'revision_v87_general_story_compact/figure3/compose_fig3_1MXE.py').write_text(code,encoding='utf-8')