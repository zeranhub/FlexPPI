from pathlib import Path
p=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v86_dense_evidence\figure3\compose_fig3_1MXE.py')
s=p.read_text(encoding='utf-8')
s=s.replace("panel('e','All sidechain states',65,72)","panel('e','Sidechain environments',65,72)")
s=s.replace('fontsize=6.1)','fontsize=6.5)')
s=s.replace('a.set_xticks([4,35,70,80,110,147])','a.set_xticks([4,40,75,110,147])')
s=s.replace("context_path=F/'context_render_manifest.json'","context_path=F/'render_manifest.json'")
s=s.replace("('cam_bound',113,'CaM-bound peptide (1MXE)')","('CaM_bound',113,'CaM-bound peptide (1MXE)')")
s=s.replace("inf=ctx['assets'][key];a=ax(125,y,55,26)","inf=ctx['assets'][key]['trim'];a=ax(125,y,55,26)")
s=s.replace("xytext=((.08,.87) if site=='W303' else (.88,.13))","xytext=((.03,.90) if site=='W303' else (.91,.52))")
s=s.replace("else:txt(152.5,105,'Structural context rendering',7,ha='center')", """  if key=='autoinhibited':
   lm=inf['landmarks']
   if 'kinase_contacts' in lm:
    a.annotate('L45 / S99',xy=lm['kinase_contacts'],xytext=(.70,.05),textcoords='axes fraction',color=BLACK,fontsize=6.4,ha='center',arrowprops=dict(arrowstyle='-',lw=.45,color=BLACK))
else:raise FileNotFoundError(context_path)""")
s=s.replace("context_render_present=context_path.exists()","context_render_present=context_path.exists()")
p.write_text(s,encoding='utf-8')