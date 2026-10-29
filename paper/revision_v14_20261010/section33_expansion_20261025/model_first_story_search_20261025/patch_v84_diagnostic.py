from pathlib import Path
p=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v84_effect_audit\compose_effect_diagnostic.py')
s=p.read_text(encoding='utf-8')
s=s.replace("txt(157,38,'Side-chain sites',7.8,ha='center',fontweight='bold');a=ax(140,43,39,15)","txt(157,37,'Side-chain sites',7.8,ha='center',fontweight='bold');a=ax(140,41,39,21)")
old="""if v<1.1:
            offset=(5,0) if s=='b' else (-5,0)
            ha='left' if s=='b' else 'right';va='center'
        else:
            above=v==max(q.full0_a,q.full0_b)
            offset=(0,3 if above else -4);ha='center';va='bottom' if above else 'top'
        a.annotate(f'{v:.2f}',(j+dx,v),xytext=offset,textcoords='offset points',ha=ha,va=va,fontsize=6.9,color=BLACK)"""
new="""above=v==max(q.full0_a,q.full0_b)
        a.annotate(f'{v:.2f}',(j+dx,v),xytext=(0,3 if above else -4),textcoords='offset points',ha='center',va='bottom' if above else 'top',fontsize=6.9,color=BLACK)"""
assert old in s
p.write_text(s.replace(old,new),encoding='utf-8')