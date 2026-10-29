from pathlib import Path
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'mel5_simplify_v52_20261022'
s=(OLD/'build_v52.py').read_text(encoding='utf-8')
s=s.replace('Remove the rejected local molecular close-ups and synchronize the v52 manuscript.','Integrate the residue-pair contact map and aligned model tracks into v53.')
s=s.replace("SOURCE = BASE / 'paper_v51_20261022_中文修订_MEL5跨界面残基机制与图4重绘_三线表.docx'", "SOURCE = BASE / 'paper_v52_20261022_中文修订_MEL5图4去除局部结构图_三线表.docx'")
s=s.replace("OUT = BASE / 'paper_v52_20261022_中文修订_MEL5图4去除局部结构图_三线表.docx'", "OUT = BASE / 'paper_v53_20261022_中文修订_MEL5残基接触与双动态联合图_三线表.docx'")
s=s.replace('fig4_mel5_compact.png','fig4_mel5_contact_map.png')
s=s.replace('2f99757f7063cc6506b25720dbe20a97197ac9d1e55eef0926d8837187f75e8d','bcf8a06760532e63b330341f315fc075acc8640f4b69e46c968e2538a74b61eb')
start=s.index('new = {}');end=s.index('changes=[]',start)
patch="""new = {}
method='局部原子间距由来源于3HG1与4JFF的MD准备起始坐标计算，分别测量Q72 OE1–β53主链N及Q72 NE2–βGly54主链O的距离。'
assert old[52].count(method)==1
new[52]=old[52].replace(method,'局部接触分析沿用来源于3HG1与4JFF的MD准备起始坐标，计算HLA第65–76位与TCR β48–56位的最短重原子间距。在任一态与该TCR片段距离≤6 Å且Hχ有效的六个HLA位置上，保留与九个TCR位置的全部54对距离；残基对距≤4 Å定义为近接触。图4a按这些位置对齐两态模型Hχ与RMSF。')
new[92]='为解释这两个位于界面两侧的模型响应，进一步将局部接触关系与两端的动态预测对齐（图4a）。亲和成熟后的βCDR2包含Gly52→Pro52、Ile53→Phe53等替换，原研究指出Ile53→Phe53带来了额外MHC接触[28]。在包含六个HLA邻接位点与TCR β48–56位的局部窗口中，Q72保留了与β51位的近接触，并新增与β52、β53、β54及β56位的近接触；Q72与β53、βGly54的最短重原子间距分别由5.46和8.27 Å缩短至2.94和3.12 Å，而与β55位的距离由3.94变为4.76 Å。这一接触重排将Q72侧链与βCDR2及相邻骨架联系到同一局部网络。模型同时预测到对侧Q72侧链状态集中和该受体片段的主链活动收紧，为界面两侧共同参与的局部稳定化提供动态解释。'
new[95]='图4：FlexPPI将MEL5亲和成熟后的界面接触与两侧动态响应对应到具体残基。（a）中部每格为一对HLA–TCR残基，橙色圆与紫色方块分别表示MEL5和α24β17的最短重原子距离：实心为≤4 Å近接触，空心为4–6 Å，小点为>6 Å。完整显示六个局部Hχ有效HLA位置与TCR β48–56位的全部54对距离；左侧Hχ预测与HLA行对齐，上方RMSF预测与TCR列对齐。突出Q72行及β53/β54列，Q72行数字为两态距离（Å），TCR位点下方标注氨基酸替换。局部位点按方法中的几何规则确定，复合物图仅作定位。（b）六个完整CDR中全部45个对齐位点的两态模型RMSF，识别环之间分隔，突出βCDR2及Gly54。颜色与形状在两部分保持一致，编号统一到MEL5；结构距离取自已有MD准备起始坐标，模型输入为氨基酸类型及N、Cα、C坐标。'
"""
s=s[:start]+patch+s[end:]
s=s.replace('Figure 4 MEL5 compact paired model readouts','Figure 4 MEL5 residue contact map with aligned dynamic predictions')
s=s.replace('Two panels: complex locator, published SPR affinity and compact paired Q72 Hchi / betaGly54 RMSF values; all45CDR model RMSF curves. Both local molecular close-ups removed.','Two panels: complete54pair HLA-TCR contact map with row-aligned Hchi and column-aligned RMSF predictions; all45CDR absolute model RMSF. Colored paired markers encode closest heavy-atom distance; no local molecular close-ups or text cards.')
s=s.replace("['5.85','2.94','10.58','3.12','1.792','1.375','1.118','0.747','Q72','βGly54','23个']", "['5.46','2.94','8.27','3.12','1.792','1.375','1.118','0.747','Q72','βGly54','23个']")
s=s.replace('prose_v52.json','prose_v53.json').replace('v52_content.txt','v53_content.txt')
(ROOT/'build_v53.py').write_text(s,encoding='utf-8')
for source,target in [('export_v52_word.ps1','export_v53_word.ps1'),('render_v52.py','render_v53.py')]:
    text=(OLD/source).read_text(encoding='utf-8-sig')
    text=text.replace('paper_v52_20261022_中文修订_MEL5图4去除局部结构图_三线表.docx','paper_v53_20261022_中文修订_MEL5残基接触与双动态联合图_三线表.docx')
    text=text.replace('v52_word.pdf','v53_word.pdf').replace('v52_page_map.json','v53_page_map.json')
    (ROOT/target).write_text(text,encoding='utf-8-sig' if target.endswith('.ps1') else 'utf-8')
print('v53 manuscript builder and render helpers prepared')
