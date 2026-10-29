from pathlib import Path
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'mel5_interface_v51_20261022'
s=(OLD/'build_v51.py').read_text(encoding='utf-8')
s=s.replace('Integrate the authorized MEL5 two-panel interface revision into v50.','Remove the rejected local molecular close-ups and synchronize the v52 manuscript.')
s=s.replace("SOURCE = BASE / 'paper_v50_20261022_中文修订_MEL5残基接触网络与图4精简_三线表.docx'", "SOURCE = BASE / 'paper_v51_20261022_中文修订_MEL5跨界面残基机制与图4重绘_三线表.docx'")
s=s.replace("OUT = BASE / 'paper_v51_20261022_中文修订_MEL5跨界面残基机制与图4重绘_三线表.docx'", "OUT = BASE / 'paper_v52_20261022_中文修订_MEL5图4去除局部结构图_三线表.docx'")
s=s.replace('fig4_mel5_interface.png','fig4_mel5_compact.png')
s=s.replace('9c05035e9af5ba5a7c2291f29ef939be7bc5e2d676d3e30c8d6c3e59cc36aaf8','2f99757f7063cc6506b25720dbe20a97197ac9d1e55eef0926d8837187f75e8d')
start=s.index('new = {');end=s.index('changes=[]',start)
patch="""new = {}
opening='将两个复合物以HLA对齐并从同一视角观察，可直接比较受体改造前后的局部配对关系（图4a）。'
assert old[92].count(opening)==1
new[92]=old[92].replace(opening,'为解释这两个位于界面两侧的模型响应，进一步检查了βCDR2附近的局部接触。')
method='跨界面局部结构沿用来源于3HG1与4JFF的MD准备起始坐标，以HLA第1–180位共同Cα原子对齐两个复合物，并用相同视角和尺度显示。图中标注Q72 OE1–β53主链N和Q72 NE2–βGly54主链O的原子间距。'
assert old[52].count(method)==1
new[52]=old[52].replace(method,'局部原子间距由来源于3HG1与4JFF的MD准备起始坐标计算，分别测量Q72 OE1–β53主链N及Q72 NE2–βGly54主链O的距离。')
new[95]='图4：FlexPPI定位MEL5亲和成熟后的界面双动态响应。（a）复合物概览、文献[28]的SPR亲和力，以及界面两侧HLA Q72的Hχ与TCR βGly54的RMSF两态模型读数。两种TCR结合相同的HLA-A*02:01与改良MART-1十肽ELAGIGILTV。（b）六个完整CDR中全部45个对齐位点的两态模型RMSF，识别环之间分隔，突出βCDR2及Gly54。橙色为MEL5，紫色为α24β17；编号统一到MEL5，配对与共同编号规则见方法。'
"""
s=s[:start]+patch+s[end:]
s=s.replace('Figure 4 MEL5 paired cross-interface residue contacts','Figure 4 MEL5 compact paired model readouts')
s=s.replace('Two panels: shared-view and shared-scale MEL5/alpha24beta17 interface structures, Q72-beta53/beta54 distances and direct FlexPPI readouts; all45CDR model RMSF curves. No panel c or d.','Two panels: complex locator, published SPR affinity and compact paired Q72 Hchi / betaGly54 RMSF values; all45CDR model RMSF curves. Both local molecular close-ups removed.')
s=s.replace("case='\\n'.join(new[i] for i in range(89,94))", "case='\\n'.join(txt(pars[i]) for i in range(89,94))")
s=s.replace('prose_v51.json','prose_v52.json').replace('v51_content.txt','v52_content.txt')
(ROOT/'build_v52.py').write_text(s,encoding='utf-8')
for source,target in [('export_v51_word.ps1','export_v52_word.ps1'),('render_v51.py','render_v52.py')]:
    text=(OLD/source).read_text(encoding='utf-8-sig')
    text=text.replace('paper_v51_20261022_中文修订_MEL5跨界面残基机制与图4重绘_三线表.docx','paper_v52_20261022_中文修订_MEL5图4去除局部结构图_三线表.docx')
    text=text.replace('v51_word.pdf','v52_word.pdf').replace('v51_page_map.json','v52_page_map.json')
    (ROOT/target).write_text(text,encoding='utf-8-sig' if target.endswith('.ps1') else 'utf-8')
print('v52 builder and render helpers prepared')
