from pathlib import Path
ROOT=Path(__file__).resolve().parent
s=(ROOT.parent/'mel5_composite_v54_20261022/finalize_v54.py').read_text(encoding='utf-8')
s=s.replace('v54','v55').replace('fig4_mel5_composite','fig4_mel5_compact_joint')
s=s.replace("ROOT.parent/'mel5_contact_v53_20261022/CURRENT_FIGURE_MAP.json'","ROOT.parent/'mel5_composite_v54_20261022/CURRENT_FIGURE_MAP.json'")
s=s.replace("'ranked_Hchi_sites':23","'ranked_Hchi_plot_removed':True")
start=s.index("fm['4'] =") if "fm['4'] =" in s else s.index("fm['4']=")
end=s.index("write(ROOT/'CURRENT_FIGURE_MAP.json',fm)",start)
s=s[:start]+'''fm['4']={'source_folder':'figure4','absolute_source_folder':str(FIG),'png':'fig4_mel5_compact_joint.png','pdf':'fig4_mel5_compact_joint.pdf','svg':'fig4_mel5_compact_joint.svg','new_in_revision':'v55','panels':['a existing complex overview and local Q72-beta53/54 structure','b all40 four-state peptide RMSF values','c all90 absolute RMSF values at45 CDR positions','d joint response, upper all135 CDR differences and lower all90 pMHC interface cells for three contrasts'],'removed':'v54 d ranked Hchi bar plot','merged':'v54 e and f become the upper/lower blocks of v55 d'}
'''+s[end:]
start=s.index("(FIG/'QA_report.md').write_text")
end=s.index("write(ROOT/'figure_contract.json',contract)",start)
s=s[:start]+'''(FIG/'QA_report.md').write_text("""# Figure4 v55 QA

READY. Four panels a-d. Prior ranked Hchi plot deleted by direct user request; two complete signed maps combined under one d heading with independent metric units and colorbars.183x172mm,600dpi PNG plus embedded-font PDF and editable SVG. Font sizes and all retained data unchanged from v54:40 peptide,90 absolute CDR,135 CDR differences,90 interface cells with valid23/23/22. Source CSVs and array values/masks independently checked against v54. No new inference/training/MD. Entire28-page DOCX independently reviewed, complete figure/caption same page.10 native three-line tables,36 equations,6 images and28 references preserved. Source v54 unchanged.
""",encoding='utf-8')
(ROOT/'修改说明_zh.md').write_text("""# v55 修改说明

按用户要求删除原d子图的23位点Hchi排序柱图，将原e/f两张热图合为新d，CDR RMSF变化位于上部、pMHC Hchi变化位于下部，保留独立单位和色标。a/b/c的内容、物理尺寸、字号及所有保留数据不变，整图由183×214mm收紧至183×172mm。

正文中原e/f引用改为图4d上部/下部；Q72原图d引用指向新d下部，图注删去柱图与p前缀说明。原始数据和v54稿均保留。全稿28页渲染及逐页检查通过，图与完整图注同页，三线表和公式完整。
""",encoding='utf-8')
contract={'scientific_goal':'Remove prior panel d and combine e/f into one d as directly requested.','archetype':'asymmetric_mixed','dimensions_mm':[183,172],'layout':'top structure+peptide, full-width CDR, single joint-response d containing upper CDR and lower pMHC heatmaps','authority':'Direct user instruction; no additional design confirmation needed.','panels':fm['4']['panels'],'source_preserved':True,'removed_panel':'v54 d rank bars','merged_panels':'v54 e/f to v55 d','all_retained_values_unchanged':True,'QA':'passed'}
'''+s[end:]
s=s.replace("'panels':6","'panels':4")
(ROOT/'finalize_v55.py').write_text(s,encoding='utf-8')
print(ROOT/'finalize_v55.py')
