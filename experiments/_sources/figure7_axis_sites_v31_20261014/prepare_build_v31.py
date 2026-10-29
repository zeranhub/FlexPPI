from pathlib import Path
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'rich_angle_figure_v30_20261014'
s=(OLD/'build_v30.py').read_text(encoding='utf8')
source='paper_v30_20261014_中文修订_全残基逐角占据与自由能验证_三线表.docx'
out='paper_v31_20261014_中文修订_图7简化与六谷氨酰胺轴标注_三线表.docx'
s=s.replace("SOURCE=BASE/'paper_v29_20261014_中文修订_全残基与同类型预测及自由能曲线_三线表.docx'","SOURCE=BASE/"+repr(source))
s=s.replace("OUT=BASE/'"+source+"'","OUT=BASE/"+repr(out))
s=s.replace("STEM='fig7_full_residue_and_chi_occupancy'","STEM='fig7_Hchi_profile_six_Gln_axis_and_occupancy'")
s=s.replace("'dc39f5dacad17bbb284621083c39834186e87d53a3a73a49d2cf5992006fa904'","'bc85e75a3d6945e123e98910404de27c08d2d30c14bd114694934119ae696d47'")
start=s.index("p=find('侧链实例使用泛素");end=s.index("caption=find('图7：')")
s=s[:start]+'''p=find('侧链实例使用泛素–EDD UBA复合物')
setp(p,p.text.replace('图7按残基编号对齐泛素链全部66个Hχ有效位置的冻结预测、主评价参考及逐位误差，并展示全部151个可定义χ角的低自由能角区总宽度。','图7按残基编号展开泛素链全部66个Hχ有效位置的冻结预测与主评价参考，并在横轴标出全部六个谷氨酰胺位点。'))
'''+s[end:]
s=s.replace('width=Inches(6.0)','width=Inches(6.2)')
s=s.replace("QA/'v30_content.txt'","QA/'v31_content.txt'")
start=s.index('oldroot=ROOT.parent/');end=s.index('print(json.dumps',start)
s=s[:start]+'''oldroot=ROOT.parent/'rich_angle_figure_v30_20261014'
for source,target in [('export_v30_word.ps1','export_v31_word.ps1'),('render_v30.py','render_v31.py')]:
    text=(oldroot/source).read_text(encoding='utf8').replace(SOURCE.name,OUT.name).replace('v30_word','v31_word').replace('render_v30','render_v31').replace('v30_page_map','v31_page_map')
    (ROOT/target).write_text(text,encoding='utf8')
'''+s[end:]
(ROOT/'build_v31.py').write_text(s,encoding='utf8')
print(ROOT/'build_v31.py')
