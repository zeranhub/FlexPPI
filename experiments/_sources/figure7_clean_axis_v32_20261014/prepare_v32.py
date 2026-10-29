from pathlib import Path
ROOT=Path(__file__).resolve().parent;FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
OLD=ROOT.parent/'figure7_axis_sites_v31_20261014'
s=(OLD/'figures/plot_Hchi_profile_six_Gln_axis_and_occupancy.py').read_text(encoding='utf8')
s=s.replace("STEM='fig7_Hchi_profile_six_Gln_axis_and_occupancy'","STEM='fig7_clean_axis_six_Gln_and_occupancy'")
s=s.replace('pa.set_xticks([1,10,20,30,50,60,72])','pa.set_xticks([1,10,20,30,40,50,60,72])')
start=s.index('axis_labels=[]');end=s.index("text(3,116,'b'",start)
s=s[:start]+'''axis_labels=[];gln_texts=[]
for site in [2,31,40,41,49,62]:
    c=identity(site)
    label_x={40:38.,41:43.}.get(site,float(site))
    pa.plot([site,site],[0,-.055],transform=pa.get_xaxis_transform(),color=c,lw=.8,clip_on=False,zorder=5)
    pa.plot([site,label_x],[-.235,-.27],transform=pa.get_xaxis_transform(),color=c,lw=.55,clip_on=False,zorder=5)
    y=-.34
    t=pa.text(label_x,y,'Q'+str(site),transform=pa.get_xaxis_transform(),ha='center',va='center',fontsize=7.5,color=c,clip_on=False)
    gln_texts.append(t)
    axis_labels.append(dict(residue_id=site,label='Q'+str(site),anchor_x_data=float(site),label_x_data=label_x,y_axes=y))
text(95,121,'泛素残基编号',fontsize=7.5,ha='center')

'''+s[end:]
s=s.replace('assert not outside,outside', '''assert not outside,outside
label_boxes=[t.get_window_extent(renderer) for t in gln_texts]
assert all(not a.overlaps(b) for i,a in enumerate(label_boxes) for b in label_boxes[i+1:])
numeric_boxes=[t.get_window_extent(renderer) for t in pa.get_xticklabels()]
assert all(not a.overlaps(b) for a in label_boxes for b in numeric_boxes)
assert len(set(item['y_axes'] for item in axis_labels))==1''')
s=s.replace("'all6_Gln_axis_labels':axis_labels,","'all6_Gln_axis_labels':axis_labels,'numeric_axis_ticks':[1,10,20,30,40,50,60,72],'single_Gln_label_row':True,'axis_labels_no_overlap':True,")
(FIG/'plot_clean_axis_six_Gln_and_occupancy.py').write_text(s,encoding='utf8')
(ROOT/'manuscript_QA').mkdir(exist_ok=True)
s=(OLD/'build_v31.py').read_text(encoding='utf8')
src='paper_v31_20261014_中文修订_图7简化与六谷氨酰胺轴标注_三线表.docx'
out='paper_v32_20261014_中文修订_图7横轴标注整理_三线表.docx'
s=s.replace("SOURCE=BASE/'paper_v30_20261014_中文修订_全残基逐角占据与自由能验证_三线表.docx'","SOURCE=BASE/"+repr(src)).replace("OUT=BASE/"+repr(src),"OUT=BASE/"+repr(out))
s=s.replace("STEM='fig7_Hchi_profile_six_Gln_axis_and_occupancy'","STEM='fig7_clean_axis_six_Gln_and_occupancy'")
s=s.replace('bc85e75a3d6945e123e98910404de27c08d2d30c14bd114694934119ae696d47','a9011e448b28be019d4adb936c84e2ef45cd8196b0705f615b96817c0efdc2b7')
start=s.index("p=find('侧链实例使用泛素");end=s.index("caption=find('图7：')",start)
s=s[:start]+s[end:]
s=s.replace("caption=find('图7：');setp(caption,(ROOT/'figures'/(STEM+'.caption_zh.md')).read_text(encoding='utf8').strip())","caption=find('图7：')\nassert caption.text==(ROOT/'figures'/(STEM+'.caption_zh.md')).read_text(encoding='utf8').strip()")
s=s.replace('v31_content.txt','v32_content.txt')
s=s.replace("oldroot=ROOT.parent/'rich_angle_figure_v30_20261014'","oldroot=ROOT.parent/'figure7_axis_sites_v31_20261014'")
s=s.replace("[('export_v30_word.ps1','export_v31_word.ps1'),('render_v30.py','render_v31.py')]","[('export_v31_word.ps1','export_v32_word.ps1'),('render_v31.py','render_v32.py')]")
s=s.replace(".replace('v30_word','v31_word').replace('render_v30','render_v31').replace('v30_page_map','v31_page_map')",".replace('v31_word','v32_word').replace('render_v31','render_v32').replace('v31_page_map','v32_page_map')")
(ROOT/'build_v32.py').write_text(s,encoding='utf8')
print(ROOT)
