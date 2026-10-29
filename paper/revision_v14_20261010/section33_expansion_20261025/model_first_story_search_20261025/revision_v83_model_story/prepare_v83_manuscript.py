from pathlib import Path
import json,shutil
r=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');old=r/'revision_v82_shared_view';v=r/'revision_v83_model_story';m=v/'manuscript';m.mkdir(exist_ok=True);(v/'render').mkdir(exist_ok=True)
shutil.copy2(old/'manuscript/v79_paragraph_inventory.json',m/'v79_paragraph_inventory.json')
for name in ['render_docx_word_wrapper.py','export_word_readonly.ps1']:shutil.copy2(old/'render'/name,v/'render'/name)
plan=json.loads((old/'manuscript/v82_text_plan.json').read_text(encoding='utf-8'));plan['output_docx']='E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v83_20261025_中文修订_图3模型柔性再分配与关键区域.docx'
plan['replacements']['68']='3.3 CaM两种识别构型的预测柔性再分配'
plan['replacements']['69']='钙调蛋白（CaM）通过两个叶瓣及连接区适应不同靶蛋白的识别环境。NHE1与IQCG复合物分别源自2YGG与4LZX，展示同一CaM序列在不同伙伴及叶瓣排列下的识别构型[29,30]。本例将DynaRepo两条记录中成熟CaM的4–148号逐位对应，用同一个冻结的第0折完整模型预测145个共同主链位置及122个有效侧链位置。两条记录均在该折训练中，本节将其用于结构解释；独立性能由3.1–3.2节留出评价建立。图3集中展示模型预测的两态绝对读数及其结构位置，MD参考和变化准确度核查移至图S3与补充材料D。'
plan['replacements']['70']='模型预测的主链变化与侧链变化呈现不同的空间模式。145个共同位置的平均RMSF由NHE1态的1.333降至IQCG态的0.954 Å，其中135个位点降低、10个升高；122个侧链的平均Hχ则由1.732变为1.750 nats，分别有60个升高、62个降低（图3c）。因此，较广泛的主链波动降低并不伴随侧链状态占据的一致集中。下面依次以完整接触分组、全部同类型Met及具体位点，解释这种位置依赖的再分配（图3d–f）。'
plan['replacements']['73']='输入结构中的接触重排提供了区域层面的解释。以跨链Cα距离≤8 Å定义伙伴接触，将122个侧链完整分为仅NHE1接触、仅IQCG接触、两态均接触及两态均不接触四组，分别含8、21、3和90个位点（图3d）。仅IQCG接触组的预测平均Hχ由1.857降至1.634 nats，两态均接触组由1.811降至1.589 nats；相反，仅NHE1接触组由1.891升至1.944 nats，两态均不接触组由1.686升至1.765 nats。相同CaM序列中，新接触区域与其他区域因而可呈相反的侧链状态变化，而非整条蛋白同步增减。'
plan['replacements']['74']='全部九个Met进一步排除了残基类型差异，并呈现多种主链与侧链组合（图3b、e、f）。M72的预测RMSF由0.877小幅降至0.810 Å，Hχ却由1.862升至2.112 nats，表现为主链变化较小而侧链占据多样性增加；M109的RMSF由1.396降至0.715 Å，Hχ由1.971降至1.759 nats，两项读数同时降低。M144则呈另一种组合：RMSF由1.345降至1.060 Å，Hχ由2.300升至2.621 nats。其相邻M145的Hχ反而由2.370降至2.159 nats，说明相邻且同类型的侧链也可发生不同方向的状态调整。'
plan['replacements']['75']='图3c按关键位点前后各四个残基框出68–76、105–113及140–148号邻域，并标明具体Met位置。68–76号的8个有效侧链Hχ全部升高，区域平均Hχ由1.966增至2.173 nats，而9个主链位置的平均RMSF仅由0.943变为0.962 Å，将主链与侧链的不同响应从单个位点延伸到局部片段。这里Hχ表示固定角区下多角联合状态的占据多样性，不等同于转动更快；RMSF沿用整链拟合口径，包含叶间运动贡献。伙伴、构型与原始钙条件同时不同，因此不将预测差异单独归因为更换伙伴或去钙。'
plan['insert_before_158'][2]+=' 图S3保留全部145个主链和122个侧链的MD参考及绝对值、变化量核查。完整模型的侧链变化排序相关为0.630；在MD绝对变化至少0.1 nats的70个位点中，升高与降低分别判对32/40和26/30，共58/70。全部九个Met的侧链变化排序相关为0.883，五个MD绝对变化至少0.1 nats的位置均判对方向；M36与M124的方向或幅度偏差仍保留在源表，不用于主图的准确度叙述。'
plan['notes']=['主图全改为冻结Full预测两态比较；MD及准确度移图S3。','结构组分灰/紫/绿与数值两态蓝/朱橙分离；d–f在c上方。','c只保留两行模型全序列，并以预先选定位点±4残基框邻域；f直接显示M72/M109双读数及原值。']
(m/'v83_text_plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
caption='图3：CaM两种识别构型的预测柔性再分配。所有定量面板均为同一冻结Full模型的原值，蓝色圆点/曲线表示NHE1，朱橙方点/曲线表示IQCG。（a）共同参考视角下的输入结构；CaM为灰色、伙伴为紫色、关键Met为绿色。显示时以N叶螺旋I/IV的Cα（7–19、65–75号）对整个复合物刚体叠合，保留真实构型差异。（b）M72与M109的局部环境，各按位点前后五残基主链叠合，同一位点两态同比例。（c）145个主链及122个有效侧链位置的两态预测；框为M72、M109、M144前后各四残基的邻域，标记为对应Met，无有效Hχ处留空。（d）按两态输入跨链Cα距离≤8 Å定义的四个完整接触组，括号为侧链数，点为组均值。（e）全部九个Met的Hχ，加粗及黑边突出三个关键位置。（f）M72、M109的RMSF与Hχ配对原值，连线对应同一位点，数字保留两位小数。两条记录均用于第0折训练，本图用于结构解释；MD参考与预测准确度见图S3。'
(v/'figure3/caption_zh.md').write_text(caption+'\n',encoding='utf-8')
s=(old/'manuscript/build_v82_manuscript.py').read_text(encoding='utf-8').replace('v82_text_plan.json','v83_text_plan.json').replace('v82_source_audit.json','v83_source_audit.json').replace('v82_build_audit.json','v83_build_audit.json').replace('v82_final_text.txt','v83_final_text.txt')
s=s.replace("        paragraph_setting(paragraphs[73], 'pageBreakBefore', val=1)","        for case_index in (73, 74):\n            paragraph_setting(paragraphs[case_index], 'keepLines', val=1)")
s=s.replace("for name in ('fig3', 'fig5', 'fig_s2', 'caption3', 'caption5'):","for name in ('fig3', 'fig5', 'fig_s2', 'fig_s3', 'caption3', 'caption5', 'caption_s3'):")
s=s.replace("    parser.add_argument('--output')","    parser.add_argument('--fig-s3'); parser.add_argument('--caption-s3')\n    parser.add_argument('--output')")
needle="        assert protected_34 == [xml_signature(p) for p in paragraphs[76:83]]"
addition="""        # Preserve the reference comparisons in a separate supplemental figure.
        s3_path = Path(args.fig_s3)
        s3_paragraph = copy.deepcopy(s2_paragraph)
        s3_paragraph.xpath('.//a:blip', namespaces=NS)[0].set(q('r', 'embed'), 'rIdSupplementV83CaM')
        for node in s3_paragraph.xpath('.//wp:docPr|.//pic:cNvPr', namespaces=NS):
            node.set('id', str(max_id + 2)); node.set('name', 'Supplementary Figure S3 CaM reference validation')
        s3_rel = etree.SubElement(rels, q('pr', 'Relationship'))
        s3_rel.set('Id', 'rIdSupplementV83CaM')
        s3_rel.set('Type', 'http://schemas.openxmlformats.org/officeDocument/2006/relationships/image')
        s3_rel.set('Target', 'media/figureS3_CaM_v83.png')
        modified_parts['word/media/figureS3_CaM_v83.png'] = s3_path.read_bytes()
        figure_report['Figure S3'] = set_image_geometry(s3_paragraph, s3_path, 'Supplementary CaM MD references and frozen-model prediction validation', width_mm=164)
        s3_caption = new_text_paragraph(paragraphs[155], read_caption(args.caption_s3))
        keep_figure_with_caption(s3_paragraph, s3_caption)
        paragraphs[167].addprevious(s3_paragraph)
        paragraphs[167].addprevious(s3_caption)
"""
assert needle in s;s=s.replace(needle,addition+needle);(m/'build_v83_manuscript.py').write_text(s,encoding='utf-8')
# Give key-site labels enough horizontal room.
for name in ['model_story_body.py','compose_fig3.py']:
    p=v/'figure3'/name;t=p.read_text(encoding='utf-8').replace("('rmsf',139,(.5,1.62)","('rmsf',137,(.5,1.62)").replace("('hchi',169,(1.45,2.35)","('hchi',166,(1.45,2.35)").replace('a=axis(x,94,13,26);ys=[]','a=axis(x,94,17 if metric==\'rmsf\' else 16,26);ys=[]');p.write_text(t,encoding='utf-8')
print(json.dumps({'caption_characters':len(caption),'output':plan['output_docx'],'section33_rewritten_around_model':True},ensure_ascii=False))