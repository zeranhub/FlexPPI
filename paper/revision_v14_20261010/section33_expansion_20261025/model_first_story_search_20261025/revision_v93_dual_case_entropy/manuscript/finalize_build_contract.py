from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
cfg=json.loads((HERE/'final_text_contract.json').read_text(encoding='utf-8'))
cfg['paragraph_replacements']=json.loads((HERE/'global_paragraph_replacements_proposed.json').read_text(encoding='utf-8'))
cfg['paragraph_replacements']['163']='全169个RMSF有效位置的MAE为Full 0.189132 Å、Intra 0.793611 Å；完整25残基靶肽为0.114859和2.383203 Å。靶肽21个有效侧链的Hχ MAE分别为0.254418和0.422994 nats；完整体系145个有效侧链的Hχ MAE为0.193119和0.228446 nats，原值Spearman相关为0.928208和0.932082。图S3保留全部靶肽位置及全体系的参考与预测，Full在N308的Hχ单点误差大于Intra，W303也存在高估，不能据此声称所有局部位置或统计口径均由Full占优。1MXE没有激酶催化域，也没有匹配的游离靶肽轨迹；1A06仅作为已有自抑制结构的文献背景，未作模型预测或MD。正文对调控构型的功能解释结合文献[29,30]与本例柔性分布提出，尚未直接计量片段回接触或激活过程，也不将Intra视作真实游离态动力学。所提供的1MXE输入结构与对应蛋白轨迹不含Ca离子，本例不用于量化钙效应。'
cfg['supplement4_caption']=Path(cfg.pop('supplement4_caption_file')).read_text(encoding='utf-8').strip()
cfg['figure3_width_mm']=164
cfg['supplement4_width_mm']=164
if len(cfg['section33']['before_figure'])==5 and not cfg['section33'].get('after_figure'):
    cfg['section33']['after_figure']=cfg['section33']['before_figure'][3:]
    cfg['section33']['before_figure']=cfg['section33']['before_figure'][:3]
    cfg['layout_note']='Move the two ubiquitin paragraphs after the full-page Figure 3 to avoid a nearly empty page before the figure; wording unchanged.'
cfg['output_docx']=str(BASE/'paper_v93_20261026_中文修订_双案例与构象占据整合.docx')
(HERE/'revision_content.json').write_text(json.dumps(cfg,ensure_ascii=False,indent=2),encoding='utf-8')
print('Final revision contract prepared; no Word file written.')
