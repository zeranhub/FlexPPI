from pathlib import Path
import hashlib,json
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
DOC=BASE/'paper_v39_20261017_中文修订_读者表述与英文图表_三线表.docx'
EXPECTED='61c88902753e1b0a826f2f3750a3dbdb40e1e57a8dbf08c5d60fa8d559721f48'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(DOC)==EXPECTED
RENDER=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\render_v39')
notes={
1:'标题、摘要、关键词及引言逐行检查；无越界、丢字、异常重叠或页码冲突。',
2:'引言续文与2.1方法正文完整；段落、标题及页码正常。页下留白由下一页图1整体排版产生，不构成内容丢失。',
3:'图1所有英文标题、箭头、图例、轴标签完整；图注与正文互不重叠，图片及图注均在版心内。',
4:'原生公式(1)、(2)逐一检查；根号、求和上下标、分数、联合熵下标及公式编号完整，无截断或行内遮挡。',
5:'原生公式(3)、(4)及其行内符号逐一检查；指数、上下标、括号、分数和公式编号完整；长公式仍在版心内，正文无挤压。',
6:'2.5与2.6所有段落、英文方法名、数值与行内根号完整；标题无孤立，页底最后一行留有正常页码间距。',
7:'2.7方法段落、温度/时间/长度单位及图5/图6引用完整；上标负号及侧链角下标可辨，无文本或公式重叠。',
8:'2.8原生公式(5)、求和条件、分数及F*上标完整；后续密集统计定义排版正常，末段自然跨页，无裁切。',
}
pages=[]
for n in range(1,9):
    p=RENDER/f'page-{n}.png'
    pages.append({'page':n,'png':str(p),'png_sha256':sha(p),'inspection':'Original-resolution individual PNG viewed at 100%; no contact sheet used.','verdict':'PASS','necessary_issues':[],'notes_zh':notes[n]})
out={'document':str(DOC),'document_sha256':EXPECTED,'status':'PASS','reviewer':'english_figures_167','scope':'Individual visual review of pages 1–8, including English Figure1 and native equations (1)–(5).','all_8_pages_individually_inspected':True,'document_modified':False,'pages':pages,'necessary_issues':[]}
(Path(__file__).parent/'visual_pages_01_08.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
print('PASS: pages 1–8; no necessary visual fixes.')
