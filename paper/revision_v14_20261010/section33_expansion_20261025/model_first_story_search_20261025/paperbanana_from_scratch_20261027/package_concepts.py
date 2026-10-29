from pathlib import Path
from PIL import Image
from datetime import datetime, timezone
from hashlib import sha256
import json
import shutil
import zipfile

root = Path(__file__).resolve().parent
run = root / 'generation_01'
base = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
source_image = root.parent / 'revision_v98_structure_fill' / 'figure3' / 'fig3_dual_case_entropy.png'
word = base / 'paper_v98_20261027_中文修订_结构区域放大.docx'
assert sha256(source_image.read_bytes()).hexdigest() == '75733ef569e662e0861ccf15a499e945700033b5e4ec31d601674029b35379a1'
assert sha256(word.read_bytes()).hexdigest() == 'eb1bc676be591345f29ab45c2b4f8f8ef9b94ea2cba6c64d74d0781d41ad4e92'
items = []
for number, label in [(1,'A'),(2,'B')]:
    image_path = run / f'candidate_{number}_final.png'
    metadata_path = run / f'candidate_{number}_pipeline.json'
    assert image_path.is_file() and metadata_path.is_file()
    metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
    with Image.open(image_path) as image:
        dimensions = image.size
        image.verify()
    final_field = metadata.get('eval_image_field')
    assert final_field and metadata.get(final_field), 'No author final image field'
    release = base / f'figure3_paperbanana_20261027_从头生成方案{label}.png'
    shutil.copy2(image_path, release)
    items.append({'candidate':label,'filename':release.name,'output':str(release),'dimensions':dimensions,
        'sha256':sha256(release.read_bytes()).hexdigest(),'final_field':final_field,
        'reference_ids':metadata.get('top10_references',[]),
        'critic_fields':[k for k in metadata if 'critic_suggestions' in k]})

html = '''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>PaperBanana 从头生成图3</title><style>body{font:17px/1.7 Arial,"Microsoft YaHei",sans-serif;color:#181818;background:white;margin:28px}h1{font-size:24px;margin:0 0 12px}p{max-width:1150px}main{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:28px;align-items:start}figure{margin:0}figcaption{font-size:20px;font-weight:700;margin-bottom:12px}img{width:100%;height:auto;display:block}a{color:#2166ac}@media(max-width:760px){main{grid-template-columns:1fr}}</style></head><body><h1>PaperBanana：从研究内容重新生成图3</h1><p>完整流程：参考检索 → 布局规划 → 风格设计 → 生成 → 最多两轮审阅修正。两套方案均未输入旧组图，按相同宽度展示，点击打开原始图片。这里展示的是独立生成的设计概念，科学成图需用真实结构和源数据脚本实现。现有论文未替换。</p><main>CARDS</main></body></html>'''
cards = ''.join(f'<figure><figcaption>方案{item["candidate"]}</figcaption><a href="{item["filename"]}"><img src="{item["filename"]}" alt="从头生成方案{item["candidate"]}"></a></figure>' for item in items)
gallery = base / 'figure3_paperbanana_20261027_从头生成对照.html'
gallery.write_text(html.replace('CARDS',cards),encoding='utf-8')
report = {'status':'complete','workflow':'author demo_full: Retriever -> Planner -> Stylist -> Visualizer -> Critic loop',
    'old_figure_supplied':False,'polish_agent_executed':False,'reference_pool_count':298,
    'max_critic_rounds':2,'requested_resolution':'4K','candidates':items,
    'output_kind':'generative layout concepts; numerical/structural fidelity not certified',
    'existing_v98_and_word_unchanged':True,'gallery':str(gallery),
    'finished_at':datetime.now(timezone.utc).isoformat()}
(root / 'run_status.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
readme = '''# PaperBanana 从头生成图3

本次实际执行作者的 demo_full 完整流程：Retriever、Planner、Stylist、Visualizer、Critic。输入由研究内容、853个冻结绘图数值的参考清单、科学解释边界、配色及期刊风格组成，没有输入旧完整组图，也没有调用现有图片润色流程。参考库为作者官方 PaperBananaBench 的298个diagram例子；自动Retriever按原逻辑从前200例中选择，两个候选共享参考检索，之后独立规划与生成，每候选最多两轮Critic。

run_status.json记录两个最终输出及哈希；generation_01/provenance.json记录原生入口和模型；candidate_N_pipeline.json保存各阶段描述、参考id、批评意见、图像文件索引。视觉概念由PaperBanana生成，图中的分子形状和数据标记不能直接作为科学证据；后续正式成图须换为真实结构渲染与原生源数据绘图。本次不修改Word或v98原图。

run_from_scratch.py是实验包装器，调用作者原生Agent和主处理器。它与已有安全包装器的区别仅为独立work_dir/期刊风格、原始参考图的绝对路径、两个并发独立设计意图、来源记录。作者Agent和内置提示词未修改。密钥通过本机配置读取，不包含在结果包。

scientific_content.txt为实际科学内容输入，visual_intent.txt为公共图意图，numeric_inventory.json为853数值清单。input_manifest.json保存源文件哈希。data/PaperBananaBench/diagram/ref.json仅用于本机运行，引用本机官方参考图绝对路径，结果包不附带第三方参考图片。图片预览按真实纵横比，不拉伸。
'''
(root / 'README_从头生成.md').write_text(readme,encoding='utf-8')
zip_path = base / 'figure3_paperbanana_20261027_从头生成结果包.zip'
with zipfile.ZipFile(zip_path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
    archive.write(gallery,gallery.name)
    for item in items:
        archive.write(item['output'],item['filename'])
    for filename in ['README_从头生成.md','run_status.json','scientific_content.txt','visual_intent.txt',
                     'numeric_inventory.json','input_manifest.json','prepare_from_scratch.py','run_from_scratch.py']:
        archive.write(root/filename,filename)
    archive.write(root/'style_guides'/'neurips2025_diagram_style_guide.md','style_guides/neurips2025_diagram_style_guide.md')
    # Full stage images remain in the local run; the portable review package
    # contains the two final concepts plus all textual stage records.
    for path in sorted(run.rglob('*')):
        if path.is_file() and path.suffix == '.json':
            archive.write(path,path.relative_to(root).as_posix())
    review = root / 'visual_review_zh.md'
    if review.is_file():
        archive.write(review,review.name)
with zipfile.ZipFile(zip_path) as archive:
    assert archive.testzip() is None
    assert not any('credential' in item.filename.lower() or item.filename.endswith('.key') for item in archive.infolist())
print(json.dumps({'gallery':str(gallery),'outputs':items,'zip':str(zip_path),'zip_bytes':zip_path.stat().st_size},ensure_ascii=False,indent=2))
