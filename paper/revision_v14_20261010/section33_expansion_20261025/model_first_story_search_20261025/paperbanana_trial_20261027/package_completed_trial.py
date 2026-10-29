from pathlib import Path
from PIL import Image
import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone

base = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
trial = Path(__file__).resolve().parent
run = trial / 'polish_v98_01'
source = trial.parent / 'revision_v98_structure_fill' / 'figure3' / 'fig3_dual_case_entropy.png'
word = base / 'paper_v98_20261027_中文修订_结构区域放大.docx'
generated = run / 'paperbanana_polished_v98.png'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

assert digest(source) == '75733ef569e662e0861ccf15a499e945700033b5e4ec31d601674029b35379a1'
assert digest(word) == 'eb1bc676be591345f29ab45c2b4f8f8ef9b94ea2cba6c64d74d0781d41ad4e92'
assert digest(generated) == 'b9a096d6690cd50bff347c51800f135b06733d2c8262747e673e234365ffdee3'
with Image.open(source) as im:
    source_size = im.size
    im.verify()
with Image.open(generated) as im:
    generated_size = im.size
    im.verify()
height_change = (generated_size[1] / generated_size[0]) / (source_size[1] / source_size[0]) - 1

preview = trial / 'comparison'
preview.mkdir(exist_ok=True)
shutil.copy2(source, preview / 'original_v98.png')
shutil.copy2(generated, preview / 'paperbanana_polished_v98.png')
release = base / 'figure3_paperbanana_20261027_润色对照稿.png'
shutil.copy2(generated, release)

html_template = '''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>图3：v98 与 PaperBanana 对照</title>
<style>body{font:17px/1.7 Arial,"Microsoft YaHei",sans-serif;margin:28px;color:#151515;background:#fff}h1{font-size:24px;margin:0 0 12px}p{max-width:1100px}main{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:28px;align-items:start}figure{margin:0}figcaption{font-weight:700;font-size:20px;margin:0 0 12px}img{display:block;width:100%;height:auto}a{color:#2166ac}@media(max-width:760px){main{grid-template-columns:1fr}}</style></head>
<body><h1>图3：v98 与 PaperBanana 实际生成稿</h1>
<p>两图按相同宽度展示，点击可打开原始 PNG。PaperBanana 使用作者原生 PolishAgent 的“视觉建议 → 图像润色”流程。生成稿仅作设计对照；布局未明显改善，结构和数据标记存在重绘差异，未替换论文。</p>
<main><figure><figcaption>原图 v98</figcaption><a href="ORIGINAL"><img src="ORIGINAL" alt="原图 v98"></a></figure>
<figure><figcaption>PaperBanana 润色稿</figcaption><a href="GENERATED"><img src="GENERATED" alt="PaperBanana 润色稿"></a></figure></main></body></html>'''
(preview / 'comparison.html').write_text(html_template.replace('ORIGINAL','original_v98.png').replace('GENERATED','paperbanana_polished_v98.png'), encoding='utf-8')
source_release = base / 'figure3_v98_20261027_结构区域放大.png'
assert digest(source_release) == digest(source)
comparison_release = base / 'figure3_paperbanana_20261027_原图对照.html'
comparison_release.write_text(html_template.replace('ORIGINAL',source_release.name).replace('GENERATED',release.name), encoding='utf-8')

report = f'''# PaperBanana 图3实际试跑结果

已于 2026-10-27 完成真实调用并生成图片。使用本机安装的作者仓库 PaperBanana，提交 f07f5f8a04d3025709d9ee3ec41d3d325a2e8529；入口为 agents.polish_agent.PolishAgent。此次是原生统计图润色流程（视觉模型提出建议，再由图像模型参照原图重绘），没有运行完整检索、规划、评论多智能体流水线。

输入为 v98 图3 PNG，另附针对现有六个面板的布局与科学约束。主模型为 gemini-3.1-pro-preview，图像模型为 gemini-3.1-flash-image-preview。密钥由本机配置读取，没有写入本结果包。

## 视觉结论

没有明显优于 v98。结构的阴影和点线的视觉权重增强，但六个面板的组织、标题层级和主要留白基本延续原图，尚未解决结构区空和整体不够紧凑的问题。生成稿为 {generated_size[0]}×{generated_size[1]}，原图为 {source_size[0]}×{source_size[1]}；按相同宽度展示，生成稿高约 {height_change:.1%}，并未缩紧版面。

## 科学图使用限制

这次输入是图片，生成过程没有直接使用 CSV 或 PDB 来重新绘制。c 曲线峰谷、散点位置以及结构丝带轮廓可见重新描画；不能将其称为数据保真的美容处理，也不能声称原有 853 个绘图数值被逐点保持。部分希腊字母的字形也需核对，尤其相关系数 rho 与拉丁 p 的区分。原生 CSV、脚本和结构才是正式成图依据。此次生成 PNG 仅保留作布局试样，未替换 Word 或 v98 原图。

## 文件与验证

- comparison/comparison.html：同宽并排对照；点击图片打开原始 PNG。
- polish_v98_01/paperbanana_polished_v98.png：实际生成稿。
- polish_v98_01/suggestions.txt：实际模型返回的十条建议。
- polish_v98_01/provenance.json：调用入口、模型、时间及输入哈希。
- run_paperbanana_polish.py：实际调用脚本；密钥不在脚本内。

输入 PNG SHA256：{digest(source)}

输出 PNG SHA256：{digest(generated)}

Word SHA256：{digest(word)}

包装时重新验证了输入 PNG 与 Word 的哈希，两者均未变化。未运行额外图像生成或模型预测。
'''
(trial / 'README_使用入口.md').write_text(report, encoding='utf-8')
(trial / 'comparison_review_zh.md').write_text(report, encoding='utf-8')
status = {
    'status': 'complete',
    'image_generated': True,
    'workflow': 'author native PolishAgent: VLM suggestions -> image refinement',
    'output': str(generated),
    'release_output': str(release),
    'comparison_html': str(comparison_release),
    'source_size': source_size,
    'generated_size': generated_size,
    'same_width_height_change_fraction': height_change,
    'source_sha256': digest(source),
    'output_sha256': digest(generated),
    'word_sha256': digest(word),
    'source_and_word_unchanged': True,
    'manuscript_replaced': False,
    'visual_verdict': 'No clear improvement over v98; keep as design trial only.',
    'scientific_fidelity': 'Raster image redraw; numerical and structural fidelity not certified.',
    'packaged_at': datetime.now(timezone.utc).isoformat(),
}
(trial / 'run_status.json').write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding='utf-8')

zip_path = base / 'figure3_paperbanana_20261027_实际结果与对照.zip'
files = [trial / 'README_使用入口.md', trial / 'comparison_review_zh.md', trial / 'run_status.json',
         trial / 'run_paperbanana_polish.py', trial / 'input_brief_zh.md', trial / 'paperbanana_content_en.txt',
         trial / 'input_source_manifest.json', run / 'status.json', run / 'provenance.json',
         run / 'suggestions.txt', run / 'result_metadata.json', generated,
         run / 'style_guides' / 'neurips2025_plot_style_guide.md',
         preview / 'comparison.html', preview / 'original_v98.png', preview / 'paperbanana_polished_v98.png']
for path in files:
    assert path.is_file(), path
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
    for path in files:
        zf.write(path, path.relative_to(trial).as_posix())
with zipfile.ZipFile(zip_path) as zf:
    assert zf.testzip() is None
    assert not any('credential' in item.filename.lower() or item.filename.endswith('.key') for item in zf.infolist())
print(json.dumps({'png': str(release), 'comparison': str(comparison_release), 'zip': str(zip_path),
                  'source_size': source_size, 'generated_size': generated_size,
                  'same_width_height_change_fraction': height_change, 'word_unchanged': True,
                  'zip_bytes': zip_path.stat().st_size}, ensure_ascii=False, indent=2))
