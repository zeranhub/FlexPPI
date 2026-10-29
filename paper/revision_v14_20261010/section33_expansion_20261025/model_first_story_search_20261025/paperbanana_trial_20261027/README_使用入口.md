# PaperBanana 图3实际试跑结果

已于 2026-10-27 完成真实调用并生成图片。使用本机安装的作者仓库 PaperBanana，提交 f07f5f8a04d3025709d9ee3ec41d3d325a2e8529；入口为 agents.polish_agent.PolishAgent。此次是原生统计图润色流程（视觉模型提出建议，再由图像模型参照原图重绘），没有运行完整检索、规划、评论多智能体流水线。

输入为 v98 图3 PNG，另附针对现有六个面板的布局与科学约束。主模型为 gemini-3.1-pro-preview，图像模型为 gemini-3.1-flash-image-preview。密钥由本机配置读取，没有写入本结果包。

## 视觉结论

没有明显优于 v98。结构的阴影和点线的视觉权重增强，但六个面板的组织、标题层级和主要留白基本延续原图，尚未解决结构区空和整体不够紧凑的问题。生成稿为 1856×2304，原图为 3171×3543；按相同宽度展示，生成稿高约 11.1%，并未缩紧版面。

## 科学图使用限制

这次输入是图片，生成过程没有直接使用 CSV 或 PDB 来重新绘制。c 曲线峰谷、散点位置以及结构丝带轮廓可见重新描画；不能将其称为数据保真的美容处理，也不能声称原有 853 个绘图数值被逐点保持。部分希腊字母的字形也需核对，尤其相关系数 rho 与拉丁 p 的区分。原生 CSV、脚本和结构才是正式成图依据。此次生成 PNG 仅保留作布局试样，未替换 Word 或 v98 原图。

## 文件与验证

- comparison/comparison.html：同宽并排对照；点击图片打开原始 PNG。
- polish_v98_01/paperbanana_polished_v98.png：实际生成稿。
- polish_v98_01/suggestions.txt：实际模型返回的十条建议。
- polish_v98_01/provenance.json：调用入口、模型、时间及输入哈希。
- run_paperbanana_polish.py：实际调用脚本；密钥不在脚本内。

输入 PNG SHA256：75733ef569e662e0861ccf15a499e945700033b5e4ec31d601674029b35379a1

输出 PNG SHA256：b9a096d6690cd50bff347c51800f135b06733d2c8262747e673e234365ffdee3

Word SHA256：eb1bc676be591345f29ab45c2b4f8f8ef9b94ea2cba6c64d74d0781d41ad4e92

包装时重新验证了输入 PNG 与 Word 的哈希，两者均未变化。未运行额外图像生成或模型预测。
