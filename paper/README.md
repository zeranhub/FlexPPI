# FlexPPI 论文目录

- `paper.tex`：论文主文件，负责组织章节。
- `preamble.tex`：宏包和统一排版设置。
- `sections/`：日常撰写的正文文件。
- `tables/`：复杂表格文件。
- `figures/`：论文图片。
- `references.bib`：参考文献数据库。

在 VS Code 中打开任意 `sections/*.tex` 文件都可以写正文；这些文件顶部的
`% !TeX root = ../paper.tex` 会让 LaTeX Workshop 编译完整论文。

编译时优先选择 `latexmk (xelatex)`，然后使用 LaTeX Workshop 的 PDF 查看命令
预览。当前本机的 TeX Live 2026 `latexmk` 启动包装器存在 Lua 路径错误；若遇到
同样错误，可在本目录连续执行两次下列命令完成正文、公式编号和交叉引用编译：

```powershell
xelatex -interaction=nonstopmode -halt-on-error paper.tex
xelatex -interaction=nonstopmode -halt-on-error paper.tex
```
