# CaM–CaMKI（1MXE）非 TCR 案例包

**用途**：回答“KIX 只定位了活跃区域，模型有没有真正使用伙伴信息”，并给出可替换图 3 a–d 的候选案例。  
**结论**：主推 CaM–CaMKI（1MXE / A001K）；若必须使用两个折叠蛋白，备选 FSH–FSHR（1XWD）。Cdt1–Geminin（2ZXX）更适合做界面幅值与侧链的补充，不适合做“位点级判别”的主证明。

## 先看什么

1. `README_案例建议.md`：完整论证、候选比较、图 3 改法、局限与审稿风险。
2. `fig_cam_camki_case_v4.png`：最新两态组图（a 游离→结合、b 两态曲线、c 界面误差、d 两态均值、e 对照）。
3. `caption_zh.md`：图注草稿与面板数据对应。
4. `draft_section33_zh.md`：3.3 节替换段落草稿。

## 关键数字

| 项目 | MD | 完整模型 | 完整模型去掉伙伴 | 链内模型 |
|---|---:|---:|---:|---:|
| 肽段平均 RMSF（Å） | 0.789 | 0.856 | 4.219 | 2.003 |
| 肽段 MAE（Å） | — | 0.099 | 3.430 | 1.214 |
| 肽段排序相关 | — | 0.806 | 0.808 | 0.505 |

- 主五折 OOF：整体 RMSF MAE 0.265 / 0.883，界面 0.157 / 0.819，肽段 0.099 / 1.214（Full / Intra）。
- 家族留出 OOF：界面 RMSF ρ 0.793 / −0.372，肽段 MAE 0.115 / 2.383。
- 伙伴移除对照：1MXE +3.36 Å；1XWD +0.08 Å；3QT2 +0.04 Å；2ZXX +0.23 Å；短链 3WWQ/2HWN/3P8M/6A7V 为 +1.19 到 +4.73 Å（必须如实披露）。

## 文件

| 文件 | 内容 |
|---|---|
| `cam_camki_per_residue.csv` | 169 个位点的 MD/Full/Intra/去伙伴 RMSF 与 χ，含界面标记和距离 |
| `case_metrics.csv` | 整体/界面/界面外/肽段/CaM 的 MAE、ρ、增益，主五折与家族留出 |
| `peptide_ablation_summary.csv` | 肽段在四种条件下的均值、MAE、ρ |
| `partner_removal_controls.csv` | 1MXE 与七个对照体系的伙伴移除 Δmean |
| `cam_structure.png` | 1MXE 的 PyMOL cartoon（链 A 灰、链 B 橙） |
| `render_structure_pymol.py` / `crop_structure.py` | 结构面板的 PyMOL 渲染与方形裁切 |
| `render_structure.py` / `3Dmol-min.js` | 早期 3Dmol 结构渲染，保留作 fallback |
| `fig_cam_camki_case_v4.{png,pdf,svg}` | 最新两态组图（a 游离→结合、b 两态曲线、c 界面误差、d 两态均值、e 对照） |
| `fig_cam_camki_case_v3.{png,pdf,svg}` | 早期结合态版本（a 只有结合态结构） |
| `make_case_figure_v4.py` | v4 两态组图脚本 |
| `make_case_figure_v3.py` | v3 结合态组图脚本 |
| `fig_cam_camki_case.{png,pdf,svg}` | 早期 CA 示意版本，保留作对照 |
| `make_case_figure.py` | 早期 CA 示意组图脚本 |
| `reproduce_case.py` | 复现主模型预测与伙伴移除消融 |
| `reproduction_check.json` | 复现结果：169 行全匹配，最大绝对差 2.5×10⁻⁶ |
| `qa_checks.py` / `QA_report.md` | 表内数字、界面计数、图值与消融表的独立核对 |
| `case_provenance.json` | 模型、checkpoint、数据路径、口径与局限 |
| `caption_zh.md` / `draft_section33_zh.md` | 图注与正文草稿 |

## 复现

```powershell
$env:PYTHONPATH='E:\Project\FlexPPI\src'
& 'E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\runtime_py312\Scripts\python.exe' `
  'E:\Project\FlexPPI\paper\revision_v14_20261010\non_tcr_case_cam_camki_20261025\reproduce_case.py'
& 'E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\runtime_py312\Scripts\python.exe' `
  'E:\Project\FlexPPI\paper\revision_v14_20261010\non_tcr_case_cam_camki_20261025\qa_checks.py'
& 'E:\Project\FlexPPI\tmp\mamba_root\envs\pymol\Scripts\pymol.bat' -cq `
  'E:\Project\FlexPPI\paper\revision_v14_20261010\non_tcr_case_cam_camki_20261025\render_structure_pymol.py'
& 'E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\runtime_py312\Scripts\python.exe' `
  'E:\Project\FlexPPI\paper\revision_v14_20261010\non_tcr_case_cam_camki_20261025\crop_structure.py'
& 'E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\runtime_py312\Scripts\python.exe' `
  'E:\Project\FlexPPI\paper\revision_v14_20261010\non_tcr_case_cam_camki_20261025\make_case_figure_v4.py'
```

主 checkpoint：`experiments/core_corrected_pipeline/trained_checkpoints/full/fold_4/epoch_040.pt`，SHA-256 `1cb429fc847dc436e542fe7926cbc4ac266a74cd705c8228415e25981a8fef74`。
