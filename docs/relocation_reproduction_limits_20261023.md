# 目录迁移后的复现范围与冻结校验限制

核验日期：2026-10-23。核验对象是 `E:/Project/FlexPPI` 中按用户指定 v58 Word 整理后的项目。本报告记录只读验收和已经确认的复现限制；没有运行训练、模型推理、分子动力学或科学统计，也没有修改科研脚本、冻结计划、原有哈希记录或既有数值结果。十二个实验入口的 `--check` 通过，只表示登记的脚本及关键本地输入存在，不能写成“全部实验复现通过”。

## 最终论文与图稿的独立核验

最终论文为 `paper/revision_v14_20261010/paper_v58_20261022_中文修订_图2至5配色统一_三线表.docx`，实测 SHA-256 为 `2058acaf9dae807bef357338a1017541f8a4a4fc85583a2e745c36ac3f38cb2a`，与 `paper/figures_v58/manifest.json` 和清理计划记录完全一致。独立读取 Word 内嵌图片后，图 1、2、3、4、5、S1 的对应图像均与图稿包登记的 Word 图像一致；S1 另外保存的不同编码 PNG，经解码比较全部像素一致。六套图的 18 个接受导出文件、151 个源文件、6 个可编辑文件，以及图稿根目录登记的 3 个文件，都通过存在性及 SHA-256 核对。这一结论不依赖对历史生产任务重新计算。

`experiments/steps.json` 的十二个标题均能在实际 Word 正文或补充材料标题中找到，包含正文 3.1—3.5 与补充 A—G。所有十二个 Python 入口以禁止生成 `.pyc` 的方式运行 `--check`，均返回成功；没有执行 `--step`。迁移记录 `docs/path_migration_20261023.json` 中的 119 个改动脚本全部可以解析，当前文件哈希均与记录的 `after_sha256` 相等。语法和文件存在性检查不能检出下文的冻结守卫冲突。

## 已确认的冻结代码冲突

下表路径均相对于 `experiments/_sources/`。这些问题来自代码文件中的目录字符串发生变化，而旧冻结记录继续保存最初执行时的文件哈希。哈希不一致不表示已有实验数值被改动；它表示当前迁移版代码不能继续声称与原来冻结的执行文件逐字节相同。下列守卫仍然保留，未被放宽。

| 生产脚本或校验程序 | 已确认冲突及相关冻结记录 |
| --- | --- |
| `label_sensitivity/run_label_sensitivity.py` | 第 110、326 行比较当前脚本与 `label_sensitivity/sensitivity_plan.json` 的 `analysis_code_sha256`。冻结值以 `bb1b02959a29` 开头，迁移后实际值以 `1f7ff8b8c360` 开头；迁移记录明确登记了这次变化。现有冻结计划不能原地继续执行。 |
| `energy_landscape_v21_20261012/trajectory_analysis/run_angle_distribution_analysis.py` | 第 47 行将上游 `run_label_sensitivity.py` 的当前哈希与旧 `sensitivity_plan.json` 的代码哈希比较，因此会受到上一项变化影响，即使本脚本自己的哈希仍与 `analysis_plan.json` 一致。第 82 行是本脚本自己的冻结检查，本次未发现这一项不一致。 |
| `model_argument_v24_20261013/experiment_predictive_occupancy/run_predictive_occupancy.py` | 第 114 行比较当前脚本与同目录 `analysis_plan.json` 的 `analysis_code_sha256`。冻结值以 `f6a527c7eba2` 开头，迁移后实际值以 `4e6a89d1a46d` 开头，现有冻结计划不能原地继续执行。 |
| `free_energy_curves_v27_20261014/case_curves/run_model_selected_ubiquitin.py` | 第 93—98 行保留原有执行修复的独立留痕：`analysis_plan.json` 冻结最初版本，`execution_repair.json` 记录后来的原有 QA 修复版本。第 97 行要求当前脚本匹配修复记录中的 `executed_code_sha256`。该记录值以 `af729da3a55b` 开头，迁移后实际值以 `5f5c356d10b6` 开头，因此原有修复记录也不能放行当前迁移版。不能把这项冲突简化为最初计划哈希不同，因为原项目本来就存在合法、明确记录的 QA 修复。 |
| `rich_angle_figure_v30_20261014/experiment/prepare_all_ubiquitin_chi.py` | 第 102 行比较当前脚本与同目录 `analysis_plan.json` 的代码哈希。冻结值以 `508230b2e243` 开头，迁移后实际值以 `b608a8d3edea` 开头。这是保留的逐角生产来源中的同类限制，未将其当成十二个入口的独立新增实验。 |

前四项直接涉及登记入口中的共同窗口标签、角度缓存、细箱构象占据和泛素模型选择曲线；部分生产程序会在到达代码哈希守卫之前先检查外部轨迹和其它冻结输入。这里没有调用它们来触发重计算，只通过读取实际守卫、冻结文件和当前文件哈希确认冲突。上述脚本的完整迁移前后哈希可以逐项对照 `docs/path_migration_20261023.json`，本报告没有替换这些记录。

## 冻结记录中原路径的限制

`energy_landscape_v21_20261012/trajectory_analysis/analysis_plan.json` 仍按来源记录保存 `paper/revision_v14_20261010/experiment_completion/label_sensitivity/run_label_sensitivity.py` 等迁移前路径。其独立校验程序 `verify_angle_distribution_results.py` 第 41—42 行直接读取 `dependency_hashes[*].path` 并检查哈希，第 45—46 行直接读取并校验 `inherited_plan_path`，第 43—44 行还读取原输出表路径。保留原路径是保留历史来源，不是当前已存在文件的保证；对应材料已移动到 `experiments/_sources/`。即使将路径在运行时映射到新位置，上游敏感性脚本的代码哈希仍存在前述冲突，因此单纯改路径或换一个输出文件夹并不能保证整个历史冻结链继续通过。

与这些已确认冲突不同，实测 `run_angle_distribution_analysis.py`、`free_energy_model_v26_20261013/experiment/run_low_Fstar_width.py` 和 `free_energy_curves_v27_20261014/population_curves/run_population_width_curves.py` 自身的文件哈希仍匹配各自现存计划。不能把所有历史代码笼统描述为哈希失效，也不能由某一个自身哈希匹配推导其上游路径、输入或全部依赖都可用。

## 整理验收与后续使用边界

迁移后的 `experiments/_sources/`、`experiments/_case_validation_sources/`、`experiments/_legacy_paper_data/` 均存在。验收快照分别记录到 31,790、104、83 个普通文件；原 `paper/figures/case_validation_v13_20261009/` 和 `paper/data/` 已不存在。原 `paper/revision_v14_20261010/experiment_completion/` 仅剩被占用而尚未移动的案例伴随 Word：`opposite_affinity_rmsf_cases_20261022/局部柔性下降与相反亲和结果_案例分析_三线表_20261022.docx`，353,191 字节；不得将它误认为仍未迁移的实验代码目录。

本次目录整理保留原冻结计划、原始哈希和已完成数值产物，不改变科研守卫，也不引入迁移哈希等价放行。当前可以确认的是最终 Word 内容未变、最终六套图及其源文件完整、十二个标题入口与论文对应、登记路径存在。未来需要重新执行旧生产任务时，应作为独立复现工作，先明确外部轨迹、模型、科学环境和目录映射，再建立新的、独立记录的冻结链及输出目录；应继续保留旧计划和旧结果，不能覆盖或悄悄改写旧哈希来获得“通过”。最终六图的独立重绘包及已保存的冻结训练载荷不因此被宣称无效，但本次也没有把图稿核验扩大为全部科学实验的重新验证。

旧稿和其它拟删除中间产物的删除仍被自动审批策略阻止，不能宣称已经删除。本报告只新增验收记录，不执行清理。
