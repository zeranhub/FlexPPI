# Figure 2 design contract

Core conclusion: FlexPPI gives useful numerical readouts and prioritizes highly dynamic sites beyond residue-type means, with descriptive gains represented across all three input sources.

The revision replaces the sparse 20-type gain points with one additional discrimination task and a compact source summary. The two top density panels retain all 1,786,309 valid observations and show four original-prediction conditional trends. They answer the continuous-agreement question without duplicating four density charts. The three PR panels distinguish RMSF retrieval, Hchi retrieval and Hchi retrieval within residue type. The bottom 4 x 9 AP matrix answers whether gains occur across ATLAS, DynaRepo and ITrace within the existing out-of-fold evaluation. It is not a leave-source-out experiment.

Layout: asymmetric mixed quantitative composition, 183 x 190 mm, top two 65 mm squares, middle three 46 x 38 mm PR panels and bottom 36-cell single-metric matrix. Source font floor is 9 pt and effective font floor at the 167.64 mm manuscript width is 8.2446 pt. Semantic colors and line styles share one key across all panels. This preserves useful reading space while reducing three formerly repeated legend strips and deleting the old sparse type points.

Asset scan: `KernelDensity/plot_CorrelationDensity1.R`, `AUROC/plot_AUROC2.py` and `heatmap/plot_composition.py`, including companion PNGs, were read and viewed. These assets use unweighted R KDE, ROC and raw-count heatmaps respectively; their statistical/data structures are incompatible with the record-weighted mass, macro PR and exact AP source matrix. Parameters are inherited rather than claiming native production execution. `inherited_v44_plot.py` records the complete prior plotting implementation. The three required style baseline blocks are retained verbatim before panel logic.

Readability risks addressed: top color-key tick/title proximity was removed by shortening both agreement headings; all three PR panels use one shared method key and a compact AP row; source headings have original record counts; all matrix cells carry exact three-decimal AP numbers and use one 0..1 scale. Source differences are descriptive. Different task positive definitions remain explicit in the caption and statistical report.

No downsampling, clipping, model retraining, inference, MD, score calibration, threshold optimization or post hoc selection occurs. Undefined empty reference bins stay undefined and exported; no interpolation bridges them.
