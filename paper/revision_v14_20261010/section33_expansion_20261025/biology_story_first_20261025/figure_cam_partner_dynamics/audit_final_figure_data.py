"""Check actual Matplotlib data artists against source tables and raw NMR-STAR.

This does not render molecular structures or change the figure composition.
It evaluates only the numeric panel functions with an in-memory Matplotlib
figure. Their matrix-writing function is intercepted and checked in memory.
"""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.container import ErrorbarContainer
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    data_module = load_module("cam_data_contract_audit", ROOT / "build_figure_data.py")
    plot_module = load_module("cam_plot_numeric_audit", ROOT / "compose_cam_figure.py")
    raw, provenance = data_module.raw_qa()
    contract = json.loads((ROOT / "numeric_audit.json").read_text(encoding="utf-8"))
    contract_outputs = []
    for entry in contract["outputs"]:
        path = Path(entry["path"])
        assert path.is_file()
        assert sha(path) == entry["sha256"], path
        contract_outputs.append(dict(path=str(path), sha256=sha(path), matches_data_contract=True))
    assert len(contract_outputs) == 13

    backbone = pd.read_csv(ROOT / "source_data/backbone_87_NH_wide.csv")
    methyl = pd.read_csv(ROOT / "source_data/sidechain_44_methyl_wide.csv")
    key = pd.read_csv(ROOT / "source_data/key_sites_matched_NH_methyl_long.csv")
    delta = pd.read_csv(ROOT / "source_data/sidechain_44_methyl_delta_smMLCK_minus_nNOS.csv")
    expected_partners = ["CaMKKalpha", "eNOS", "nNOS", "smMLCK"]
    expected_keys = [("L18", "CD1"), ("L39", "CD2"), ("L105", "CD1"), ("M109", "CE"), ("M144", "CE")]
    assert plot_module.PARTNERS == expected_partners
    assert plot_module.KEYS == expected_keys
    assert len(backbone) == 87 and len(methyl) == 44 and len(delta) == 44

    comparisons = 0
    for frame in (backbone, methyl):
        for row in frame.to_dict("records"):
            identity = (row["observable"], int(row["seq_id"]), row["amino_acid"], row["atom"])
            for partner in expected_partners:
                true_value, true_error = raw[identity][partner][0]
                assert row[partner + "_s2"] == true_value
                assert row[partner + "_fit_error"] == true_error
                comparisons += 1
    assert comparisons == 524

    # Check the actual image arrays produced by the two heatmap functions.
    matrix_writes = []
    real_savetxt = np.savetxt
    np.savetxt = lambda path, values, **kwargs: matrix_writes.append((str(path), np.array(values)))
    try:
        fig = plt.figure(figsize=(8, 10))
        plot_module.heatmaps(fig, backbone, methyl)
        numeric_images = [a.images[0] for a in fig.axes if a.images]
        assert len(numeric_images) == 2
        for artist, frame in zip(numeric_images, (backbone, methyl)):
            expected = frame[[p + "_s2" for p in expected_partners]].to_numpy().T
            np.testing.assert_array_equal(np.asarray(artist.get_array()), expected)
            assert tuple(artist.get_clim()) == (0, 1)
        assert len(matrix_writes) == 2
        for (filename, values), frame in zip(matrix_writes, (backbone, methyl)):
            expected = frame[[p + "_s2" for p in expected_partners]].to_numpy().T
            np.testing.assert_array_equal(values, expected)
            if Path(filename).exists():
                np.testing.assert_array_equal(np.loadtxt(filename, delimiter=","), expected)
        plt.close(fig)
    finally:
        np.savetxt = real_savetxt

    # Verify each selected point and the lower/upper errorbar segments read by
    # Matplotlib, not only a reimplementation of the plotting arithmetic.
    fig = plt.figure(figsize=(8, 10))
    plot_module.local_data(fig, key)
    selected_points = []
    assert len(fig.axes) == 2
    for axis, observable in zip(fig.axes, ("backbone_NH", "methyl_axis")):
        containers = [item for item in axis.containers if isinstance(item, ErrorbarContainer)]
        assert len(containers) == 10
        assert tuple(axis.get_ylim()) == (0, 1.05)
        assert [label.get_text() for label in axis.get_xticklabels()] == [r for r, a in expected_keys]
        for index, container in enumerate(containers):
            site_index, partner_index = divmod(index, 2)
            residue, atom = expected_keys[site_index]
            partner = ("nNOS", "smMLCK")[partner_index]
            position = int(residue[1:])
            aa = residue[0]
            raw_key = (observable, position, aa, "N" if observable == "backbone_NH" else atom)
            true_value, true_error = raw[raw_key][partner][0]
            value = float(container.lines[0].get_ydata()[0])
            segment = container.lines[2][0].get_segments()[0]
            expected_limits = np.array([true_value - true_error, true_value + true_error])
            assert value == true_value
            np.testing.assert_allclose(segment[:, 1], expected_limits, atol=1e-14, rtol=0)
            selected_points.append(dict(panel="e" if observable == "backbone_NH" else "f", observable=observable,
                                        residue=residue, atom=raw_key[3], partner=partner, s2=value,
                                        fit_error=true_error, actual_errorbar_y_limits=segment[:, 1].tolist()))
    plt.close(fig)

    # The bar artist heights must equal all 44 raw endpoint differences, with
    # the same probe order as the complete methyl heatmap.
    fig = plt.figure(figsize=(8, 10))
    plot_module.all_changes(fig, delta)
    assert len(fig.axes) == 1
    bar_heights = np.array([p.get_height() for p in fig.axes[0].patches])
    expected_delta = methyl.smMLCK_s2.to_numpy() - methyl.nNOS_s2.to_numpy()
    assert list(delta.probe_id) == list(methyl.probe_id)
    np.testing.assert_allclose(bar_heights, expected_delta, atol=1e-14, rtol=0)
    assert len(bar_heights) == 44
    assert int((bar_heights > 0).sum()) == 21
    assert int((bar_heights < 0).sum()) == 23
    footer = " ".join(t.get_text() for t in fig.texts)
    assert "21 increase / 23 decrease" in footer
    assert "0.546 → 0.549" in footer
    plt.close(fig)

    zh = (ROOT / "caption_zh.md").read_text(encoding="utf-8")
    en = (ROOT / "caption_en.md").read_text(encoding="utf-8")
    for required in ("**a，", "**b，", "**c、d，", "**e、f，", "**g，", "1NIW", "SeMet", "第一模型", "L105-CD2 未在 f", "21 根上升、23 根下降", "0.5459", "0.5486", "拟合误差", "置信区间", "尚未验证"):
        assert required in zh, required
    for required in ("**a,", "**b,", "**c,d,", "**e,f,", "**g,", "SeMet", "not validated FlexPPI output", "21 positive and 23 negative"):
        assert required in en, required

    # Structure metadata is recorded but atom mapping/rendering QA remains a
    # separate structural audit; numeric agreement does not validate a camera.
    structure_files = [p for p in (ROOT / "structure_data").glob("*") if p.suffix in {".json", ".csv"}]
    report = dict(status="passed", audit_scope="Numeric panel artists c-g and captions; structure geometry/image layout audited separately",
                  figure_script=dict(path=str(ROOT / "compose_cam_figure.py"), sha256=sha(ROOT / "compose_cam_figure.py")),
                  source_table_contract_count=13, source_tables=contract_outputs,
                  independent_raw_partner_probe_comparisons=comparisons,
                  raw_NMR_STAR_sources=provenance,
                  plotted_heatmap_shapes=[[4, 87], [4, 44]], heatmaps_shared_scale=[0, 1],
                  all_heatmap_values_match_raw_and_CSV=True,
                  selected_point_count=20, selected_points_and_actual_errorbars_match_raw=True,
                  selected_points=selected_points,
                  all_44_delta_bar_heights_match_raw_endpoint_differences=True,
                  delta_definition="smMLCK S2 minus nNOS S2", positive_bar_count=21, negative_bar_count=23,
                  nNOS_mean_methyl_s2=float(methyl.nNOS_s2.mean()), smMLCK_mean_methyl_s2=float(methyl.smMLCK_s2.mean()),
                  mean_delta_s2=float(expected_delta.mean()), median_absolute_delta_s2=float(np.median(np.abs(expected_delta))),
                  main_key_site_methyl_atoms=dict(expected_keys),
                  L105_CD2_retained_in_whole_set=True,
                  caption_panel_mapping="a four structures; b L39/L105/M144 local pairs; c 87 NH heatmap; d 44 methyl heatmap; e five NH sites; f five methyl sites; g all 44 methyl differences",
                  captions=[dict(path=str(ROOT / f), sha256=sha(ROOT / f)) for f in ("caption_zh.md", "caption_en.md")],
                  structure_metadata_files=[dict(path=str(p), sha256=sha(p)) for p in structure_files],
                  inferential_tests=[], model_predictions_used=[],
                  limitations=["Deposited fit errors are not confidence intervals or complete experimental uncertainties",
                               "44 methyl probes cover 34 residues, not 44 independent biological samples",
                               "Selected sites explain the whole-set pattern and are not a separate population",
                               "NMR S2 is not RMSF or Hchi; similar means do not establish equivalence",
                               "Static structural snapshots and construct differences do not prove a dynamical mechanism",
                               "This does not validate current FlexPPI predictions"])
    (ROOT / "final_numeric_and_caption_audit.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(dict(status="passed", raw_probe_comparisons=comparisons, actual_plotted_points=20,
                          actual_bar_count=44, source_table_contract_count=13, captions="a-g mapping updated")))


if __name__ == "__main__":
    main()
