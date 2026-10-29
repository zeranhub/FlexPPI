"""Task-local EditaPlot/Origin renderer for the redesigned quantitative Figure 3.

This is intentionally a task script, not a claim that the public generic CLI
supports exact typography overrides. Default execution only validates the data
and writes a reviewable renderer plan. COM requires --execute, a source-bound
human approval record, and the successful EditaPlot Origin smoke report.

Official API references consulted (2026-10-27):
https://docs.originlab.com/originpro/classoriginpro_1_1graph_1_1GLayer.html
https://docs.originlab.com/originpro/classoriginpro_1_1graph_1_1Plot.html
https://docs.originlab.com/originpro/classoriginpro_1_1base_1_1Label.html
The physical units, native line styles, text placement and readback follow
EditaPlot's scientific_renderer.py, evidence_renderer.py and verify_utils.py.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPOSITORY = Path(r"C:\Users\DELL\.agent-reach\editaplot")
sys.path.insert(0, str(REPOSITORY / "runtime" / "src"))

FULL = "#2166AC"
INTRA = "#9A633D"
MD = "#222222"
MUTED = "#B8BCC1"
STYLE = {"font_family": "Arial", "title_pt": 10.0, "panel_pt": 11.0,
         "axis_pt": 9.0, "tick_pt": 9.0, "legend_pt": 9.0,
         "frame_pt": 0.6, "line_pt": 1.0, "full_color": FULL,
         "intra_color": INTRA, "full_line_style": 0, "intra_line_style": 1,
         "layer_factor": 1.0}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2,
                               allow_nan=False), encoding="utf-8")


def emit(event: str, **payload: Any) -> None:
    print(json.dumps({"event": event, **payload}, ensure_ascii=False), flush=True)


@dataclass
class Series:
    label: str
    x: list[float]
    y: list[float]
    kind: str = "s"
    color: str = FULL
    line_style: int = 0
    line_pt: float = 1.0
    symbol_pt: float = 3.2
    symbol_kind: int = 2
    symbol_interior: int = 0  # LabTalk -kf enumeration: 0 solid, 1 open.
    transparency: float = 0.0
    role: str = "measurement"


@dataclass
class Text:
    name: str
    text: str
    x_percent: float
    y_percent: float
    size_pt: float = 9.0
    bold: bool = False
    color: str = "#20262B"


@dataclass
class Panel:
    id: str
    letter: str
    title: str
    x_title: str
    y_title: str
    source_files: list[str]
    series: list[Series]
    texts: list[Text] = field(default_factory=list)
    xlim: tuple[float, float] | None = None
    ylim: tuple[float, float] | None = None
    width_cm: float = 8.6
    height_cm: float = 6.7
    left: float = 18.0
    top: float = 17.0
    width: float = 76.0
    height: float = 62.0
    y_tick_labels: bool = True
    legend: bool = True
    core_band: bool = False


def floats(values: Any) -> list[float]:
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError("Non-finite plotted numeric values")
    return array.tolist()


def scatter(label: str, x: Any, y: Any, **options: Any) -> Series:
    return Series(label, floats(x), floats(y), **options)


def source_path(name: str) -> Path:
    # Derived paired-analysis tables live beside the scientific proposal;
    # immutable copied measurement tables live in source_data.
    base = ROOT if name in {"f_paired_differences.csv", "f_paired_summary.csv"} else ROOT / "source_data"
    return base / name


def source_table(name: str) -> pd.DataFrame:
    return pd.read_csv(source_path(name))


def prepare_panels() -> list[Panel]:
    """Read source values; derive only declared coordinate/display helpers."""
    panels: list[Panel] = []
    cam_name = "cam_all25.csv"
    cam = source_table(cam_name)
    if len(cam) != 25:
        raise ValueError("Panel b requires every one of the 25 peptide positions")
    b_series = [
        scatter("MD", cam.CaMKI_native_position, cam.target_rmsf,
                kind="y", color=MD, line_style=2, symbol_interior=1),
        scatter("Full", cam.CaMKI_native_position, cam.pred_rmsf,
                kind="l", color=FULL, line_style=0),
        scatter("Intra", cam.CaMKI_native_position, cam.intra_rmsf,
                kind="l", color=INTRA, line_style=1),
    ]
    panels.append(Panel("b_rmsf_all25", "b", "CaMKI backbone flexibility",
                        "CaMKI residue", "RMSF (Å)", [cam_name], b_series,
                        width_cm=8.9, xlim=(293.5, 318.5), ylim=(0,8.2), core_band=True,
                        texts=[Text("MAEFull", "MAE: Full 0.115 Å", 58, 25),
                               Text("MAEIntra", "Intra 2.383 Å", 58, 32),
                               Text("Core", "W303–M316", 61, 44)]))
    panels.append(Panel("b_detail_md_full", "b", "MD + Full: all 25 residues",
                        "CaMKI residue", "RMSF (Å)", [cam_name], b_series[:2],
                        width_cm=8.9, xlim=(293.5, 318.5), ylim=(0.0, 2.3), core_band=True,
                        texts=[Text("Core", "W303–M316", 62, 26)]))

    c_name = "ubiquitin_72_position_profile_with_original_mask.csv"
    profile = source_table(c_name)
    mask = profile.valid_hchi.astype(str).str.lower().eq("true")
    valid = profile.loc[mask]
    if len(valid) != 66:
        raise ValueError("Panel c must retain exactly the 66 original valid positions")
    parity = scatter("Identity", [0, 5], [0, 5], kind="l", color=MUTED,
                     line_pt=0.6, line_style=1, role="mathematical_identity")
    points = scatter("Full", valid.MD_Hchi_nats, valid.Full_Hchi_nats,
                     symbol_pt=3.5)
    lo = -.10  # Margin keeps valid zero-valued prolines fully visible.
    hi = float(max(valid.MD_Hchi_nats.max(), valid.Full_Hchi_nats.max())) + .15
    parity.x = [lo, hi]
    parity.y = [lo, hi]
    panels.append(Panel("c_ubiquitin_parity", "c", "Ubiquitin side-chain entropy",
                        "MD Hχ (nats)", "Full Hχ (nats)", [c_name], [parity, points],
                        xlim=(lo, hi), ylim=(lo, hi), width_cm=8.9,
                        height_cm=8.9, left=22, top=17, width=68, height=68,
                        legend=False,
                        texts=[Text("SampleCount", "n = 66", 38, 23),
                               Text("MAE", "MAE = 0.080", 41, 29),
                               Text("Rho", "ρ = 0.979", 39, 35)]))
    c=panels[-1]
    for residue, label_x, label_y in [(40, 42, 43), (41, 40, 60), (62, 77, 52)]:
        row=valid.loc[valid.residue_id.eq(residue)].iloc[0]
        c.texts.append(Text(f"Q{residue}",f"Q{residue}",label_x,label_y))
        lx=lo+(label_x-c.left)/c.width*(hi-lo)
        ly=hi-(label_y-c.top)/c.height*(hi-lo)
        c.series.insert(1,scatter(f"Q{residue} leader",[lx,row.MD_Hchi_nats],[ly,row.Full_Hchi_nats],
            kind="l",color=MUTED,line_pt=.5,role="residue_annotation_leader"))

    d_name = "ubiquitin_all6_Gln_Hchi_SASA_plot.csv"
    gln = source_table(d_name).sort_values("MD_Hchi_nats", kind="stable")
    if len(gln) != 6:
        raise ValueError("Panel d must contain all six supplied Gln residues")
    ys = np.arange(len(gln), 0, -1, dtype=float)
    pair_series: list[Series] = []
    for index, row in enumerate(gln.itertuples(index=False)):
        pair_series.append(scatter(f"Pair_{row.residue}",
                                   [row.MD_Hchi_nats, row.Full_Hchi_nats],
                                   [ys[index], ys[index]], kind="l", color=MUTED,
                                   line_pt=.65, role="paired_connector"))
    pair_series += [scatter("MD", gln.MD_Hchi_nats, ys, color=MD, symbol_pt=4,
                            symbol_interior=1),
                    scatter("Full", gln.Full_Hchi_nats, ys, symbol_pt=4)]
    d = Panel("d_gln_pairs", "d", "Same-type side-chain environments",
              "Joint Hχ (nats)", "", [d_name], pair_series, width_cm=8.9,
              left=16, top=22, width=59, height=57, ylim=(.3, 6.7),
              y_tick_labels=False, legend=False)
    d.texts.extend([Text("MDKey", "○ MD", 33, 17, color=MD),
                    Text("FullKey", "● Full", 54, 17, color=FULL),
                    Text("SASAHeading", "SC SASA (Å²)", 87, 20, size_pt=9)])
    for index, row in enumerate(gln.itertuples(index=False)):
        y_page = d.top + d.height * (6.7 - ys[index]) / 6.4
        d.texts += [Text(f"Residue_{row.residue}", str(row.residue), 7, y_page),
                    Text(f"SASA_{row.residue}", f"{row.sidechain_SASA_A2:.1f}",
                         87, y_page)]
    panels.append(d)

    e_name = "defined_angles_long.csv"
    angles = source_table(e_name)
    corr = source_table("per_chi_correlations.csv")
    for number in range(1, 5):
        chi = f"chi{number}"
        part = angles.loc[angles.chi.eq(chi)]
        xcol = "predicted_Hchi_nats_AA_centered"
        ycol = "fine_entropy_four_origins_mean_nats_AA_centered"
        rho = corr.loc[corr.chi.eq(chi) & corr.predictor.eq("predicted_joint_Hchi")
                       & corr.metric.eq("fine_entropy"), "rho_AA_centered"]
        note = f"n = {len(part)}"
        if len(rho) == 1:
            note += f"; ρ = {float(rho.iloc[0]):.2f}"
        normal=part.loc[~part.residue.isin(["Q41","Q2"])]
        e_series=[scatter("X zero",[0,0],[-1.05,.80],kind="l",color=MUTED,line_pt=.5,role="zero_guide"),
                  scatter("Y zero",[-.75,.60],[0,0],kind="l",color=MUTED,line_pt=.5,role="zero_guide"),
                  scatter("Full", normal[xcol], normal[ycol],symbol_pt=3.2,transparency=20)]
        for residue,shape in [("Q41",5),("Q2",1)]:
            special=part.loc[part.residue.eq(residue)]
            if len(special):
                e_series.append(scatter(residue,special[xcol],special[ycol],color=MD,
                    symbol_kind=shape,symbol_pt=4.8,symbol_interior=1))
        panels.append(Panel(f"e_chi{number}", "e" if number == 1 else "",
                            f"χ{number}: within-type association",
                            "Joint Full Hχ\nAA-centered (nats)",
                            "MD single-χ entropy\nAA-centered (nats)",
                            [e_name, "per_chi_correlations.csv"],
                            e_series,
                            width_cm=8.9, height_cm=8.2, left=25, width=68,height=58,
                            xlim=(-.75,.60),ylim=(-1.05,.80),
                            legend=False, texts=[Text("Association", note, 51, 22),
                                Text("QKey", "◇ Q41    □ Q2" if number<4 else "No Gln χ4",57,96)]))

    f_name = "f_paired_differences.csv"
    fs_name = "f_paired_summary.csv"
    differences = source_table(f_name)
    summaries = source_table(fs_name)
    for metric, label in [("fine_entropy", "MD per-χ entropy"),
                          ("lowF_width", "Low-F* angular width")]:
        part = differences.loc[differences.metric_key.eq(metric)]
        summary = summaries.loc[summaries.metric_key.eq(metric)]
        if len(part) != 83 or len(summary) != 1:
            raise ValueError(f"Panel f {metric} needs 83 matched records and one summary")
        row = summary.iloc[0]
        # The deterministic ordinate only separates marks. It is not a value,
        # estimate, rank, error or extra observation; record identities remain.
        jitter_y = .78 + .18 * np.sin(np.arange(83) * 2.399963229728653)
        delta = part.delta_rho_full_hchi_minus_full_rmsf
        series = [scatter("Zero", [0, 0], [.05, 1.13], kind="l", color=MUTED,
                          line_style=1, line_pt=.6, role="mathematical_zero")]
        for dataset,shape in [("atlas",2),("dynarepo",1),("itrace",3)]:
            keep=part.dataset.eq(dataset).to_numpy()
            series.append(scatter(dataset,delta.loc[keep],jitter_y[keep],symbol_pt=3.0,
                transparency=55,symbol_kind=shape))
        series += [scatter("Bootstrap 95% CI", [row.bootstrap_95_low, row.bootstrap_95_high],
                          [.27, .27], kind="l", color=MD, line_pt=1.2,
                          role="supplied_bootstrap_interval"),
                  scatter("Paired mean", [row["mean"]], [.27], symbol_pt=5.5,
                          role="supplied_paired_mean")]
        xmin, xmax = -.18, .42  # Shared physical comparison range retains all 166 values.
        panels.append(Panel(f"f_{metric}", "f" if metric == "fine_entropy" else "",
                            label, "Δρ = Full Hχ − Full RMSF", "",
                            [f_name, fs_name], series, width_cm=8.9, height_cm=6.7,
                            left=13, width=82, top=17, height=55,
                            xlim=(xmin, xmax), ylim=(0, 1.2), y_tick_labels=False,
                            legend=False,
                            texts=[Text("RecordCount", "83 paired records", 55, 23),
                                   Text("IntervalMeaning", f"{row['mean']:.3f} [{row.bootstrap_95_low:.3f}, {row.bootstrap_95_high:.3f}]", 56, 66),
                                   Text("SourceKey", "○ ATLAS  □ DynaRepo  △ ITrace", 52, 94)]))
    return panels


def renderer_plan(panels: list[Panel]) -> dict[str, Any]:
    names = sorted({name for panel in panels for name in panel.source_files})
    return {"renderer": "task-local EditaPlot native Origin renderer v1",
            "style": STYLE, "source_sha256": {name: sha256(source_path(name))
                                               for name in names},
            "panels": [asdict(panel) for panel in panels],
            "display_helpers": ["identity line", "paired connectors",
                                "source-order row ordinates", "deterministic strip jitter"],
            "scientific_analysis_performed_by_renderer": False,
            "approval_status": "requires_external_user_approval_record"}


def plan_digest(plan: dict[str, Any]) -> str:
    payload = json.dumps(plan, sort_keys=True, ensure_ascii=False,
                         separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def require_approval(path: Path, digest: str, plan: dict[str, Any]) -> dict[str, Any]:
    approval = json.loads(path.read_text(encoding="utf-8"))
    if approval.get("confirmed") is not True or not approval.get("user_confirmation"):
        raise RuntimeError("A recorded human confirmation is required before COM")
    if approval.get("renderer_plan_sha256") != digest:
        raise RuntimeError("User approval does not match the current renderer plan")
    if approval.get("source_sha256") != plan["source_sha256"]:
        raise RuntimeError("Source hashes changed after approval")
    if set(approval.get("approved_panel_ids", [])) != {p["id"] for p in plan["panels"]}:
        raise RuntimeError("Approval does not cover this exact panel set")
    scientific = ROOT / "scientific_design_proposal.json"
    if approval.get("scientific_proposal_sha256") != sha256(scientific):
        raise RuntimeError("Approval does not bind the reviewed scientific proposal")
    return approval


def require_smoke(path: Path) -> dict[str, Any]:
    smoke = json.loads(path.read_text(encoding="utf-8"))
    if smoke.get("ok") is not True and smoke.get("status") not in {"passed", "success"}:
        raise RuntimeError("Successful EditaPlot Origin smoke report required")
    return smoke


def add_page_text(op: Any, graph: Any, layer: Any, text: Text) -> Any:
    # Reuse the verified editable page-label route, including clipping checks.
    from origin_sciplot.origin_backend.evidence_renderer import _add_shap_page_label
    label, _ = _add_shap_page_label(
        op, graph, layer, name=text.name, text=text.text,
        center_x_percent=text.x_percent, center_y_percent=text.y_percent,
        size_pt=text.size_pt, font_family="Arial", bold=text.bold, color=text.color)
    return label


def style_axis(op: Any, layer: Any, panel: Panel) -> None:
    font = int(round(op.lt_float("font(Arial)")))
    for axis in ("x", "y"):
        layer.set_int(f"{axis}.showAxes", 0 if axis=='y' and panel.id.startswith('f_') else 1)
        layer.set_int(f"{axis}.showGrids", 0)
        layer.set_int(f"{axis}.ticks", 2)
        layer.set_int(f"{axis}.showLabels", 1 if axis == "x" or panel.y_tick_labels else 0)
        layer.set_int(f"{axis}.showlabel", 1 if axis == "x" or panel.y_tick_labels else 0)
        layer.set_int(f"{axis}.label.show", 1 if axis == "x" or panel.y_tick_labels else 0)
        layer.set_int(f"{axis}.label.type", 1)
        layer.set_int(f"{axis}.minorTicks", 0)
        layer.set_int(f"{axis}.atZero", 0)
        layer.set_float(f"{axis}.label.rotate", 0)
        layer.set_float(f"{axis}.label.pt", 9)
        layer.set_int(f"{axis}.label.font", font)
        layer.set_int(f"{axis}.label.color", 1)
        layer.set_float(f"{axis}.thickness", .6)
        layer.set_float(f"{axis}.tickthickness", .6)
        layer.set_float(f"{axis}.mtickthickness", .6)
        layer.set_float(f"{axis}.ticklength", 3)
    if panel.xlim:
        layer.set_xlim(*panel.xlim)
    if panel.ylim:
        layer.set_ylim(*panel.ylim)
    tick_positions=None
    if panel.id.startswith('b_'):
        tick_positions='294 300 306 312 318'
    elif panel.id.startswith('c_'):
        tick_positions='0 1 2 3 4'
    elif panel.id.startswith('d_'):
        tick_positions='1.8 2.2 2.6 3.0'
    elif panel.id.startswith('e_'):
        tick_positions='-0.5 0 0.5'
    elif panel.id.startswith('f_'):
        tick_positions='-0.1 0 0.1 0.2 0.3 0.4'
    if tick_positions:
        layer.set_str('x.ticksbydata',tick_positions)
    layer.axis("x").title = panel.x_title
    layer.axis("y").title = panel.y_title
    for axis, name in [("x", "xb"), ("y", "yl")]:
        label = layer.label(name)
        if label is None:
            continue
        label.text = panel.x_title if axis == "x" else panel.y_title
        label.set_float("fsize", 9)
        label.set_int("font", font)
        label.set_int("bold", 0)
        label.set_int("color", 1)
        label.set_int("attach", 1)
        label.set_float("rotate", 90 if axis == "y" else 0)
        op.lt_exec("doc -uw;")
        page_width = float(op.lt_float("page.width"))
        page_height = float(op.lt_float("page.height"))
        label_width = float(label.get_float("width"))
        label_height = float(label.get_float("height"))
        if axis == "x":
            label.set_float("left", page_width * (panel.left + panel.width / 2) / 100 - label_width / 2)
            label.set_float("top", page_height * .855 - label_height / 2)
        else:
            label.set_float("left", page_width * .045 - label_width / 2)
            label.set_float("top", page_height * (panel.top + panel.height / 2) / 100 - label_height / 2)
    # Lock page/axis physical geometry and layer factor after all label binding.
    layer.set_int("fixed", 1)
    layer.set_float("factor", 1.0)
    op.lt_exec("doc -uw;")


def native_axis_state(op: Any, graph: Any, layer: Any, panel: Panel) -> dict[str, Any]:
    from origin_sciplot.origin_backend.verify_utils import (
        read_layer_geometry_percent, read_page_geometry_cm)
    font = int(round(op.lt_float("font(Arial)")))
    state: dict[str, Any] = {"page": read_page_geometry_cm(op, graph),
                             "layer": read_layer_geometry_percent(op, layer)}
    for axis in ("x", "y"):
        values = {prop: float(layer.get_float(f"{axis}.{prop}"))
                  for prop in ["from", "to", "label.pt", "label.font", "thickness",
                               "tickthickness", "label.rotate", "showAxes", "showLabels"]}
        if abs(values["label.pt"] - 9) > .05 or int(round(values["label.font"])) != font:
            raise RuntimeError(f"{panel.id}: {axis} typography readback differs")
        for prop in ("thickness", "tickthickness"):
            if abs(values[prop] - .6) > .05:
                raise RuntimeError(f"{panel.id}: {axis}.{prop} is not 0.6 pt")
        expected_show=0 if axis=='y' and panel.id.startswith('f_') else 1
        if int(round(values['showAxes']))!=expected_show or int(round(values['showLabels']))!=(1 if axis=='x' or panel.y_tick_labels else 0):
            raise RuntimeError(f"{panel.id}: axis line or label visibility incorrect: {values}")
        state[axis] = values
    if abs(float(layer.get_float("factor")) - 1) > .001:
        raise RuntimeError("Origin layer scale factor changed physical typography")
    if abs(state['page']['width_cm']-panel.width_cm)>.03 or abs(state['page']['height_cm']-panel.height_cm)>.03:
        raise RuntimeError("Saved page changed physical dimensions")
    if panel.core_band:
        core={key:float(op.lt_float(f"CoreBand.{key}")) for key in ['x','dx','fillcolor','order']}
        if abs(core['x']-309.5)>.05 or abs(core['dx']-14)>.05:
            raise RuntimeError(f"Core interval changed from native W303-M316: {core}")
        if int(core['fillcolor'])!=int(op.ocolor('#ECECEC')) or int(core['order'])!=256:
            raise RuntimeError("Core background color/drawing order changed")
        state['core_band']=core
    return state


def readback(op: Any, graph: Any, layer: Any, panel: Panel,
             plots: list[Any], labels: dict[str, Any]) -> dict[str, Any]:
    from origin_sciplot.origin_backend.verify_utils import (
        verify_text_fonts, verify_text_sizes, verify_plot_color,
        verify_plot_line_widths, verify_symbol_style)
    text_sizes = {"FigureTitle": 10.0}
    if panel.letter:
        text_sizes["PanelLetter"] = 11.0
    if panel.x_title:
        labels["XAxisTitle"] = layer.label("xb")
        text_sizes["XAxisTitle"] = 9.0
    if panel.y_title:
        labels["YAxisTitle"] = layer.label("yl")
        text_sizes["YAxisTitle"] = 9.0
    text_sizes.update({text.name: text.size_pt for text in panel.texts})
    for key in labels:
        if key.startswith("Legend"):
            text_sizes[key] = 9.0
    result = {"axis_geometry": native_axis_state(op, graph, layer, panel),
              "text_sizes": verify_text_sizes(labels, text_sizes),
              "text_fonts": verify_text_fonts(op, labels, "Arial"), "plots": []}
    for index, (spec, plot) in enumerate(zip(panel.series, plots, strict=True)):
        record = {"label": spec.label, "kind": spec.kind,
                  "source_values": {"x": spec.x, "y": spec.y},
                  "origin_range": plot.lt_range(),
                  "color": verify_plot_color(op, plot, spec.color,
                                               variable_name=f"__fig3_color_{index}")}
        if spec.kind in {"l", "y"}:
            record["width"] = verify_plot_line_widths(op, {spec.label: plot}, spec.line_pt)
            variable = f"__fig3_dash_{index}"
            plot.layer.LT_execute(f"{{range rr={plot.lt_range()};get rr -d {variable};}}")
            value = float(op.lt_float(variable))
            if abs(value - spec.line_style) > .05:
                raise RuntimeError(f"{panel.id}: {spec.label} line style changed")
            record["line_style"] = value
        if spec.kind in {"s", "y"}:
            record["symbol"] = verify_symbol_style(
                op, plot, expected_size_pt=spec.symbol_pt,
                expected_edge_percent=50.0, expected_symbol_interior=spec.symbol_interior,
                expected_symbol_kind=spec.symbol_kind)
        result["plots"].append(record)
    return result


def create_panel(op: Any, panel: Panel, output: Path) -> dict[str, Any]:
    from origin_sciplot.scientific_visual import AdaptiveOriginStyle
    from origin_sciplot.origin_backend.verify_utils import verify_page_and_layer, NativePlotReference
    from origin_sciplot.origin_backend.export_utils import export_graph
    op.new(asksave=False)  # Only inside the task-owned isolated Origin session.
    for index, source in enumerate(panel.source_files):
        sheet = op.new_sheet("w", f"Original source {index + 1}")
        sheet.from_df(source_table(source))
    plot_sheet = op.new_sheet("w", "Plot coordinates - declared display helpers")
    plot_sheet.name = "PlotData"
    graph = op.new_graph(panel.id, template="Line")
    graph.name = "Figure3"
    graph.set_int("background", op.ocolor("#FFFFFF"))
    layer = graph[0]
    style = AdaptiveOriginStyle(
        profile_name="FlexPPI Figure3 exact user specification",
        page_width_cm=panel.width_cm, page_height_cm=panel.height_cm,
        layer_left_percent=panel.left, layer_top_percent=panel.top,
        layer_width_percent=panel.width, layer_height_percent=panel.height,
        axis_title_size_pt=9, tick_label_size_pt=9, legend_size_pt=9,
        plot_line_width_pt=1, frame_line_width_pt=.6, layer_factor=1)
    verify_page_and_layer(graph, layer, origin=op, style=style)
    plots: list[Any] = []
    for index, item in enumerate(panel.series):
        plot_sheet.from_list(index * 2, item.x, lname=f"{item.label} X", axis="X")
        plot_sheet.from_list(index * 2 + 1, item.y, lname=item.label, axis="Y")
        plot = layer.add_plot(plot_sheet, index * 2 + 1, index * 2, type=item.kind)
        if plot is None:
            raise RuntimeError(f"Origin could not add {item.label}")
        plot.color = item.color
        plot.set_cmd(f"-w {round(item.line_pt * 500)}")
        if item.kind in {"l", "y"}:
            plot.set_cmd(f"-d {item.line_style}")
        if item.kind in {"s", "y"}:
            plot.symbol_kind = item.symbol_kind
            plot.symbol_size = item.symbol_pt
            plot.set_cmd(f"-kf {item.symbol_interior}", "-kh 50")
        if item.transparency:
            plot.transparency = item.transparency
        plots.append(NativePlotReference(f"[{graph.name}]1!{index+1}",layer.obj))
    layer.rescale()
    style_axis(op, layer, panel)
    if panel.core_band:
        y0,y1=panel.ylim
        command=(f'draw -n CoreBand -b;CoreBand.attach=2;CoreBand.dx=14;CoreBand.dy={y1-y0};'
                 f'CoreBand.x=309.5;CoreBand.y={(y0+y1)/2};'
                 'CoreBand.fillcolor=color("#ECECEC");CoreBand.color=color("#ECECEC");'
                 'CoreBand.order=256;doc -uw;')
        if not op.lt_exec(command):
            raise RuntimeError("Could not create editable structural core interval")
    for name in ("legend", "Legend"):
        obj = layer.label(name)
        if obj is not None:
            obj.remove()
    labels = {"FigureTitle": add_page_text(
        op, graph, layer, Text("FigureTitle", panel.title, 55, 7, 10))}
    if panel.letter:
        labels["PanelLetter"] = add_page_text(
            op, graph, layer, Text("PanelLetter", panel.letter, 4, 7, 11, True))
    for text in panel.texts:
        labels[text.name] = add_page_text(op, graph, layer, text)
    if panel.legend:
        for index, item in enumerate(panel.series):
            if item.role != "measurement":
                continue
            name = f"Legend{index}"
            glyph = "○" if item.label=="MD" else ("—" if item.line_style == 0 else "- -")
            legend_x,legend_y=(33,24+8*index) if panel.id.startswith('b_') else (32+24*index,94)
            labels[name] = add_page_text(op, graph, layer, Text(
                name, f"{glyph} {item.label}", legend_x,legend_y,9,color=item.color))
    op.lt_exec("doc -uw;")
    style_axis(op, layer, panel)
    before = readback(op, graph, layer, panel, plots, labels)
    workbook_columns = plot_sheet.to_df().columns.tolist()
    plot_sheet_range = plot_sheet.lt_range()
    saved_label_names = {name: label.name for name, label in labels.items()}
    numeric_before = [plot_sheet.to_list(i) for i in range(len(panel.series) * 2)]
    for index,spec in enumerate(panel.series):
        for offset,expected in enumerate([spec.x,spec.y]):
            actual=np.asarray(numeric_before[2*index+offset],dtype=float)
            if not np.array_equal(actual[:len(expected)],np.asarray(expected,dtype=float)) or np.any(np.isfinite(actual[len(expected):])):
                raise RuntimeError(f"Origin worksheet differs from source: {spec.label}")
    project_path = output / "result.opju"
    if not op.save(str(project_path)) or project_path.stat().st_size <= 0:
        raise RuntimeError("Origin failed to save editable OPJU")
    exports = export_graph(op, graph, output / "result.png", output / "result.pdf",
                           output / "result.tif", raster_width=2600)
    # Independent saved-project reopen, retaining the same owned COM session.
    op.open(str(project_path), readonly=False, asksave=False)
    reopened = op.find_graph("Figure3")
    if reopened is None:
        raise RuntimeError("Saved OPJU did not reopen its editable graph")
    reopened_layer = reopened[0]
    if len(reopened_layer.plot_list()) != len(panel.series):
        raise RuntimeError("Saved OPJU plot count changed")
    saved_plots=[NativePlotReference(f"[{reopened.name}]1!{i+1}",reopened_layer.obj) for i in range(len(panel.series))]
    saved_labels = {name: reopened_layer.label(label_name) for name, label_name in saved_label_names.items()}
    after = readback(op, reopened, reopened_layer, panel, saved_plots, saved_labels)
    saved_sheet = op.find_sheet("w", plot_sheet_range)
    # Names of source workbooks depend on their count; locate by worksheet label.
    if saved_sheet is None:
        for book in op.pages("w"):
            for candidate in book:
                if candidate.name == "PlotData":
                    saved_sheet = candidate
                    break
    if saved_sheet is None:
        raise RuntimeError("Saved project lost the editable plot-coordinate worksheet")
    for index, expected in enumerate(numeric_before):
        actual = saved_sheet.to_list(index)
        expected_array = np.asarray(expected, dtype=float)
        actual_array = np.asarray(actual, dtype=float)
        if expected_array.shape != actual_array.shape or not np.array_equal(
                expected_array, actual_array, equal_nan=True):
            raise RuntimeError("Saved project numeric plot data changed")
    reopened.save_fig(str(output / "reopened.png"), type="png", width=2600)
    report = {"panel_id": panel.id, "style": STYLE, "before_save": before,
              "after_independent_reopen": after, "exports": exports,
              "editable_numeric_data_exact": True,
              "plot_sheet_columns": workbook_columns,
              "visual_qa_status": "pending_image_inspection"}
    write_json(output / "origin_verify_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--approval-manifest", type=Path)
    parser.add_argument("--smoke-report", type=Path)
    parser.add_argument("--keep-open", action="store_true")
    parser.add_argument("--panel", action="append", help="Render selected panel IDs from the complete approved plan")
    parser.add_argument("--output-subdir", default="origin_panels")
    args = parser.parse_args()
    panels = prepare_panels()
    plan = renderer_plan(panels)
    digest = plan_digest(plan)
    write_json(ROOT / "origin_renderer_plan.json", {**plan, "renderer_plan_sha256": digest})
    emit("renderer_plan_ready", renderer_plan_sha256=digest,
         panels=[panel.id for panel in panels], origin_called=False)
    if not args.execute:
        return 0
    if args.approval_manifest is None or args.smoke_report is None:
        parser.error("--execute requires --approval-manifest and --smoke-report")
    approval = require_approval(args.approval_manifest, digest, plan)
    smoke = require_smoke(args.smoke_report)
    from origin_sciplot.origin_backend.execution_context import require_interactive_origin_context
    from origin_sciplot.origin_backend.job_queue import origin_job_slot
    from origin_sciplot.origin_backend.session import OriginSession
    require_interactive_origin_context()
    delivery = ROOT / args.output_subdir
    if delivery.resolve().parent!=ROOT.resolve():
        raise ValueError("Output must remain directly within the task folder")
    delivery.mkdir(exist_ok=True)
    with origin_job_slot(job_kind="FlexPPI Figure3 custom renderer", wait_report_interval=30,
                         on_wait=lambda elapsed: emit("origin_job_queue", wait_seconds=elapsed)):
        with OriginSession(keep_open=args.keep_open) as session:
            op = session.op
            if op is None:
                raise RuntimeError("Origin session did not initialize")
            reports = []
            for panel in panels:
                if args.panel and panel.id not in args.panel:
                    continue
                output = delivery / panel.id
                output.mkdir(exist_ok=True)
                if any(output.iterdir()):
                    raise RuntimeError(f"Refusing to overwrite previous panel artifacts: {output}")
                emit("drawing_panel", panel_id=panel.id)
                reports.append(create_panel(op, panel, output))
            write_json(delivery / "delivery_manifest.json", {
                "renderer_plan_sha256": digest, "environment": session.environment.to_dict(),
                "approval_record": approval, "smoke_report": smoke,
                "panels": [{"id": p.id, "directory": str(delivery / p.id)} for p in panels if not args.panel or p.id in args.panel],
                "object_readback_passed": True, "visual_qa_status": "pending_image_inspection"})
    emit("origin_export_complete", panel_count=len(reports),
         output_directory=str(delivery), visual_qa_status="pending_image_inspection")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
