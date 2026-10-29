"""Prepared main OOF figure entry; no fits/forward/SSH/synthetic figures.

--preflight checks real mature baseline schema, then requires an explicit
CURRENT complete analysis audit. Pending returns exit 2 and emits QA only.
Formal render consumes independently audited caches and same-mask CSVs.
"""
# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})
# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"
# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

import argparse
from datetime import datetime, timezone
import gzip
import importlib.util
import json
from pathlib import Path
import sys
import uuid
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE
PIPE = HERE.parent
sys.path.insert(0, str(BASE))
from collect_completed_neural_exports import require_unlinked

ANALYZER = BASE / "analyze_corrected_main_comparison.py"
spec = importlib.util.spec_from_file_location("current_main_comparison", ANALYZER)
analyzer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analyzer)
METHODS = ["full", "intrachain", "static_sequence_complex", "aa_only", "backflip2_official"]
NAMES = dict(zip(METHODS, ["Full", "Intrachain", "Matched static", "AA only", "Released BF2"]))
COLORS = dict(zip(METHODS, [CATEGORICAL[0], CATEGORICAL[3], CATEGORICAL[2], CATEGORICAL[5], CATEGORICAL[4]]))
METRICS = ["rmsf", "chi_entropy"]
SOURCE_NAMES = {"atlas": "ATLAS", "dynarepo": "DynaRepo", "itrace": "ITrace"}
SOURCE_COLORS = dict(zip(SOURCE_NAMES, CATEGORICAL_EXTENDED[6:9]))
CSV_NAMES = ["per_record_metrics.csv", "all_3063_record_points.csv", "source_fold_record_and_system_macros.csv",
             "per_source_biological_system_metrics.csv", "five_fold_mean_sample_sd.csv",
             "paired_record_improvements.csv", "paired_common_support_macros.csv"]
STATUS = "complete; all actual OOF records compared on identical targets and masks"


def bind(path, evidence, expected=None):
    require_unlinked(Path(path).absolute())
    addition={}
    digest=analyzer.bind(path, addition, expected)
    for key,value in addition.items():
        analyzer.require(key not in evidence or evidence[key]==value,"Conflicting previously bound SHA: "+key)
        evidence[key]=value
    return digest


def read_frame(path):
    return pd.read_csv(path, keep_default_na=False, na_values=[""])


def compare_frame(a, b, keys, description):
    analyzer.require(not a.duplicated(keys).any() and not b.duplicated(keys).any(), "Duplicate " + description)
    x, y = a.set_index(keys).sort_index(), b.set_index(keys).sort_index()
    analyzer.require(x.index.equals(y.index) and set(x.columns) == set(y.columns), "Table membership/schema differs: " + description)
    for key in x:
        if pd.api.types.is_numeric_dtype(x[key]) or pd.api.types.is_bool_dtype(x[key]):
            analyzer.require(np.allclose(x[key].to_numpy(float), y[key].to_numpy(float), equal_nan=True, rtol=1e-9, atol=1e-9), "Recomputed table differs: " + description + "/" + key)
        else:
            analyzer.require(x[key].fillna("").equals(y[key].fillna("")), "Text values differ: " + description + "/" + key)


def baseline_schema(evidence):
    """Real mature baseline schema only, never final neural performance."""
    pre = BASE / "qa/corrected_main_comparison_preflight.json"
    bind(pre, evidence)
    report = analyzer.read_json(pre)
    analyzer.require(report["analyzer_sha256"] == bind(ANALYZER, evidence) and report["status"].startswith("prepared;") and
                     report["actual_records"] == 3063 and report["valid_residues"] == [961321, 824988] and
                     report["final_analysis_complete"] is False and report["final_neural_scores_emitted"] is False,
                     "Mature baseline preflight is not current prepared evidence")
    for name in ["training_schema.json", "training_features.npz", "static_sequence_complex_oof.npz", "aa_only_oof.npz"]:
        path = PIPE / "static_corrected" / name
        bind(path, evidence, report["input_and_baseline_sha256"][str(path.resolve())])
    schema = analyzer.read_json(PIPE / "static_corrected/training_schema.json")
    analyzer.require(len(schema["records"]) == 3063 and len({r["sample_id"] for r in schema["records"]}) == 3063, "Baseline schema membership differs")
    source_counts=pd.Series([r["dataset"] for r in schema["records"]]).value_counts().to_dict()
    fold_counts=pd.Series([r["fold"] for r in schema["records"]]).value_counts().sort_index().tolist()
    analyzer.require(source_counts=={"atlas":1929,"dynarepo":890,"itrace":244} and fold_counts==[612,612,611,614,614],
                     "Baseline source or held-out fold census differs")
    with np.load(PIPE / "static_corrected/training_features.npz", allow_pickle=False) as z:
        target, valid = z["target"], z["valid"]
        analyzer.require(target.shape == valid.shape == (961321, 2) and valid.dtype.kind == "b" and
                         valid.sum(0).tolist() == [961321, 824988], "Baseline common mask schema differs")
    for method in ["static_sequence_complex", "aa_only"]:
        with np.load(PIPE / "static_corrected" / (method + "_oof.npz"), allow_pickle=False) as z:
            analyzer.require(z["prediction"].shape == (961321, 2) and np.isfinite(z["prediction"][valid]).all(), "Invalid mature baseline prediction")
    return dict(actual_baseline_records=3063, actual_baseline_residues=961321, effective_counts=[961321, 824988],
                source_records=source_counts,heldout_records=fold_counts,
                baseline_only_schema_checked=True, final_neural_accuracy_available=False)


def completed_context(audit_path, export_root, weight_root, evidence):
    analyzer.require(audit_path is not None, "Explicit complete --analysis-audit is absent; final neural OOF unavailable")
    path = Path(audit_path).absolute()
    require_unlinked(path)
    analyzer.require(path.is_file(), "Completed analysis_audit.json absent")
    bind(path, evidence)
    audit = analyzer.read_json(path)
    analyzer.require(audit["status"] == STATUS and audit["final_analysis_complete"] is True and audit["final_neural_scores_emitted"] is True,
                     "Analysis is prepared/pending, not a completed final OOF analysis")
    analyzer.require(audit["analyzer_sha256"] == bind(ANALYZER, evidence) and audit["actual_records"] == 3063 and
                     audit["per_record_metric_rows"] == 3063*9 and audit["new_fits"] == audit["new_forward_calls"] == audit["bootstrap_replicates"] == 0,
                     "Analysis version or complete record support differs")
    for source, digest in audit["input_and_consumed_output_sha256"].items():
        bind(source, evidence, digest)
    neural_path = Path(audit["neural_audit_path"])
    bind(neural_path, evidence, audit["neural_audit_sha256"])
    neural_evidence={}
    neural, neural_hashes = analyzer.neural_gate(neural_path, export_root, weight_root, neural_evidence)
    analyzer.merge_evidence(evidence,neural_evidence)
    analyzer.require(set(audit["generated_csv_sha256"]) == set(CSV_NAMES), "Complete seven analysis CSVs required")
    frames = {}
    for name in CSV_NAMES:
        p = path.parent / name
        bind(p, evidence, audit["generated_csv_sha256"][name])
        frames[name] = read_frame(p)
    per = frames["per_record_metrics.csv"]
    analyzer.require(len(per) == 3063*9 and not per.duplicated(["sample_id", "variant", "metric"]).any(), "Incomplete record metric census")
    points = frames["all_3063_record_points.csv"]
    analyzer.require(len(points) == 3063 and points.sample_id.is_unique, "Incomplete 3063 point table")
    wanted = {(v,m) for v in METHODS for m in METRICS if not (v == "backflip2_official" and m == "chi_entropy")}
    analyzer.require(set(zip(per.variant,per.metric)) == wanted, "Missing method/target, or invented BF2 Hchi")
    for (v,m), group in per.groupby(["variant","metric"]):
        analyzer.require(len(group) == 3063 and set(group.sample_id) == set(points.sample_id), "Method lost actual records")
    macro, systems, five = analyzer.summarize_records(per)
    pairs = analyzer.paired_records(per)
    pairmacro = analyzer.summarize_pairs(pairs)
    compare_frame(macro,frames[CSV_NAMES[2]],["variant","metric","dataset","fold","estimand"],"macro")
    compare_frame(systems,frames[CSV_NAMES[3]],["variant","metric","dataset","biological_system"],"system")
    compare_frame(five,frames[CSV_NAMES[4]],["variant","metric","dataset","estimand"],"fivefold")
    compare_frame(pairs,frames[CSV_NAMES[5]],["sample_id","dataset","fold","biological_system","metric","comparator"],"pairs")
    compare_frame(pairmacro,frames[CSV_NAMES[6]],["comparator","metric","dataset","fold","score","estimand"],"common-support macro")
    schema = analyzer.read_json(PIPE / "static_corrected/training_schema.json")
    analyzer.require(set(points.sample_id)=={i["sample_id"] for i in schema["records"]}, "Point/schema identity differs")
    compare_frame(analyzer.all_record_points(per, schema["records"]), points, ["sample_id"], "all-record point table")
    analyzer.require((per["count"] > 0).all() and np.isfinite(per.mae).all(),
                     "Each real record must retain valid target residues and finite MAE")
    return dict(audit=audit,neural=neural,frames=frames,per=per,pairs=pairs,schema=schema,output_hashes=neural_hashes,
                analysis_directory=path.parent)


def sidechain_rows(context, export_root, evidence):
    """Each record/AA gets equal weight after within-record valid-residue MAE."""
    sys.path.insert(0,str(PIPE / "training_payload/src"))
    from flexppi.data.structures import parse_pdb
    from flexppi.constants import RESTYPES, THREE_TO_ONE, CHI_ATOMS
    arrays = {}
    with np.load(PIPE/"static_corrected/training_features.npz",allow_pickle=False) as z:
        target,valid=z["target"].copy(),z["valid"].copy()
    for method in ["static_sequence_complex","aa_only"]:
        with np.load(PIPE/"static_corrected"/(method+"_oof.npz"),allow_pickle=False) as z:arrays[method]=z["prediction"].copy()
    bfpath=PIPE/"backflip2_corrected_scoring/training_aligned.npz"
    bind(bfpath,evidence,context["audit"]["input_and_consumed_output_sha256"][str(bfpath.resolve())])
    with np.load(bfpath,allow_pickle=False) as z:
        analyzer.exact_array(z["target"],target[:,0],"BF/common target");analyzer.exact_array(z["valid"],valid[:,0],"BF/common mask")
        arrays["backflip2_official"]=z["prediction"].copy()
    qmap={one:len(CHI_ATOMS.get(three,[])) for three,one in THREE_TO_ONE.items() if three in CHI_ATOMS or three in ["ALA","GLY"]}
    out=[]
    records=analyzer.read_lines(PIPE/"training_payload/splits/all_ready.jsonl");lookup={r["sample_id"]:r for r in records}
    for info in context["schema"]["records"]:
        span=slice(info["start"],info["stop"]);r=lookup[info["sample_id"]]
        ip=analyzer.local_input(r,"structure_path");bind(ip,evidence,info["input_sha256"]);s=parse_pdb(ip)
        aa=np.asarray(list(s.sequence));y,v=target[span],valid[span]
        prediction={m:values[span] for m,values in arrays.items()}
        for method in ["full","intrachain"]:
            p=export_root/"evaluation"/method/(info["sample_id"]+".npz");bind(p,evidence,context["output_hashes"][str(p.resolve())])
            with np.load(p,allow_pickle=False) as z:
                analyzer.exact_array(z["target"][:,:2],y,"AA/neural target");analyzer.exact_array(z["valid"][:,:2],v,"AA/neural mask")
                analyzer.exact_array(z["amino_acid"],aa,"AA/neural identity")
                prediction[method]=z["mean"].copy()
        for letter in RESTYPES:
            for method in METHODS:
                for m,metric in enumerate(METRICS):
                    if method=="backflip2_official" and m==1:continue
                    mask=(aa==letter)&v[:,m];p=prediction[method];p=p[:,m] if p.ndim==2 else p
                    scored=analyzer.statistics(y[mask,m],p[mask])
                    out.append(dict(sample_id=info["sample_id"],dataset=info["dataset"],fold=info["fold"],
                        biological_system=info["biological_system"],amino_acid=letter,canonical_chi_count=qmap[letter],
                        variant=method,metric=metric,**scored))
    frame=pd.DataFrame(out)
    analyzer.require(len(frame)==3063*20*9 and not frame.duplicated(["sample_id","amino_acid","variant","metric"]).any(),"AA record census differs")
    for (method,metric),part in frame.groupby(["variant","metric"]):
        analyzer.require(int(part["count"].sum())==[961321,824988][METRICS.index(metric)],
                         "AA groups omitted valid residues: "+method+"/"+metric)
    macros=[]
    for (method,metric,aa),part in frame.groupby(["variant","metric","amino_acid"]):
        for source,sub in [("all",part)]+list(part.groupby("dataset")):
            rows=sub[sub["count"]>0];systems=rows.groupby(["dataset","biological_system"]).mae.mean()
            macros.append(dict(variant=method,metric=metric,amino_acid=aa,canonical_chi_count=int(part.canonical_chi_count.iloc[0]),
                dataset=source,records=len(sub),records_valid=len(rows),valid_residues=int(sub["count"].sum()),
                record_macro_mae=analyzer.finite_mean(sub.mae),source_system_macro_mae=analyzer.finite_mean(systems),
                spearman_records_defined=int(np.isfinite(sub.spearman).sum()),pearson_records_defined=int(np.isfinite(sub.pearson).sum())))
    keys = ["sample_id", "dataset", "fold", "biological_system", "amino_acid", "canonical_chi_count", "metric"]
    full = frame[frame.variant == "full"][keys + ["count", "mae"]]
    pairs = []
    for method in METHODS[1:]:
        part = full.merge(frame[frame.variant == method][keys + ["count", "mae"]], on=keys,
                          suffixes=("_full", "_comparator"), validate="one_to_one")
        analyzer.require(len(part) == 3063 * 20 * (1 if method == "backflip2_official" else 2), "AA paired census differs")
        analyzer.require((part.count_full == part.count_comparator).all(), "AA paired valid support differs")
        part["comparator"] = method
        part["mae_full_defined"] = np.isfinite(part.mae_full)
        part["mae_comparator_defined"] = np.isfinite(part.mae_comparator)
        part["mae_both_defined"] = part.mae_full_defined & part.mae_comparator_defined
        part["mae_improvement"] = (part.mae_comparator - part.mae_full).where(part.mae_both_defined, np.nan)
        pairs.append(part)
    paired = pd.concat(pairs, ignore_index=True)
    paired_macros = []
    for (method, metric, aa), part in paired.groupby(["comparator", "metric", "amino_acid"]):
        for source, sub in [("all", part)] + list(part.groupby("dataset")):
            common = sub[sub.mae_both_defined]
            system = common.groupby(["dataset", "biological_system"])[["mae_full", "mae_comparator", "mae_improvement"]].mean()
            paired_macros.append(dict(comparator=method, metric=metric, amino_acid=aa, dataset=source,
                canonical_chi_count=int(part.canonical_chi_count.iloc[0]), records=len(sub),
                full_records_defined=int(sub.mae_full_defined.sum()), comparator_records_defined=int(sub.mae_comparator_defined.sum()),
                paired_records_defined=len(common), paired_systems_defined=len(system),
                valid_residues=int(sub.count_full.sum()), record_macro_improvement=analyzer.finite_mean(common.mae_improvement),
                source_system_macro_improvement=analyzer.finite_mean(system.mae_improvement),
                full_record_mean_common_support=analyzer.finite_mean(common.mae_full),
                comparator_record_mean_common_support=analyzer.finite_mean(common.mae_comparator),
                negative_records=int((common.mae_improvement < 0).sum()), zero_records=int((common.mae_improvement == 0).sum()),
                positive_records=int((common.mae_improvement > 0).sum())))
    analyzer.require(len(paired) == 3063 * 20 * 7, "Complete AA pair census differs")
    return frame, pd.DataFrame(macros), paired, pd.DataFrame(paired_macros)


def correlation_shared_support(per):
    """Joint finite support across 5/4 methods, separate for each correlation."""
    census, summaries=[] ,[]
    identity=["sample_id","dataset","fold","biological_system"]
    for metric in METRICS:
        methods=[m for m in METHODS if not (metric=="chi_entropy" and m=="backflip2_official")]
        part=per[per.metric==metric]
        infos=part[identity].drop_duplicates().set_index("sample_id")
        analyzer.require(len(infos)==3063 and infos.index.is_unique,"Shared correlation identity census differs")
        for score in ["spearman","pearson"]:
            wide=part.pivot(index="sample_id",columns="variant",values=score).reindex(columns=methods)
            analyzer.require(len(wide)==3063 and set(wide.index)==set(infos.index),"Shared correlation membership differs")
            finite=np.isfinite(wide.to_numpy(float))
            common=finite.all(axis=1)
            common_ids=set(wide.index[common])
            for i,(sample_id,values) in enumerate(wide.iterrows()):
                row=dict(sample_id=sample_id,metric=metric,score=score,**infos.loc[sample_id].to_dict(),
                         required_methods="|".join(methods),required_method_count=len(methods),
                         own_methods_defined=int(finite[i].sum()),common_defined=bool(common[i]))
                for j,method in enumerate(methods):
                    row[method+"__value"]=float(values[method]);row[method+"__defined"]=bool(finite[i,j])
                census.append(row)
            for method in methods:
                original=part[part.variant==method]
                for source in ["all",*SOURCE_NAMES]:
                    src=original if source=="all" else original[original.dataset==source]
                    for fold in [-1,0,1,2,3,4]:
                        own=src if fold==-1 else src[src.fold==fold]
                        joint=own[own.sample_id.isin(common_ids)]
                        system=joint.groupby(["dataset","biological_system"])[score].mean()
                        for estimand,units in [("record_macro",joint[score]),("source_system_macro",system)]:
                            summaries.append(dict(variant=method,metric=metric,score=score,dataset=source,fold=fold,estimand=estimand,
                                required_methods="|".join(methods),records_original=len(own),
                                records_own_defined=int(np.isfinite(own[score]).sum()),records_common_defined=len(joint),
                                records_common_undefined=len(own)-len(joint),
                                biological_systems_original=int(own.groupby(["dataset","biological_system"]).ngroups),
                                biological_systems_common=len(system),units_common_defined=int(np.isfinite(units).sum()),
                                mean_common_support=analyzer.finite_mean(units)))
    return pd.DataFrame(census),pd.DataFrame(summaries)


def render(context, side, side_pairs, output):
    """No call until complete audits and source hashes have passed."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    mpl.rcParams.update({"xtick.labelsize":7.5,"ytick.labelsize":7.5,"legend.fontsize":7.5})
    per,pairs=context["per"],context["pairs"]
    macro=context["frames"]["source_fold_record_and_system_macros.csv"]
    figure_specs=[]
    def export_fixed(fig, stem):
        # Explicit full-canvas bounds preserve 183 mm and actual point sizes;
        # a tight crop would change the width before manuscript placement.
        texts=[t for t in fig.findobj(mpl.text.Text) if t.get_text()]
        floor=min(t.get_fontsize() for t in texts)
        analyzer.require(floor >= 7.5, "Rendered text below 7.5 pt")
        width,height=fig.get_size_inches()*25.4
        analyzer.require(np.isclose(width,183), "Canvas differs from 183 mm")
        for suffix in ["pdf","png","svg"]:
            fig.savefig(output/(stem+"."+suffix),bbox_inches=fig.bbox_inches,dpi=300)
        figure_specs.append(dict(stem=stem,width_mm=float(width),height_mm=float(height),minimum_text_pt=float(floor),
                                 text_objects=len(texts),png_dpi=300,pdf_fonttype=42,svg_editable_text=True))
    fig_accuracy,axes_accuracy=plt.subplots(2,2,figsize=(183/25.4,145/25.4))
    fig_partner,axes_partner=plt.subplots(2,2,figsize=(183/25.4,145/25.4))
    for f in [fig_accuracy,fig_partner]:
        f.subplots_adjust(left=.13,right=.98,bottom=.15,top=.83,hspace=.85,wspace=.52)
    axes=np.vstack([axes_accuracy,axes_partner])
    fig_accuracy.suptitle("Corrected OOF: common labels and held-out members",y=.99)
    fig_accuracy.legend(handles=[Line2D([],[],color=COLORS[m],label=NAMES[m],linewidth=1.3) for m in METHODS],
               loc="upper center",bbox_to_anchor=(.52,.969),ncol=3)
    fig_partner.suptitle("Corrected OOF: partner ablation and source sensitivity",y=.99)
    fig_partner.legend(handles=[Line2D([],[],color=SOURCE_COLORS[s],marker="o",label=SOURCE_NAMES[s],linewidth=0,markersize=4) for s in SOURCE_NAMES],
                       loc="upper center",bbox_to_anchor=(.52,.969),ncol=3)
    common_census,common_macros=correlation_shared_support(per)
    common_census.to_csv(output/"all_3063_correlation_joint_support.csv",index=False,na_rep="")
    common_macros.to_csv(output/"correlation_joint_support_record_and_system_macros.csv",index=False,na_rep="")
    traces,points,correlation_points,support=[] ,[],[],[]
    for col,metric in enumerate(METRICS):
        unit="Å" if metric=="rmsf" else "nats"
        for method in METHODS:
            if method=="backflip2_official" and col==1:continue
            group=per[(per.variant==method)&(per.metric==metric)]
            good=group[np.isfinite(group.mae)].sort_values("mae",kind="mergesort")
            analyzer.require(len(good)==3063,"ECDF must retain all real record MAE values")
            x=good.mae.to_numpy();ecdf=np.arange(1,len(x)+1)/len(x)
            axes[0,col].step(x,ecdf,where="post",color=COLORS[method],linewidth=1)
            traces.extend(dict(sample_id=r.sample_id,variant=method,metric=metric,mae=float(r.mae),ecdf=float(y)) for r,y in zip(good.itertuples(),ecdf))
            support.append(dict(panel="accuracy_ecdf",variant=method,metric=metric,records=len(group),
                mae_records_defined=len(good),mae_records_undefined=len(group)-len(good),
                spearman_records_defined=int(np.isfinite(group.spearman).sum()),
                spearman_records_undefined=int((~np.isfinite(group.spearman)).sum()),
                pearson_records_defined=int(np.isfinite(group.pearson).sum()),
                pearson_records_undefined=int((~np.isfinite(group.pearson)).sum())))
        axes[0,col].set_xlabel("Within-record MAE ("+unit+")");axes[0,col].set_ylabel("Fraction of finite records")
        axes[0,col].set_title(("RMSF" if col==0 else "Hχ")+": all-record accuracy (n = 3063)",loc="left")
        if col==0:
            from matplotlib.ticker import FuncFormatter, NullFormatter
            axes[0,col].set_xscale("symlog",linthresh=.1)
            # Decimal labels preserve axis scale and avoid 5.25 pt exponents.
            axes[0,col].xaxis.set_major_formatter(FuncFormatter(lambda value, position: format(value, "g")))
            axes[0,col].xaxis.set_minor_formatter(NullFormatter())
        for i,method in enumerate(METHODS):
            if method=="backflip2_official" and col==1:continue
            group=common_macros[(common_macros.variant==method)&(common_macros.metric==metric)&
                                (common_macros.dataset=="all")&(common_macros.score=="spearman")]
            folds=group[(group.fold>=0)&(group.estimand=="record_macro")].sort_values("fold")
            analyzer.require(len(folds)==5,"Correlation diagnostic requires all five held-out folds")
            for xx,r in zip(i+np.linspace(-.16,.16,5),folds.itertuples()):
                axes[1,col].scatter(xx,r.mean_common_support,s=8,color=COLORS[method],alpha=.55,edgecolors="none")
                correlation_points.append(dict(variant=method,metric=metric,fold=int(r.fold),estimand="record_macro",
                    x=float(xx),spearman=float(r.mean_common_support),records=int(r.records_original),
                    records_own_defined=int(r.records_own_defined),records_common_defined=int(r.records_common_defined),
                    units_common_defined=int(r.units_common_defined)))
            for estimand,marker,shift in [("record_macro","o",-.055),("source_system_macro","D",.055)]:
                row=group[(group.fold==-1)&(group.estimand==estimand)]
                analyzer.require(len(row)==1,"Correlation macro support missing")
                r=row.iloc[0];axes[1,col].scatter(i+shift,r.mean_common_support,s=22,color=COLORS[method],marker=marker,edgecolors=BLACK,linewidths=.4)
                correlation_points.append(dict(variant=method,metric=metric,fold=-1,estimand=estimand,
                    x=float(i+shift),spearman=float(r.mean_common_support),records=int(r.records_original),
                    records_own_defined=int(r.records_own_defined),records_common_defined=int(r.records_common_defined),
                    units_common_defined=int(r.units_common_defined)))
        axes[1,col].set_xticks(range(5 if col==0 else 4),["Full","Intra","Static","AA","BF2"][:5 if col==0 else 4])
        joint_n=int(common_macros[(common_macros.variant=="full")&(common_macros.metric==metric)&
                    (common_macros.dataset=="all")&(common_macros.score=="spearman")&(common_macros.fold==-1)&
                    (common_macros.estimand=="record_macro")].records_common_defined.iloc[0])
        axes[1,col].set_ylabel("Spearman ρ, same finite records")
        axes[1,col].set_title(("RMSF" if col==0 else "Hχ")+": joint-finite correlation (n = "+str(joint_n)+")",loc="left")
        part=pairs[(pairs.comparator=="intrachain")&(pairs.metric==metric)].sort_values(["dataset","biological_system","sample_id"],kind="mergesort")
        analyzer.require(len(part)==3063 and np.isfinite(part.mae_improvement).all(),"Partner panel must retain all 3063 signed differences")
        offset=0
        for source in SOURCE_NAMES:
            group=part[part.dataset==source];x=np.arange(offset,offset+len(group));y=group.mae_improvement.to_numpy()
            axes[2,col].scatter(x,y,s=3,color=SOURCE_COLORS[source],alpha=.5,linewidths=0,rasterized=True)
            points.extend(dict(sample_id=r.sample_id,dataset=source,metric=metric,x=int(xx),mae_improvement=float(r.mae_improvement)) for r,xx in zip(group.itertuples(),x))
            support.append(dict(panel="partner_difference",variant="full_vs_intrachain",metric=metric,dataset=source,
                records=len(group),paired_records_defined=int(np.isfinite(y).sum()),paired_records_undefined=int((~np.isfinite(y)).sum()),
                negative_records=int((y<0).sum()),zero_records=int((y==0).sum()),positive_records=int((y>0).sum())))
            offset+=len(group)
        axes[2,col].axhline(0,color=GREY,linewidth=.7,linestyle="--")
        axes[2,col].set_xlabel("All 3063 records, source / system order");axes[2,col].set_ylabel("Intrachain − full MAE ("+unit+")")
        axes[2,col].set_title(("RMSF" if col==0 else "Hχ")+": partner-input ablation",loc="left")
        for i,method in enumerate(METHODS):
            if method=="backflip2_official" and col==1:continue
            for j,source in enumerate(SOURCE_NAMES):
                group=macro[(macro.variant==method)&(macro.metric==metric)&(macro.dataset==source)&(macro.fold==-1)]
                for estimand,marker,shift in [("record_macro","o",-.07),("source_system_macro","D",.07)]:
                    row=group[group.estimand==estimand];analyzer.require(len(row)==1,"Source macro row missing")
                    axes[3,col].scatter(row.mae.to_numpy(),[i+(j-1)*.22+shift],s=15,color=SOURCE_COLORS[source],marker=marker,edgecolors="none")
        axes[3,col].set_yticks(range(5 if col==0 else 4),[NAMES[m] for m in METHODS[:5 if col==0 else 4]])
        axes[3,col].invert_yaxis();axes[3,col].set_xlabel("Source macro MAE ("+unit+")")
        axes[3,col].set_title(("RMSF" if col==0 else "Hχ")+": source / averaging sensitivity",loc="left")
    for f,axs in [(fig_accuracy,axes_accuracy),(fig_partner,axes_partner)]:
        for k,ax in enumerate(axs.flat):
            b=ax.get_position();f.text(.01 if k%2==0 else .505,b.y1+.008,chr(97+k),fontsize=10,fontweight="bold")
    fig_accuracy.text(.13,.025,"MAE retains every record. Correlations use the joint finite intersection of all five / four methods.\nCircles: record macro; diamonds: first within source-specific system, then macro; small dots: five folds.\nSupport census preserves own-defined and excluded records; no CI or p values.",fontsize=7.5)
    fig_partner.text(.13,.025,"Each record uses its held-out fold; all 3063 negative / zero / positive ablation differences remain.\nPositive difference favors full. Circles: record macro; diamonds: first within source-specific system.\nAblation is descriptive and does not establish causal partner effects. No CI or p values.",fontsize=7.5)
    export_fixed(fig_accuracy,"fig_main_corrected_oof_accuracy");plt.close(fig_accuracy)
    export_fixed(fig_partner,"fig_main_corrected_oof_partner_source");plt.close(fig_partner)
    pd.DataFrame(traces).to_csv(output/"plot_all_record_mae_ecdf.csv",index=False)
    pd.DataFrame(points).to_csv(output/"plot_all_record_partner_mae.csv",index=False)
    pd.DataFrame(correlation_points).to_csv(output/"plot_correlation_macro_and_fold_points.csv",index=False,na_rep="")
    pd.DataFrame(support).to_csv(output/"plot_support_counts.csv",index=False,na_rep="")
    aa_macro=side
    fig,axes=plt.subplots(2,2,figsize=(183/25.4,165/25.4))
    fig.subplots_adjust(left=.13,right=.98,bottom=.13,top=.83,hspace=.72,wspace=.52)
    fig.suptitle("Side-chain diagnostics: every amino-acid / record group",y=.99)
    fig.legend(handles=[Line2D([],[],color=COLORS[m],marker="o",label=NAMES[m],linewidth=.8,markersize=3) for m in METHODS],loc="upper center",bbox_to_anchor=(.5,.966),ncol=3)
    aa_order=list("AGCSTVNDHILFPWYQEMRK")
    for col,metric in enumerate(METRICS):
        for method in METHODS:
            if method=="backflip2_official" and col==1:continue
            rows=aa_macro[(aa_macro.variant==method)&(aa_macro.metric==metric)&(aa_macro.dataset=="all")].set_index("amino_acid").reindex(aa_order)
            axes[0,col].plot(range(20),rows.record_macro_mae,color=COLORS[method],marker="o",markersize=2.5,linewidth=.9)
        for method,style in [("intrachain","-"),("static_sequence_complex","--")]:
            rows=side_pairs[(side_pairs.comparator==method)&(side_pairs.metric==metric)&(side_pairs.dataset=="all")].set_index("amino_acid").reindex(aa_order)
            axes[1,col].plot(range(20),rows.record_macro_improvement,color=COLORS[method],linestyle=style,marker="o",markersize=2.5,linewidth=.9)
        axes[1,col].axhline(0,color=GREY,linewidth=.7,linestyle=":")
        axes[0,col].set_ylabel("Mean within-record / AA MAE ("+("Å" if col==0 else "nats")+")")
        axes[1,col].set_ylabel("Comparator − full MAE ("+("Å" if col==0 else "nats")+")")
        axes[0,col].set_title(("RMSF" if col==0 else "Hχ")+": AA-specific absolute MAE",loc="left")
        axes[1,col].set_title(("RMSF" if col==0 else "Hχ")+": paired AA-specific improvement",loc="left")
        for row in range(2):
            axes[row,col].set_xticks(range(20),aa_order)
            axes[row,col].set_xlabel("AA in canonical χ-count order")
    for k,ax in enumerate(axes.flat):
        b=ax.get_position();fig.text(.01 if k%2==0 else .505,b.y1+.01,chr(97+k),fontsize=10,fontweight="bold")
    fig.text(.13,.022,"Within-record / AA means precede equal-record macros; improvement is paired on identical valid residues.\nA/G Hχ remains NaN. Intrachain: solid orange; static: dashed green in c–d. Full support census is in CSV.\nχ count affects target scale; curves do not establish intrinsic chemical difficulty. No CI or p values.",fontsize=7.5)
    export_fixed(fig,"fig_diagnostic_corrected_oof_amino_acids");plt.close(fig)
    return figure_specs


def refuse_report(report, error, status, stage):
    """Reject completion without discarding actual output counts or paths."""
    report.update(status=status,error_type=type(error).__name__,error=str(error),
        refusal_stage=stage,final_neural_results_complete=False,final_figures_complete=False,
        complete_input_gate=False,output_directory_complete=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight",action="store_true")
    parser.add_argument("--analysis-audit",type=Path)
    parser.add_argument("--export-root",type=Path,default=PIPE/"neural_exports")
    parser.add_argument("--weight-root",type=Path,default=PIPE/"trained_checkpoints")
    args=parser.parse_args();evidence={}
    report=dict(status="prepared; no final figures",checked_utc=datetime.now(timezone.utc).isoformat(),
        plot_script_sha256=analyzer.sha(__file__),final_figures_complete=False,final_neural_results_complete=False,
        complete_input_gate=False,output_directory_complete=False,
        figures_generated=0,new_fits=0,new_forward_calls=0,remote_operations=0,
        runtime_versions=dict(python=sys.version.split()[0],numpy=np.__version__,pandas=pd.__version__,matplotlib=mpl.__version__))
    code=0
    try:
        bind(Path(__file__),evidence,report["plot_script_sha256"])
        bind(HERE/"prepare_and_plot_main_oof.py",evidence,"0c2e488d07dbf5c841d602493aabd99d4f51e87fe078e1d1deafdcaa429e30a5")
        report["presentation_change"]="Decimal symlog labels only; no data, masks, axes, method, support or statistics changes"
        bind(BASE/"scripts/collect_completed_neural_exports.py",evidence)
        require_unlinked(args.export_root.absolute());require_unlinked(args.weight_root.absolute())
        report.update(baseline_schema(evidence))
        context=completed_context(args.analysis_audit,args.export_root.resolve(),args.weight_root.resolve(),evidence)
        analyzer.check_snapshot(evidence)
        if args.preflight:
            report.update(status="ready entry; preflight only; no figures generated",complete_input_gate=True)
        else:
            side,side_macro,side_pairs,side_pair_macro=sidechain_rows(context,args.export_root.resolve(),evidence)
            analyzer.check_snapshot(evidence)
            output=HERE/"completed"/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8])
            require_unlinked(output);output.mkdir(parents=True,exist_ok=False)
            report.update(output_directory=str(output),output_directory_complete=False)
            side.to_csv(output/"all_record_amino_acid_metrics.csv.gz",index=False,na_rep="",compression="gzip")
            side_macro.to_csv(output/"amino_acid_support_and_macro_metrics.csv",index=False,na_rep="")
            side_pairs.to_csv(output/"all_record_amino_acid_paired_mae.csv.gz",index=False,na_rep="",compression="gzip")
            side_pair_macro.to_csv(output/"amino_acid_paired_mae_support_and_macros.csv",index=False,na_rep="")
            for name in CSV_NAMES:(output/name).write_bytes((context["analysis_directory"]/name).read_bytes())
            figure_specs=render(context,side_macro,side_pair_macro,output)
            analyzer.check_snapshot(evidence)
            report.update(status="rendered complete data; actual visual review still required",complete_input_gate=True,
                figures_generated=3,final_neural_results_complete=True,output_directory=str(output),
                output_sha256={p.name:analyzer.sha(p) for p in output.iterdir() if p.is_file()},
                all_records=3063,per_record_metric_rows=27567,all_record_amino_acid_rows=len(side),
                all_record_amino_acid_paired_mae_rows=len(side_pairs),valid_counts=[961321,824988],
                font_floor_pt=7.5,figure_width_mm=183,figure_specs=figure_specs,visual_review_complete=False,
                input_sha256=evidence,output_directory_complete=True)
            # Completion refers to actual outputs only after a separate visual
            # review; raw render never labels final_figures_complete=True.
            (output/"render_audit.json").write_text(json.dumps(report,indent=2,allow_nan=False),encoding="utf-8")
    except Exception as error:
        refuse_report(report,error,"pending/refused; complete predictions unavailable or contract mismatch","execution")
        code=2
    finally:
        report["input_sha256"]=evidence
        report["completed_utc"]=datetime.now(timezone.utc).isoformat()
        # Even a refused entry must report the actual code hash and detect a
        # file changed while it was checking the mature baseline.
        try:
            analyzer.check_snapshot(evidence)
        except Exception as error:
            refuse_report(report,error,"refused; bound source changed during execution","final_snapshot")
            code=2
        try:
            (HERE/("entry_contract_audit_typography_v2_"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8]+".json")).write_text(json.dumps(report,indent=2,allow_nan=False),encoding="utf-8")
        except Exception as error:
            refuse_report(report,error,"refused; final entry audit could not be serialized or written","entry_audit_write")
            code=2
        print(json.dumps({k:report[k] for k in ["status","figures_generated","final_figures_complete","final_neural_results_complete"]}))
    return code


if __name__=="__main__":
    raise SystemExit(main())
