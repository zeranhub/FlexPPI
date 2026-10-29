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


"""Plot every full-cohort completed method, with fixed rows across strata.

The journal reference is AlphaMissense Fig2A-C; this adapts its aligned method
rows and independent benchmark columns to points/intervals, not truncated bars.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Asset confirmation: visual adaptation of the existing benchmark dot/interval
# script and its inspected PNG. The first63 lines are copied byte-for-byte.
# Archetype: quantitative_grid. This module never runs inference or statistics.
import argparse
import struct
from matplotlib.transforms import Bbox

ROOT = Path(__file__).resolve().parent
BENCHMARK_ROOT = ROOT.parent
OUT = ROOT/'figures'
GROUPS = ['protein','protein_complex']
TITLES = {'protein':'Proteins','protein_complex':'Protein complexes'}
ORDER = ['Full','Intra','backflip-2.1 (assembly)','backflip-2.1 (chains)',
    'backflip-1.0 (assembly)','backflip-1.0 (chains)',
    'backflip-1.0-seq (assembly)','backflip-1.0-seq (chains)',
    'DynaProt-M','FastProtFlex','ANM','GNM','Static','AA']
DISPLAY = {'Full':'FlexPPI','Intra':'FlexPPI · intra',
    'backflip-2.1 (assembly)':'BackFlip-2.1 · assembly',
    'backflip-2.1 (chains)':'BackFlip-2.1 · chains',
    'backflip-1.0 (assembly)':'BackFlip-1 · assembly',
    'backflip-1.0 (chains)':'BackFlip-1 · chains',
    'backflip-1.0-seq (assembly)':'BackFlip-1 + seq · assembly',
    'backflip-1.0-seq (chains)':'BackFlip-1 + seq · chains',
    'Static':'Static features','AA':'Residue-type mean'}
BLUE = '#31688E'; OCHRE = '#B78038'; TEAL = '#348980'; INK = '#253846'
RELATIVE = {'FastProtFlex','ANM','GNM'}
METRICS = ['spearman','top10_ap','mae']
ECDF_METHODS = ['Full','backflip-2.1 (assembly)','backflip-2.1 (chains)',
    'DynaProt-M','FastProtFlex','ANM']
# Named styles keep every existing method stable when other entries are added.
ECDF_LINE_STYLES = {'Full':'-','backflip-2.1 (assembly)':'-',
    'backflip-2.1 (chains)':'--','DynaProt-M':':','FastProtFlex':'-.',
    'ANM':(0,(5,2,1,2,1,2))}
MAIN_SIZE_MM = [183,141]
ECDF_SIZE_MM = [183,78]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024**2), b''):
            result.update(block)
    return result.hexdigest()


def visual(name, meta):
    color = BLUE if name in ('Full','Intra') else TEAL if name in ('ANM','GNM','Static','AA') else OCHRE
    context = meta[name]['input_context']
    if name=='Full':
        return color,'s',3.5
    if context in ('independent_chain','within_chain','within_chain_encoder',
                   'independent_chain_predictions_reassembled_per_record'):
        return color,'D',2.9
    if 'sequence' in context:
        return color,'^',3.1
    return color,'o',3.0


def membership_file(value):
    path = Path(value)
    path = path if path.is_absolute() else BENCHMARK_ROOT/path
    require(path.is_file() and path.resolve().is_relative_to(BENCHMARK_ROOT.resolve()),
            f'Frozen membership missing/outside benchmark root: {path}')
    return path.resolve()


def load_inputs():
    paths = {name:ROOT/'evaluation'/name for name in
        ('group_summary.csv','per_record_metrics.csv','coverage_summary.csv','scoring_audit.json')}
    exposure = ROOT/'references/training_exposure_protocol_sources_sha256.json'
    require(all(p.is_file() for p in paths.values()) and exposure.is_file(),
            'Real completed two-group aggregation and frozen training-exposure audit are required before rendering')
    audit_bytes = paths['scoring_audit.json'].read_bytes()
    audit = json.loads(audit_bytes.decode('utf-8-sig'))
    require(len(audit['complete_methods'])==14 and set(audit['complete_methods'])==set(ORDER),
            'Exactly the14 completed method/configuration entries are required; partial/heavy methods cannot enter this plot')
    meta = audit['methods']
    require(set(ORDER)<=set(meta), 'Method metadata is incomplete')
    require({name for name in ORDER if not meta[name]['absolute_angstrom']}==RELATIVE,
            'Relative-output method metadata changed; never assign Angstrom MAE automatically')
    counts = {key:int(value) for key,value in audit['group_counts'].items()}
    residues = {key:int(value) for key,value in audit['group_residues'].items()}
    require(set(counts)==set(residues)==set(GROUPS) and all(n>0 for n in counts.values()) and
            sum(counts.values())==3063 and sum(residues.values())==961321,
            'The real exhaustive disjoint protein/complex cohort must total3063 records and961321 positions')
    membership = membership_file(audit['membership_path'])
    require(digest(membership)==audit['membership_sha256'], 'Membership SHA differs from the frozen scoring input')
    hashes = {str(p.resolve()):digest(p) for p in [*paths.values(),exposure,membership,Path(__file__)]}
    require(hashes[str(paths['scoring_audit.json'].resolve())]==hashlib.sha256(audit_bytes).hexdigest() and
            hashes[str(membership)]==audit['membership_sha256'],
            'Scoring audit or frozen membership changed during its first verified read')
    members = pd.read_csv(membership,keep_default_na=False,dtype={'sample_id':str})
    member_group_column = 'group' if 'group' in members.columns else 'benchmark_group'
    if member_group_column=='benchmark_group' and 'benchmark_group' in members.columns:
        members = members.rename(columns={'benchmark_group':'group'})
    elif {'group','benchmark_group'}<=set(members.columns):
        require((members.group==members.benchmark_group).all(), 'Membership group aliases disagree')
    require({'sample_id','group'}<=set(members.columns) and len(members)==3063 and members.sample_id.is_unique,
            'Frozen membership must have3063 unique sample_id rows and an explicit group column')
    require(set(members.group)==set(GROUPS) and members.groupby('group').size().to_dict()==counts,
            'Frozen membership group counts differ from the scoring audit')
    identity = members.set_index('sample_id')['group'].to_dict()
    summary = pd.read_csv(paths['group_summary.csv'])
    per_all = pd.read_csv(paths['per_record_metrics.csv'],dtype={'sample_id':str})
    coverage = pd.read_csv(paths['coverage_summary.csv'])
    require({'method','group','metric','mean','ci_low','ci_high','n_records','n_families'}<=set(summary.columns),
            'Summary schema differs from the declared scorer output')
    require({'method','sample_id','group','family','n','spearman'}<=set(per_all.columns), 'Per-record schema differs')
    excluded_methods = sorted(set(per_all.method)-set(ORDER))
    per = per_all.loc[per_all.method.isin(ORDER)].copy()
    require(set(per.method)==set(ORDER), 'Per-record results omit a completed method')
    families = {}
    for name in ORDER:
        block = per.loc[per.method==name]
        require(len(block)==3063 and block.sample_id.is_unique and set(block.sample_id)==set(identity),
                f'Incomplete/duplicate record identities: {name}')
        require(all(identity[row.sample_id]==row.group for row in block.itertuples()),
                f'{name}: per-record group labels differ from frozen membership')
        require(block.groupby('group').size().to_dict()==counts and
                np.isfinite(block['n']).all() and (block['n']>0).all() and
                (block['n']==np.floor(block['n'])).all() and
                block.groupby('group')['n'].sum().astype(int).to_dict()==residues,
                f'{name}: record/residue coverage differs between groups')
        current = block.groupby('group')['family'].nunique().to_dict()
        require(not families or current==families, 'Family support differs between compared methods')
        families = current
    coverage = coverage.loc[coverage.method.isin(ORDER)&coverage.group.isin(GROUPS)]
    require(len(coverage)==28 and not coverage.duplicated(['method','group']).any(),
            'Expected exactly14 methods by2 coverage rows')
    for row in coverage.itertuples():
        require(int(row.records)==int(row.complete_records)==counts[row.group] and
                int(row.expected_residues)==int(row.covered_residues)==residues[row.group],
                f'Incomplete coverage: {row.method}/{row.group}')
    point_rows,relative_cells = [],[]
    for group in GROUPS:
        for metric in METRICS:
            for name in ORDER:
                selected = summary.loc[(summary.method==name)&(summary.group==group)&(summary.metric==metric)]
                if metric=='mae' and name in RELATIVE:
                    require(selected.empty, f'Unexpected physical MAE for relative output: {name}')
                    relative_cells.append(dict(method=name,group=group,metric=metric,label='relative score'))
                    continue
                require(len(selected)==1, f'Missing/duplicate summary cell: {name}/{group}/{metric}')
                row = selected.iloc[0]
                require(int(row.n_records)==counts[group] and int(row.n_families)==families[group] and
                        np.isfinite(row['mean']), f'Invalid summary support/value: {name}/{group}/{metric}')
                finite = np.isfinite([row.ci_low,row.ci_high])
                require(finite.all() or (not finite.any() and families[group]<2),
                        f'Incomplete family confidence interval: {name}/{group}/{metric}')
                if finite.all():
                    require(row.ci_low<=row.ci_high, 'Confidence interval endpoints reversed')
                value_range = [-1,1] if metric=='spearman' else [0,1] if metric=='top10_ap' else [0,np.inf]
                actual = np.asarray([row['mean'],row.ci_low,row.ci_high],dtype=float)
                actual = actual[np.isfinite(actual)]
                require((actual>=value_range[0]-1e-12).all() and (actual<=value_range[1]+1e-12).all(),
                        f'Metric outside its physical/statistical range: {name}/{group}/{metric}')
                point_rows.append(row.to_dict())
    require(len(point_rows)==78 and len(relative_cells)==6, 'Expected78 measured points and6 relative-output annotations')
    for name in ECDF_METHODS:
        values = per.loc[per.method==name,'spearman'].to_numpy()
        require(len(values)==3063 and np.isfinite(values).all() and (values>=-1).all() and (values<=1).all(),
                f'ECDF cannot drop observations or hide invalid correlations: {name}')
    require(all(digest(Path(p))==h for p,h in hashes.items()), 'Input changed while its support was being checked')
    return dict(summary=summary,per=per,audit=audit,meta=meta,counts=counts,residues=residues,
        families=families,membership=membership,hashes=hashes,point_rows=point_rows,
        relative_cells=relative_cells,exposure=exposure,excluded_methods=excluded_methods,
        member_group_column=member_group_column)


def save_exact_canvas(fig, stem):
    # Preserve the mandated original export helper above, but use the explicit
    # canvas bounds here so tight-cropping cannot change the183mm journal width.
    box = Bbox.from_bounds(0,0,*fig.get_size_inches())
    for suffix in ('pdf','png','svg'):
        fig.savefig(stem.with_suffix('.'+suffix),bbox_inches=box,dpi=300)


def draw_main(data):
    points = pd.DataFrame(data['point_rows'])
    rho = points.loc[points.metric=='spearman',['mean','ci_low','ci_high']].to_numpy()
    rank_left = max(-1.,min(-.1,float(np.floor(np.nanmin(rho)*10)/10)))
    mae = points.loc[points.metric=='mae',['mean','ci_high']].to_numpy()
    max_mae = float(np.ceil(max(1.,float(np.nanmax(mae))*1.10)*4)/4)
    rho_ticks = ([rank_left] if rank_left<=-.5 else [])+[0,.5,1]
    specs = [('spearman','Spearman ρ ↑',(rank_left,1),rho_ticks),
        ('top10_ap','High-RMSF average precision ↑',(0,1),[0,.5,1]),
        ('mae','RMSF MAE (Å) ↓',(0,max_mae),[0,.5,1.0])]
    fig,axes = plt.subplots(2,3,figsize=(MAIN_SIZE_MM[0]/25.4,MAIN_SIZE_MM[1]/25.4),sharey=True)
    fig.subplots_adjust(left=.257,right=.985,bottom=.15,top=.92,wspace=.22,hspace=.42)
    for ri,group in enumerate(GROUPS):
        position = axes[ri,0].get_position()
        header_y = position.y1+.032
        fig.text(.015,header_y,TITLES[group],fontsize=7,fontweight='bold',color=INK)
        fig.text(.985,header_y,f'n = {data["counts"][group]:,} records · {data["families"][group]:,} families',
                 fontsize=6.5,ha='right',color=INK)
        for ci,(metric,label,limits,ticks) in enumerate(specs):
            ax = axes[ri,ci]
            for yi,name in enumerate(ORDER):
                if metric=='mae' and name in RELATIVE:
                    ax.text(.97,yi,'relative score',transform=ax.get_yaxis_transform(),ha='right',va='center',
                            fontsize=6.3,color='#6C797F')
                    continue
                row = points.loc[(points.method==name)&(points.group==group)&(points.metric==metric)].iloc[0]
                color,marker,size = visual(name,data['meta'])
                if np.isfinite(row.ci_low) and np.isfinite(row.ci_high):
                    ax.plot([row.ci_low,row.ci_high],[yi,yi],color=color,lw=.9,solid_capstyle='round',zorder=2)
                ax.plot(row['mean'],yi,marker=marker,ms=size,color=color,mec='white',mew=.3,ls='',zorder=3)
            ax.axhspan(-.45,.45,color=BLUE,alpha=.055,lw=0,zorder=0)
            ax.set_ylim(len(ORDER)-.45,-.65)
            ax.set_yticks(range(len(ORDER)),[DISPLAY.get(name,name) for name in ORDER],fontsize=6.5)
            ax.tick_params(axis='y',length=0,pad=5,labelleft=ci==0)
            ax.tick_params(axis='x',labelsize=6.3,length=2.5,pad=2)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_color('#80909A')
            ax.set_xlim(*limits);ax.set_xticks(ticks)
            ax.set_xlabel(label,fontsize=6.7,labelpad=4)
            ax.text(-.06,1.0,chr(97+ri*3+ci),transform=ax.transAxes,
                    fontsize=8.5,fontweight='bold',color=INK)
    handles = [Line2D([],[],color=BLUE,marker='o',lw=0,ms=3.2,label='FlexPPI configurations'),
        Line2D([],[],color=OCHRE,marker='o',lw=0,ms=3.2,label='Published learned models'),
        Line2D([],[],color=TEAL,marker='o',lw=0,ms=3.2,label='Graph / simple controls')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.007),ncol=3,fontsize=6.3,
               handletextpad=.45,columnspacing=1.6)
    fig.text(.015,.066,'14 completed entries · same 3,063 records · mixed training exposure',fontsize=6.3,color=INK)
    fig.text(.015,.043,'Points: mean across records; intervals: 95% family-cluster bootstrap.',fontsize=6.3,color=INK)
    save_exact_canvas(fig,OUT/'benchmark_current')
    plt.close(fig)
    points.to_csv(OUT/'benchmark_current_source.csv',index=False)
    return dict(spearman_mean_x_limits=[rank_left,1],absolute_mae_x_limits=[0,max_mae])


def draw_ecdf(data):
    per = data['per']
    selected_values=per.loc[per.method.isin(ECDF_METHODS),'spearman'].to_numpy()
    ecdf_left=max(-1.,min(-.2,float(np.floor(selected_values.min()*10)/10)))
    fig,axes = plt.subplots(1,2,figsize=(ECDF_SIZE_MM[0]/25.4,ECDF_SIZE_MM[1]/25.4))
    fig.subplots_adjust(left=.10,right=.985,bottom=.31,top=.87,wspace=.23)
    for ci,group in enumerate(GROUPS):
        ax = axes[ci]
        for name in ECDF_METHODS:
            values = np.sort(per.loc[(per.method==name)&(per.group==group),'spearman'].to_numpy())
            require(len(values)==data['counts'][group], 'ECDF group support changed')
            # Endpoints extend the empirical CDF over its domain; no records are added.
            xs = np.r_[ecdf_left,values,1.]
            ys = np.r_[0.,np.arange(1,len(values)+1)/len(values),1.]
            ax.step(xs,ys,where='post',color=visual(name,data['meta'])[0],lw=1.05,
                    ls=ECDF_LINE_STYLES[name],label=DISPLAY.get(name,name))
        ax.set(xlim=(ecdf_left,1),ylim=(0,1),
            xticks=([ecdf_left] if ecdf_left<=-.5 else [])+[0,.5,1],yticks=[0,.5,1],
            xlabel='Per-record Spearman ρ',title=TITLES[group]+f' (n = {data["counts"][group]:,})')
        if ci==0:
            ax.set_ylabel('Fraction of records',fontsize=6.7)
        ax.tick_params(labelsize=6.3,length=2.5,pad=2)
        ax.title.set_fontsize(7);ax.xaxis.label.set_fontsize(6.7)
        ax.text(-.12,1.10,chr(97+ci),transform=ax.transAxes,fontsize=8.5,fontweight='bold',color=INK)
    handles,labels = axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.066),ncol=3,fontsize=6.3,
               handlelength=3.5,columnspacing=1.6,labelspacing=.7)
    fig.text(.015,.018,'All records · six predeclared methods · mixed training exposure',fontsize=6.3,color=INK)
    save_exact_canvas(fig,OUT/'benchmark_distributions_current')
    plt.close(fig)
    return dict(ecdf_x_limits=[ecdf_left,1])


def main():
    parser = argparse.ArgumentParser(description='Render the frozen real protein/complex14-entry benchmark only.')
    parser.add_argument('--check-only',action='store_true',help='Validate real inputs without rendering or writing outputs')
    args = parser.parse_args()
    data = load_inputs()
    if args.check_only:
        print(json.dumps(dict(status='real_inputs_valid_no_render_written',methods=ORDER,
            primary_groups=data['counts'],membership_sha256=data['audit']['membership_sha256'])),flush=True)
        return
    OUT.mkdir(exist_ok=True)
    bounds = {**draw_main(data),**draw_ecdf(data)}
    require(all(digest(Path(p))==h for p,h in data['hashes'].items()), 'Input/source changed during rendering; output is not accepted')
    artifacts = [OUT/(stem+'.'+suffix) for stem in ('benchmark_current','benchmark_distributions_current')
                 for suffix in ('pdf','png','svg')]+[OUT/'benchmark_current_source.csv']
    png_sizes = {}
    for stem in ('benchmark_current','benchmark_distributions_current'):
        with (OUT/(stem+'.png')).open('rb') as handle:
            header = handle.read(24)
        require(header[:8]==b'\x89PNG\r\n\x1a\n', 'PNG export header invalid')
        png_sizes[stem+'.png'] = list(struct.unpack('>II',header[16:24]))
    report = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
        plot_script_sha256=digest(__file__),summary_sha256=digest(ROOT/'evaluation/group_summary.csv'),
        per_record_metrics_sha256=digest(ROOT/'evaluation/per_record_metrics.csv'),
        coverage_summary_sha256=digest(ROOT/'evaluation/coverage_summary.csv'),
        scoring_audit_sha256=digest(ROOT/'evaluation/scoring_audit.json'),
        membership_path=str(data['membership']),membership_sha256=digest(data['membership']),
        membership_group_column=data['member_group_column'],
        training_exposure_audit_path=str(data['exposure']),training_exposure_audit_sha256=digest(data['exposure']),
        input_sha256=data['hashes'],methods=ORDER,records=3063,valid_rmsf_positions=961321,
        primary_groups=data['counts'],group_residues=data['residues'],group_families=data['families'],
        group_definition='Proteins include homomers; protein complexes contain distinct biological protein partners. Membership is the separately frozen audited input.',
        style_reference='Science AlphaMissense Fig2A-C: fixed named method rows with aligned benchmark metrics; adapted to a two-row dot/interval quantitative grid.',
        layout=dict(main_rows=2,main_columns=3,main_row_groups=GROUPS,main_column_metrics=METRICS,
            method_rows_per_axes=14,main_canvas_mm=MAIN_SIZE_MM,ecdf_panels=2,ecdf_canvas_mm=ECDF_SIZE_MM,
            exact_canvas_export=True,export_png_dpi=300,pdf_fonttype=42,svg_fonttype='none',font_size_floor_pt=6.3),
        main_panel_mapping={chr(97+ri*3+ci):dict(group=group,metric=metric)
            for ri,group in enumerate(GROUPS) for ci,metric in enumerate(METRICS)},
        semantic_colors=dict(FlexPPI=BLUE,learned_external=OCHRE,graph_simple=TEAL),
        relative_output_methods=sorted(RELATIVE),relative_output_annotations=data['relative_cells'],
        numeric_points_in_main=78,all_observations_in_ecdf=True,ecdf_methods=ECDF_METHODS,
        ecdf_line_styles=ECDF_LINE_STYLES,ecdf_ANM_style_name='dash-dot-dot',
        ecdf_observation_counts={name:data['counts'] for name in ECDF_METHODS},
        confidence_intervals='Record-weighted means and 95% family-cluster bootstrap intervals, read from the real aggregate; no statistical recomputation by this script.',
        training_exposure='Mixed training conditions; own OOF does not imply sequence-family heldout. Family bootstrap is not training isolation.',
        absolute_mae_missing='FastProtFlex/ANM/GNM relative outputs have no Angstrom MAE; no fitted rescaling.',
        no_score_based_selection=True,unfinished_methods_in_primary_figure=False,
        excluded_noncomplete_methods=data['excluded_methods'],progress_snapshot=True,
        progress_scope='All14 completed entries on the full frozen cohort; deferred heavy methods remain outside this snapshot.',
        visual_review_status='pending_actual_review_of_current_output_hashes',
        source_data='benchmark_current_source.csv contains78 real summary cells;6 unit-inapplicable annotations are explicit in this audit.',
        artifact_sha256={p.name:digest(p) for p in artifacts},png_pixel_dimensions=png_sizes,**bounds)
    (OUT/'plot_audit.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    print('Rendered14 complete method/configuration entries in2 audited biological groups',flush=True)


if __name__=='__main__':
    main()
