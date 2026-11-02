# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) existing native PDF -> unchanged crop including b's independent method legend
# (c) original complete PR vertices / AUROC/plot_AUROC1.py -> param inherit
# (d) frozen external PR and ROC paired as one panel -> param inherit
# (e) interface-stratified native PDF -> native PDF reuse, d relabeled e; centered group legend
# (S3a,b) original complete PR vertices -> param inherit, original precision scales
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

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

from pathlib import Path
import hashlib
import json
import shutil
import xml.etree.ElementTree as ET
import fitz
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import Bbox
from matplotlib.ticker import FuncFormatter, NullLocator

OUT = Path(__file__).resolve().parent
ARCH = OUT / 'source_archive'
ORIGINAL = Path(r'E:\Project\FlexPPI\paper\figure2_compact_layouts_20261023')
REV = OUT.parent
EXTERNAL = REV / 'figure2_external_curves_20261024'
BOTTOM = REV / 'figure2_panel_d_20261024'
W, TOP_H, ROW_H = 183.0, 86.0, 65.0
BOTTOM_H = 52.0 * W / 150.0
H = TOP_H + ROW_H + BOTTOM_H
PT = 72.0 / 25.4
INTERNAL = ['Full', 'Intra', 'Static', 'AA']
INT_COL = dict(Full='#2166AC', Intra='#B87936', Static='#3A8D8C', AA='#999999')
INT_LS = dict(Full='-', Intra='--', Static='-.', AA=':')
INT_LW = dict(Full=1.25, Intra=1.1, Static=1.1, AA=1.0)
INT_Z = dict(Full=5, Intra=6, Static=4, AA=3)
EXT = [
    ('Full', 'FlexPPI', '#2166AC', '-', None),
    ('backflip-2.1 (assembly)', 'BackFlip-2.1', '#B35806', (0,(4,1.6)), 'o'),
    ('backflip-1.0 (assembly)', 'BackFlip-1.0', '#D6604D', (0,(4,1.4,1,1.4)), 's'),
    ('DynaProt-M', 'DynaProt-M', '#1B7837', (0,(1,1.5)), '^'),
    ('FastProtFlex', 'FastProtFlex', '#762A83', (0,(6,1.6)), 'D'),
    ('ANM', 'ANM', '#666666', (0,(2,1.4)), 'v'),
    ('GNM', 'GNM', '#999999', (0,(4,1.4,1,1.4,1,1.4)), 'P'),
]
SCALES = {
    'rmsf': ([0,.5,1], [0,.25,1], [0,.5,1]),
    'chi_entropy': ([0,.5,.9,1], [0,.10,.45,1], [0,.9,1]),
    'within_type_hchi': ([0,.45,.75,1], [0,.12,.90,1], [0,.6,1]),
}
mpl.rcParams.update({'font.size':9, 'axes.labelsize':9, 'axes.titlesize':10,
    'xtick.labelsize':8, 'ytick.labelsize':8, 'legend.fontsize':8,
    'axes.unicode_minus':False})

def archive_sources():
    # Self-contained archive: validated locally without requiring the original disk tree.
    if (OUT/'source_manifest.json').exists() and ARCH.exists():
        frozen=json.loads((OUT/'source_manifest.json').read_text(encoding='utf-8'))
        for item in frozen:
            src=ARCH/item['file']
            assert hashlib.sha256(src.read_bytes()).hexdigest()==item['sha256'],item['file']
        return
    ARCH.mkdir(exist_ok=True)
    pairs = {
        'original_fig2_compact.pdf': ORIGINAL/'fig2_compact_merged.pdf',
        'original_fig2_panel_d.pdf': BOTTOM/'fig2_panel_d.pdf',
        'original_fig2_panel_d_source.csv': BOTTOM/'fig2_panel_d_source.csv',
        'internal_PR_vertices.csv': ORIGINAL/'complete_PR_vertices.csv',
        'internal_summary.json': ORIGINAL/'source_archive'/'summary.json',
        'external_mean_curves.csv': EXTERNAL/'mean_curves.csv',
        'external_summary_metrics.csv': EXTERNAL/'summary_metrics.csv',
        'external_calculation_audit.json': EXTERNAL/'calculation_audit.json',
        'external_per_record_metrics.csv': EXTERNAL/'per_record_metrics.csv',
        'external_cohort_records.csv': EXTERNAL/'cohort_records.csv',
    }
    manifest=[]
    for name, src in pairs.items():
        dst=ARCH/name
        if not dst.exists(): shutil.copy2(src,dst)
        a=hashlib.sha256(src.read_bytes()).hexdigest()
        b=hashlib.sha256(dst.read_bytes()).hexdigest()
        assert a==b, (name,'source changed since first archive')
        manifest.append(dict(file=name,source=str(src),sha256=b))
    (OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')

def text(fig, height, x, top, s, **kw):
    return fig.text(x/W, 1-top/height, s, **kw)

def axis(fig, height, x, top, width, h):
    return fig.add_axes([x/W,1-(top+h)/height,width/W,h/height])

def piecewise(values, knots, heights):
    v=np.asarray(values,dtype=float)
    a=np.asarray(knots,dtype=float);b=np.asarray(heights,dtype=float)
    out=np.interp(v,a,b)
    return np.where(v<a[0],b[0]+(v-a[0])*(b[1]-b[0])/(a[1]-a[0]),
                    np.where(v>a[-1],b[-1]+(v-a[-1])*(b[-1]-b[-2])/(a[-1]-a[-2]),out))

def original_pr(fig, height, task, x, width, top, h, label, heading_top, metric_top):
    ax=axis(fig,height,x,top,width,h)
    part=vertices[vertices.task.eq(task)]
    means=summary['within_type'] if task=='within_type_hchi' else summary['metrics'][task]
    for m in INTERNAL:
        q=part[part.method.eq(m)]
        assert len(q)==102 and q.final_endpoint.sum()==1
        assert q.recall.iloc[-1]==1
        # The exact 101 mean-curve vertices plus the original all-selected endpoint.
        # No resampling, smoothing, recalibration or metric recomputation.
        ax.step(q.recall,q.precision,where='post',color=INT_COL[m],ls=INT_LS[m],
                lw=INT_LW[m],zorder=INT_Z[m])
    baseline=float(means['Full']['prevalence_mean']) if task!='within_type_hchi' else float(means['Full']['mean_record_type_prevalence'])
    ax.axhline(baseline,color='#AAAAAA',ls=(0,(2,2)),lw=.65,zorder=0)
    knots,heights,ticks=SCALES[task]
    forward=lambda v:piecewise(v,knots,heights)
    inverse=lambda v:piecewise(v,heights,knots)
    ax.set_yscale('function',functions=(forward,inverse))
    ax.set(xlim=(0,1),ylim=(0,float(inverse(1.025))),xlabel='Recall',ylabel='Precision')
    ax.set_xticks([0,.5,1]);ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v,pos:format(v,'g')))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.tick_params(pad=2,labelsize=9);ax.xaxis.labelpad=2.5
    if label:text(fig,height,x-11,heading_top,label,fontsize=11,fontweight='bold',va='bottom')
    heading={'rmsf':'RMSF','chi_entropy':'Hχ','within_type_hchi':'Hχ within type'}[task]
    text(fig,height,x+width/2,heading_top,heading,fontsize=10,ha='center',va='bottom')
    # User removed the entire AP/numeric annotation strip from the main figure.
    if task!='within_type_hchi':
        text(fig,height,x-10,metric_top,'AP',fontsize=8,ha='left',va='center',color='#444444')
        for k,m in enumerate(INTERNAL):
            text(fig,height,x+(k+.5)*width/4,metric_top,f'{means[m]["mean_ap"]:.3f}',
                 fontsize=8,ha='center',va='center',color=INT_COL[m])
    return ax

def scalar(method, metric):
    row=metrics[(metrics.method.eq(method))&(metrics.metric.eq(metric))]
    assert len(row)==1
    return float(row.iloc[0]['mean'])

def external_panel(fig,height,kind,x,label):
    ax=axis(fig,height,x,10,42,40)
    for method,name,color,ls,marker in reversed(EXT):
        part=curves[(curves.method.eq(method))&(curves.curve_type.eq(kind))].sort_values('x')
        xx,yy=part.x.to_numpy(),part['mean'].to_numpy()
        assert len(xx)==1001 and np.isfinite(yy).all() and yy.min()>=0 and yy.max()<=1
        assert np.allclose(xx,np.linspace(0,1,1001),atol=1e-15)
        if kind=='ROC':assert np.diff(yy).min()>=-1e-12 and yy[-1]==1
        ax.plot(xx,yy,color=color,ls=ls,lw=1.35 if method=='Full' else .95,
                marker=marker,markevery=[180,380,580,780] if marker else None,
                markersize=2.25,markerfacecolor='white',markeredgewidth=.55,
                zorder=10 if method=='Full' else 4)
        if kind=='PR':
            ax.plot([1,1],[yy[-1],prevalence],color=color,ls=ls,lw=.65,zorder=3,clip_on=False)
        else:
            ax.plot([0,0],[0,yy[0]],color=color,ls=ls,lw=.65,zorder=3,clip_on=False)
    if kind=='PR':
        ax.axhline(prevalence,color='#AAAAAA',ls=(0,(2,2)),lw=.65,zorder=0)
        ax.set(xlabel='Recall',ylabel='Precision')
    else:
        ax.plot([0,1],[0,1],color='#AAAAAA',ls=(0,(2,2)),lw=.65,zorder=0)
        ax.set(xlabel='False-positive rate',ylabel='True-positive rate')
    ax.set(xlim=(0,1),ylim=(0,1.025))
    ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
    ax.tick_params(pad=2,length=2.5);ax.xaxis.labelpad=2.5;ax.yaxis.labelpad=2.5
    if label:text(fig,height,x-11,5,label,fontsize=11,fontweight='bold',va='bottom')
    text(fig,height,x+21,5,'RMSF '+kind,fontsize=9.5,ha='center',va='bottom')
    if kind=='ROC':
        handles=[Line2D([],[],color=c,ls=ls,lw=1.35 if m=='Full' else .95,
            marker=mark,ms=2.25,mfc='white',mew=.55,label=name)
            for m,name,c,ls,mark in EXT]
        legend=ax.legend(handles=handles,loc='lower right',bbox_to_anchor=(.99,.01),
            fontsize=7.5,ncol=1,handlelength=1.35,handletextpad=.45,
            labelspacing=.10,borderpad=.1,borderaxespad=0,frameon=False)
        legend.get_texts()[0].set_fontweight('bold')
        LEGENDS['external']=legend
    return ax

def shared_external_table(fig,height):
    # One shared, aligned table gives both scalar metrics for seven curves.
    # Keep all line keys outside the plotting areas; do not repeat 7-item legends.
    text(fig,height,100,65,'Method',fontsize=7.8,ha='left',va='center')
    text(fig,height,154,65,'AP',fontsize=7.8,ha='center',va='center')
    text(fig,height,174,65,'AUROC',fontsize=7.8,ha='center',va='center')
    for i,(m,name,color,ls,marker) in enumerate(EXT):
        y=68+i*3.1
        line=Line2D([88/W,96/W],[1-y/height]*2,transform=fig.transFigure,
            color=color,ls=ls,lw=1.35 if m=='Full' else .95,marker=None)
        fig.add_artist(line)
        if marker:fig.add_artist(Line2D([92/W],[1-y/height],transform=fig.transFigure,
            color=color,marker=marker,ms=2.25,mfc='white',mew=.55,ls='none'))
        weight='bold' if m=='Full' else 'normal'
        text(fig,height,100,y,name,fontsize=7.5,ha='left',va='center',fontweight=weight)
        text(fig,height,154,y,f'{scalar(m,"AP"):.3f}',fontsize=7.5,ha='center',va='center',fontweight=weight)
        text(fig,height,174,y,f'{scalar(m,"AUROC"):.3f}',fontsize=7.5,ha='center',va='center',fontweight=weight)

def save_plot(fig,stem,height):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer();outside=[]
    for t in fig.findobj(mpl.text.Text):
        if not t.get_visible() or not t.get_text().strip():continue
        b=t.get_window_extent(renderer)
        if b.x0<-.75 or b.y0<-.75 or b.x1>fig.bbox.width+.75 or b.y1>fig.bbox.height+.75:
            outside.append(t.get_text())
    assert not outside,outside
    box=Bbox.from_bounds(0,0,W/25.4,height/25.4)
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{stem}.{ext}',bbox_inches=box,dpi=450,facecolor='white')
    plt.close(fig)
    return dict(stem=stem,height_mm=height,outside_text=outside)

def compose_row():
    fig=plt.figure(figsize=(W/25.4,ROW_H/25.4))
    ax=original_pr(fig,ROW_H,'within_type_hchi',15,42,10,40,'c',5,11)
    handles=[Line2D([],[],color=INT_COL[m],ls=INT_LS[m],lw=1.4,label=m) for m in INTERNAL]
    LEGENDS['internal']=ax.legend(handles=handles,loc='upper right',bbox_to_anchor=(.99,1.0),
        fontsize=9,ncol=1,handlelength=2,handletextpad=.4,labelspacing=.10,
        borderpad=.15,borderaxespad=0,frameon=False)
    # Cohort labels and explanatory text are documented in the caption instead.
    # Pair PR/ROC more tightly than the c-to-d separation to show one group.
    external_panel(fig,ROW_H,'PR',81,'d')
    external_panel(fig,ROW_H,'ROC',137,'')
    # AP/AUROC values are archived unchanged and moved to supplementary Table S10.
    fig.canvas.draw()
    for kind,legend in LEGENDS.items():
        box=legend.get_window_extent(fig.canvas.get_renderer())
        ax=legend.axes
        assert (box.x0>=ax.bbox.x0-1e-6 and box.y0>=ax.bbox.y0-1e-6
            and box.x1<=ax.bbox.x1+1e-6 and box.y1<=ax.bbox.y1+1e-6),(kind,box)
    return save_plot(fig,'fig2_second_row',ROW_H)

def compose_supplement():
    height=89.
    fig=plt.figure(figsize=(W/25.4,height/25.4))
    original_pr(fig,height,'rmsf',25,52,18,48,'a',6,12)
    original_pr(fig,height,'chi_entropy',119,52,18,48,'b',6,12)
    handles=[Line2D([],[],color=INT_COL[m],ls=INT_LS[m],lw=INT_LW[m],label=m) for m in INTERNAL]
    fig.legend(handles=handles,bbox_to_anchor=(.5,1-81/height),loc='center',
        ncol=4,fontsize=9,columnspacing=1.5,handlelength=2)
    return save_plot(fig,'figS3_internal_PR',height)

def native_assembly():
    top=fitz.open(ARCH/'original_fig2_compact.pdf')
    bottom=fitz.open(ARCH/'original_fig2_panel_d.pdf')
    row=fitz.open(OUT/'fig2_second_row.pdf')
    # Keep a/b completely unchanged, including b's independent 4-method key.
    # Remove only the precisely identified old 3-key bar legend drawing block.
    bcontents=bottom[0].get_contents()
    assert len(bcontents)==1
    bstream=bottom[0].read_contents()
    anchor=(b'Q q 0.7254901961 0.8078431373 0.8901960784 rg 0 j 0.25 w 0.168627451\n'
        b'0.2901960784 0.3882352941 RG 0.7254901961 0.8078431373 0.8901960784 rg\n\n'
        b'62.638113 137.536575 m')
    start=bstream.rfind(anchor)
    assert start>0 and bstream.count(anchor)==1
    block=bstream[start:]
    assert block.count(b'\nBT\n')==3 and block.count(b'\nB\n')==3
    assert block.endswith(b'Q\nQ\n')
    bottom.update_stream(bcontents[0],bstream[:start]+b'Q\n')
    # Rename only the isolated panel letter in a temporary native PDF.
    span=next(s for b in bottom[0].get_text('dict')['blocks'] if 'lines' in b
        for l in b['lines'] for s in l['spans'] if s['text']=='d')
    r=fitz.Rect(span['bbox']);r.x0-=.25;r.x1+=.25
    # The font's theoretical descender bbox overlaps the 1.0 tick's ascender bbox;
    # d has no descender ink. Keep the isolated redaction safely above that tick.
    r.y1=span['origin'][1]+.75
    bottom[0].add_redact_annot(r,fill=(1,1,1))
    bottom[0].apply_redactions(images=0,graphics=0,text=0)
    bottom[0].insert_font(fontname='relayoutArialBold',fontfile=r'C:\Windows\Fonts\arialbd.ttf')
    bottom[0].insert_text(span['origin'],'e',fontname='relayoutArialBold',fontsize=span['size'],
        color=(31/255,44/255,51/255))
    # Reposition a concise shared group key, centered just above all three bars.
    legfig=plt.figure(figsize=(W/25.4,BOTTOM_H/25.4))
    handles=[Patch(facecolor=col,edgecolor='#2B4A63',lw=.25*1.22,label=name)
        for name,col in [('1 interface','#B9CEE3'),('2–3 interfaces','#5E8FBF'),('≥4 interfaces','#2060A8')]]
    legfig.legend(handles=handles,loc='center',bbox_to_anchor=(.5,1-9/BOTTOM_H),ncol=3,
        fontsize=6.1*1.22,handlelength=.85,handletextpad=.4,columnspacing=2.4,frameon=False)
    box=Bbox.from_bounds(0,0,W/25.4,BOTTOM_H/25.4)
    legfig.savefig(OUT/'fig2_bottom_legend.pdf',bbox_inches=box,transparent=True)
    plt.close(legfig)
    legend=fitz.open(OUT/'fig2_bottom_legend.pdf')
    doc=fitz.open();p=doc.new_page(width=W*PT,height=H*PT)
    top_clip=fitz.Rect(0,0,W*PT,TOP_H*PT)
    p.show_pdf_page(top_clip,top,0,clip=top_clip)
    p.show_pdf_page(fitz.Rect(0,TOP_H*PT,W*PT,(TOP_H+ROW_H)*PT),row,0)
    p.show_pdf_page(fitz.Rect(0,(TOP_H+ROW_H)*PT,W*PT,H*PT),bottom,0)
    p.show_pdf_page(fitz.Rect(0,(TOP_H+ROW_H)*PT,W*PT,H*PT),legend,0)
    doc.subset_fonts()
    doc.save(OUT/'fig2_grouped.pdf',garbage=4,deflate=True)
    doc.close();top.close();row.close();bottom.close();legend.close()
    d=fitz.open(OUT/'fig2_grouped.pdf');p=d[0]
    p.get_pixmap(dpi=450,alpha=False).save(OUT/'fig2_grouped.png')
    p.get_pixmap(dpi=160,alpha=False).save(OUT/'fig2_grouped_preview.png')
    svg=p.get_svg_image(text_as_path=False)
    (OUT/'fig2_grouped.svg').write_text(svg,encoding='utf-8')
    images=p.get_images(full=True)
    assert len(images)==0,images
    assert '<image' not in svg
    assert abs(p.rect.width/PT-W)<1e-4 and abs(p.rect.height/PT-H)<1e-4
    visible=p.get_text(clip=p.rect)
    for letter in 'abcde': assert '\n'+letter+'\n' in '\n'+visible+'\n',letter
    fonts=p.get_fonts(full=True)
    d.close()
    return dict(width_mm=W,height_mm=H,raster_images_in_pdf=len(images),raster_images_in_svg=0,
        top_native_crop_mm=[0,0,W,TOP_H],bottom_scaling=W/150,
        font_names=sorted(set(f[3] for f in fonts)),all_panels_native_vector=True)

def qa_statistics(assembled,layouts):
    counts={}
    for task in ['rmsf','chi_entropy','within_type_hchi']:
        q=vertices[vertices.task.eq(task)]
        assert len(q)==408 and set(q.method)==set(INTERNAL)
        counts[task]=len(q)
    assert audit['n_records']==615 and audit['n_families']==261 and audit['n_residues']==260340
    assert len(curves)==14014 and len(metrics)==14
    assert set(metrics.n_records)=={615}
    assert abs(scalar('Full','AP')-.6589512895302075)<1e-14
    assert abs(scalar('Full','AUROC')-.9173297776735228)<1e-14
    assert metrics['mean'].between(0,1).all()
    report=dict(status='code_data_passed_visual_review_pending',native_assembly=assembled,
        layout=layouts,main_letters=list('abcde'),second_row=['c: within-type Hchi PR','d: paired external RMSF PR and ROC'],
        internal_vertices_preserved=counts,external_mean_vertices=14014,
        data_recomputed=False,models_retrained=False,source_data_dropped=False,
        original_top_a_b_preserved_with_independent_b_legend=True,original_bottom_three_metrics_preserved=True,
        grouped_external_curves=True,bottom_group_label='e',bottom_legend='3 short names centered above all bars at local y=9mm',
        removed_annotations=['c AP annotation and four scalar values','external cohort subtitles',
            'external AP/AUROC table','c Internal controls and record count','bottom sample/family count suffixes'],
        concise_external_titles=['RMSF PR','RMSF ROC'],
        main_internal_precision_scale='Original continuous piecewise-linear precision scale, identical knots and tick values',
        external_precision_and_ROC_scales='Linear 0-1; no smoothing',
        supplement='S3 unchanged; external AP/AUROC moved to Table S10 by document editor',
        external_metric_policy='Exact per-record AP/AUROC averaged equally across records; not area under displayed mean curve',
        known_limitations=['Original internal PR precision-axis transformation is retained at user request; external PR/ROC axes are linear.',
            'A and B retain frozen native PDF content; numerical prediction evaluation is not rerun.',
            'No new significance tests are added; source confidence intervals remain in archived tables.'])
    (OUT/'figure_QA.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    stats='''# Figure 2 relayout: statistics and reproducibility

The figure uses the frozen original data. No model is retrained, no data point is dropped, and no model metric is recalculated. Main-panel a and b are reused directly from the original native PDF at the same 183 mm physical width, including b's independent four-method legend. The external PR and ROC axes share one panel d; PR moves 5 mm closer to ROC while all axis sizes, data, fonts and scales are retained. The three interface-stratified bar charts are reused directly from their original native PDF with the same 183/150 scaling used by the manuscript composite; the bottom group is labeled e. Its concise three-item group legend is centered 9 mm below the top of the bottom group, closer to the bar axes and separate from the e letter. Counts remain in the caption. Internal c and external d-ROC legends occupy blank plot regions. Unnecessary AP/cohort/count annotations remain removed, and external axis titles are RMSF PR and RMSF ROC.

Main c: within-type Hchi retrieval, 3,063 records. Four controls (Full, Intra, Static, AA) retain all 102 original PR vertices each, including the all-selected endpoint. The original continuous piecewise-linear precision-axis transformation is retained with knots [0,0.45,0.75,1], display heights [0,0.12,0.90,1], ticks [0,0.6,1]. AP is the frozen equal-record mean (Full 0.623, Intra 0.609, Static 0.555, AA 0.146); curves are not integrated to replace AP. Each eligible record/type is evaluated before aggregation as documented by the original calculation. No confidence band is added.

Main d (paired PR and ROC): seven methods, identical 615 DynaRepo protein complexes, 261 families and 260,340 residues. Positives are residues at or above each record's true RMSF 90th percentile; tied values are retained. Curves are equal-record means on the frozen 1,001-point common recall/FPR grid, without smoothing. PR retains the first threshold reaching recall one and a separate all-selected endpoint; ROC retains empirical vertical-start segments and the corrected interpolation that respects ordered duplicate FPR vertices. The internal ROC legend contains method names only. AP/AUROC numbers and their original 95% intervals are moved to supplementary Table S10 by the document editor. They remain independently averaged exact per-record metrics, not the area of the display curve. No curve confidence band is shown. Frozen scalar confidence intervals in the source summary use 2,000 family-cluster bootstrap draws, seed 20261024, percentile 95% intervals.

Main e: all original seven-method grouped bars for Spearman rho, high-RMSF AP and standardized MAE remain. Strata contain 414 records/196 families (one interface), 144/58 (two to three), 57/34 (at least four). Means are equal-record averages. Error bars are original 95% family-cluster bootstrap percentile intervals from 2,000 draws, seed 20261024. Family memberships can span strata. No new hypothesis test or multiplicity adjustment is introduced.

Supplement S3a/b: original RMSF and overall Hchi internal-control PR panels. Each panel preserves all 408 vertices, all four methods and frozen AP. Precision axes retain the original task-specific continuous piecewise-linear transformations: RMSF knots [0,0.5,1] to display heights [0,0.25,1]; overall Hchi knots [0,0.5,0.9,1] to [0,0.10,0.45,1]. These unequal axis spacings must be disclosed in the caption; AP remains the original raw-threshold calculation.

Traceability: source_manifest.json records absolute original paths and SHA-256 hashes; source_archive contains the exact frozen input files needed to reproduce every newly rendered curve and reuse the native panels. Reproduction: python grouped_layout_fig2.py. The delivered PDF and SVG contain vector paths and text, zero raster image objects. PNG is a 450 dpi preview. Document insertion and wording are handled separately. Supplement S3 is unchanged and is not regenerated.
'''
    (OUT/'statistics_and_reproducibility.md').write_text(stats,encoding='utf-8')
    return report

if __name__=='__main__':
    archive_sources()
    vertices=pd.read_csv(ARCH/'internal_PR_vertices.csv',float_precision='round_trip')
    summary=json.loads((ARCH/'internal_summary.json').read_text(encoding='utf-8'))
    curves=pd.read_csv(ARCH/'external_mean_curves.csv',float_precision='round_trip')
    curves=curves[curves.scope.eq('all')]
    metrics=pd.read_csv(ARCH/'external_summary_metrics.csv',float_precision='round_trip')
    metrics=metrics[metrics.scope.eq('all')]
    audit=json.loads((ARCH/'external_calculation_audit.json').read_text(encoding='utf-8'))
    prevalence=float(audit['all_selected_precision'])
    LEGENDS={}
    layouts=[compose_row()]
    assembled=native_assembly()
    report=qa_statistics(assembled,layouts)
    print(json.dumps({'outputs':str(OUT),'native_assembly':assembled,'layouts':layouts},indent=2))
