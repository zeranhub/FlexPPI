# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Molecular context -> existing manuscript PyMOL coordinate views -> param inherit
# (b) Absolute peptide profiles -> LineTrend/plot_trend.py -> param inherit
# (c) All same-AA interface entropy pairs -> PairedBoxScatter visual logic -> param inherit
# (d) Complete-region MD comparison -> atlas_style.delta_comparison -> param inherit
# RULE: no native-run panels; every defined source position retained, no performance selection.
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
import sys,json
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parent.parent
sys.path.insert(0,str(ROOT.parent/'figure_design'))
from atlas_style import *
import pandas as pd

CAT=json.loads((ROOT.parent/'catalog/candidate_inventory.json').read_text(encoding='utf8'))
TITLES={
 'C18':'Influenza M1 recognition by different TCRs',
 'C19':'Recognition of the p53 R175H neoantigen',
 'C20':'HCMV NLV recognition by C7 and C25',
 'C21':'A6 recognition of Tax variants and cross-reactive antigens',
 'C22':'HBV peptide variants presented by HLA-E',
 'C23':'gp100 peptide anchor substitutions',
 'C24':'TCR55 recognition of HIV Pol and Pep20',
 'C25':'MAGE-A3 and Titin cross-recognition',
 'C26':'HIV SL9 variants recognized by the 868 TCR',
 'C27':'1E6 recognition of self and altered peptide ligands',
 'C28':'Public and private TCR structural comparison',
 'C29':'P14 recognition of viral and self-derived peptides',
 'C30':'MART-1 recognition by DMF5 and α24β17',
}
MANIFEST=[]

def read_pair(case,a,b):
    return pd.read_csv(ROOT/case['case_id']/'data'/f'{a}_{b}_all_chain_alignment.csv',dtype={'id_a':str,'id_b':str})
def select(frame,metric):
    if metric=='rmsf':return frame[frame.role.eq('peptide')&frame.paired&frame.rmsf_valid].copy()
    return frame[frame.role.isin(['mhc','peptide'])&frame.interface_union&frame.paired&frame.hchi_valid&frame.same_aa].copy()

def main():
    for case in CAT['cases']:
        cid=case['case_id']
        if cid not in TITLES:continue
        out=ROOT/cid/'figure';out.mkdir(exist_ok=True)
        a,b=case['contrasts'][0]
        f=read_pair(case,a,b)
        p=f[f.role.eq('peptide')].copy()
        q=select(f,'hchi')
        q['interface_label']=[('H' if role=='mhc' else 'P')+rid for role,rid in zip(q.role,q.id_a)]
        assert len(p) and len(q)
        fig=make_canvas(200)
        heading(fig,'a',TITLES[cid],.975)
        assets=[]
        for path,rect,title,col in [
            (ROOT/cid/'structures'/f'{a}_overview.png',[.055,.779,.27,.17],a.upper(),FIRST),
            (ROOT/cid/'structures'/f'{b}_overview.png',[.355,.779,.27,.17],b.upper(),SECOND),
            (ROOT/cid/'structures/all_state_peptide_superposition.png',[.675,.800,.27,.11],'Aligned peptides',GRAY),
        ]:
            _,r=add_molecule(fig,path,rect,title,col);assets.append(r)
        fig.text(.68,.779,f"All {len(case['pdbs'])} available structures",color=GRAY,fontsize=9)
        fig.text(.065,.751,f'Primary comparison: {a.upper()} → {b.upper()}',fontsize=9)
        fig.text(.065,.730,'MHC aligned; existing starting coordinates',fontsize=9,color=GRAY)
        heading(fig,'b','Model: RMSF across the complete aligned peptide',.690)
        legend(fig,[a.upper(),b.upper()],anchor=(.965,.718))
        ax=fig.add_axes([.11,.525,.84,.133])
        x=np.arange(len(p))
        for suffix,col,marker,face in [('a',FIRST,'o','white'),('b',SECOND,'s',SECOND)]:
            y=p['rmsf_full_'+suffix].to_numpy(float)
            ax.plot(x,y,color=col,marker=marker,ms=3.5,mfc=face,mec=col,lw=1.1)
        ymax=float(np.nanmax(p[['rmsf_full_a','rmsf_full_b']].to_numpy(float)))
        ytop=max(.4,np.ceil(ymax*1.25*5)/5)
        ax.set(ylabel='RMSF (Å)',ylim=(0,ytop),xlim=(-.4,len(p)-.6),xticks=x)
        lbl=[]
        for row in p.itertuples():
            rid=row.id_a if not row.gap_a else row.id_b
            if row.same_aa:label=f'{row.aa_a}{rid}'
            else:label=f'{row.aa_a if not row.gap_a else "–"}/{row.aa_b if not row.gap_b else "–"}{rid}'
            lbl.append(label)
        ax.set_xticklabels(lbl,rotation=0 if len(p)<14 else 45)
        ax.set_xlabel('Aligned peptide position · altered chemistry marked by /',labelpad=6)
        for tick,row in zip(ax.get_xticklabels(),p.itertuples()):
            if not row.same_aa:tick.set_color(GRAY)
        heading(fig,'c','Model: entropy across the complete comparable pMHC interface',.449)
        ax=fig.add_axes([.11,.306,.84,.111])
        hchi_records=paired_residue_values(ax,q,'hchi_full_a','hchi_full_b',q.interface_label.tolist(),rotate=60)
        fig.text(.53,.240,f'{len(q)} same-amino-acid sites · H: HLA/MHC; P: peptide',ha='center',fontsize=9,color=GRAY)
        heading(fig,'d','Existing MD: complete-region comparisons for all listed contrasts',.209)
        ds={};plotted=[]
        for metric,rect,unit,title in [('rmsf',[.115,.058,.335,.108],'Å','Peptide backbone'),('hchi',[.61,.058,.335,.108],'nats','Same-AA pMHC interface')]:
            data=[]
            for pa,pb in case['contrasts']:
                g=select(read_pair(case,pa,pb),metric)
                g['contrast']=pa+'>'+pb
                data.append(g)
            g=pd.concat(data,ignore_index=True)
            ax=fig.add_axes(rect)
            ds[metric]=delta_comparison(ax,g[metric+'_md_delta'],g[metric+'_full_delta'],unit=unit,title=title+f' · n={len(g)}')
            ds[metric]['contrasts']=case['contrasts'];ds[metric]['no_position_selection']=True
            g.to_csv(out/f'panel_d_{metric}_all_contrasts.csv',index=False,encoding='utf-8-sig')
            plotted.append(g)
        p.to_csv(out/'panel_b_complete_aligned_peptide.csv',index=False,encoding='utf-8-sig')
        q.to_csv(out/'panel_c_complete_comparable_pMHC_interface.csv',index=False,encoding='utf-8-sig')
        stem=out/f'{cid}_case_figure'
        audit=export_audited(fig,stem,extra=dict(case_id=cid,primary_pair=[a,b],tier=case['tier'],assets=assets,
             RMSF_rows=len(p),Hchi_rows=len(q),hchi_all_values=hchi_records,MD_comparison=ds,
             source_scope='Complete predefined peptide and same-AA union pMHC interface; gaps/undefined chi retained in full source CSV',
             no_invented_mechanism=True,all_data_positions_retained=True))
        color_order=[a,b]+[pdb for pdb in case['pdbs'] if pdb not in {a,b}]
        color_names=['橙','紫','蓝','绿','红','灰','浅蓝','棕','浅紫']
        color_key='、'.join(pdb.upper()+'为'+color_names[i] for i,pdb in enumerate(color_order))
        caption=(f"{cid} {' / '.join(z.upper() for z in case['pdbs'])}。a，真实既有MD起始结构的MHC配准视图，左、中为{a.upper()}和{b.upper()}，右为本系列全部可用肽结构叠合；叠合配色：{color_key}。b，主对照完整肽段的两状态模型RMSF，斜杠标注氨基酸改变，缺失对应保留为断点。c，主对照两状态8 Å界面并集中全部可自然对应、氨基酸相同且χ有定义的pMHC侧链模型Hχ；H/P分别表示MHC/肽，完整无χ或替换位置记录于源表。d，本系列所有列明生物对照在完整肽主链及同AA pMHC界面侧链上的模型变化与既有MD变化；两者变化均为后态减前态。点为残基位置，并非独立实验重复；ITrace各状态3次MD，未对残基或帧作显著性检验。具体支持与局限见正文。")
        (out/'caption_zh.txt').write_text(caption,encoding='utf8')
        (out/'statistics_and_reproducibility.md').write_text(f"# {cid} statistics and reproducibility\n\nFixed existing joint Full fold0; ATLAS1866 + DynaRepo564 training records, ITrace0. No new inference, training or MD. Model curves are a single deterministic existing prediction with no uncertainty bars. Model Hchi unit nats and RMSF unit angstrom. Existing MD labels are means over 3 independent trajectories per ITrace state; source sample variances and nrep are in complete alignment CSV. Panel d pools all aligned residue comparisons across all listed contrasts; residue positions are not independent experimental samples. No significance tests, regression fits or multiple comparison procedures were performed. MAE is descriptive change error, not proof that MD is physical truth. All positions in the predefined scope are shown; unpaired gaps, altered chemistry and undefined chi are explicitly retained in source data. Structural images are real existing topology coordinates, MHC-aligned with no outlier rejection.\n\nSource audit: ../data_audit.json and data/*_source_audit.json.\n",encoding='utf8')
        MANIFEST.append(dict(case_id=cid,png=str(stem)+'.png',pdf=str(stem)+'.pdf',svg=str(stem)+'.svg',audit=str(stem)+'.audit.json',caption=str(out/'caption_zh.txt')))
        plt.close(fig)
        print(cid,'figure exported',len(p),len(q),flush=True)
    (ROOT/'general_figure_manifest.json').write_text(json.dumps(dict(status='code_checks_passed_visual_review_pending',figures=MANIFEST),ensure_ascii=False,indent=2),encoding='utf8')

if __name__=='__main__':main()
